use std::collections::HashMap;
use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::utils::DynError;

const AFFINITY_FILE_BASENAME: &str = "yi-corr-affinity.txt";

fn home_dir() -> Option<PathBuf> {
    if let Some(v) = env::var_os("HOME") {
        return Some(PathBuf::from(v));
    }
    if let Some(v) = env::var_os("USERPROFILE") {
        return Some(PathBuf::from(v));
    }
    match (env::var_os("HOMEDRIVE"), env::var_os("HOMEPATH")) {
        (Some(d), Some(p)) => Some(PathBuf::from(d).join(p)),
        _ => None,
    }
}

fn cargo_tmp_dir() -> Option<PathBuf> {
    if let Some(v) = env::var_os("CARGO_HOME") {
        return Some(PathBuf::from(v).join("tmp"));
    }
    home_dir().map(|h| h.join(".cargo").join("tmp"))
}

fn affinity_file_path() -> Option<PathBuf> {
    cargo_tmp_dir().map(|d| d.join(AFFINITY_FILE_BASENAME))
}

fn ensure_parent_dir(path: &Path) {
    if let Some(parent) = path.parent() {
        let _ = fs::create_dir_all(parent);
    }
}

fn parse_affinity_spec(spec: &str) -> Result<Vec<usize>, DynError> {
    let mut out = Vec::<usize>::new();
    for token in spec.split(|c: char| c == ',' || c.is_ascii_whitespace()) {
        let t = token.trim();
        if t.is_empty() {
            continue;
        }
        if let Some((a, b)) = t.split_once('-') {
            let start: usize = a.trim().parse()?;
            let end: usize = b.trim().parse()?;
            if start > end {
                return Err(format!("invalid core range: {start}-{end}").into());
            }
            out.extend(start..=end);
        } else {
            out.push(t.parse::<usize>()?);
        }
    }
    if out.is_empty() {
        return Err("affinity set is empty".into());
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

fn parse_affinity_file(path: &Path) -> Result<Vec<Vec<usize>>, DynError> {
    let text = fs::read_to_string(path)?;
    let mut specs = Vec::<Vec<usize>>::new();
    for line in text.lines() {
        let content = line.split('#').next().unwrap_or("").trim();
        if content.is_empty() {
            continue;
        }
        specs.push(parse_affinity_spec(content)?);
    }
    Ok(specs)
}

struct SlotLockGuard {
    path: PathBuf,
    _file: File,
}

impl Drop for SlotLockGuard {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

fn try_claim_slot_lock(lock_path: PathBuf) -> Option<SlotLockGuard> {
    let mut f = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&lock_path)
        .ok()?;
    let _ = writeln!(f, "pid={}", std::process::id());
    Some(SlotLockGuard {
        path: lock_path,
        _file: f,
    })
}

fn map_requested_to_core_ids(requested: &[usize]) -> Option<Vec<core_affinity::CoreId>> {
    let available = core_affinity::get_core_ids()?;
    if available.is_empty() {
        return None;
    }
    let mut by_id = HashMap::<usize, core_affinity::CoreId>::new();
    for c in available {
        by_id.insert(c.id, c);
    }
    let mut out = Vec::<core_affinity::CoreId>::new();
    for &id in requested {
        if let Some(c) = by_id.get(&id) {
            out.push(*c);
        }
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

pub struct AffinityRuntime {
    worker_cores: Option<Arc<Vec<core_affinity::CoreId>>>,
    _slot_lock: Option<SlotLockGuard>,
    info: Option<String>,
}

impl AffinityRuntime {
    pub fn from_default_file() -> Result<Self, DynError> {
        let Some(path) = affinity_file_path() else {
            return Ok(Self {
                worker_cores: None,
                _slot_lock: None,
                info: None,
            });
        };
        ensure_parent_dir(&path);
        if !path.exists() {
            return Ok(Self {
                worker_cores: None,
                _slot_lock: None,
                info: None,
            });
        }
        let groups = parse_affinity_file(&path)?;
        if groups.is_empty() {
            return Ok(Self {
                worker_cores: None,
                _slot_lock: None,
                info: Some(format!(
                    "Affinity file exists but has no valid entries: {}",
                    path.display()
                )),
            });
        }

        let mut slot_idx = 0usize;
        let mut lock_guard = None;
        if groups.len() > 1 {
            let lock_dir = path.parent().unwrap_or_else(|| Path::new("."));
            for idx in 0..groups.len() {
                let lock_path = lock_dir.join(format!("yi-corr-affinity.slot{idx}.lock"));
                if let Some(g) = try_claim_slot_lock(lock_path) {
                    slot_idx = idx;
                    lock_guard = Some(g);
                    break;
                }
            }
            if lock_guard.is_none() {
                slot_idx = 0;
            }
        }

        let selected = &groups[slot_idx];
        let mapped = map_requested_to_core_ids(selected);
        let info = if mapped.is_some() {
            let lock_msg = if groups.len() > 1 {
                if lock_guard.is_some() {
                    format!(" slot={slot_idx}")
                } else {
                    " slot=fallback(0)".to_string()
                }
            } else {
                "".to_string()
            };
            Some(format!(
                "CPU affinity enabled via {}{} cores={}",
                path.display(),
                lock_msg,
                selected
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ))
        } else {
            Some(format!(
                "CPU affinity config found in {}, but requested cores are unavailable on this host: {}",
                path.display(),
                selected
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ))
        };

        Ok(Self {
            worker_cores: mapped.map(Arc::new),
            _slot_lock: lock_guard,
            info,
        })
    }

    pub fn worker_cores(&self) -> Option<Arc<Vec<core_affinity::CoreId>>> {
        self.worker_cores.clone()
    }

    pub fn info(&self) -> Option<&str> {
        self.info.as_deref()
    }
}

pub fn reader_core_from_env() -> Result<Option<core_affinity::CoreId>, DynError> {
    let Some(raw) = env::var_os("YI_READER_CORE") else {
        return Ok(None);
    };
    let raw = raw.to_string_lossy();
    let requested = parse_affinity_spec(&raw)?;
    Ok(map_requested_to_core_ids(&requested).and_then(|v| v.first().copied()))
}

pub fn set_current_thread_core(core: core_affinity::CoreId) -> bool {
    core_affinity::set_for_current(core)
}

pub struct CpuAllocation {
    pub workers: Vec<core_affinity::CoreId>,
    pub io: core_affinity::CoreId,
    pub total: usize,
}

pub fn allocate_cpus(
    available: &[core_affinity::CoreId],
    requested: usize,
    preferred_io: Option<core_affinity::CoreId>,
) -> Result<CpuAllocation, DynError> {
    if requested == 0 || available.is_empty() {
        return Err("CPU allocation requires at least one available CPU".into());
    }
    let total = requested.min(available.len());
    let io = preferred_io.unwrap_or(available[total - 1]);
    if !available.iter().any(|c| c.id == io.id) {
        return Err("YI_READER_CORE is outside the allowed CPU affinity set".into());
    }
    let workers = if total == 1 {
        vec![io]
    } else {
        available
            .iter()
            .copied()
            .filter(|c| c.id != io.id)
            .take(total - 1)
            .collect()
    };
    Ok(CpuAllocation { workers, io, total })
}

// Pin only the file-output section. Restoring the mask before workflows spawn
// child processes prevents them from inheriting a one-CPU restriction.
#[cfg(target_os = "linux")]
pub struct IoAffinityGuard(libc::cpu_set_t);

#[cfg(target_os = "linux")]
impl IoAffinityGuard {
    pub fn enter(core: Option<core_affinity::CoreId>) -> Result<Option<Self>, DynError> {
        let Some(core) = core else {
            return Ok(None);
        };
        let mut previous = unsafe { std::mem::zeroed::<libc::cpu_set_t>() };
        let result = unsafe {
            libc::sched_getaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &mut previous)
        };
        if result != 0 {
            return Err(std::io::Error::last_os_error().into());
        }
        if !set_current_thread_core(core) {
            return Err("failed to pin output thread to I/O CPU".into());
        }
        Ok(Some(Self(previous)))
    }
}

#[cfg(target_os = "linux")]
impl Drop for IoAffinityGuard {
    fn drop(&mut self) {
        unsafe {
            libc::sched_setaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &self.0);
        }
    }
}

#[cfg(not(target_os = "linux"))]
pub struct IoAffinityGuard;
#[cfg(not(target_os = "linux"))]
impl IoAffinityGuard {
    pub fn enter(_core: Option<core_affinity::CoreId>) -> Result<Option<Self>, DynError> {
        Ok(None)
    }
}

#[cfg(test)]
mod cpu_allocation_tests {
    use super::*;
    #[test]
    fn reserves_one_of_requested_cpus_and_respects_allowed_mask() {
        let available: Vec<_> = [2, 4, 8, 10]
            .into_iter()
            .map(|id| core_affinity::CoreId { id })
            .collect();
        for n in 1..=6 {
            let plan = allocate_cpus(&available, n, None).unwrap();
            assert_eq!(plan.total, n.min(4));
            assert_eq!(plan.workers.len(), n.min(4).saturating_sub(1).max(1));
            if n > 1 {
                assert!(plan.workers.iter().all(|c| c.id != plan.io.id));
            }
        }
        let plan = allocate_cpus(&available, 3, Some(available[0])).unwrap();
        assert_eq!(plan.io.id, 2);
        assert_eq!(
            plan.workers.iter().map(|c| c.id).collect::<Vec<_>>(),
            [4, 8]
        );
        assert!(allocate_cpus(&available, 0, None).is_err());
        assert!(allocate_cpus(&[], 1, None).is_err());
        assert!(allocate_cpus(&available, 2, Some(core_affinity::CoreId { id: 0 })).is_err());
    }
}
