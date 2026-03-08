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
    if out.is_empty() { None } else { Some(out) }
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
