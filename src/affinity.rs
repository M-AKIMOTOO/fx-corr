use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use crate::utils::DynError;

const AFFINITY_FILE_BASENAME: &str = "yi-corr-affinity.txt";
const CPU_RESERVATION_STATE: &str = ".yi-corr";
const CPU_RESERVATION_LOCK: &str = ".yi-corr.lock";
const CPU_RESERVATION_HEADER: &str = "yi-corr-cpu-reservations-v1";
static NEXT_RESERVATION_ID: AtomicU64 = AtomicU64::new(1);

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

impl CpuAllocation {
    pub fn used_cores(&self) -> Vec<core_affinity::CoreId> {
        let mut cores = self.workers.clone();
        cores.push(self.io);
        cores.sort_unstable_by_key(|c| c.id);
        cores.dedup_by_key(|c| c.id);
        cores
    }
}

struct CpuTopology {
    groups: Vec<Vec<usize>>,
}

impl CpuTopology {
    fn discover(universe: &[core_affinity::CoreId]) -> Self {
        let mut groups = Vec::new();
        for cpu in universe {
            #[cfg(target_os = "linux")]
            let siblings = fs::read_to_string(format!(
                "/sys/devices/system/cpu/cpu{}/topology/thread_siblings_list",
                cpu.id
            ))
            .ok()
            .and_then(|text| parse_affinity_spec(&text).ok());
            #[cfg(not(target_os = "linux"))]
            let siblings: Option<Vec<usize>> = None;
            let group = siblings
                .filter(|g| g.contains(&cpu.id))
                .unwrap_or_else(|| vec![cpu.id]);
            if !groups.contains(&group) {
                groups.push(group);
            }
        }
        groups.sort_unstable();
        Self { groups }
    }
}

// Keep live assignments fixed. Reserve all SMT siblings of each selected
// physical core, while only the N used logical CPUs enter the thread masks.
fn isolated_allocation(
    available: &[core_affinity::CoreId],
    universe: &[core_affinity::CoreId],
    topology: &CpuTopology,
    requested: usize,
    preferred_io: Option<core_affinity::CoreId>,
    occupied: &HashSet<usize>,
) -> Result<(CpuAllocation, Vec<usize>), DynError> {
    let universe_ids: HashSet<_> = universe.iter().map(|c| c.id).collect();
    let allowed: HashSet<_> = available
        .iter()
        .map(|c| c.id)
        .filter(|id| universe_ids.contains(id))
        .collect();
    let target = requested.min(allowed.len());
    if target == 0 {
        return Err("CPU reservation requires a nonzero request and available CPUs".into());
    }
    let mut blocked = occupied.clone();
    for group in &topology.groups {
        if group.iter().any(|id| occupied.contains(id)) {
            blocked.extend(group);
        }
    }
    let free_count = universe_ids
        .iter()
        .filter(|id| !blocked.contains(id))
        .count();
    let budget = free_count.saturating_sub(2);
    let mut candidates: Vec<_> = topology
        .groups
        .iter()
        .filter_map(|group| {
            if group.iter().any(|id| blocked.contains(id)) {
                return None;
            }
            let usable: Vec<_> = group
                .iter()
                .copied()
                .filter(|id| allowed.contains(id))
                .collect();
            if usable.is_empty() {
                None
            } else {
                Some((group, usable))
            }
        })
        .collect();
    if let Some(io) = preferred_io {
        let index = candidates.iter().position(|(_, ids)| ids.contains(&io.id))
            .ok_or_else(|| format!("YI_READER_CORE={} is outside the allowed set or shares a reserved physical core", io.id))?;
        candidates.swap(0, index);
    }
    let mut selected = Vec::new();
    let mut claimed_groups = Vec::new();
    let mut cost = 0;
    // First spread workers over physical cores, then use their SMT siblings
    // only if the requested logical count needs them.
    for (group, ids) in &candidates {
        let extra = group.iter().filter(|id| universe_ids.contains(id)).count();
        if cost + extra > budget {
            continue;
        }
        let id = preferred_io
            .filter(|io| ids.contains(&io.id))
            .map_or(ids[0], |io| io.id);
        selected.push(id);
        claimed_groups.push((*group, ids));
        cost += extra;
        if selected.len() == target {
            break;
        }
    }
    if preferred_io.is_some_and(|io| !selected.contains(&io.id)) {
        return Err("cannot reserve YI_READER_CORE while leaving two visible CPUs free".into());
    }
    if selected.len() < target {
        for (_, ids) in &claimed_groups {
            for &id in *ids {
                if selected.len() == target {
                    break;
                }
                if !selected.contains(&id) {
                    selected.push(id);
                }
            }
            if selected.len() == target {
                break;
            }
        }
    }
    if selected.len() < target {
        return Err(format!("cannot reserve {target} logical CPUs without sharing physical cores with running yi-corr processes and leaving two visible CPUs free (only {} can be assigned)", selected.len()).into());
    }
    let io = preferred_io.or_else(|| {
        selected
            .iter()
            .rev()
            .find(|id| {
                topology.groups.iter().any(|g| {
                    g.contains(id) && selected.iter().filter(|c| g.contains(c)).count() == 1
                })
            })
            .map(|id| core_affinity::CoreId { id: *id })
    });
    let cores = selected
        .iter()
        .map(|id| core_affinity::CoreId { id: *id })
        .collect::<Vec<_>>();
    let allocation = allocate_cpus(&cores, target, io)?;
    let mut reserved: Vec<_> = claimed_groups
        .iter()
        .flat_map(|(g, _)| g.iter().copied())
        .collect();
    reserved.sort_unstable();
    reserved.dedup();
    Ok((allocation, reserved))
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

#[derive(Clone, Debug)]
struct CpuReservationEntry {
    token: String,
    pid: u32,
    start_time: String,
    cpus: Vec<usize>,
}

#[cfg(unix)]
struct CpuReservationFileLock(File);

#[cfg(unix)]
impl CpuReservationFileLock {
    fn acquire(path: &Path) -> Result<Self, DynError> {
        use std::os::fd::AsRawFd;

        let file = OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .open(path)?;
        if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX) } != 0 {
            return Err(std::io::Error::last_os_error().into());
        }
        Ok(Self(file))
    }
}

#[cfg(unix)]
impl Drop for CpuReservationFileLock {
    fn drop(&mut self) {
        use std::os::fd::AsRawFd;
        unsafe {
            libc::flock(self.0.as_raw_fd(), libc::LOCK_UN);
        }
    }
}

/// A live claim in `$HOME/.yi-corr`.
///
/// The registry is protected by an advisory lock for each update. The entry
/// remains present for this guard's lifetime and is reclaimed on Drop. On
/// Linux, entries left by a kill -9 or power loss are recognized by PID plus
/// `/proc` start time and removed by the next claimant.
#[cfg(unix)]
pub struct CpuReservationGuard {
    registry_path: PathBuf,
    lock_path: PathBuf,
    token: String,
    cpus: Vec<usize>,
}

#[cfg(unix)]
impl CpuReservationGuard {
    pub fn info(&self) -> String {
        format!(
            "registry={} reserved CPU IDs (including SMT siblings)={}",
            self.registry_path.display(),
            self.cpus
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

#[cfg(unix)]
impl Drop for CpuReservationGuard {
    fn drop(&mut self) {
        if let Err(error) =
            release_cpu_reservation(&self.registry_path, &self.lock_path, &self.token)
        {
            eprintln!("[warn] could not release CPU reservation: {error}");
        }
    }
}

#[cfg(not(unix))]
pub struct CpuReservationGuard;

#[cfg(not(unix))]
impl CpuReservationGuard {
    pub fn info(&self) -> String {
        "unavailable".to_string()
    }
}

/// Reserve a requested CPU set across concurrent yi-corr processes.
///
/// `available` is the CPU set this invocation may use (including any fixed
/// affinity group); `universe` is the process-visible host CPU set. At least
/// two IDs in `universe` are kept outside active reservations.
pub fn reserve_cpus(
    available: &[core_affinity::CoreId],
    universe: &[core_affinity::CoreId],
    requested: usize,
    preferred_io: Option<core_affinity::CoreId>,
) -> Result<(CpuAllocation, CpuReservationGuard), DynError> {
    #[cfg(unix)]
    {
        let home = home_dir().ok_or("cannot determine HOME for CPU reservation registry")?;
        reserve_cpus_in(
            &home.join(CPU_RESERVATION_STATE),
            &home.join(CPU_RESERVATION_LOCK),
            available,
            universe,
            requested,
            preferred_io,
        )
    }
    #[cfg(not(unix))]
    {
        let _ = (available, universe, requested, preferred_io);
        Err("cross-process --cpu reservations are currently supported on Unix".into())
    }
}

#[cfg(unix)]
fn reserve_cpus_in(
    state_path: &Path,
    lock_path: &Path,
    available: &[core_affinity::CoreId],
    universe: &[core_affinity::CoreId],
    requested: usize,
    preferred_io: Option<core_affinity::CoreId>,
) -> Result<(CpuAllocation, CpuReservationGuard), DynError> {
    reserve_cpus_with_topology_in(
        state_path,
        lock_path,
        available,
        universe,
        requested,
        preferred_io,
        &CpuTopology::discover(universe),
    )
}

#[cfg(unix)]
fn reserve_cpus_with_topology_in(
    state_path: &Path,
    lock_path: &Path,
    available: &[core_affinity::CoreId],
    universe: &[core_affinity::CoreId],
    requested: usize,
    preferred_io: Option<core_affinity::CoreId>,
    topology: &CpuTopology,
) -> Result<(CpuAllocation, CpuReservationGuard), DynError> {
    if requested == 0 || available.is_empty() || universe.is_empty() {
        return Err("CPU reservation requires a nonzero request and available CPUs".into());
    }
    let _lock = CpuReservationFileLock::acquire(&lock_path)?;

    let current_pid = std::process::id();
    let current_start = process_start_time(current_pid)?;
    let mut entries = read_cpu_reservations(state_path)?;
    entries.retain(reservation_is_live);

    let mut occupied = HashSet::<usize>::new();
    for entry in &entries {
        occupied.extend(entry.cpus.iter().copied());
    }
    let (allocation, reserved) = isolated_allocation(
        available,
        universe,
        topology,
        requested,
        preferred_io,
        &occupied,
    )?;

    let sequence = NEXT_RESERVATION_ID.fetch_add(1, Ordering::Relaxed);
    let token = format!("{current_pid}-{current_start}-{sequence}");
    entries.push(CpuReservationEntry {
        token: token.clone(),
        pid: current_pid,
        start_time: current_start,
        cpus: reserved.clone(),
    });
    write_cpu_reservations(state_path, &entries, current_pid, sequence)?;

    Ok((
        allocation,
        CpuReservationGuard {
            registry_path: state_path.to_path_buf(),
            lock_path: lock_path.to_path_buf(),
            token,
            cpus: reserved,
        },
    ))
}

#[cfg(unix)]
fn read_cpu_reservations(path: &Path) -> Result<Vec<CpuReservationEntry>, DynError> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines.next().transpose()?.unwrap_or_default();
    if header != CPU_RESERVATION_HEADER {
        return Err(format!(
            "unrecognized CPU reservation registry format in {}",
            path.display()
        )
        .into());
    }
    let mut entries = Vec::new();
    for (index, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 4 {
            return Err(format!(
                "invalid CPU reservation at {}:{}",
                path.display(),
                index + 2
            )
            .into());
        }
        let pid = fields[1].parse::<u32>()?;
        let cpus = fields[3]
            .split(',')
            .map(str::parse::<usize>)
            .collect::<Result<Vec<_>, _>>()?;
        if fields[0].is_empty() || fields[2].is_empty() || cpus.is_empty() {
            return Err(format!(
                "invalid CPU reservation at {}:{}",
                path.display(),
                index + 2
            )
            .into());
        }
        entries.push(CpuReservationEntry {
            token: fields[0].to_string(),
            pid,
            start_time: fields[2].to_string(),
            cpus,
        });
    }
    Ok(entries)
}

#[cfg(unix)]
fn write_cpu_reservations(
    path: &Path,
    entries: &[CpuReservationEntry],
    pid: u32,
    sequence: u64,
) -> Result<(), DynError> {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|e| format!("system time error: {e}"))?
        .as_nanos();
    let temporary = path.with_file_name(format!(
        ".{CPU_RESERVATION_STATE}.{pid}.{sequence}.{now}.tmp"
    ));
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&temporary)?;
    writeln!(file, "{CPU_RESERVATION_HEADER}")?;
    for entry in entries {
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            entry.token,
            entry.pid,
            entry.start_time,
            entry
                .cpus
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )?;
    }
    file.sync_all()?;
    fs::rename(&temporary, path)?;
    Ok(())
}

#[cfg(unix)]
fn release_cpu_reservation(
    state_path: &Path,
    lock_path: &Path,
    token: &str,
) -> Result<(), DynError> {
    let _lock = CpuReservationFileLock::acquire(lock_path)?;
    let mut entries = read_cpu_reservations(state_path)?;
    entries.retain(|entry| entry.token != token && reservation_is_live(entry));
    if entries.is_empty() {
        match fs::remove_file(state_path) {
            Ok(()) => Ok(()),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(error) => Err(error.into()),
        }
    } else {
        let sequence = NEXT_RESERVATION_ID.fetch_add(1, Ordering::Relaxed);
        write_cpu_reservations(state_path, &entries, std::process::id(), sequence)
    }
}

#[cfg(target_os = "linux")]
fn process_start_time(pid: u32) -> Result<String, DynError> {
    fn parse(stat: &str) -> Option<&str> {
        // The command name (field 2) is parenthesized and may contain spaces.
        let close = stat.rfind(')')?;
        stat.get(close + 1..)?.split_whitespace().nth(19)
    }
    let stat = fs::read_to_string(format!("/proc/{pid}/stat"))?;
    parse(&stat)
        .map(str::to_string)
        .ok_or_else(|| format!("cannot parse process start time for PID {pid}").into())
}

#[cfg(not(target_os = "linux"))]
fn process_start_time(_pid: u32) -> Result<String, DynError> {
    Ok("unknown".to_string())
}

#[cfg(target_os = "linux")]
fn reservation_is_live(entry: &CpuReservationEntry) -> bool {
    match process_start_time(entry.pid) {
        Ok(start_time) => start_time == entry.start_time,
        Err(error) => {
            let not_found = error
                .downcast_ref::<std::io::Error>()
                .map(|e| e.kind() == std::io::ErrorKind::NotFound)
                .unwrap_or(false);
            !not_found
        }
    }
}

#[cfg(not(target_os = "linux"))]
fn reservation_is_live(entry: &CpuReservationEntry) -> bool {
    let result = unsafe { libc::kill(entry.pid as libc::pid_t, 0) };
    result == 0 || std::io::Error::last_os_error().raw_os_error() != Some(libc::ESRCH)
}

// Pin the calling thread only, and restore its previous mask on scope exit.
// The outer guard bounds auxiliary threads; a nested output guard temporarily
// selects I/O without leaking a one-CPU mask into later workflow children.
#[cfg(target_os = "linux")]
pub struct ThreadAffinityGuard(libc::cpu_set_t);

#[cfg(target_os = "linux")]
impl ThreadAffinityGuard {
    pub fn enter(core: Option<core_affinity::CoreId>) -> Result<Option<Self>, DynError> {
        let Some(core) = core else {
            return Ok(None);
        };
        Self::enter_cores(&[core])
    }

    pub fn enter_cores(cores: &[core_affinity::CoreId]) -> Result<Option<Self>, DynError> {
        if cores.is_empty() {
            return Err("cannot set an empty thread CPU mask".into());
        }
        let mut previous = unsafe { std::mem::zeroed::<libc::cpu_set_t>() };
        let result = unsafe {
            libc::sched_getaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &mut previous)
        };
        if result != 0 {
            return Err(std::io::Error::last_os_error().into());
        }
        let mut mask = unsafe { std::mem::zeroed::<libc::cpu_set_t>() };
        for core in cores {
            if core.id >= libc::CPU_SETSIZE as usize
                || !unsafe { libc::CPU_ISSET(core.id, &previous) }
            {
                return Err(
                    format!("CPU {} is outside this thread's inherited mask", core.id).into(),
                );
            }
            unsafe {
                libc::CPU_SET(core.id, &mut mask);
            }
        }
        if unsafe { libc::sched_setaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &mask) } != 0
        {
            return Err(std::io::Error::last_os_error().into());
        }
        Ok(Some(Self(previous)))
    }
}

#[cfg(target_os = "linux")]
impl Drop for ThreadAffinityGuard {
    fn drop(&mut self) {
        unsafe {
            libc::sched_setaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &self.0);
        }
    }
}

#[cfg(not(target_os = "linux"))]
pub struct ThreadAffinityGuard;
#[cfg(not(target_os = "linux"))]
impl ThreadAffinityGuard {
    pub fn enter(_core: Option<core_affinity::CoreId>) -> Result<Option<Self>, DynError> {
        Ok(None)
    }
    pub fn enter_cores(_cores: &[core_affinity::CoreId]) -> Result<Option<Self>, DynError> {
        Ok(None)
    }
}

pub type IoAffinityGuard = ThreadAffinityGuard;

#[cfg(test)]
mod cpu_allocation_tests {
    use super::*;

    // Original registry tests use synthetic single-threaded cores, independent
    // of the machine executing the test.
    #[cfg(unix)]
    fn reserve_cpus_in(
        state: &Path,
        lock: &Path,
        available: &[core_affinity::CoreId],
        universe: &[core_affinity::CoreId],
        requested: usize,
        io: Option<core_affinity::CoreId>,
    ) -> Result<(CpuAllocation, CpuReservationGuard), DynError> {
        reserve_cpus_with_topology_in(
            state,
            lock,
            available,
            universe,
            requested,
            io,
            &CpuTopology {
                groups: universe.iter().map(|c| vec![c.id]).collect(),
            },
        )
    }

    #[test]
    fn different_processes_never_select_smt_siblings_even_with_partial_masks() {
        let cpus: Vec<_> = (0..12).map(|id| core_affinity::CoreId { id }).collect();
        let topology = CpuTopology {
            groups: (0..6).map(|i| vec![2 * i, 2 * i + 1]).collect(),
        };
        let (first, reserved) =
            isolated_allocation(&cpus, &cpus, &topology, 2, None, &HashSet::new()).unwrap();
        assert_eq!(
            first.used_cores().iter().map(|c| c.id).collect::<Vec<_>>(),
            [0, 2]
        );
        assert_eq!(reserved, [0, 1, 2, 3]);
        let occupied = reserved.into_iter().collect();
        let (second, reserved2) =
            isolated_allocation(&cpus, &cpus, &topology, 2, None, &occupied).unwrap();
        assert_eq!(
            second.used_cores().iter().map(|c| c.id).collect::<Vec<_>>(),
            [4, 6]
        );
        assert_eq!(reserved2, [4, 5, 6, 7]);
        assert_eq!(
            first.used_cores().iter().map(|c| c.id).collect::<Vec<_>>(),
            [0, 2]
        );
        let partial: Vec<_> = cpus.iter().copied().filter(|c| c.id % 2 == 1).collect();
        let (third, _) =
            isolated_allocation(&partial, &partial, &topology, 2, None, &HashSet::from([0]))
                .unwrap();
        assert_eq!(
            third.used_cores().iter().map(|c| c.id).collect::<Vec<_>>(),
            [3, 5]
        );
        assert!(isolated_allocation(
            &cpus,
            &cpus,
            &topology,
            2,
            Some(core_affinity::CoreId { id: 1 }),
            &occupied
        )
        .is_err());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn thread_mask_is_local_restored_and_inherited_by_own_children() {
        let initial = core_affinity::get_core_ids().unwrap();
        let selected: Vec<_> = initial.iter().copied().take(2).collect();
        let (tx, rx) = std::sync::mpsc::channel();
        let observer = std::thread::spawn(move || {
            rx.recv().unwrap();
            core_affinity::get_core_ids().unwrap()
        });
        let guard = ThreadAffinityGuard::enter_cores(&selected).unwrap();
        assert_eq!(core_affinity::get_core_ids().unwrap(), selected);
        tx.send(()).unwrap();
        assert_eq!(observer.join().unwrap(), initial);
        assert_eq!(
            std::thread::spawn(|| core_affinity::get_core_ids().unwrap())
                .join()
                .unwrap(),
            selected
        );
        let io = ThreadAffinityGuard::enter(Some(selected[0])).unwrap();
        assert_eq!(core_affinity::get_core_ids().unwrap(), [selected[0]]);
        drop(io);
        assert_eq!(core_affinity::get_core_ids().unwrap(), selected);
        drop(guard);
        assert_eq!(core_affinity::get_core_ids().unwrap(), initial);
    }
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

    #[cfg(unix)]
    #[test]
    fn concurrent_reservations_are_disjoint_keep_two_cpus_free_and_release() {
        let sequence = NEXT_RESERVATION_ID.fetch_add(1, Ordering::Relaxed);
        let directory = std::env::temp_dir().join(format!(
            "yi-corr-reservation-test-{}-{sequence}",
            std::process::id()
        ));
        fs::create_dir(&directory).unwrap();
        let registry = directory.join(CPU_RESERVATION_STATE);
        let lock = directory.join(CPU_RESERVATION_LOCK);
        let cpus: Vec<_> = (0..32).map(|id| core_affinity::CoreId { id }).collect();

        let (first, first_guard) =
            reserve_cpus_in(&registry, &lock, &cpus, &cpus, 6, None).unwrap();
        let (second, second_guard) =
            reserve_cpus_in(&registry, &lock, &cpus, &cpus, 6, None).unwrap();
        let (third, third_guard) =
            reserve_cpus_in(&registry, &lock, &cpus, &cpus, 12, None).unwrap();
        assert_eq!(first.io.id, 5);
        assert_eq!(second.io.id, 11);
        assert_eq!(third.io.id, 23);
        assert_eq!(
            first.workers.iter().map(|core| core.id).collect::<Vec<_>>(),
            [0, 1, 2, 3, 4]
        );
        assert_eq!(
            second
                .workers
                .iter()
                .map(|core| core.id)
                .collect::<Vec<_>>(),
            [6, 7, 8, 9, 10]
        );
        assert_eq!(
            third.workers.iter().map(|core| core.id).collect::<Vec<_>>(),
            (12..23).collect::<Vec<_>>()
        );
        assert!(reserve_cpus_in(&registry, &lock, &cpus, &cpus, 8, None).is_err());

        drop(second_guard);
        let (reused, reused_guard) =
            reserve_cpus_in(&registry, &lock, &cpus, &cpus, 6, None).unwrap();
        assert_eq!(reused.io.id, 11);
        assert_eq!(
            reused
                .workers
                .iter()
                .map(|core| core.id)
                .collect::<Vec<_>>(),
            [6, 7, 8, 9, 10]
        );

        drop((first_guard, third_guard, reused_guard));
        let entries = read_cpu_reservations(&registry).unwrap();
        assert!(entries.is_empty());
        let (maximum, maximum_guard) =
            reserve_cpus_in(&registry, &lock, &cpus, &cpus, 30, None).unwrap();
        assert_eq!(maximum.total, 30);
        assert!(reserve_cpus_in(&registry, &lock, &cpus, &cpus, 1, None).is_err());
        drop(maximum_guard);
        fs::remove_dir_all(directory).unwrap();
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn stale_pid_reservations_are_removed_before_allocation() {
        let sequence = NEXT_RESERVATION_ID.fetch_add(1, Ordering::Relaxed);
        let directory = std::env::temp_dir().join(format!(
            "yi-corr-stale-reservation-test-{}-{sequence}",
            std::process::id()
        ));
        fs::create_dir(&directory).unwrap();
        let registry = directory.join(CPU_RESERVATION_STATE);
        let lock = directory.join(CPU_RESERVATION_LOCK);
        let cpus: Vec<_> = (0..8).map(|id| core_affinity::CoreId { id }).collect();
        write_cpu_reservations(
            &registry,
            &[CpuReservationEntry {
                token: "dead-process".to_string(),
                pid: u32::MAX,
                start_time: "0".to_string(),
                cpus: (0..6).collect(),
            }],
            std::process::id(),
            sequence,
        )
        .unwrap();

        let (allocation, guard) = reserve_cpus_in(&registry, &lock, &cpus, &cpus, 6, None).unwrap();
        assert_eq!(allocation.io.id, 5);
        let entries = read_cpu_reservations(&registry).unwrap();
        assert_eq!(entries.len(), 1);
        assert_ne!(entries[0].token, "dead-process");
        drop(guard);
        fs::remove_dir_all(directory).unwrap();
    }
}
