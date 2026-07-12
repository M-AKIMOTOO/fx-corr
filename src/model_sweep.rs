use crate::args::Args;
use crate::utils::DynError;
use std::ffi::{OsStr, OsString};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

pub(crate) const SWEEP_CHILD_ENV: &str = "YI_MODEL_SWEEP_CHILD";

const LOCK_FILE: &str = ".model_sweep.lock";
const OK_MARKER: &str = ".model_sweep_ok";
const CASE_LOG: &str = "model_sweep_stdout_stderr.log";
const MANIFEST_FILE: &str = "model_sweep_status.tsv";

#[derive(Clone, Copy)]
struct SweepCase {
    key: &'static str,
    source_mode: Option<&'static str>,
    delay_mode: Option<&'static str>,
    force_zero_eop: bool,
    description: &'static str,
}

const CASES: [SweepCase; 5] = [
    SweepCase {
        key: "01_active",
        source_mode: None,
        delay_mode: None,
        force_zero_eop: false,
        description: "unchanged caller/XML/environment model",
    },
    SweepCase {
        key: "02_mean_anchored_eop_zero",
        source_mode: Some("mean-gast"),
        delay_mode: Some("anchored"),
        force_zero_eop: true,
        description: "legacy source/delay model with DUT1=xp=yp=0",
    },
    SweepCase {
        key: "03_pnm_anchored_eop_zero",
        source_mode: Some("pnm-gast"),
        delay_mode: Some("anchored"),
        force_zero_eop: true,
        description: "PNM06A+GAST source model with DUT1=xp=yp=0",
    },
    SweepCase {
        key: "04_pnm_vlbi_minus_eop_zero",
        source_mode: Some("pnm-gast"),
        delay_mode: Some("vlbi-minus"),
        force_zero_eop: true,
        description: "first-order IERS-sign VLBI model with DUT1=xp=yp=0",
    },
    SweepCase {
        key: "05_pnm_vlbi_minus_active_eop",
        source_mode: Some("pnm-gast"),
        delay_mode: Some("vlbi-minus"),
        force_zero_eop: false,
        description: "first-order IERS-sign VLBI model with caller/XML EOP",
    },
];

const REMOVED_VALUE_OPTIONS: [&str; 3] = ["cor-directory", "output", "cor"];
const REMOVED_FLAGS: [&str; 4] = ["model-sweep", "model-diagnostics", "stdout", "compact-logs"];

fn env_flag_value(value: Option<&OsStr>) -> bool {
    value.is_some_and(|value| {
        let value = value.to_string_lossy().trim().to_ascii_lowercase();
        !(value.is_empty() || value == "0" || value == "false" || value == "off" || value == "no")
    })
}

fn recursion_guard(model_sweep_requested: bool, marker: Option<&OsStr>) -> Result<(), String> {
    if model_sweep_requested && env_flag_value(marker) {
        Err(format!(
            "refusing recursive --model-sweep in a sweep child ({SWEEP_CHILD_ENV} is set)"
        ))
    } else {
        Ok(())
    }
}

pub(crate) fn ensure_not_recursive_model_sweep(
    model_sweep_requested: bool,
) -> Result<(), DynError> {
    recursion_guard(
        model_sweep_requested,
        std::env::var_os(SWEEP_CHILD_ENV).as_deref(),
    )
    .map_err(Into::into)
}

fn inferred_long_match(candidate: &str, names: &[&str]) -> bool {
    !candidate.is_empty() && names.iter().any(|name| name.starts_with(candidate))
}

fn split_long_option(argument: &str) -> Option<(&str, bool)> {
    let body = argument.strip_prefix("--")?;
    Some(match body.split_once('=') {
        Some((name, _)) => (name, true),
        None => (body, false),
    })
}

fn forwarded_arguments_from<I, T>(arguments: I) -> Vec<OsString>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString>,
{
    let mut output = Vec::new();
    let mut skip_next = false;
    for argument in arguments {
        let argument = argument.into();
        if skip_next {
            skip_next = false;
            continue;
        }
        let text = argument.to_string_lossy();
        if let Some((name, attached_value)) = split_long_option(&text) {
            if inferred_long_match(name, &REMOVED_VALUE_OPTIONS) {
                skip_next = !attached_value;
                continue;
            }
            if inferred_long_match(name, &REMOVED_FLAGS) {
                continue;
            }
        }
        output.push(argument);
    }
    output
}

fn forwarded_arguments() -> Vec<OsString> {
    forwarded_arguments_from(std::env::args_os().skip(1))
}

fn configure_case(command: &mut Command, case: SweepCase) {
    command.env(SWEEP_CHILD_ENV, "1");
    if let Some(source_mode) = case.source_mode {
        command.env("YI_SOURCE_VECTOR_MODE", source_mode);
    }
    if let Some(delay_mode) = case.delay_mode {
        command.env("YI_GEOM_DELAY_MODE", delay_mode);
    }
    if case.force_zero_eop {
        command
            .env("YI_EOP_MODE", "none")
            .env_remove("YI_EOP_FILE")
            .env("YI_DUT1_S", "0")
            .env("YI_XP_ARCSEC", "0")
            .env("YI_YP_ARCSEC", "0");
        // TT-UTC is a time-scale conversion, not an EOP term to zero. Preserve
        // the caller/XML/YI_TT_UTC_S value for a controlled comparison.
    }
}

fn hash_bytes(hash: &mut u64, bytes: &[u8]) {
    const FNV_PRIME: u64 = 1_099_511_628_211;
    for byte in bytes {
        *hash ^= u64::from(*byte);
        *hash = hash.wrapping_mul(FNV_PRIME);
    }
    // Field separator so concatenated values cannot alias different argv/env splits.
    *hash ^= 0xff;
    *hash = hash.wrapping_mul(FNV_PRIME);
}

fn hash_os_str(hash: &mut u64, value: &OsStr) {
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        hash_bytes(hash, value.as_bytes());
    }
    #[cfg(not(unix))]
    {
        hash_bytes(hash, value.to_string_lossy().as_bytes());
    }
}

fn hash_file_state(
    hash: &mut u64,
    path: &Path,
    metadata: &fs::Metadata,
    include_full_contents: bool,
) -> Result<(), DynError> {
    hash_os_str(hash, path.as_os_str());
    hash_bytes(hash, &metadata.len().to_le_bytes());
    hash_bytes(hash, &[u8::from(metadata.is_file())]);
    if let Ok(modified) = metadata.modified() {
        if let Ok(duration) = modified.duration_since(UNIX_EPOCH) {
            hash_bytes(hash, &duration.as_secs().to_le_bytes());
            hash_bytes(hash, &duration.subsec_nanos().to_le_bytes());
        }
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        hash_bytes(hash, &metadata.dev().to_le_bytes());
        hash_bytes(hash, &metadata.ino().to_le_bytes());
        hash_bytes(hash, &metadata.mtime().to_le_bytes());
        hash_bytes(hash, &metadata.mtime_nsec().to_le_bytes());
    }
    if !metadata.is_file() {
        return Ok(());
    }
    if include_full_contents && metadata.len() <= 32 * 1024 * 1024 {
        hash_bytes(hash, &fs::read(path)?);
        return Ok(());
    }

    let mut file = File::open(path)?;
    let mut edge = [0u8; 4096];
    let first = file.read(&mut edge)?;
    hash_bytes(hash, &edge[..first]);
    if metadata.len() > edge.len() as u64 {
        file.seek(SeekFrom::End(-(edge.len() as i64)))?;
        let last = file.read(&mut edge)?;
        hash_bytes(hash, &edge[..last]);
    }
    Ok(())
}

fn is_raw_input_candidate(path: &Path) -> bool {
    path.extension()
        .and_then(|value| value.to_str())
        .is_some_and(|extension| {
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "raw" | "vdif" | "dat"
            )
        })
}

fn hash_input_path(
    hash: &mut u64,
    path: &Path,
    include_full_contents: bool,
) -> Result<(), DynError> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    let canonical = fs::canonicalize(&absolute).unwrap_or(absolute);
    hash_os_str(hash, canonical.as_os_str());
    let metadata = match fs::metadata(&canonical) {
        Ok(metadata) => metadata,
        Err(error) => {
            hash_bytes(hash, format!("missing:{error}").as_bytes());
            return Ok(());
        }
    };
    if metadata.is_file() {
        return hash_file_state(hash, &canonical, &metadata, include_full_contents);
    }
    if !metadata.is_dir() {
        return Ok(());
    }

    let mut inputs = fs::read_dir(&canonical)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|entry| is_raw_input_candidate(entry))
        .collect::<Vec<_>>();
    inputs.sort();
    for input in inputs {
        let metadata = fs::metadata(&input)?;
        hash_file_state(hash, &input, &metadata, false)?;
    }
    Ok(())
}

fn hash_cli_input_state(hash: &mut u64, args: &Args) -> Result<(), DynError> {
    let cwd = std::env::current_dir()?;
    hash_os_str(hash, fs::canonicalize(&cwd).unwrap_or(cwd).as_os_str());
    if let Some(path) = args.schedule.as_deref() {
        hash_input_path(hash, path, true)?;
    }
    if let Some(schedule) = args.schedule.as_ref() {
        if let Ok(parsed) = crate::parse_ifile_cached(schedule) {
            if let Some(eop_file) = parsed.eop_file.as_deref() {
                hash_input_path(hash, eop_file, true)?;
            }
        }
    }
    let default_eop = Path::new("data/eop/finals2000A.data");
    if default_eop.exists() {
        hash_input_path(hash, default_eop, true)?;
    }
    if args.ant1.is_none() || args.ant2.is_none() {
        if let Some(path) = args.raw_directory.as_deref() {
            hash_input_path(hash, path, false)?;
        }
    }
    for path in [args.ant1.as_deref(), args.ant2.as_deref()]
        .into_iter()
        .flatten()
    {
        hash_input_path(hash, path, false)?;
        if path.is_relative() {
            if let Some(raw_directory) = args.raw_directory.as_deref() {
                hash_input_path(hash, &raw_directory.join(path), false)?;
            }
        }
    }
    Ok(())
}

fn base_run_signature(
    executable: &Path,
    forwarded: &[OsString],
    args: &Args,
) -> Result<u64, DynError> {
    const FNV_OFFSET: u64 = 14_695_981_039_346_656_037;
    let mut hash = FNV_OFFSET;
    let metadata = executable.metadata()?;
    hash_file_state(&mut hash, executable, &metadata, true)?;
    hash_cli_input_state(&mut hash, args)?;
    for argument in forwarded {
        hash_os_str(&mut hash, argument);
    }
    let mut yi_environment = std::env::vars_os()
        .filter(|(key, _)| key.to_string_lossy().starts_with("YI_"))
        .collect::<Vec<_>>();
    yi_environment.sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
    for (key, value) in yi_environment {
        hash_os_str(&mut hash, &key);
        hash_os_str(&mut hash, &value);
        if key.to_string_lossy().ends_with("_FILE") {
            hash_input_path(&mut hash, Path::new(&value), true)?;
        }
    }
    Ok(hash)
}

fn case_run_signature(base_signature: u64, case: SweepCase) -> String {
    let mut hash = base_signature;
    hash_bytes(&mut hash, case.key.as_bytes());
    hash_bytes(
        &mut hash,
        case.source_mode.unwrap_or("caller-active").as_bytes(),
    );
    hash_bytes(
        &mut hash,
        case.delay_mode.unwrap_or("caller-active").as_bytes(),
    );
    hash_bytes(&mut hash, &[u8::from(case.force_zero_eop)]);
    format!("{hash:016x}")
}

fn unix_timestamp() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| format!("{}.{:03}", duration.as_secs(), duration.subsec_millis()))
        .unwrap_or_else(|_| "0.000".to_string())
}

fn tsv_field(value: impl AsRef<str>) -> String {
    value
        .as_ref()
        .replace(['\t', '\r', '\n'], " ")
        .trim()
        .to_string()
}

fn open_manifest(path: &Path, root: &Path) -> Result<File, DynError> {
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .read(true)
        .open(path)?;
    if file.metadata()?.len() == 0 {
        writeln!(file, "# yi-corr unattended model sweep")?;
        writeln!(file, "# root={}", root.display())?;
        writeln!(file, "timestamp\tevent\tcase\tstatus\texit_code\telapsed_s\tdescription\toutput\tlog\tdetail")?;
    } else {
        let contents = fs::read_to_string(path)?;
        let root_header = format!("# root={}", root.display());
        let schema = "timestamp\tevent\tcase\tstatus\texit_code\telapsed_s\tdescription\toutput\tlog\tdetail";
        if !contents.starts_with("# yi-corr unattended model sweep\n")
            || !contents.lines().any(|line| line == root_header)
            || !contents.lines().any(|line| line == schema)
        {
            return Err(format!(
                "existing model sweep manifest has an incompatible header: {}",
                path.display()
            )
            .into());
        }
        if !contents.ends_with('\n') {
            writeln!(file)?;
            writeln!(
                file,
                "# recovered_partial_manifest_line_at_unix={}",
                unix_timestamp()
            )?;
        }
    }
    file.flush()?;
    file.sync_data()?;
    Ok(file)
}
#[allow(clippy::too_many_arguments)]
fn record_manifest(
    manifest: &mut File,
    event: &str,
    case: &str,
    status: &str,
    exit_code: &str,
    elapsed_s: f64,
    description: &str,
    output: &Path,
    log: Option<&Path>,
    detail: &str,
) -> std::io::Result<()> {
    writeln!(
        manifest,
        "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}",
        unix_timestamp(),
        tsv_field(event),
        tsv_field(case),
        tsv_field(status),
        tsv_field(exit_code),
        elapsed_s,
        tsv_field(description),
        tsv_field(output.to_string_lossy()),
        log.map(|path| tsv_field(path.to_string_lossy()))
            .unwrap_or_else(|| "-".to_string()),
        tsv_field(detail),
    )?;
    manifest.flush()?;
    manifest.sync_data()
}

#[allow(clippy::too_many_arguments)]
fn record_or_warn(
    manifest: &mut File,
    manifest_error: &mut Option<String>,
    event: &str,
    case: &str,
    status: &str,
    exit_code: &str,
    elapsed_s: f64,
    description: &str,
    output: &Path,
    log: Option<&Path>,
    detail: &str,
) {
    if let Err(error) = record_manifest(
        manifest,
        event,
        case,
        status,
        exit_code,
        elapsed_s,
        description,
        output,
        log,
        detail,
    ) {
        if manifest_error.is_none() {
            *manifest_error = Some(error.to_string());
        }
        crate::emit_stdout_fmt(
            format_args!("[model-sweep] manifest write failed; cases continue: {error}"),
            true,
        );
    }
}

struct SweepLock {
    file: File,
}

impl SweepLock {
    #[cfg(unix)]
    fn acquire(root: &Path) -> Result<Self, DynError> {
        use std::os::fd::AsRawFd;
        use std::os::unix::fs::OpenOptionsExt;

        let path = root.join(LOCK_FILE);
        let file = OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .custom_flags(libc::O_NOFOLLOW)
            .open(&path)?;
        if !file.metadata()?.is_file() {
            return Err(
                format!("model sweep lock is not a regular file: {}", path.display()).into(),
            );
        }
        let fd = file.as_raw_fd();
        if unsafe { libc::flock(fd, libc::LOCK_EX | libc::LOCK_NB) } != 0 {
            let error = std::io::Error::last_os_error();
            if error.kind() == std::io::ErrorKind::WouldBlock {
                return Err(
                    format!("model sweep is already running (lock={})", path.display()).into(),
                );
            }
            return Err(format!("cannot lock {}: {error}", path.display()).into());
        }

        let fd_flags = unsafe { libc::fcntl(fd, libc::F_GETFD) };
        if fd_flags < 0
            || unsafe { libc::fcntl(fd, libc::F_SETFD, fd_flags & !libc::FD_CLOEXEC) } < 0
        {
            let error = std::io::Error::last_os_error();
            unsafe {
                libc::flock(fd, libc::LOCK_UN);
            }
            return Err(format!(
                "cannot make sweep lock inheritable for child protection: {error}"
            )
            .into());
        }
        Ok(Self { file })
    }

    #[cfg(not(unix))]
    fn acquire(_root: &Path) -> Result<Self, DynError> {
        Err("--model-sweep requires Unix flock support".into())
    }
}

impl Drop for SweepLock {
    fn drop(&mut self) {
        #[cfg(unix)]
        {
            use std::os::fd::AsRawFd;
            unsafe {
                libc::flock(self.file.as_raw_fd(), libc::LOCK_UN);
            }
        }
    }
}
fn validate_cor_file(path: &Path, metadata: &fs::Metadata) -> Result<i32, String> {
    let mut file =
        File::open(path).map_err(|error| format!("cannot open {}: {error}", path.display()))?;
    let mut header = [0u8; 32];
    file.read_exact(&mut header)
        .map_err(|error| format!("short .cor header {}: {error}", path.display()))?;
    if header[..4] != [0x83, 0xf9, 0xa2, 0x3e] {
        return Err(format!("invalid .cor magic: {}", path.display()));
    }
    let fft_point = i32::from_le_bytes(header[24..28].try_into().unwrap());
    let sectors = i32::from_le_bytes(header[28..32].try_into().unwrap());
    if fft_point <= 0 || fft_point % 2 != 0 || sectors <= 0 {
        return Err(format!(
            "invalid .cor header values: {} fft_point={} sectors={}",
            path.display(),
            fft_point,
            sectors
        ));
    }
    let sector_bytes = 128u64
        .checked_add((fft_point as u64 / 2).saturating_mul(8))
        .ok_or_else(|| format!(".cor size overflow: {}", path.display()))?;
    let expected_size = 256u64
        .checked_add((sectors as u64).saturating_mul(sector_bytes))
        .ok_or_else(|| format!(".cor size overflow: {}", path.display()))?;
    if metadata.len() != expected_size {
        return Err(format!(
            "incomplete .cor file: {} actual_bytes={} expected_bytes={} sectors={}",
            path.display(),
            metadata.len(),
            expected_size,
            sectors
        ));
    }
    Ok(sectors)
}

#[derive(Clone, Copy, Debug)]
struct ArtifactSummary {
    cor_files: usize,
    diagnostic_files: usize,
    delay_model_files: usize,
    cor_sector_total: u64,
    fingerprint: u64,
}

fn file_is_fresh(metadata: &fs::Metadata, minimum_modified: Option<SystemTime>) -> bool {
    let Some(minimum) = minimum_modified else {
        return true;
    };
    let minimum = minimum
        .checked_sub(Duration::from_secs(2))
        .unwrap_or(UNIX_EPOCH);
    metadata
        .modified()
        .map(|modified| modified >= minimum)
        .unwrap_or(false)
}

fn validate_artifacts(
    output_dir: &Path,
    log_path: &Path,
    minimum_modified: Option<SystemTime>,
) -> Result<ArtifactSummary, String> {
    let log_metadata = log_path
        .metadata()
        .map_err(|error| format!("missing case log {}: {error}", log_path.display()))?;
    if log_metadata.len() == 0 || !file_is_fresh(&log_metadata, minimum_modified) {
        return Err(format!(
            "case log is empty or stale: {}",
            log_path.display()
        ));
    }

    let mut cor_files = 0usize;
    let mut diagnostic_files = 0usize;
    let mut delay_model_files = 0usize;
    let mut cor_sector_total = 0u64;
    let mut fingerprint = 14_695_981_039_346_656_037u64;
    let mut entries = fs::read_dir(output_dir)
        .map_err(|error| format!("cannot inspect {}: {error}", output_dir.display()))?
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| error.to_string())?;
    entries.sort_by_key(|entry| entry.file_name());

    for entry in entries {
        let path = entry.path();
        let name = entry.file_name();
        let name = name.to_string_lossy();
        let is_cor = name.ends_with(".cor");
        let is_diagnostic = name.starts_with("model_diagnostics_") && name.ends_with(".tsv");
        let is_delay_model = name.starts_with("delay_model_") && name.ends_with(".tsv");
        if !is_cor && !is_diagnostic && !is_delay_model {
            continue;
        }
        let metadata = entry.metadata().map_err(|error| error.to_string())?;
        if !metadata.is_file() || metadata.len() == 0 {
            return Err(format!("empty or non-file artifact: {}", path.display()));
        }
        if !file_is_fresh(&metadata, minimum_modified) {
            return Err(format!(
                "stale artifact from another attempt: {}",
                path.display()
            ));
        }
        hash_file_state(&mut fingerprint, &path, &metadata, false)
            .map_err(|error| error.to_string())?;
        if is_cor {
            let sectors = validate_cor_file(&path, &metadata)?;
            cor_files += 1;
            cor_sector_total = cor_sector_total
                .checked_add(sectors as u64)
                .ok_or_else(|| "cor sector count overflow".to_string())?;
        } else if is_diagnostic {
            diagnostic_files += 1;
        } else if is_delay_model {
            delay_model_files += 1;
        }
    }
    if cor_files < 3 {
        return Err(format!(
            "expected at least three complete .cor products, found {cor_files}"
        ));
    }
    if diagnostic_files == 0 {
        return Err("no non-empty model_diagnostics_*.tsv product found".to_string());
    }
    if delay_model_files == 0 {
        return Err("no non-empty delay_model_*.tsv product found".to_string());
    }
    Ok(ArtifactSummary {
        cor_files,
        diagnostic_files,
        delay_model_files,
        cor_sector_total,
        fingerprint,
    })
}

fn marker_matches_case(marker_path: &Path, case: SweepCase, signature: &str) -> bool {
    let mut contents = String::new();
    File::open(marker_path)
        .and_then(|mut file| file.read_to_string(&mut contents))
        .is_ok()
        && contents.lines().any(|line| line == "status=ok")
        && contents
            .lines()
            .any(|line| line == format!("case={}", case.key))
        && contents
            .lines()
            .any(|line| line == format!("run_signature={signature}"))
}

fn marker_matches_artifacts(marker_path: &Path, artifacts: ArtifactSummary) -> bool {
    let mut contents = String::new();
    File::open(marker_path)
        .and_then(|mut file| file.read_to_string(&mut contents))
        .is_ok()
        && contents
            .lines()
            .any(|line| line == format!("cor_files={}", artifacts.cor_files))
        && contents
            .lines()
            .any(|line| line == format!("diagnostic_files={}", artifacts.diagnostic_files))
        && contents
            .lines()
            .any(|line| line == format!("delay_model_files={}", artifacts.delay_model_files))
        && contents
            .lines()
            .any(|line| line == format!("cor_sector_total={}", artifacts.cor_sector_total))
        && contents
            .lines()
            .any(|line| line == format!("artifact_fingerprint={:016x}", artifacts.fingerprint))
}

fn write_ok_marker(
    output_dir: &Path,
    case: SweepCase,
    signature: &str,
    artifacts: ArtifactSummary,
) -> Result<(), DynError> {
    let marker_path = output_dir.join(OK_MARKER);
    let temporary_path = output_dir.join(format!(".model_sweep_ok.{}.tmp", std::process::id()));
    let mut file = File::create(&temporary_path)?;
    writeln!(file, "status=ok")?;
    writeln!(file, "case={}", case.key)?;
    writeln!(file, "run_signature={signature}")?;
    writeln!(file, "cor_files={}", artifacts.cor_files)?;
    writeln!(file, "diagnostic_files={}", artifacts.diagnostic_files)?;
    writeln!(file, "delay_model_files={}", artifacts.delay_model_files)?;
    writeln!(file, "cor_sector_total={}", artifacts.cor_sector_total)?;
    writeln!(file, "artifact_fingerprint={:016x}", artifacts.fingerprint)?;
    writeln!(file, "completed_unix={}", unix_timestamp())?;
    file.flush()?;
    file.sync_all()?;
    fs::rename(&temporary_path, &marker_path)?;
    Ok(())
}

fn append_case_log(path: &Path, message: &str) {
    if let Ok(mut file) = OpenOptions::new().create(true).append(true).open(path) {
        let _ = writeln!(file, "{message}");
        let _ = file.flush();
    }
}

fn clear_previous_case_products(output_dir: &Path) -> Result<(), String> {
    let entries = fs::read_dir(output_dir)
        .map_err(|error| format!("cannot inspect {}: {error}", output_dir.display()))?;
    for entry in entries {
        let entry = entry.map_err(|error| error.to_string())?;
        let metadata = entry.metadata().map_err(|error| error.to_string())?;
        if !metadata.is_file() {
            continue;
        }
        let name = entry.file_name();
        let name = name.to_string_lossy();
        let generated = name.ends_with(".cor")
            || (name.starts_with("model_diagnostics_") && name.ends_with(".tsv"))
            || (name.starts_with("delay_model_") && name.ends_with(".tsv"))
            || name.starts_with(".model_sweep_ok.");
        if generated {
            fs::remove_file(entry.path())
                .map_err(|error| format!("cannot remove {}: {error}", entry.path().display()))?;
        }
    }
    Ok(())
}

fn record_case_setup_error(
    manifest: &mut File,
    manifest_error: &mut Option<String>,
    failures: &mut Vec<String>,
    case: SweepCase,
    output_dir: &Path,
    log_path: &Path,
    error: &str,
) {
    failures.push(format!("{}({error})", case.key));
    record_or_warn(
        manifest,
        manifest_error,
        "case-end",
        case.key,
        "setup-error",
        "-1",
        0.0,
        case.description,
        output_dir,
        Some(log_path),
        error,
    );
}

pub(crate) fn run_unattended_model_sweep(args: &Args) -> Result<(), DynError> {
    ensure_not_recursive_model_sweep(true)?;
    let root = args
        .cor_directory
        .as_ref()
        .ok_or("--model-sweep requires --cor-directory")?
        .join("model_sweep");
    fs::create_dir_all(&root)?;
    let _lock = SweepLock::acquire(&root)?;

    let manifest_path = root.join(MANIFEST_FILE);
    let mut manifest = open_manifest(&manifest_path, &root)?;
    let mut manifest_error = None;
    let executable = std::env::current_exe()?;
    let forwarded = forwarded_arguments();
    let fixed_base_signature = base_run_signature(&executable, &forwarded, args)?;
    let run_signatures = CASES
        .iter()
        .copied()
        .map(|case| case_run_signature(fixed_base_signature, case))
        .collect::<Vec<_>>();
    let sweep_started = Instant::now();
    let mut failures = Vec::new();

    record_or_warn(
        &mut manifest,
        &mut manifest_error,
        "sweep-start",
        "-",
        "running",
        "-",
        0.0,
        "unattended five-case model sweep",
        &root,
        None,
        &format!(
            "pid={} executable={} base_signature={fixed_base_signature:016x} argv={forwarded:?}",
            std::process::id(),
            executable.display()
        ),
    );
    crate::emit_stdout_fmt(
        format_args!(
            "[model-sweep] unattended sweep start: cases={} root={} manifest={}",
            CASES.len(),
            root.display(),
            manifest_path.display()
        ),
        true,
    );
    crate::emit_stdout_fmt(
        format_args!(
            "[model-sweep] raw reading and FFT are repeated once per case; completed, validated cases are resumed without rerunning"
        ),
        true,
    );

    for (index, case) in CASES.iter().copied().enumerate() {
        let output_dir = root.join(case.key);
        let log_path = output_dir.join(CASE_LOG);
        let marker_path = output_dir.join(OK_MARKER);
        let run_signature = run_signatures[index].clone();
        match base_run_signature(&executable, &forwarded, args) {
            Ok(current) if current == fixed_base_signature => {}
            Ok(current) => {
                record_case_setup_error(
                    &mut manifest,
                    &mut manifest_error,
                    &mut failures,
                    case,
                    &output_dir,
                    &log_path,
                    &format!(
                        "input/executable state changed during sweep: expected={fixed_base_signature:016x} current={current:016x}"
                    ),
                );
                continue;
            }
            Err(error) => {
                record_case_setup_error(
                    &mut manifest,
                    &mut manifest_error,
                    &mut failures,
                    case,
                    &output_dir,
                    &log_path,
                    &format!("cannot verify fixed run signature: {error}"),
                );
                continue;
            }
        }

        if let Err(error) = fs::create_dir_all(&output_dir) {
            record_case_setup_error(
                &mut manifest,
                &mut manifest_error,
                &mut failures,
                case,
                &output_dir,
                &log_path,
                &format!("cannot create output directory: {error}"),
            );
            continue;
        }

        if marker_matches_case(&marker_path, case, &run_signature) {
            match validate_artifacts(&output_dir, &log_path, None).and_then(|artifacts| {
                if marker_matches_artifacts(&marker_path, artifacts) {
                    Ok(artifacts)
                } else {
                    Err("completion marker does not match current artifact fingerprint".to_string())
                }
            }) {
                Ok(artifacts) => {
                    let detail =
                        format!(
                        "validated resume: cor_files={} diagnostic_files={} delay_model_files={}",
                        artifacts.cor_files, artifacts.diagnostic_files, artifacts.delay_model_files
                    );
                    record_or_warn(
                        &mut manifest,
                        &mut manifest_error,
                        "case-resume",
                        case.key,
                        "ok",
                        "0",
                        0.0,
                        case.description,
                        &output_dir,
                        Some(&log_path),
                        &detail,
                    );
                    crate::emit_stdout_fmt(
                        format_args!(
                            "[model-sweep] case {}/{} resume: {} ({detail})",
                            index + 1,
                            CASES.len(),
                            case.key
                        ),
                        true,
                    );
                    continue;
                }
                Err(error) => {
                    record_or_warn(
                        &mut manifest,
                        &mut manifest_error,
                        "case-resume",
                        case.key,
                        "invalid",
                        "-",
                        0.0,
                        case.description,
                        &output_dir,
                        Some(&log_path),
                        &format!("rerunning: {error}"),
                    );
                    if let Err(remove_error) = fs::remove_file(&marker_path) {
                        record_case_setup_error(
                            &mut manifest,
                            &mut manifest_error,
                            &mut failures,
                            case,
                            &output_dir,
                            &log_path,
                            &format!("cannot remove invalid completion marker: {remove_error}"),
                        );
                        continue;
                    }
                }
            }
        } else if marker_path.exists() {
            if let Err(error) = fs::remove_file(&marker_path) {
                record_case_setup_error(
                    &mut manifest,
                    &mut manifest_error,
                    &mut failures,
                    case,
                    &output_dir,
                    &log_path,
                    &format!("cannot remove mismatched completion marker: {error}"),
                );
                continue;
            }
        }

        if let Err(error) = clear_previous_case_products(&output_dir) {
            record_case_setup_error(
                &mut manifest,
                &mut manifest_error,
                &mut failures,
                case,
                &output_dir,
                &log_path,
                &format!("cannot clear products from previous attempt: {error}"),
            );
            continue;
        }
        crate::emit_stdout_fmt(
            format_args!(
                "[model-sweep] case {}/{} start: {} ({}) output={} log={}",
                index + 1,
                CASES.len(),
                case.key,
                case.description,
                output_dir.display(),
                log_path.display()
            ),
            true,
        );
        record_or_warn(
            &mut manifest,
            &mut manifest_error,
            "case-start",
            case.key,
            "running",
            "-",
            0.0,
            case.description,
            &output_dir,
            Some(&log_path),
            &format!("child process starting run_signature={run_signature}"),
        );

        let case_started_wall = SystemTime::now();
        let started = Instant::now();
        let mut log_file = match OpenOptions::new().create(true).append(true).open(&log_path) {
            Ok(file) => file,
            Err(error) => {
                record_case_setup_error(
                    &mut manifest,
                    &mut manifest_error,
                    &mut failures,
                    case,
                    &output_dir,
                    &log_path,
                    &format!("cannot create case log: {error}"),
                );
                continue;
            }
        };
        let _ = writeln!(
            log_file,
            "# yi-corr model sweep case={} description={} run_signature={} started_unix={}",
            case.key,
            case.description,
            run_signature,
            unix_timestamp()
        );
        let _ = log_file.flush();
        let stderr_file = match log_file.try_clone() {
            Ok(file) => file,
            Err(error) => {
                record_case_setup_error(
                    &mut manifest,
                    &mut manifest_error,
                    &mut failures,
                    case,
                    &output_dir,
                    &log_path,
                    &format!("cannot clone case log: {error}"),
                );
                continue;
            }
        };

        let mut command = Command::new(&executable);
        command
            .args(&forwarded)
            .arg("--cor-directory")
            .arg(&output_dir)
            .arg("--model-diagnostics")
            .stdout(Stdio::from(log_file))
            .stderr(Stdio::from(stderr_file));
        configure_case(&mut command, case);
        let status = command.status();
        let elapsed_s = started.elapsed().as_secs_f64();

        let post_signature_error =
            match base_run_signature(&executable, &forwarded, args) {
                Ok(current) if current == fixed_base_signature => None,
                Ok(current) => Some(format!(
                    "input/executable state changed while case ran: expected={fixed_base_signature:016x} current={current:016x}"
                )),
                Err(error) => Some(format!(
                    "cannot verify input/executable state after case: {error}"
                )),
            };
        if let Some(error) = post_signature_error {
            let child_exit = status
                .as_ref()
                .ok()
                .and_then(|child_status| child_status.code())
                .unwrap_or(-1);
            failures.push(format!("{}(state changed)", case.key));
            append_case_log(
                &log_path,
                &format!("# sweep state validation failed: {error}"),
            );
            record_or_warn(
                &mut manifest,
                &mut manifest_error,
                "case-end",
                case.key,
                "state-changed",
                &child_exit.to_string(),
                elapsed_s,
                case.description,
                &output_dir,
                Some(&log_path),
                &error,
            );
            crate::emit_stdout_fmt(
                format_args!(
                    "[model-sweep] case {}/{} invalidated because inputs or executable changed: {}",
                    index + 1,
                    CASES.len(),
                    case.key
                ),
                true,
            );
            continue;
        }
        match status {
            Ok(status) if status.success() => {
                append_case_log(
                    &log_path,
                    &format!(
                        "# model sweep child exit=0 elapsed_s={elapsed_s:.3} completed_unix={}",
                        unix_timestamp()
                    ),
                );
                match validate_artifacts(&output_dir, &log_path, Some(case_started_wall)) {
                    Ok(artifacts) => {
                        match write_ok_marker(&output_dir, case, &run_signature, artifacts) {
                            Ok(()) => {
                                let detail = format!(
                                    "cor_files={} diagnostic_files={} delay_model_files={}",
                                    artifacts.cor_files,
                                    artifacts.diagnostic_files,
                                    artifacts.delay_model_files
                                );
                                record_or_warn(
                                    &mut manifest,
                                    &mut manifest_error,
                                    "case-end",
                                    case.key,
                                    "ok",
                                    &status.code().unwrap_or(0).to_string(),
                                    elapsed_s,
                                    case.description,
                                    &output_dir,
                                    Some(&log_path),
                                    &detail,
                                );
                                crate::emit_stdout_fmt(
                                    format_args!(
                                    "[model-sweep] case {}/{} complete: {} elapsed={:.3}s {detail}",
                                    index + 1,
                                    CASES.len(),
                                    case.key,
                                    elapsed_s
                                ),
                                    true,
                                );
                            }
                            Err(error) => {
                                failures.push(format!("{}(ok marker: {error})", case.key));
                                record_or_warn(
                                    &mut manifest,
                                    &mut manifest_error,
                                    "case-end",
                                    case.key,
                                    "marker-error",
                                    "0",
                                    elapsed_s,
                                    case.description,
                                    &output_dir,
                                    Some(&log_path),
                                    &error.to_string(),
                                );
                            }
                        }
                    }
                    Err(error) => {
                        failures.push(format!("{}(artifacts: {error})", case.key));
                        append_case_log(
                            &log_path,
                            &format!("# artifact validation failed: {error}"),
                        );
                        record_or_warn(
                            &mut manifest,
                            &mut manifest_error,
                            "case-end",
                            case.key,
                            "artifact-error",
                            "0",
                            elapsed_s,
                            case.description,
                            &output_dir,
                            Some(&log_path),
                            &error,
                        );
                        crate::emit_stdout_fmt(
                            format_args!(
                                "[model-sweep] case {}/{} produced incomplete artifacts but sweep continues: {} error={error}",
                                index + 1,
                                CASES.len(),
                                case.key
                            ),
                            true,
                        );
                    }
                }
            }
            Ok(status) => {
                let code = status.code().unwrap_or(-1);
                failures.push(format!("{}(exit={code})", case.key));
                append_case_log(
                    &log_path,
                    &format!("# model sweep child failed exit={code} elapsed_s={elapsed_s:.3}"),
                );
                record_or_warn(
                    &mut manifest,
                    &mut manifest_error,
                    "case-end",
                    case.key,
                    "failed",
                    &code.to_string(),
                    elapsed_s,
                    case.description,
                    &output_dir,
                    Some(&log_path),
                    "child returned non-zero; sweep continues",
                );
                crate::emit_stdout_fmt(
                    format_args!(
                        "[model-sweep] case {}/{} failed but sweep continues: {} exit={} elapsed={:.3}s",
                        index + 1,
                        CASES.len(),
                        case.key,
                        code,
                        elapsed_s
                    ),
                    true,
                );
            }
            Err(error) => {
                failures.push(format!("{}({error})", case.key));
                append_case_log(&log_path, &format!("# child spawn error: {error}"));
                record_or_warn(
                    &mut manifest,
                    &mut manifest_error,
                    "case-end",
                    case.key,
                    "spawn-error",
                    "-1",
                    elapsed_s,
                    case.description,
                    &output_dir,
                    Some(&log_path),
                    &error.to_string(),
                );
                crate::emit_stdout_fmt(
                    format_args!(
                        "[model-sweep] case {}/{} could not start but sweep continues: {} error={}",
                        index + 1,
                        CASES.len(),
                        case.key,
                        error
                    ),
                    true,
                );
            }
        }
    }

    if let Some(error) = manifest_error.as_ref() {
        failures.push(format!("manifest({error})"));
    }
    let final_status = if failures.is_empty() { "ok" } else { "failed" };
    record_or_warn(
        &mut manifest,
        &mut manifest_error,
        "sweep-end",
        "-",
        final_status,
        if failures.is_empty() { "0" } else { "1" },
        sweep_started.elapsed().as_secs_f64(),
        "unattended five-case model sweep",
        &root,
        None,
        &format!("failures={}", failures.join(", ")),
    );
    if failures.is_empty() {
        if let Some(error) = manifest_error.as_ref() {
            failures.push(format!("manifest({error})"));
        }
    }
    crate::emit_stdout_fmt(
        format_args!(
            "[model-sweep] sweep finished: failures={} manifest={}",
            failures.len(),
            manifest_path.display()
        ),
        true,
    );
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "model sweep completed with failures: {}",
            failures.join(", ")
        )
        .into())
    }
}

#[cfg(test)]
mod tests {
    use super::{env_flag_value, forwarded_arguments_from, recursion_guard};
    use std::ffi::{OsStr, OsString};

    fn strings(values: &[&str]) -> Vec<OsString> {
        values.iter().map(OsString::from).collect()
    }

    #[test]
    fn inferred_sweep_diagnostic_and_output_arguments_are_removed() {
        let input = strings(&[
            "--model-s",
            "--model-d",
            "--std",
            "--compa",
            "--cor-dir",
            "old-a",
            "--out=old-b",
            "--cor",
            "old-c",
            "--schedule",
            "scan.xml",
            "--length=3600",
        ]);
        assert_eq!(
            forwarded_arguments_from(input),
            strings(&["--schedule", "scan.xml", "--length=3600"])
        );
    }

    #[test]
    fn process_index_is_preserved_in_separate_and_attached_forms() {
        let input = strings(&[
            "--process-index",
            "7",
            "--process-i=8",
            "--cor-directory",
            "old",
        ]);
        assert_eq!(
            forwarded_arguments_from(input),
            strings(&["--process-index", "7", "--process-i=8"])
        );
    }

    #[test]
    fn unrelated_model_and_cpu_options_are_preserved() {
        let input = strings(&[
            "--model-time-offset",
            "0.5",
            "--cpu",
            "12",
            "--clock-test-value",
        ]);
        assert_eq!(forwarded_arguments_from(input.clone()), input);
    }

    #[test]
    fn recursive_sweep_is_rejected_only_for_marked_children() {
        assert!(recursion_guard(true, Some(OsStr::new("1"))).is_err());
        assert!(recursion_guard(true, Some(OsStr::new("true"))).is_err());
        assert!(recursion_guard(true, Some(OsStr::new("0"))).is_ok());
        assert!(recursion_guard(false, Some(OsStr::new("1"))).is_ok());
        assert!(recursion_guard(true, None).is_ok());
    }

    #[test]
    fn child_marker_uses_normal_boolean_spelling() {
        assert!(env_flag_value(Some(OsStr::new("yes"))));
        assert!(!env_flag_value(Some(OsStr::new("off"))));
        assert!(!env_flag_value(None));
    }
}
