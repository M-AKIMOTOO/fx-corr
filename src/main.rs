mod acf;
mod affinity;
mod args;
mod cor;
mod geom;
mod ifile;
mod plot;
mod utils;
mod xcf;
mod xml;

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::{read_to_string, File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex, OnceLock};
use std::thread;

use clap::Parser;
use num_complex::Complex;
use rayon::prelude::*;

use acf::finalize_auto_spectrum;
use args::{parse_levels, parse_shuffle, resolve_per_antenna_config, resolve_weight, DEFAULT_SHUFFLE_IN};
use cor::{epoch_to_yyyydddhhmmss, unix_seconds_to_yyyydddhhmmss, CorHeaderConfig, CorStation, CorWriter};
use plot::{plot_multi_series_f64_x, plot_series_f64_x, plot_series_with_x, BLUE, GREEN, RED};
use utils::{apply_delay_and_rate_regular_bins, apply_delay_and_rate_regular_bins_range, apply_integer_sample_shift_zerofill, build_decode_plan, decode_block_into_with_plan, quantise_frame, DynError, FftHelper, FftScratch};
use xcf::finalize_cross_spectrum;

const DEFAULT_COARSE_DELAY_S: f64 = 0.0;
// GICO3-compatible one-sided spectrum normalization factor:
// .cor stores only positive-frequency bins for real-input FFT outputs.
const COR_ONE_SIDED_POWER_FACTOR: f64 = 0.5;

fn carrier_phase_from_delay(f_hz: f64, tau_s: f64) -> Complex<f32> {
    if f_hz == 0.0 {
        Complex::new(1.0_f32, 0.0_f32)
    } else {
        Complex::from_polar(1.0_f32, (-2.0 * std::f64::consts::PI * f_hz * tau_s) as f32)
    }
}

fn delay_seconds_at_time(d_s: f64, r_sps: f64, a_sps2: f64, t_s: f64) -> f64 {
    d_s + r_sps * t_s + 0.5 * a_sps2 * t_s.powi(2)
}

fn split_delay_to_integer_and_fractional(delay_seconds: f64, fs_hz: f64) -> (i64, f64) {
    let integer_samples = (delay_seconds * fs_hz).round() as i64;
    let fractional_seconds = delay_seconds - (integer_samples as f64 / fs_hz);
    (integer_samples, fractional_seconds)
}

#[derive(Clone, Copy)]
struct FrameDelayEntry {
    int1: i64,
    int2: i64,
    frac1: f64,
    frac2: f64,
    fr_lo1: Complex<f32>,
    fr_lo2: Complex<f32>,
}

struct DelayEvalConfig {
    frame_dt: f64,
    fs: f64,
    time_offset_s: f64,
    geom_delay_table_1s: Option<Arc<[f64]>>,
    coarse_delay_s: f64,
    delay_user_samples: f64,
    extra_delay_rate_sps: f64,
    clock1_delay_s: f64,
    clock1_rate_sps: f64,
    clock2_delay_s: f64,
    clock2_rate_sps: f64,
    d_seek: f64,
    net_d1_base: f64,
    total_rate1_base: f64,
    total_accel1_base: f64,
    net_d2_base: f64,
    total_rate2_base: f64,
    total_accel2_base: f64,
    lo1_hz: f64,
    lo2_hz: f64,
}

fn compute_frame_delay_entry(frame_idx: usize, cfg: &DelayEvalConfig) -> FrameDelayEntry {
    let t_mid_local = (frame_idx as f64 + 0.5) * cfg.frame_dt;
    let t_mid = t_mid_local + cfg.time_offset_s;
    let (tau1, tau2) = if let Some(gd_table) = cfg.geom_delay_table_1s.as_deref() {
        let sec_idx = t_mid_local.floor().max(0.0) as usize;
        let i0 = sec_idx.min(gd_table.len().saturating_sub(1));
        let i1 = (i0 + 1).min(gd_table.len().saturating_sub(1));
        let frac = if i1 == i0 {
            0.0
        } else {
            (t_mid_local - i0 as f64).clamp(0.0, 1.0)
        };
        let gd_t = gd_table[i0] * (1.0 - frac) + gd_table[i1] * frac;
        let net_d_rel_no_clock_t = gd_t
            + cfg.coarse_delay_s
            + cfg.delay_user_samples / cfg.fs
            + cfg.extra_delay_rate_sps * t_mid;
        let clock1_t = cfg.clock1_delay_s + cfg.clock1_rate_sps * t_mid;
        let clock2_t = cfg.clock2_delay_s + cfg.clock2_rate_sps * t_mid;
        // Apply only relative clock delay/rate to avoid large common-mode shifts.
        // Large absolute per-station clocks can exceed small FFT frame length and
        // zero-fill whole frames even though the baseline-relative delay is small.
        let rel_clock_t = clock2_t - clock1_t;
        (net_d_rel_no_clock_t + rel_clock_t - cfg.d_seek, 0.0)
    } else {
        (
            delay_seconds_at_time(
                cfg.net_d1_base,
                cfg.total_rate1_base,
                cfg.total_accel1_base,
                t_mid,
            ),
            delay_seconds_at_time(
                cfg.net_d2_base,
                cfg.total_rate2_base,
                cfg.total_accel2_base,
                t_mid,
            ),
        )
    };
    let (int1, frac1) = split_delay_to_integer_and_fractional(tau1, cfg.fs);
    let (int2, frac2) = split_delay_to_integer_and_fractional(tau2, cfg.fs);
    FrameDelayEntry {
        int1,
        int2,
        frac1,
        frac2,
        fr_lo1: carrier_phase_from_delay(cfg.lo1_hz, tau1),
        fr_lo2: carrier_phase_from_delay(cfg.lo2_hz, tau2),
    }
}

#[derive(Clone, Copy, Debug)]
enum OutputGrid { Ant1, Ant2 }

#[derive(Clone, Copy, Debug)]
enum RunMode {
    Acf,
    Xcf,
    PhasedArray,
    Corr,
}

impl RunMode {
    fn label(self) -> &'static str {
        match self {
            RunMode::Acf => "yi_acf",
            RunMode::Xcf => "yi-xcf",
            RunMode::PhasedArray => "yi-phasedarray",
            RunMode::Corr => "yi-corr",
        }
    }
}

fn detect_run_mode_from_argv0() -> RunMode {
    let stem = std::env::args()
        .next()
        .and_then(|p| PathBuf::from(p).file_stem().map(|s| s.to_string_lossy().to_string()))
        .unwrap_or_default()
        .to_ascii_lowercase();
    match stem.as_str() {
        "yi_acf" | "yi-acf" => RunMode::Acf,
        "yi-xcf" | "yi_xcf" => RunMode::Xcf,
        "yi_phasedarray" | "yi-phasedarray" => RunMode::PhasedArray,
        "yi-corr" | "yi_corr" => RunMode::Corr,
        _ => RunMode::Corr,
    }
}

fn parse_cache_size_bytes(text: &str) -> Option<u64> {
    let s = text.trim();
    if s.is_empty() {
        return None;
    }
    let (num, mul) = if let Some(v) = s.strip_suffix('K') {
        (v, 1024u64)
    } else if let Some(v) = s.strip_suffix('M') {
        (v, 1024u64 * 1024u64)
    } else if let Some(v) = s.strip_suffix('G') {
        (v, 1024u64 * 1024u64 * 1024u64)
    } else {
        (s, 1u64)
    };
    num.trim().parse::<u64>().ok().map(|v| v.saturating_mul(mul))
}

fn detect_l3_cache_bytes() -> Option<u64> {
    let paths = [
        "/sys/devices/system/cpu/cpu0/cache/index3/size",
        "/sys/devices/system/cpu/cpu0/cache/index2/size",
    ];
    for p in paths {
        if let Ok(txt) = read_to_string(p) {
            if let Some(bytes) = parse_cache_size_bytes(&txt) {
                if bytes > 0 {
                    return Some(bytes);
                }
            }
        }
    }
    None
}

fn auto_chunk_frames(cpu_threads: usize, bytes_per_frame_pair: usize) -> usize {
    if bytes_per_frame_pair == 0 {
        return 4096;
    }
    let l3_bytes = detect_l3_cache_bytes().unwrap_or(32 * 1024 * 1024);
    let scale = (cpu_threads as u64).clamp(1, 16);
    let target_bytes = (l3_bytes.saturating_mul(scale) / 4).clamp(8 * 1024 * 1024, 128 * 1024 * 1024);
    let mut frames = (target_bytes / bytes_per_frame_pair as u64) as usize;
    frames = frames.clamp(1024, 32768);
    ((frames / 256).max(1)) * 256
}

fn auto_pipeline_depth(cpu_threads: usize) -> usize {
    (cpu_threads / 4).clamp(2, 16)
}

fn schedule_stdout_txt_path(schedule: &PathBuf) -> PathBuf {
    let mut os: OsString = schedule.as_os_str().to_os_string();
    os.push(".stdout.txt");
    PathBuf::from(os)
}

fn parse_ifile_cached(path: &PathBuf) -> Result<Arc<ifile::IFileData>, DynError> {
    static IF_CACHE: OnceLock<Mutex<HashMap<PathBuf, Arc<ifile::IFileData>>>> = OnceLock::new();
    let cache = IF_CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Ok(guard) = cache.lock() {
        if let Some(found) = guard.get(path) {
            return Ok(Arc::clone(found));
        }
    }
    let parsed = Arc::new(ifile::parse_ifile(path)?);
    let mut guard = cache
        .lock()
        .map_err(|_| std::io::Error::other("ifile cache lock poisoned while inserting"))?;
    guard.insert(path.clone(), Arc::clone(&parsed));
    Ok(parsed)
}

fn chrono_like_now_string() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let Ok(now) = SystemTime::now().duration_since(UNIX_EPOCH) else {
        return "unix:0".to_string();
    };
    format!("unix:{}", now.as_secs())
}

#[derive(Clone, Copy, Debug)]
struct BandAlignment { shift_bins: isize, a1s: usize, a1e: usize, a2s: usize }

fn compute_band_alignment(fft: usize, fs: f64, c1: f64, c2: f64, bw1: f64, bw2: f64) -> Result<BandAlignment, DynError> {
    let df = fs / fft as f64; let h = fft / 2 + 1;
    let lim = if fft % 2 == 0 { h - 1 } else { h };
    let v1 = (((bw1 * 1e6) / df).floor() as usize).saturating_add(1).min(lim);
    let v2 = (((bw2 * 1e6) / df).floor() as usize).saturating_add(1).min(lim);
    let sr = (c1 * 1e6 - c2 * 1e6) / df; let sb = sr.round() as isize;
    let a1s = 0isize.max(-sb) as usize; let a1e = (v1 as isize).min(v2 as isize - sb) as usize;
    if a1e <= a1s { return Err("No band overlap".into()); }
    Ok(BandAlignment { shift_bins: sb, a1s, a1e, a2s: (a1s as isize + sb) as usize })
}

fn format_bit_codes(b: usize) -> String { if b == 0 { "n/a".into() } else { (0..(1<<b).min(1024)).map(|c| format!("{c:0b$b}", b=b)).collect::<Vec<_>>().join(" ") } }
fn format_level_map(b: usize, l: &[f64]) -> String { if b == 0 { "n/a".into() } else { (0..(1<<b).min(l.len())).map(|c| format!("{c:0b$b}->{:.1}", l[c], b=b)).collect::<Vec<_>>().join(", ") } }
fn format_shuffle_compact(values: &[usize]) -> String { values.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",") }

fn quantization_mean_power(levels: &[f64], label: &str) -> Result<f64, DynError> {
    if levels.is_empty() {
        return Err(format!("{label} levels are empty").into());
    }
    let p = levels.iter().map(|v| v * v).sum::<f64>() / levels.len() as f64;
    if !p.is_finite() || p <= 0.0 {
        return Err(format!("{label} mean-square level is invalid: {p}").into());
    }
    Ok(p)
}

fn normalize_level_args(level_args: &[String], bit1: usize, bit2: usize) -> Result<Vec<String>, DynError> {
    if level_args.is_empty() {
        return Ok(level_args.to_vec());
    }
    if level_args
        .iter()
        .any(|v| v.contains("ant1:") || v.contains("ant2:"))
    {
        return Ok(level_args.to_vec());
    }
    if bit1 != bit2 {
        return Ok(level_args.to_vec());
    }
    let expected = 1usize << bit1;
    let all_numeric = level_args.iter().all(|v| v.parse::<f64>().is_ok());
    if all_numeric {
        if level_args.len() != expected {
            return Err(format!(
                "--level without ant1:/ant2: must provide exactly {expected} values for --bin {bit1}"
            )
            .into());
        }
        return Ok(vec![level_args.join(",")]);
    }
    Ok(level_args.to_vec())
}

fn normalize_shuffle_args(shuffle_args: &[String]) -> Result<Vec<String>, DynError> {
    if shuffle_args.is_empty() {
        return Ok(shuffle_args.to_vec());
    }
    if shuffle_args
        .iter()
        .any(|v| v.contains("ant1:") || v.contains("ant2:"))
    {
        return Ok(shuffle_args.to_vec());
    }
    let all_numeric = shuffle_args.iter().all(|v| v.parse::<usize>().is_ok());
    if all_numeric {
        if shuffle_args.len() != 32 {
            return Err("--shuffle without ant1:/ant2: must provide exactly 32 values".into());
        }
        return Ok(vec![shuffle_args.join(",")]);
    }
    Ok(shuffle_args.to_vec())
}

fn seek_forward_samples(r: &mut BufReader<File>, s: u64, b: usize, len: u64, sought: &mut isize, label: &str, reason: &str) -> Result<(), DynError> {
    let bytes = (s * b as u64 + 7) / 8;
    if bytes >= len { return Err(format!("Seek for {} exceeds file", label).into()); }
    r.seek(SeekFrom::Current(bytes as i64))?;
    *sought += ((bytes * 8) / b as u64) as isize;
    println!("[info] Seeking {} forward by {} bytes for {}", label, bytes, reason);
    Ok(())
}

fn read_with_padding(reader: &mut BufReader<File>, buf: &mut [u8]) -> Result<usize, DynError> {
    let mut total = 0;
    while total < buf.len() {
        match reader.read(&mut buf[total..]) {
            Ok(0) => { buf[total..].fill(0); return Ok(total); }
            Ok(n) => total += n,
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    Ok(total)
}

#[cfg(any(target_os = "linux", target_os = "android"))]
fn advise_sequential_readahead(file: &File) {
    use std::os::fd::AsRawFd;
    let fd = file.as_raw_fd();
    unsafe {
        let _ = libc::posix_fadvise(fd, 0, 0, libc::POSIX_FADV_SEQUENTIAL);
        let _ = libc::posix_fadvise(fd, 0, 0, libc::POSIX_FADV_WILLNEED);
    }
}

#[cfg(not(any(target_os = "linux", target_os = "android")))]
fn advise_sequential_readahead(_file: &File) {}

fn apply_integer_delay_with_forward_seek(ds: i64, b1: usize, b2: usize, res: &str, r1: &mut BufReader<File>, r2: &mut BufReader<File>, l1: u64, l2: u64, s1: &mut isize, s2: &mut isize) -> Result<(), DynError> {
    if ds == 0 { return Ok(()); }
    if ds > 0 { seek_forward_samples(r2, ds as u64, b2, l2, s2, "ant2", res)?; }
    else { seek_forward_samples(r1, (-ds) as u64, b1, l1, s1, "ant1", res)?; }
    Ok(())
}

fn resolve_input_paths(args: &args::Args, pe: &str, meta: Option<&ifile::IFileData>) -> Result<(PathBuf, PathBuf, String, String), DynError> {
    if let (Some(a1), Some(a2)) = (&args.ant1, &args.ant2) { let (_, tag) = epoch_to_yyyydddhhmmss(pe)?; return Ok((a1.clone(), a2.clone(), pe.to_string(), tag)); }
    let data_dir = args.raw_directory.clone().unwrap_or_else(|| PathBuf::from("."));
    let (_, tag) = epoch_to_yyyydddhhmmss(pe)?;
    let mut candidates = vec![("YAMAGU32", "YAMAGU34")];
    if let Some(m) = meta { if let (Some(n1), Some(n2)) = (&m.ant1_station_name, &m.ant2_station_name) { candidates.insert(0, (n1, n2)); } }
    for (p1, p2) in candidates {
        let a1 = data_dir.join(format!("{}_{}.raw", p1, tag)); let a2 = data_dir.join(format!("{}_{}.raw", p2, tag));
        if a1.exists() && a2.exists() { println!("[info] Auto-resolved inputs: {} / {}", a1.display(), a2.display()); return Ok((a1, a2, pe.to_string(), tag)); }
    }
    Err("Input files not found".into())
}

fn resolve_output_layout(args: &args::Args, tag: &str, run_mode: RunMode) -> Result<(PathBuf, PathBuf), DynError> {
    let stem = format!("YAMAGU66_{tag}");
    let dir = if matches!(run_mode, RunMode::PhasedArray) {
        std::env::current_dir()?
    } else {
        let d = args
            .cor_directory
            .clone()
            .ok_or("--cor-directory is required for yi-acf/yi-xcf/yi-corr")?;
        std::fs::create_dir_all(&d)?;
        d
    };
    let path = dir.join(format!("{stem}.raw"));
    Ok((dir, path))
}

fn fmt_opt<T: std::fmt::Display>(v: Option<T>) -> String {
    v.map(|x| x.to_string()).unwrap_or_else(|| "-".to_string())
}

fn fmt_opt_f64(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.9e}")).unwrap_or_else(|| "-".to_string())
}

fn fmt_opt_mhz(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.3} MHz")).unwrap_or_else(|| "-".to_string())
}

fn fmt_mhz(v: f64) -> String { format!("{v:.3} MHz") }
fn fmt_hz_to_mhz(v_hz: f64) -> String { format!("{:.3} MHz", v_hz / 1e6) }
fn sanitize_file_token(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for c in raw.chars() {
        if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
            out.push(c);
        } else {
            out.push('_');
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        "unknown".to_string()
    } else {
        trimmed.to_string()
    }
}

fn fmt_opt_range_mhz(center_mhz: Option<f64>, bw_mhz: Option<f64>) -> String {
    match (center_mhz, bw_mhz) {
        (Some(c), Some(bw)) => {
            let low = c - 0.5 * bw;
            let high = c + 0.5 * bw;
            format!("{low:.3}..{high:.3} MHz")
        }
        _ => "-".to_string(),
    }
}

struct ScheduleGlobalView<'a> {
    source_frame: &'a str,
    process_skip_sec: Option<f64>,
    process_length_sec: Option<f64>,
    geom_delay_s: f64,
    geom_rate_sps: f64,
    geom_accel_sps2: f64,
    geom_rate_hz_at_obs: f64,
    rel_clock_delay_s: f64,
    rel_clock_rate_sps: f64,
    coarse_delay_s: f64,
    coarse_delay_samples: f64,
    res_delay_input_samples: f64,
    read_align_delay_samples: f64,
    read_align_delay_s: f64,
    delay_model: &'a str,
    delay_rate_hz: f64,
    delay_rate_geom_hz: f64,
    delay_rate_clock_hz: f64,
    delay_rate_rot_hz: f64,
    delay_rate_user_hz: f64,
    obsfreq_mhz: f64,
    bw_mhz: f64,
    sampling_mhz: f64,
    band_align_desc: String,
    band_overlap_bins: usize,
    band_overlap_mhz: f64,
    band_bin_mhz: f64,
    rotation_shift_desc: String,
    rotation_delta_hz: f64,
    rotation_shift_hz: f64,
    rotation_residual_hz: f64,
    output_grid: &'a str,
}

struct ScheduleAntennaView<'a> {
    name: &'a str,
    key: &'a str,
    path: &'a PathBuf,
    size_bytes: u64,
    est_obs_sec: f64,
    ecef_m: [f64; 3],
    bit: usize,
    levels: &'a [f64],
    shuffle_ext: &'a [usize],
    sideband: &'a str,
    rotation_hz: f64,
    clock_delay_s: f64,
    clock_rate_sps: f64,
    center_mhz: Option<f64>,
    bw_mhz: Option<f64>,
}

fn print_kv_table(title: &str, rows: &[(String, String)]) {
    let mut w1 = "Field".len();
    let mut w2 = "Value".len();
    for (k, v) in rows {
        w1 = w1.max(k.len());
        w2 = w2.max(v.len());
    }
    let border = format!("+{}+{}+", "-".repeat(w1 + 2), "-".repeat(w2 + 2));
    println!("{title}");
    println!("{border}");
    println!("| {:<w1$} | {:<w2$} |", "Field", "Value", w1 = w1, w2 = w2);
    println!("{border}");
    for (k, v) in rows {
        println!("| {:<w1$} | {:<w2$} |", k, v, w1 = w1, w2 = w2);
    }
    println!("{border}");
}

fn print_ant_table(header_l: &str, header_r: &str, rows: &[(String, String, String)]) {
    let mut w0 = "Parameter".len();
    let mut w1 = header_l.len();
    let mut w2 = header_r.len();
    for (p, l, r) in rows {
        w0 = w0.max(p.len());
        w1 = w1.max(l.len());
        w2 = w2.max(r.len());
    }
    let border = format!(
        "+{}+{}+{}+",
        "-".repeat(w0 + 2),
        "-".repeat(w1 + 2),
        "-".repeat(w2 + 2)
    );
    println!("Antenna Parameters");
    println!("{border}");
    println!(
        "| {:<w0$} | {:<w1$} | {:<w2$} |",
        "Parameter",
        header_l,
        header_r,
        w0 = w0,
        w1 = w1,
        w2 = w2
    );
    println!("{border}");
    for (p, l, r) in rows {
        println!(
            "| {:<w0$} | {:<w1$} | {:<w2$} |",
            p,
            l,
            r,
            w0 = w0,
            w1 = w1,
            w2 = w2
        );
    }
    println!("{border}");
}

fn print_schedule_summary(
    path: &PathBuf,
    d: &ifile::IFileData,
    gv: &ScheduleGlobalView<'_>,
    a1: &ScheduleAntennaView<'_>,
    a2: &ScheduleAntennaView<'_>,
    include_global: bool,
) {
    println!("[info] Parsed schedule XML parameters:");
    if include_global {
        let global_rows = vec![
            ("file".to_string(), path.display().to_string()),
            (
                "source".to_string(),
                d.source.as_deref().unwrap_or("-").to_string(),
            ),
            (
                "stream label".to_string(),
                d.stream_label.as_deref().unwrap_or("-").to_string(),
            ),
            ("source frame".to_string(), gv.source_frame.to_string()),
            ("ra/dec (J2000)".to_string(), format!("{} / {}", d.ra, d.dec)),
            (
                "process epochs (all)".to_string(),
                if d.process_epochs.is_empty() {
                    "-".to_string()
                } else {
                    d.process_epochs.join(", ")
                },
            ),
            (
                "process skip (xml) [s]".to_string(),
                gv.process_skip_sec
                    .map(|v| format!("{v:.3}"))
                    .unwrap_or_else(|| "-".to_string()),
            ),
            (
                "process length (xml) [s]".to_string(),
                gv.process_length_sec
                    .map(|v| format!("{v:.3}"))
                    .unwrap_or_else(|| "-".to_string()),
            ),
            ("delay convention".to_string(), "apply delay correction as ant2-ant1".to_string()),
            (
                "geom delay [s]".to_string(),
                format!("{:.6e}", gv.geom_delay_s),
            ),
            (
                "geom rate [s/s]".to_string(),
                format!("{:.6e}", gv.geom_rate_sps),
            ),
            (
                "geom accel [s/s^2]".to_string(),
                format!("{:.6e}", gv.geom_accel_sps2),
            ),
            (
                "geom rate @obsfreq [Hz]".to_string(),
                format!("{:.6}", gv.geom_rate_hz_at_obs),
            ),
            (
                "clock relative delay [s]".to_string(),
                format!("{:.6e}", gv.rel_clock_delay_s),
            ),
            (
                "clock relative rate [s/s]".to_string(),
                format!("{:.6e}", gv.rel_clock_rate_sps),
            ),
            (
                "coarse delay fixed [s]".to_string(),
                format!("{:.6e}", gv.coarse_delay_s),
            ),
            (
                "coarse delay fixed [samples]".to_string(),
                format!("{:.3}", gv.coarse_delay_samples),
            ),
            (
                "res-delay input [samples]".to_string(),
                format!("{:.3}", gv.res_delay_input_samples),
            ),
            (
                "read-align delay [samples]".to_string(),
                format!("{:.3}", gv.read_align_delay_samples),
            ),
            (
                "read-align delay [s]".to_string(),
                format!("{:.6e}", gv.read_align_delay_s),
            ),
            ("delay model".to_string(), gv.delay_model.to_string()),
            (
                "delay-rate total [Hz]".to_string(),
                format!("{:.6}", gv.delay_rate_hz),
            ),
            (
                "delay-rate terms [Hz]".to_string(),
                format!(
                    "geom {:.6} + clock {:.6} + rot-res {:.6} + user {:.6}",
                    gv.delay_rate_geom_hz, gv.delay_rate_clock_hz, gv.delay_rate_rot_hz, gv.delay_rate_user_hz
                ),
            ),
            ("band-align".to_string(), gv.band_align_desc.clone()),
            (
                "band-overlap".to_string(),
                format!(
                    "{} bins ({:.3} MHz, {:.3} MHz/bin)",
                    gv.band_overlap_bins, gv.band_overlap_mhz, gv.band_bin_mhz
                ),
            ),
            ("rotation-shift".to_string(), gv.rotation_shift_desc.clone()),
            (
                "rotation-fringestop [Hz]".to_string(),
                format!(
                    "delta {:.6}, shift {:.6}, residual {:.6}",
                    gv.rotation_delta_hz, gv.rotation_shift_hz, gv.rotation_residual_hz
                ),
            ),
            ("output-grid".to_string(), gv.output_grid.to_string()),
            ("stream fft".to_string(), fmt_opt(d.fft)),
            ("stream sampling (xml) [Hz]".to_string(), fmt_opt_f64(d.sampling_hz)),
            ("stream sampling (runtime)".to_string(), fmt_mhz(gv.sampling_mhz)),
            ("stream obsfreq".to_string(), fmt_opt_mhz(d.obsfreq_mhz)),
            ("stream obsfreq (runtime)".to_string(), fmt_mhz(gv.obsfreq_mhz)),
            ("stream bw (runtime)".to_string(), fmt_mhz(gv.bw_mhz)),
        ];
        print_kv_table("Global Parameters", &global_rows);
    }

    let ant1_label = format!("Ant1 [{}:{}]", a1.key, a1.name);
    let ant2_label = format!("Ant2 [{}:{}]", a2.key, a2.name);
    let ant_rows = vec![
        (
            "input raw".to_string(),
            a1.path.display().to_string(),
            a2.path.display().to_string(),
        ),
        (
            "input size [bytes]".to_string(),
            a1.size_bytes.to_string(),
            a2.size_bytes.to_string(),
        ),
        (
            "input obs time [s]".to_string(),
            format!("{:.2}", a1.est_obs_sec),
            format!("{:.2}", a2.est_obs_sec),
        ),
        (
            "ecef [m]".to_string(),
            format!("[{:.3}, {:.3}, {:.3}]", a1.ecef_m[0], a1.ecef_m[1], a1.ecef_m[2]),
            format!("[{:.3}, {:.3}, {:.3}]", a2.ecef_m[0], a2.ecef_m[1], a2.ecef_m[2]),
        ),
        (
            "bit / bit-code".to_string(),
            format!("{} / {}", a1.bit, format_bit_codes(a1.bit)),
            format!("{} / {}", a2.bit, format_bit_codes(a2.bit)),
        ),
        (
            "level / level-map".to_string(),
            format!("{:?} / {}", a1.levels, format_level_map(a1.bit, a1.levels)),
            format!("{:?} / {}", a2.levels, format_level_map(a2.bit, a2.levels)),
        ),
        (
            "shuffle-in".to_string(),
            format_shuffle_compact(a1.shuffle_ext),
            format_shuffle_compact(a2.shuffle_ext),
        ),
        (
            "sideband".to_string(),
            a1.sideband.to_string(),
            a2.sideband.to_string(),
        ),
        (
            "rotation".to_string(),
            fmt_hz_to_mhz(a1.rotation_hz),
            fmt_hz_to_mhz(a2.rotation_hz),
        ),
        (
            "clock delay/rate [s, s/s]".to_string(),
            format!("{:.6e} / {:.6e}", a1.clock_delay_s, a1.clock_rate_sps),
            format!("{:.6e} / {:.6e}", a2.clock_delay_s, a2.clock_rate_sps),
        ),
        (
            "freq range [MHz]".to_string(),
            fmt_opt_range_mhz(a1.center_mhz, a1.bw_mhz),
            fmt_opt_range_mhz(a2.center_mhz, a2.bw_mhz),
        ),
    ];
    print_ant_table(&ant1_label, &ant2_label, &ant_rows);
}

fn main() -> Result<(), DynError> {
    let args = args::Args::parse();
    if args.mkxml {
        let out = PathBuf::from("example.xml");
        xml::write_example_xml(&out)?;
        println!("[info] Wrote example schedule XML: {}", out.display());
        return Ok(());
    }
    let run_mode = detect_run_mode_from_argv0();
    match run_mode {
        RunMode::PhasedArray => {
            if args.schedule.is_none() {
                return Err("--schedule is required for yi-phasedarray".into());
            }
            if args.raw_directory.is_none() {
                return Err("--raw-directory is required for yi-phasedarray".into());
            }
        }
        RunMode::Acf | RunMode::Xcf | RunMode::Corr => {
            if args.schedule.is_none() {
                return Err("--schedule is required for yi-acf/yi-xcf/yi-corr".into());
            }
            if args.raw_directory.is_none() {
                return Err("--raw-directory is required for yi-acf/yi-xcf/yi-corr".into());
            }
            if args.cor_directory.is_none() {
                return Err("--cor-directory is required for yi-acf/yi-xcf/yi-corr".into());
            }
        }
    }
    let cpu_auto = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(2)
        .saturating_sub(2)
        .max(1);
    let cpu_threads = match args.cpu {
        Some(0) => return Err("--cpu must be >= 1".into()),
        Some(v) => v,
        None => cpu_auto,
    };
    let affinity_runtime = affinity::AffinityRuntime::from_default_file()?;
    if let Some(msg) = affinity_runtime.info() {
        println!("[info] {msg}");
    }
    let mut tp_builder = rayon::ThreadPoolBuilder::new().num_threads(cpu_threads);
    if let Some(core_ids) = affinity_runtime.worker_cores() {
        tp_builder = tp_builder.start_handler(move |thread_idx| {
            let c = core_ids[thread_idx % core_ids.len()];
            let _ = core_affinity::set_for_current(c);
        });
    }
    tp_builder
        .build_global()
        .map_err(|e| format!("failed to configure rayon thread pool: {e}"))?;
    let if_d = if let Some(p) = &args.schedule {
        let is_xml = p
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("xml"))
            .unwrap_or(false);
        if !is_xml {
            return Err("--schedule must point to a .xml file".into());
        }
        Some(parse_ifile_cached(p)?)
    } else {
        None
    };
    if let Some(meta) = if_d.as_ref() {
        let run_all_processes = meta.processes.len() > 1
            && args.process_index.is_none()
            && args.epoch.is_none()
            && args.ra.is_none()
            && args.dec.is_none()
            && args.length.is_none()
            && args.skip == 0.0;
        if run_all_processes {
            println!(
                "[info] Processing all scans sequentially: {} process entries",
                meta.processes.len()
            );
            for (idx, p) in meta.processes.iter().enumerate() {
                println!(
                    "[info] process {}/{}: epoch={} skip={} length={}",
                    idx + 1,
                    meta.processes.len(),
                    p.epoch,
                    p.skip_sec,
                    p.length_sec
                        .map(|v| format!("{v:.3}"))
                        .unwrap_or_else(|| "-".to_string())
                );
                let mut run_args = args.clone();
                run_args.process_index = Some(idx);
                run_args.compact_logs = idx > 0;
                let process_started = std::time::Instant::now();
                let run_result = run_once(run_args, run_mode, cpu_threads);
                let elapsed_sec = process_started.elapsed().as_secs_f64();
                match run_result {
                    Ok(()) => {
                        println!(
                            "[info] process {}/{} elapsed: {:.3}s",
                            idx + 1,
                            meta.processes.len(),
                            elapsed_sec
                        );
                    }
                    Err(e) => {
                        println!(
                            "[error] process {}/{} failed after {:.3}s",
                            idx + 1,
                            meta.processes.len(),
                            elapsed_sec
                        );
                        return Err(e);
                    }
                }
            }
            return Ok(());
        }
    }
    run_once(args, run_mode, cpu_threads)
}

fn run_once(args: args::Args, run_mode: RunMode, cpu_threads: usize) -> Result<(), DynError> {
    let fringe_enabled = matches!(run_mode, RunMode::PhasedArray) && args.fringe;
    if args.fringe && !matches!(run_mode, RunMode::PhasedArray) {
        println!("[info] --fringe is ignored for this mode");
    }
    let do_xcf = fringe_enabled;
    let do_synth = matches!(run_mode, RunMode::Acf | RunMode::Xcf | RunMode::PhasedArray | RunMode::Corr);
    let write_raw = matches!(run_mode, RunMode::PhasedArray);
    let write_acf_cor = matches!(run_mode, RunMode::Acf | RunMode::Corr);
    let write_xcf_cor = matches!(run_mode, RunMode::Xcf | RunMode::Corr);
    let write_phased_cor = matches!(run_mode, RunMode::PhasedArray);
    let plot_phased = matches!(run_mode, RunMode::PhasedArray);
    let if_d = if let Some(p) = &args.schedule {
        let is_xml = p
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("xml"))
            .unwrap_or(false);
        if !is_xml {
            return Err("--schedule must point to a .xml file".into());
        }
        Some(parse_ifile_cached(p)?)
    } else {
        None
    };
    let selected_process = if let Some(meta) = if_d.as_ref() {
        if let Some(idx) = args.process_index {
            Some(
                meta.processes
                    .get(idx)
                    .cloned()
                    .ok_or_else(|| format!("--process-index {} out of range (0..{})", idx, meta.processes.len().saturating_sub(1)))?,
            )
        } else {
            meta.processes.first().cloned()
        }
    } else {
        None
    };
    if let (Some(idx), Some(p)) = (args.process_index, selected_process.as_ref()) {
        println!(
            "[info] Selected process index {}: epoch={} skip={} length={}",
            idx,
            p.epoch,
            p.skip_sec,
            p.length_sec
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "-".to_string())
        );
    }
    let schedule_stdout_path = args.schedule.as_ref().map(schedule_stdout_txt_path);
    let fft_len = if args.fft == args::DEFAULT_FFT { if_d.as_ref().and_then(|d| d.fft).unwrap_or(args.fft) } else { args.fft };
    
    let bit_a = if args.bin.is_empty() { vec!["2".into()] } else { args.bin.clone() };
    let (bit1, bit2) = resolve_per_antenna_config(&bit_a, if_d.as_ref().and_then(|d| d.ant1_bit).unwrap_or(2), |s| Ok(s.parse()?))?;
    let bit_out = bit1.max(bit2);
    let level_src = args.level.clone();
    let level_args = normalize_level_args(&level_src, bit1, bit2)?;
    let (lv1_s, lv2_s) = resolve_per_antenna_config(&level_args, if_d.as_ref().and_then(|d| d.ant1_level.clone()).unwrap_or("-1.5,-0.5,0.5,1.5".into()), |s| Ok(s.to_string()))?;
    let (levels1, levels2) = (Arc::new(parse_levels(bit1, &lv1_s)?), Arc::new(parse_levels(bit2, &lv2_s)?));
    let level_power1 = quantization_mean_power(levels1.as_ref(), "ant1")?;
    let level_power2 = quantization_mean_power(levels2.as_ref(), "ant2")?;
    let ds_s = DEFAULT_SHUFFLE_IN.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
    let shuffle_src = args.shuffle_in.clone();
    let shuffle_args = normalize_shuffle_args(&shuffle_src)?;
    let (sh1_s, sh2_s) = resolve_per_antenna_config(&shuffle_args, if_d.as_ref().and_then(|d| d.ant1_shuffle.clone()).unwrap_or(ds_s), |s| Ok(s.to_string()))?;
    let (sh1, sh2) = (Arc::new(parse_shuffle(&sh1_s)?), Arc::new(parse_shuffle(&sh2_s)?));
    let sh1_ext: Vec<usize> = sh1_s
        .split(',')
        .map(|v| v.trim().parse::<usize>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("invalid ant1 shuffle display map: {e}"))?;
    let sh2_ext: Vec<usize> = sh2_s
        .split(',')
        .map(|v| v.trim().parse::<usize>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("invalid ant2 shuffle display map: {e}"))?;
    
    let sb_a = if args.sideband.is_empty() { vec!["LSB".into()] } else { args.sideband.clone() };
    let (sb1_s, sb2_s) = resolve_per_antenna_config(&sb_a, if_d.as_ref().and_then(|d| d.ant1_sideband.clone()).unwrap_or("LSB".into()), |s| Ok(s.to_uppercase()))?;
    let (lsb1, lsb2) = (sb1_s == "LSB", sb2_s == "LSB");
    let output_lsb = lsb1 && lsb2;

    let (tsys1, tsys2) = resolve_per_antenna_config(&args.tsys, 1.0, |s| Ok(s.parse()?))?;
    let (dia1, dia2) = resolve_per_antenna_config(&args.diameter, 0.0, |s| Ok(s.parse()?))?;
    let (eta1, eta2) = resolve_per_antenna_config(&args.eta, 0.65, |s| Ok(s.parse()?))?;
    let (gain1, gain2) = resolve_per_antenna_config(&args.gain, 1.0, |s| Ok(s.parse()?))?;
    let (sefd1, sefd2) = resolve_per_antenna_config(&args.sefd, 0.0, |s| Ok(s.parse()?))?;

    let fs = if_d.as_ref().and_then(|d| d.sampling_hz).unwrap_or(args.sampling * 1e6);
    let obs_mhz = if_d.as_ref().and_then(|d| d.obsfreq_mhz).unwrap_or(args.obsfreq);
    let ep_i = args
        .epoch
        .clone()
        .or_else(|| selected_process.as_ref().map(|p| p.epoch.clone()))
        .or_else(|| if_d.as_ref().and_then(|d| d.epoch.clone()))
        .unwrap_or("2000".into());
    let (a1p, a2p, _, _unused_tag) =
        resolve_input_paths(&args, &ep_i, if_d.as_ref().map(|v| &**v))?;
    let (c_unix_base, _) = epoch_to_yyyydddhhmmss(&ep_i)?;
    let ant1_ecef = if_d.as_ref().and_then(|d| d.ant1_ecef_m).unwrap_or(geom::YAMAGU32_ECEF);
    let ant2_ecef = if_d.as_ref().and_then(|d| d.ant2_ecef_m).unwrap_or(geom::YAMAGU34_ECEF);

    let mut gdi: Option<GDI> = None; struct GDI { ra: f64, dec: f64, mjd: f64 }
    let (mut gd0, mut gr0, mut ga0) = (0.0, 0.0, 0.0);
    let ra_in = args
        .ra
        .clone()
        .or_else(|| selected_process.as_ref().and_then(|p| p.ra.clone()))
        .or_else(|| if_d.as_ref().map(|d| d.ra.clone()));
    let dec_in = args
        .dec
        .clone()
        .or_else(|| selected_process.as_ref().and_then(|p| p.dec.clone()))
        .or_else(|| if_d.as_ref().map(|d| d.dec.clone()));
    if let (Some(ra_s), Some(dec_s)) = (ra_in, dec_in) {
        let (ra_raw, dec_raw, mjd) = (
            geom::parse_ra(&ra_s)?,
            geom::parse_dec(&dec_s)?,
            geom::parse_epoch_to_mjd(&ep_i)?,
        );
        // Input sky coordinates are always interpreted as J2000.
        let (ra, dec) = geom::precess_j2000_to_mean_of_date(ra_raw, dec_raw, mjd);
        let (_, _, gd, gr, ga) = geom::calculate_geometric_delay_and_derivatives(ant1_ecef, ant2_ecef, ra, dec, mjd);
        gdi = Some(GDI { ra, dec, mjd }); (gd0, gr0, ga0) = (gd, gr, ga);
    }
    let clock_delay_rel_legacy_s = if_d.as_ref().and_then(|d| d.clock_delay_s).unwrap_or(0.0);
    let clock_rate_rel_legacy_sps = if_d.as_ref().and_then(|d| d.clock_rate_sps).unwrap_or(0.0);
    let clock1_delay_s = if_d.as_ref().and_then(|d| d.ant1_clock_delay_s).unwrap_or(0.0);
    let clock2_delay_s = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_delay_s)
        .unwrap_or(clock1_delay_s + clock_delay_rel_legacy_s);
    let clock1_rate_sps = if_d.as_ref().and_then(|d| d.ant1_clock_rate_sps).unwrap_or(0.0);
    let clock2_rate_sps = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_rate_sps)
        .unwrap_or(clock1_rate_sps + clock_rate_rel_legacy_sps);
    let clock_delay_s = clock2_delay_s - clock1_delay_s;
    let clock_rate_sps = clock2_rate_sps - clock1_rate_sps;
    let coarse_delay_s = args.coarse.unwrap_or(DEFAULT_COARSE_DELAY_S);
    let delay_user_samples = args.delay + args.resdelay;
    let rate_user_hz = args.rate + args.resrate;
    let net_d_rel_no_clock0 = gd0 + coarse_delay_s + delay_user_samples / fs;
    let net_d0 = net_d_rel_no_clock0 + clock_delay_s;

    let mut rot1 = if_d.as_ref().and_then(|d| d.ant1_rotation_hz).unwrap_or(0.0);
    let mut rot2 = if_d.as_ref().and_then(|d| d.ant2_rotation_hz).unwrap_or(0.0);
    for entry in &args.rotation {
        if entry.contains("ant1:") || entry.contains("ant2:") {
            let parts: Vec<&str> = entry
                .split(|c: char| c == ',' || c == ' ')
                .filter(|s| !s.is_empty())
                .collect();
            for part in parts {
                if let Some(v) = part.strip_prefix("ant1:") {
                    rot1 = v.trim().parse::<f64>()?;
                } else if let Some(v) = part.strip_prefix("ant2:") {
                    rot2 = v.trim().parse::<f64>()?;
                }
            }
        } else {
            let v = entry.trim().parse::<f64>()?;
            rot1 = v;
            rot2 = v;
        }
    }
    let bw = fs / 2e6;
    let a1_data_low = obs_mhz + rot1/1e6; let a2_data_low = obs_mhz + rot2/1e6;
    let lo1_hz = a1_data_low * 1e6;
    let lo2_hz = a2_data_low * 1e6;
    let ba = compute_band_alignment(fft_len, fs, a1_data_low + 0.5*bw, a2_data_low + 0.5*bw, bw, bw)?;
    let out_grid = if rot1.abs() <= rot2.abs() { OutputGrid::Ant1 } else { OutputGrid::Ant2 };
    
    let bpf1 = (fft_len * bit1 + 7) / 8; let bpf2 = (fft_len * bit2 + 7) / 8;
    let bpf_o = (fft_len * bit_out + 7) / 8;
    let bytes_per_frame_pair = bpf1 + bpf2;
    let io_chunk_frames = match args.chunk_frames {
        Some(0) => return Err("--chunk-frames must be >= 1".into()),
        Some(v) => v,
        None => auto_chunk_frames(cpu_threads, bytes_per_frame_pair),
    };
    let io_pipeline_depth = match args.pipeline_depth {
        Some(0) => return Err("--pipeline-depth must be >= 1".into()),
        Some(v) => v,
        None => auto_pipeline_depth(cpu_threads),
    };
    let f1_m = std::fs::metadata(&a1p)?; let f2_m = std::fs::metadata(&a2p)?;

    if args.skip < 0.0 {
        return Err("--skip must be >= 0".into());
    }
    if let Some(length_sec) = args.length {
        if length_sec < 0.0 {
            return Err("--length must be >= 0".into());
        }
    }
    let xml_skip_sec = selected_process
        .as_ref()
        .map(|p| p.skip_sec)
        .or_else(|| if_d.as_ref().and_then(|d| d.process_skip_sec))
        .unwrap_or(0.0);
    let xml_length_sec = selected_process
        .as_ref()
        .and_then(|p| p.length_sec)
        .or_else(|| if_d.as_ref().and_then(|d| d.process_length_sec));
    if xml_skip_sec < 0.0 {
        return Err("XML process/skip must be >= 0".into());
    }
    if let Some(length_sec) = xml_length_sec {
        if length_sec < 0.0 {
            return Err("XML process/length must be >= 0".into());
        }
    }
    let cli_skip_sec = args.skip;
    let total_skip_sec = xml_skip_sec + cli_skip_sec;
    let c_unix = c_unix_base + total_skip_sec.round() as i64;
    let c_tag = unix_seconds_to_yyyydddhhmmss(c_unix)?;
    let (o_dir, o_path) = resolve_output_layout(&args, &c_tag, run_mode)?;
    let process_window_sec = args.length.or(xml_length_sec);

    let total_samples1 = (f1_m.len() * 8) / bit1 as u64;
    let total_samples2 = (f2_m.len() * 8) / bit2 as u64;
    let init_delay_samples = (net_d0 * fs).round() as i64;
    let init_seek_s1 = if init_delay_samples < 0 {
        (-init_delay_samples) as u64
    } else {
        0
    };
    let init_seek_s2 = if init_delay_samples > 0 {
        init_delay_samples as u64
    } else {
        0
    };
    let total_skip_samples = (total_skip_sec * fs).round() as u64;
    let start_s1 = init_seek_s1.saturating_add(total_skip_samples);
    let start_s2 = init_seek_s2.saturating_add(total_skip_samples);
    let avail_samples1 = total_samples1.saturating_sub(start_s1);
    let avail_samples2 = total_samples2.saturating_sub(start_s2);
    let avail_samples = avail_samples1.min(avail_samples2);
    let available_sec = avail_samples as f64 / fs;
    let requested_sec = if let Some(window_sec) = process_window_sec {
        (window_sec - total_skip_sec).max(0.0)
    } else {
        available_sec
    };
    let total_sec = requested_sec.min(available_sec);
    // Keep only complete FFT frames.
    let total_f = (total_sec * fs / fft_len as f64).floor() as usize;

    let (w1, _, _, _) = resolve_weight(tsys1, gain1, if sefd1 > 0.0 { Some(sefd1) } else { None }, if dia1 > 0.0 { Some(dia1) } else { None }, eta1, "A1")?;
    let (w2, _, _, _) = resolve_weight(tsys2, gain2, if sefd2 > 0.0 { Some(sefd2) } else { None }, if dia2 > 0.0 { Some(dia2) } else { None }, eta2, "A2")?;
    let a1_name = if_d.as_ref().and_then(|d| d.ant1_station_name.as_deref()).unwrap_or("YAMAGU32");
    let a2_name = if_d.as_ref().and_then(|d| d.ant2_station_name.as_deref()).unwrap_or("YAMAGU34");
    let a1_key_opt = if_d.as_ref().and_then(|d| d.ant1_station_key.as_deref());
    let a2_key_opt = if_d.as_ref().and_then(|d| d.ant2_station_key.as_deref());
    let a1_key = a1_key_opt.unwrap_or("-");
    let a2_key = a2_key_opt.unwrap_or("-");
    let a1_code = a1_key_opt.and_then(|k| k.as_bytes().first().copied()).or_else(|| a1_name.as_bytes().first().copied()).unwrap_or(b'1');
    let a2_code = a2_key_opt.and_then(|k| k.as_bytes().first().copied()).or_else(|| a2_name.as_bytes().first().copied()).unwrap_or(b'2');
    let schedule_mode = if_d.is_some();
    let correction_sign = -1.0f64; // Correct ant1 toward ant2 while geometric model is (ant2 - ant1).
    let geometric_rate_hz = correction_sign * gr0 * obs_mhz * 1e6;
    let clock_rate_hz = clock_rate_sps * obs_mhz * 1e6;
    let df_hz = fs / fft_len as f64;
    let rotation_delta_hz = rot1 - rot2; // ant1 - ant2
    let rotation_shift_hz = ba.shift_bins as f64 * df_hz;
    let rotation_residual_hz = rotation_delta_hz - rotation_shift_hz;
    let rotation_fringe_hz = -rotation_residual_hz; // correction applied on ant1 toward ant2
    let total_rate_hz = geometric_rate_hz + clock_rate_hz + rotation_fringe_hz + rate_user_hz;
    let overlap_bins = ba.a1e - ba.a1s;
    let bin_mhz = fs / fft_len as f64 / 1e6;
    let overlap_mhz = overlap_bins as f64 * bin_mhz;
    let band_align_desc = format!(
        "{}->{} shift {} bins, overlap {}[{}..{})",
        a2_name, a1_name, ba.shift_bins, a1_name, ba.a1s, ba.a1e
    );
    let rotation_shift_desc = format!(
        "target {}, grid {}-ref, shift {} bins",
        a2_name, a1_name, ba.shift_bins
    );
    let delay_model = "per-frame midpoint delay + integer/fractional correction";

    if schedule_mode {
        if let (Some(sch), Some(meta)) = (&args.schedule, &if_d) {
            let gv = ScheduleGlobalView {
                source_frame: "J2000 (fixed)",
                process_skip_sec: Some(xml_skip_sec),
                process_length_sec: xml_length_sec,
                geom_delay_s: gd0,
                geom_rate_sps: gr0,
                geom_accel_sps2: ga0,
                geom_rate_hz_at_obs: gr0 * obs_mhz * 1e6,
                rel_clock_delay_s: clock_delay_s,
                rel_clock_rate_sps: clock_rate_sps,
                coarse_delay_s,
                coarse_delay_samples: coarse_delay_s * fs,
                res_delay_input_samples: delay_user_samples,
                read_align_delay_samples: net_d0 * fs,
                read_align_delay_s: net_d0,
                delay_model,
                delay_rate_hz: total_rate_hz,
                delay_rate_geom_hz: geometric_rate_hz,
                delay_rate_clock_hz: clock_rate_hz,
                delay_rate_rot_hz: rotation_fringe_hz,
                delay_rate_user_hz: rate_user_hz,
                obsfreq_mhz: obs_mhz,
                bw_mhz: bw,
                sampling_mhz: fs / 1e6,
                band_align_desc: band_align_desc.clone(),
                band_overlap_bins: overlap_bins,
                band_overlap_mhz: overlap_mhz,
                band_bin_mhz: bin_mhz,
                rotation_shift_desc: rotation_shift_desc.clone(),
                rotation_delta_hz,
                rotation_shift_hz,
                rotation_residual_hz,
                output_grid: match out_grid { OutputGrid::Ant1 => a1_name, OutputGrid::Ant2 => a2_name },
            };
            let a1v = ScheduleAntennaView {
                name: a1_name,
                key: a1_key,
                path: &a1p,
                size_bytes: f1_m.len(),
                est_obs_sec: (f1_m.len() * 8) as f64 / bit1 as f64 / fs,
                ecef_m: ant1_ecef,
                bit: bit1,
                levels: &levels1,
                shuffle_ext: &sh1_ext,
                sideband: &sb1_s,
                rotation_hz: rot1,
                clock_delay_s: clock1_delay_s,
                clock_rate_sps: clock1_rate_sps,
                center_mhz: meta.ant1_center_mhz,
                bw_mhz: meta.ant1_bw_mhz,
            };
            let a2v = ScheduleAntennaView {
                name: a2_name,
                key: a2_key,
                path: &a2p,
                size_bytes: f2_m.len(),
                est_obs_sec: (f2_m.len() * 8) as f64 / bit2 as f64 / fs,
                ecef_m: ant2_ecef,
                bit: bit2,
                levels: &levels2,
                shuffle_ext: &sh2_ext,
                sideband: &sb2_s,
                rotation_hz: rot2,
                clock_delay_s: clock2_delay_s,
                clock_rate_sps: clock2_rate_sps,
                center_mhz: meta.ant2_center_mhz,
                bw_mhz: meta.ant2_bw_mhz,
            };
            print_schedule_summary(sch, meta, &gv, &a1v, &a2v, !args.compact_logs);
        }
    }

    let is_phased_mode = matches!(run_mode, RunMode::PhasedArray);
    if args.compact_logs {
        println!(
            "[info] process summary: mode={} epoch={} fft={} frames={} length={:.6}s",
            run_mode.label(),
            ep_i,
            fft_len,
            total_f,
            total_f as f64 * fft_len as f64 / fs
        );
    } else if is_phased_mode {
        println!("Starting phased array processing with the following arguments:");
    } else {
        println!("Starting correlation processing with the following arguments:");
    }
    if !args.compact_logs {
        println!("--------------------------------------------------");
    }
    if let Some(p) = schedule_stdout_path.as_ref() {
        let mut outf = OpenOptions::new()
            .create(true)
            .append(true)
            .open(p)?;
        writeln!(
            outf,
            "[{}] mode={} schedule={} process_index={} epoch={} fft={} cpu={} skip={:.6}s length_req={:.6}s length_proc={:.6}s cor_dir={}",
            chrono_like_now_string(),
            run_mode.label(),
            args.schedule
                .as_ref()
                .map(|v| v.display().to_string())
                .unwrap_or_else(|| "-".to_string()),
            args.process_index
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            ep_i,
            fft_len,
            cpu_threads,
            total_skip_sec,
            requested_sec,
            total_f as f64 * fft_len as f64 / fs,
            o_dir.display()
        )?;
        writeln!(
            outf,
            "  delay: geom={:.6e}s rate={:.6e}s/s accel={:.6e}s/s^2 total_rate_hz={:.6} read_align={:.3}samples",
            gd0, gr0, ga0, total_rate_hz, net_d0 * fs
        )?;
        writeln!(
            outf,
            "  band: overlap_bins={} overlap_mhz={:.3} bin_mhz={:.3} align_shift_bins={}",
            overlap_bins, overlap_mhz, bin_mhz, ba.shift_bins
        )?;
    }
    if !args.compact_logs {
        println!("  mode:       {}", run_mode.label());
    }
    if !schedule_mode {
        println!("  {}:  {} (Size: {} bytes, Estimated Obs Time: {:.2}s)", a1_name, a1p.display(), f1_m.len(), (f1_m.len()*8) as f64 / bit1 as f64 / fs);
        println!("  {}:  {} (Size: {} bytes, Estimated Obs Time: {:.2}s)", a2_name, a2p.display(), f2_m.len(), (f2_m.len()*8) as f64 / bit2 as f64 / fs);
    }
    if !schedule_mode {
    if let Some(info) = &gdi {
        println!("  ra/dec:     {:.6} / {:.6}", info.ra.to_degrees(), info.dec.to_degrees());
        println!("  source-frame: J2000 (fixed)");
        println!("  epoch:      {}", ep_i);
        println!("  geom-delay: {:.6e} s ({} - {})", gd0, a2_name, a1_name);
        println!("  geom-rate:  {:.6e} s/s ({} - {}) => {:.6} Hz @ obsfreq", gr0, a2_name, a1_name, gr0 * obs_mhz * 1e6);
        println!("  geom-accel: {:.6e} s/s^2 ({} - {})", ga0, a2_name, a1_name);
    }
    println!("  coarse-delay fixed: {:.6e} s (relative pre-align) => applied {:.3} samples", coarse_delay_s, coarse_delay_s*fs);
    println!("  clock-delay {}: {:.6e} s => applied {:.3} samples", a1_name, clock1_delay_s, clock1_delay_s * fs);
    println!("  clock-delay {}: {:.6e} s => applied {:.3} samples", a2_name, clock2_delay_s, clock2_delay_s * fs);
    println!("  clock-rate  {}: {:.6e} s/s", a1_name, clock1_rate_sps);
    println!("  clock-rate  {}: {:.6e} s/s", a2_name, clock2_rate_sps);
    }
    if !schedule_mode {
        println!("  res-delay input: {} samples (relative pre-align)", delay_user_samples);
        println!("  read-align delay: {:.3} samples ({:.3e} s)", net_d0 * fs, net_d0);
        println!("  delay-model: {delay_model}");
        println!(
            "  delay-rate: {:.6} Hz (geom {:.6} + clock {:.6} + rot-res {:.6} + user {:.6})",
            total_rate_hz, geometric_rate_hz, clock_rate_hz, rotation_fringe_hz, rate_user_hz
        );
    }
    if !schedule_mode {
        println!("  obsfreq:    {:.3} MHz (Reference)", obs_mhz);
        println!(
            "  rotation:   {:.3} MHz ({}), {:.3} MHz ({})",
            rot1 / 1e6,
            a1_name,
            rot2 / 1e6,
            a2_name
        );
        println!(
            "  phase-freq: {:.3} MHz ({}) / {:.3} MHz ({}) (delay/rate correction carriers)",
            a1_data_low, a1_name, a2_data_low, a2_name
        );
        println!("  bw:         {:.3} MHz", bw);
        println!(
            "  {}-param: key={} sideband={} obsfreq_mhz={:.3} rotation_mhz={:.3} ref_band_low_mhz={:.3} data_band_low_mhz={:.3} data_band_center_mhz={:.3} bw_mhz={:.3}",
            a1_name, a1_key, sb1_s, obs_mhz, rot1 / 1e6, obs_mhz, a1_data_low, a1_data_low + 0.5 * bw, bw
        );
        println!(
            "  {}-param: key={} sideband={} obsfreq_mhz={:.3} rotation_mhz={:.3} ref_band_low_mhz={:.3} data_band_low_mhz={:.3} data_band_center_mhz={:.3} bw_mhz={:.3}",
            a2_name, a2_key, sb2_s, obs_mhz, rot2 / 1e6, obs_mhz, a2_data_low, a2_data_low + 0.5 * bw, bw
        );
        println!("  sampling:   {:.0} Hz ({:.3} MHz)", fs, fs/1e6);
    }
    if !args.compact_logs {
        println!("  samples/s:  {:.0}", fs);
        println!("  frames/s:   {:.6}", fs/fft_len as f64);
        println!("  fft:        {}", fft_len);
        println!("  debug:      {}", args.debug);
    }
    if !schedule_mode {
        println!("  bit:        {}={} {}={} -> out={}", a1_name, bit1, a2_name, bit2, bit_out);
        println!("  bit-code:   {}=({}) {}=({})", a1_name, format_bit_codes(bit1), a2_name, format_bit_codes(bit2));
        println!("  level:      {}={:?}", a1_name, levels1);
        println!("  level:      {}={:?}", a2_name, levels2);
        println!("  level-map:  {}={}", a1_name, format_level_map(bit1, &levels1));
        println!("  level-map:  {}={}", a2_name, format_level_map(bit2, &levels2));
        println!(
            "  cor-normalization: inv = 1 / (0.5 * P * nf * fft^2), P11={:.6}, P22={:.6}, P12={:.6}",
            level_power1,
            level_power2,
            (level_power1 * level_power2).sqrt()
        );
        println!("  shuffle-in: {}={:?}", a1_name, sh1_ext);
        println!("  shuffle-in: {}={:?}", a2_name, sh2_ext);
        println!("  sideband:   {}={} {}={}", a1_name, sb1_s, a2_name, sb2_s);
        println!("  sideband-normalize: {}={} {}={} (internal USB domain)", a1_name, if lsb1 { "LSB->USB" } else { "USB" }, a2_name, if lsb2 { "LSB->USB" } else { "USB" });
        println!("  sideband-output: {} (YAMAGU66 raw)", if output_lsb { "LSB" } else { "USB" });
        println!("  obs-band:   {:.3} .. {:.3} MHz", obs_mhz, obs_mhz + bw);
    }
    if !schedule_mode {
        println!("  band-align: {band_align_desc}");
        println!("  band-overlap: {} bins ({:.3} MHz, {:.3} MHz/bin)", overlap_bins, overlap_mhz, bin_mhz);
        println!("  rotation-shift: {rotation_shift_desc}");
        println!("  rotation-fringestop: delta_hz={:.6} shift_hz={:.6} residual_hz={:.6}", rotation_delta_hz, rotation_shift_hz, rotation_residual_hz);
        println!("  output-grid: {}", match out_grid { OutputGrid::Ant1 => a1_name, OutputGrid::Ant2 => a2_name });
        println!("  ant-fixed:  {}={:?}, {}={:?}", a1_name, ant1_ecef, a2_name, ant2_ecef);
    }
    if !args.compact_logs {
        println!("  cpu:        {} (compute threads: {})", cpu_threads, rayon::current_num_threads());
        println!(
            "  io:         chunk={} frames (pair-bytes={}), pipeline={} chunks",
            io_chunk_frames, bytes_per_frame_pair, io_pipeline_depth
        );
        println!(
            "  skip:       {:.6}s (xml {:.6}s + cli {:.6}s, {} samples)",
            total_skip_sec, xml_skip_sec, cli_skip_sec, total_skip_samples
        );
        if let Some(window_sec) = process_window_sec {
            println!(
                "  process-window: {:.6}s from epoch (post-skip request {:.6}s)",
                window_sec, requested_sec
            );
        }
    }
    let processed_sec = total_f as f64 * fft_len as f64 / fs;
    if !args.compact_logs {
        println!("  length:     {:.6}s requested, {:.6}s processed", requested_sec, processed_sec);
        if is_phased_mode {
            println!("  tsys:       {} ({}), {} ({})", tsys1, a1_name, tsys2, a2_name);
            println!("  eta:        {} ({}), {} ({})", eta1, a1_name, eta2, a2_name);
            println!("  gain:       {} ({}), {} ({})", gain1, a1_name, gain2, a2_name);
            println!("  weight:     {:.6} ({}), {:.6} ({})", w1, a1_name, w2, a2_name);
            println!("  results:    {}", o_dir.display());
            println!("  output:     {}", o_path.display());
        } else {
            println!("  cor-dir:    {}", o_dir.display());
        }
        println!("--------------------------------------------------");
    } else if !is_phased_mode {
        println!("[info] cor-dir: {}", o_dir.display());
    }

    let fr1 = File::open(&a1p)?;
    let fr2 = File::open(&a2p)?;
    advise_sequential_readahead(&fr1);
    advise_sequential_readahead(&fr2);
    let (mut r1, mut r2) = (BufReader::new(fr1), BufReader::new(fr2));
    let (mut s1_s, mut s2_s) = (0, 0); apply_integer_delay_with_forward_seek((net_d0 * fs).round() as i64, bit1, bit2, "init", &mut r1, &mut r2, f1_m.len(), f2_m.len(), &mut s1_s, &mut s2_s)?;
    if total_skip_samples > 0 {
        seek_forward_samples(&mut r1, total_skip_samples, bit1, f1_m.len(), &mut s1_s, "ant1", "skip")?;
        seek_forward_samples(&mut r2, total_skip_samples, bit2, f2_m.len(), &mut s2_s, "ant2", "skip")?;
    }
    let d_seek = (s2_s - s1_s) as f64 / fs; let helper = Arc::new(FftHelper::new(fft_len));
    let frame_dt = fft_len as f64 / fs;
    // Keep integer/fractional shifts on baseline-relative delay only.
    let net_d1_base = net_d_rel_no_clock0 + clock_delay_s - d_seek;
    let net_d2_base = 0.0;
    let rate_rel_no_clock_base = (geometric_rate_hz + rotation_fringe_hz + rate_user_hz) / (obs_mhz * 1e6);
    let total_rate1_base = rate_rel_no_clock_base + clock_rate_sps;
    let total_rate2_base = 0.0;
    let total_accel_base = correction_sign * ga0;
    let total_accel1_base = total_accel_base;
    let total_accel2_base = 0.0;
    let geom_delay_table_1s: Option<Arc<[f64]>> = gdi.as_ref().map(|v| {
        let total_duration_sec = total_f as f64 * frame_dt;
        // Include both interval endpoints [0, total_duration] so 1-second interpolation
        // remains valid through the last second of the process window.
        let n_update_secs = total_duration_sec.ceil().max(1.0) as usize;
        let n_points = n_update_secs + 1;
        let mut table = Vec::with_capacity(n_points);
        for sec in 0..n_points {
            let t_abs_sec = total_skip_sec + sec as f64;
            let mjd_t = v.mjd + t_abs_sec / 86400.0;
            let (_, _, gd_t, _, _) =
                geom::calculate_geometric_delay_and_derivatives(ant1_ecef, ant2_ecef, v.ra, v.dec, mjd_t);
            table.push(gd_t);
        }
        Arc::from(table.into_boxed_slice())
    });
    let extra_delay_rate_sps = (rotation_fringe_hz + rate_user_hz) / (obs_mhz * 1e6);
    if geom_delay_table_1s.is_some() {
        println!(
            "[info] Delay model refinement: 1-second geometric-delay table enabled (process scope, {} entries)",
            geom_delay_table_1s.as_ref().map(|t| t.len()).unwrap_or(0)
        );
    }
    let delay_cfg = DelayEvalConfig {
        frame_dt,
        fs,
        time_offset_s: total_skip_sec,
        geom_delay_table_1s,
        coarse_delay_s,
        delay_user_samples,
        extra_delay_rate_sps,
        clock1_delay_s,
        clock1_rate_sps,
        clock2_delay_s,
        clock2_rate_sps,
        d_seek,
        net_d1_base,
        total_rate1_base,
        total_accel1_base,
        net_d2_base,
        total_rate2_base,
        total_accel2_base,
        lo1_hz,
        lo2_hz,
    };
    println!("[info] Delay model cache: per-sector streaming evaluation");
    let source_name = selected_process
        .as_ref()
        .and_then(|p| p.object.clone())
        .or_else(|| if_d.as_ref().and_then(|d| d.source.clone()))
        .unwrap_or_else(|| "UNKNOWN".into());
    let cor_label_raw = if matches!(run_mode, RunMode::PhasedArray) {
        "phasedarray".to_string()
    } else {
        if_d
            .as_ref()
            .and_then(|d| d.stream_label.clone())
            .unwrap_or_else(|| "phasedarray".to_string())
    };
    let cor_label = sanitize_file_token(&cor_label_raw);
    let ant1_name_file = sanitize_file_token(a1_name);
    let ant2_name_file = sanitize_file_token(a2_name);

    if do_xcf {
        let mut ac = vec![Complex::new(0.0_f64, 0.0_f64); fft_len/2+1];
        println!(
            "[info] XCF pipeline: chunk={} frames, depth={} chunks",
            io_chunk_frames, io_pipeline_depth
        );
        let (tx, rx) = mpsc::sync_channel::<Result<Vec<Vec<u8>>, String>>(io_pipeline_depth);
        let (r1_p, r2_p, s1_c, s2_c) = (a1p.clone(), a2p.clone(), s1_s, s2_s);
        thread::spawn(move || {
            let mut rd1 = match File::open(&r1_p) {
                Ok(f) => {
                    advise_sequential_readahead(&f);
                    BufReader::new(f)
                }
                Err(e) => {
                    let _ = tx.send(Err(format!("failed to open {}: {e}", r1_p.display())));
                    return;
                }
            };
            let mut rd2 = match File::open(&r2_p) {
                Ok(f) => {
                    advise_sequential_readahead(&f);
                    BufReader::new(f)
                }
                Err(e) => {
                    let _ = tx.send(Err(format!("failed to open {}: {e}", r2_p.display())));
                    return;
                }
            };
            if let Err(e) = rd1.seek(SeekFrom::Start(s1_c as u64 * bit1 as u64 / 8)) {
                let _ = tx.send(Err(format!("failed to seek ant1 input: {e}")));
                return;
            }
            if let Err(e) = rd2.seek(SeekFrom::Start(s2_c as u64 * bit2 as u64 / 8)) {
                let _ = tx.send(Err(format!("failed to seek ant2 input: {e}")));
                return;
            }
            let mut nf = 0;
            while nf < total_f {
                let n = (total_f - nf).min(io_chunk_frames);
                let (mut b1, mut b2) = (vec![0u8; n * bpf1], vec![0u8; n * bpf2]);
                if let Err(e) = read_with_padding(&mut rd1, &mut b1) {
                    let _ = tx.send(Err(format!("failed reading ant1 input: {e}")));
                    return;
                }
                if let Err(e) = read_with_padding(&mut rd2, &mut b2) {
                    let _ = tx.send(Err(format!("failed reading ant2 input: {e}")));
                    return;
                }
                if tx.send(Ok(vec![b1, b2])).is_err() {
                    return;
                }
                nf += n;
            }
        });
        let (dp1, dp2) = (
            Arc::new(build_decode_plan(bit1, sh1.as_ref(), levels1.as_ref())?),
            Arc::new(build_decode_plan(bit2, sh2.as_ref(), levels2.as_ref())?),
        );
        let mut processed = 0;
        let mut dropped_xcf_frames = 0usize;
        let inv_fft2 = 1.0 / (fft_len as f64).powi(2);
        for bufs_res in rx {
            let bufs = bufs_res.map_err(|e| std::io::Error::other(format!("xcf reader error: {e}")))?;
            let (raw1, raw2) = (&bufs[0], &bufs[1]); let nf = raw1.len() / bpf1;
            let frame_delays: Vec<FrameDelayEntry> = (0..nf)
                .map(|i| compute_frame_delay_entry(processed + i, &delay_cfg))
                .collect();
            let frame_errors = AtomicUsize::new(0);
            let bc = raw1.par_chunks(bpf1).zip(raw2.par_chunks(bpf2)).enumerate().fold(|| vec![Complex::new(0.0_f64, 0.0_f64); fft_len/2+1], |mut c, (i, (r1, r2))| {
                let d = frame_delays[i];
                let (mut f1, mut f2, mut s1, mut s2) = (vec![0.0_f32; fft_len], vec![0.0_f32; fft_len], vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1], vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1]);
                if decode_block_into_with_plan(r1, fft_len, &dp1, &mut f1, lsb1).is_err() {
                    frame_errors.fetch_add(1, Ordering::Relaxed);
                    return c;
                }
                if decode_block_into_with_plan(r2, fft_len, &dp2, &mut f2, lsb2).is_err() {
                    frame_errors.fetch_add(1, Ordering::Relaxed);
                    return c;
                }
                let int_shift1 = d.int1;
                let int_shift2 = d.int2;
                let frac_delay1 = d.frac1;
                let frac_delay2 = d.frac2;
                apply_integer_sample_shift_zerofill(&mut f1, int_shift1);
                apply_integer_sample_shift_zerofill(&mut f2, int_shift2);
                if helper.forward_r2c_process(&mut f1, &mut s1).is_err() {
                    frame_errors.fetch_add(1, Ordering::Relaxed);
                    return c;
                }
                if helper.forward_r2c_process(&mut f2, &mut s2).is_err() {
                    frame_errors.fetch_add(1, Ordering::Relaxed);
                    return c;
                }
                
                // Use observing reference frequency in .cor header (not upper band edge)
                // so residual-rate definition matches frinZ delay/rate search.
                let fr_mix = d.fr_lo1 * d.fr_lo2.conj();
                apply_delay_and_rate_regular_bins(&mut s1, fft_len, fs/fft_len as f64, frac_delay1, 0.0, 0.0, 0.0, false);
                apply_delay_and_rate_regular_bins(&mut s2, fft_len, fs/fft_len as f64, frac_delay2, 0.0, 0.0, 0.0, false);
                
                for k in 0..(ba.a1e - ba.a1s) {
                    let i1 = ba.a1s + k;
                    let i2 = ba.a2s + k;
                    let v = (s1[i1] * s2[i2].conj()) * fr_mix;
                    match out_grid {
                        OutputGrid::Ant1 => c[i1] += Complex::new(v.re as f64, v.im as f64),
                        OutputGrid::Ant2 => c[i2] += Complex::new(v.re as f64, v.im as f64),
                    }
                }
                c
            }).reduce(|| vec![Complex::new(0.0_f64, 0.0_f64); fft_len/2+1], |mut c1, c2| { for k in 0..c1.len() { c1[k] += c2[k]; } c1 });
            dropped_xcf_frames += frame_errors.load(Ordering::Relaxed);
            for k in 0..ac.len() { ac[k] += bc[k] * inv_fft2; }
            processed += nf;
            print!("\rCorrelating ({}/{})", processed, total_f);
            let _ = std::io::stdout().flush();
        }
        println!();
        if dropped_xcf_frames > 0 {
            println!(
                "[warn] XCF skipped {} frame(s) due to decode/FFT errors",
                dropped_xcf_frames
            );
        }
        let xcf_res = finalize_cross_spectrum(
            &mut ac,
            helper.as_ref(),
            fft_len,
            fft_len / 2 + 1,
            &(-(fft_len as i32 / 2)..fft_len as i32 / 2).collect::<Vec<_>>(),
            fs,
            obs_mhz,
            true,
            &o_dir,
            false,
            false,
        )?;
        if let Some(delay_s) = xcf_res.delay_seconds_from_phase {
            println!("[info] Delay estimate from phase slope: {:.9e} s", delay_s);
        }
    }

    if total_f > 0 && do_synth {
        if is_phased_mode {
            println!("[info] Synthesising output...");
        } else {
            println!("[info] Integrating correlation sectors...");
        }
        let mut wr: Option<BufWriter<File>> = if write_raw {
            Some(BufWriter::new(File::create(&o_path)?))
        } else {
            None
        };
        let mut dbg_writer: Option<BufWriter<File>> = if args.debug {
            let dbg_path = o_dir.join(format!("YAMAGU66_{}_debug.log", c_tag));
            println!("[info] Debug log: {}", dbg_path.display());
            let mut w = BufWriter::new(File::create(&dbg_path)?);
            writeln!(w, "# phased_array debug log")?;
            writeln!(w, "# schedule_mode={schedule_mode} fft={fft_len} fs_hz={fs:.3}")?;
            writeln!(w, "# all frames are logged")?;
            Some(w)
        } else {
            None
        };
        let frame_sec = fft_len as f64 / fs;
        let total_duration_sec = total_f as f64 * frame_sec;
        // "second sectors" should track actual covered duration without
        // creating an extra tail sector from tiny floating-point excess.
        let mut sector_count = total_duration_sec.round().max(1.0) as usize;
        sector_count = sector_count.min(total_f.max(1));
        let base = total_f / sector_count;
        let extra = total_f % sector_count;
        let mut sec_counts = Vec::with_capacity(sector_count);
        for si in 0..sector_count {
            let nf = base + if si < extra { 1 } else { 0 };
            if nf > 0 {
                sec_counts.push(nf);
            }
        }
        let dp1 = Arc::new(build_decode_plan(bit1, sh1.as_ref(), levels1.as_ref())?);
        let dp2 = Arc::new(build_decode_plan(bit2, sh2.as_ref(), levels2.as_ref())?);
        let mut emitted = 0;
        let cor_h_freq_hz = obs_mhz * 1e6;
        let mut cw_ph = if write_phased_cor {
            Some(CorWriter::create(&o_dir.join(format!("YAMAGU66_YAMAGU66_{}_{}.cor", c_tag, cor_label)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: "YAMAGU66", code: b'M', ecef_m: ant1_ecef }, CorStation { name: "YAMAGU66", code: b'M', ecef_m: ant1_ecef })?)
        } else {
            None
        };
        let mut cw_11 = if write_acf_cor {
            Some(CorWriter::create(&o_dir.join(format!("{}_{}_{}_{}.cor", ant1_name_file, ant1_name_file, c_tag, cor_label)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: a1_name, code: a1_code, ecef_m: ant1_ecef }, CorStation { name: a1_name, code: a1_code, ecef_m: ant1_ecef })?)
        } else {
            None
        };
        let mut cw_12 = if write_xcf_cor {
            Some(CorWriter::create(&o_dir.join(format!("{}_{}_{}_{}.cor", ant1_name_file, ant2_name_file, c_tag, cor_label)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: a1_name, code: a1_code, ecef_m: ant1_ecef }, CorStation { name: a2_name, code: a2_code, ecef_m: ant2_ecef })?)
        } else {
            None
        };
        let mut cw_22 = if write_acf_cor {
            Some(CorWriter::create(&o_dir.join(format!("{}_{}_{}_{}.cor", ant2_name_file, ant2_name_file, c_tag, cor_label)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: a2_name, code: a2_code, ecef_m: ant2_ecef }, CorStation { name: a2_name, code: a2_code, ecef_m: ant2_ecef })?)
        } else {
            None
        };
        
        let prefetch_depth = io_pipeline_depth.max(2);
        println!(
            "[info] I/O prefetch: process-window reader enabled (pipeline={} chunks)",
            prefetch_depth
        );
        let sec_counts_for_read = sec_counts.clone();
        let (tx_sec, rx_sec) =
            mpsc::sync_channel::<Result<(Vec<u8>, Vec<u8>), String>>(prefetch_depth);
        let (a1_read, a2_read) = (a1p.clone(), a2p.clone());
        let reader_handle = thread::spawn(move || {
            let fpr1 = match File::open(&a1_read) {
                Ok(f) => f,
                Err(e) => {
                    let _ = tx_sec.send(Err(format!(
                        "failed to open {}: {e}",
                        a1_read.display()
                    )));
                    return;
                }
            };
            let fpr2 = match File::open(&a2_read) {
                Ok(f) => f,
                Err(e) => {
                    let _ = tx_sec.send(Err(format!(
                        "failed to open {}: {e}",
                        a2_read.display()
                    )));
                    return;
                }
            };
            advise_sequential_readahead(&fpr1);
            advise_sequential_readahead(&fpr2);
            let (mut pr1, mut pr2) = (BufReader::new(fpr1), BufReader::new(fpr2));
            if let Err(e) = pr1.seek(SeekFrom::Start(s1_s as u64 * bit1 as u64 / 8)) {
                let _ = tx_sec.send(Err(format!("failed to seek ant1 input: {e}")));
                return;
            }
            if let Err(e) = pr2.seek(SeekFrom::Start(s2_s as u64 * bit2 as u64 / 8)) {
                let _ = tx_sec.send(Err(format!("failed to seek ant2 input: {e}")));
                return;
            }
            for nf in sec_counts_for_read {
                let mut b1 = vec![0u8; nf * bpf1];
                let mut b2 = vec![0u8; nf * bpf2];
                if let Err(e) = read_with_padding(&mut pr1, &mut b1) {
                    let _ = tx_sec.send(Err(format!("failed reading ant1 input: {e}")));
                    return;
                }
                if let Err(e) = read_with_padding(&mut pr2, &mut b2) {
                    let _ = tx_sec.send(Err(format!("failed reading ant2 input: {e}")));
                    return;
                }
                if tx_sec.send(Ok((b1, b2))).is_err() {
                    return;
                }
            }
        });
        let need_phased_products = write_raw || write_phased_cor || plot_phased;
        let need_acf_products = write_acf_cor;
        let need_xcf_products = write_xcf_cor;
        let acf_overlap_only = matches!(run_mode, RunMode::Acf)
            && need_acf_products
            && !need_xcf_products
            && !need_phased_products;
        let mut acc_ph_total = if need_phased_products {
            Some(vec![0.0; fft_len / 2 + 1])
        } else {
            None
        };
        let mut acc_11_total = vec![0.0; fft_len/2+1];
        let mut acc_22_total = vec![0.0; fft_len/2+1];
        for (si, &nf) in sec_counts.iter().enumerate() {
            let frame_delays: Vec<FrameDelayEntry> = (0..nf)
                .map(|i| compute_frame_delay_entry(emitted + i, &delay_cfg))
                .collect();
            if let Some(dw) = dbg_writer.as_mut() {
                writeln!(dw, "[sec {}] frames={} start_frame={}", si + 1, nf, emitted)?;
                for i in 0..nf {
                    let frame_idx = emitted + i;
                    let d = frame_delays[i];
                    let t_mid = (frame_idx as f64 + 0.5) * frame_dt;
                    let int_shift1 = d.int1;
                    let int_shift2 = d.int2;
                    let frac_delay1 = d.frac1;
                    let frac_delay2 = d.frac2;
                    let tau1 = int_shift1 as f64 / fs + frac_delay1;
                    let tau2 = int_shift2 as f64 / fs + frac_delay2;
                    writeln!(
                        dw,
                        "  frame={} t_mid={:.9} tau1={:.9e} tau2={:.9e} int1={} int2={} frac1={:.9e} frac2={:.9e}",
                        emitted + i, t_mid, tau1, tau2, int_shift1, int_shift2, frac_delay1, frac_delay2
                    )?;
                }
            }
            let (raw1_vec, raw2_vec) = match rx_sec.recv() {
                Ok(Ok(v)) => v,
                Ok(Err(e)) => {
                    return Err(std::io::Error::other(format!(
                        "synth reader error: {e}"
                    ))
                    .into())
                }
                Err(e) => {
                    return Err(std::io::Error::other(format!(
                        "synth reader channel error: {e}"
                    ))
                    .into())
                }
            };
            let raw1: &[u8] = &raw1_vec;
            let raw2: &[u8] = &raw2_vec;
            let sector_failures = AtomicUsize::new(0);
            let process_frame = |i: usize, out_f: Option<&mut [u8]>| -> Option<(Vec<f64>, Vec<f64>, Vec<Complex<f64>>, Vec<f64>)> {
                let d = frame_delays[i];
                let (mut f1, mut f2, mut s1, mut s2) = (vec![0.0_f32; fft_len], vec![0.0_f32; fft_len], vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1], vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1]);
                if decode_block_into_with_plan(&raw1[i*bpf1..(i+1)*bpf1], fft_len, &dp1, &mut f1, lsb1).is_err() {
                    if let Some(out_f) = out_f { out_f.fill(0); }
                    sector_failures.fetch_add(1, Ordering::Relaxed);
                    return None;
                }
                if decode_block_into_with_plan(&raw2[i*bpf2..(i+1)*bpf2], fft_len, &dp2, &mut f2, lsb2).is_err() {
                    if let Some(out_f) = out_f { out_f.fill(0); }
                    sector_failures.fetch_add(1, Ordering::Relaxed);
                    return None;
                }
                let int_shift1 = d.int1;
                let int_shift2 = d.int2;
                let frac_delay1 = d.frac1;
                let frac_delay2 = d.frac2;
                apply_integer_sample_shift_zerofill(&mut f1, int_shift1);
                apply_integer_sample_shift_zerofill(&mut f2, int_shift2);
                if helper.forward_r2c_process(&mut f1, &mut s1).is_err() {
                    if let Some(out_f) = out_f { out_f.fill(0); }
                    sector_failures.fetch_add(1, Ordering::Relaxed);
                    return None;
                }
                if helper.forward_r2c_process(&mut f2, &mut s2).is_err() {
                    if let Some(out_f) = out_f { out_f.fill(0); }
                    sector_failures.fetch_add(1, Ordering::Relaxed);
                    return None;
                }
                // Keep the same reference as cross-correlation path.
                let fr_lo1 = d.fr_lo1;
                let fr_lo2 = d.fr_lo2;
                apply_delay_and_rate_regular_bins(&mut s1, fft_len, fs/fft_len as f64, frac_delay1, 0.0, 0.0, 0.0, false);
                apply_delay_and_rate_regular_bins(&mut s2, fft_len, fs/fft_len as f64, frac_delay2, 0.0, 0.0, 0.0, false);

                let mut phased_pow = vec![0.0; fft_len/2+1];
                let mut p11 = vec![0.0; fft_len/2+1];
                let mut p12 = vec![Complex::new(0.0_f64, 0.0_f64); fft_len/2+1];
                let mut p22 = vec![0.0; fft_len/2+1];

                let mut s1_aligned = vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1];
                let mut s2_aligned = vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1];
                match out_grid {
                    OutputGrid::Ant1 => {
                        if need_xcf_products || need_acf_products || need_phased_products {
                            for k in 0..(ba.a1e - ba.a1s) {
                                let i1 = ba.a1s + k;
                                let i2 = ba.a2s + k;
                                s2_aligned[i1] = s2[i2] * fr_lo2;
                            }
                        }
                        let s1c: Vec<Complex<f32>> = if need_xcf_products || need_acf_products || need_phased_products {
                            s1.iter().map(|z| *z * fr_lo1).collect::<Vec<_>>()
                        } else {
                            Vec::new()
                        };
                        if need_phased_products {
                            let mut cb = vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1];
                            for k in 0..cb.len() { cb[k] = s1c[k] * (w1 as f32) + s2_aligned[k] * (w2 as f32); }
                            phased_pow = cb.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                            if let Some(out_f) = out_f {
                                cb[0].im = 0.0_f32;
                                if fft_len % 2 == 0 { cb[fft_len/2].im = 0.0_f32; }
                                let mut out_t = vec![0.0_f32; fft_len];
                                if helper.inverse_c2r_process(&mut cb, &mut out_t).is_err() {
                                    out_f.fill(0);
                                    sector_failures.fetch_add(1, Ordering::Relaxed);
                                    return None;
                                }
                                if output_lsb {
                                    for odd in out_t.iter_mut().skip(1).step_by(2) { *odd = -*odd; }
                                }
                                let mut tmp_enc = Vec::new();
                                if quantise_frame(&out_t, bit_out, &levels1, sh1.as_ref(), &mut tmp_enc).is_err() {
                                    out_f.fill(0);
                                    sector_failures.fetch_add(1, Ordering::Relaxed);
                                    return None;
                                }
                                out_f.copy_from_slice(&tmp_enc);
                            }
                        }
                        if need_acf_products {
                            p11 = s1c.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                            p22 = s2_aligned.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                        }
                        if need_xcf_products {
                            p12 = s1c.iter().zip(s2_aligned.iter()).map(|(z1, z2)| {
                                let v = *z1 * z2.conj();
                                Complex::new(v.re as f64, v.im as f64)
                            }).collect::<Vec<_>>();
                        }
                    }
                    OutputGrid::Ant2 => {
                        if need_xcf_products || need_acf_products || need_phased_products {
                            for k in 0..(ba.a1e - ba.a1s) {
                                let i1 = ba.a1s + k;
                                let i2 = ba.a2s + k;
                                s1_aligned[i2] = s1[i1] * fr_lo1;
                            }
                        }
                        let s2c: Vec<Complex<f32>> = if need_xcf_products || need_acf_products || need_phased_products {
                            s2.iter().map(|z| *z * fr_lo2).collect::<Vec<_>>()
                        } else {
                            Vec::new()
                        };
                        if need_phased_products {
                            let mut cb = vec![Complex::new(0.0_f32, 0.0_f32); fft_len/2+1];
                            for k in 0..cb.len() { cb[k] = s1_aligned[k] * (w1 as f32) + s2c[k] * (w2 as f32); }
                            phased_pow = cb.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                            if let Some(out_f) = out_f {
                                cb[0].im = 0.0_f32;
                                if fft_len % 2 == 0 { cb[fft_len/2].im = 0.0_f32; }
                                let mut out_t = vec![0.0_f32; fft_len];
                                if helper.inverse_c2r_process(&mut cb, &mut out_t).is_err() {
                                    out_f.fill(0);
                                    sector_failures.fetch_add(1, Ordering::Relaxed);
                                    return None;
                                }
                                if output_lsb {
                                    for odd in out_t.iter_mut().skip(1).step_by(2) { *odd = -*odd; }
                                }
                                let mut tmp_enc = Vec::new();
                                if quantise_frame(&out_t, bit_out, &levels1, sh1.as_ref(), &mut tmp_enc).is_err() {
                                    out_f.fill(0);
                                    sector_failures.fetch_add(1, Ordering::Relaxed);
                                    return None;
                                }
                                out_f.copy_from_slice(&tmp_enc);
                            }
                        }
                        if need_acf_products {
                            p11 = s1_aligned.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                            p22 = s2c.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                        }
                        if need_xcf_products {
                            p12 = s1_aligned.iter().zip(s2c.iter()).map(|(z1, z2)| {
                                let v = *z1 * z2.conj();
                                Complex::new(v.re as f64, v.im as f64)
                            }).collect::<Vec<_>>();
                        }
                    }
                }
                Some((phased_pow, p11, p12, p22))
            };
            let zero_acc = || {
                (
                    if need_phased_products {
                        vec![0.0; fft_len / 2 + 1]
                    } else {
                        Vec::new()
                    },
                    vec![0.0; fft_len / 2 + 1],
                    vec![Complex::new(0.0, 0.0); fft_len / 2 + 1],
                    vec![0.0; fft_len / 2 + 1],
                )
            };
            let reduce_acc = |mut acc1: (Vec<f64>, Vec<f64>, Vec<Complex<f64>>, Vec<f64>), acc2: (Vec<f64>, Vec<f64>, Vec<Complex<f64>>, Vec<f64>)| {
                if need_phased_products {
                    for k in 0..acc1.0.len() {
                        acc1.0[k] += acc2.0[k];
                    }
                }
                for k in 0..acc1.1.len() {
                    acc1.1[k] += acc2.1[k];
                    acc1.2[k] += acc2.2[k];
                    acc1.3[k] += acc2.3[k];
                }
                acc1
            };
            let (batch_ph, batch_11, batch_12, batch_22) = if write_raw {
                let mut enc = vec![0u8; nf * bpf_o];
                let acc = enc
                    .par_chunks_mut(bpf_o)
                    .enumerate()
                    .map(|(i, out_f)| process_frame(i, Some(out_f)))
                    .fold(zero_acc, |mut acc, res| {
                        if let Some((p_ph, p_11, p_12, p_22)) = res {
                            for k in 0..acc.0.len() {
                                acc.0[k] += p_ph[k];
                                acc.1[k] += p_11[k];
                                acc.2[k] += p_12[k];
                                acc.3[k] += p_22[k];
                            }
                        }
                        acc
                    })
                    .reduce(zero_acc, reduce_acc);
                if let Some(w) = wr.as_mut() {
                    w.write_all(&enc)?;
                }
                acc
            } else {
                struct ThreadAccum {
                    acc_ph: Vec<f64>,
                    acc_11: Vec<f64>,
                    acc_12: Vec<Complex<f64>>,
                    acc_22: Vec<f64>,
                    f1: Vec<f32>,
                    f2: Vec<f32>,
                    s1: Vec<Complex<f32>>,
                    s2: Vec<Complex<f32>>,
                    fft_scratch: FftScratch,
                }
                let half = fft_len / 2 + 1;
                let init = || ThreadAccum {
                    acc_ph: if need_phased_products {
                        vec![0.0; half]
                    } else {
                        Vec::new()
                    },
                    acc_11: vec![0.0; half],
                    acc_12: vec![Complex::new(0.0_f64, 0.0_f64); half],
                    acc_22: vec![0.0; half],
                    f1: vec![0.0_f32; fft_len],
                    f2: vec![0.0_f32; fft_len],
                    s1: vec![Complex::new(0.0_f32, 0.0_f32); half],
                    s2: vec![Complex::new(0.0_f32, 0.0_f32); half],
                    fft_scratch: helper.make_scratch(),
                };
                let frames_per_job = (nf / (cpu_threads.saturating_mul(8)).max(1)).clamp(128, 2048);
                let chunk_starts: Vec<usize> = (0..nf).step_by(frames_per_job).collect();
                let mut out = chunk_starts
                    .into_par_iter()
                    .map(|start| {
                        let mut st = init();
                        let end = (start + frames_per_job).min(nf);
                        for i in start..end {
                        let d = frame_delays[i];
                        if decode_block_into_with_plan(&raw1[i*bpf1..(i+1)*bpf1], fft_len, &dp1, &mut st.f1, lsb1).is_err() {
                            sector_failures.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        if decode_block_into_with_plan(&raw2[i*bpf2..(i+1)*bpf2], fft_len, &dp2, &mut st.f2, lsb2).is_err() {
                            sector_failures.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        let int_shift1 = d.int1;
                        let int_shift2 = d.int2;
                        let frac_delay1 = d.frac1;
                        let frac_delay2 = d.frac2;
                        apply_integer_sample_shift_zerofill(&mut st.f1, int_shift1);
                        apply_integer_sample_shift_zerofill(&mut st.f2, int_shift2);
                        if helper
                            .forward_r2c_process_with_scratch(
                                &mut st.f1,
                                &mut st.s1,
                                &mut st.fft_scratch,
                            )
                            .is_err()
                        {
                            sector_failures.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        if helper
                            .forward_r2c_process_with_scratch(
                                &mut st.f2,
                                &mut st.s2,
                                &mut st.fft_scratch,
                            )
                            .is_err()
                        {
                            sector_failures.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        let fr_lo1 = d.fr_lo1;
                        let fr_lo2 = d.fr_lo2;
                        let fr_mix = fr_lo1 * fr_lo2.conj();
                        let overlap_len = ba.a1e - ba.a1s;

                        match out_grid {
                            OutputGrid::Ant1 => {
                                if need_xcf_products {
                                    apply_delay_and_rate_regular_bins_range(
                                        &mut st.s1,
                                        fft_len,
                                        df_hz,
                                        frac_delay1,
                                        0.0,
                                        0.0,
                                        0.0,
                                        false,
                                        ba.a1s,
                                        ba.a1e,
                                    );
                                    apply_delay_and_rate_regular_bins_range(
                                        &mut st.s2,
                                        fft_len,
                                        df_hz,
                                        frac_delay2,
                                        0.0,
                                        0.0,
                                        0.0,
                                        false,
                                        ba.a2s,
                                        ba.a2s + overlap_len,
                                    );
                                }
                                if need_acf_products {
                                    if acf_overlap_only {
                                        for k in 0..overlap_len {
                                            let i1 = ba.a1s + k;
                                            let i2 = ba.a2s + k;
                                            st.acc_11[i1] += st.s1[i1].norm_sqr() as f64;
                                            st.acc_22[i1] += st.s2[i2].norm_sqr() as f64;
                                        }
                                    } else {
                                        for k in 0..half {
                                            st.acc_11[k] += st.s1[k].norm_sqr() as f64;
                                        }
                                        if need_xcf_products {
                                            for k in 0..overlap_len {
                                                let i1 = ba.a1s + k;
                                                let i2 = ba.a2s + k;
                                                st.acc_22[i1] += st.s2[i2].norm_sqr() as f64;
                                                let v = (st.s1[i1] * st.s2[i2].conj()) * fr_mix;
                                                st.acc_12[i1] += Complex::new(v.re as f64, v.im as f64);
                                            }
                                        } else {
                                            for k in 0..overlap_len {
                                                let i1 = ba.a1s + k;
                                                let i2 = ba.a2s + k;
                                                st.acc_22[i1] += st.s2[i2].norm_sqr() as f64;
                                            }
                                        }
                                    }
                                }
                                if need_xcf_products {
                                    if !(need_acf_products && !acf_overlap_only) {
                                        for k in 0..overlap_len {
                                            let i1 = ba.a1s + k;
                                            let i2 = ba.a2s + k;
                                            let v = (st.s1[i1] * st.s2[i2].conj()) * fr_mix;
                                            st.acc_12[i1] += Complex::new(v.re as f64, v.im as f64);
                                        }
                                    }
                                }
                            }
                            OutputGrid::Ant2 => {
                                if need_xcf_products {
                                    apply_delay_and_rate_regular_bins_range(
                                        &mut st.s2,
                                        fft_len,
                                        df_hz,
                                        frac_delay2,
                                        0.0,
                                        0.0,
                                        0.0,
                                        false,
                                        ba.a2s,
                                        ba.a2s + (ba.a1e - ba.a1s),
                                    );
                                    apply_delay_and_rate_regular_bins_range(
                                        &mut st.s1,
                                        fft_len,
                                        df_hz,
                                        frac_delay1,
                                        0.0,
                                        0.0,
                                        0.0,
                                        false,
                                        ba.a1s,
                                        ba.a1e,
                                    );
                                }
                                if need_acf_products {
                                    if acf_overlap_only {
                                        for k in 0..overlap_len {
                                            let i1 = ba.a1s + k;
                                            let i2 = ba.a2s + k;
                                            st.acc_11[i2] += st.s1[i1].norm_sqr() as f64;
                                            st.acc_22[i2] += st.s2[i2].norm_sqr() as f64;
                                        }
                                    } else {
                                        for k in 0..half {
                                            st.acc_22[k] += st.s2[k].norm_sqr() as f64;
                                        }
                                        if need_xcf_products {
                                            for k in 0..overlap_len {
                                                let i1 = ba.a1s + k;
                                                let i2 = ba.a2s + k;
                                                st.acc_11[i2] += st.s1[i1].norm_sqr() as f64;
                                                let v = (st.s1[i1] * st.s2[i2].conj()) * fr_mix;
                                                st.acc_12[i2] += Complex::new(v.re as f64, v.im as f64);
                                            }
                                        } else {
                                            for k in 0..overlap_len {
                                                let i1 = ba.a1s + k;
                                                let i2 = ba.a2s + k;
                                                st.acc_11[i2] += st.s1[i1].norm_sqr() as f64;
                                            }
                                        }
                                    }
                                }
                                if need_xcf_products {
                                    if !(need_acf_products && !acf_overlap_only) {
                                        for k in 0..overlap_len {
                                            let i1 = ba.a1s + k;
                                            let i2 = ba.a2s + k;
                                            let v = (st.s1[i1] * st.s2[i2].conj()) * fr_mix;
                                            st.acc_12[i2] += Complex::new(v.re as f64, v.im as f64);
                                        }
                                    }
                                }
                            }
                        }
                        }
                        st
                    })
                    .reduce(init, |mut a, b| {
                        if need_phased_products {
                            for k in 0..half {
                                a.acc_ph[k] += b.acc_ph[k];
                            }
                        }
                        for k in 0..half {
                            a.acc_11[k] += b.acc_11[k];
                            a.acc_12[k] += b.acc_12[k];
                            a.acc_22[k] += b.acc_22[k];
                        }
                        a
                    });
                (
                    std::mem::take(&mut out.acc_ph),
                    std::mem::take(&mut out.acc_11),
                    std::mem::take(&mut out.acc_12),
                    std::mem::take(&mut out.acc_22),
                )
            };
            let sec_failed = sector_failures.load(Ordering::Relaxed);
            if sec_failed > 0 {
                println!(
                    "[warn] Sector {} skipped {} frame(s) due to decode/FFT/quantize errors",
                    si + 1,
                    sec_failed
                );
            }
            emitted += nf;
            if is_phased_mode {
                print!(
                    "\r[info] Synthesised second {}/{} ({:.2}%)",
                    si + 1,
                    sec_counts.len(),
                    (emitted as f64 / total_f as f64) * 100.0
                );
            } else {
                print!(
                    "\r[info] Processed second {}/{} ({:.2}%)",
                    si + 1,
                    sec_counts.len(),
                    (emitted as f64 / total_f as f64) * 100.0
                );
            }
            std::io::stdout().flush()?;
            if let Some(acc_ph) = acc_ph_total.as_mut() {
                for k in 0..acc_ph.len() {
                    acc_ph[k] += batch_ph[k];
                }
            }
            for k in 0..acc_11_total.len() {
                acc_11_total[k] += batch_11[k];
                acc_22_total[k] += batch_22[k];
            }
            let norm_base = nf as f64 * (fft_len as f64).powi(2);
            let inv_11 = 1.0 / (COR_ONE_SIDED_POWER_FACTOR * level_power1 * norm_base);
            let inv_22 = 1.0 / (COR_ONE_SIDED_POWER_FACTOR * level_power2 * norm_base);
            let inv_12 = 1.0
                / (COR_ONE_SIDED_POWER_FACTOR * (level_power1 * level_power2).sqrt() * norm_base);
            let phased_power_ref = (w1 * w1 * level_power1 + w2 * w2 * level_power2).max(1e-12);
            let inv_ph = 1.0 / (COR_ONE_SIDED_POWER_FACTOR * phased_power_ref * norm_base);
            if let Some(w) = cw_ph.as_mut() {
                let s_ph: Vec<Complex<f32>> = batch_ph
                    .iter()
                    .take(fft_len / 2)
                    .map(|&v| Complex::new((v * inv_ph) as f32, 0.0))
                    .collect();
                w.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_ph)?;
            }
            if let Some(w) = cw_11.as_mut() {
                let s_11: Vec<Complex<f32>> = batch_11
                    .iter()
                    .take(fft_len / 2)
                    .map(|&v| Complex::new((v * inv_11) as f32, 0.0))
                    .collect();
                w.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_11)?;
            }
            if let Some(w) = cw_12.as_mut() {
                let s_12: Vec<Complex<f32>> = batch_12
                    .iter()
                    .take(fft_len / 2)
                    .map(|v| Complex::new((v.re * inv_12) as f32, (v.im * inv_12) as f32))
                    .collect();
                w.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_12)?;
            }
            if let Some(w) = cw_22.as_mut() {
                let s_22: Vec<Complex<f32>> = batch_22
                    .iter()
                    .take(fft_len / 2)
                    .map(|&v| Complex::new((v * inv_22) as f32, 0.0))
                    .collect();
                w.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_22)?;
            }
        }
        drop(rx_sec);
        if reader_handle.join().is_err() {
            return Err("synth reader thread panicked".into());
        }
        if !sec_counts.is_empty() {
            println!();
        }
        if let Some(w) = wr.as_mut() {
            w.flush()?;
        }
        if let Some(w) = cw_ph.take() {
            w.finalize()?;
        }
        if let Some(w) = cw_11.take() {
            w.finalize()?;
        }
        if let Some(w) = cw_12.take() {
            w.finalize()?;
        }
        if let Some(w) = cw_22.take() {
            w.finalize()?;
        }
        if let Some(mut dw) = dbg_writer {
            dw.flush()?;
        }

        if let Some(p) = schedule_stdout_path.as_ref() {
            let mut outf = OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)?;
            writeln!(
                outf,
                "  done: emitted_frames={} seconds={} output_tag={} label={}",
                emitted,
                sec_counts.len(),
                c_tag,
                cor_label
            )?;
        }

        if plot_phased {
            println!("[info] Generating phased-array plots...");
            let power_norm = (emitted as f64 * fft_len as f64).max(1.0);
            let acc_ph = acc_ph_total
                .as_mut()
                .ok_or("internal error: phased accumulator not available")?;
            let phased_auto_mag = finalize_auto_spectrum(acc_ph, power_norm)?;
            let a11_auto_mag = finalize_auto_spectrum(&mut acc_11_total, power_norm)?;
            let a22_auto_mag = finalize_auto_spectrum(&mut acc_22_total, power_norm)?;
            // Plotting convention: use [0 .. fft/2-1] for even fft (exclude Nyquist bin).
            let spec_bins = if fft_len % 2 == 0 {
                fft_len / 2
            } else {
                phased_auto_mag.len()
            };
            let phased_plot = &phased_auto_mag[..spec_bins];
            let a11_plot = &a11_auto_mag[..spec_bins];
            let a22_plot = &a22_auto_mag[..spec_bins];

            let df = df_hz / 1e6;
            let freqs_obs_mhz: Vec<f64> = (0..spec_bins).map(|i| obs_mhz + (i as f64 * df)).collect();

            plot_series_f64_x(&freqs_obs_mhz, phased_plot, "Phased Auto-Spectrum (ObsRef)", &o_dir.join(format!("YAMAGU66_{}_phased_auto_spectrum.png", c_tag)).to_string_lossy(), "Frequency (MHz)", "Power", None, "Auto-Spectrum")?;
            
            let amp_ph: Vec<f64> = phased_plot.iter().map(|&v| v.sqrt()).collect();
            let amp_11: Vec<f64> = a11_plot.iter().map(|&v| v.sqrt()).collect();
            let amp_22: Vec<f64> = a22_plot.iter().map(|&v| v.sqrt()).collect();
            let l1 = format!("{a1_name} (Ref)");
            let l2 = format!("{a2_name} (Target)");
            plot_multi_series_f64_x(&freqs_obs_mhz, &[(&amp_ph, &BLUE, "YAMAGU66 (Phased)"), (&amp_11, &GREEN, &l1), (&amp_22, &RED, &l2)], "Phased Spectrum Amplitude (ObsRef)", &o_dir.join(format!("YAMAGU66_{}_phased_spectrum_amplitude.png", c_tag)).to_string_lossy(), "Frequency (MHz)", "Amplitude", None)?;

            let mut full_spec = vec![Complex::new(0.0_f32, 0.0_f32); fft_len];
            for (i, &v) in phased_plot.iter().enumerate() {
                let vf = v as f32;
                full_spec[i] = Complex::new(vf, 0.0_f32);
                if i > 0 && i < fft_len/2 { full_spec[fft_len - i] = Complex::new(vf, 0.0_f32); }
            }
            helper.inverse_c2c(&mut full_spec)?;
            let acf_mag: Vec<f64> = full_spec.iter().map(|c| c.re as f64).collect();
            let acf_shifted: Vec<f64> = acf_mag.iter().cycle().skip(fft_len/2).take(fft_len).copied().collect();
            let lags: Vec<i32> = (-(fft_len as i32 / 2)..(fft_len as i32 / 2)).collect();
            plot_series_with_x(&lags, &[(&acf_shifted, &BLUE)], "Phased Autocorrelation", &o_dir.join(format!("YAMAGU66_{}_phased_autocorrelation.png", c_tag)).to_string_lossy(), "Lag (samples)", "ACF", None, None)?;
        }
    }
    Ok(())
}
