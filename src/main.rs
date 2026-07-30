mod acf;
mod affinity;
mod args;
mod cor;
mod eop;
mod fringe;
mod geom;
mod ifile;
mod model_diag;
mod model_diag_output;
mod model_sweep;
mod plot;
mod pulsar;
mod utils;
mod xcf;
mod xml;

use std::collections::{HashMap, HashSet};
use std::fs::{read_to_string, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{mpsc, Arc, Mutex, OnceLock};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use clap::Parser;
use num_complex::Complex;
use rayon::prelude::*;

use acf::finalize_auto_spectrum;
use args::{
    parse_levels, parse_shuffle, resolve_per_antenna_config,
    resolve_per_antenna_config_with_defaults, resolve_weight, DEFAULT_SHUFFLE_IN,
};
use cor::{
    epoch_to_yyyydddhhmmss, unix_seconds_to_yyyydddhhmmss, CorHeaderConfig, CorStation, CorWriter,
};
use plot::{plot_multi_series_f64_x, plot_series_f64_x, plot_series_with_x, BLUE, GREEN, RED};
use pulsar::{FoldAccum, FoldProduct, PulsarRuntime};
use utils::{
    build_decode_plan, decode_block_into_with_plan, quantise_frame, DecodePlan, DynError,
    FftHelper, FftScratch,
};
use xcf::finalize_cross_spectrum;

const DEFAULT_COARSE_DELAY_S: f64 = 0.0;
// GICO3-compatible one-sided spectrum normalization factor:
// .cor stores only positive-frequency bins for real-input FFT outputs.
const COR_ONE_SIDED_POWER_FACTOR: f64 = 0.5;
static STDOUT_LOG_WRITER: OnceLock<Mutex<Option<BufWriter<File>>>> = OnceLock::new();
static STDOUT_LOG_LOCK: OnceLock<Mutex<()>> = OnceLock::new();

fn stdout_log_writer() -> &'static Mutex<Option<BufWriter<File>>> {
    STDOUT_LOG_WRITER.get_or_init(|| Mutex::new(None))
}

fn stdout_log_lock() -> &'static Mutex<()> {
    STDOUT_LOG_LOCK.get_or_init(|| Mutex::new(()))
}

fn emit_stdout_fmt(args: std::fmt::Arguments<'_>, newline: bool) {
    let text = args.to_string();
    let _guard = stdout_log_lock().lock().ok();
    {
        let mut out = std::io::stdout().lock();
        let _ = out.write_all(text.as_bytes());
        if newline {
            let _ = out.write_all(b"\n");
        }
        let _ = out.flush();
    }
    if let Ok(mut writer_guard) = stdout_log_writer().lock() {
        if let Some(writer) = writer_guard.as_mut() {
            let _ = writer.write_all(text.as_bytes());
            if newline {
                let _ = writer.write_all(b"\n");
            }
            let _ = writer.flush();
        }
    }
}

macro_rules! print {
    ($($arg:tt)*) => {{
        crate::emit_stdout_fmt(format_args!($($arg)*), false);
    }};
}

macro_rules! println {
    () => {{
        crate::emit_stdout_fmt(format_args!(""), true);
    }};
    ($($arg:tt)*) => {{
        crate::emit_stdout_fmt(format_args!($($arg)*), true);
    }};
}

fn init_stdout_log_for_yi_corr(run_mode: RunMode, enabled: bool) -> Result<(), DynError> {
    if !enabled || !matches!(run_mode, RunMode::Corr) {
        return Ok(());
    }
    let now_sec = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("system time error: {e}"))?
        .as_secs() as i64;
    let tag = unix_seconds_to_yyyydddhhmmss(now_sec)?;
    let log_dir = std::env::current_dir()?.join("stdout");
    std::fs::create_dir_all(&log_dir)?;
    let log_path = log_dir.join(format!("stdout_{}.log", tag));
    let mut writer = BufWriter::new(File::create(&log_path)?);
    writeln!(writer, "# yi-corr stdout log started at {}", tag)?;
    writer.flush()?;
    {
        let mut slot = stdout_log_writer()
            .lock()
            .map_err(|_| "stdout log lock poisoned")?;
        *slot = Some(writer);
    }
    println!("[info] stdout log dir: {}", log_dir.display());
    println!("[info] stdout log: {}", log_path.display());
    Ok(())
}

fn carrier_phase_from_delay(f_hz: f64, tau_s: f64) -> Complex<f32> {
    if f_hz == 0.0 {
        Complex::new(1.0_f32, 0.0_f32)
    } else {
        let two_pi = 2.0 * std::f64::consts::PI;
        let phase = -two_pi * f_hz * tau_s;
        // Reduce phase before f32 conversion to keep precision for large absolute delays.
        let wrapped = ((phase + std::f64::consts::PI).rem_euclid(two_pi)) - std::f64::consts::PI;
        Complex::from_polar(1.0_f32, wrapped as f32)
    }
}

fn delay_seconds_at_time(
    d_s: f64,
    r_sps: f64,
    a_sps2: f64,
    j_sps3: f64,
    s_sps4: f64,
    t_s: f64,
) -> f64 {
    let t2 = t_s * t_s;
    let t3 = t2 * t_s;
    let t4 = t2 * t2;
    d_s + r_sps * t_s + 0.5 * a_sps2 * t2 + (1.0 / 6.0) * j_sps3 * t3 + (1.0 / 24.0) * s_sps4 * t4
}

fn split_delay_to_integer_and_fractional(delay_seconds: f64, fs_hz: f64) -> (i64, f64) {
    let integer_samples = (delay_seconds * fs_hz).round() as i64;
    let fractional_seconds = delay_seconds - (integer_samples as f64 / fs_hz);
    (integer_samples, fractional_seconds)
}

fn complete_fft_frame_count(duration_s: f64, fs_hz: f64, fft_len: usize) -> usize {
    let frames = duration_s * fs_hz / fft_len as f64;
    let nearest = frames.round();
    if (frames - nearest).abs() <= 1e-6 {
        nearest.max(0.0) as usize
    } else {
        frames.floor().max(0.0) as usize
    }
}

fn display_whole_seconds(duration_s: f64) -> i64 {
    let nearest = duration_s.round();
    if (duration_s - nearest).abs() <= 1e-5 {
        nearest as i64
    } else {
        duration_s.floor() as i64
    }
}

fn synth_timing_sample_stride() -> usize {
    static STRIDE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *STRIDE.get_or_init(|| {
        std::env::var("YI_TIMING_SAMPLE_STRIDE")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&v| v > 0)
            .unwrap_or(4096)
    })
}

#[inline]
fn map_real_fft_bin_to_xml_grid(
    src: &[Complex<f32>],
    fft_len: usize,
    raw_idx: isize,
) -> Complex<f32> {
    let n = fft_len as isize;
    let half = (fft_len / 2) as isize;
    let idx = raw_idx.rem_euclid(n);
    if idx <= half {
        src.get(idx as usize)
            .copied()
            .unwrap_or_else(|| Complex::new(0.0, 0.0))
    } else {
        src.get((n - idx) as usize)
            .map(|v| v.conj())
            .unwrap_or_else(|| Complex::new(0.0, 0.0))
    }
}

#[derive(Clone, Copy)]
struct RealFftGridMap {
    src_idx: usize,
    conjugate: bool,
}

fn build_real_fft_xml_grid_map(
    fft_len: usize,
    rotation_bins: isize,
    extra_raw_offset: isize,
    bins: usize,
) -> Vec<RealFftGridMap> {
    let n = fft_len as isize;
    let half = (fft_len / 2) as isize;
    (0..bins)
        .map(|xml_bin| {
            let raw_idx = xml_bin as isize - rotation_bins + extra_raw_offset;
            let idx = raw_idx.rem_euclid(n);
            if idx <= half {
                RealFftGridMap {
                    src_idx: idx as usize,
                    conjugate: false,
                }
            } else {
                RealFftGridMap {
                    src_idx: (n - idx) as usize,
                    conjugate: true,
                }
            }
        })
        .collect()
}

#[inline]
fn mapped_real_fft_bin(
    src: &[Complex<f32>],
    map: &[RealFftGridMap],
    xml_bin: usize,
) -> Complex<f32> {
    let m = map[xml_bin];
    let z = src[m.src_idx];
    if m.conjugate {
        z.conj()
    } else {
        z
    }
}

/// Return the first source bin when an XML-grid range maps to one forward,
/// non-conjugated, contiguous FFT range. This common case lets the accumulation
/// hot loop read spectra directly instead of loading and branching on two map
/// entries for every frequency bin of every frame.
fn contiguous_real_fft_src_start(
    map: &[RealFftGridMap],
    xml_start: usize,
    len: usize,
) -> Option<usize> {
    let xml_end = xml_start.checked_add(len)?;
    let range = map.get(xml_start..xml_end)?;
    let first = range.first()?;
    if first.conjugate {
        return None;
    }
    let src_start = first.src_idx;
    if range
        .iter()
        .enumerate()
        .all(|(offset, entry)| !entry.conjugate && entry.src_idx == src_start + offset)
    {
        Some(src_start)
    } else {
        None
    }
}

#[inline(always)]
fn accumulate_direct_acf_xcf(
    spectrum1: &[Complex<f32>],
    spectrum2: &[Complex<f32>],
    acc11: &mut [f64],
    acc12: &mut [Complex<f64>],
    acc22: &mut [f64],
    fr_mix: Complex<f32>,
    mut phase_corr: Complex<f32>,
    phase_step: Complex<f32>,
) {
    debug_assert_eq!(spectrum1.len(), spectrum2.len());
    debug_assert_eq!(spectrum1.len(), acc11.len());
    debug_assert_eq!(spectrum1.len(), acc12.len());
    debug_assert_eq!(spectrum1.len(), acc22.len());
    for k in 0..spectrum1.len() {
        let z1 = spectrum1[k];
        let z2 = spectrum2[k];
        acc11[k] += z1.norm_sqr() as f64;
        acc22[k] += z2.norm_sqr() as f64;
        let raw_xcf = z1 * z2.conj();
        let value = raw_xcf * fr_mix * phase_corr;
        acc12[k] += Complex::new(value.re as f64, value.im as f64);
        phase_corr *= phase_step;
    }
}

#[cfg(test)]
mod correlation_hot_path_tests {
    use super::*;

    #[test]
    fn contiguous_grid_map_detects_direct_and_wrapped_ranges() {
        let direct = build_real_fft_xml_grid_map(16, 0, 2, 9);
        assert_eq!(contiguous_real_fft_src_start(&direct, 1, 5), Some(3));

        let wrapped = build_real_fft_xml_grid_map(16, 1, 0, 9);
        assert_eq!(contiguous_real_fft_src_start(&wrapped, 0, 5), None);
        assert_eq!(contiguous_real_fft_src_start(&wrapped, 1, 5), Some(0));
    }

    #[test]
    fn direct_accumulator_is_bit_identical_to_mapped_loop() {
        let spectrum1 = vec![
            Complex::new(1.0, -0.5),
            Complex::new(-2.0, 0.25),
            Complex::new(0.75, 1.5),
            Complex::new(-0.125, -3.0),
        ];
        let spectrum2 = vec![
            Complex::new(-0.25, 2.0),
            Complex::new(1.5, -1.0),
            Complex::new(-0.5, -0.75),
            Complex::new(2.25, 0.125),
        ];
        let map = (0..spectrum1.len())
            .map(|src_idx| RealFftGridMap {
                src_idx,
                conjugate: false,
            })
            .collect::<Vec<_>>();
        let fr_mix = Complex::from_polar(1.0_f32, 0.37);
        let phase0 = Complex::from_polar(1.0_f32, -0.21);
        let phase_step = Complex::from_polar(1.0_f32, 0.013);

        let mut direct11 = vec![0.0; spectrum1.len()];
        let mut direct12 = vec![Complex::new(0.0, 0.0); spectrum1.len()];
        let mut direct22 = vec![0.0; spectrum1.len()];
        accumulate_direct_acf_xcf(
            &spectrum1,
            &spectrum2,
            &mut direct11,
            &mut direct12,
            &mut direct22,
            fr_mix,
            phase0,
            phase_step,
        );

        let mut mapped11 = vec![0.0; spectrum1.len()];
        let mut mapped12 = vec![Complex::new(0.0, 0.0); spectrum1.len()];
        let mut mapped22 = vec![0.0; spectrum1.len()];
        let mut phase_corr = phase0;
        for k in 0..spectrum1.len() {
            let z1 = mapped_real_fft_bin(&spectrum1, &map, k);
            let z2 = mapped_real_fft_bin(&spectrum2, &map, k);
            mapped11[k] += z1.norm_sqr() as f64;
            mapped22[k] += z2.norm_sqr() as f64;
            let raw_xcf = z1 * z2.conj();
            let value = raw_xcf * fr_mix * phase_corr;
            mapped12[k] += Complex::new(value.re as f64, value.im as f64);
            phase_corr *= phase_step;
        }

        assert_eq!(direct11, mapped11);
        assert_eq!(direct12, mapped12);
        assert_eq!(direct22, mapped22);
    }

    #[test]
    fn beam_weights_preserve_ratio_and_output_noise_power() {
        let (w1, w2, scale) = normalized_beam_weights(2.0, 1.0, 1.25, 1.25, 1.25).unwrap();
        assert!((w1 / w2 - 2.0).abs() < 1.0e-12);
        assert!(scale > 0.0);
        let power = w1 * w1 * 1.25 + w2 * w2 * 1.25;
        assert!((power - 1.25).abs() < 1.0e-12);

        let (tiny1, tiny2, _) = normalized_beam_weights(2.0e-6, 1.0e-6, 1.25, 1.25, 1.25).unwrap();
        assert!((tiny1 - w1).abs() < 1.0e-12);
        assert!((tiny2 - w2).abs() < 1.0e-12);
    }

    #[test]
    fn vsrec_phased_raw_round_trips_through_decoder() {
        let levels = vec![-1.5, -0.5, 0.5, 1.5];
        let samples = vec![
            -1.5_f32, -0.5, 0.5, 1.5, 1.5, 0.5, -0.5, -1.5, -1.5, 0.5, -0.5, 1.5, 0.5, -1.5, 1.5,
            -0.5,
        ];
        let vsrec = parse_shuffle(
            "24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7",
        )
        .unwrap();
        let mut packed = Vec::new();
        quantise_frame(&samples, 2, &levels, &vsrec, &mut packed).unwrap();
        assert_eq!(packed.len(), 4);

        let plan = build_decode_plan(2, &vsrec, &levels).unwrap();
        let mut decoded = vec![0.0_f32; samples.len()];
        decode_block_into_with_plan(&packed, samples.len(), &plan, &mut decoded, false, false)
            .unwrap();
        assert_eq!(decoded, samples);

        // LSB input is normalized to USB internally by negating odd samples.
        // Applying the same involution before output packing must reproduce
        // the original native VSREC bytes exactly.
        let mut internal_usb = vec![0.0_f32; samples.len()];
        decode_block_into_with_plan(
            &packed,
            samples.len(),
            &plan,
            &mut internal_usb,
            true,
            false,
        )
        .unwrap();
        for odd in internal_usb.iter_mut().skip(1).step_by(2) {
            *odd = -*odd;
        }
        let mut restored = Vec::new();
        quantise_frame(&internal_usb, 2, &levels, &vsrec, &mut restored).unwrap();
        assert_eq!(restored, packed);
    }

    #[test]
    fn native_word_aligned_seek_preserves_vsrec_sample_order() {
        let levels = vec![-1.5, -0.5, 0.5, 1.5];
        let vsrec = parse_shuffle(
            "24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7",
        )
        .unwrap();
        let samples = (0..64)
            .map(|i| levels[(i * 3 + i / 7) % levels.len()] as f32)
            .collect::<Vec<_>>();
        let mut packed = Vec::new();
        quantise_frame(&samples, 2, &levels, &vsrec, &mut packed).unwrap();

        let plan = build_decode_plan(2, &vsrec, &levels).unwrap();
        let mut decoded = vec![0.0_f32; 32];
        let mut scratch = DecodeWindowScratch::new();
        // Requested sample 13 is represented by a native-word read at sample
        // zero plus an integer delay shift. No physical bitstream repacking is
        // allowed before the VSREC word shuffle.
        decode_shifted_frame_from_chunk(
            &packed[..12],
            0,
            0,
            32,
            2,
            16,
            &plan,
            false,
            -13,
            &mut decoded,
            &mut scratch,
        )
        .unwrap();
        assert_eq!(decoded, samples[13..45]);
    }
}

#[inline]
fn station_grid_origin_offset_hz(name: &str) -> f64 {
    // Empirical backend/grid-origin convention correction in physical
    // frequency units, not in FFT-bin units.
    //
    // G9.62 maser validation with fs=1024 MHz and FFT=1048576 showed that
    // HITACH32 is displaced by about -1 high-resolution bin relative to
    // YAMAGU32/YAMAGU34 and GICO3.  One such bin is:
    //
    //   1024e6 / 1048576 = 976.5625 Hz
    //
    // Therefore the station/backend convention is represented as a small
    // frequency-origin offset.  It must not be applied as a fixed -1 bin for
    // all FFT lengths.  For example, FFT=8192 has df=125 kHz, where this
    // correction rounds to 0 bin.
    match name {
        "HITACH32" => -976.5625,
        _ => 0.0,
    }
}

#[inline]
fn station_grid_origin_offset_bins(name: &str, fs_hz: f64, fft_len: usize) -> isize {
    let df_hz = fs_hz / fft_len as f64;
    if !df_hz.is_finite() || df_hz <= 0.0 {
        return 0;
    }
    (station_grid_origin_offset_hz(name) / df_hz).round() as isize
}

fn ant2_grid_extra_offset() -> isize {
    static OFFSET: std::sync::OnceLock<isize> = std::sync::OnceLock::new();
    *OFFSET.get_or_init(|| {
        std::env::var("YI_ANT2_GRID_OFFSET")
            .ok()
            .and_then(|v| v.trim().parse::<isize>().ok())
            .unwrap_or(0)
    })
}

fn shift_real_fft_to_xml_grid_with_extra_offset(
    src: &[Complex<f32>],
    dst: &mut [Complex<f32>],
    fft_len: usize,
    rotation_bins: isize,
    extra_raw_offset: isize,
) {
    for (out_idx, out) in dst.iter_mut().enumerate() {
        let raw_idx = out_idx as isize - rotation_bins + extra_raw_offset;
        *out = map_real_fft_bin_to_xml_grid(src, fft_len, raw_idx);
    }
}
#[inline]
fn maybe_dump_xcf_debug(
    args_debug: bool,
    frame_global: usize,
    local_frame: usize,
    sector_idx: usize,
    output_grid: OutputGrid,
    bin: usize,
    raw_xcf: Complex<f32>,
    fr_mix: Complex<f32>,
    phase_corr: Complex<f32>,
    value: Complex<f32>,
    d: FrameDelayEntry,
    fs: f64,
) {
    if !args_debug || frame_global >= 40 {
        return;
    }
    let dbg_bin = std::env::var("YI_XCFDBG2_BIN")
        .ok()
        .and_then(|v| v.trim().parse::<usize>().ok())
        .unwrap_or(69866);
    if bin != dbg_bin {
        return;
    }
    let v_fr = raw_xcf * fr_mix;
    let v_pc = raw_xcf * phase_corr;
    let grid_label = match output_grid {
        OutputGrid::Ant1 => "ant1",
        OutputGrid::Ant2 => "ant2",
    };
    eprintln!(
        "[xcfdbg2] grid={} frame={} sec={} local_frame={} bin={} raw_abs={:.9e} raw_ph={:+.3} fr_ph={:+.3} pc_ph={:+.3} all_ph={:+.3} fr_mix_ph={:+.3} phase_corr_ph={:+.3} int=({},{}) frac_sample=({:+.6},{:+.6}) tau_sample=({:+.6},{:+.6})",
        grid_label,
        frame_global,
        sector_idx,
        local_frame,
        bin,
        raw_xcf.norm(),
        raw_xcf.arg().to_degrees(),
        v_fr.arg().to_degrees(),
        v_pc.arg().to_degrees(),
        value.arg().to_degrees(),
        fr_mix.arg().to_degrees(),
        phase_corr.arg().to_degrees(),
        d.int1,
        d.int2,
        d.frac1 * fs,
        d.frac2 * fs,
        d.tau1_s * fs,
        d.tau2_s * fs,
    );
}

#[inline]
fn fft_peak_dbg_enabled() -> bool {
    matches!(
        std::env::var("YI_FFTPEAKDBG")
            .unwrap_or_else(|_| "0".to_string())
            .trim()
            .to_ascii_lowercase()
            .as_str(),
        "1" | "true" | "yes" | "on"
    )
}

#[inline]
fn fft_peak_dbg_max_frames() -> usize {
    std::env::var("YI_FFTPEAKDBG_MAX")
        .ok()
        .and_then(|v| v.trim().parse::<usize>().ok())
        .unwrap_or(4)
}

#[inline]
fn fft_peak_dbg_range(default_start: usize, default_end: usize) -> (usize, usize) {
    if let Ok(v) = std::env::var("YI_FFTPEAKDBG_BINS") {
        let t = v.trim();
        if let Some((a, b)) = t.split_once(':') {
            if let (Ok(start), Ok(end)) = (a.trim().parse::<usize>(), b.trim().parse::<usize>()) {
                return (start, end.max(start + 1));
            }
        }
    }
    (default_start, default_end)
}

#[inline]
fn complex_peak_in_range(
    xs: &[Complex<f32>],
    start: usize,
    end: usize,
) -> Option<(usize, f64, f64)> {
    let end = end.min(xs.len());
    let start = start.min(end);
    if start >= end {
        return None;
    }

    let mut best_i = start;
    let mut best_abs = -1.0_f64;

    for i in start..end {
        let z = xs[i];
        let a = ((z.re as f64) * (z.re as f64) + (z.im as f64) * (z.im as f64)).sqrt();
        if a > best_abs {
            best_abs = a;
            best_i = i;
        }
    }

    let z = xs[best_i];
    let phase_deg = (z.im as f64).atan2(z.re as f64).to_degrees();
    Some((best_i, best_abs, phase_deg))
}

#[inline]
fn acf_peak_dbg_enabled() -> bool {
    matches!(
        std::env::var("YI_ACFPEAKDBG")
            .unwrap_or_else(|_| "0".to_string())
            .trim()
            .to_ascii_lowercase()
            .as_str(),
        "1" | "true" | "yes" | "on"
    )
}

#[inline]
fn acf_peak_dbg_range(default_start: usize, default_end: usize) -> (usize, usize) {
    if let Ok(v) = std::env::var("YI_ACFPEAKDBG_BINS") {
        let t = v.trim();
        if let Some((a, b)) = t.split_once(':') {
            if let (Ok(start), Ok(end)) = (a.trim().parse::<usize>(), b.trim().parse::<usize>()) {
                return (start, end.max(start + 1));
            }
        }
    }
    (default_start, default_end)
}

#[inline]
fn real_peak_in_range(xs: &[f64], start: usize, end: usize) -> Option<(usize, f64)> {
    let end = end.min(xs.len());
    let start = start.min(end);
    if start >= end {
        return None;
    }

    let mut best_i = start;
    let mut best_v = f64::NEG_INFINITY;
    for k in start..end {
        if xs[k] > best_v {
            best_v = xs[k];
            best_i = k;
        }
    }

    Some((best_i, best_v))
}

#[inline]
fn print_acf_peak_dbg(
    sector_index: usize,
    ant_name: &str,
    acc: &[f64],
    rotation_bins: isize,
    extra_offset: isize,
    fs_hz: f64,
    fft_len: usize,
) {
    let (b0, b1) = acf_peak_dbg_range(69840, 69890);

    // extra_offset is the diagnostic/env offset passed from the call site.
    // station_offset is the built-in backend/grid-origin correction.
    let station_offset = station_grid_origin_offset_bins(ant_name, fs_hz, fft_len);
    let env_offset = extra_offset;
    let total_offset = station_offset + env_offset;

    if let Some((pk, val)) = real_peak_in_range(acc, b0, b1) {
        eprintln!(
            "[acfpeakdbg] sector={} ant={} bins={}:{} rot={} station_offset={} env_offset={} total_offset={} peak={} pow={:.9e}",
            sector_index,
            ant_name,
            b0,
            b1,
            rotation_bins,
            station_offset,
            env_offset,
            total_offset,
            pk,
            val
        );
    } else {
        eprintln!(
            "[acfpeakdbg] sector={} ant={} bins={}:{} rot={} station_offset={} env_offset={} total_offset={} no_peak",
            sector_index,
            ant_name,
            b0,
            b1,
            rotation_bins,
            station_offset,
            env_offset,
            total_offset
        );
    }
}

fn print_fft_peak_dbg(
    tag: &str,
    frame_local: usize,
    ant_name: &str,
    raw: &[Complex<f32>],
    grid: &[Complex<f32>],
    rotation_bins: isize,
    extra_offset: isize,
) {
    let (b0, b1) = fft_peak_dbg_range(69840, 69890);

    let raw_peak = complex_peak_in_range(raw, b0, b1);
    let grid_peak = complex_peak_in_range(grid, b0, b1);

    match (raw_peak, grid_peak) {
        (Some((ri, ra, rp)), Some((gi, ga, gp))) => {
            eprintln!(
                "[fftpeakdbg] tag={} frame={} ant={} bins={}:{} rot={} extra={} raw_peak={} raw_abs={:.9e} raw_ph={:+.3} grid_peak={} grid_abs={:.9e} grid_ph={:+.3}",
                tag,
                frame_local,
                ant_name,
                b0,
                b1,
                rotation_bins,
                extra_offset,
                ri,
                ra,
                rp,
                gi,
                ga,
                gp,
            );
        }
        _ => {
            eprintln!(
                "[fftpeakdbg] tag={} frame={} ant={} bins={}:{} rot={} extra={} no_peak",
                tag, frame_local, ant_name, b0, b1, rotation_bins, extra_offset,
            );
        }
    }
}

#[inline]
fn ant_lsb_override(original: bool, var: &str) -> bool {
    match std::env::var(var)
        .unwrap_or_else(|_| "".to_string())
        .trim()
        .to_ascii_lowercase()
        .as_str()
    {
        "0" | "false" | "off" | "usb" => false,
        "1" | "true" | "on" | "lsb" => true,
        _ => original,
    }
}

fn fx_integer_delay_enabled() -> bool {
    // Production default: ON.
    //
    // yi-corr uses the same FX delay-correction path for short baselines,
    // long VLBI baselines, continuum sources, and maser sources.
    //
    // Normal path:
    //   integer sample delay -> raw/decode FFT-window shift
    //   fractional delay     -> frequency-domain phase slope
    //
    // Set YI_FX_INT_DELAY=0/false/no/off only when intentionally reproducing
    // the old legacy behavior for debugging.
    std::env::var("YI_FX_INT_DELAY")
        .ok()
        .map(|v| {
            !matches!(
                v.trim().to_ascii_lowercase().as_str(),
                "" | "0" | "false" | "no" | "off" | "legacy"
            )
        })
        .unwrap_or(true)
}

fn fx_read_offset_samples() -> i64 {
    std::env::var("YI_FX_READ_OFFSET")
        .ok()
        .and_then(|v| v.trim().parse::<i64>().ok())
        .unwrap_or(0)
}

fn env_f64_override(name: &str, default: f64) -> Result<f64, DynError> {
    match std::env::var(name) {
        Ok(raw) => {
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                Ok(default)
            } else {
                trimmed
                    .parse::<f64>()
                    .map_err(|e| format!("invalid {name}='{trimmed}': {e}").into())
            }
        }
        Err(_) => Ok(default),
    }
}

fn fx_delay_phase_offset_samples() -> f64 {
    std::env::var("YI_FX_DELAY_PHASE_OFFSET")
        .ok()
        .and_then(|v| v.trim().parse::<f64>().ok())
        .unwrap_or(0.0)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FringeStopFreqMode {
    Low,
    Center,
    High,
    Obs,
    Baseband,
    SidebandEdge,
}

impl FringeStopFreqMode {
    fn from_env() -> Result<Self, DynError> {
        let raw = std::env::var("YI_FRINGE_STOP_FREQ_MODE").unwrap_or_else(|_| "low".to_string());
        match raw.trim().to_ascii_lowercase().as_str() {
            "" | "low" | "data-low" | "current" | "default" => Ok(Self::Low),
            "center" | "centre" | "data-center" => Ok(Self::Center),
            "high" | "data-high" => Ok(Self::High),
            "obs" | "obsfreq" | "xml" | "ref" | "reference" => Ok(Self::Obs),
            "baseband" | "zero" | "none" => Ok(Self::Baseband),
            "sideband-edge" | "sb-edge" => Ok(Self::SidebandEdge),
            other => Err(format!(
                "invalid YI_FRINGE_STOP_FREQ_MODE='{other}' (use low, center, high, obs, baseband, or sideband-edge)"
            )
            .into()),
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Low => "low",
            Self::Center => "center",
            Self::High => "high",
            Self::Obs => "obsfreq",
            Self::Baseband => "baseband",
            Self::SidebandEdge => "sideband-edge",
        }
    }

    fn carrier_mhz(self, data_low_mhz: f64, bw_mhz: f64, obs_mhz: f64, is_lsb: bool) -> f64 {
        match self {
            Self::Low => data_low_mhz,
            Self::Center => data_low_mhz + 0.5 * bw_mhz,
            Self::High => data_low_mhz + bw_mhz,
            Self::Obs => obs_mhz,
            Self::Baseband => 0.0,
            Self::SidebandEdge => {
                if is_lsb {
                    data_low_mhz + bw_mhz
                } else {
                    data_low_mhz
                }
            }
        }
    }
}

fn xcf_phase_start_and_step(
    df_hz: f64,
    phase_delay1_s: f64,
    phase_delay2_s: f64,
    start_bin1: isize,
    start_bin2: isize,
) -> (Complex<f32>, Complex<f32>) {
    // Equivalent to rotating s1/s2 spectra independently and then forming
    // s1 * conj(s2).
    //
    // Apply the residual sub-sample delay left after integer sample tracking.
    // The large read-align component is already represented by the decoded
    // sample window and must not be applied again as a baseband phase slope.
    let two_pi_df = -2.0_f64 * std::f64::consts::PI * df_hz;
    let step1 = two_pi_df * phase_delay1_s;
    let step2 = two_pi_df * phase_delay2_s;
    let phase_start_angle = step1 * start_bin1 as f64 - step2 * start_bin2 as f64;
    let phase_step_angle = step1 - step2;
    let phase_start = Complex::from_polar(1.0_f32, phase_start_angle as f32);
    let phase_step = Complex::from_polar(1.0_f32, phase_step_angle as f32);
    (phase_start, phase_step)
}

fn antenna_phase_start_and_step(
    df_hz: f64,
    phase_delay_s: f64,
    start_bin: isize,
) -> (Complex<f32>, Complex<f32>) {
    let step_angle = -2.0_f64 * std::f64::consts::PI * df_hz * phase_delay_s;
    let phase_start = Complex::from_polar(1.0_f32, (step_angle * start_bin as f64) as f32);
    let phase_step = Complex::from_polar(1.0_f32, step_angle as f32);
    (phase_start, phase_step)
}

#[derive(Clone, Copy)]
struct FrameDelayEntry {
    t_mid_s: f64,
    full_rel_s: f64,
    residual_s: f64,
    tau1_s: f64,
    tau2_s: f64,
    int1: i64,
    int2: i64,
    frac1: f64,
    frac2: f64,
    fr_lo1: Complex<f32>,
    fr_lo2: Complex<f32>,
}

#[derive(Clone, Copy)]
struct GeomDelaySample {
    delay_s: f64,
    rate_sps: f64,
    accel_sps2: f64,
    jerk_sps3: f64,
    snap_sps4: f64,
}

struct DelayEvalConfig {
    frame_dt: f64,
    model_time_offset_s: f64,
    fs: f64,
    time_offset_s: f64,
    geom_delay_table_1s: Option<Arc<[GeomDelaySample]>>,
    geom_poly_order: u8,
    coarse_delay_s: f64,
    delay_user_samples: f64,
    extra_delay_rate_sps: f64,
    extra_delay_accel_sps2: f64,
    clock1_delay_s: f64,
    clock1_rate_sps: f64,
    clock1_accel_sps2: f64,
    clock1_jerk_sps3: f64,
    clock1_snap_sps4: f64,
    clock2_delay_s: f64,
    clock2_rate_sps: f64,
    clock2_accel_sps2: f64,
    clock2_jerk_sps3: f64,
    clock2_snap_sps4: f64,
    d_seek: f64,
    residual_on_ant2: bool,
    // Opt-in experiment: use process-start read alignment and update the
    // time-domain integer delay cumulatively whenever the residual exceeds
    // one sample.  Default is false to preserve the stable fixed-process
    // midpoint read-align behavior.
    fx_start_cumulative_seek: bool,
    net_d1_base: f64,
    total_rate1_base: f64,
    total_accel1_base: f64,
    total_jerk1_base: f64,
    total_snap1_base: f64,
    net_d2_base: f64,
    total_rate2_base: f64,
    total_accel2_base: f64,
    total_jerk2_base: f64,
    total_snap2_base: f64,
    lo1_hz: f64,
    lo2_hz: f64,
    fx_integer_delay: bool,
}

impl DelayEvalConfig {
    #[inline]
    fn residual_s_for_start_seek(&self, _frame_idx: usize, full_rel_s: f64, d_seek_s: f64) -> f64 {
        (full_rel_s - d_seek_s) * self.fs
    }
}

fn compute_frame_delay_entry(
    frame_idx: usize,
    cfg: &DelayEvalConfig,
    d_seek_s: f64,
) -> FrameDelayEntry {
    let t_mid_local = (frame_idx as f64 + 0.5) * cfg.frame_dt + cfg.model_time_offset_s;
    let t_mid = t_mid_local + cfg.time_offset_s;
    let (tau1, tau2, tau1_for_fringe, tau2_for_fringe, full_rel_dbg, residual_dbg) =
        if let Some(gd_table) = cfg.geom_delay_table_1s.as_deref() {
            let sec_idx = t_mid_local.floor().max(0.0) as usize;
            let i0 = sec_idx.min(gd_table.len().saturating_sub(1));
            let dt = t_mid_local - i0 as f64;
            let gd0 = gd_table[i0];
            let dt2 = dt * dt;
            let dt3 = dt2 * dt;
            let dt4 = dt2 * dt2;
            let gd_t = gd0.delay_s
                + gd0.rate_sps * dt
                + if cfg.geom_poly_order >= 2 {
                    0.5 * gd0.accel_sps2 * dt2
                } else {
                    0.0
                }
                + if cfg.geom_poly_order >= 3 {
                    (1.0 / 6.0) * gd0.jerk_sps3 * dt3
                } else {
                    0.0
                }
                + if cfg.geom_poly_order >= 4 {
                    (1.0 / 24.0) * gd0.snap_sps4 * dt4
                } else {
                    0.0
                };
            let net_d_rel_no_clock_t = gd_t
                + cfg.coarse_delay_s
                + cfg.delay_user_samples / cfg.fs
                + cfg.extra_delay_rate_sps * t_mid
                + 0.5 * cfg.extra_delay_accel_sps2 * t_mid * t_mid;
            let t2 = t_mid * t_mid;
            let t3 = t2 * t_mid;
            let t4 = t2 * t2;
            let clock1_t = cfg.clock1_delay_s
                + cfg.clock1_rate_sps * t_mid
                + 0.5 * cfg.clock1_accel_sps2 * t2
                + (1.0 / 6.0) * cfg.clock1_jerk_sps3 * t3
                + (1.0 / 24.0) * cfg.clock1_snap_sps4 * t4;
            let clock2_t = cfg.clock2_delay_s
                + cfg.clock2_rate_sps * t_mid
                + 0.5 * cfg.clock2_accel_sps2 * t2
                + (1.0 / 6.0) * cfg.clock2_jerk_sps3 * t3
                + (1.0 / 24.0) * cfg.clock2_snap_sps4 * t4;
            // Apply only relative clock delay/rate to avoid large common-mode shifts.
            // Large absolute per-station clocks can exceed small FFT frame length and
            // zero-fill whole frames even though the baseline-relative delay is small.
            let rel_clock_t = clock2_t - clock1_t;
            let full_rel_t = net_d_rel_no_clock_t + rel_clock_t;
            let residual_t = full_rel_t - d_seek_s;
            if cfg.residual_on_ant2 {
                (0.0, -residual_t, 0.0, -full_rel_t, full_rel_t, residual_t)
            } else {
                (residual_t, 0.0, full_rel_t, 0.0, full_rel_t, residual_t)
            }
        } else {
            let tau1_res = delay_seconds_at_time(
                cfg.net_d1_base,
                cfg.total_rate1_base,
                cfg.total_accel1_base,
                cfg.total_jerk1_base,
                cfg.total_snap1_base,
                t_mid,
            );
            let tau2_res = delay_seconds_at_time(
                cfg.net_d2_base,
                cfg.total_rate2_base,
                cfg.total_accel2_base,
                cfg.total_jerk2_base,
                cfg.total_snap2_base,
                t_mid,
            );
            if cfg.residual_on_ant2 {
                let residual_dbg = -tau2_res;
                let full_rel_dbg = residual_dbg + d_seek_s;
                (
                    tau1_res + d_seek_s,
                    tau2_res,
                    0.0,
                    tau2_res - d_seek_s,
                    full_rel_dbg,
                    residual_dbg,
                )
            } else {
                let residual_dbg = tau1_res;
                let full_rel_dbg = residual_dbg + d_seek_s;
                (
                    tau1_res,
                    tau2_res,
                    // Keep fringe phase tied to the pre-seek absolute delay so changing
                    // read-start alignment between scans does not reset phase origin.
                    tau1_res,
                    tau2_res,
                    full_rel_dbg,
                    residual_dbg,
                )
            }
        };
    // Split the delay correction into a time-domain integer sample shift and
    // a frequency-domain fractional phase slope.
    //
    // Correct invariant for every FFT frame:
    //
    //     tau*_s == int*/fs + frac*
    //     frac* in roughly [-0.5, +0.5] sample
    //
    // The previous geometric-delay-table path deliberately forced int=(0,0)
    // and left the full residual delay in frac.  That is not fractional delay;
    // on long baselines it can be thousands of samples and the raw data are not
    // time-domain delay-tracked.
    //
    // Set YI_FX_INT_DELAY=0 only for the old legacy/debug behavior.
    let (int1, frac1, int2, frac2) = if cfg.fx_start_cumulative_seek {
        if cfg.residual_on_ant2 {
            let residual_samples = cfg.residual_s_for_start_seek(frame_idx, full_rel_dbg, d_seek_s);
            let seek_delta_samples = residual_samples.round() as i64;
            let frac2_s = -(residual_samples - seek_delta_samples as f64) / cfg.fs;
            (0, 0.0, seek_delta_samples, frac2_s)
        } else {
            let residual_samples = cfg.residual_s_for_start_seek(frame_idx, full_rel_dbg, d_seek_s);
            let seek_delta_samples = residual_samples.round() as i64;
            let frac1_s = (residual_samples - seek_delta_samples as f64) / cfg.fs;
            (seek_delta_samples, frac1_s, 0, 0.0)
        }
    } else if cfg.fx_integer_delay {
        let (i1, f1) = split_delay_to_integer_and_fractional(tau1, cfg.fs);
        let (i2, f2) = split_delay_to_integer_and_fractional(tau2, cfg.fs);
        (i1, f1, i2, f2)
    } else if cfg.geom_delay_table_1s.is_some() {
        // Legacy fixed-process mode, enabled only by YI_FX_INT_DELAY=0.
        // This keeps the full residual delay in the post-FFT phase slope.
        (0, tau1, 0, tau2)
    } else {
        let (i1, f1) = split_delay_to_integer_and_fractional(tau1, cfg.fs);
        let (i2, f2) = split_delay_to_integer_and_fractional(tau2, cfg.fs);
        (i1, f1, i2, f2)
    };

    FrameDelayEntry {
        t_mid_s: t_mid,
        full_rel_s: full_rel_dbg,
        residual_s: residual_dbg,
        tau1_s: tau1,
        tau2_s: tau2,
        int1,
        int2,
        frac1,
        frac2,
        // Keep the carrier / fringe phase origin tied to the continuous
        // pre-read-align delay branch.  The baseband XCF phase slope uses
        // only the residual fractional delay left after integer sample tracking.
        fr_lo1: carrier_phase_from_delay(cfg.lo1_hz, tau1_for_fringe),
        fr_lo2: carrier_phase_from_delay(cfg.lo2_hz, tau2_for_fringe),
    }
}

fn print_delay_debug_samples(
    label: &str,
    start_frame: usize,
    frame_delays: &[FrameDelayEntry],
    fs: f64,
) {
    if frame_delays.is_empty() {
        return;
    }

    // Print a compact diagnostic set:
    //   - start / middle / end of this sector
    //   - the first frame of each UTC elapsed second covered by this sector
    //   - the frames immediately around any integer-delay transition
    // This is much lighter than dumping every FFT frame, but it catches the
    // ±0.5 sample boundary where round(delay*fs) changes d.int1/d.int2.
    let mut sample_indices = Vec::<usize>::new();
    let mut push_idx = |idx: usize| {
        if idx < frame_delays.len() {
            sample_indices.push(idx);
        }
    };

    push_idx(0);
    push_idx(frame_delays.len() / 2);
    push_idx(frame_delays.len() - 1);

    let mut last_whole_sec: Option<i64> = None;
    for (i, d) in frame_delays.iter().enumerate() {
        let whole_sec = d.t_mid_s.floor() as i64;
        if last_whole_sec != Some(whole_sec) {
            push_idx(i);
            last_whole_sec = Some(whole_sec);
        }

        if i > 0 {
            let p = frame_delays[i - 1];
            if d.int1 != p.int1 || d.int2 != p.int2 {
                push_idx(i - 1);
                push_idx(i);
                push_idx(i + 1);
            }
        }
    }

    sample_indices.sort_unstable();
    sample_indices.dedup();
    println!(
        "[debug] {label}: frames={} start_frame={} sampled={}",
        frame_delays.len(),
        start_frame,
        sample_indices.len()
    );
    for i in sample_indices {
        let d = frame_delays[i];
        println!(
            "[debug]   frame={} t_mid={:.9e} full_rel={:.9} residual={:.9} tau1={:.9} tau2={:.9} int=({},{}) frac=({:.9},{:.9}) sample",
            start_frame + i,
            d.t_mid_s,
            d.full_rel_s * fs,
            d.residual_s * fs,
            d.tau1_s * fs,
            d.tau2_s * fs,
            d.int1,
            d.int2,
            d.frac1 * fs,
            d.frac2 * fs
        );
    }
}

#[derive(Clone, Copy, Debug)]
enum OutputGrid {
    Ant1,
    Ant2,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum NormalCorrKernel {
    Ant1Grid,
    Ant2Grid,
}

impl NormalCorrKernel {
    fn from_output_grid(out_grid: OutputGrid) -> Self {
        match out_grid {
            OutputGrid::Ant1 => Self::Ant1Grid,
            OutputGrid::Ant2 => Self::Ant2Grid,
        }
    }

    fn output_grid(self) -> OutputGrid {
        match self {
            Self::Ant1Grid => OutputGrid::Ant1,
            Self::Ant2Grid => OutputGrid::Ant2,
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum RunMode {
    PhasedArray,
    Corr,
}

impl RunMode {
    fn label(self) -> &'static str {
        match self {
            RunMode::PhasedArray => "yi-phasedarray",
            RunMode::Corr => "yi-corr",
        }
    }
}

fn detect_run_mode_from_argv0() -> RunMode {
    let stem = std::env::args()
        .next()
        .and_then(|p| {
            PathBuf::from(p)
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
        })
        .unwrap_or_default()
        .to_ascii_lowercase();
    match stem.as_str() {
        "yi_phasedarray" | "yi-phasedarray" => RunMode::PhasedArray,
        "yi-corr" | "yi_corr" => RunMode::Corr,
        s if s.starts_with("yi-corr-v") => RunMode::Corr,
        _ => RunMode::Corr,
    }
}

fn build_logical_cpus() -> Option<usize> {
    option_env!("YI_BUILD_LOGICAL_CPUS")
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&v| v > 0)
}

fn runtime_physical_cores() -> Option<usize> {
    let available = core_affinity::get_core_ids()?;
    let mut physical = HashSet::<(usize, usize)>::new();
    for cpu in available {
        let topology = PathBuf::from(format!("/sys/devices/system/cpu/cpu{}/topology", cpu.id));
        let package = read_to_string(topology.join("physical_package_id"))
            .ok()?
            .trim()
            .parse::<usize>()
            .ok()?;
        let core = read_to_string(topology.join("core_id"))
            .ok()?
            .trim()
            .parse::<usize>()
            .ok()?;
        physical.insert((package, core));
    }
    (!physical.is_empty()).then_some(physical.len())
}

fn build_l3_cache_bytes() -> Option<u64> {
    option_env!("YI_BUILD_L3_CACHE_BYTES")
        .and_then(|v| v.parse::<u64>().ok())
        .filter(|&v| v > 0)
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
    num.trim()
        .parse::<u64>()
        .ok()
        .map(|v| v.saturating_mul(mul))
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
    let l3_bytes = build_l3_cache_bytes()
        .or_else(detect_l3_cache_bytes)
        .unwrap_or(32 * 1024 * 1024);
    let scale = (cpu_threads as u64).clamp(1, 16);
    let target_bytes =
        (l3_bytes.saturating_mul(scale) / 4).clamp(8 * 1024 * 1024, 128 * 1024 * 1024);
    let mut frames = (target_bytes / bytes_per_frame_pair as u64) as usize;
    frames = frames.clamp(1024, 32768);
    ((frames / 256).max(1)) * 256
}

fn auto_pipeline_depth(cpu_threads: usize) -> usize {
    (cpu_threads / 4).clamp(2, 16)
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

#[derive(Clone, Copy, Debug)]
struct BandAlignment {
    shift_bins: isize,
    a1s: usize,
    a1e: usize,
    a2s: usize,
}

fn compute_band_alignment(
    fft: usize,
    fs: f64,
    c1: f64,
    c2: f64,
    bw1: f64,
    bw2: f64,
) -> Result<BandAlignment, DynError> {
    let df = fs / fft as f64;
    let h = fft / 2 + 1;
    let lim = if fft % 2 == 0 { h - 1 } else { h };
    let v1 = (((bw1 * 1e6) / df).floor() as usize)
        .saturating_add(1)
        .min(lim);
    let v2 = (((bw2 * 1e6) / df).floor() as usize)
        .saturating_add(1)
        .min(lim);
    let sr = (c1 * 1e6 - c2 * 1e6) / df;
    let sb = sr.round() as isize;
    let a1s = 0isize.max(-sb) as usize;
    let a1e = (v1 as isize).min(v2 as isize - sb) as usize;
    if a1e <= a1s {
        return Err("No band overlap".into());
    }
    Ok(BandAlignment {
        shift_bins: sb,
        a1s,
        a1e,
        a2s: (a1s as isize + sb) as usize,
    })
}

fn format_bit_codes(b: usize) -> String {
    if b == 0 {
        "n/a".into()
    } else {
        (0..(1 << b).min(1024))
            .map(|c| format!("{c:0b$b}", b = b))
            .collect::<Vec<_>>()
            .join(" ")
    }
}
fn format_level_map(b: usize, l: &[f64]) -> String {
    if b == 0 {
        "n/a".into()
    } else {
        (0..(1 << b).min(l.len()))
            .map(|c| format!("{c:0b$b}->{:.1}", l[c], b = b))
            .collect::<Vec<_>>()
            .join(", ")
    }
}
fn format_shuffle_compact(values: &[usize]) -> String {
    values
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

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

fn normalized_beam_weights(
    weight1: f64,
    weight2: f64,
    input_power1: f64,
    input_power2: f64,
    output_power: f64,
) -> Result<(f64, f64, f64), DynError> {
    if !weight1.is_finite() || !weight2.is_finite() || weight1 <= 0.0 || weight2 <= 0.0 {
        return Err("phased-array weights must be positive finite values".into());
    }
    if !input_power1.is_finite()
        || !input_power2.is_finite()
        || !output_power.is_finite()
        || input_power1 <= 0.0
        || input_power2 <= 0.0
        || output_power <= 0.0
    {
        return Err("phased-array quantizer powers must be positive finite values".into());
    }
    let unscaled_power = weight1 * weight1 * input_power1 + weight2 * weight2 * input_power2;
    let voltage_scale = (output_power / unscaled_power).sqrt();
    Ok((
        weight1 * voltage_scale,
        weight2 * voltage_scale,
        voltage_scale,
    ))
}

fn normalize_level_args(
    level_args: &[String],
    bit1: usize,
    bit2: usize,
) -> Result<Vec<String>, DynError> {
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

fn sample_to_total_bits(sample_offset: u64, bits_per_sample: usize) -> Result<u64, DynError> {
    if bits_per_sample == 0 {
        return Err("bits_per_sample must be >= 1".into());
    }
    sample_offset
        .checked_mul(bits_per_sample as u64)
        .ok_or("sample offset overflow while converting to bit offset".into())
}

fn div_floor_i128(v: i128, d: i128) -> i128 {
    let q = v / d;
    let r = v % d;
    if r != 0 && ((r > 0) != (d > 0)) {
        q - 1
    } else {
        q
    }
}

struct DecodeWindowScratch {
    raw: Vec<u8>,
    samples: Vec<f32>,
}

impl DecodeWindowScratch {
    fn new() -> Self {
        Self {
            raw: Vec::new(),
            samples: Vec::new(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn decode_shifted_frame_from_chunk(
    raw_chunk: &[u8],
    chunk_abs_start_sample: u64,
    frame_in_chunk: usize,
    fft_len: usize,
    bit: usize,
    samples_per_word: u64,
    plan: &DecodePlan,
    lsb_to_usb: bool,
    int_shift: i64,
    out: &mut [f32],
    scratch: &mut DecodeWindowScratch,
) -> Result<(), DynError> {
    if out.len() != fft_len {
        return Err("decode output length does not match fft length".into());
    }
    out.fill(0.0);
    let chunk_samples = (raw_chunk.len() * 8 / bit) as i128;
    let frame_start = frame_in_chunk as i128 * fft_len as i128;
    let needed_start = frame_start - int_shift as i128;
    let needed_end = needed_start + fft_len as i128;
    let copy_start = needed_start.max(0);
    let copy_end = needed_end.min(chunk_samples);
    if copy_end <= copy_start {
        return Ok(());
    }

    let out_start = (copy_start - needed_start) as usize;
    let copy_len = (copy_end - copy_start) as usize;
    let spw = samples_per_word as i128;
    let aligned_start = div_floor_i128(copy_start, spw) * spw;
    let offset_in_decoded = (copy_start - aligned_start) as usize;
    let decode_samples_needed = offset_in_decoded + copy_len;
    let decode_samples = ((decode_samples_needed as u64 + samples_per_word - 1) / samples_per_word
        * samples_per_word) as usize;
    let raw_byte_start = (aligned_start as usize * bit) / 8;
    let raw_bytes_needed = (decode_samples * bit + 7) / 8;

    scratch.raw.resize(raw_bytes_needed, 0);
    if raw_byte_start < raw_chunk.len() {
        let available = (raw_chunk.len() - raw_byte_start).min(raw_bytes_needed);
        scratch.raw[..available]
            .copy_from_slice(&raw_chunk[raw_byte_start..raw_byte_start + available]);
        scratch.raw[available..].fill(0);
    } else {
        scratch.raw.fill(0);
    }

    scratch.samples.resize(decode_samples, 0.0);
    let aligned_abs = chunk_abs_start_sample as u128 + aligned_start as u128;
    let first_sample_odd = (aligned_abs & 1) != 0;
    decode_block_into_with_plan(
        &scratch.raw,
        decode_samples,
        plan,
        &mut scratch.samples,
        lsb_to_usb,
        first_sample_odd,
    )?;
    out[out_start..out_start + copy_len]
        .copy_from_slice(&scratch.samples[offset_in_decoded..offset_in_decoded + copy_len]);
    Ok(())
}

fn read_with_padding(reader: &mut BufReader<File>, buf: &mut [u8]) -> Result<usize, DynError> {
    let mut total = 0;
    while total < buf.len() {
        match reader.read(&mut buf[total..]) {
            Ok(0) => {
                buf[total..].fill(0);
                return Ok(total);
            }
            Ok(n) => total += n,
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    Ok(total)
}

struct PackedSampleReader {
    reader: BufReader<File>,
    bit_offset: u8,
    prev_byte: Option<u8>,
    next_start_byte: Option<u64>,
}

impl PackedSampleReader {
    fn open(path: &PathBuf, start_byte: u64, bit_offset: u8) -> Result<Self, DynError> {
        if bit_offset >= 8 {
            return Err("bit_offset must be in 0..8".into());
        }
        let f = File::open(path)?;
        advise_sequential_readahead(&f);
        let mut reader = BufReader::new(f);
        reader.seek(SeekFrom::Start(start_byte))?;
        let prev_byte = if bit_offset == 0 {
            None
        } else {
            let mut b = [0u8; 1];
            let _ = read_with_padding(&mut reader, &mut b)?;
            Some(b[0])
        };
        Ok(Self {
            reader,
            bit_offset,
            prev_byte,
            next_start_byte: Some(start_byte),
        })
    }

    fn seek_to(&mut self, start_byte: u64, bit_offset: u8) -> Result<(), DynError> {
        if bit_offset >= 8 {
            return Err("bit_offset must be in 0..8".into());
        }
        if self.next_start_byte == Some(start_byte) && self.bit_offset == bit_offset {
            return Ok(());
        }
        self.reader.seek(SeekFrom::Start(start_byte))?;
        self.bit_offset = bit_offset;
        self.prev_byte = if bit_offset == 0 {
            None
        } else {
            let mut b = [0u8; 1];
            let _ = read_with_padding(&mut self.reader, &mut b)?;
            Some(b[0])
        };
        self.next_start_byte = Some(start_byte);
        Ok(())
    }

    fn read_packed_with_padding(&mut self, out: &mut [u8]) -> Result<(), DynError> {
        if self.bit_offset == 0 {
            let _ = read_with_padding(&mut self.reader, out)?;
            if let Some(next) = self.next_start_byte.as_mut() {
                *next += out.len() as u64;
            }
            return Ok(());
        }
        let _ = read_with_padding(&mut self.reader, out)?;
        let shift = self.bit_offset;
        let mut prev = self.prev_byte.unwrap_or(0);
        for dst in out.iter_mut() {
            let next = *dst;
            *dst = (prev >> shift) | (next << (8 - shift));
            prev = next;
        }
        self.prev_byte = Some(prev);
        if let Some(next) = self.next_start_byte.as_mut() {
            *next += out.len() as u64;
        }
        Ok(())
    }
}

type PairedSectorRead = Result<(usize, Vec<u8>, Vec<u8>), String>;
type StationSectorRead = Result<(usize, Vec<u8>), String>;

struct SynthReaderPipeline {
    receiver: Option<mpsc::Receiver<PairedSectorRead>>,
    handles: Vec<(&'static str, thread::JoinHandle<()>)>,
}

impl SynthReaderPipeline {
    fn new(
        receiver: mpsc::Receiver<PairedSectorRead>,
        handles: Vec<(&'static str, thread::JoinHandle<()>)>,
    ) -> Self {
        Self {
            receiver: Some(receiver),
            handles,
        }
    }

    fn recv(&self) -> Result<PairedSectorRead, mpsc::RecvError> {
        self.receiver
            .as_ref()
            .expect("synth reader receiver is available")
            .recv()
    }

    fn finish(mut self) -> Result<(), DynError> {
        self.receiver.take();
        let mut panicked = Vec::new();
        for (name, handle) in self.handles.drain(..) {
            if handle.join().is_err() {
                panicked.push(name);
            }
        }
        if panicked.is_empty() {
            Ok(())
        } else {
            Err(format!("synth reader thread(s) panicked: {}", panicked.join(", ")).into())
        }
    }
}

impl Drop for SynthReaderPipeline {
    fn drop(&mut self) {
        // Dropping the ready receiver first releases a producer blocked on a
        // bounded send.  In USB mode that also lets the pairer drop both
        // rendezvous receivers before the station-reader joins below.
        self.receiver.take();
        for (_, handle) in self.handles.drain(..) {
            let _ = handle.join();
        }
    }
}

fn run_usb_station_reader(
    path: PathBuf,
    station_label: &'static str,
    bit_depth: usize,
    sector_reads: Vec<(u64, u64)>,
    reader_core: Option<core_affinity::CoreId>,
    sender: mpsc::SyncSender<StationSectorRead>,
) {
    if let Some(core) = reader_core {
        let _ = affinity::set_current_thread_core(core);
    }
    let mut reader = match PackedSampleReader::open(&path, 0, 0) {
        Ok(v) => v,
        Err(e) => {
            let _ = sender.send(Err(format!(
                "failed to open {} for USB {} sector reader: {e}",
                path.display(),
                station_label
            )));
            return;
        }
    };
    for (sector_idx, (start_sample, sample_count)) in sector_reads.into_iter().enumerate() {
        let mut buf = vec![0u8; ((sample_count as usize * bit_depth) + 7) / 8];
        let start_bits = start_sample * bit_depth as u64;
        if let Err(e) = reader.seek_to(start_bits / 8, (start_bits % 8) as u8) {
            let _ = sender.send(Err(format!(
                "failed to seek {} USB {} input at sample {} (sector {}): {e}",
                path.display(),
                station_label,
                start_sample,
                sector_idx + 1
            )));
            return;
        }
        if let Err(e) = reader.read_packed_with_padding(&mut buf) {
            let _ = sender.send(Err(format!(
                "failed reading {} USB {} input at sample {} (sector {}): {e}",
                path.display(),
                station_label,
                start_sample,
                sector_idx + 1
            )));
            return;
        }
        if sender.send(Ok((sector_idx, buf))).is_err() {
            return;
        }
    }
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

fn resolve_input_paths(
    args: &args::Args,
    pe: &str,
    meta: Option<&ifile::IFileData>,
) -> Result<(PathBuf, PathBuf, String, String), DynError> {
    if let (Some(a1), Some(a2)) = (&args.ant1, &args.ant2) {
        let (_, tag) = epoch_to_yyyydddhhmmss(pe)?;
        return Ok((a1.clone(), a2.clone(), pe.to_string(), tag));
    }
    let data_dir = args
        .raw_directory
        .clone()
        .unwrap_or_else(|| PathBuf::from("."));
    let (_, tag) = epoch_to_yyyydddhhmmss(pe)?;
    let mut candidates = vec![("YAMAGU32", "YAMAGU34")];
    if let Some(m) = meta {
        if let (Some(n1), Some(n2)) = (&m.ant1_station_name, &m.ant2_station_name) {
            candidates.insert(0, (n1, n2));
        }
    }
    for (p1, p2) in candidates {
        let a1 = data_dir.join(format!("{}_{}.raw", p1, tag));
        let a2 = data_dir.join(format!("{}_{}.raw", p2, tag));
        if a1.exists() && a2.exists() {
            println!(
                "[info] Auto-resolved inputs: {} / {}",
                a1.display(),
                a2.display()
            );
            return Ok((a1, a2, pe.to_string(), tag));
        }
    }
    Err("Input files not found".into())
}

fn resolve_output_layout(
    args: &args::Args,
    tag: &str,
    run_mode: RunMode,
) -> Result<(PathBuf, PathBuf), DynError> {
    let phased_name = sanitize_file_token(&args.phased_name);
    let stem = format!("{phased_name}_{tag}");
    let dir = if matches!(run_mode, RunMode::PhasedArray) {
        let d = args
            .cor_directory
            .clone()
            .unwrap_or(std::env::current_dir()?);
        std::fs::create_dir_all(&d)?;
        d
    } else {
        let d = args
            .cor_directory
            .clone()
            .ok_or("--cor-directory is required for yi-corr")?;
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
    v.map(|x| format!("{x:.9e}"))
        .unwrap_or_else(|| "-".to_string())
}

fn fmt_opt_mhz(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.3} MHz"))
        .unwrap_or_else(|| "-".to_string())
}

fn fmt_mhz(v: f64) -> String {
    format!("{v:.3} MHz")
}
fn fmt_hz_to_mhz(v_hz: f64) -> String {
    format!("{:.3} MHz", v_hz / 1e6)
}
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
    geom_eval_offset_s: f64,
    read_align_ref_offset_s: f64,
    geom_rate_hz_at_obs: f64,
    rel_clock_delay_s: f64,
    rel_clock_rate_sps: f64,
    clock_epoch_mode: &'a str,
    coarse_delay_s: f64,
    coarse_delay_samples: f64,
    res_delay_input_samples: f64,
    read_align_delay_samples: f64,
    read_align_delay_s: f64,
    delay_model: &'a str,
    model_time_offset_s: f64,
    dut1_s: f64,
    tt_utc_s: f64,
    xp_arcsec: f64,
    yp_arcsec: f64,
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
    clock_epoch: Option<&'a str>,
    clock_epoch_dt_s: Option<f64>,
    clock_delay_raw_s: f64,
    clock_rate_raw_sps: f64,
    clock_accel_raw_sps2: f64,
    clock_jerk_raw_sps3: f64,
    clock_snap_raw_sps4: f64,
    clock_delay_s: f64,
    clock_rate_sps: f64,
    clock_accel_sps2: f64,
    clock_jerk_sps3: f64,
    clock_snap_sps4: f64,
    center_mhz: Option<f64>,
    bw_mhz: Option<f64>,
}

fn print_kv_table(title: &str, rows: &[(String, String)]) {
    let key_width = rows.iter().map(|(k, _)| k.len()).max().unwrap_or(0);
    println!("{title}");
    for (key, value) in rows {
        println!("  {key:<key_width$}: {value}");
    }
}

fn abbreviate_middle(s: &str, max_len: usize) -> String {
    if s.chars().count() <= max_len {
        return s.to_string();
    }
    if max_len <= 3 {
        return "...".chars().take(max_len).collect();
    }
    let keep = max_len - 3;
    let left = keep / 2;
    let right = keep - left;
    let prefix: String = s.chars().take(left).collect();
    let suffix: String = s
        .chars()
        .rev()
        .take(right)
        .collect::<Vec<_>>()
        .into_iter()
        .rev()
        .collect();
    format!("{prefix}...{suffix}")
}

fn path_file_display(path: &PathBuf) -> String {
    path.file_name()
        .and_then(|s| s.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| path.display().to_string())
}

fn print_ant_table(header_l: &str, header_r: &str, rows: &[(String, String, String)]) {
    let key_width = rows.iter().map(|(p, _, _)| p.len()).max().unwrap_or(0);
    println!("Antenna Parameters");
    println!("  {header_l}");
    for (param, left, _) in rows {
        println!("    {param:<key_width$}: {left}");
    }
    println!("  {header_r}");
    for (param, _, right) in rows {
        println!("    {param:<key_width$}: {right}");
    }
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
        let baselines = d
            .processes
            .iter()
            .filter_map(|p| {
                let a = p.ant1_station_key.as_deref()?;
                let b = p.ant2_station_key.as_deref()?;
                Some(format!("{a}{b}"))
            })
            .collect::<Vec<_>>();
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
            (
                "ra/dec (J2000)".to_string(),
                format!("{} / {}", d.ra, d.dec),
            ),
            (
                "process epochs (all)".to_string(),
                if d.process_epochs.is_empty() {
                    "-".to_string()
                } else {
                    d.process_epochs.join(", ")
                },
            ),
            (
                "process skip/length [s]".to_string(),
                format!(
                    "{} / {}",
                    gv.process_skip_sec
                        .map(|v| format!("{v:.3}"))
                        .unwrap_or_else(|| "-".to_string()),
                    gv.process_length_sec
                        .map(|v| format!("{v:.3}"))
                        .unwrap_or_else(|| "-".to_string())
                ),
            ),
            (
                "baselines".to_string(),
                if baselines.is_empty() {
                    format!("{}{}", a1.key, a2.key)
                } else {
                    baselines.join(", ")
                },
            ),
            ("stream fft".to_string(), fmt_opt(d.fft)),
            ("stream output [s]".to_string(), fmt_opt_f64(d.output_sec)),
            ("stream inband".to_string(), fmt_opt(d.inband)),
            (
                "stream output rate".to_string(),
                d.output_sec
                    .map(|v| format!("{v:.9} Hz ({:.9} s)", 1.0 / v))
                    .unwrap_or_else(|| "1.000000000 Hz (1.000000000 s)".to_string()),
            ),
            (
                "stream sampling".to_string(),
                format!(
                    "xml {} / runtime {}",
                    fmt_opt_f64(d.sampling_hz),
                    fmt_mhz(gv.sampling_mhz)
                ),
            ),
            (
                "stream obsfreq".to_string(),
                format!(
                    "xml {} / runtime {}",
                    fmt_opt_mhz(d.obsfreq_mhz),
                    fmt_mhz(gv.obsfreq_mhz)
                ),
            ),
            ("stream bw (runtime)".to_string(), fmt_mhz(gv.bw_mhz)),
            (
                "delay convention".to_string(),
                "read-align ant2-ant1; residual on delayed antenna".to_string(),
            ),
        ];
        print_kv_table("Schedule / Stream Parameters", &global_rows);
    }

    let baseline_rows = vec![
        (
            "baseline".to_string(),
            format!("{}{} ({} - {})", a1.key, a2.key, a1.name, a2.name),
        ),
        (
            "geom delay/rate".to_string(),
            format!("{:.6e} s / {:.6e} s/s", gv.geom_delay_s, gv.geom_rate_sps),
        ),
        (
            "geom accel".to_string(),
            format!("{:.6e} s/s^2", gv.geom_accel_sps2),
        ),
        (
            "geom eval offset".to_string(),
            format!(
                "start {:+.3} s, read-align {:+.3} s",
                gv.geom_eval_offset_s, gv.read_align_ref_offset_s
            ),
        ),
        (
            "geom rate @obsfreq".to_string(),
            format!("{:.6} Hz", gv.geom_rate_hz_at_obs),
        ),
        (
            "clock relative delay/rate".to_string(),
            format!(
                "{:.6e} s / {:.6e} s/s",
                gv.rel_clock_delay_s, gv.rel_clock_rate_sps
            ),
        ),
        ("clock epoch mode".to_string(), gv.clock_epoch_mode.to_string()),
        (
            "coarse/res-delay".to_string(),
            format!(
                "{:.6e} s ({:.3} sample) / {:.3} sample",
                gv.coarse_delay_s, gv.coarse_delay_samples, gv.res_delay_input_samples
            ),
        ),
        (
            "read-align delay".to_string(),
            format!(
                "{:.3} sample ({:.6e} s)",
                gv.read_align_delay_samples, gv.read_align_delay_s
            ),
        ),
        (
            "delay model".to_string(),
            format!(
                "{}; model-time offset {:+.6} s",
                gv.delay_model, gv.model_time_offset_s
            ),
        ),
        (
            "earth orientation".to_string(),
            format!(
                "DUT1 {:+.6} s, TT-UTC {:+.6} s, xp {:+.6} arcsec, yp {:+.6} arcsec",
                gv.dut1_s, gv.tt_utc_s, gv.xp_arcsec, gv.yp_arcsec
            ),
        ),
        (
            "delay-rate terms [Hz]".to_string(),
            format!(
                "total {:.6} = geom {:.6} + clock {:.6} + rot-res {:.6} + user {:.6}",
                gv.delay_rate_hz,
                gv.delay_rate_geom_hz,
                gv.delay_rate_clock_hz,
                gv.delay_rate_rot_hz,
                gv.delay_rate_user_hz
            ),
        ),
        (
            "band-overlap".to_string(),
            format!(
                "{} bins ({:.3} MHz, {:.3} MHz/bin)",
                gv.band_overlap_bins, gv.band_overlap_mhz, gv.band_bin_mhz
            ),
        ),
        (
            "rotation-shift".to_string(),
            abbreviate_middle(&gv.rotation_shift_desc, 100),
        ),
        (
            "rotation-fringestop [Hz]".to_string(),
            format!(
                "delta {:.6}, shift {:.6}, residual {:.6}",
                gv.rotation_delta_hz, gv.rotation_shift_hz, gv.rotation_residual_hz
            ),
        ),
        (
            "band-align".to_string(),
            abbreviate_middle(&gv.band_align_desc, 120),
        ),
        ("output-grid".to_string(), gv.output_grid.to_string()),
        (
            format!("band-param {}", a1.name),
            format!(
                "key={} sideband={} obsfreq_mhz={:.9} rotation_mhz={:.9} ref_band_low_mhz={:.9} data_band_low_mhz={:.9} data_band_center_mhz={:.9} bw_mhz={:.9}",
                a1.key,
                a1.sideband,
                gv.obsfreq_mhz,
                a1.rotation_hz / 1.0e6,
                gv.obsfreq_mhz,
                gv.obsfreq_mhz + a1.rotation_hz / 1.0e6,
                gv.obsfreq_mhz + a1.rotation_hz / 1.0e6 + 0.5 * gv.bw_mhz,
                gv.bw_mhz
            ),
        ),
        (
            format!("band-param {}", a2.name),
            format!(
                "key={} sideband={} obsfreq_mhz={:.9} rotation_mhz={:.9} ref_band_low_mhz={:.9} data_band_low_mhz={:.9} data_band_center_mhz={:.9} bw_mhz={:.9}",
                a2.key,
                a2.sideband,
                gv.obsfreq_mhz,
                a2.rotation_hz / 1.0e6,
                gv.obsfreq_mhz,
                gv.obsfreq_mhz + a2.rotation_hz / 1.0e6,
                gv.obsfreq_mhz + a2.rotation_hz / 1.0e6 + 0.5 * gv.bw_mhz,
                gv.bw_mhz
            ),
        ),
    ];
    print_kv_table("Baseline Parameters", &baseline_rows);

    let ant1_label = format!("Ant1 [{}:{}]", a1.key, a1.name);
    let ant2_label = format!("Ant2 [{}:{}]", a2.key, a2.name);
    let ant_rows = vec![
        (
            "input raw".to_string(),
            path_file_display(a1.path),
            path_file_display(a2.path),
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
            format!(
                "[{:.3}, {:.3}, {:.3}]",
                a1.ecef_m[0], a1.ecef_m[1], a1.ecef_m[2]
            ),
            format!(
                "[{:.3}, {:.3}, {:.3}]",
                a2.ecef_m[0], a2.ecef_m[1], a2.ecef_m[2]
            ),
        ),
        (
            "bit / bit-code".to_string(),
            format!("{} / {}", a1.bit, format_bit_codes(a1.bit)),
            format!("{} / {}", a2.bit, format_bit_codes(a2.bit)),
        ),
        (
            "level / level-map".to_string(),
            abbreviate_middle(
                &format!("{:?} / {}", a1.levels, format_level_map(a1.bit, a1.levels)),
                72,
            ),
            abbreviate_middle(
                &format!("{:?} / {}", a2.levels, format_level_map(a2.bit, a2.levels)),
                72,
            ),
        ),
        (
            "shuffle-in".to_string(),
            abbreviate_middle(&format_shuffle_compact(a1.shuffle_ext), 56),
            abbreviate_middle(&format_shuffle_compact(a2.shuffle_ext), 56),
        ),
        (
            "sideband / rotation".to_string(),
            format!("{} / {}", a1.sideband, fmt_hz_to_mhz(a1.rotation_hz)),
            format!("{} / {}", a2.sideband, fmt_hz_to_mhz(a2.rotation_hz)),
        ),
        (
            "clock epoch".to_string(),
            a1.clock_epoch
                .map(|ep| {
                    format!(
                        "{} (dt to process {:+.3} s)",
                        ep,
                        a1.clock_epoch_dt_s.unwrap_or(0.0)
                    )
                })
                .unwrap_or_else(|| "-".to_string()),
            a2.clock_epoch
                .map(|ep| {
                    format!(
                        "{} (dt to process {:+.3} s)",
                        ep,
                        a2.clock_epoch_dt_s.unwrap_or(0.0)
                    )
                })
                .unwrap_or_else(|| "-".to_string()),
        ),
        (
            "clock XML delay/rate/accel/jerk/snap".to_string(),
            format!(
                "{:.6e} / {:.6e} / {:.6e} / {:.6e} / {:.6e}",
                a1.clock_delay_raw_s,
                a1.clock_rate_raw_sps,
                a1.clock_accel_raw_sps2,
                a1.clock_jerk_raw_sps3,
                a1.clock_snap_raw_sps4
            ),
            format!(
                "{:.6e} / {:.6e} / {:.6e} / {:.6e} / {:.6e}",
                a2.clock_delay_raw_s,
                a2.clock_rate_raw_sps,
                a2.clock_accel_raw_sps2,
                a2.clock_jerk_raw_sps3,
                a2.clock_snap_raw_sps4
            ),
        ),
        (
            "clock effective delay/rate/accel/jerk/snap".to_string(),
            format!(
                "{:.6e} / {:.6e} / {:.6e} / {:.6e} / {:.6e}",
                a1.clock_delay_s,
                a1.clock_rate_sps,
                a1.clock_accel_sps2,
                a1.clock_jerk_sps3,
                a1.clock_snap_sps4
            ),
            format!(
                "{:.6e} / {:.6e} / {:.6e} / {:.6e} / {:.6e}",
                a2.clock_delay_s,
                a2.clock_rate_sps,
                a2.clock_accel_sps2,
                a2.clock_jerk_sps3,
                a2.clock_snap_sps4
            ),
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
        RunMode::Corr => {
            if args.schedule.is_none() {
                return Err("--schedule is required for yi-corr".into());
            }
            if args.raw_directory.is_none() {
                return Err("--raw-directory is required for yi-corr".into());
            }
            if args.cor_directory.is_none() {
                return Err("--cor-directory is required for yi-corr".into());
            }
        }
    }
    if args.model_sweep {
        model_sweep::ensure_not_recursive_model_sweep(true)?;
        if !matches!(run_mode, RunMode::Corr) {
            return Err("--model-sweep is supported only by yi-corr".into());
        }
        return model_sweep::run_unattended_model_sweep(&args);
    }
    init_stdout_log_for_yi_corr(run_mode, args.stdout)?;
    let physical_cores = runtime_physical_cores();
    let cpu_auto = physical_cores
        .map(|n| n.saturating_sub(1).max(1))
        .unwrap_or_else(|| {
            build_logical_cpus()
                .or_else(|| std::thread::available_parallelism().ok().map(|n| n.get()))
                .unwrap_or(2)
                .saturating_sub(2)
                .max(1)
        });
    let cpu_threads = match args.cpu {
        Some(0) => return Err("--cpu must be >= 1".into()),
        Some(v) => v,
        None => cpu_auto,
    };
    if args.cpu.is_none() {
        if let Some(physical) = physical_cores {
            println!(
                "[info] CPU auto-tuning: physical-cores={} reader-reserve=1 compute-threads={}",
                physical, cpu_threads
            );
        }
    }
    let affinity_runtime = affinity::AffinityRuntime::from_default_file()?;
    if let Some(msg) = affinity_runtime.info() {
        println!("[info] {msg}");
    }
    let reader_core = affinity::reader_core_from_env()?;
    if let Some(core) = reader_core {
        println!(
            "[info] I/O reader thread affinity enabled via YI_READER_CORE: core={}",
            core.id
        );
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
                let run_result = run_once(run_args, run_mode, cpu_threads, reader_core);
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
    run_once(args, run_mode, cpu_threads, reader_core)
}

fn run_once(
    args: args::Args,
    run_mode: RunMode,
    cpu_threads: usize,
    reader_core: Option<core_affinity::CoreId>,
) -> Result<(), DynError> {
    let fringe_interval_s = args.fringe;
    if let Some(v) = fringe_interval_s {
        if !v.is_finite() || v <= 0.0 {
            return Err("--fringe interval must be a positive finite number of seconds".into());
        }
        if !matches!(run_mode, RunMode::Corr) {
            return Err("--fringe requires yi-corr because it uses XCF spectra".into());
        }
    }
    let do_xcf = false;
    let do_synth = true;
    let write_raw = matches!(run_mode, RunMode::PhasedArray);
    let write_acf_cor = matches!(run_mode, RunMode::Corr);
    let write_xcf_cor = matches!(run_mode, RunMode::Corr);
    let write_phased_cor = matches!(run_mode, RunMode::PhasedArray) && args.phased_diagnostics;
    let plot_phased = matches!(run_mode, RunMode::PhasedArray) && args.phased_diagnostics;
    let if_d = if let Some(p) = &args.schedule {
        let is_xml = p
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("xml"))
            .unwrap_or(false);
        if !is_xml {
            return Err("--schedule must point to a .xml file".into());
        }
        Some(Arc::new(ifile::parse_ifile_for_process(
            p,
            args.process_index,
        )?))
    } else {
        None
    };
    let selected_process = if let Some(meta) = if_d.as_ref() {
        if let Some(idx) = args.process_index {
            Some(meta.processes.get(idx).cloned().ok_or_else(|| {
                format!(
                    "--process-index {} out of range (0..{})",
                    idx,
                    meta.processes.len().saturating_sub(1)
                )
            })?)
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
    let fft_len = if args.fft == args::DEFAULT_FFT {
        if_d.as_ref().and_then(|d| d.fft).unwrap_or(args.fft)
    } else {
        args.fft
    };

    let bit_a = args.bin.clone();
    let bit1_def = if_d.as_ref().and_then(|d| d.ant1_bit).unwrap_or(2);
    let bit2_def = if_d.as_ref().and_then(|d| d.ant2_bit).unwrap_or(2);
    let (bit1, bit2): (usize, usize) =
        resolve_per_antenna_config_with_defaults(&bit_a, bit1_def, bit2_def, |s: &str| {
            Ok(s.parse()?)
        })?;
    let level_src = args.level.clone();
    let level_args = normalize_level_args(&level_src, bit1, bit2)?;
    let lv1_def = if_d
        .as_ref()
        .and_then(|d| d.ant1_level.clone())
        .unwrap_or("-1.5,-0.5,0.5,1.5".into());
    let lv2_def = if_d
        .as_ref()
        .and_then(|d| d.ant2_level.clone())
        .unwrap_or("-1.5,-0.5,0.5,1.5".into());
    let (lv1_s, lv2_s): (String, String) =
        resolve_per_antenna_config_with_defaults(&level_args, lv1_def, lv2_def, |s: &str| {
            Ok(s.to_string())
        })?;
    let (levels1, levels2) = (
        Arc::new(parse_levels(bit1, &lv1_s)?),
        Arc::new(parse_levels(bit2, &lv2_s)?),
    );
    let level_power1 = quantization_mean_power(levels1.as_ref(), "ant1")?;
    let level_power2 = quantization_mean_power(levels2.as_ref(), "ant2")?;
    let ds_s = DEFAULT_SHUFFLE_IN
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let shuffle_src = args.shuffle_in.clone();
    let shuffle_args = normalize_shuffle_args(&shuffle_src)?;
    let sh1_def = if_d
        .as_ref()
        .and_then(|d| d.ant1_shuffle.clone())
        .unwrap_or(ds_s.clone());
    let sh2_def = if_d
        .as_ref()
        .and_then(|d| d.ant2_shuffle.clone())
        .unwrap_or(ds_s);
    let (sh1_s, sh2_s): (String, String) =
        resolve_per_antenna_config_with_defaults(&shuffle_args, sh1_def, sh2_def, |s: &str| {
            Ok(s.to_string())
        })?;
    let (sh1, sh2) = (
        Arc::new(parse_shuffle(&sh1_s)?),
        Arc::new(parse_shuffle(&sh2_s)?),
    );
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

    let sb_a = args.sideband.clone();
    let sb1_def = if_d
        .as_ref()
        .and_then(|d| d.ant1_sideband.clone())
        .unwrap_or("USB".into());
    let sb2_def = if_d
        .as_ref()
        .and_then(|d| d.ant2_sideband.clone())
        .unwrap_or("USB".into());
    let (sb1_s, sb2_s): (String, String) =
        resolve_per_antenna_config_with_defaults(&sb_a, sb1_def, sb2_def, |s: &str| {
            Ok(s.to_uppercase())
        })?;
    let (lsb1_raw, lsb2_raw) = (sb1_s.as_str() == "LSB", sb2_s.as_str() == "LSB");
    let lsb1 = ant_lsb_override(lsb1_raw, "YI_ANT1_LSB_TO_USB");
    let lsb2 = ant_lsb_override(lsb2_raw, "YI_ANT2_LSB_TO_USB");
    let (tsys1, tsys2) = resolve_per_antenna_config(&args.tsys, 1.0, |s| Ok(s.parse()?))?;
    let (dia1, dia2) = resolve_per_antenna_config(&args.diameter, 0.0, |s| Ok(s.parse()?))?;
    let (eta1, eta2) = resolve_per_antenna_config(&args.eta, 0.65, |s| Ok(s.parse()?))?;
    let (gain1, gain2) = resolve_per_antenna_config(&args.gain, 1.0, |s| Ok(s.parse()?))?;
    let (sefd1, sefd2) = resolve_per_antenna_config(&args.sefd, 0.0, |s| Ok(s.parse()?))?;

    let fs = if_d
        .as_ref()
        .and_then(|d| d.sampling_hz)
        .unwrap_or(args.sampling * 1e6);
    let obs_mhz = if_d
        .as_ref()
        .and_then(|d| d.obsfreq_mhz)
        .unwrap_or(args.obsfreq);
    let ep_i = args
        .epoch
        .clone()
        .or_else(|| selected_process.as_ref().map(|p| p.epoch.clone()))
        .or_else(|| if_d.as_ref().and_then(|d| d.epoch.clone()))
        .unwrap_or("2000".into());
    let (a1p, a2p, _, _unused_tag) =
        resolve_input_paths(&args, &ep_i, if_d.as_ref().map(|v| &**v))?;
    let (c_unix_base, process_epoch_tag) = epoch_to_yyyydddhhmmss(&ep_i)?;
    let ant1_ecef = if_d
        .as_ref()
        .and_then(|d| d.ant1_ecef_m)
        .unwrap_or(geom::YAMAGU32_ECEF);
    let ant2_ecef = if_d
        .as_ref()
        .and_then(|d| d.ant2_ecef_m)
        .unwrap_or(geom::YAMAGU34_ECEF);

    let schedule_param_source_xml = if_d.is_some();
    let model_time_offset_s = if schedule_param_source_xml {
        if_d.as_ref()
            .and_then(|d| d.model_time_offset_s)
            .unwrap_or(args::DEFAULT_MODEL_TIME_OFFSET_S)
    } else {
        args.model_time_offset
    };
    let dut1_s = if schedule_param_source_xml {
        if_d.as_ref().and_then(|d| d.dut1_s).unwrap_or(0.0)
    } else {
        args.dut1
    };
    let tt_utc_s = if schedule_param_source_xml {
        if_d.as_ref().and_then(|d| d.tt_utc_s).unwrap_or(69.184)
    } else {
        args.tt_utc
    };
    let xp_arcsec = if schedule_param_source_xml {
        if_d.as_ref().and_then(|d| d.xp_arcsec).unwrap_or(0.0)
    } else {
        args.xp
    };
    let yp_arcsec = if schedule_param_source_xml {
        if_d.as_ref().and_then(|d| d.yp_arcsec).unwrap_or(0.0)
    } else {
        args.yp
    };

    let explicit_xml_eop = if schedule_param_source_xml {
        if_d.as_ref().is_some_and(|d| {
            d.dut1_s.is_some()
                || d.tt_utc_s.is_some()
                || d.xp_arcsec.is_some()
                || d.yp_arcsec.is_some()
        })
    } else {
        false
    };
    let explicit_env_eop = ["YI_DUT1_S", "YI_TT_UTC_S", "YI_XP_ARCSEC", "YI_YP_ARCSEC"]
        .iter()
        .all(|name| std::env::var_os(name).is_some());
    let eop_mode = std::env::var("YI_EOP_MODE")
        .unwrap_or_else(|_| "auto".to_string())
        .to_ascii_lowercase();
    if !matches!(eop_mode.as_str(), "auto" | "strict" | "none") {
        return Err(format!("invalid YI_EOP_MODE='{eop_mode}' (use auto, strict, or none)").into());
    }
    let mut dut1_s = dut1_s;
    let mut tt_utc_s = tt_utc_s;
    let mut xp_arcsec = xp_arcsec;
    let mut yp_arcsec = yp_arcsec;
    if !explicit_xml_eop && !explicit_env_eop && eop_mode != "none" {
        let eop_path = if_d
            .as_ref()
            .and_then(|d| d.eop_file.clone())
            .or_else(|| std::env::var_os("YI_EOP_FILE").map(PathBuf::from))
            .or_else(|| {
                let default = PathBuf::from("data/eop/finals2000A.data");
                default.exists().then_some(default)
            });
        if let Some(path) = eop_path {
            match eop::interpolate_finals2000a(&path, geom::parse_epoch_to_mjd(&ep_i)?) {
                Ok(v) => {
                    dut1_s = v.dut1_s;
                    tt_utc_s = v.tt_utc_s;
                    xp_arcsec = v.xp_arcsec;
                    yp_arcsec = v.yp_arcsec;
                    println!(
                        "[info] EOP file: {} MJD {:.1}..{:.1} {}..{} -> DUT1={:+.6}s TT-UTC={:+.6}s xp={:+.6}arcsec yp={:+.6}arcsec",
                        v.source,
                        v.mjd0,
                        v.mjd1,
                        v.kind0.as_str(),
                        v.kind1.as_str(),
                        dut1_s,
                        tt_utc_s,
                        xp_arcsec,
                        yp_arcsec
                    );
                }
                Err(e) if eop_mode == "strict" => return Err(e),
                Err(e) => {
                    println!(
                        "[warn] EOP file unavailable for this epoch ({e}); using fallback DUT1/TT-UTC/xp/yp"
                    );
                }
            }
        } else if eop_mode == "strict" {
            return Err("YI_EOP_MODE=strict requires XML <eop file=...>, YI_EOP_FILE, or data/eop/finals2000A.data".into());
        } else {
            println!(
                "[warn] no EOP file specified/found; using fallback DUT1/TT-UTC/xp/yp (YI_EOP_MODE=auto)"
            );
        }
    }

    let dut1_s = env_f64_override("YI_DUT1_S", dut1_s)?;
    let tt_utc_s = env_f64_override("YI_TT_UTC_S", tt_utc_s)?;
    let xp_arcsec = env_f64_override("YI_XP_ARCSEC", xp_arcsec)?;
    let yp_arcsec = env_f64_override("YI_YP_ARCSEC", yp_arcsec)?;

    let mut gdi: Option<GDI> = None;
    struct GDI {
        ra: f64,
        dec: f64,
        ra_raw: f64,
        dec_raw: f64,
        ra_header: f64,
        dec_header: f64,
        mjd: f64,
    }
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
    let earth_orientation = geom::EarthOrientation {
        dut1_s,
        tt_minus_utc_s: tt_utc_s,
        xp_arcsec,
        yp_arcsec,
    };
    let geom_delay_mode = match std::env::var("YI_GEOM_DELAY_MODE") {
        Ok(v) => match v.trim().to_ascii_lowercase().as_str() {
            "" | "anchored" | "mixed" | "default" => geom::GeometricDelayMode::Anchored,
            "bary" | "barycentric" | "absolute-barycentric" => {
                geom::GeometricDelayMode::Barycentric
            }
            "vlbi-minus" | "barycentric-minus" | "first-order-minus" => {
                geom::GeometricDelayMode::VlbiMinus
            }
            "vlbi-plus" | "barycentric-plus" | "first-order-plus" | "vlbi" => {
                geom::GeometricDelayMode::VlbiPlus
            }
            "geo" | "geocentric" | "absolute-geocentric" => geom::GeometricDelayMode::Geocentric,
            other => {
                return Err(format!(
                    "invalid YI_GEOM_DELAY_MODE='{other}' (use anchored, barycentric, vlbi-minus, vlbi-plus, or geocentric)"
                )
                .into())
            }
        },
        Err(_) => geom::GeometricDelayMode::Anchored,
    };
    let source_vector_mode = match std::env::var("YI_SOURCE_VECTOR_MODE") {
        Ok(v) => match v.trim().to_ascii_lowercase().as_str() {
            "" | "mean" | "mean-gast" | "precess-gast" | "default" => {
                geom::SourceVectorMode::MeanGast
            }
            "pnm" | "pnm-gast" | "true-gast" => geom::SourceVectorMode::PnmGast,
            "pnm-era" | "era" | "c2t" | "c2t-era" => geom::SourceVectorMode::PnmEra,
            other => {
                return Err(format!(
                    "invalid YI_SOURCE_VECTOR_MODE='{other}' (use mean-gast, pnm-gast, or pnm-era)"
                )
                .into())
            }
        },
        Err(_) => geom::SourceVectorMode::MeanGast,
    };
    if let (Some(ra_s), Some(dec_s)) = (ra_in, dec_in) {
        let (ra_raw, dec_raw, mjd) = (
            geom::parse_ra(&ra_s)?,
            geom::parse_dec(&dec_s)?,
            geom::parse_epoch_to_mjd(&ep_i)?,
        );
        // Input sky coordinates are J2000; precess to date for the delay model.
        // Keep the original J2000 coordinates for .cor headers separately.
        let mjd_tt = mjd + tt_utc_s / 86400.0;
        let (ra, dec) = geom::precess_j2000_to_mean_of_date(ra_raw, dec_raw, mjd_tt);
        let (ra_model, dec_model) = match source_vector_mode {
            geom::SourceVectorMode::MeanGast => (ra, dec),
            geom::SourceVectorMode::PnmGast | geom::SourceVectorMode::PnmEra => (ra_raw, dec_raw),
        };
        let (_, _, gd, gr, ga) = geom::calculate_geometric_delay_and_derivatives_full_with_eop(
            ant1_ecef,
            ant2_ecef,
            ra_model,
            dec_model,
            mjd,
            mjd,
            earth_orientation,
            geom_delay_mode,
            source_vector_mode,
        );
        gdi = Some(GDI {
            ra,
            dec,
            ra_raw,
            dec_raw,
            ra_header: ra_raw,
            dec_header: dec_raw,
            mjd,
        });
        (gd0, gr0, ga0) = (gd, gr, ga);
    }
    let clock_delay_rel_legacy_s = if_d.as_ref().and_then(|d| d.clock_delay_s).unwrap_or(0.0);
    let clock_rate_rel_legacy_sps = if_d.as_ref().and_then(|d| d.clock_rate_sps).unwrap_or(0.0);
    let clock_accel_rel_legacy_sps2 = if_d
        .as_ref()
        .and_then(|d| d.clock_accel_sps2)
        .unwrap_or(0.0);
    let clock_jerk_rel_legacy_sps3 = if_d.as_ref().and_then(|d| d.clock_jerk_sps3).unwrap_or(0.0);
    let clock_snap_rel_legacy_sps4 = if_d.as_ref().and_then(|d| d.clock_snap_sps4).unwrap_or(0.0);
    let clock1_delay_base_s = if_d
        .as_ref()
        .and_then(|d| d.ant1_clock_delay_s)
        .unwrap_or(0.0);
    let clock2_delay_base_s = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_delay_s)
        .unwrap_or(clock1_delay_base_s + clock_delay_rel_legacy_s);
    let clock1_rate_sps = if_d
        .as_ref()
        .and_then(|d| d.ant1_clock_rate_sps)
        .unwrap_or(0.0);
    let clock2_rate_sps = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_rate_sps)
        .unwrap_or(clock1_rate_sps + clock_rate_rel_legacy_sps);
    let clock1_rate_raw_sps = clock1_rate_sps;
    let clock2_rate_raw_sps = clock2_rate_sps;
    let clock1_accel_sps2_from_meta = if_d
        .as_ref()
        .and_then(|d| d.ant1_clock_accel_sps2)
        .unwrap_or(0.0);
    let clock2_accel_sps2_from_meta = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_accel_sps2)
        .unwrap_or(clock1_accel_sps2_from_meta + clock_accel_rel_legacy_sps2);
    let clock1_jerk_sps3_from_meta = if_d
        .as_ref()
        .and_then(|d| d.ant1_clock_jerk_sps3)
        .unwrap_or(0.0);
    let clock2_jerk_sps3_from_meta = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_jerk_sps3)
        .unwrap_or(clock1_jerk_sps3_from_meta + clock_jerk_rel_legacy_sps3);
    let clock1_snap_sps4_from_meta = if_d
        .as_ref()
        .and_then(|d| d.ant1_clock_snap_sps4)
        .unwrap_or(0.0);
    let clock2_snap_sps4_from_meta = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_snap_sps4)
        .unwrap_or(clock1_snap_sps4_from_meta + clock_snap_rel_legacy_sps4);
    let clock1_accel_raw_sps2 = clock1_accel_sps2_from_meta;
    let clock2_accel_raw_sps2 = clock2_accel_sps2_from_meta;
    let clock1_jerk_raw_sps3 = clock1_jerk_sps3_from_meta;
    let clock2_jerk_raw_sps3 = clock2_jerk_sps3_from_meta;
    let clock1_snap_raw_sps4 = clock1_snap_sps4_from_meta;
    let clock2_snap_raw_sps4 = clock2_snap_sps4_from_meta;
    let clock1_epoch = if_d
        .as_ref()
        .and_then(|d| d.ant1_clock_epoch.as_deref())
        .map(|s| s.to_string());
    let clock2_epoch = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_epoch.as_deref())
        .map(|s| s.to_string());
    let parse_clock_epoch_offset = |station: &str,
                                    ep: Option<&str>|
     -> Result<Option<f64>, DynError> {
        let Some(ep) = ep else {
            return Ok(None);
        };
        let ep_trim = ep.trim();
        match ep_trim.to_ascii_lowercase().as_str() {
            "" | "process" | "process-epoch" | "scan" | "scan-epoch" | "current" => Ok(Some(0.0)),
            _ => {
                let (clock_unix, _) = epoch_to_yyyydddhhmmss(ep_trim)
                    .map_err(|e| format!("invalid {station} clock epoch '{}': {}", ep_trim, e))?;
                Ok(Some((c_unix_base - clock_unix) as f64))
            }
        }
    };
    let clock1_epoch_offset_s = parse_clock_epoch_offset("ant1", clock1_epoch.as_deref())?;
    let clock2_epoch_offset_s = parse_clock_epoch_offset("ant2", clock2_epoch.as_deref())?;
    let clock_epoch_mode = match std::env::var("YI_CLOCK_EPOCH_MODE") {
        Ok(v) => match v.trim().to_ascii_lowercase().as_str() {
            "" | "process" | "process-epoch" | "frinz" | "residual" => "process",
            "clock" | "clock-epoch" | "propagate" | "absolute" => "clock",
            other => {
                return Err(
                    format!("invalid YI_CLOCK_EPOCH_MODE='{other}' (use process or clock)").into(),
                )
            }
        },
        Err(_) => "clock",
    };
    let propagate_clock_epoch = clock_epoch_mode == "clock";

    // XML clock terms are defined at their own clock epoch.  Propagate them to
    // the process epoch before per-frame evaluation.  Use
    // YI_CLOCK_EPOCH_MODE=process only for empirical frinZ residual corrections
    // that have already been solved at the process epoch.
    let clock1_dt0 = if propagate_clock_epoch {
        clock1_epoch_offset_s.unwrap_or(0.0)
    } else {
        0.0
    };
    let clock2_dt0 = if propagate_clock_epoch {
        clock2_epoch_offset_s.unwrap_or(0.0)
    } else {
        0.0
    };
    let clock1_dt0_2 = clock1_dt0 * clock1_dt0;
    let clock1_dt0_3 = clock1_dt0_2 * clock1_dt0;
    let clock1_dt0_4 = clock1_dt0_2 * clock1_dt0_2;
    let clock2_dt0_2 = clock2_dt0 * clock2_dt0;
    let clock2_dt0_3 = clock2_dt0_2 * clock2_dt0;
    let clock2_dt0_4 = clock2_dt0_2 * clock2_dt0_2;
    let clock1_delay_s = clock1_delay_base_s
        + clock1_rate_sps * clock1_dt0
        + 0.5 * clock1_accel_sps2_from_meta * clock1_dt0_2
        + (1.0 / 6.0) * clock1_jerk_sps3_from_meta * clock1_dt0_3
        + (1.0 / 24.0) * clock1_snap_sps4_from_meta * clock1_dt0_4;
    let clock2_delay_s = clock2_delay_base_s
        + clock2_rate_sps * clock2_dt0
        + 0.5 * clock2_accel_sps2_from_meta * clock2_dt0_2
        + (1.0 / 6.0) * clock2_jerk_sps3_from_meta * clock2_dt0_3
        + (1.0 / 24.0) * clock2_snap_sps4_from_meta * clock2_dt0_4;
    let clock1_rate_sps = clock1_rate_sps
        + clock1_accel_sps2_from_meta * clock1_dt0
        + 0.5 * clock1_jerk_sps3_from_meta * clock1_dt0_2
        + (1.0 / 6.0) * clock1_snap_sps4_from_meta * clock1_dt0_3;
    let clock2_rate_sps = clock2_rate_sps
        + clock2_accel_sps2_from_meta * clock2_dt0
        + 0.5 * clock2_jerk_sps3_from_meta * clock2_dt0_2
        + (1.0 / 6.0) * clock2_snap_sps4_from_meta * clock2_dt0_3;
    let clock1_accel_sps2_from_meta = clock1_accel_sps2_from_meta
        + clock1_jerk_sps3_from_meta * clock1_dt0
        + 0.5 * clock1_snap_sps4_from_meta * clock1_dt0_2;
    let clock2_accel_sps2_from_meta = clock2_accel_sps2_from_meta
        + clock2_jerk_sps3_from_meta * clock2_dt0
        + 0.5 * clock2_snap_sps4_from_meta * clock2_dt0_2;
    let clock1_jerk_sps3_from_meta =
        clock1_jerk_sps3_from_meta + clock1_snap_sps4_from_meta * clock1_dt0;
    let clock2_jerk_sps3_from_meta =
        clock2_jerk_sps3_from_meta + clock2_snap_sps4_from_meta * clock2_dt0;
    // Optional relative clock polynomial offsets for engineering tests.
    // Production clock acceleration should be supplied in XML / ifile metadata.
    // These offsets are applied with t=0 at the process epoch and
    // t=skip+local_time in per-frame evaluation.
    let parse_env_f64 = |name: &str, default: f64| -> Result<f64, DynError> {
        match std::env::var(name) {
            Ok(v) => v
                .trim()
                .parse::<f64>()
                .map_err(|e| format!("invalid {name}='{v}': {e}").into()),
            Err(_) => Ok(default),
        }
    };
    let rel_clock_delay_offset_s = parse_env_f64("YI_REL_CLOCK_DELAY_S", 0.0)?;
    let rel_clock_rate_offset_sps = parse_env_f64("YI_REL_CLOCK_RATE_SPS", 0.0)?;
    let rel_clock_accel_sps2 = parse_env_f64("YI_REL_CLOCK_ACCEL_SPS2", 0.0)?;
    let rel_clock_jerk_sps3 = parse_env_f64("YI_REL_CLOCK_JERK_SPS3", 0.0)?;
    let rel_clock_snap_sps4 = parse_env_f64("YI_REL_CLOCK_SNAP_SPS4", 0.0)?;
    let clock1_accel_sps2 =
        clock1_accel_sps2_from_meta + parse_env_f64("YI_CLOCK1_ACCEL_SPS2", 0.0)?;
    let clock2_accel_sps2 = clock2_accel_sps2_from_meta
        + parse_env_f64("YI_CLOCK2_ACCEL_SPS2", 0.0)?
        + rel_clock_accel_sps2;
    let clock1_jerk_sps3 = clock1_jerk_sps3_from_meta + parse_env_f64("YI_CLOCK1_JERK_SPS3", 0.0)?;
    let clock2_jerk_sps3 = clock2_jerk_sps3_from_meta
        + parse_env_f64("YI_CLOCK2_JERK_SPS3", 0.0)?
        + rel_clock_jerk_sps3;
    let clock1_snap_sps4 = clock1_snap_sps4_from_meta + parse_env_f64("YI_CLOCK1_SNAP_SPS4", 0.0)?;
    let clock2_snap_sps4 = clock2_snap_sps4_from_meta
        + parse_env_f64("YI_CLOCK2_SNAP_SPS4", 0.0)?
        + rel_clock_snap_sps4;
    let clock1_delay_s = clock1_delay_s;
    let clock2_delay_s = clock2_delay_s + rel_clock_delay_offset_s;
    let clock1_rate_sps = clock1_rate_sps;
    let clock2_rate_sps = clock2_rate_sps + rel_clock_rate_offset_sps;
    let clock_delay_s = clock2_delay_s - clock1_delay_s;
    let clock_rate_sps = clock2_rate_sps - clock1_rate_sps;
    let clock_accel_sps2 = clock2_accel_sps2 - clock1_accel_sps2;
    let clock_jerk_sps3 = clock2_jerk_sps3 - clock1_jerk_sps3;
    let clock_snap_sps4 = clock2_snap_sps4 - clock1_snap_sps4;
    let coarse_delay_s = args.coarse.unwrap_or(DEFAULT_COARSE_DELAY_S);
    let delay_user_samples = args.delay + args.resdelay;
    let rate_user_hz = args.rate + args.resrate;
    let accel_user_hzps = args.resacel;
    let accel_user_sps2 = accel_user_hzps / (obs_mhz * 1e6);
    let rot1_regular = if_d
        .as_ref()
        .and_then(|d| d.ant1_rotation_hz)
        .unwrap_or(0.0);
    let rot2_regular = if_d
        .as_ref()
        .and_then(|d| d.ant2_rotation_hz)
        .unwrap_or(0.0);
    let rot1_hermitian = if_d.as_ref().and_then(|d| d.ant1_rotation2_hz);
    let rot2_hermitian = if_d.as_ref().and_then(|d| d.ant2_rotation2_hz);
    let rotation_hermitian_mode = rot1_hermitian.is_some() || rot2_hermitian.is_some();
    let mut rot1 = rot1_hermitian.unwrap_or(rot1_regular);
    let mut rot2 = rot2_hermitian.unwrap_or(rot2_regular);
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
    let a1_data_low = obs_mhz + rot1 / 1e6;
    let a2_data_low = obs_mhz + rot2 / 1e6;
    let raw_ba = compute_band_alignment(
        fft_len,
        fs,
        a1_data_low + 0.5 * bw,
        a2_data_low + 0.5 * bw,
        bw,
        bw,
    )?;
    let rotation_bins1 = (rot1 / (fs / fft_len as f64)).round() as isize;
    let rotation_bins2 = (rot2 / (fs / fft_len as f64)).round() as isize;
    let a1_corr_low = a1_data_low - rot1 / 1e6;
    let a2_corr_low = a2_data_low - rot2 / 1e6;
    let fringe_stop_freq_mode = FringeStopFreqMode::from_env()?;
    let lo1_mhz = fringe_stop_freq_mode.carrier_mhz(a1_data_low, bw, obs_mhz, lsb1_raw);
    let lo2_mhz = fringe_stop_freq_mode.carrier_mhz(a2_data_low, bw, obs_mhz, lsb2_raw);
    let lo1_hz = lo1_mhz * 1e6;
    let lo2_hz = lo2_mhz * 1e6;
    let ba = if rotation_hermitian_mode {
        compute_band_alignment(
            fft_len,
            fs,
            a1_corr_low + 0.5 * bw,
            a2_corr_low + 0.5 * bw,
            bw,
            bw,
        )?
    } else {
        let half_bins = fft_len / 2;
        let valid_start1 = rotation_bins1.max(0) as usize;
        let valid_start2 = rotation_bins2.max(0) as usize;
        let valid_end1 =
            (rotation_bins1 + half_bins as isize).clamp(0, half_bins as isize) as usize;
        let valid_end2 =
            (rotation_bins2 + half_bins as isize).clamp(0, half_bins as isize) as usize;
        let start = valid_start1.max(valid_start2);
        let end = valid_end1.min(valid_end2);
        if end <= start {
            return Err("No real signal band overlap after rotation".into());
        }
        BandAlignment {
            shift_bins: raw_ba.shift_bins,
            a1s: start,
            a1e: end,
            a2s: start,
        }
    };
    let out_grid = if rot1.abs() <= rot2.abs() {
        OutputGrid::Ant1
    } else {
        OutputGrid::Ant2
    };
    // The phased raw is returned to the native packed-data representation of
    // the selected output-grid antenna. Bit depth, levels, shuffle, sideband,
    // and rotation are one format bundle and must not be mixed between
    // stations.
    let (
        bit_out,
        levels_out,
        level_power_out,
        output_shuffle,
        output_shuffle_ext,
        output_sideband,
        output_lsb_restore,
        output_rotation_hz,
    ) = match out_grid {
        OutputGrid::Ant1 => (
            bit1,
            Arc::clone(&levels1),
            level_power1,
            sh1.as_ref().clone(),
            sh1_ext.clone(),
            sb1_s.clone(),
            lsb1,
            rot1,
        ),
        OutputGrid::Ant2 => (
            bit2,
            Arc::clone(&levels2),
            level_power2,
            sh2.as_ref().clone(),
            sh2_ext.clone(),
            sb2_s.clone(),
            lsb2,
            rot2,
        ),
    };

    let bpf1 = (fft_len * bit1 + 7) / 8;
    let bpf2 = (fft_len * bit2 + 7) / 8;
    let bpf_o = (fft_len * bit_out + 7) / 8;
    let bytes_per_frame_pair = bpf1 + bpf2;
    let io_chunk_frames = if args.razoku5bay {
        let forced = ((0.5_f64 * fs) / fft_len as f64).round().max(1.0) as usize;
        if args.chunk_frames.is_some() {
            println!(
                "[info] --razoku5bay is enabled: overriding --chunk-frames with {}",
                forced
            );
        } else {
            println!(
                "[info] --razoku5bay is enabled: forcing chunk size to {} frames (~0.5 s)",
                forced
            );
        }
        forced
    } else {
        match args.chunk_frames {
            Some(0) => return Err("--chunk-frames must be >= 1".into()),
            Some(v) => v,
            None => auto_chunk_frames(cpu_threads, bytes_per_frame_pair),
        }
    };
    let io_pipeline_depth = match args.pipeline_depth {
        Some(0) => return Err("--pipeline-depth must be >= 1".into()),
        Some(v) => v,
        None => auto_pipeline_depth(cpu_threads),
    };
    let f1_m = std::fs::metadata(&a1p)?;
    let f2_m = std::fs::metadata(&a2p)?;

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

    let eval_geom_at_elapsed = |elapsed_s: f64| -> Result<(f64, f64, f64), DynError> {
        if let Some(v) = gdi.as_ref() {
            let mjd_t = v.mjd + elapsed_s / 86400.0;
            let mjd_tt_t = mjd_t + tt_utc_s / 86400.0;
            let (ra_t, dec_t) = match source_vector_mode {
                geom::SourceVectorMode::MeanGast => {
                    geom::precess_j2000_to_mean_of_date(v.ra_raw, v.dec_raw, mjd_tt_t)
                }
                geom::SourceVectorMode::PnmGast | geom::SourceVectorMode::PnmEra => {
                    (v.ra_raw, v.dec_raw)
                }
            };
            let (_, _, gd_t, gr_t, ga_t) =
                geom::calculate_geometric_delay_and_derivatives_full_with_eop(
                    ant1_ecef,
                    ant2_ecef,
                    ra_t,
                    dec_t,
                    mjd_t,
                    v.mjd,
                    earth_orientation,
                    geom_delay_mode,
                    source_vector_mode,
                );
            Ok((gd_t, gr_t, ga_t))
        } else {
            Ok((
                gd0 + gr0 * elapsed_s + 0.5 * ga0 * elapsed_s * elapsed_s,
                gr0 + ga0 * elapsed_s,
                ga0,
            ))
        }
    };
    let (geom_delay_at_start_s, geom_rate_at_start_sps, geom_accel_at_start_sps2) =
        eval_geom_at_elapsed(total_skip_sec)?;
    let net_d_rel_no_clock0 = geom_delay_at_start_s + coarse_delay_s + delay_user_samples / fs;
    let net_d0 = net_d_rel_no_clock0 + clock_delay_s;

    let total_samples1 = (f1_m.len() * 8) / bit1 as u64;
    let total_samples2 = (f2_m.len() * 8) / bit2 as u64;

    // Choose a fixed read-align delay for the whole processed window.
    //
    // Older paths used round(net_d0*fs), i.e. the delay at the process start.
    // That removed sector-to-sector read-start jumps, but it could leave the
    // whole scan on the neighboring one-sample delay branch when the model delay
    // drifted across a half-sample point shortly after the start.  This was seen
    // as yi-corr residual delays being about one sample larger than KLH/frinZ4
    // on the 2025/302 08:15 NRAO530 test.
    //
    // Use the midpoint of the requested processing window as the read-align
    // reference instead.  The read start is still fixed for the whole process,
    // so the 0.5-sample sector-boundary problem does not come back, but the
    // chosen branch minimizes the residual delay over the scan rather than only
    // at its first frame.
    let requested_sec_for_align = process_window_sec
        .map(|window_sec| (window_sec - total_skip_sec).max(0.0))
        .unwrap_or(0.0);
    let fx_start_cumulative_seek = std::env::var("YI_FX_START_SEEK")
        .map(|v| {
            let v = v.trim().to_ascii_lowercase();
            !(v.is_empty() || v == "0" || v == "false" || v == "off" || v == "no")
        })
        .unwrap_or(false);
    let read_align_ref_local_s = if fx_start_cumulative_seek {
        total_skip_sec
    } else {
        total_skip_sec + 0.5 * requested_sec_for_align
    };
    let read_align_ref_label = if fx_start_cumulative_seek {
        "start"
    } else {
        "midpoint"
    };
    if fx_start_cumulative_seek {
        println!("[info] FX integer seek mode: start-based cumulative sample seek enabled (YI_FX_START_SEEK=1)");
    }
    let df_hz_for_align = fs / fft_len as f64;
    let rotation_shift_hz_for_align = raw_ba.shift_bins as f64 * df_hz_for_align;
    let rotation_fringe_hz_for_align = -((rot1 - rot2) - rotation_shift_hz_for_align);
    let extra_delay_rate_sps_for_align =
        (rotation_fringe_hz_for_align + rate_user_hz) / (obs_mhz * 1e6);
    let (geom_delay_at_align_s, _geom_rate_at_align_sps, _geom_accel_at_align_sps2) =
        eval_geom_at_elapsed(read_align_ref_local_s)?;
    let align_t2 = read_align_ref_local_s * read_align_ref_local_s;
    let align_t3 = align_t2 * read_align_ref_local_s;
    let align_t4 = align_t2 * align_t2;
    let clock1_align_s = clock1_delay_s
        + clock1_rate_sps * read_align_ref_local_s
        + 0.5 * clock1_accel_sps2 * align_t2
        + (1.0 / 6.0) * clock1_jerk_sps3 * align_t3
        + (1.0 / 24.0) * clock1_snap_sps4 * align_t4;
    let clock2_align_s = clock2_delay_s
        + clock2_rate_sps * read_align_ref_local_s
        + 0.5 * clock2_accel_sps2 * align_t2
        + (1.0 / 6.0) * clock2_jerk_sps3 * align_t3
        + (1.0 / 24.0) * clock2_snap_sps4 * align_t4;
    let read_align_delay_s = geom_delay_at_align_s
        + coarse_delay_s
        + delay_user_samples / fs
        + extra_delay_rate_sps_for_align * read_align_ref_local_s
        + 0.5 * accel_user_sps2 * read_align_ref_local_s * read_align_ref_local_s
        + (clock2_align_s - clock1_align_s);
    let read_align_delay_samples = read_align_delay_s * fs;
    // Directed read-align branch selection.
    //
    // The read-align integer is also the one-sample branch used to decode the
    // packed raw stream. Nearest rounding can select the adjacent low-SNR branch
    // even when the later XCF phase correction is continuous. Select the integer
    // sample on the delayed side of the continuous model delay: ceil for positive
    // relative delay, floor for negative relative delay.
    //
    // The remaining residual is corrected by the XCF frequency-domain phase
    // slope and is not split into per-frame integer shifts in the sector
    // read-align path.
    let init_delay_samples = if read_align_delay_samples >= 0.0 {
        read_align_delay_samples.ceil() as i64
    } else {
        read_align_delay_samples.floor() as i64
    };
    let read_align_residual_samples = read_align_delay_samples - init_delay_samples as f64;
    println!(
        "[info] Read-align reference: {} t={:.6}s, delay={:+.6} sample, fixed integer={} sample, residual={:+.6} sample (directed branch)",
        read_align_ref_label,
        read_align_ref_local_s,
        read_align_delay_samples,
        init_delay_samples,
        read_align_residual_samples
    );
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
    let total_f = complete_fft_frame_count(total_sec, fs, fft_len);

    let (weight1_raw, _, _, _) = resolve_weight(
        tsys1,
        gain1,
        if sefd1 > 0.0 { Some(sefd1) } else { None },
        if dia1 > 0.0 { Some(dia1) } else { None },
        eta1,
        "A1",
    )?;
    let (weight2_raw, _, _, _) = resolve_weight(
        tsys2,
        gain2,
        if sefd2 > 0.0 { Some(sefd2) } else { None },
        if dia2 > 0.0 { Some(dia2) } else { None },
        eta2,
        "A2",
    )?;
    // Absolute SEFD weights can be very small. Normalize their common scale so
    // independent receiver noise occupies the output quantizer power while
    // preserving the scientifically meaningful weight ratio.
    let (w1, w2, beam_voltage_scale) = normalized_beam_weights(
        weight1_raw,
        weight2_raw,
        level_power1,
        level_power2,
        level_power_out,
    )?;
    let a1_name = if_d
        .as_ref()
        .and_then(|d| d.ant1_station_name.as_deref())
        .unwrap_or("YAMAGU32");
    let a2_name = if_d
        .as_ref()
        .and_then(|d| d.ant2_station_name.as_deref())
        .unwrap_or("YAMAGU34");
    let output_format_name = match out_grid {
        OutputGrid::Ant1 => a1_name,
        OutputGrid::Ant2 => a2_name,
    };
    let phased_name = sanitize_file_token(&args.phased_name);
    let a1_key_opt = if_d.as_ref().and_then(|d| d.ant1_station_key.as_deref());
    let a2_key_opt = if_d.as_ref().and_then(|d| d.ant2_station_key.as_deref());
    let a1_key = a1_key_opt.unwrap_or("-");
    let a2_key = a2_key_opt.unwrap_or("-");
    let a1_code = a1_key_opt
        .and_then(|k| k.as_bytes().first().copied())
        .or_else(|| a1_name.as_bytes().first().copied())
        .unwrap_or(b'1');
    let a2_code = a2_key_opt
        .and_then(|k| k.as_bytes().first().copied())
        .or_else(|| a2_name.as_bytes().first().copied())
        .unwrap_or(b'2');
    let schedule_mode = if_d.is_some();
    let correction_sign = -1.0f64; // Correct ant1 toward ant2 while geometric model is (ant2 - ant1).
    let geometric_rate_hz = correction_sign * geom_rate_at_start_sps * obs_mhz * 1e6;
    let clock_rate_hz = clock_rate_sps * obs_mhz * 1e6;
    let df_hz = fs / fft_len as f64;
    let rotation_delta_hz = rot1 - rot2; // ant1 - ant2
    let rotation_shift_hz = raw_ba.shift_bins as f64 * df_hz;
    let rotation_residual_hz = rotation_delta_hz - rotation_shift_hz;
    let rotation_fringe_hz = -rotation_residual_hz; // correction applied on ant1 toward ant2
    let total_rate_hz = geometric_rate_hz + clock_rate_hz + rotation_fringe_hz + rate_user_hz;
    let overlap_bins = ba.a1e - ba.a1s;
    let bin_mhz = fs / fft_len as f64 / 1e6;
    let overlap_mhz = overlap_bins as f64 * bin_mhz;
    let rotation_mode_desc = if rotation_hermitian_mode {
        "rotation2 Hermitian-completion mode"
    } else {
        "rotation real-overlap mode"
    };
    let band_align_desc = format!(
        "{}: frequency-shift after FFT: {} raw shift {} bins -> XML grid; {} raw shift 0 bins -> XML grid; integrate XML-grid[{}..{})",
        rotation_mode_desc, a2_name, raw_ba.shift_bins, a1_name, ba.a1s, ba.a1e
    );
    let rotation_shift_desc = format!(
        "FFT frequency-shift {} by {} bins onto XML-frequency grid",
        a2_name, raw_ba.shift_bins
    );
    let delay_model = "per-frame midpoint delay + integer/fractional correction";

    if schedule_mode {
        if let (Some(sch), Some(meta)) = (&args.schedule, &if_d) {
            let gv = ScheduleGlobalView {
                source_frame: "J2000 catalog; delay model precessed to date",
                process_skip_sec: Some(xml_skip_sec),
                process_length_sec: xml_length_sec,
                geom_delay_s: geom_delay_at_start_s,
                geom_rate_sps: geom_rate_at_start_sps,
                geom_accel_sps2: geom_accel_at_start_sps2,
                geom_eval_offset_s: total_skip_sec,
                read_align_ref_offset_s: read_align_ref_local_s,
                geom_rate_hz_at_obs: geom_rate_at_start_sps * obs_mhz * 1e6,
                rel_clock_delay_s: clock_delay_s,
                rel_clock_rate_sps: clock_rate_sps,
                clock_epoch_mode,
                coarse_delay_s,
                coarse_delay_samples: coarse_delay_s * fs,
                res_delay_input_samples: delay_user_samples,
                read_align_delay_samples: read_align_delay_s * fs,
                read_align_delay_s,
                delay_model,
                model_time_offset_s: model_time_offset_s,
                dut1_s,
                tt_utc_s: tt_utc_s,
                xp_arcsec,
                yp_arcsec,
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
                output_grid: match out_grid {
                    OutputGrid::Ant1 => a1_name,
                    OutputGrid::Ant2 => a2_name,
                },
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
                clock_epoch: clock1_epoch.as_deref(),
                clock_epoch_dt_s: clock1_epoch_offset_s,
                clock_delay_raw_s: clock1_delay_base_s,
                clock_rate_raw_sps: clock1_rate_raw_sps,
                clock_accel_raw_sps2: clock1_accel_raw_sps2,
                clock_jerk_raw_sps3: clock1_jerk_raw_sps3,
                clock_snap_raw_sps4: clock1_snap_raw_sps4,
                clock_delay_s: clock1_delay_s,
                clock_rate_sps: clock1_rate_sps,
                clock_accel_sps2: clock1_accel_sps2,
                clock_jerk_sps3: clock1_jerk_sps3,
                clock_snap_sps4: clock1_snap_sps4,
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
                clock_epoch: clock2_epoch.as_deref(),
                clock_epoch_dt_s: clock2_epoch_offset_s,
                clock_delay_raw_s: clock2_delay_base_s,
                clock_rate_raw_sps: clock2_rate_raw_sps,
                clock_accel_raw_sps2: clock2_accel_raw_sps2,
                clock_jerk_raw_sps3: clock2_jerk_raw_sps3,
                clock_snap_raw_sps4: clock2_snap_raw_sps4,
                clock_delay_s: clock2_delay_s,
                clock_rate_sps: clock2_rate_sps,
                clock_accel_sps2: clock2_accel_sps2,
                clock_jerk_sps3: clock2_jerk_sps3,
                clock_snap_sps4: clock2_snap_sps4,
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
    if !args.compact_logs {
        println!("  mode:       {}", run_mode.label());
    }
    if !schedule_mode {
        println!(
            "  {}:  {} (Size: {} bytes, Estimated Obs Time: {:.2}s)",
            a1_name,
            a1p.display(),
            f1_m.len(),
            (f1_m.len() * 8) as f64 / bit1 as f64 / fs
        );
        println!(
            "  {}:  {} (Size: {} bytes, Estimated Obs Time: {:.2}s)",
            a2_name,
            a2p.display(),
            f2_m.len(),
            (f2_m.len() * 8) as f64 / bit2 as f64 / fs
        );
    }
    if !schedule_mode {
        if let Some(info) = &gdi {
            println!(
                "  ra/dec:     {:.6} / {:.6}",
                info.ra.to_degrees(),
                info.dec.to_degrees()
            );
            println!("  source-frame: J2000 catalog; delay model precessed to date");
            println!("  epoch:      {}", ep_i);
            println!(
                "  geom-eval:  start t={:+.3}s, read-align t={:+.3}s",
                total_skip_sec, read_align_ref_local_s
            );
            println!(
                "  geom-delay: {:.6e} s ({} - {})",
                geom_delay_at_start_s, a2_name, a1_name
            );
            println!(
                "  geom-rate:  {:.6e} s/s ({} - {}) => {:.6} Hz @ obsfreq",
                geom_rate_at_start_sps,
                a2_name,
                a1_name,
                geom_rate_at_start_sps * obs_mhz * 1e6
            );
            println!(
                "  geom-accel: {:.6e} s/s^2 ({} - {})",
                geom_accel_at_start_sps2, a2_name, a1_name
            );
        }
        println!(
            "  coarse-delay fixed: {:.6e} s (relative pre-align) => applied {:.3} samples",
            coarse_delay_s,
            coarse_delay_s * fs
        );
        println!(
            "  clock-delay {}: {:.6e} s => applied {:.3} samples",
            a1_name,
            clock1_delay_s,
            clock1_delay_s * fs
        );
        println!(
            "  clock-delay {}: {:.6e} s => applied {:.3} samples",
            a2_name,
            clock2_delay_s,
            clock2_delay_s * fs
        );
        println!("  clock-rate  {}: {:.6e} s/s", a1_name, clock1_rate_sps);
        println!("  clock-rate  {}: {:.6e} s/s", a2_name, clock2_rate_sps);
        if let Some(ep) = clock1_epoch.as_deref() {
            println!(
                "  clock-epoch {}: {} (to process epoch: {:+.3} s)",
                a1_name,
                ep,
                clock1_epoch_offset_s.unwrap_or(0.0)
            );
        }
        if let Some(ep) = clock2_epoch.as_deref() {
            println!(
                "  clock-epoch {}: {} (to process epoch: {:+.3} s)",
                a2_name,
                ep,
                clock2_epoch_offset_s.unwrap_or(0.0)
            );
        }
    }
    if !schedule_mode {
        println!(
            "  res-delay input: {} samples (relative pre-align)",
            delay_user_samples
        );
        println!(
            "  read-align delay: {:.3} samples ({:.3e} s)",
            read_align_delay_s * fs,
            read_align_delay_s
        );
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
            "  phase-freq: {:.3} MHz ({}) / {:.3} MHz ({}) mode={} (pre-FFT raw delay/rate carriers)",
            lo1_mhz,
            a1_name,
            lo2_mhz,
            a2_name,
            fringe_stop_freq_mode.label()
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
        println!("  sampling:   {:.0} Hz ({:.3} MHz)", fs, fs / 1e6);
    }
    if !args.compact_logs {
        println!("  samples/s:  {:.0}", fs);
        println!("  frames/s:   {:.6}", fs / fft_len as f64);
        println!("  fft:        {}", fft_len);
        println!("  debug:      {}", args.debug);
    }
    if !schedule_mode {
        println!(
            "  bit:        {}={} {}={} -> out={}",
            a1_name, bit1, a2_name, bit2, bit_out
        );
        println!(
            "  bit-code:   {}=({}) {}=({})",
            a1_name,
            format_bit_codes(bit1),
            a2_name,
            format_bit_codes(bit2)
        );
        println!("  level:      {}={:?}", a1_name, levels1);
        println!("  level:      {}={:?}", a2_name, levels2);
        println!(
            "  level-map:  {}={}",
            a1_name,
            format_level_map(bit1, &levels1)
        );
        println!(
            "  level-map:  {}={}",
            a2_name,
            format_level_map(bit2, &levels2)
        );
        println!(
            "  cor-normalization: inv = 1 / (0.5 * P * nf * fft^2), P11={:.6}, P22={:.6}, P12={:.6}",
            level_power1,
            level_power2,
            (level_power1 * level_power2).sqrt()
        );
        println!("  shuffle-in: {}={:?}", a1_name, sh1_ext);
        println!("  shuffle-in: {}={:?}", a2_name, sh2_ext);
        println!("  sideband:   {}={} {}={}", a1_name, sb1_s, a2_name, sb2_s);
        println!(
            "  sideband-normalize: {}={} {}={} (internal USB domain)",
            a1_name,
            if lsb1 { "LSB->USB" } else { "USB" },
            a2_name,
            if lsb2 { "LSB->USB" } else { "USB" }
        );
        println!(
            "  sideband-output: {} ({} raw; inverse normalization={})",
            output_sideband,
            phased_name,
            if output_lsb_restore {
                "odd-sample sign restore"
            } else {
                "none"
            }
        );
        println!(
            "  phased-output-format: {} bit, {}-compatible packed order, level={:?}",
            bit_out, output_format_name, levels_out
        );
        println!(
            "  shuffle-out: {}={:?} (from {})",
            phased_name, output_shuffle_ext, output_format_name
        );
        println!(
            "  phased-diagnostics: {}",
            if args.phased_diagnostics {
                "enabled"
            } else {
                "disabled (raw only)"
            }
        );
        println!("  obs-band:   {:.3} .. {:.3} MHz", obs_mhz, obs_mhz + bw);
    }
    if !schedule_mode {
        println!("  band-align: {band_align_desc}");
        println!(
            "  band-overlap: {} bins ({:.3} MHz, {:.3} MHz/bin)",
            overlap_bins, overlap_mhz, bin_mhz
        );
        println!("  rotation-shift: {rotation_shift_desc}");
        println!(
            "  rotation-fringestop: delta_hz={:.6} shift_hz={:.6} residual_hz={:.6}",
            rotation_delta_hz, rotation_shift_hz, rotation_residual_hz
        );
        println!(
            "  output-grid: {}",
            match out_grid {
                OutputGrid::Ant1 => a1_name,
                OutputGrid::Ant2 => a2_name,
            }
        );
        println!(
            "  ant-fixed:  {}={:?}, {}={:?}",
            a1_name, ant1_ecef, a2_name, ant2_ecef
        );
    }
    if !args.compact_logs {
        println!(
            "  cpu:        {} (compute threads: {})",
            cpu_threads,
            rayon::current_num_threads()
        );
        let build_cpu = build_logical_cpus()
            .map(|v| v.to_string())
            .unwrap_or_else(|| "n/a".to_string());
        let build_l3 = build_l3_cache_bytes()
            .map(|b| format!("{:.1} MiB", b as f64 / (1024.0 * 1024.0)))
            .unwrap_or_else(|| "n/a".to_string());
        println!(
            "  build-host: logical-cpu={} l3-cache={}",
            build_cpu, build_l3
        );
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
        println!(
            "  length:     {:.6}s requested, {:.6}s processed",
            requested_sec, processed_sec
        );
        if is_phased_mode {
            println!(
                "  tsys:       {} ({}), {} ({})",
                tsys1, a1_name, tsys2, a2_name
            );
            println!(
                "  eta:        {} ({}), {} ({})",
                eta1, a1_name, eta2, a2_name
            );
            println!(
                "  gain:       {} ({}), {} ({})",
                gain1, a1_name, gain2, a2_name
            );
            println!(
                "  weight-raw: {:.9e} ({}), {:.9e} ({})",
                weight1_raw, a1_name, weight2_raw, a2_name
            );
            println!(
                "  weight-effective: {:.9e} ({}), {:.9e} ({}) voltage-scale={:.9e}",
                w1, a1_name, w2, a2_name, beam_voltage_scale
            );
            println!("  results:    {}", o_dir.display());
            println!("  output:     {}", o_path.display());
        } else {
            println!("  cor-dir:    {}", o_dir.display());
        }
        println!("--------------------------------------------------");
    } else if !is_phased_mode {
        println!("[info] cor-dir: {}", o_dir.display());
    }

    let fx_read_offset_samples = fx_read_offset_samples();
    let start_s1_samples = start_s1;
    let start_s2_nominal_samples = start_s2;
    let start_s2_with_offset = start_s2_nominal_samples as i128 + fx_read_offset_samples as i128;
    if start_s2_with_offset < 0 {
        return Err(format!(
            "YI_FX_READ_OFFSET_SAMPLE={} makes ant2 input start negative (nominal={})",
            fx_read_offset_samples, start_s2_nominal_samples
        )
        .into());
    }
    let start_s2_samples = start_s2_with_offset as u64;
    if fx_read_offset_samples != 0 {
        println!(
            "[info] FX user read offset: ant2 raw input start {:+} sample(s): nominal {} -> {}",
            fx_read_offset_samples, start_s2_nominal_samples, start_s2_samples
        );
        println!(
            "[info] FX user read offset: phase-delay correction is NOT changed by this offset"
        );
    }
    let samples_per_word1 = 32_u64
        .checked_div(bit1 as u64)
        .ok_or("invalid bit1 for word alignment")?;
    let samples_per_word2 = 32_u64
        .checked_div(bit2 as u64)
        .ok_or("invalid bit2 for word alignment")?;
    if samples_per_word1 == 0 || samples_per_word2 == 0 {
        return Err("invalid samples-per-word while decoding input".into());
    }
    // VSREC shuffle is defined on the original physical 32-bit word boundary.
    // Never bit-shift/repack the physical stream before applying that shuffle:
    // doing so redefines the word boundary and creates four spectral images.
    // Read from the preceding native word and leave the remaining sample
    // displacement to the normal integer/fractional delay correction.
    let start_s1_actual_samples = start_s1_samples / samples_per_word1 * samples_per_word1;
    let start_s2_actual_samples = start_s2_samples / samples_per_word2 * samples_per_word2;
    let start_s1_residual_samples = start_s1_samples as i128 - start_s1_actual_samples as i128;
    let start_s2_residual_samples = start_s2_samples as i128 - start_s2_actual_samples as i128;
    let start_s1_bits = sample_to_total_bits(start_s1_actual_samples, bit1)?;
    let start_s2_bits = sample_to_total_bits(start_s2_actual_samples, bit2)?;
    let start_s1_byte = start_s1_bits / 8;
    let start_s2_byte = start_s2_bits / 8;
    let start_s1_bit = (start_s1_bits % 8) as u8;
    let start_s2_bit = (start_s2_bits % 8) as u8;
    println!(
        "[info] Input start ant1: desired sample {} -> native word byte {} + bit {} (actual sample {}, delay residual {} sample)",
        start_s1_samples, start_s1_byte, start_s1_bit, start_s1_actual_samples, start_s1_residual_samples
    );
    println!(
        "[info] Input start ant2: desired sample {} -> native word byte {} + bit {} (actual sample {}, delay residual {} sample)",
        start_s2_samples, start_s2_byte, start_s2_bit, start_s2_actual_samples, start_s2_residual_samples
    );
    if start_s1_residual_samples != 0 || start_s2_residual_samples != 0 {
        println!(
            "[info] Bit-exact seek residual: ant1={} sample(s), ant2={} sample(s)",
            start_s1_residual_samples, start_s2_residual_samples
        );
    }
    let d_seek_samples = start_s2_actual_samples as i128 - start_s1_actual_samples as i128;
    let d_seek = d_seek_samples as f64 / fs;

    let helper = Arc::new(FftHelper::new(fft_len));
    let frame_dt = fft_len as f64 / fs;
    // Keep integer/fractional shifts on baseline-relative delay only.
    let residual_on_ant2 = init_delay_samples >= 0;
    println!(
        "[info] Residual delay correction target: {}",
        if residual_on_ant2 { "ant2" } else { "ant1" }
    );
    println!(
        "[info] Fringe-stop carrier frequency: mode={} {}={:.6} MHz {}={:.6} MHz",
        fringe_stop_freq_mode.label(),
        a1_name,
        lo1_mhz,
        a2_name,
        lo2_mhz
    );
    let residual_base = net_d_rel_no_clock0 + clock_delay_s - d_seek;
    println!(
        "[info] Delay residual after read-align: {:+.6} samples ({:+.9e} s)",
        residual_base * fs,
        residual_base
    );
    let net_d1_base = if residual_on_ant2 { 0.0 } else { residual_base };
    let net_d2_base = if residual_on_ant2 {
        -residual_base
    } else {
        0.0
    };
    let rate_rel_no_clock_base =
        (geometric_rate_hz + rotation_fringe_hz + rate_user_hz) / (obs_mhz * 1e6);
    let residual_rate_base = rate_rel_no_clock_base + clock_rate_sps;
    let total_rate1_base = if residual_on_ant2 {
        0.0
    } else {
        residual_rate_base
    };
    let total_rate2_base = if residual_on_ant2 {
        -residual_rate_base
    } else {
        0.0
    };
    let total_accel_base =
        correction_sign * geom_accel_at_start_sps2 + clock_accel_sps2 + accel_user_sps2;
    let total_jerk_base = clock_jerk_sps3;
    let total_snap_base = clock_snap_sps4;
    let total_accel1_base = if residual_on_ant2 {
        0.0
    } else {
        total_accel_base
    };
    let total_accel2_base = if residual_on_ant2 {
        -total_accel_base
    } else {
        0.0
    };
    let total_jerk1_base = if residual_on_ant2 {
        0.0
    } else {
        total_jerk_base
    };
    let total_jerk2_base = if residual_on_ant2 {
        -total_jerk_base
    } else {
        0.0
    };
    let total_snap1_base = if residual_on_ant2 {
        0.0
    } else {
        total_snap_base
    };
    let total_snap2_base = if residual_on_ant2 {
        -total_snap_base
    } else {
        0.0
    };
    let geom_delay_table_1s: Option<Arc<[GeomDelaySample]>> = gdi.as_ref().map(|v| {
        let total_duration_sec = total_f as f64 * frame_dt + model_time_offset_s.abs();

        // Need extra guard points because rate/accel/jerk/snap are obtained
        // by finite differences from the 1-s geometric-delay grid.
        let n_update_secs = total_duration_sec.ceil().max(1.0) as usize;
        let n_points = n_update_secs + 5;

        // delay_grid[i] is evaluated at:
        //   t_abs = total_skip_sec + i [s]
        //
        // The grid spacing is exactly 1 s, so finite-difference denominators
        // are unity in SI units.
        // Keep guard points on both sides so rate/accel/jerk/snap are
        // evaluated by a local quartic fit even at the process edges.
        let guard = 10isize;
        let mut delay_grid = Vec::with_capacity(n_points + 2 * guard as usize);
        for sec in -guard..(n_points as isize + guard) {
            let t_abs_sec = total_skip_sec + sec as f64;
            let mjd_t = v.mjd + t_abs_sec / 86400.0;
            let mjd_tt_t = mjd_t + tt_utc_s / 86400.0;
            let (ra_t, dec_t) = match source_vector_mode {
                geom::SourceVectorMode::MeanGast => {
                    geom::precess_j2000_to_mean_of_date(v.ra_raw, v.dec_raw, mjd_tt_t)
                }
                geom::SourceVectorMode::PnmGast | geom::SourceVectorMode::PnmEra => {
                    (v.ra_raw, v.dec_raw)
                }
            };
            let (_, _, gd_t, _gr_t, _ga_t) =
                geom::calculate_geometric_delay_and_derivatives_full_with_eop(
                    ant1_ecef,
                    ant2_ecef,
                    ra_t,
                    dec_t,
                    mjd_t,
                    v.mjd,
                    earth_orientation,
                    geom_delay_mode,
                    source_vector_mode,
                );
            delay_grid.push(gd_t);
        }

        let sample = |i: usize, offset: isize| -> f64 {
            let idx = i as isize + guard + offset;
            delay_grid[idx as usize]
        };

        let fit_radius = guard as usize;
        let mut table = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let center_idx = i + guard as usize;
            let (rate_sps, accel_sps2, jerk_sps3, snap_sps4) =
                model_diag::local_quartic_derivatives(&delay_grid, center_idx, fit_radius);
            table.push(GeomDelaySample {
                delay_s: sample(i, 0),
                rate_sps,
                accel_sps2,
                jerk_sps3,
                snap_sps4,
            });
        }

        Arc::from(table.into_boxed_slice())
    });
    let extra_delay_rate_sps = (rotation_fringe_hz + rate_user_hz) / (obs_mhz * 1e6);
    let extra_delay_accel_sps2 = accel_user_sps2;
    let geom_poly_order = match std::env::var("YI_GEOM_POLY_ORDER") {
        Ok(v) => {
            let order = v
                .trim()
                .parse::<u8>()
                .map_err(|e| format!("invalid YI_GEOM_POLY_ORDER='{v}': {e}"))?;
            if !(1..=4).contains(&order) {
                return Err(format!(
                    "invalid YI_GEOM_POLY_ORDER='{v}': expected 1..4 (1=rate, 2=accel, 3=jerk, 4=snap)"
                )
                .into());
            }
            order
        }
        Err(_) => 4,
    };
    if geom_delay_table_1s.is_some() {
        println!(
            "[info] Delay model refinement: 1-second geometric-delay/rate/accel/jerk/snap table enabled (per-frame eval order {}, process scope, {} entries, model-time offset {:+.3} s)",
            geom_poly_order,
            geom_delay_table_1s.as_ref().map(|t| t.len()).unwrap_or(0),
            model_time_offset_s
        );
        println!(
            "[info] Earth orientation model: DUT1={:+.6}s TT-UTC={:+.6}s xp={:+.6}arcsec yp={:+.6}arcsec source-vector={} delay-mode={}",
            dut1_s,
            tt_utc_s,
            xp_arcsec,
            yp_arcsec,
            source_vector_mode.label(),
            geom_delay_mode.label()
        );
        println!(
            "[info] Delay model time-offset equivalent: {:+.6} s -> geom-rate {:+.6} sample",
            model_time_offset_s,
            model_time_offset_s * geom_rate_at_start_sps * fs
        );
    }
    let delay_cfg = DelayEvalConfig {
        frame_dt,
        model_time_offset_s: model_time_offset_s,
        fs,
        time_offset_s: total_skip_sec,
        geom_delay_table_1s,
        geom_poly_order,
        coarse_delay_s,
        delay_user_samples,
        extra_delay_rate_sps,
        extra_delay_accel_sps2,
        clock1_delay_s,
        clock1_rate_sps,
        clock1_accel_sps2,
        clock1_jerk_sps3,
        clock1_snap_sps4,
        clock2_delay_s,
        clock2_rate_sps,
        clock2_accel_sps2,
        clock2_jerk_sps3,
        clock2_snap_sps4,
        d_seek,
        residual_on_ant2,
        fx_start_cumulative_seek,
        net_d1_base,
        total_rate1_base,
        total_accel1_base,
        total_jerk1_base,
        total_snap1_base,
        net_d2_base,
        total_rate2_base,
        total_accel2_base,
        total_jerk2_base,
        total_snap2_base,
        lo1_hz,
        lo2_hz,
        fx_integer_delay: fx_integer_delay_enabled(),
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
        if_d.as_ref()
            .and_then(|d| d.stream_label.clone())
            .unwrap_or_else(|| "phasedarray".to_string())
    };
    let cor_label = sanitize_file_token(&cor_label_raw);
    let ant1_name_file = sanitize_file_token(a1_name);
    let ant2_name_file = sanitize_file_token(a2_name);

    if args.model_diagnostics {
        if let Some(source) = gdi.as_ref() {
            let diagnostic_path = o_dir.join(format!(
                "model_diagnostics_{}_{}_{}.tsv",
                ant1_name_file, ant2_name_file, c_tag
            ));
            let carrier_hz = if residual_on_ant2 { lo2_hz } else { lo1_hz };
            let geometry = model_diag::GeometryContext {
                ant1_ecef_m: ant1_ecef,
                ant2_ecef_m: ant2_ecef,
                ra_j2000_rad: source.ra_raw,
                dec_j2000_rad: source.dec_raw,
                reference_mjd_utc: source.mjd,
            };
            let clock_coefficients = model_diag::DelayDerivatives {
                delay_s: clock_delay_s,
                rate_sps: clock_rate_sps,
                accel_sps2: clock_accel_sps2,
                jerk_sps3: clock_jerk_sps3,
                snap_sps4: clock_snap_sps4,
            };
            let other_coefficients = model_diag::DelayDerivatives {
                delay_s: coarse_delay_s + delay_user_samples / fs,
                rate_sps: extra_delay_rate_sps,
                accel_sps2: extra_delay_accel_sps2,
                jerk_sps3: 0.0,
                snap_sps4: 0.0,
            };
            let residual_target = if residual_on_ant2 { a2_name } else { a1_name };
            let diagnostic_input = model_diag_output::DiagnosticInput {
                ant1_name: a1_name.to_string(),
                ant2_name: a2_name.to_string(),
                epoch: ep_i.clone(),
                geometry,
                active_source_mode: source_vector_mode,
                active_delay_mode: geom_delay_mode,
                active_eop: earth_orientation,
                start_elapsed_process_s: total_skip_sec + model_time_offset_s,
                total_skip_s: total_skip_sec,
                model_time_offset_s,
                duration_s: processed_sec,
                frame_dt_s: frame_dt,
                total_frames: total_f,
                fs_hz: fs,
                obs_hz: obs_mhz * 1.0e6,
                carrier_hz,
                clock_epoch_coefficients: clock_coefficients,
                other_epoch_coefficients: other_coefficients,
                read_align_reference_s: read_align_ref_local_s,
                d_seek_s: d_seek,
                residual_target: residual_target.to_string(),
            };
            let baseline = [
                ant2_ecef[0] - ant1_ecef[0],
                ant2_ecef[1] - ant1_ecef[1],
                ant2_ecef[2] - ant1_ecef[2],
            ];
            let baseline_length_m =
                (baseline[0] * baseline[0] + baseline[1] * baseline[1] + baseline[2] * baseline[2])
                    .sqrt();
            println!(
                "[geom-diag:begin] one_run=true raw_fft_runs=1 full_model_variants=16 per_second_tsv={} duration_s={:.9}",
                diagnostic_path.display(),
                processed_sec
            );
            println!(
                "[geom-diag:context] epoch={} reference_mjd_utc={:.12} total_skip_s={:+.17e} model_time_offset_s={:+.17e} frame_dt_s={:.17e} frames={} fs_hz={:.9} obs_hz={:.9} carrier_hz={:.9}",
                ep_i,
                source.mjd,
                total_skip_sec,
                model_time_offset_s,
                frame_dt,
                total_f,
                fs,
                obs_mhz * 1.0e6,
                carrier_hz
            );
            println!(
                "[geom-diag:source] ra_j2000_rad={:+.17e} dec_j2000_rad={:+.17e} ra_j2000_deg={:+.12} dec_j2000_deg={:+.12} source_mode={} delay_mode={} poly_order={}",
                source.ra_raw,
                source.dec_raw,
                source.ra_raw.to_degrees(),
                source.dec_raw.to_degrees(),
                source_vector_mode.label(),
                geom_delay_mode.label(),
                geom_poly_order
            );
            println!(
                "[geom-diag:baseline] ant1={} xyz_m={:+.9},{:+.9},{:+.9} ant2={} xyz_m={:+.9},{:+.9},{:+.9} baseline_2minus1_m={:+.9},{:+.9},{:+.9} length_m={:.9}",
                a1_name,
                ant1_ecef[0], ant1_ecef[1], ant1_ecef[2],
                a2_name,
                ant2_ecef[0], ant2_ecef[1], ant2_ecef[2],
                baseline[0], baseline[1], baseline[2], baseline_length_m
            );
            println!(
                "[geom-diag:eop] dut1_s={:+.17e} tt_minus_utc_s={:+.17e} xp_arcsec={:+.17e} yp_arcsec={:+.17e}",
                dut1_s, tt_utc_s, xp_arcsec, yp_arcsec
            );
            println!(
                "[geom-diag:clock-raw] ant1_epoch={} delay_s={:+.17e} rate_sps={:+.17e} accel_sps2={:+.17e} jerk_sps3={:+.17e} snap_sps4={:+.17e} ant2_epoch={} delay_s={:+.17e} rate_sps={:+.17e} accel_sps2={:+.17e} jerk_sps3={:+.17e} snap_sps4={:+.17e}",
                clock1_epoch.as_deref().unwrap_or("process"),
                clock1_delay_base_s,
                clock1_rate_raw_sps,
                clock1_accel_raw_sps2,
                clock1_jerk_raw_sps3,
                clock1_snap_raw_sps4,
                clock2_epoch.as_deref().unwrap_or("process"),
                clock2_delay_base_s,
                clock2_rate_raw_sps,
                clock2_accel_raw_sps2,
                clock2_jerk_raw_sps3,
                clock2_snap_raw_sps4
            );
            println!(
                "[geom-diag:clock-effective] epoch_mode={} ant1_delay_s={:+.17e} ant1_rate_sps={:+.17e} ant1_accel_sps2={:+.17e} ant1_jerk_sps3={:+.17e} ant1_snap_sps4={:+.17e} ant2_delay_s={:+.17e} ant2_rate_sps={:+.17e} ant2_accel_sps2={:+.17e} ant2_jerk_sps3={:+.17e} ant2_snap_sps4={:+.17e}",
                clock_epoch_mode,
                clock1_delay_s,
                clock1_rate_sps,
                clock1_accel_sps2,
                clock1_jerk_sps3,
                clock1_snap_sps4,
                clock2_delay_s,
                clock2_rate_sps,
                clock2_accel_sps2,
                clock2_jerk_sps3,
                clock2_snap_sps4
            );
            println!(
                "[geom-diag:nongeom] coarse_delay_s={:+.17e} user_delay_sample={:+.17e} rotation_residual_hz={:+.17e} user_rate_hz={:+.17e} extra_delay_rate_sps={:+.17e} user_accel_hz_s={:+.17e} extra_delay_accel_sps2={:+.17e} read_align_ref_s={:+.17e} d_seek_s={:+.17e} residual_target={}",
                coarse_delay_s,
                delay_user_samples,
                rotation_fringe_hz,
                rate_user_hz,
                extra_delay_rate_sps,
                accel_user_hzps,
                extra_delay_accel_sps2,
                read_align_ref_local_s,
                d_seek,
                residual_target
            );
            println!(
                "[geom-diag:phase-convention] correction_cycles=-frequency_hz*delay_s rate_hz=-f*tau1 accel_hz_s=-f*tau2 jerk_hz_s2=-f*tau3 polynomial_deg=c0+c1*dt+c2*dt2+c3*dt3+c4*dt4"
            );
            println!(
                "[geom-diag:not-modeled] troposphere ionosphere solid-earth-tide ocean-loading pole-tide station-velocity axis-offset gravitational-delay"
            );
            println!(
                "[geom-diag:warning] pnm-era-engineering is a sign-control path, not standard CIO-plus-ERA"
            );
            let report =
                model_diag_output::generate_model_diagnostics(&diagnostic_input, &diagnostic_path)?;
            println!(
                "[geom-diag:time-origin] scan_elapsed_s=0_is_first_data_timestamp frame_midpoint_is_model_evaluation_only data_and_model_times_are_reported_separately"
            );
            println!(
                "[geom-diag:derivative-stencil] {}",
                model_diag::diagnostic_derivative_stencil_label()
            );
            println!(
                "[geom-diag:file] path={} rows={} variants={} correlation_numerics_changed=false",
                report.path.display(),
                report.rows,
                report.variants
            );
            for checkpoint in &report.checkpoints {
                let total_phase_obs = checkpoint.total.phase_derivatives(obs_mhz * 1.0e6);
                let total_phase_carrier = checkpoint.total.phase_derivatives(carrier_hz);
                let phase_coeff = checkpoint.total.phase_polynomial_deg(obs_mhz * 1.0e6);
                let applied =
                    compute_frame_delay_entry(checkpoint.frame_idx, &delay_cfg, delay_cfg.d_seek);
                let applied_minus_exact_s = applied.full_rel_s - checkpoint.total.delay_s;
                println!(
                    "[geom-diag:checkpoint] point={} frame_idx={} scan_elapsed_s={:.9} data_elapsed_process_s={:.9} data_mjd_utc={:.12} model_elapsed_process_s={:.9} model_mjd_utc={:.12} geom_sample={:+.17e} geom_rate_sample_s={:+.17e} geom_accel_sample_s2={:+.17e} geom_jerk_sample_s3={:+.17e} geom_snap_sample_s4={:+.17e}",
                    checkpoint.label,
                    checkpoint.frame_idx,
                    checkpoint.scan_elapsed_s,
                    checkpoint.data_elapsed_process_s,
                    checkpoint.data_mjd_utc,
                    checkpoint.model_elapsed_process_s,
                    checkpoint.model_mjd_utc,
                    checkpoint.geom.delay_s * fs,
                    checkpoint.geom.rate_sps * fs,
                    checkpoint.geom.accel_sps2 * fs,
                    checkpoint.geom.jerk_sps3 * fs,
                    checkpoint.geom.snap_sps4 * fs
                );
                println!(
                    "[geom-diag:components] point={} clock_delay_sample={:+.17e} clock_rate_sample_s={:+.17e} clock_accel_sample_s2={:+.17e} clock_jerk_sample_s3={:+.17e} clock_snap_sample_s4={:+.17e} other_delay_sample={:+.17e} other_rate_sample_s={:+.17e} other_accel_sample_s2={:+.17e}",
                    checkpoint.label,
                    checkpoint.clock.delay_s * fs,
                    checkpoint.clock.rate_sps * fs,
                    checkpoint.clock.accel_sps2 * fs,
                    checkpoint.clock.jerk_sps3 * fs,
                    checkpoint.clock.snap_sps4 * fs,
                    checkpoint.other.delay_s * fs,
                    checkpoint.other.rate_sps * fs,
                    checkpoint.other.accel_sps2 * fs
                );
                println!(
                    "[geom-diag:total] point={} delay_sample={:+.17e} rate_sample_s={:+.17e} accel_sample_s2={:+.17e} jerk_sample_s3={:+.17e} snap_sample_s4={:+.17e} obs_phase_rate_hz={:+.17e} obs_phase_accel_hz_s={:+.17e} obs_phase_jerk_hz_s2={:+.17e} carrier_phase_rate_hz={:+.17e} carrier_phase_accel_hz_s={:+.17e} carrier_phase_jerk_hz_s2={:+.17e} c1_deg_s={:+.17e} c2_deg_s2={:+.17e} c3_deg_s3={:+.17e} c4_deg_s4={:+.17e}",
                    checkpoint.label,
                    checkpoint.total.delay_s * fs,
                    checkpoint.total.rate_sps * fs,
                    checkpoint.total.accel_sps2 * fs,
                    checkpoint.total.jerk_sps3 * fs,
                    checkpoint.total.snap_sps4 * fs,
                    total_phase_obs.rate_sps,
                    total_phase_obs.accel_sps2,
                    total_phase_obs.jerk_sps3,
                    total_phase_carrier.rate_sps,
                    total_phase_carrier.accel_sps2,
                    total_phase_carrier.jerk_sps3,
                    phase_coeff[1], phase_coeff[2], phase_coeff[3], phase_coeff[4]
                );
                println!(
                    "[geom-diag:applied] point={} full_rel_sample={:+.17e} residual_sample={:+.17e} int1={} int2={} frac1_sample={:+.17e} frac2_sample={:+.17e} applied_minus_exact_sample={:+.17e} applied_minus_exact_phase_deg={:+.17e}",
                    checkpoint.label,
                    applied.full_rel_s * fs,
                    applied.residual_s * fs,
                    applied.int1,
                    applied.int2,
                    applied.frac1 * fs,
                    applied.frac2 * fs,
                    applied_minus_exact_s * fs,
                    -360.0 * obs_mhz * 1.0e6 * applied_minus_exact_s
                );
                let csc_delta_valid = checkpoint.ant1_csc_el_valid && checkpoint.ant2_csc_el_valid;
                let csc_delta = if csc_delta_valid {
                    checkpoint.ant2_csc_el - checkpoint.ant1_csc_el
                } else {
                    0.0
                };
                println!(
                    "[geom-diag:azel] point={} ant1={} az_deg={:.9} el_deg={:.9} csc_el={:+.17e} csc_valid={} ant2={} az_deg={:.9} el_deg={:.9} csc_el={:+.17e} csc_valid={} csc_delta_2minus1={:+.17e} csc_delta_valid={} u_lambda={:+.17e} v_lambda={:+.17e} w_lambda={:+.17e} troposphere_mode=not-modeled",
                    checkpoint.label,
                    a1_name,
                    checkpoint.ant1_az_deg,
                    checkpoint.ant1_el_deg,
                    checkpoint.ant1_csc_el,
                    checkpoint.ant1_csc_el_valid,
                    a2_name,
                    checkpoint.ant2_az_deg,
                    checkpoint.ant2_el_deg,
                    checkpoint.ant2_csc_el,
                    checkpoint.ant2_csc_el_valid,
                    csc_delta,
                    csc_delta_valid,
                    checkpoint.u_lambda,
                    checkpoint.v_lambda,
                    checkpoint.w_lambda
                );
                for variant in &checkpoint.variants {
                    let geom_phase = variant.geom.phase_derivatives(obs_mhz * 1.0e6);
                    let delta_phase = variant.delta.phase_derivatives(obs_mhz * 1.0e6);
                    println!(
                        "[geom-diag:variant] point={} model={} geom_delay_sample={:+.17e} geom_rate_hz={:+.17e} geom_accel_hz_s={:+.17e} geom_jerk_hz_s2={:+.17e} delta_delay_sample={:+.17e} delta_rate_hz={:+.17e} delta_accel_hz_s={:+.17e} delta_jerk_hz_s2={:+.17e} delta_phase_from_start_deg={:+.17e} delta_phase_after_linear_deg={:+.17e} delta_phase_after_quadratic_deg={:+.17e}",
                        checkpoint.label,
                        variant.key,
                        variant.geom.delay_s * fs,
                        geom_phase.rate_sps,
                        geom_phase.accel_sps2,
                        geom_phase.jerk_sps3,
                        variant.delta.delay_s * fs,
                        delta_phase.rate_sps,
                        delta_phase.accel_sps2,
                        delta_phase.jerk_sps3,
                        variant.delta_phase_from_start_deg,
                        variant.delta_phase_after_linear_deg,
                        variant.delta_phase_after_quadratic_deg
                    );
                }
            }
            println!(
                "[geom-diag:end] collect={} and stdout log; no additional correlation is required for parameter comparison",
                report.path.display()
            );
        } else {
            return Err(
                "--model-diagnostics requires source coordinates; no diagnostic TSV was written"
                    .into(),
            );
        }
    }

    let dump_delay_model = args.model_diagnostics
        || std::env::var("YI_DUMP_DELAY_MODEL")
            .map(|v| {
                let v = v.trim().to_ascii_lowercase();
                !(v.is_empty() || v == "0" || v == "false" || v == "off" || v == "no")
            })
            .unwrap_or(false);
    if dump_delay_model {
        std::fs::create_dir_all(&o_dir)?;
        let dump_path = o_dir.join(format!(
            "delay_model_{}_{}_{}.tsv",
            ant1_name_file, ant2_name_file, c_tag
        ));
        let mut w = BufWriter::new(File::create(&dump_path)?);
        writeln!(
            w,
            "# t_sec\tfull_rel_sample\tresidual_sample\tmodel_rate_sample_s\tmodel_accel_sample_s2\tmodel_jerk_sample_s3\tmodel_snap_sample_s4\tu_lambda\tv_lambda\tw_lambda\tphase_dra_cycles_per_rad\tphase_ddec_cycles_per_rad\tphase_dbx_cycles_per_m\tphase_dby_cycles_per_m\tphase_dbz_cycles_per_m\ttau1_sample\ttau2_sample\tint1\tint2\tfrac1_sample\tfrac2_sample\tfringe_phase_cycles_xmlfreq\tfringe_phase_deg_wrapped\tcarrier_hz\tcarrier_phase_cycles\tcarrier_phase_deg_wrapped"
        )?;
        let dump_secs = processed_sec.ceil().max(1.0) as usize;
        for sec in 0..=dump_secs {
            let frame_idx = ((sec as f64) / frame_dt).round() as usize;
            let frame_idx = frame_idx.min(total_f.saturating_sub(1));
            let d = compute_frame_delay_entry(frame_idx, &delay_cfg, delay_cfg.d_seek);
            let model_sec = (d.t_mid_s - total_skip_sec).floor().max(0.0) as usize;
            let (model_rate, model_accel, model_jerk, model_snap) = if let Some(table) =
                delay_cfg.geom_delay_table_1s.as_deref()
            {
                let idx = model_sec.min(table.len().saturating_sub(1));
                let g = table[idx];
                (
                    (g.rate_sps
                        + delay_cfg.extra_delay_rate_sps
                        + delay_cfg.extra_delay_accel_sps2 * d.t_mid_s
                        + (delay_cfg.clock2_rate_sps - delay_cfg.clock1_rate_sps)
                        + (delay_cfg.clock2_accel_sps2 - delay_cfg.clock1_accel_sps2) * d.t_mid_s
                        + 0.5
                            * (delay_cfg.clock2_jerk_sps3 - delay_cfg.clock1_jerk_sps3)
                            * d.t_mid_s
                            * d.t_mid_s
                        + (1.0 / 6.0)
                            * (delay_cfg.clock2_snap_sps4 - delay_cfg.clock1_snap_sps4)
                            * d.t_mid_s
                            * d.t_mid_s
                            * d.t_mid_s)
                        * fs,
                    (g.accel_sps2
                        + delay_cfg.extra_delay_accel_sps2
                        + (delay_cfg.clock2_accel_sps2 - delay_cfg.clock1_accel_sps2)
                        + (delay_cfg.clock2_jerk_sps3 - delay_cfg.clock1_jerk_sps3) * d.t_mid_s
                        + 0.5
                            * (delay_cfg.clock2_snap_sps4 - delay_cfg.clock1_snap_sps4)
                            * d.t_mid_s
                            * d.t_mid_s)
                        * fs,
                    (g.jerk_sps3
                        + (delay_cfg.clock2_jerk_sps3 - delay_cfg.clock1_jerk_sps3)
                        + (delay_cfg.clock2_snap_sps4 - delay_cfg.clock1_snap_sps4) * d.t_mid_s)
                        * fs,
                    (g.snap_sps4 + (delay_cfg.clock2_snap_sps4 - delay_cfg.clock1_snap_sps4)) * fs,
                )
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
            let basis = if let Some(v) = gdi.as_ref() {
                let mjd_t = v.mjd + d.t_mid_s / 86400.0;
                let mjd_tt_t = mjd_t + tt_utc_s / 86400.0;
                let (ra_t, dec_t) = match source_vector_mode {
                    geom::SourceVectorMode::MeanGast => {
                        geom::precess_j2000_to_mean_of_date(v.ra_raw, v.dec_raw, mjd_tt_t)
                    }
                    geom::SourceVectorMode::PnmGast | geom::SourceVectorMode::PnmEra => {
                        (v.ra_raw, v.dec_raw)
                    }
                };
                geom::baseline_phase_basis(
                    ant1_ecef,
                    ant2_ecef,
                    ra_t,
                    dec_t,
                    mjd_t,
                    earth_orientation,
                    source_vector_mode,
                    obs_mhz * 1.0e6,
                )
            } else {
                geom::BaselinePhaseBasis {
                    u_lambda: 0.0,
                    v_lambda: 0.0,
                    w_lambda: 0.0,
                    phase_dra_cycles_per_rad: 0.0,
                    phase_ddec_cycles_per_rad: 0.0,
                    phase_dbx_cycles_per_m: 0.0,
                    phase_dby_cycles_per_m: 0.0,
                    phase_dbz_cycles_per_m: 0.0,
                }
            };
            let phase_cycles = -obs_mhz * 1.0e6 * d.full_rel_s;
            let wrapped = ((phase_cycles.rem_euclid(1.0) + 0.5).rem_euclid(1.0) - 0.5) * 360.0;
            let carrier_hz = if residual_on_ant2 { lo2_hz } else { lo1_hz };
            let carrier_phase_cycles = -carrier_hz * d.full_rel_s;
            let carrier_wrapped =
                ((carrier_phase_cycles.rem_euclid(1.0) + 0.5).rem_euclid(1.0) - 0.5) * 360.0;
            writeln!(
                w,
                "{:.6}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.6}\t{:.6}\t{:.6}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{}\t{}\t{:.12}\t{:.12}\t{:.12}\t{:.6}\t{:.6}\t{:.12}\t{:.6}",
                sec as f64,
                d.full_rel_s * fs,
                d.residual_s * fs,
                model_rate,
                model_accel,
                model_jerk,
                model_snap,
                basis.u_lambda,
                basis.v_lambda,
                basis.w_lambda,
                basis.phase_dra_cycles_per_rad,
                basis.phase_ddec_cycles_per_rad,
                basis.phase_dbx_cycles_per_m,
                basis.phase_dby_cycles_per_m,
                basis.phase_dbz_cycles_per_m,
                d.tau1_s * fs,
                d.tau2_s * fs,
                d.int1,
                d.int2,
                d.frac1 * fs,
                d.frac2 * fs,
                phase_cycles,
                wrapped,
                carrier_hz,
                carrier_phase_cycles,
                carrier_wrapped,
            )?;
        }
        w.flush()?;
        println!("[info] Delay model dump: {}", dump_path.display());
    }

    if do_xcf {
        let mut ac = vec![Complex::new(0.0_f64, 0.0_f64); fft_len / 2 + 1];
        println!(
            "[info] XCF pipeline: chunk={} frames, depth={} chunks",
            io_chunk_frames, io_pipeline_depth
        );
        let (tx, rx) = mpsc::sync_channel::<Result<Vec<Vec<u8>>, String>>(io_pipeline_depth);
        let xcf_produced_chunks = Arc::new(AtomicUsize::new(0));
        let xcf_produced_bytes = Arc::new(AtomicU64::new(0));
        let xcf_consumed_chunks = Arc::new(AtomicUsize::new(0));
        let (r1_p, r2_p, s1_b, s2_b) = (a1p.clone(), a2p.clone(), start_s1_byte, start_s2_byte);
        let (s1_bit, s2_bit) = (start_s1_bit, start_s2_bit);
        let xcf_produced_chunks_rd = Arc::clone(&xcf_produced_chunks);
        let xcf_produced_bytes_rd = Arc::clone(&xcf_produced_bytes);
        thread::spawn(move || {
            if let Some(core) = reader_core {
                let _ = affinity::set_current_thread_core(core);
            }
            let mut rd1 = match PackedSampleReader::open(&r1_p, s1_b, s1_bit) {
                Ok(v) => v,
                Err(e) => {
                    let _ = tx.send(Err(format!("failed to open/seek {}: {e}", r1_p.display())));
                    return;
                }
            };
            let mut rd2 = match PackedSampleReader::open(&r2_p, s2_b, s2_bit) {
                Ok(v) => v,
                Err(e) => {
                    let _ = tx.send(Err(format!("failed to open/seek {}: {e}", r2_p.display())));
                    return;
                }
            };
            let mut nf = 0;
            while nf < total_f {
                let n = (total_f - nf).min(io_chunk_frames);
                let (mut b1, mut b2) = (vec![0u8; n * bpf1], vec![0u8; n * bpf2]);
                if let Err(e) = rd1.read_packed_with_padding(&mut b1) {
                    let _ = tx.send(Err(format!("failed reading ant1 input: {e}")));
                    return;
                }
                if let Err(e) = rd2.read_packed_with_padding(&mut b2) {
                    let _ = tx.send(Err(format!("failed reading ant2 input: {e}")));
                    return;
                }
                let chunk_bytes = (b1.len() + b2.len()) as u64;
                if tx.send(Ok(vec![b1, b2])).is_err() {
                    return;
                }
                xcf_produced_bytes_rd.fetch_add(chunk_bytes, Ordering::Relaxed);
                xcf_produced_chunks_rd.fetch_add(1, Ordering::Relaxed);
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
        let xcf_stats_start = Instant::now();
        let mut xcf_last_report_at = xcf_stats_start;
        let mut xcf_read_bytes_total: u64 = 0;
        let mut xcf_last_report_bytes: u64 = 0;
        let mut xcf_queue_hwm = 0usize;
        for bufs_res in rx {
            let bufs =
                bufs_res.map_err(|e| std::io::Error::other(format!("xcf reader error: {e}")))?;
            xcf_consumed_chunks.fetch_add(1, Ordering::Relaxed);
            let (raw1, raw2) = (&bufs[0], &bufs[1]);
            let nf = raw1.len() / bpf1;
            xcf_read_bytes_total += (raw1.len() + raw2.len()) as u64;
            let produced = xcf_produced_chunks.load(Ordering::Relaxed);
            let consumed = xcf_consumed_chunks.load(Ordering::Relaxed);
            let queue_fill = produced.saturating_sub(consumed);
            if queue_fill > xcf_queue_hwm {
                xcf_queue_hwm = queue_fill;
            }
            let frame_delays: Vec<FrameDelayEntry> = (0..nf)
                .map(|i| compute_frame_delay_entry(processed + i, &delay_cfg, delay_cfg.d_seek))
                .collect();
            if args.debug {
                print_delay_debug_samples("xcf chunk", processed, &frame_delays, fs);
            }
            let frame_errors = AtomicUsize::new(0);
            struct XcfThreadAccum {
                acc: Vec<Complex<f64>>,
                f1: Vec<f32>,
                f2: Vec<f32>,
                s1: Vec<Complex<f32>>,
                s2: Vec<Complex<f32>>,
                fft_scratch: FftScratch,
                dw1: DecodeWindowScratch,
                dw2: DecodeWindowScratch,
            }
            let half = fft_len / 2 + 1;
            let init = || XcfThreadAccum {
                acc: vec![Complex::new(0.0_f64, 0.0_f64); half],
                f1: vec![0.0_f32; fft_len],
                f2: vec![0.0_f32; fft_len],
                s1: vec![Complex::new(0.0_f32, 0.0_f32); half],
                s2: vec![Complex::new(0.0_f32, 0.0_f32); half],
                fft_scratch: helper.make_scratch(),
                dw1: DecodeWindowScratch::new(),
                dw2: DecodeWindowScratch::new(),
            };
            let frames_per_job = (nf / (cpu_threads.saturating_mul(8)).max(1)).clamp(128, 2048);
            let chunk_starts: Vec<usize> = (0..nf).step_by(frames_per_job).collect();
            let chunk_abs_start1 = start_s1_actual_samples + processed as u64 * fft_len as u64;
            let chunk_abs_start2 = start_s2_actual_samples + processed as u64 * fft_len as u64;
            let mut out = chunk_starts
                .into_par_iter()
                .map(|start| {
                    let mut st = init();
                    let end = (start + frames_per_job).min(nf);
                    for i in start..end {
                        let d = frame_delays[i];
                        if decode_shifted_frame_from_chunk(
                            raw1,
                            chunk_abs_start1,
                            i,
                            fft_len,
                            bit1,
                            samples_per_word1,
                            &dp1,
                            lsb1,
                            d.int1,
                            &mut st.f1,
                            &mut st.dw1,
                        )
                        .is_err()
                        {
                            frame_errors.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        if decode_shifted_frame_from_chunk(
                            raw2,
                            chunk_abs_start2,
                            i,
                            fft_len,
                            bit2,
                            samples_per_word2,
                            &dp2,
                            lsb2,
                            d.int2,
                            &mut st.f2,
                            &mut st.dw2,
                        )
                        .is_err()
                        {
                            frame_errors.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }
                        if helper
                            .forward_r2c_process_with_scratch(
                                &mut st.f1,
                                &mut st.s1,
                                &mut st.fft_scratch,
                            )
                            .is_err()
                        {
                            frame_errors.fetch_add(1, Ordering::Relaxed);
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
                            frame_errors.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }

                        // Use observing reference frequency in .cor header (not upper band edge)
                        // so residual-rate definition matches frinZ delay/rate search.
                        let fr_mix = d.fr_lo1 * d.fr_lo2.conj();
                        let mut g1 = vec![Complex::new(0.0_f32, 0.0_f32); fft_len / 2 + 1];
                        let mut g2 = vec![Complex::new(0.0_f32, 0.0_f32); fft_len / 2 + 1];
                        shift_real_fft_to_xml_grid_with_extra_offset(
                            &st.s1,
                            &mut g1,
                            fft_len,
                            rotation_bins1,
                            station_grid_origin_offset_bins(a1_name, fs, fft_len),
                        );
                        shift_real_fft_to_xml_grid_with_extra_offset(
                            &st.s2,
                            &mut g2,
                            fft_len,
                            rotation_bins2,
                            station_grid_origin_offset_bins(a2_name, fs, fft_len)
                                + ant2_grid_extra_offset(),
                        );

                        if fft_peak_dbg_enabled() && i < fft_peak_dbg_max_frames() {
                            print_fft_peak_dbg(
                                "raw_to_grid",
                                i,
                                a1_name,
                                &st.s1,
                                &g1,
                                rotation_bins1,
                                0,
                            );
                            print_fft_peak_dbg(
                                "raw_to_grid",
                                i,
                                a2_name,
                                &st.s2,
                                &g2,
                                rotation_bins2,
                                station_grid_origin_offset_bins(a2_name, fs, fft_len)
                                    + ant2_grid_extra_offset(),
                            );
                        }
                        let phase_bin_start1 = ba.a1s as isize - rotation_bins1;
                        let phase_bin_start2 = ba.a2s as isize - rotation_bins2;
                        let phase_delay1_s = d.frac1;
                        let phase_delay2_s = d.frac2;
                        let (mut phase_corr, phase_step) = xcf_phase_start_and_step(
                            df_hz,
                            phase_delay1_s,
                            phase_delay2_s,
                            phase_bin_start1,
                            phase_bin_start2,
                        );
                        for k in 0..(ba.a1e - ba.a1s) {
                            let i1 = ba.a1s + k;
                            let i2 = ba.a2s + k;
                            let v = (g1[i1] * g2[i2].conj()) * fr_mix * phase_corr;
                            match out_grid {
                                OutputGrid::Ant1 => {
                                    st.acc[i1] += Complex::new(v.re as f64, v.im as f64)
                                }
                                OutputGrid::Ant2 => {
                                    st.acc[i2] += Complex::new(v.re as f64, v.im as f64)
                                }
                            }
                            phase_corr *= phase_step;
                        }
                    }
                    st
                })
                .reduce(init, |mut a, b| {
                    for k in 0..half {
                        a.acc[k] += b.acc[k];
                    }
                    a
                });
            let bc = std::mem::take(&mut out.acc);
            dropped_xcf_frames += frame_errors.load(Ordering::Relaxed);
            for k in 0..ac.len() {
                ac[k] += bc[k] * inv_fft2;
            }
            processed += nf;
            print!("\rCorrelating ({}/{})", processed, total_f);
            let _ = std::io::stdout().flush();
            if !args.compact_logs && xcf_last_report_at.elapsed() >= Duration::from_secs(2) {
                let elapsed = xcf_stats_start.elapsed().as_secs_f64().max(1e-9);
                let dt = xcf_last_report_at.elapsed().as_secs_f64().max(1e-9);
                let avg_mib_s = (xcf_read_bytes_total as f64 / (1024.0 * 1024.0)) / elapsed;
                let inst_mib_s = ((xcf_read_bytes_total - xcf_last_report_bytes) as f64
                    / (1024.0 * 1024.0))
                    / dt;
                println!(
                    "\n[info] XCF I/O reader: avg={:.1} MiB/s inst={:.1} MiB/s queue={}/{} hwm={}/{}",
                    avg_mib_s,
                    inst_mib_s,
                    queue_fill,
                    io_pipeline_depth,
                    xcf_queue_hwm,
                    io_pipeline_depth
                );
                xcf_last_report_at = Instant::now();
                xcf_last_report_bytes = xcf_read_bytes_total;
            }
        }
        println!();
        if !args.compact_logs {
            let elapsed = xcf_stats_start.elapsed().as_secs_f64().max(1e-9);
            let avg_mib_s = (xcf_read_bytes_total as f64 / (1024.0 * 1024.0)) / elapsed;
            let produced = xcf_produced_chunks.load(Ordering::Relaxed);
            let consumed = xcf_consumed_chunks.load(Ordering::Relaxed);
            let queue_fill = produced.saturating_sub(consumed);
            println!(
                "[info] XCF I/O reader summary: avg={:.1} MiB/s total={:.1} MiB queue_end={}/{} hwm={}/{}",
                avg_mib_s,
                xcf_read_bytes_total as f64 / (1024.0 * 1024.0),
                queue_fill,
                io_pipeline_depth,
                xcf_queue_hwm,
                io_pipeline_depth
            );
        }
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
        let raw_partial_path = o_path.with_extension("raw.part");
        let mut wr: Option<BufWriter<File>> = if write_raw {
            Some(BufWriter::new(File::create(&raw_partial_path)?))
        } else {
            None
        };
        let mut dbg_writer: Option<BufWriter<File>> = if args.debug {
            let debug_base_dir = args
                .schedule
                .as_ref()
                .and_then(|p| p.parent().map(PathBuf::from))
                .unwrap_or_else(|| o_dir.clone());
            let debug_dir = debug_base_dir.join("debug_yi-corr");
            std::fs::create_dir_all(&debug_dir)?;
            let dbg_path = debug_dir.join(format!("debug_{}.log", c_tag));
            println!("[info] Debug log: {}", dbg_path.display());
            let mut w = BufWriter::new(File::create(&dbg_path)?);
            let processed_sec = total_f as f64 * fft_len as f64 / fs;
            writeln!(w, "# phased_array debug log")?;
            writeln!(
                w,
                "# schedule_mode={schedule_mode} fft={fft_len} fs_hz={fs:.3}"
            )?;
            writeln!(
                w,
                "# process_window: start={} end_floor={} end_offset_s={:.6} length_req_s={:.6} length_proc_s={:.6}",
                c_tag,
                unix_seconds_to_yyyydddhhmmss(c_unix + display_whole_seconds(processed_sec))?,
                processed_sec,
                requested_sec,
                processed_sec
            )?;
            writeln!(
                w,
                "# geom_time_reference: process_epoch={} t_abs=process_epoch+(xml_skip+cli_skip)+frame_mid",
                ep_i
            )?;
            writeln!(
                w,
                "# phase_convention: x12 = X1 * conj(X2), fringe_mix = fr_lo1 * conj(fr_lo2)"
            )?;
            writeln!(
                w,
                "# delay_convention: model=(ant2-ant1), read-align by baseline-relative delay, residual applied on delayed antenna"
            )?;
            writeln!(
                w,
                "# read_align_components_s: geom={:.9e} coarse={:.9e} user={:.9e} clock={:.9e} total={:.9e}",
                geom_delay_at_start_s,
                coarse_delay_s,
                delay_user_samples / fs,
                clock_delay_s,
                net_d0
            )?;
            writeln!(
                w,
                "# read_align_components_samples: geom={:.6} coarse={:.6} user={:.6} clock={:.6} total={:.6}",
                geom_delay_at_start_s * fs,
                coarse_delay_s * fs,
                delay_user_samples,
                clock_delay_s * fs,
                net_d0 * fs
            )?;
            writeln!(
                w,
                "# delay_rate_terms_hz: total={:.9} geom={:.9} clock={:.9} rot_res={:.9} user={:.9}",
                total_rate_hz,
                geometric_rate_hz,
                clock_rate_hz,
                rotation_fringe_hz,
                rate_user_hz
            )?;
            writeln!(
                w,
                "# delay_accel_terms_hz_s: total={:.9e} geom={:.9e} clock={:.9e} user={:.9e}",
                (geom_accel_at_start_sps2 + clock_accel_sps2 + accel_user_sps2) * obs_mhz * 1e6,
                geom_accel_at_start_sps2 * obs_mhz * 1e6,
                clock_accel_sps2 * obs_mhz * 1e6,
                accel_user_hzps
            )?;
            writeln!(
                w,
                "# xcf_phase_bins: a1=[{}..{}) a2=[{}..{}) rotation_bins=({}, {}) phase_start_raw=({}, {}) phase_start_xml=({}, {}) overlap_bins={}",
                ba.a1s,
                ba.a1e,
                ba.a2s,
                ba.a2s + ba.a1e.saturating_sub(ba.a1s),
                rotation_bins1,
                rotation_bins2,
                ba.a1s as isize - rotation_bins1,
                ba.a2s as isize - rotation_bins2,
                ba.a1s,
                ba.a2s,
                ba.a1e.saturating_sub(ba.a1s)
            )?;
            writeln!(
                w,
                "# clock_epoch: ant1={} dt_s={:+.3} ant2={} dt_s={:+.3}",
                clock1_epoch.as_deref().unwrap_or("-"),
                clock1_epoch_offset_s.unwrap_or(0.0),
                clock2_epoch.as_deref().unwrap_or("-"),
                clock2_epoch_offset_s.unwrap_or(0.0)
            )?;
            writeln!(
                w,
                "# seek_ant1: desired_samples={} actual_samples={} residual_samples={} byte_offset={} bit_offset={} byte_floor_samples={} byte_rounding_loss_samples={}",
                start_s1_samples,
                start_s1_actual_samples,
                start_s1_residual_samples,
                start_s1_byte,
                start_s1_bit,
                (start_s1_byte * 8) / bit1 as u64,
                start_s1_residual_samples.unsigned_abs()
            )?;
            writeln!(
                w,
                "# seek_ant2: desired_samples={} actual_samples={} residual_samples={} byte_offset={} bit_offset={} byte_floor_samples={} byte_rounding_loss_samples={}",
                start_s2_samples,
                start_s2_actual_samples,
                start_s2_residual_samples,
                start_s2_byte,
                start_s2_bit,
                (start_s2_byte * 8) / bit2 as u64,
                start_s2_residual_samples.unsigned_abs()
            )?;
            writeln!(w, "# frame_logging: all frames")?;
            Some(w)
        } else {
            None
        };
        let frame_sec = fft_len as f64 / fs;
        let output_hz = if_d.as_ref().and_then(|d| d.output_sec).unwrap_or(1.0);

        if !output_hz.is_finite() || output_hz <= 0.0 {
            return Err(format!(
                "stream/output must be positive finite rate [Hz], got {}",
                output_hz
            )
            .into());
        }

        let output_sec = 1.0 / output_hz;
        let frames_per_sector = (output_sec / frame_sec).round().max(1.0) as usize;
        let actual_output_sec = frames_per_sector as f64 * frame_sec;
        let mut sec_counts = Vec::new();
        let mut remaining = total_f;
        while remaining > 0 {
            let nf = remaining.min(frames_per_sector);
            sec_counts.push(nf);
            remaining -= nf;
        }
        println!(
            "[info] Output integration: XML stream/output requested {:.9} Hz, actual {:.9} s, frames/sector {}, sectors {}",
            output_sec,
            actual_output_sec,
            frames_per_sector,
            sec_counts.len()
        );
        let dp1 = Arc::new(build_decode_plan(bit1, sh1.as_ref(), levels1.as_ref())?);
        let dp2 = Arc::new(build_decode_plan(bit2, sh2.as_ref(), levels2.as_ref())?);
        let mut emitted = 0;
        let cor_h_freq_hz = obs_mhz * 1e6;
        let inband = if_d.as_ref().and_then(|d| d.inband).unwrap_or(1);
        if inband == 0 || !inband.is_power_of_two() {
            return Err(format!(
                "stream/inband must be a positive power of 2, got {}",
                inband
            )
            .into());
        }
        let cor_bins_total = fft_len / 2;
        if cor_bins_total % inband != 0 {
            return Err(format!(
                "stream/inband={} does not evenly divide {} positive-frequency bins",
                inband, cor_bins_total
            )
            .into());
        }
        let inband_bins = cor_bins_total / inband;
        let inband_fft_len = fft_len / inband;
        if inband_fft_len < 2 || inband_fft_len % 2 != 0 {
            return Err(format!(
                "stream/inband={} leaves invalid .cor FFT length {}",
                inband, inband_fft_len
            )
            .into());
        }
        let inband_fs = fs / inband as f64;
        let inband_bw_mhz = bw / inband as f64;
        let bw_tag = if (bw - bw.round()).abs() < 1.0e-6 {
            format!("{:.0}", bw.round())
        } else {
            format!("{:.3}", bw).replace(".", "p")
        };
        if inband > 1 {
            println!(
                "[info] In-band output: split {} MHz into {} x {:.6} MHz .cor files ({} bins each)",
                bw_tag, inband, inband_bw_mhz, inband_bins
            );
        }
        let pulsar_runtime = if let Some(cfg) = if_d.as_ref().and_then(|d| d.pulsar.clone()) {
            if write_raw {
                return Err("XML <pulsar> folding is only supported for .cor spectral products, not phased raw output".into());
            }
            if fringe_interval_s.is_some() {
                return Err(
                    "--fringe QL output is not yet supported together with XML <pulsar> folding"
                        .into(),
                );
            }
            let rt = PulsarRuntime::new(cfg, &ep_i, obs_mhz, df_hz, cor_bins_total)?;
            println!(
                "[info] Pulsar folding: name={} bins={} period={:.12e}s pdot={:.6e} pddot={:.6e} dm={:.6e} dedisperse={} time-correction={}",
                rt.config().name.as_deref().unwrap_or("-"),
                rt.bins(),
                rt.config().period_s,
                rt.config().pdot,
                rt.config().pddot,
                rt.config().dm,
                rt.config().dedisperse,
                rt.config().time_correction
            );
            Some(rt)
        } else {
            None
        };
        let pulsar_bin_count = pulsar_runtime.as_ref().map(|rt| rt.bins()).unwrap_or(1);
        let pbin_tag = |pbin: usize| {
            format!(
                "pbin{:0width$}",
                pbin,
                width = pulsar_bin_count.saturating_sub(1).to_string().len().max(2)
            )
        };
        let make_cor_path = |left: &str, right: &str, ch: usize, pbin: Option<usize>| {
            let pbin_suffix = pbin
                .map(|b| format!(".{}", pbin_tag(b)))
                .unwrap_or_default();
            if inband == 1 {
                o_dir.join(format!(
                    "{}_{}_{}_{}{}.cor",
                    left, right, c_tag, cor_label, pbin_suffix
                ))
            } else {
                o_dir.join(format!(
                    "{}_{}_{}_{}.ch{}bw{}{}.cor",
                    left,
                    right,
                    c_tag,
                    cor_label,
                    ch + 1,
                    bw_tag,
                    pbin_suffix
                ))
            }
        };
        let make_cor_header = |ch: usize| CorHeaderConfig {
            sampling_speed_hz: if inband == 1 {
                fs.round() as i32
            } else {
                inband_fs.round() as i32
            },
            observing_frequency_hz: if inband == 1 {
                cor_h_freq_hz
            } else {
                (obs_mhz + ch as f64 * inband_bw_mhz) * 1.0e6
            },
            fft_point: if inband == 1 {
                fft_len as i32
            } else {
                inband_fft_len as i32
            },
            number_of_sector_hint: sec_counts.len() as i32,
            clock_reference_unix_sec: c_unix,
            source_name: source_name.clone(),
            source_ra_rad: gdi.as_ref().map(|v| v.ra_header).unwrap_or(0.0),
            source_dec_rad: gdi.as_ref().map(|v| v.dec_header).unwrap_or(0.0),
        };
        let create_cor_writers = |enabled: bool,
                                  left_file: &str,
                                  right_file: &str,
                                  st1: CorStation<'_>,
                                  st2: CorStation<'_>|
         -> Result<Vec<CorWriter>, DynError> {
            if !enabled {
                return Ok(Vec::new());
            }
            let mut writers = Vec::with_capacity(inband * pulsar_bin_count);
            for pbin in 0..pulsar_bin_count {
                for ch in 0..inband {
                    writers.push(CorWriter::create(
                        &make_cor_path(
                            left_file,
                            right_file,
                            ch,
                            if pulsar_bin_count > 1 {
                                Some(pbin)
                            } else {
                                None
                            },
                        ),
                        &make_cor_header(ch),
                        st1,
                        st2,
                    )?);
                }
            }
            Ok(writers)
        };
        let (phased_reference_name, phased_reference_ecef) = if residual_on_ant2 {
            (a1_name, ant1_ecef)
        } else {
            (a2_name, ant2_ecef)
        };
        let phased_code = phased_name.as_bytes().first().copied().unwrap_or(b'M');
        let mut cw_ph = create_cor_writers(
            write_phased_cor,
            &phased_name,
            &phased_name,
            CorStation {
                name: &phased_name,
                code: phased_code,
                ecef_m: phased_reference_ecef,
            },
            CorStation {
                name: &phased_name,
                code: phased_code,
                ecef_m: phased_reference_ecef,
            },
        )?;
        let mut cw_11 = create_cor_writers(
            write_acf_cor,
            &ant1_name_file,
            &ant1_name_file,
            CorStation {
                name: a1_name,
                code: a1_code,
                ecef_m: ant1_ecef,
            },
            CorStation {
                name: a1_name,
                code: a1_code,
                ecef_m: ant1_ecef,
            },
        )?;
        let mut cw_12 = create_cor_writers(
            write_xcf_cor,
            &ant1_name_file,
            &ant2_name_file,
            CorStation {
                name: a1_name,
                code: a1_code,
                ecef_m: ant1_ecef,
            },
            CorStation {
                name: a2_name,
                code: a2_code,
                ecef_m: ant2_ecef,
            },
        )?;
        let mut cw_22 = create_cor_writers(
            write_acf_cor,
            &ant2_name_file,
            &ant2_name_file,
            CorStation {
                name: a2_name,
                code: a2_code,
                ecef_m: ant2_ecef,
            },
            CorStation {
                name: a2_name,
                code: a2_code,
                ecef_m: ant2_ecef,
            },
        )?;
        let mut cor_write_buf_ph = vec![Complex::new(0.0_f32, 0.0_f32); inband_bins];
        let mut cor_write_buf_11 = vec![Complex::new(0.0_f32, 0.0_f32); inband_bins];
        let mut cor_write_buf_12 = vec![Complex::new(0.0_f32, 0.0_f32); inband_bins];
        let mut cor_write_buf_22 = vec![Complex::new(0.0_f32, 0.0_f32); inband_bins];

        let ql_fringe_dir = fringe_interval_s.map(|_| {
            let schedule_base = args
                .schedule
                .as_ref()
                .and_then(|p| p.file_stem())
                .and_then(|s| s.to_str())
                .map(sanitize_file_token)
                .unwrap_or_else(|| "direct".to_string());
            o_dir
                .join("fringe_yicorr_ql")
                .join(schedule_base)
                .join(&process_epoch_tag)
        });
        let mut ql_fringe_acc = ql_fringe_dir
            .as_ref()
            .map(|_| vec![Complex::new(0.0_f64, 0.0_f64); cor_bins_total]);
        let mut ql_fringe_frames = 0usize;
        let mut ql_fringe_start_offset_s = 0.0_f64;
        let mut ql_fringe_index = 0usize;

        // auto_pipeline_depth() is already at least two.  Keep an explicit
        // --pipeline-depth value unchanged, including a requested depth of one.
        let prefetch_depth = io_pipeline_depth;
        let paired_ready_capacity = if args.usb {
            prefetch_depth.saturating_sub(1).max(1)
        } else {
            prefetch_depth
        };
        println!(
            "[info] I/O prefetch: process-window reader enabled (pipeline={} chunks, paired-ready={} chunks)",
            prefetch_depth, paired_ready_capacity
        );
        if args.usb {
            println!(
                "[info] Input reader mode: USB-attached storage concurrent readers (one per antenna; unrelated to USB signal sideband)"
            );
        }
        let read_align_first = compute_frame_delay_entry(0, &delay_cfg, d_seek);
        let read_align_last_frame = total_f.saturating_sub(1);
        let read_align_last = compute_frame_delay_entry(read_align_last_frame, &delay_cfg, d_seek);
        let read_align_drift_samples =
            (read_align_last.full_rel_s - read_align_first.full_rel_s).abs() * fs;
        let sector_readalign_threshold_samples = ((fft_len as f64) * 0.25).max(128.0);
        let use_sector_readalign = read_align_drift_samples > sector_readalign_threshold_samples;
        println!(
            "[info] Read-align mode: {}{} (full_rel drift {:+.3} sample over process, threshold {:.1} sample)",
            if use_sector_readalign { "adaptive-sector" } else { "fixed-process" },
            if fx_start_cumulative_seek { "; start-cumulative integer seek" } else { "" },
            (read_align_last.full_rel_s - read_align_first.full_rel_s) * fs,
            sector_readalign_threshold_samples
        );
        if fx_integer_delay_enabled() {
            println!(
                "[info] FX delay correction: integer sample tracking enabled; phase slope uses fractional delay only"
            );
        } else {
            println!(
                "[warn] FX delay correction: legacy mode, full residual delay is left in post-FFT phase slope (YI_FX_INT_DELAY=0)"
            );
        }

        let sec_counts_for_read = sec_counts.clone();
        let (tx_sec, rx_sec) = mpsc::sync_channel::<PairedSectorRead>(paired_ready_capacity);
        let synth_produced_chunks = Arc::new(AtomicUsize::new(0));
        let synth_produced_bytes = Arc::new(AtomicU64::new(0));
        let synth_consumed_chunks = Arc::new(AtomicUsize::new(0));
        let (a1_read, a2_read) = (a1p.clone(), a2p.clone());
        if fx_integer_delay_enabled() {
            println!(
                "[info] FX delay correction: fixed read-align + per-frame integer sample shift enabled (YI_FX_INT_DELAY=1)"
            );
        } else {
            println!(
                "[warn] FX delay correction: legacy mode, full residual delay is left in post-FFT phase slope (YI_FX_INT_DELAY=0)"
            );
        }
        let fx_delay_phase_offset_samples = fx_delay_phase_offset_samples();
        if fx_delay_phase_offset_samples != 0.0 {
            println!(
                "[info] FX user delay offset: {:+.6} sample added to post-FFT phase-delay correction only",
                fx_delay_phase_offset_samples
            );
        }
        let mut sector_read_starts = Vec::with_capacity(sec_counts_for_read.len());
        let mut sector_read_sample_counts = Vec::with_capacity(sec_counts_for_read.len());
        let mut sector_d_seeks = Vec::with_capacity(sec_counts_for_read.len());
        let mut sector_start_frame = 0usize;
        for &nf in &sec_counts_for_read {
            // Track the sector-level read anchor using the delay at the sector midpoint.
            // The remaining per-frame residual is split into integer + fractional
            // parts relative to this sector d_seek.
            let sector_ref_frame = sector_start_frame + nf / 2;
            let sector_ref =
                compute_frame_delay_entry(sector_ref_frame, &delay_cfg, delay_cfg.d_seek);

            // Keep the physical read branch and the delay-model branch
            // consistent.  In fixed-process mode, both must remain tied to the
            // initial process d_seek.  Only adaptive-sector read-align is
            // allowed to move the raw/decode window and update the sector
            // d_seek.  Otherwise the read branch silently changes at the
            // +/-0.5 sample boundaries even though the log says fixed-process.
            let sector_d_seek_samples = if use_sector_readalign {
                (sector_ref.full_rel_s * fs).round() as i128
            } else {
                d_seek_samples
            };
            let delta_samples = if use_sector_readalign {
                sector_d_seek_samples - d_seek_samples
            } else {
                0
            };

            let common_samples = sector_start_frame as i128 * fft_len as i128;
            let mut sector_s1 = start_s1_actual_samples as i128 + common_samples;
            let mut sector_s2 = start_s2_actual_samples as i128 + common_samples;
            if use_sector_readalign {
                if residual_on_ant2 {
                    sector_s2 += delta_samples;
                } else {
                    sector_s1 -= delta_samples;
                }
            }
            if sector_s1 < 0 || sector_s2 < 0 {
                return Err(format!(
                    "negative sector read start: ant1={} ant2={} sector_frame={}",
                    sector_s1, sector_s2, sector_start_frame
                )
                .into());
            }
            let logical_s1 = sector_s1 as u64;
            let logical_s2 = sector_s2 as u64;

            // Preserve the native 32-bit shuffle boundary for every adaptive
            // sector as well. The delay model uses the actual aligned read
            // separation below, so the discarded prefix remains corrected.
            let read_s1 = logical_s1 / samples_per_word1 * samples_per_word1;
            let read_s2 = logical_s2 / samples_per_word2 * samples_per_word2;
            let payload_samples = nf as u64 * fft_len as u64;
            // One trailing native word supplies the positive look-ahead used
            // when the residual integer delay selects samples after the
            // aligned read origin, including the final frame of a sector.
            let read_n1 = payload_samples + samples_per_word1;
            let read_n2 = payload_samples + samples_per_word2;
            sector_read_starts.push((read_s1, read_s2));
            sector_read_sample_counts.push((read_n1, read_n2));
            let actual_sector_d_seek_samples = read_s2 as i128 - read_s1 as i128;
            sector_d_seeks.push(actual_sector_d_seek_samples as f64 / fs);
            sector_start_frame += nf;
        }
        let sector_read_starts_rd = sector_read_starts.clone();
        let sector_read_sample_counts_rd = sector_read_sample_counts.clone();
        let synth_produced_chunks_rd = Arc::clone(&synth_produced_chunks);
        let synth_produced_bytes_rd = Arc::clone(&synth_produced_bytes);
        let reader_pipeline = if args.usb {
            let sector_count = sec_counts_for_read.len();
            let ant1_sector_reads = sector_read_starts_rd
                .iter()
                .zip(sector_read_sample_counts_rd.iter())
                .map(|(&(start1, _), &(count1, _))| (start1, count1))
                .collect::<Vec<_>>();
            let ant2_sector_reads = sector_read_starts_rd
                .iter()
                .zip(sector_read_sample_counts_rd.iter())
                .map(|(&(_, start2), &(_, count2))| (start2, count2))
                .collect::<Vec<_>>();

            // Rendezvous channels prevent either station reader from queuing
            // an additional sector independently of the pairer.
            let (tx_ant1, rx_ant1) = mpsc::sync_channel::<StationSectorRead>(0);
            let (tx_ant2, rx_ant2) = mpsc::sync_channel::<StationSectorRead>(0);
            let ant1_handle = thread::spawn(move || {
                run_usb_station_reader(
                    a1_read,
                    "ant1",
                    bit1,
                    ant1_sector_reads,
                    reader_core,
                    tx_ant1,
                );
            });
            let ant2_handle = thread::spawn(move || {
                run_usb_station_reader(a2_read, "ant2", bit2, ant2_sector_reads, None, tx_ant2);
            });
            let pairer_handle = thread::spawn(move || {
                for expected_sector in 0..sector_count {
                    let (sector1, b1) = match rx_ant1.recv() {
                        Ok(Ok(v)) => v,
                        Ok(Err(e)) => {
                            let _ = tx_sec.send(Err(e));
                            return;
                        }
                        Err(e) => {
                            let _ = tx_sec.send(Err(format!(
                                "USB ant1 reader channel failed at sector {}: {e}",
                                expected_sector + 1
                            )));
                            return;
                        }
                    };
                    let (sector2, b2) = match rx_ant2.recv() {
                        Ok(Ok(v)) => v,
                        Ok(Err(e)) => {
                            let _ = tx_sec.send(Err(e));
                            return;
                        }
                        Err(e) => {
                            let _ = tx_sec.send(Err(format!(
                                "USB ant2 reader channel failed at sector {}: {e}",
                                expected_sector + 1
                            )));
                            return;
                        }
                    };
                    if sector1 != expected_sector || sector2 != expected_sector {
                        let _ = tx_sec.send(Err(format!(
                            "USB reader sector pairing mismatch: expected={} ant1={} ant2={}",
                            expected_sector + 1,
                            sector1 + 1,
                            sector2 + 1
                        )));
                        return;
                    }
                    let chunk_bytes = (b1.len() + b2.len()) as u64;
                    if tx_sec.send(Ok((expected_sector, b1, b2))).is_err() {
                        return;
                    }
                    synth_produced_bytes_rd.fetch_add(chunk_bytes, Ordering::Relaxed);
                    synth_produced_chunks_rd.fetch_add(1, Ordering::Relaxed);
                }
            });
            SynthReaderPipeline::new(
                rx_sec,
                vec![
                    ("USB sector pairer", pairer_handle),
                    ("USB ant1 reader", ant1_handle),
                    ("USB ant2 reader", ant2_handle),
                ],
            )
        } else {
            let reader_handle = thread::spawn(move || {
                if let Some(core) = reader_core {
                    let _ = affinity::set_current_thread_core(core);
                }
                let mut pr1 = match PackedSampleReader::open(&a1_read, 0, 0) {
                    Ok(v) => v,
                    Err(e) => {
                        let _ = tx_sec.send(Err(format!(
                            "failed to open {} for sector reader: {e}",
                            a1_read.display()
                        )));
                        return;
                    }
                };
                let mut pr2 = match PackedSampleReader::open(&a2_read, 0, 0) {
                    Ok(v) => v,
                    Err(e) => {
                        let _ = tx_sec.send(Err(format!(
                            "failed to open {} for sector reader: {e}",
                            a2_read.display()
                        )));
                        return;
                    }
                };
                for (sector_idx, _nf) in sec_counts_for_read.into_iter().enumerate() {
                    let (read_n1, read_n2) = sector_read_sample_counts_rd[sector_idx];
                    let mut b1 = vec![0u8; ((read_n1 as usize * bit1) + 7) / 8];
                    let mut b2 = vec![0u8; ((read_n2 as usize * bit2) + 7) / 8];
                    let (sector_s1, sector_s2) = sector_read_starts_rd[sector_idx];
                    let bits1 = sector_s1 * bit1 as u64;
                    let bits2 = sector_s2 * bit2 as u64;
                    if let Err(e) = pr1.seek_to(bits1 / 8, (bits1 % 8) as u8) {
                        let _ = tx_sec.send(Err(format!(
                            "failed to seek {} at sample {}: {e}",
                            a1_read.display(),
                            sector_s1
                        )));
                        return;
                    }
                    if let Err(e) = pr2.seek_to(bits2 / 8, (bits2 % 8) as u8) {
                        let _ = tx_sec.send(Err(format!(
                            "failed to seek {} at sample {}: {e}",
                            a2_read.display(),
                            sector_s2
                        )));
                        return;
                    }
                    if let Err(e) = pr1.read_packed_with_padding(&mut b1) {
                        let _ = tx_sec.send(Err(format!(
                            "failed reading ant1 input at sample {}: {e}",
                            sector_s1
                        )));
                        return;
                    }
                    if let Err(e) = pr2.read_packed_with_padding(&mut b2) {
                        let _ = tx_sec.send(Err(format!(
                            "failed reading ant2 input at sample {}: {e}",
                            sector_s2
                        )));
                        return;
                    }
                    let chunk_bytes = (b1.len() + b2.len()) as u64;
                    if tx_sec.send(Ok((sector_idx, b1, b2))).is_err() {
                        return;
                    }
                    synth_produced_bytes_rd.fetch_add(chunk_bytes, Ordering::Relaxed);
                    synth_produced_chunks_rd.fetch_add(1, Ordering::Relaxed);
                }
            });
            SynthReaderPipeline::new(rx_sec, vec![("serial sector reader", reader_handle)])
        };
        let need_phased_products = write_raw || write_phased_cor || plot_phased;
        let need_acf_products = write_acf_cor;
        let need_xcf_products = write_xcf_cor;
        let fft_peak_debug = fft_peak_dbg_enabled();
        let normal_corr_kernel = NormalCorrKernel::from_output_grid(out_grid);
        let need_grid_buffers = pulsar_runtime.is_some() || fft_peak_debug;
        let overlap_len = ba.a1e - ba.a1s;
        // Frequency-grid maps are constant for the whole process.  Building
        // them inside every output sector is especially expensive for very
        // large FFTs, where a map pair can occupy many MiB.
        let normal_grid_maps = if write_raw {
            None
        } else {
            let half = fft_len / 2 + 1;
            let station_offset1 = station_grid_origin_offset_bins(a1_name, fs, fft_len);
            let station_offset2 =
                station_grid_origin_offset_bins(a2_name, fs, fft_len) + ant2_grid_extra_offset();
            let grid_map1 =
                build_real_fft_xml_grid_map(fft_len, rotation_bins1, station_offset1, half);
            let grid_map2 =
                build_real_fft_xml_grid_map(fft_len, rotation_bins2, station_offset2, half);
            let direct_src1 = contiguous_real_fft_src_start(&grid_map1, ba.a1s, overlap_len);
            let direct_src2 = contiguous_real_fft_src_start(&grid_map2, ba.a2s, overlap_len);
            Some((grid_map1, grid_map2, direct_src1, direct_src2))
        };
        if let Some((_, _, direct_src1, direct_src2)) = normal_grid_maps.as_ref() {
            match (direct_src1, direct_src2) {
                (Some(src1), Some(src2)) => println!(
                    "[info] Correlation grid path: direct contiguous FFT bins (src starts {}, {})",
                    src1, src2
                ),
                _ => println!("[info] Correlation grid path: general mapped/Hermitian FFT bins"),
            }
        }
        let acf_overlap_only = false;
        let mut acc_ph_total = if need_phased_products {
            Some(vec![0.0; fft_len / 2 + 1])
        } else {
            None
        };
        let mut acc_11_total = vec![0.0; fft_len / 2 + 1];
        let mut acc_22_total = vec![0.0; fft_len / 2 + 1];
        if need_xcf_products {
            println!(
                "[info] XCF phase bins: a1=[{}..{}) a2=[{}..{}) rotation_bins=({}, {}) phase_start_used_raw=({}, {}) phase_start_xml=({}, {}) overlap_bins={}",
                ba.a1s,
                ba.a1e,
                ba.a2s,
                ba.a2s + ba.a1e.saturating_sub(ba.a1s),
                rotation_bins1,
                rotation_bins2,
                ba.a1s as isize - rotation_bins1,
                ba.a2s as isize - rotation_bins2,
                ba.a1s,
                ba.a2s,
                ba.a1e.saturating_sub(ba.a1s)
            );
        }
        let synth_stats_start = Instant::now();
        let mut synth_read_bytes_total: u64 = 0;
        let mut synth_queue_hwm = 0usize;
        let mut timing_delay_s = 0.0_f64;
        let mut timing_recv_wait_s = 0.0_f64;
        let mut timing_compute_s = 0.0_f64;
        let mut timing_output_s = 0.0_f64;
        let timing_sample_stride = synth_timing_sample_stride();
        let mut timing_sample_decode_s = 0.0_f64;
        let mut timing_sample_fft_s = 0.0_f64;
        let mut timing_sample_accum_s = 0.0_f64;
        let mut timing_sample_frames = 0usize;
        for (si, &nf) in sec_counts.iter().enumerate() {
            let sector_d_seek = sector_d_seeks[si];
            let (sector_sample_start1, sector_sample_start2) = sector_read_starts[si];
            let timing_delay_start = Instant::now();
            let frame_delays: Vec<FrameDelayEntry> = (0..nf)
                .map(|i| compute_frame_delay_entry(emitted + i, &delay_cfg, sector_d_seek))
                .collect();
            timing_delay_s += timing_delay_start.elapsed().as_secs_f64();
            if args.debug && need_xcf_products {
                print_delay_debug_samples(
                    &format!("delay sector {}", si + 1),
                    emitted,
                    &frame_delays,
                    fs,
                );
                if !frame_delays.is_empty() {
                    let mid_idx = frame_delays.len() / 2;
                    let sample_refs = [
                        ("start", 0usize),
                        ("mid", mid_idx),
                        ("end", frame_delays.len() - 1),
                    ];
                    for (pos, idx) in sample_refs {
                        let d = frame_delays[idx];
                        let fr_mix = d.fr_lo1 * d.fr_lo2.conj();
                        let fr_mix_deg = (fr_mix.im as f64).atan2(fr_mix.re as f64).to_degrees();
                        println!(
                            "[info] XCF phase reference sector {} {}: frame={} t_mid={:.9e}s f_ref_xml={:.6}MHz fr_mix={:+.6}deg raw_phase_bins=({}, {}) xml_phase_bins=({}, {})",
                            si + 1,
                            pos,
                            emitted + idx,
                            d.t_mid_s,
                            obs_mhz,
                            fr_mix_deg,
                            ba.a1s as isize - rotation_bins1,
                            ba.a2s as isize - rotation_bins2,
                            ba.a1s,
                            ba.a2s
                        );
                    }
                }
            }
            if args.debug && !need_xcf_products {
                print_delay_debug_samples(
                    &format!("sector {}", si + 1),
                    emitted,
                    &frame_delays,
                    fs,
                );
            }
            if let Some(dw) = dbg_writer.as_mut() {
                writeln!(
                    dw,
                    "[sector {}] frames={} start_frame={}",
                    si + 1,
                    nf,
                    emitted
                )?;
                for i in 0..nf {
                    let d = frame_delays[i];
                    writeln!(
                        dw,
                        "  frame={} t_mid={:.9e} full_rel_samples={:.9} residual_samples={:.9} tau1_samples={:.9} tau2_samples={:.9} int1={} int2={} frac1_samples={:.9} frac2_samples={:.9}",
                        emitted + i,
                        d.t_mid_s,
                        d.full_rel_s * fs,
                        d.residual_s * fs,
                        d.tau1_s * fs,
                        d.tau2_s * fs,
                        d.int1,
                        d.int2,
                        d.frac1 * fs,
                        d.frac2 * fs
                    )?;
                }
            }
            let timing_recv_start = Instant::now();
            let (received_sector, raw1_vec, raw2_vec) = match reader_pipeline.recv() {
                Ok(Ok(v)) => v,
                Ok(Err(e)) => {
                    return Err(std::io::Error::other(format!("synth reader error: {e}")).into())
                }
                Err(e) => {
                    return Err(
                        std::io::Error::other(format!("synth reader channel error: {e}")).into(),
                    )
                }
            };
            if received_sector != si {
                return Err(std::io::Error::other(format!(
                    "synth reader sector order mismatch: expected={} received={}",
                    si + 1,
                    received_sector + 1
                ))
                .into());
            }
            timing_recv_wait_s += timing_recv_start.elapsed().as_secs_f64();
            synth_consumed_chunks.fetch_add(1, Ordering::Relaxed);
            let raw1: &[u8] = &raw1_vec;
            let raw2: &[u8] = &raw2_vec;
            synth_read_bytes_total += (raw1.len() + raw2.len()) as u64;
            let produced = synth_produced_chunks.load(Ordering::Relaxed);
            let consumed = synth_consumed_chunks.load(Ordering::Relaxed);
            let queue_fill = produced.saturating_sub(consumed);
            if queue_fill > synth_queue_hwm {
                synth_queue_hwm = queue_fill;
            }
            let consumed = synth_consumed_chunks.load(Ordering::Relaxed);
            let queue_fill = produced.saturating_sub(consumed);
            if queue_fill > synth_queue_hwm {
                synth_queue_hwm = queue_fill;
            }
            let sector_failures = AtomicUsize::new(0);
            let process_frame =
                |i: usize,
                 out_f: Option<&mut [u8]>|
                 -> Option<(Vec<f64>, Vec<f64>, Vec<Complex<f64>>, Vec<f64>)> {
                    let d = frame_delays[i];
                    let (mut f1, mut f2, mut s1, mut s2) = (
                        vec![0.0_f32; fft_len],
                        vec![0.0_f32; fft_len],
                        vec![Complex::new(0.0_f32, 0.0_f32); fft_len / 2 + 1],
                        vec![Complex::new(0.0_f32, 0.0_f32); fft_len / 2 + 1],
                    );
                    let mut dw1 = DecodeWindowScratch::new();
                    let mut dw2 = DecodeWindowScratch::new();
                    if decode_shifted_frame_from_chunk(
                        raw1,
                        sector_sample_start1,
                        i,
                        fft_len,
                        bit1,
                        samples_per_word1,
                        &dp1,
                        lsb1,
                        d.int1,
                        &mut f1,
                        &mut dw1,
                    )
                    .is_err()
                    {
                        if let Some(out_f) = out_f {
                            out_f.fill(0);
                        }
                        sector_failures.fetch_add(1, Ordering::Relaxed);
                        return None;
                    }
                    if decode_shifted_frame_from_chunk(
                        raw2,
                        sector_sample_start2,
                        i,
                        fft_len,
                        bit2,
                        samples_per_word2,
                        &dp2,
                        lsb2,
                        d.int2,
                        &mut f2,
                        &mut dw2,
                    )
                    .is_err()
                    {
                        if let Some(out_f) = out_f {
                            out_f.fill(0);
                        }
                        sector_failures.fetch_add(1, Ordering::Relaxed);
                        return None;
                    }
                    // The integer-delay sign convention is opposite to the sample-window
                    // displacement used by the decoded FFT frame.
                    if helper.forward_r2c_process(&mut f1, &mut s1).is_err() {
                        if let Some(out_f) = out_f {
                            out_f.fill(0);
                        }
                        sector_failures.fetch_add(1, Ordering::Relaxed);
                        return None;
                    }
                    if helper.forward_r2c_process(&mut f2, &mut s2).is_err() {
                        if let Some(out_f) = out_f {
                            out_f.fill(0);
                        }
                        sector_failures.fetch_add(1, Ordering::Relaxed);
                        return None;
                    }
                    // Keep the same reference as cross-correlation path.
                    let fr_lo1 = d.fr_lo1;
                    let fr_lo2 = d.fr_lo2;
                    let half = fft_len / 2 + 1;
                    let mut g1 = vec![Complex::new(0.0_f32, 0.0_f32); half];
                    let mut g2 = vec![Complex::new(0.0_f32, 0.0_f32); half];
                    shift_real_fft_to_xml_grid_with_extra_offset(
                        &s1,
                        &mut g1,
                        fft_len,
                        rotation_bins1,
                        station_grid_origin_offset_bins(a1_name, fs, fft_len),
                    );
                    shift_real_fft_to_xml_grid_with_extra_offset(
                        &s2,
                        &mut g2,
                        fft_len,
                        rotation_bins2,
                        station_grid_origin_offset_bins(a2_name, fs, fft_len)
                            + ant2_grid_extra_offset(),
                    );

                    let mut phased_pow = vec![0.0; half];
                    let mut p11 = vec![0.0; half];
                    let mut p12 = vec![Complex::new(0.0_f64, 0.0_f64); half];
                    let mut p22 = vec![0.0; half];

                    let mut s1_aligned = vec![Complex::new(0.0_f32, 0.0_f32); half];
                    let mut s2_aligned = vec![Complex::new(0.0_f32, 0.0_f32); half];
                    match out_grid {
                        OutputGrid::Ant1 => {
                            if need_xcf_products || need_acf_products || need_phased_products {
                                let (mut phase1, step1) = antenna_phase_start_and_step(
                                    df_hz,
                                    d.frac1,
                                    ba.a1s as isize - rotation_bins1,
                                );
                                let (mut phase2, step2) = antenna_phase_start_and_step(
                                    df_hz,
                                    d.frac2,
                                    ba.a2s as isize - rotation_bins2,
                                );
                                for k in 0..(ba.a1e - ba.a1s) {
                                    let i1 = ba.a1s + k;
                                    let i2 = ba.a2s + k;
                                    s1_aligned[i1] = g1[i1] * fr_lo1 * phase1;
                                    s2_aligned[i1] = g2[i2] * fr_lo2 * phase2;
                                    phase1 *= step1;
                                    phase2 *= step2;
                                }
                            }
                            let s1c = &s1_aligned;
                            if need_phased_products {
                                let mut cb = vec![Complex::new(0.0_f32, 0.0_f32); fft_len / 2 + 1];
                                for k in 0..cb.len() {
                                    cb[k] = s1c[k] * (w1 as f32) + s2_aligned[k] * (w2 as f32);
                                }
                                phased_pow =
                                    cb.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                                if let Some(out_f) = out_f {
                                    cb[0].im = 0.0_f32;
                                    if fft_len % 2 == 0 {
                                        cb[fft_len / 2].im = 0.0_f32;
                                    }
                                    let mut out_t = vec![0.0_f32; fft_len];
                                    if helper.inverse_c2r_process(&mut cb, &mut out_t).is_err() {
                                        out_f.fill(0);
                                        sector_failures.fetch_add(1, Ordering::Relaxed);
                                        return None;
                                    }
                                    if output_lsb_restore {
                                        for odd in out_t.iter_mut().skip(1).step_by(2) {
                                            *odd = -*odd;
                                        }
                                    }
                                    let mut tmp_enc = Vec::new();
                                    if quantise_frame(
                                        &out_t,
                                        bit_out,
                                        &levels_out,
                                        &output_shuffle,
                                        &mut tmp_enc,
                                    )
                                    .is_err()
                                    {
                                        out_f.fill(0);
                                        sector_failures.fetch_add(1, Ordering::Relaxed);
                                        return None;
                                    }
                                    out_f.copy_from_slice(&tmp_enc);
                                }
                            }
                            if need_acf_products {
                                p11 = s1c.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                                p22 = s2_aligned
                                    .iter()
                                    .map(|c| c.norm_sqr() as f64)
                                    .collect::<Vec<_>>();
                            }
                            if need_xcf_products {
                                p12 = s1c
                                    .iter()
                                    .zip(s2_aligned.iter())
                                    .map(|(z1, z2)| {
                                        let v = *z1 * z2.conj();
                                        Complex::new(v.re as f64, v.im as f64)
                                    })
                                    .collect::<Vec<_>>();
                            }
                        }
                        OutputGrid::Ant2 => {
                            if need_xcf_products || need_acf_products || need_phased_products {
                                let (mut phase1, step1) = antenna_phase_start_and_step(
                                    df_hz,
                                    d.frac1,
                                    ba.a1s as isize - rotation_bins1,
                                );
                                let (mut phase2, step2) = antenna_phase_start_and_step(
                                    df_hz,
                                    d.frac2,
                                    ba.a2s as isize - rotation_bins2,
                                );
                                for k in 0..(ba.a1e - ba.a1s) {
                                    let i1 = ba.a1s + k;
                                    let i2 = ba.a2s + k;
                                    s1_aligned[i2] = g1[i1] * fr_lo1 * phase1;
                                    s2_aligned[i2] = g2[i2] * fr_lo2 * phase2;
                                    phase1 *= step1;
                                    phase2 *= step2;
                                }
                            }
                            let s2c = &s2_aligned;
                            if need_phased_products {
                                let mut cb = vec![Complex::new(0.0_f32, 0.0_f32); fft_len / 2 + 1];
                                for k in 0..cb.len() {
                                    cb[k] = s1_aligned[k] * (w1 as f32) + s2c[k] * (w2 as f32);
                                }
                                phased_pow =
                                    cb.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                                if let Some(out_f) = out_f {
                                    cb[0].im = 0.0_f32;
                                    if fft_len % 2 == 0 {
                                        cb[fft_len / 2].im = 0.0_f32;
                                    }
                                    let mut out_t = vec![0.0_f32; fft_len];
                                    if helper.inverse_c2r_process(&mut cb, &mut out_t).is_err() {
                                        out_f.fill(0);
                                        sector_failures.fetch_add(1, Ordering::Relaxed);
                                        return None;
                                    }
                                    if output_lsb_restore {
                                        for odd in out_t.iter_mut().skip(1).step_by(2) {
                                            *odd = -*odd;
                                        }
                                    }
                                    let mut tmp_enc = Vec::new();
                                    if quantise_frame(
                                        &out_t,
                                        bit_out,
                                        &levels_out,
                                        &output_shuffle,
                                        &mut tmp_enc,
                                    )
                                    .is_err()
                                    {
                                        out_f.fill(0);
                                        sector_failures.fetch_add(1, Ordering::Relaxed);
                                        return None;
                                    }
                                    out_f.copy_from_slice(&tmp_enc);
                                }
                            }
                            if need_acf_products {
                                p11 = s1_aligned
                                    .iter()
                                    .map(|c| c.norm_sqr() as f64)
                                    .collect::<Vec<_>>();
                                p22 = s2c.iter().map(|c| c.norm_sqr() as f64).collect::<Vec<_>>();
                            }
                            if need_xcf_products {
                                p12 = s1_aligned
                                    .iter()
                                    .zip(s2c.iter())
                                    .map(|(z1, z2)| {
                                        let v = *z1 * z2.conj();
                                        Complex::new(v.re as f64, v.im as f64)
                                    })
                                    .collect::<Vec<_>>();
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
            let reduce_acc =
                |mut acc1: (Vec<f64>, Vec<f64>, Vec<Complex<f64>>, Vec<f64>),
                 acc2: (Vec<f64>, Vec<f64>, Vec<Complex<f64>>, Vec<f64>)| {
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
            let timing_compute_start = Instant::now();
            let (batch_ph, batch_11, batch_12, batch_22, batch_fold) = if write_raw {
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
                (acc.0, acc.1, acc.2, acc.3, None)
            } else {
                struct ThreadAccum {
                    acc_ph: Vec<f64>,
                    acc_11: Vec<f64>,
                    acc_12: Vec<Complex<f64>>,
                    acc_22: Vec<f64>,
                    fold: Option<FoldAccum>,
                    f1: Vec<f32>,
                    f2: Vec<f32>,
                    s1: Vec<Complex<f32>>,
                    s2: Vec<Complex<f32>>,
                    g1: Vec<Complex<f32>>,
                    g2: Vec<Complex<f32>>,
                    fft_scratch: FftScratch,
                    dw1: DecodeWindowScratch,
                    dw2: DecodeWindowScratch,
                    timing_decode_s: f64,
                    timing_fft_s: f64,
                    timing_accum_s: f64,
                    timing_samples: usize,
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
                    fold: pulsar_runtime
                        .as_ref()
                        .map(|rt| FoldAccum::new(rt.bins(), half, need_phased_products)),
                    f1: vec![0.0_f32; fft_len],
                    f2: vec![0.0_f32; fft_len],
                    s1: vec![Complex::new(0.0_f32, 0.0_f32); half],
                    s2: vec![Complex::new(0.0_f32, 0.0_f32); half],
                    g1: if need_grid_buffers {
                        vec![Complex::new(0.0_f32, 0.0_f32); half + 1]
                    } else {
                        Vec::new()
                    },
                    g2: if need_grid_buffers {
                        vec![Complex::new(0.0_f32, 0.0_f32); half + 1]
                    } else {
                        Vec::new()
                    },
                    fft_scratch: helper.make_scratch(),
                    dw1: DecodeWindowScratch::new(),
                    dw2: DecodeWindowScratch::new(),
                    timing_decode_s: 0.0,
                    timing_fft_s: 0.0,
                    timing_accum_s: 0.0,
                    timing_samples: 0,
                };
                let frames_per_job = (nf / (cpu_threads.saturating_mul(8)).max(1)).clamp(128, 2048);
                let chunk_starts: Vec<usize> = (0..nf).step_by(frames_per_job).collect();
                let chunk_abs_start1 = sector_sample_start1;
                let chunk_abs_start2 = sector_sample_start2;
                let (grid_map1, grid_map2, direct_src1, direct_src2) = normal_grid_maps
                    .as_ref()
                    .expect("normal correlation grid maps");
                let direct_src1 = *direct_src1;
                let direct_src2 = *direct_src2;
                let station_offset1 = station_grid_origin_offset_bins(a1_name, fs, fft_len);
                let station_offset2 = station_grid_origin_offset_bins(a2_name, fs, fft_len)
                    + ant2_grid_extra_offset();
                let mut out = chunk_starts
                    .into_par_iter()
                    .map(|start| {
                        let mut st = init();
                        let end = (start + frames_per_job).min(nf);
                        for i in start..end {
                            let d = frame_delays[i];
                            let timing_frame_sampled = i % timing_sample_stride == 0;
                            let timing_decode_start = timing_frame_sampled.then(Instant::now);
                            if decode_shifted_frame_from_chunk(
                                raw1,
                                chunk_abs_start1,
                                i,
                                fft_len,
                                bit1,
                                samples_per_word1,
                                &dp1,
                                lsb1,
                                d.int1,
                                &mut st.f1,
                                &mut st.dw1,
                            )
                            .is_err()
                            {
                                sector_failures.fetch_add(1, Ordering::Relaxed);
                                continue;
                            }
                            if decode_shifted_frame_from_chunk(
                                raw2,
                                chunk_abs_start2,
                                i,
                                fft_len,
                                bit2,
                                samples_per_word2,
                                &dp2,
                                lsb2,
                                d.int2,
                                &mut st.f2,
                                &mut st.dw2,
                            )
                            .is_err()
                            {
                                sector_failures.fetch_add(1, Ordering::Relaxed);
                                continue;
                            }
                            if let Some(t) = timing_decode_start {
                                st.timing_decode_s += t.elapsed().as_secs_f64();
                            }
                            let timing_fft_start = timing_frame_sampled.then(Instant::now);
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
                            if let Some(t) = timing_fft_start {
                                st.timing_fft_s += t.elapsed().as_secs_f64();
                            }
                            let timing_accum_start = timing_frame_sampled.then(Instant::now);
                            let fr_lo1 = d.fr_lo1;
                            let fr_lo2 = d.fr_lo2;
                            let fr_mix = fr_lo1 * fr_lo2.conj();
                            if need_grid_buffers {
                                shift_real_fft_to_xml_grid_with_extra_offset(
                                    &st.s1,
                                    &mut st.g1,
                                    fft_len,
                                    rotation_bins1,
                                    station_offset1,
                                );
                                shift_real_fft_to_xml_grid_with_extra_offset(
                                    &st.s2,
                                    &mut st.g2,
                                    fft_len,
                                    rotation_bins2,
                                    station_offset2,
                                );

                                if fft_peak_debug && i < fft_peak_dbg_max_frames() {
                                    print_fft_peak_dbg(
                                        "raw_to_grid",
                                        i,
                                        a1_name,
                                        &st.s1,
                                        &st.g1,
                                        rotation_bins1,
                                        0,
                                    );
                                    print_fft_peak_dbg(
                                        "raw_to_grid",
                                        i,
                                        a2_name,
                                        &st.s2,
                                        &st.g2,
                                        rotation_bins2,
                                        station_offset2,
                                    );
                                }
                            }
                            if let (Some(rt), Some(fold)) =
                                (pulsar_runtime.as_ref(), st.fold.as_mut())
                            {
                                let t_since_process_s = d.t_mid_s - total_skip_sec;
                                match out_grid {
                                    OutputGrid::Ant1 => {
                                        let phase_delay1_s = d.frac1;
                                        let phase_delay2_s = d.frac2;
                                        let (mut phase_corr, phase_step) = xcf_phase_start_and_step(
                                            df_hz,
                                            phase_delay1_s,
                                            phase_delay2_s,
                                            ba.a1s as isize - rotation_bins1,
                                            ba.a2s as isize - rotation_bins2,
                                        );
                                        for k in 0..overlap_len {
                                            let i1 = ba.a1s + k;
                                            let i2 = ba.a2s + k;
                                            let raw_xcf = st.g1[i1] * st.g2[i2].conj();
                                            let v = raw_xcf * fr_mix * phase_corr;
                                            fold.add_values(
                                                rt,
                                                t_since_process_s,
                                                i1,
                                                None,
                                                st.g1[i1].norm_sqr() as f64,
                                                Complex::new(v.re as f64, v.im as f64),
                                                st.g2[i2].norm_sqr() as f64,
                                            );
                                            phase_corr *= phase_step;
                                        }
                                    }
                                    OutputGrid::Ant2 => {
                                        let phase_delay1_s = d.frac1;
                                        let phase_delay2_s = d.frac2;
                                        let (mut phase_corr, phase_step) = xcf_phase_start_and_step(
                                            df_hz,
                                            phase_delay1_s,
                                            phase_delay2_s,
                                            ba.a1s as isize - rotation_bins1,
                                            ba.a2s as isize - rotation_bins2,
                                        );
                                        for k in 0..overlap_len {
                                            let i1 = ba.a1s + k;
                                            let i2 = ba.a2s + k;
                                            let raw_xcf = st.g1[i1] * st.g2[i2].conj();
                                            let v = raw_xcf * fr_mix * phase_corr;
                                            fold.add_values(
                                                rt,
                                                t_since_process_s,
                                                i2,
                                                None,
                                                st.g1[i1].norm_sqr() as f64,
                                                Complex::new(v.re as f64, v.im as f64),
                                                st.g2[i2].norm_sqr() as f64,
                                            );
                                            phase_corr *= phase_step;
                                        }
                                    }
                                }
                                if let Some(t) = timing_accum_start {
                                    st.timing_accum_s += t.elapsed().as_secs_f64();
                                    st.timing_samples += 1;
                                }
                                continue;
                            }

                            // Fast path: do not materialize XML-grid spectra for normal
                            // correlation.  Read the mapped FFT bin once and update ACF/XCF in
                            // the same pass to reduce memory traffic and extra 1024-bin loops.
                            let output_start = match normal_corr_kernel {
                                NormalCorrKernel::Ant1Grid => ba.a1s,
                                NormalCorrKernel::Ant2Grid => ba.a2s,
                            };
                            let direct_accumulated = if need_acf_products
                                && need_xcf_products
                                && !args.debug
                            {
                                if let (Some(src1), Some(src2)) = (direct_src1, direct_src2) {
                                    let (phase0, phase_step) = xcf_phase_start_and_step(
                                        df_hz,
                                        d.frac1,
                                        d.frac2,
                                        ba.a1s as isize - rotation_bins1,
                                        ba.a2s as isize - rotation_bins2,
                                    );
                                    match normal_corr_kernel {
                                        NormalCorrKernel::Ant1Grid => accumulate_direct_acf_xcf(
                                            &st.s1[src1..src1 + overlap_len],
                                            &st.s2[src2..src2 + overlap_len],
                                            &mut st.acc_11[ba.a1s..ba.a1s + overlap_len],
                                            &mut st.acc_12[ba.a1s..ba.a1s + overlap_len],
                                            &mut st.acc_22[ba.a1s..ba.a1s + overlap_len],
                                            fr_mix,
                                            phase0,
                                            phase_step,
                                        ),
                                        NormalCorrKernel::Ant2Grid => accumulate_direct_acf_xcf(
                                            &st.s1[src1..src1 + overlap_len],
                                            &st.s2[src2..src2 + overlap_len],
                                            &mut st.acc_11[ba.a2s..ba.a2s + overlap_len],
                                            &mut st.acc_12[ba.a2s..ba.a2s + overlap_len],
                                            &mut st.acc_22[ba.a2s..ba.a2s + overlap_len],
                                            fr_mix,
                                            phase0,
                                            phase_step,
                                        ),
                                    }
                                    true
                                } else {
                                    false
                                }
                            } else {
                                false
                            };
                            if !direct_accumulated {
                                let mut phase_corr = None;
                                if need_xcf_products {
                                    let (phase0, step) = xcf_phase_start_and_step(
                                        df_hz,
                                        d.frac1,
                                        d.frac2,
                                        ba.a1s as isize - rotation_bins1,
                                        ba.a2s as isize - rotation_bins2,
                                    );
                                    phase_corr = Some((phase0, step));
                                }
                                for k in 0..overlap_len {
                                    let i1 = ba.a1s + k;
                                    let i2 = ba.a2s + k;
                                    let output_bin = output_start + k;
                                    let z1 = mapped_real_fft_bin(&st.s1, &grid_map1, i1);
                                    let z2 = mapped_real_fft_bin(&st.s2, &grid_map2, i2);
                                    if need_acf_products {
                                        st.acc_11[output_bin] += z1.norm_sqr() as f64;
                                        st.acc_22[output_bin] += z2.norm_sqr() as f64;
                                    }
                                    if let Some((ref mut pc, step)) = phase_corr {
                                        let raw_xcf = z1 * z2.conj();
                                        let v = raw_xcf * fr_mix * *pc;
                                        maybe_dump_xcf_debug(
                                            args.debug,
                                            emitted + i,
                                            i,
                                            si,
                                            normal_corr_kernel.output_grid(),
                                            output_bin,
                                            raw_xcf,
                                            fr_mix,
                                            *pc,
                                            v,
                                            d,
                                            fs,
                                        );
                                        st.acc_12[output_bin] +=
                                            Complex::new(v.re as f64, v.im as f64);
                                        *pc *= step;
                                    }
                                }
                            }
                            if need_acf_products && !acf_overlap_only {
                                let (spectrum, map, accumulator) = match normal_corr_kernel {
                                    NormalCorrKernel::Ant1Grid => {
                                        (&st.s1, &grid_map1, &mut st.acc_11)
                                    }
                                    NormalCorrKernel::Ant2Grid => {
                                        (&st.s2, &grid_map2, &mut st.acc_22)
                                    }
                                };
                                for k in 0..ba.a1s {
                                    let z = mapped_real_fft_bin(spectrum, map, k);
                                    accumulator[k] += z.norm_sqr() as f64;
                                }
                                for k in ba.a1e..half {
                                    let z = mapped_real_fft_bin(spectrum, map, k);
                                    accumulator[k] += z.norm_sqr() as f64;
                                }
                            }
                            if let Some(t) = timing_accum_start {
                                st.timing_accum_s += t.elapsed().as_secs_f64();
                                st.timing_samples += 1;
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
                        if let (Some(a_fold), Some(b_fold)) = (a.fold.as_mut(), b.fold) {
                            a_fold.merge(b_fold);
                        }
                        a.timing_decode_s += b.timing_decode_s;
                        a.timing_fft_s += b.timing_fft_s;
                        a.timing_accum_s += b.timing_accum_s;
                        a.timing_samples += b.timing_samples;
                        a
                    });

                timing_sample_decode_s += out.timing_decode_s;
                timing_sample_fft_s += out.timing_fft_s;
                timing_sample_accum_s += out.timing_accum_s;
                timing_sample_frames += out.timing_samples;

                if acf_peak_dbg_enabled() {
                    print_acf_peak_dbg(si, a1_name, &out.acc_11, rotation_bins1, 0, fs, fft_len);
                    print_acf_peak_dbg(
                        si,
                        a2_name,
                        &out.acc_22,
                        rotation_bins2,
                        ant2_grid_extra_offset(),
                        fs,
                        fft_len,
                    );
                }

                (
                    std::mem::take(&mut out.acc_ph),
                    std::mem::take(&mut out.acc_11),
                    std::mem::take(&mut out.acc_12),
                    std::mem::take(&mut out.acc_22),
                    out.fold.take(),
                )
            };
            timing_compute_s += timing_compute_start.elapsed().as_secs_f64();
            let timing_output_start = Instant::now();
            let sec_failed = sector_failures.load(Ordering::Relaxed);
            if sec_failed > 0 {
                println!(
                    "[warn] Sector {} skipped {} frame(s) due to decode/FFT/quantize errors",
                    si + 1,
                    sec_failed
                );
            }
            let sector_start_offset_s = emitted as f64 * frame_sec;
            let sector_integ_s = nf as f64 * frame_sec;
            emitted += nf;
            if is_phased_mode {
                print!(
                    "\r[info] Synthesised sector {}/{} ({:.2}%)",
                    si + 1,
                    sec_counts.len(),
                    (emitted as f64 / total_f as f64) * 100.0
                );
            } else {
                print!(
                    "\r[info] Processed sector {}/{} ({:.2}%)",
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
            if let (Some(interval_s), Some(acc)) = (fringe_interval_s, ql_fringe_acc.as_mut()) {
                if ql_fringe_frames == 0 {
                    ql_fringe_start_offset_s = sector_start_offset_s;
                }
                for k in 0..cor_bins_total {
                    acc[k] += batch_12[k];
                }
                ql_fringe_frames += nf;
                let ql_integ_s = ql_fringe_frames as f64 * frame_sec;
                let is_last_sector = si + 1 == sec_counts.len();
                if ql_integ_s + 0.5 * frame_sec >= interval_s || is_last_sector {
                    let ql_norm = 1.0
                        / (COR_ONE_SIDED_POWER_FACTOR
                            * (level_power1 * level_power2).sqrt()
                            * (fft_len as f64).powi(2)
                            * ql_fringe_frames.max(1) as f64);
                    let ql_spec: Vec<Complex<f64>> = acc
                        .iter()
                        .map(|v| Complex::new(v.re * ql_norm, v.im * ql_norm))
                        .collect();
                    let ql_label = format!("{}_ql{:04}", cor_label, ql_fringe_index);
                    fringe::write_ql_fringe_products(
                        ql_fringe_dir.as_ref().expect("fringe output directory"),
                        &c_tag,
                        &ql_label,
                        ql_fringe_start_offset_s,
                        ql_integ_s,
                        &ql_spec,
                        obs_mhz,
                        df_hz,
                    )?;
                    acc.fill(Complex::new(0.0_f64, 0.0_f64));
                    ql_fringe_frames = 0;
                    ql_fringe_index += 1;
                }
            }
            if let Some(fold) = batch_fold.as_ref() {
                let norm_11_per_frame =
                    1.0 / (COR_ONE_SIDED_POWER_FACTOR * level_power1 * (fft_len as f64).powi(2));
                let norm_22_per_frame =
                    1.0 / (COR_ONE_SIDED_POWER_FACTOR * level_power2 * (fft_len as f64).powi(2));
                let norm_12_per_frame = 1.0
                    / (COR_ONE_SIDED_POWER_FACTOR
                        * (level_power1 * level_power2).sqrt()
                        * (fft_len as f64).powi(2));
                let norm_ph_per_frame = 1.0
                    / (COR_ONE_SIDED_POWER_FACTOR * phased_power_ref * (fft_len as f64).powi(2));
                for pbin in 0..pulsar_bin_count {
                    for ch in 0..inband {
                        let start = ch * inband_bins;
                        let end = start + inband_bins;
                        let count_max = fold.count_max(pbin, start, end);
                        if count_max == 0 {
                            continue;
                        }
                        let effective_s = count_max as f64 * frame_sec;
                        let idx = pbin * inband + ch;
                        if let Some(w) = cw_ph.get_mut(idx) {
                            let s_ph = fold.normalized_real(
                                FoldProduct::Phased,
                                pbin,
                                start,
                                end,
                                norm_ph_per_frame,
                            );
                            w.write_sector_at(
                                c_unix,
                                sector_start_offset_s,
                                effective_s as f32,
                                &s_ph,
                            )?;
                        }
                        if let Some(w) = cw_11.get_mut(idx) {
                            let s_11 = fold.normalized_real(
                                FoldProduct::P11,
                                pbin,
                                start,
                                end,
                                norm_11_per_frame,
                            );
                            w.write_sector_at(
                                c_unix,
                                sector_start_offset_s,
                                effective_s as f32,
                                &s_11,
                            )?;
                        }
                        if let Some(w) = cw_12.get_mut(idx) {
                            let s_12 = fold.normalized_complex(pbin, start, end, norm_12_per_frame);
                            w.write_sector_at(
                                c_unix,
                                sector_start_offset_s,
                                effective_s as f32,
                                &s_12,
                            )?;
                        }
                        if let Some(w) = cw_22.get_mut(idx) {
                            let s_22 = fold.normalized_real(
                                FoldProduct::P22,
                                pbin,
                                start,
                                end,
                                norm_22_per_frame,
                            );
                            w.write_sector_at(
                                c_unix,
                                sector_start_offset_s,
                                effective_s as f32,
                                &s_22,
                            )?;
                        }
                    }
                }
            } else {
                for (ch, w) in cw_ph.iter_mut().enumerate() {
                    let start = ch * inband_bins;
                    let end = start + inband_bins;
                    for (out, &v) in cor_write_buf_ph.iter_mut().zip(&batch_ph[start..end]) {
                        *out = Complex::new((v * inv_ph) as f32, 0.0);
                    }
                    w.write_sector_at(
                        c_unix,
                        sector_start_offset_s,
                        sector_integ_s as f32,
                        &cor_write_buf_ph,
                    )?;
                }
                for (ch, w) in cw_11.iter_mut().enumerate() {
                    let start = ch * inband_bins;
                    let end = start + inband_bins;
                    for (out, &v) in cor_write_buf_11.iter_mut().zip(&batch_11[start..end]) {
                        *out = Complex::new((v * inv_11) as f32, 0.0);
                    }
                    w.write_sector_at(
                        c_unix,
                        sector_start_offset_s,
                        sector_integ_s as f32,
                        &cor_write_buf_11,
                    )?;
                }
                for (ch, w) in cw_12.iter_mut().enumerate() {
                    let start = ch * inband_bins;
                    let end = start + inband_bins;
                    for (out, v) in cor_write_buf_12.iter_mut().zip(&batch_12[start..end]) {
                        *out = Complex::new((v.re * inv_12) as f32, (v.im * inv_12) as f32);
                    }
                    w.write_sector_at(
                        c_unix,
                        sector_start_offset_s,
                        sector_integ_s as f32,
                        &cor_write_buf_12,
                    )?;
                }
                for (ch, w) in cw_22.iter_mut().enumerate() {
                    let start = ch * inband_bins;
                    let end = start + inband_bins;
                    for (out, &v) in cor_write_buf_22.iter_mut().zip(&batch_22[start..end]) {
                        *out = Complex::new((v * inv_22) as f32, 0.0);
                    }
                    w.write_sector_at(
                        c_unix,
                        sector_start_offset_s,
                        sector_integ_s as f32,
                        &cor_write_buf_22,
                    )?;
                }
            }
            timing_output_s += timing_output_start.elapsed().as_secs_f64();
        }
        reader_pipeline.finish()?;
        if !args.compact_logs {
            println!();
            let elapsed = synth_stats_start.elapsed().as_secs_f64().max(1e-9);
            let avg_mib_s = (synth_read_bytes_total as f64 / (1024.0 * 1024.0)) / elapsed;
            let produced = synth_produced_chunks.load(Ordering::Relaxed);
            let consumed = synth_consumed_chunks.load(Ordering::Relaxed);
            let queue_fill = produced.saturating_sub(consumed);
            println!(
                "[info] Synth I/O reader summary: avg={:.1} MiB/s total={:.1} MiB queue_end={}/{} hwm={}/{}",
                avg_mib_s,
                synth_read_bytes_total as f64 / (1024.0 * 1024.0),
                queue_fill,
                paired_ready_capacity,
                synth_queue_hwm,
                paired_ready_capacity
            );
            let timed = timing_delay_s + timing_recv_wait_s + timing_compute_s + timing_output_s;
            let pct = |v: f64| {
                if timed > 0.0 {
                    100.0 * v / timed
                } else {
                    0.0
                }
            };
            println!(
                "[info] Synth timing summary: delay={:.3}s ({:.1}%) recv-wait={:.3}s ({:.1}%) compute={:.3}s ({:.1}%) output={:.3}s ({:.1}%) accounted={:.3}s elapsed={:.3}s",
                timing_delay_s,
                pct(timing_delay_s),
                timing_recv_wait_s,
                pct(timing_recv_wait_s),
                timing_compute_s,
                pct(timing_compute_s),
                timing_output_s,
                pct(timing_output_s),
                timed,
                elapsed
            );
            let sample_timed = timing_sample_decode_s + timing_sample_fft_s + timing_sample_accum_s;
            if timing_sample_frames > 0 && sample_timed > 0.0 {
                let sample_pct = |v: f64| 100.0 * v / sample_timed;
                println!(
                    "[info] Compute sampled breakdown: frames={} stride={} decode={:.6}s ({:.1}%) fft={:.6}s ({:.1}%) accum={:.6}s ({:.1}%) sampled_total={:.6}s",
                    timing_sample_frames,
                    timing_sample_stride,
                    timing_sample_decode_s,
                    sample_pct(timing_sample_decode_s),
                    timing_sample_fft_s,
                    sample_pct(timing_sample_fft_s),
                    timing_sample_accum_s,
                    sample_pct(timing_sample_accum_s),
                    sample_timed
                );
            }
        }
        if let Some(mut w) = wr.take() {
            w.flush()?;
            drop(w);
            let expected_raw_bytes = total_f
                .checked_mul(bpf_o)
                .ok_or("phased raw byte count overflow")?
                as u64;
            let actual_raw_bytes = std::fs::metadata(&raw_partial_path)?.len();
            if actual_raw_bytes != expected_raw_bytes {
                return Err(format!(
                    "phased raw size mismatch: wrote {} bytes, expected {}",
                    actual_raw_bytes, expected_raw_bytes
                )
                .into());
            }
            std::fs::rename(&raw_partial_path, &o_path)?;

            let meta_path = o_path.with_extension("raw.meta");
            let mut meta = BufWriter::new(File::create(&meta_path)?);
            writeln!(meta, "format=yi-phasedarray-raw-v1")?;
            writeln!(meta, "software_version={}", env!("CARGO_PKG_VERSION"))?;
            writeln!(meta, "virtual_station={}", phased_name)?;
            writeln!(meta, "raw_file={}", o_path.display())?;
            writeln!(meta, "epoch_utc={}", ep_i)?;
            writeln!(meta, "output_tag={}", c_tag)?;
            writeln!(meta, "source={}", source_name)?;
            writeln!(meta, "reference_station={}", phased_reference_name)?;
            writeln!(meta, "native_format_station={}", output_format_name)?;
            writeln!(
                meta,
                "reference_ecef_m={:.9},{:.9},{:.9}",
                phased_reference_ecef[0], phased_reference_ecef[1], phased_reference_ecef[2]
            )?;
            writeln!(meta, "sampling_hz={:.9}", fs)?;
            writeln!(meta, "observing_frequency_mhz={:.9}", obs_mhz)?;
            writeln!(meta, "rotation_hz={:.9}", output_rotation_hz)?;
            writeln!(meta, "sideband={}", output_sideband)?;
            writeln!(meta, "bit={}", bit_out)?;
            writeln!(
                meta,
                "level={}",
                levels_out
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
            writeln!(
                meta,
                "shuffle_external={}",
                output_shuffle_ext
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
            writeln!(meta, "fft={}", fft_len)?;
            writeln!(meta, "frames={}", total_f)?;
            writeln!(meta, "samples={}", total_f * fft_len)?;
            writeln!(meta, "duration_s={:.12}", total_f as f64 * frame_sec)?;
            writeln!(meta, "input_ant1={}", a1p.display())?;
            writeln!(meta, "input_ant2={}", a2p.display())?;
            writeln!(meta, "weight_ant1={:.17e}", w1)?;
            writeln!(meta, "weight_ant2={:.17e}", w2)?;
            writeln!(meta, "voltage_scale={:.17e}", beam_voltage_scale)?;
            meta.flush()?;
            println!(
                "[info] Completed phased raw: {} ({} bytes), metadata: {}",
                o_path.display(),
                actual_raw_bytes,
                meta_path.display()
            );
        }
        for w in cw_ph {
            w.finalize()?;
        }
        for w in cw_11 {
            w.finalize()?;
        }
        for w in cw_12 {
            w.finalize()?;
        }
        for w in cw_22 {
            w.finalize()?;
        }
        if let Some(mut dw) = dbg_writer {
            dw.flush()?;
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
            let freqs_obs_mhz: Vec<f64> =
                (0..spec_bins).map(|i| obs_mhz + (i as f64 * df)).collect();

            plot_series_f64_x(
                &freqs_obs_mhz,
                phased_plot,
                "Phased Auto-Spectrum (ObsRef)",
                &o_dir
                    .join(format!(
                        "{}_{}_phased_auto_spectrum.png",
                        phased_name, c_tag
                    ))
                    .to_string_lossy(),
                "Frequency (MHz)",
                "Power",
                None,
                "Auto-Spectrum",
            )?;

            let amp_ph: Vec<f64> = phased_plot.iter().map(|&v| v.sqrt()).collect();
            let amp_11: Vec<f64> = a11_plot.iter().map(|&v| v.sqrt()).collect();
            let amp_22: Vec<f64> = a22_plot.iter().map(|&v| v.sqrt()).collect();
            let l1 = format!("{a1_name} (Ref)");
            let l2 = format!("{a2_name} (Target)");
            plot_multi_series_f64_x(
                &freqs_obs_mhz,
                &[
                    (&amp_ph, &BLUE, &format!("{} (Phased)", phased_name)),
                    (&amp_11, &GREEN, &l1),
                    (&amp_22, &RED, &l2),
                ],
                "Phased Spectrum Amplitude (ObsRef)",
                &o_dir
                    .join(format!(
                        "{}_{}_phased_spectrum_amplitude.png",
                        phased_name, c_tag
                    ))
                    .to_string_lossy(),
                "Frequency (MHz)",
                "Amplitude",
                None,
            )?;

            let mut full_spec = vec![Complex::new(0.0_f32, 0.0_f32); fft_len];
            for (i, &v) in phased_plot.iter().enumerate() {
                let vf = v as f32;
                full_spec[i] = Complex::new(vf, 0.0_f32);
                if i > 0 && i < fft_len / 2 {
                    full_spec[fft_len - i] = Complex::new(vf, 0.0_f32);
                }
            }
            helper.inverse_c2c(&mut full_spec)?;
            let acf_mag: Vec<f64> = full_spec.iter().map(|c| c.re as f64).collect();
            let acf_shifted: Vec<f64> = acf_mag
                .iter()
                .cycle()
                .skip(fft_len / 2)
                .take(fft_len)
                .copied()
                .collect();
            let lags: Vec<i32> = (-(fft_len as i32 / 2)..(fft_len as i32 / 2)).collect();
            plot_series_with_x(
                &lags,
                &[(&acf_shifted, &BLUE)],
                "Phased Autocorrelation",
                &o_dir
                    .join(format!(
                        "{}_{}_phased_autocorrelation.png",
                        phased_name, c_tag
                    ))
                    .to_string_lossy(),
                "Lag (samples)",
                "ACF",
                None,
                None,
            )?;
        }
    }
    Ok(())
}
