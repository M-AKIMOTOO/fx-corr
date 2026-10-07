//! Two simultaneous, separated RF bands. Never concatenate them onto a false
//! contiguous frequency axis. Rate is a delay derivative expressed at ref_hz.
use std::f64::consts::TAU;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;

use crate::args::Args;
use crate::utils::DynError;

#[path = "multiband_diagnostics.rs"]
mod diagnostics;

/// Resolve one CX schedule into two native schedules in memory. Shared
/// station/source/process definitions occur only once; each band supplies its
/// stream and optional overrides (for example its clock). No temporary XMLs.
pub fn configure_from_xml(args: &mut Args) -> Result<(), DynError> {
    let Some(path) = args.schedule.clone() else {
        return Ok(());
    };
    if !path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("xml"))
    {
        return Ok(());
    }
    let text = std::fs::read_to_string(&path)?;
    let doc = roxmltree::Document::parse(&text)?;
    let root = doc.root_element();
    let is_tag = |n: roxmltree::Node<'_, '_>, tag: &str| {
        n.is_element() && n.tag_name().name().eq_ignore_ascii_case(tag)
    };
    let definitions: Vec<_> = root
        .children()
        .filter(|n| is_tag(*n, "multiband"))
        .collect();
    if definitions.is_empty() {
        if args.multiband_schedule.is_none()
            && (args.multiband_raw_directory.is_some()
                || args.multiband_ant1.is_some()
                || args.multiband_ant2.is_some()
                || !args.multiband_calibrator.is_empty())
        {
            return Err(
                "multiband options require a <multiband> schedule or --multiband-schedule".into(),
            );
        }
        return Ok(());
    }
    if definitions.len() != 1
        || args.multiband_schedule.is_some()
        || root.children().any(|n| is_tag(n, "stream"))
    {
        return Err("use one <multiband> element with two bands, without a root stream or --multiband-schedule".into());
    }
    let definition = definitions[0];
    let bands: Vec<_> = definition
        .children()
        .filter(|n| is_tag(*n, "band"))
        .collect();
    if bands.len() != 2 {
        return Err("<multiband> requires exactly two <band> elements".into());
    }
    let base = args
        .raw_directory
        .clone()
        .ok_or("--raw is required for a multiband XML")?;
    let mut resolved = Vec::new();
    let mut names = Vec::new();
    for band in bands {
        let name = band
            .attribute("name")
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .ok_or("multiband band/name is required")?;
        if name == "."
            || name == ".."
            || !name
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "_-.".contains(c))
            || names.contains(&name)
        {
            return Err("multiband band names must be distinct directory-safe names".into());
        }
        names.push(name);
        if band.children().filter(|n| is_tag(*n, "stream")).count() != 1 {
            return Err("each multiband band requires exactly one <stream>".into());
        }
        let raw = band
            .children()
            .find(|n| is_tag(*n, "raw-directory"))
            .and_then(|n| n.text())
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .unwrap_or(name);
        let raw = base.join(raw);
        let mut native = String::from("<schedule>\n");
        for node in root
            .children()
            .filter(|n| n.is_element() && !is_tag(*n, "multiband"))
        {
            native.push_str(&text[node.range()]);
            native.push('\n');
        }
        for node in band
            .children()
            .filter(|n| n.is_element() && !is_tag(*n, "raw-directory"))
        {
            native.push_str(&text[node.range()]);
            native.push('\n');
        }
        native.push_str("</schedule>\n");
        let meta = crate::xml::parse_xml_schedule_text(&native, None)?;
        let frequency = meta.obsfreq_mhz.ok_or("multiband frequency missing")?;
        if !frequency.is_finite() || frequency <= 0.0 {
            return Err("invalid multiband RF frequency".into());
        }
        resolved.push((frequency, raw, Arc::new(native)));
    }
    resolved.sort_by(|a, b| a.0.total_cmp(&b.0));
    let high = resolved.pop().unwrap();
    let low = resolved.pop().unwrap();
    args.raw_directory = Some(low.1);
    args.schedule_xml = Some(low.2);
    args.multiband_schedule = Some(path);
    args.multiband_schedule_xml = Some(high.2);
    if args.multiband_raw_directory.is_none() {
        args.multiband_raw_directory = Some(high.1);
    }
    if args.multiband_calibrator.is_empty() {
        if let Some(names) = definition.attribute("calibrator") {
            args.multiband_calibrator = names
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect();
        }
    }
    println!(
        "[multiband] CX XML: band1 raw={} band2 raw={} calibrators={}",
        args.raw_directory.as_ref().unwrap().display(),
        args.multiband_raw_directory.as_ref().unwrap().display(),
        args.multiband_calibrator.join(",")
    );
    Ok(())
}

#[derive(Clone, Debug)]
pub struct Solution {
    pub start_s: f64,
    pub end_s: f64,
    pub reference_s: f64,
    pub delay_s: f64,
    pub rate_hz: f64,
    pub phase_rad: f64,
    pub coherence: f64,
    pub phase_connected: bool,
}

#[derive(Clone, Debug)]
pub struct Corrections {
    pub reference_hz: f64,
    pub band_phase_rad: [f64; 2],
    pub solutions: Vec<Solution>,
}

impl Corrections {
    pub fn phase_start_and_step(
        &self,
        band: usize,
        time_s: f64,
        frequency_hz: f64,
        df_hz: f64,
    ) -> (Complex<f64>, Complex<f64>) {
        // Solutions cover the complete scan. Boundary frames select the next
        // window; the final FFT midpoint always precedes the final end time.
        let index = self
            .solutions
            .partition_point(|s| s.end_s <= time_s)
            .min(self.solutions.len() - 1);
        let s = &self.solutions[index];
        let dt = time_s - s.reference_s;
        let delay = s.delay_s + s.rate_hz / self.reference_hz * dt;
        let cycles = (frequency_hz - self.reference_hz) * s.delay_s
            + frequency_hz / self.reference_hz * s.rate_hz * dt;
        (
            Complex::from_polar(
                1.0,
                -(TAU * cycles + s.phase_rad + self.band_phase_rad[band]),
            ),
            Complex::from_polar(1.0, -TAU * df_hz * delay),
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct Layout {
    low_hz: [f64; 2],
    df_hz: f64,
    channels: usize,
    reference_hz: f64,
}

impl Layout {
    fn validate(&self) -> Result<(), DynError> {
        if self.channels == 0
            || !self.df_hz.is_finite()
            || self.df_hz <= 0.0
            || !self.reference_hz.is_finite()
            || self.reference_hz <= 0.0
            || self.low_hz.iter().any(|f| !f.is_finite() || *f <= 0.0)
            || self.low_hz[0] + self.channels as f64 * self.df_hz > self.low_hz[1]
        {
            return Err("multiband requires positive, ordered, non-overlapping RF bands".into());
        }
        let offset = (self.low_hz[1] - self.low_hz[0]) / self.df_hz;
        if !offset.is_finite() || offset > (1 << 21) as f64 {
            return Err("multiband RF grid is too large for the delay search".into());
        }
        if (offset - offset.round()).abs() > 1e-6 {
            return Err("multiband RF separation must be an integer multiple of the FFT channel spacing; increase FFT length".into());
        }
        Ok(())
    }
    fn frequency(&self, band: usize, bin: usize) -> f64 {
        self.low_hz[band] + bin as f64 * self.df_hz
    }
}

struct Row {
    time_s: f64,
    duration_s: f64,
    visibility: [Vec<Complex<f64>>; 2],
}

#[derive(Clone, Copy)]
struct Search {
    delay_limit_s: f64,
    rate_limit_hz: f64,
    min_coherence: f64,
}

fn coherent_sums(
    layout: Layout,
    rows: &[Row],
    reference_s: f64,
    delay: f64,
    rate: f64,
    offsets: [f64; 2],
) -> [Complex<f64>; 2] {
    let mut sums = [Complex::new(0.0, 0.0); 2];
    for row in rows {
        let dt = row.time_s - reference_s;
        let slope = delay + rate / layout.reference_hz * dt;
        let step = Complex::from_polar(1.0, -TAU * layout.df_hz * slope);
        for band in 0..2 {
            let cycles = (layout.low_hz[band] - layout.reference_hz) * delay
                + layout.low_hz[band] / layout.reference_hz * rate * dt;
            let mut phase = Complex::from_polar(row.duration_s, -(TAU * cycles + offsets[band]));
            for &v in &row.visibility[band] {
                sums[band] += v * phase;
                phase *= step;
            }
        }
    }
    sums
}

fn score(sums: [Complex<f64>; 2], connected: bool) -> f64 {
    if connected {
        (sums[0] + sums[1]).norm_sqr()
    } else {
        sums[0].norm_sqr() + sums[1].norm_sqr()
    }
}

fn fit(
    layout: Layout,
    rows: &[Row],
    search: Search,
    offsets: [f64; 2],
    connected: bool,
) -> Result<(Solution, [f64; 2]), DynError> {
    layout.validate()?;
    if rows.len() < 2 {
        return Err(
            "at least two short integrations are required to estimate multiband rate".into(),
        );
    }
    let start = rows[0].time_s - 0.5 * rows[0].duration_s;
    let last = rows.last().unwrap();
    let end = last.time_s + 0.5 * last.duration_s;
    let reference_s = 0.5 * (start + end);
    if search.delay_limit_s >= 0.5 / layout.df_hz {
        return Err("delay window exceeds channel-spacing ambiguity; increase FFT length or narrow the delay window".into());
    }
    let max_rf_ratio = layout.frequency(1, layout.channels - 1) / layout.reference_hz;
    if rows
        .iter()
        .any(|r| r.duration_s * search.rate_limit_hz * max_rf_ratio > 0.25)
    {
        return Err("solution-pass integrations are too long for the rate window; reduce --multiband-solve-integration".into());
    }
    let offset = ((layout.low_hz[1] - layout.low_hz[0]) / layout.df_hz).round() as usize;
    let fft_len = (offset + layout.channels)
        .checked_next_power_of_two()
        .and_then(|n| n.checked_mul(2))
        .ok_or("multiband FFT size overflow")?;
    if fft_len > 1 << 22 {
        return Err("multiband search grid exceeds 4M bins; reduce spectral resolution for the solution pass".into());
    }
    let delay_step = 1.0 / (fft_len as f64 * layout.df_hz);
    let desired_rate_step = 1.0 / (4.0 * (end - start) * max_rf_ratio);
    let divisions = ((2.0 * search.rate_limit_hz / desired_rate_step).ceil() as usize)
        .max(2)
        .next_multiple_of(2);
    if divisions > 4096 {
        return Err("multiband rate grid exceeds 4096 steps; narrow the rate window".into());
    }
    let rate_step = 2.0 * search.rate_limit_hz / divisions as f64;
    let fft = FftPlanner::<f64>::new().plan_fft_forward(fft_len);
    // Zero-weight frequency gaps keep their actual RF positions. Workspaces
    // are recycled per Rayon job; no FFT plan/allocation per delay trial.
    let (_, mut delay, mut rate) = (0..=divisions)
        .into_par_iter()
        .map_init(
            || {
                (
                    [
                        vec![Complex::new(0.0, 0.0); fft_len],
                        vec![Complex::new(0.0, 0.0); fft_len],
                    ],
                    vec![Complex::new(0.0, 0.0); fft.get_inplace_scratch_len()],
                )
            },
            |(spectra, scratch), index| {
                let rate = -search.rate_limit_hz + index as f64 * rate_step;
                spectra[0].fill(Complex::new(0.0, 0.0));
                spectra[1].fill(Complex::new(0.0, 0.0));
                for row in rows {
                    let dt = row.time_s - reference_s;
                    let step = Complex::from_polar(
                        1.0,
                        -TAU * layout.df_hz / layout.reference_hz * rate * dt,
                    );
                    for band in 0..2 {
                        let mut phase = Complex::from_polar(
                            row.duration_s,
                            -(TAU * layout.low_hz[band] / layout.reference_hz * rate * dt
                                + offsets[band]),
                        );
                        let begin = if band == 0 { 0 } else { offset };
                        for (k, &v) in row.visibility[band].iter().enumerate() {
                            spectra[band][begin + k] += v * phase;
                            phase *= step;
                        }
                    }
                }
                fft.process_with_scratch(&mut spectra[0], scratch);
                fft.process_with_scratch(&mut spectra[1], scratch);
                let mut best = (f64::NEG_INFINITY, 0.0, rate);
                for k in 0..fft_len {
                    let lag = if k <= fft_len / 2 {
                        k as f64
                    } else {
                        k as f64 - fft_len as f64
                    };
                    let delay = lag * delay_step;
                    if delay.abs() > search.delay_limit_s {
                        continue;
                    }
                    let power = score([spectra[0][k], spectra[1][k]], connected);
                    if power > best.0 {
                        best = (power, delay, rate);
                    }
                }
                best
            },
        )
        .reduce_with(|a, b| if a.0 >= b.0 { a } else { b })
        .unwrap();
    let mut step_d = delay_step.min(search.delay_limit_s);
    let mut step_r = rate_step;
    // Refine on the actual RF coordinates and frequency-scaled rate, rather
    // than fitting a parabola across delay ambiguity sidelobes.
    for _ in 0..18 {
        let candidates: Vec<_> = (-1..=1)
            .flat_map(|i| (-1..=1).map(move |j| (i, j)))
            .collect();
        let best = candidates
            .par_iter()
            .map(|&(i, j)| {
                let d =
                    (delay + i as f64 * step_d).clamp(-search.delay_limit_s, search.delay_limit_s);
                let r =
                    (rate + j as f64 * step_r).clamp(-search.rate_limit_hz, search.rate_limit_hz);
                (
                    score(
                        coherent_sums(layout, rows, reference_s, d, r, offsets),
                        connected,
                    ),
                    d,
                    r,
                )
            })
            .reduce_with(|a, b| if a.0 >= b.0 { a } else { b })
            .unwrap();
        delay = best.1;
        rate = best.2;
        step_d *= 0.5;
        step_r *= 0.5;
    }
    let sums = coherent_sums(layout, rows, reference_s, delay, rate, offsets);
    let mut norm = [0.0; 2];
    for row in rows {
        for band in 0..2 {
            norm[band] +=
                row.visibility[band].iter().map(|v| v.norm()).sum::<f64>() * row.duration_s;
        }
    }
    let denominator = if connected {
        norm.iter().sum::<f64>()
    } else {
        norm[0].hypot(norm[1])
    };
    let coherence = score(sums, connected).sqrt() / denominator.max(f64::MIN_POSITIVE);
    if !coherence.is_finite() || denominator == 0.0 || coherence < search.min_coherence {
        return Err(format!(
            "unreliable multiband solution at {start:.6}..{end:.6}s: coherence={coherence:.6}"
        )
        .into());
    }
    if delay.abs() >= search.delay_limit_s * (1.0 - 1e-6)
        || rate.abs() >= search.rate_limit_hz * (1.0 - 1e-6)
    {
        return Err(
            "multiband solution reached a search boundary; widen the delay/rate windows".into(),
        );
    }
    let phase = if connected {
        (sums[0] + sums[1]).arg()
    } else {
        sums[0].arg()
    };
    let band_phase = if connected {
        offsets
    } else {
        [0.0, (sums[1] * sums[0].conj()).arg()]
    };
    Ok((
        Solution {
            start_s: start,
            end_s: end,
            reference_s,
            delay_s: delay,
            rate_hz: rate,
            phase_rad: phase,
            coherence,
            phase_connected: connected,
        },
        band_phase,
    ))
}

struct CorReader {
    reader: BufReader<File>,
    header: [u8; 256],
    channels: usize,
    sectors: usize,
    low_hz: f64,
    df_hz: f64,
}

fn i32_at(bytes: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}
fn f32_at(bytes: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}
fn f64_at(bytes: &[u8], offset: usize) -> f64 {
    f64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}
fn stamp(bytes: &[u8]) -> (i32, u32) {
    (
        i32_at(bytes, 0),
        u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
    )
}

impl CorReader {
    fn open(path: &Path) -> Result<Self, DynError> {
        let file = File::open(path)?;
        let length = file.metadata()?.len();
        let mut reader = BufReader::new(file);
        let mut header = [0; 256];
        reader.read_exact(&mut header)?;
        let fs = i32_at(&header, 12);
        let fft = i32_at(&header, 24);
        let sectors = i32_at(&header, 28);
        let low = f64_at(&header, 16);
        if header[..4] != [0x83, 0xf9, 0xa2, 0x3e]
            || fs <= 0
            || fft < 2
            || fft % 2 != 0
            || sectors < 1
            || !low.is_finite()
            || low <= 0.0
        {
            return Err(format!("invalid multiband .cor header: {}", path.display()).into());
        }
        let channels = fft as usize / 2;
        if length != 256 + sectors as u64 * (128 + channels as u64 * 8) {
            return Err(format!("incomplete multiband .cor file: {}", path.display()).into());
        }
        Ok(Self {
            reader,
            header,
            channels,
            sectors: sectors as usize,
            low_hz: low,
            df_hz: fs as f64 / fft as f64,
        })
    }
    fn sector(&mut self) -> Result<([u8; 128], Vec<Complex<f64>>), DynError> {
        let mut header = [0; 128];
        self.reader.read_exact(&mut header)?;
        let mut bytes = vec![0; self.channels * 8];
        self.reader.read_exact(&mut bytes)?;
        let spectrum = bytes
            .chunks_exact(8)
            .map(|v| Complex::new(f32_at(v, 0) as f64, f32_at(v, 4) as f64))
            .collect();
        Ok((header, spectrum))
    }
}

fn reader_pair(paths: &[PathBuf; 2]) -> Result<([CorReader; 2], Layout), DynError> {
    let a = CorReader::open(&paths[0])?;
    let b = CorReader::open(&paths[1])?;
    if a.channels != b.channels
        || a.df_hz != b.df_hz
        || a.sectors != b.sectors
        || a.header[32..160] != b.header[32..160]
    {
        return Err(
            "multiband .cor files differ in FFT grid, sector count, baseline geometry or source"
                .into(),
        );
    }
    let bw = a.channels as f64 * a.df_hz;
    let layout = Layout {
        low_hz: [a.low_hz, b.low_hz],
        df_hz: a.df_hz,
        channels: a.channels,
        reference_hz: 0.5 * (a.low_hz + b.low_hz + bw),
    };
    layout.validate()?;
    Ok(([a, b], layout))
}

fn estimate(
    paths: &[PathBuf; 2],
    args: &Args,
    scan_start_s: f64,
    initial_offsets: Option<[f64; 2]>,
) -> Result<Corrections, DynError> {
    let (mut readers, layout) = reader_pair(paths)?;
    if readers[0].sectors < 2 {
        return Err("scan is too short for a multiband delay/rate solution".into());
    }
    let search = Search {
        delay_limit_s: args.multiband_delay_window_ns * 1e-9,
        rate_limit_hz: args.multiband_rate_window_hz,
        min_coherence: args.multiband_min_coherence,
    };
    let mut result = Corrections {
        reference_hz: layout.reference_hz,
        band_phase_rad: initial_offsets.unwrap_or([0.0; 2]),
        solutions: Vec::new(),
    };
    let mut consumed = 0;
    let mut origin = None;
    let mut previous_end = None;
    while consumed < readers[0].sectors {
        let mut rows: Vec<Row> = Vec::new();
        loop {
            let (h1, v1) = readers[0].sector()?;
            let (h2, v2) = readers[1].sector()?;
            if h1[..16] != h2[..16] || h1[112..116] != h2[112..116] {
                return Err(
                    "multiband sectors have different time boundaries or integration durations"
                        .into(),
                );
            }
            let duration = f32_at(&h1, 112) as f64;
            if !duration.is_finite()
                || duration <= 0.0
                || v1
                    .iter()
                    .chain(&v2)
                    .any(|v| !v.re.is_finite() || !v.im.is_finite())
            {
                return Err("invalid multiband visibility/integration".into());
            }
            let (sec, nsec) = stamp(&h1);
            let (end_sec, end_nsec) = stamp(&h1[8..]);
            if nsec >= 1_000_000_000 || end_nsec >= 1_000_000_000 {
                return Err("invalid multiband nanosecond timestamp".into());
            }
            let stamped_duration = (end_sec as i64 - sec as i64) as f64
                + (end_nsec as i64 - nsec as i64) as f64 * 1e-9;
            // Native .cor end stamps are formed from the f32 duration,
            // whereas the next start comes from the f64 FFT sample grid.
            // Allow their rounding difference, including across fit windows.
            let gap = previous_end.map_or(0.0, |(psec, pnsec)| {
                (sec as i64 - psec as i64) as f64 + (nsec as i64 - pnsec as i64) as f64 * 1e-9
            });
            if (stamped_duration - duration).abs() > 1e-6 * duration.max(1.0)
                || gap.abs() > 1e-6 * duration.max(1.0)
            {
                return Err("non-contiguous or inconsistent multiband integration times".into());
            }
            previous_end = Some((end_sec, end_nsec));
            let (origin_sec, origin_nsec) = *origin.get_or_insert((sec, nsec));
            let elapsed = (sec as i64 - origin_sec as i64) as f64
                + (nsec as i64 - origin_nsec as i64) as f64 * 1e-9;
            if let Some(last) = rows.last() {
                let expected = last.time_s + 0.5 * last.duration_s;
                if (elapsed + scan_start_s - expected).abs() > 1e-6 {
                    return Err("non-contiguous multiband integration times".into());
                }
            }
            rows.push(Row {
                time_s: scan_start_s + elapsed + duration * 0.5,
                duration_s: duration,
                visibility: [v1, v2],
            });
            consumed += 1;
            let first = &rows[0];
            let last = rows.last().unwrap();
            let span = last.time_s + last.duration_s * 0.5 - first.time_s + first.duration_s * 0.5;
            if consumed == readers[0].sectors
                || (rows.len() >= 2
                    && span + 1e-6 >= args.multiband_window
                    && readers[0].sectors - consumed != 1)
            {
                break;
            }
        }
        let connected = initial_offsets.is_some()
            || args.multiband_phase_mode == "connected"
            || !result.solutions.is_empty();
        let (solution, offsets) = fit(layout, &rows, search, result.band_phase_rad, connected)?;
        if !connected {
            result.band_phase_rad = offsets;
        }
        println!("[multiband] t={:.6}..{:.6}s delay={:+.6}ns rate={:+.9}Hz@{:.6}MHz coherence={:.6} phase-connected={}",
            solution.start_s, solution.end_s, solution.delay_s * 1e9, solution.rate_hz,
            result.reference_hz / 1e6, solution.coherence, connected);
        result.solutions.push(solution);
    }
    Ok(result)
}

fn write_solutions(
    path: &Path,
    correction: &Corrections,
    time_basis: &str,
) -> Result<(), DynError> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "# format=yi-multiband-solutions-v1\n# software_version={}\n# reference_hz={:.12e}",
        env!("CARGO_PKG_VERSION"),
        correction.reference_hz
    )?;
    writeln!(w, "# time_basis={time_basis}")?;
    writeln!(w, "# observed_phase=phase+band_phase+2pi*((rf-reference)*delay+rf/reference*rate*(time-reference_time))")?;
    writeln!(
        w,
        "# band_phase_rad={:.12e},{:.12e}",
        correction.band_phase_rad[0], correction.band_phase_rad[1]
    )?;
    writeln!(w, "start_s\tend_s\treference_s\tdelay_s\trate_hz\tdelay_rate_sps\tphase_rad\tcoherence\tphase_connected")?;
    for s in &correction.solutions {
        writeln!(
            w,
            "{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.12e}\t{:.9}\t{}",
            s.start_s,
            s.end_s,
            s.reference_s,
            s.delay_s,
            s.rate_hz,
            s.rate_hz / correction.reference_hz,
            s.phase_rad,
            s.coherence,
            s.phase_connected
        )?;
    }
    w.flush()?;
    Ok(())
}

fn pack_joint(
    paths: &[PathBuf; 2],
    output: &Path,
    reference_hz: f64,
    mut inspect: impl FnMut(usize, &[u8; 128], &[Vec<Complex<f64>>; 2]) -> Result<(), DynError>,
) -> Result<(), DynError> {
    let (mut readers, layout) = reader_pair(paths)?;
    let partial = output.with_extension("mbcor.part");
    let mut w = BufWriter::new(File::create(&partial)?);
    w.write_all(b"YIMBCOR\0")?;
    w.write_all(&1u32.to_le_bytes())?;
    w.write_all(&2u32.to_le_bytes())?;
    let occupied = 2.0 * layout.channels as f64 * layout.df_hz;
    w.write_all(&occupied.to_le_bytes())?;
    w.write_all(&reference_hz.to_le_bytes())?;
    for r in &readers {
        w.write_all(&r.header)?;
    }
    for row in 0..readers[0].sectors {
        let (ha, va) = readers[0].sector()?;
        let (hb, vb) = readers[1].sector()?;
        if ha[..16] != hb[..16] || ha[112..116] != hb[112..116] {
            return Err("cannot package multiband output with unequal time grids".into());
        }
        let spectra = [va, vb];
        inspect(row, &ha, &spectra)?;
        for (header, values) in [(&ha, &spectra[0]), (&hb, &spectra[1])] {
            w.write_all(header)?;
            for v in values {
                w.write_all(&(v.re as f32).to_le_bytes())?;
                w.write_all(&(v.im as f32).to_le_bytes())?;
            }
        }
    }
    w.flush()?;
    drop(w);
    std::fs::rename(&partial, output)?;
    println!(
        "[multiband] Joint output: {} occupied bandwidth={:.6}MHz RF span={:.6}MHz",
        output.display(),
        occupied / 1e6,
        (layout.low_hz[1] + layout.channels as f64 * layout.df_hz - layout.low_hz[0]) / 1e6
    );
    Ok(())
}

fn auto_products(directory: &Path, cross: &Path) -> Result<[PathBuf; 2], DynError> {
    let x = CorReader::open(cross)?;
    let stations = [&x.header[32..48], &x.header[80..96]];
    let mut autos = [None, None];
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.extension().and_then(|s| s.to_str()) != Some("cor") {
            continue;
        }
        let r = CorReader::open(&path)?;
        if r.header[32..48] != r.header[80..96] {
            continue;
        }
        for (index, station) in stations.iter().enumerate() {
            if &r.header[32..48] == *station && autos[index].replace(path.clone()).is_some() {
                return Err(format!("duplicate ACF in {}", directory.display()).into());
            }
        }
    }
    Ok([
        autos[0].take().ok_or("missing station 1 multiband ACF")?,
        autos[1].take().ok_or("missing station 2 multiband ACF")?,
    ])
}

fn only_cross_product(directory: &Path) -> Result<PathBuf, DynError> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(directory)? {
        let path = entry?.path();
        if path.extension().and_then(|s| s.to_str()) == Some("cor") {
            let r = CorReader::open(&path)?;
            if r.header[32..48] != r.header[80..96] {
                files.push(path);
            }
        }
    }
    if files.len() != 1 {
        return Err(format!("expected exactly one XCF in {}", directory.display()).into());
    }
    Ok(files.remove(0))
}

// Calibration references use Unix seconds; target corrections use seconds
// since the target process epoch, like FrameDelayEntry. Transfer the nearest
// local linear model, including its phase evolution, never fit target noise.
fn transfer_calibration(
    calibration: &Corrections,
    epoch_s: f64,
    start_s: f64,
    end_s: f64,
    max_gap: f64,
) -> Result<Corrections, DynError> {
    if calibration.solutions.is_empty() || end_s <= start_s {
        return Err("no multiband calibration solutions or empty target scan".into());
    }
    let mut solutions = calibration.solutions.clone();
    solutions.sort_by(|a, b| a.reference_s.total_cmp(&b.reference_s));
    let mut target = Vec::new();
    for i in 0..solutions.len() {
        let mut s = solutions[i].clone();
        let left = if i == 0 {
            f64::NEG_INFINITY
        } else {
            0.5 * (solutions[i - 1].reference_s + s.reference_s) - epoch_s
        };
        let right = if i + 1 == solutions.len() {
            f64::INFINITY
        } else {
            0.5 * (solutions[i + 1].reference_s + s.reference_s) - epoch_s
        };
        s.start_s = start_s.max(left);
        s.end_s = end_s.min(right);
        s.reference_s -= epoch_s;
        if s.end_s <= s.start_s {
            continue;
        }
        let gap = (s.start_s - s.reference_s)
            .abs()
            .max((s.end_s - s.reference_s).abs());
        if gap > max_gap {
            return Err(format!("target is {gap:.3}s from the nearest multiband calibration solution; maximum is {max_gap:.3}s").into());
        }
        target.push(s);
    }
    Ok(Corrections {
        reference_hz: calibration.reference_hz,
        band_phase_rad: calibration.band_phase_rad,
        solutions: target,
    })
}

pub fn run(
    args: &Args,
    cpu_threads: usize,
    reader_core: Option<core_affinity::CoreId>,
) -> Result<(), DynError> {
    for (name, value) in [
        ("window", args.multiband_window),
        ("solve-integration", args.multiband_solve_integration),
        ("delay-window-ns", args.multiband_delay_window_ns),
        ("rate-window-hz", args.multiband_rate_window_hz),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(format!("--multiband-{name} must be finite and positive").into());
        }
    }
    if args
        .multiband_calibration_max_gap
        .is_some_and(|v| !v.is_finite() || v <= 0.0)
    {
        return Err("--multiband-calibration-max-gap must be finite and positive".into());
    }
    if !args.multiband_min_coherence.is_finite()
        || !(0.0..=1.0).contains(&args.multiband_min_coherence)
    {
        return Err("--multiband-min-coherence must be in 0..1".into());
    }
    if args.multiband_raw_directory.is_none()
        && (args.multiband_ant1.is_none() || args.multiband_ant2.is_none())
    {
        return Err("specify --multiband-raw-directory or both --multiband-ant1/--multiband-ant2 (separate RF-band files)".into());
    }
    if args.gain_phasecal
        || args.gain_uncalibrated_schedule.is_some()
        || args.phased_validation
        || args.model_sweep
        || args.fringe.is_some()
        || args.band.is_some()
    {
        return Err("multiband cannot be combined with gain/phased workflows, model sweep, fringe QL or --band".into());
    }
    let schedules = [
        args.schedule.as_ref().unwrap(),
        args.multiband_schedule.as_ref().unwrap(),
    ];
    let mut base_args = [args.clone(), args.clone()];
    base_args[1].schedule = Some(schedules[1].clone());
    base_args[1].schedule_xml = args.multiband_schedule_xml.clone();
    let metadata = [
        crate::parse_schedule_args(&base_args[0], None)?,
        crate::parse_schedule_args(&base_args[1], None)?,
    ];
    if metadata[0].processes.len() != metadata[1].processes.len() {
        return Err("multiband XMLs must contain the same simultaneous scans".into());
    }
    let indices: Vec<_> = args
        .process_index
        .map(|i| vec![i])
        .unwrap_or_else(|| (0..metadata[0].processes.len()).collect());
    if indices.is_empty() {
        return Err("multiband XMLs contain no scans".into());
    }
    let mut work_indices = Vec::new();
    if !args.multiband_calibrator.is_empty() {
        if args.epoch.is_some() {
            return Err(
                "multiband calibration transfer requires the original XML epochs; omit --epoch"
                    .into(),
            );
        }
        for name in &args.multiband_calibrator {
            if !metadata[0]
                .processes
                .iter()
                .any(|p| p.object.as_ref() == Some(name))
            {
                return Err(format!("multiband calibrator {name:?} is absent from the XML").into());
            }
        }
        work_indices.extend(
            metadata[0]
                .processes
                .iter()
                .enumerate()
                .filter(|(_, p)| {
                    p.object
                        .as_ref()
                        .is_some_and(|n| args.multiband_calibrator.contains(n))
                })
                .map(|(i, _)| i),
        );
    }
    for &index in &indices {
        if !work_indices.contains(&index) {
            work_indices.push(index);
        }
    }
    let mut calibration: Option<Corrections> = None;
    let mut calibration_layout: Option<Layout> = None;
    let mut baseline = None;
    let root = args.cor_directory.as_ref().unwrap().join("multiband");
    for index in work_indices {
        let a = metadata[0]
            .processes
            .get(index)
            .ok_or("multiband process index out of range")?;
        let b = &metadata[1].processes[index];
        if !crate::same_validation_scan(a, b) {
            return Err("multiband process epoch, skip, length and source must match".into());
        }
        let meta = [
            crate::parse_schedule_args(&base_args[0], Some(index))?,
            crate::parse_schedule_args(&base_args[1], Some(index))?,
        ];
        if meta[0].sampling_hz != meta[1].sampling_hz
            || meta[0].fft != meta[1].fft
            || meta[0].ant1_station_name != meta[1].ant1_station_name
            || meta[0].ant2_station_name != meta[1].ant2_station_name
            || meta[0].ant1_ecef_m != meta[1].ant1_ecef_m
            || meta[0].ant2_ecef_m != meta[1].ant2_ecef_m
            || meta[0].ra != meta[1].ra
            || meta[0].dec != meta[1].dec
            || meta[0].output_sec != meta[1].output_sec
            || meta
                .iter()
                .any(|m| m.inband.unwrap_or(1) != 1 || m.pulsar.is_some())
        {
            return Err("multiband requires identical baseline/source, sampling, FFT and output integration; inband=1 and no pulsar folding".into());
        }
        let current_baseline = (
            meta[0].ant1_station_name.clone(),
            meta[0].ant2_station_name.clone(),
            meta[0].ant1_ecef_m,
            meta[0].ant2_ecef_m,
        );
        if baseline.as_ref().is_some_and(|b| b != &current_baseline) {
            return Err(
                "all multiband scans must use the same ordered baseline and station coordinates"
                    .into(),
            );
        }
        baseline = Some(current_baseline);
        let directory = root.join(format!("scan{index:04}"));
        std::fs::create_dir_all(&directory)?;
        let mut band_args = base_args.clone();
        for band in 0..2 {
            let v = &mut band_args[band];
            v.schedule = Some(schedules[band].clone());
            v.process_index = Some(index);
            v.multiband_schedule = None;
            v.fringe = None;
            v.model_diagnostics = false;
            v.compact_logs = true;
            v.cor_label_override = Some("multiband".into());
            if band == 1 {
                v.raw_directory = args
                    .multiband_raw_directory
                    .clone()
                    .or_else(|| args.raw_directory.clone());
                v.ant1 = args.multiband_ant1.clone();
                v.ant2 = args.multiband_ant2.clone();
            }
        }
        let epoch = args.epoch.as_deref().unwrap_or(&a.epoch);
        let files = [
            crate::resolve_input_paths(&band_args[0], epoch, Some(&meta[0]))?,
            crate::resolve_input_paths(&band_args[1], epoch, Some(&meta[1]))?,
        ];
        if std::fs::canonicalize(&files[0].0)? == std::fs::canonicalize(&files[1].0)?
            || std::fs::canonicalize(&files[0].1)? == std::fs::canonicalize(&files[1].1)?
        {
            return Err("the two RF bands must not use the same raw file".into());
        }
        let is_calibrator = a
            .object
            .as_ref()
            .is_some_and(|n| args.multiband_calibrator.contains(n));
        let transfer = !args.multiband_calibrator.is_empty() && !is_calibrator;
        let mut solution_paths = None;
        let correction = if transfer {
            let layout = calibration_layout.ok_or("multiband calibration layout is missing")?;
            if meta.iter().enumerate().any(|(band, m)| {
                m.obsfreq_mhz.map(|f| f * 1e6) != Some(layout.low_hz[band])
                    || m.sampling_hz.map(|fs| fs as f64)
                        != Some(layout.df_hz * layout.channels as f64 * 2.0)
            }) {
                return Err("calibrator and target RF bands/sampling must match".into());
            }
            let (epoch, _) = crate::epoch_to_yyyydddhhmmss(&a.epoch)?;
            let start = a.skip_sec + args.skip;
            let end = args
                .length
                .or(a.length_sec)
                .ok_or("calibrated target scan requires XML length or --length")?;
            let c = transfer_calibration(
                calibration.as_ref().ok_or("no multiband calibration")?,
                epoch as f64,
                start,
                end,
                args.multiband_calibration_max_gap.unwrap_or(f64::INFINITY),
            )?;
            println!("[multiband] scan {} source {:?}: transfer calibrator solutions (no target fringe fit)", index + 1, a.object);
            c
        } else {
            let mut solve_paths = [PathBuf::new(), PathBuf::new()];
            for band in 0..2 {
                let mut v = band_args[band].clone();
                let out = directory.join(format!("solve-band{}", band + 1));
                v.cor_directory = Some(out.clone());
                v.xcf_only = true;
                v.integration_rate_override = Some(1.0 / args.multiband_solve_integration);
                println!(
                    "[multiband] scan {} solution pass band {}",
                    index + 1,
                    band + 1
                );
                crate::run_once(v, crate::RunMode::Corr, cpu_threads, reader_core)?;
                solve_paths[band] = only_cross_product(&out)?;
            }
            let c = estimate(
                &solve_paths,
                args,
                a.skip_sec + args.skip,
                calibration.as_ref().map(|c| c.band_phase_rad),
            )?;
            solution_paths = Some(solve_paths.clone());
            if is_calibrator {
                let (_, layout) = reader_pair(&solve_paths)?;
                if let Some(previous) = calibration_layout {
                    if previous.low_hz != layout.low_hz
                        || previous.df_hz != layout.df_hz
                        || previous.channels != layout.channels
                    {
                        return Err(
                            "all multiband calibrators must use the same RF/FFT grid".into()
                        );
                    }
                }
                calibration_layout = Some(layout);
                let (epoch, _) = crate::epoch_to_yyyydddhhmmss(&a.epoch)?;
                let pool = calibration.get_or_insert_with(|| Corrections {
                    reference_hz: c.reference_hz,
                    band_phase_rad: c.band_phase_rad,
                    solutions: Vec::new(),
                });
                for s in &c.solutions {
                    let mut s = s.clone();
                    s.start_s += epoch as f64;
                    s.end_s += epoch as f64;
                    s.reference_s += epoch as f64;
                    pool.solutions.push(s);
                }
            }
            c
        };
        let correction = Arc::new(correction);
        write_solutions(
            &directory.join("solutions.tsv"),
            &correction,
            "seconds_since_XML_process_epoch",
        )?;
        diagnostics::plot_solutions(&directory, &correction)?;
        if !indices.contains(&index) {
            continue;
        }
        let mut final_paths = [PathBuf::new(), PathBuf::new()];
        let mut auto_paths = [
            [PathBuf::new(), PathBuf::new()],
            [PathBuf::new(), PathBuf::new()],
        ];
        for band in 0..2 {
            let mut v = band_args[band].clone();
            let out = directory.join(format!("band{}", band + 1));
            v.cor_directory = Some(out.clone());
            v.multiband_correction = Some(Arc::clone(&correction));
            v.multiband_band_index = band;
            // Final multiband products always include both station powers.
            // Only the short solution pass uses the XCF-only fast path.
            v.xcf_only = false;
            println!(
                "[multiband] scan {} corrected correlation band {}",
                index + 1,
                band + 1
            );
            crate::run_once(v, crate::RunMode::Corr, cpu_threads, reader_core)?;
            final_paths[band] = only_cross_product(&out)?;
            let autos = auto_products(&out, &final_paths[band])?;
            for station in 0..2 {
                auto_paths[station][band] = autos[station].clone();
            }
        }
        diagnostics::write(
            &directory,
            &final_paths,
            &auto_paths,
            solution_paths.as_ref(),
            &correction,
            args,
            epoch,
            transfer,
        )?;
    }
    if let Some(mut c) = calibration {
        c.solutions
            .sort_by(|a, b| a.reference_s.total_cmp(&b.reference_s));
        write_solutions(&root.join("calibration.tsv"), &c, "unix_seconds")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cx_schedule_keeps_shared_scans_and_each_bands_clock_and_raw_directory() {
        use clap::Parser;
        let path = std::env::temp_dir().join(format!("cx-schedule-{}.xml", std::process::id()));
        let xml = include_str!("../examples/I26280X_all_KL_CX.xml");
        std::fs::write(&path, xml).unwrap();
        let mut args = Args::parse_from([
            "yi-corr",
            "--schedule",
            path.to_str().unwrap(),
            "--raw",
            "/observations",
        ]);
        configure_from_xml(&mut args).unwrap();
        assert_eq!(
            args.raw_directory.as_deref(),
            Some(Path::new("/observations/c"))
        );
        assert_eq!(
            args.multiband_raw_directory.as_deref(),
            Some(Path::new("/observations/x"))
        );
        assert_eq!(args.multiband_calibrator, vec!["NRAO530"]);
        let low = crate::parse_schedule_args(&args, Some(0)).unwrap();
        args.schedule_xml = args.multiband_schedule_xml.clone();
        let high = crate::parse_schedule_args(&args, Some(0)).unwrap();
        assert_eq!(low.processes.len(), 10);
        assert_eq!(high.processes.len(), 10);
        assert_eq!(low.obsfreq_mhz, Some(6600.0));
        assert_eq!(high.obsfreq_mhz, Some(8192.0));
        assert_eq!(low.output_sec, Some(1.0));
        assert_eq!(low.ant2_clock_delay_s, Some(1.707118524609375e-6));
        assert_eq!(high.ant2_clock_delay_s, Some(1.7176596209375e-6));
        for (a, b) in low.processes.iter().zip(&high.processes) {
            assert!(crate::same_validation_scan(a, b));
        }
        for bad in [
            xml.replace("name=\"x\"", "name=\"c\""),
            xml.replace("name=\"x\"", "name=\"../x\""),
            xml.replace("<band name=\"x\">", "<other>")
                .replace("</band>\n  </multiband>", "</other>\n  </multiband>"),
        ] {
            std::fs::write(&path, bad).unwrap();
            let mut args = Args::parse_from([
                "yi-corr",
                "--schedule",
                path.to_str().unwrap(),
                "--raw",
                "/observations",
            ]);
            assert!(configure_from_xml(&mut args).is_err());
        }
        std::fs::remove_file(path).unwrap();
    }

    fn fixture(offset: f64, delay: f64, rate: f64) -> (Layout, Vec<Row>) {
        let layout = Layout {
            low_hz: [6600e6, 8192e6],
            df_hz: 8e6,
            channels: 64,
            reference_hz: 7652e6,
        };
        let rows = (0..16)
            .map(|i| {
                let time = (i as f64 + 0.5) * 0.25;
                let visibility = std::array::from_fn(|band| {
                    (0..layout.channels)
                        .map(|k| {
                            let f = layout.frequency(band, k);
                            let phase = 0.3
                                + if band == 1 { offset } else { 0.0 }
                                + TAU
                                    * ((f - layout.reference_hz) * delay
                                        + f / layout.reference_hz * rate * (time - 2.0));
                            Complex::from_polar(1.0 + k as f64 / 256.0, phase)
                        })
                        .collect()
                });
                Row {
                    time_s: time,
                    duration_s: 0.25,
                    visibility,
                }
            })
            .collect();
        (layout, rows)
    }

    #[test]
    fn actual_rf_gap_and_frequency_scaled_rate_recover_injected_solution() {
        let (layout, rows) = fixture(0.0, 2.371e-9, -0.0713);
        let search = Search {
            delay_limit_s: 10e-9,
            rate_limit_hz: 0.2,
            min_coherence: 0.9,
        };
        let (s, offsets) = fit(layout, &rows, search, [0.0; 2], true).unwrap();
        assert!((s.delay_s - 2.371e-9).abs() < 1e-12, "{s:?}");
        assert!((s.rate_hz + 0.0713).abs() < 2e-5, "{s:?}");
        assert!(s.coherence > 0.99999);
        let correction = Corrections {
            reference_hz: layout.reference_hz,
            band_phase_rad: offsets,
            solutions: vec![s],
        };
        for row in &rows {
            for band in 0..2 {
                let (mut phase, step) = correction.phase_start_and_step(
                    band,
                    row.time_s,
                    layout.low_hz[band],
                    layout.df_hz,
                );
                for v in &row.visibility[band] {
                    assert!((v * phase).arg().abs() < 0.002);
                    phase *= step;
                }
            }
        }
    }

    #[test]
    fn bootstrap_fits_unknown_if_phase_then_connects_later_window() {
        let (layout, rows) = fixture(1.17, -1.913e-9, 0.0632);
        let search = Search {
            delay_limit_s: 10e-9,
            rate_limit_hz: 0.2,
            min_coherence: 0.9,
        };
        let (s, offsets) = fit(layout, &rows, search, [0.0; 2], false).unwrap();
        assert!((s.delay_s + 1.913e-9).abs() < 1e-12, "{s:?}");
        assert!((offsets[1] - 1.17).abs() < 0.01, "{offsets:?}");
        assert!(!s.phase_connected);
        let (_, later) = fixture(1.17, 3.157e-9, -0.0824);
        let (later_s, _) = fit(layout, &later, search, offsets, true).unwrap();
        assert!((later_s.delay_s - 3.157e-9).abs() < 2e-12, "{later_s:?}");
        assert!((later_s.rate_hz + 0.0824).abs() < 2e-5);
    }

    #[test]
    fn invalid_or_unidentifiable_searches_are_rejected() {
        let (layout, mut rows) = fixture(0.0, 0.0, 0.0);
        let search = Search {
            delay_limit_s: 10e-9,
            rate_limit_hz: 0.2,
            min_coherence: 0.5,
        };
        assert!(fit(layout, &rows[..1], search, [0.0; 2], true).is_err());
        let wrong = Layout {
            low_hz: [6600e6, 8193e6],
            ..layout
        };
        assert!(wrong.validate().is_err());
        for row in &mut rows {
            for values in &mut row.visibility {
                values.fill(Complex::new(0.0, 0.0));
            }
        }
        assert!(fit(layout, &rows, search, [0.0; 2], true).is_err());
    }

    #[test]
    fn calibration_transfer_preserves_phase_and_rate_across_epochs() {
        let (layout, rows) = fixture(1.17, 2.371e-9, -0.0713);
        let (s, offsets) = fit(
            layout,
            &rows,
            Search {
                delay_limit_s: 10e-9,
                rate_limit_hz: 0.2,
                min_coherence: 0.9,
            },
            [0.0; 2],
            false,
        )
        .unwrap();
        let local = Corrections {
            reference_hz: layout.reference_hz,
            band_phase_rad: offsets,
            solutions: vec![s.clone()],
        };
        let epoch = 946684800.0;
        let mut absolute = s;
        absolute.start_s += epoch;
        absolute.end_s += epoch;
        absolute.reference_s += epoch;
        let calibration = Corrections {
            solutions: vec![absolute],
            ..local.clone()
        };
        let target = transfer_calibration(&calibration, epoch + 10.0, 0.0, 4.0, 30.0).unwrap();
        for band in 0..2 {
            let a = local.phase_start_and_step(band, 12.5, layout.low_hz[band], layout.df_hz);
            let b = target.phase_start_and_step(band, 2.5, layout.low_hz[band], layout.df_hz);
            assert!((a.0 - b.0).norm() < 1e-12);
            assert!((a.1 - b.1).norm() < 1e-12);
        }
        assert!(transfer_calibration(&calibration, epoch + 10.0, 0.0, 4.0, 5.0).is_err());
    }

    #[test]
    fn native_f32_integration_timestamps_allow_rounding_across_windows() {
        use crate::cor::{CorHeaderConfig, CorStation, CorWriter};
        use clap::Parser;
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory =
            std::env::temp_dir().join(format!("multiband-rounding-{}-{nonce}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let paths = [directory.join("low.cor"), directory.join("high.cor")];
        let station = CorStation {
            name: "ANT1",
            code: b'A',
            ecef_m: [0.0; 3],
        };
        let (layout, _) = fixture(0.0, 0.0, 0.0);
        for band in 0..2 {
            let cfg = CorHeaderConfig {
                sampling_speed_hz: 1024000000,
                observing_frequency_hz: layout.low_hz[band],
                fft_point: 128,
                number_of_sector_hint: 16,
                clock_reference_unix_sec: 946684800,
                source_name: "CAL".into(),
                source_ra_rad: 0.0,
                source_dec_rad: 0.0,
            };
            let mut writer = CorWriter::create(
                &paths[band],
                &cfg,
                station,
                CorStation {
                    name: "ANT2",
                    code: b'B',
                    ecef_m: [110.0, 0.0, 0.0],
                },
            )
            .unwrap();
            for i in 0..16 {
                let spectrum = (0..64)
                    .map(|k| {
                        Complex::<f32>::from_polar(
                            1.0,
                            (0.3 + TAU
                                * ((layout.frequency(band, k) - layout.reference_hz) * 2.371e-9
                                    + layout.frequency(band, k) / layout.reference_hz
                                        * -0.0713
                                        * ((i as f64 + 0.5) * 0.05 - 0.4)))
                                as f32,
                        )
                    })
                    .collect::<Vec<_>>();
                writer
                    .write_sector_at(946684800, i as f64 * 0.05, 0.05, &spectrum)
                    .unwrap();
            }
            writer.finalize().unwrap();
        }
        let args = Args::parse_from([
            "yi-corr",
            "--multiband-window",
            "0.2",
            "--multiband-delay-window-ns",
            "10",
        ]);
        let c = estimate(&paths, &args, 0.0, Some([0.0; 2])).unwrap();
        assert!(c.solutions.len() >= 3);
        for s in &c.solutions {
            let expected = 2.371e-9 - 0.0713 / layout.reference_hz * (s.reference_s - 0.4);
            assert!((s.delay_s - expected).abs() < 2e-12, "{s:?}");
            assert!((s.rate_hz + 0.0713).abs() < 1e-4, "{s:?}");
        }
        std::fs::remove_dir_all(directory).unwrap();
    }
}
