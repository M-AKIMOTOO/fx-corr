//! Streaming QA from native COR products; no additional RAW pass or target fit.
use super::*;
use plotters::coord::Shift;
use plotters::prelude::*;
use std::path::Path;

const MAX_PLOT_POINTS: usize = 2048;
const COLORS: [RGBColor; 4] = [BLUE, RED, GREEN, RGBColor(140, 80, 180)];
type Points = Vec<(f64, f64)>;

fn plot_stride(length: usize) -> usize {
    length
        .saturating_sub(1)
        .div_ceil(MAX_PLOT_POINTS - 1)
        .max(1)
}

fn padded_range(values: impl Iterator<Item = f64>) -> std::ops::Range<f64> {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for x in values.filter(|x| x.is_finite()) {
        lo = lo.min(x);
        hi = hi.max(x);
    }
    if !lo.is_finite() {
        return -1.0..1.0;
    }
    let pad = if hi > lo {
        (hi - lo) * 0.06
    } else {
        lo.abs().max(1.0) * 0.06
    };
    (lo - pad)..(hi + pad)
}

fn panel(
    area: &DrawingArea<BitMapBackend<'_>, Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[(&str, Points)],
    phase: bool,
) -> Result<(), DynError> {
    let xr = padded_range(series.iter().flat_map(|(_, p)| p.iter().map(|v| v.0)));
    let yr = if phase {
        -180.0..180.0
    } else {
        padded_range(series.iter().flat_map(|(_, p)| p.iter().map(|v| v.1)))
    };
    let mut chart = ChartBuilder::on(area)
        .caption(title, ("sans-serif", 20))
        .margin(12)
        .x_label_area_size(42)
        .y_label_area_size(88)
        .build_cartesian_2d(xr, yr)?;
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .label_style(("sans-serif", 14))
        .axis_desc_style(("sans-serif", 16))
        .draw()?;
    // Scatter preserves phase wrapping and never joins the unobserved RF gap.
    for (index, (label, points)) in series.iter().enumerate() {
        let color = COLORS[index % COLORS.len()];
        chart
            .draw_series(
                points
                    .iter()
                    .filter(|(x, y)| x.is_finite() && y.is_finite())
                    .map(|p| Circle::new(*p, 2, color.filled())),
            )?
            .label(*label)
            .legend(move |(x, y)| Circle::new((x + 8, y), 4, color.filled()));
    }
    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.9))
        .border_style(BLACK)
        .draw()?;
    Ok(())
}

fn figure(
    path: &Path,
    panels: &[(&str, &str, &str, Vec<(&str, Points)>, bool)],
) -> Result<(), DynError> {
    let root = BitMapBackend::new(path, (1280, 380 * panels.len() as u32)).into_drawing_area();
    root.fill(&WHITE)?;
    for (area, (title, x, y, data, phase)) in
        root.split_evenly((panels.len(), 1)).iter().zip(panels)
    {
        panel(area, title, x, y, data, *phase)?;
    }
    root.present()?;
    Ok(())
}

pub(super) fn plot_solutions(directory: &Path, correction: &Corrections) -> Result<(), DynError> {
    let points = |value: fn(&Solution) -> f64| -> Vec<(&str, Points)> {
        vec![(
            "model reference",
            correction
                .solutions
                .iter()
                .step_by(plot_stride(correction.solutions.len()))
                .map(|s| (s.reference_s, value(s)))
                .collect(),
        )]
    };
    figure(
        &directory.join("solutions.png"),
        &[
            (
                "Common residual delay",
                "Model reference time from XML epoch [s]",
                "Delay [ns]",
                points(|s| s.delay_s * 1e9),
                false,
            ),
            (
                "Common residual rate at reference RF",
                "Model reference time from XML epoch [s]",
                "Rate [Hz]",
                points(|s| s.rate_hz),
                false,
            ),
            (
                "Common phase at reference time / RF",
                "Model reference time from XML epoch [s]",
                "Phase [deg]",
                points(|s| s.phase_rad.to_degrees()),
                true,
            ),
            (
                "Calibrator fit coherence (not detection S/N)",
                "Model reference time from XML epoch [s]",
                "Coherence",
                points(|s| s.coherence),
                false,
            ),
        ],
    )
}

struct Summary {
    layout: Layout,
    spectra: [Vec<Complex<f64>>; 2],
    exposure: [Vec<f64>; 2],
    means: [Complex<f64>; 2],
    weights: [f64; 2],
    duration: f64,
    rows: usize,
    epoch: f64,
    amplitudes: [Points; 3],
    phases: [Points; 3],
    relative_phase: Points,
    applied: [Points; 4],
}

impl Summary {
    fn new(layout: Layout) -> Self {
        Self {
            layout,
            spectra: std::array::from_fn(|_| vec![Complex::new(0.0, 0.0); layout.channels]),
            exposure: std::array::from_fn(|_| vec![0.0; layout.channels]),
            means: [Complex::new(0.0, 0.0); 2],
            weights: [0.0; 2],
            duration: 0.0,
            rows: 0,
            epoch: 0.0,
            amplitudes: Default::default(),
            phases: Default::default(),
            relative_phase: Vec::new(),
            applied: Default::default(),
        }
    }

    fn add(
        &mut self,
        header: &[u8; 128],
        spectra: &[Vec<Complex<f64>>; 2],
        keep: bool,
        text: &mut impl Write,
    ) -> Result<(), DynError> {
        let (sec, nsec) = stamp(header);
        let duration = f32_at(header, 112) as f64;
        if nsec >= 1_000_000_000 || !duration.is_finite() || duration <= 0.0 {
            return Err("invalid multiband diagnostic timestamp/integration".into());
        }
        let unix = sec as f64 + nsec as f64 * 1e-9;
        if self.rows == 0 {
            self.epoch = unix;
        }
        let elapsed = unix - self.epoch + duration * 0.5;
        let mut means = [Complex::new(0.0, 0.0); 3];
        let mut counts = [0; 2];
        for band in 0..2 {
            for (k, &v) in spectra[band].iter().enumerate() {
                if !v.re.is_finite() || !v.im.is_finite() {
                    return Err("nonfinite multiband QA spectrum".into());
                }
                // Zero bins have no valid common-band XCF weight.
                if v == Complex::new(0.0, 0.0) {
                    continue;
                }
                self.spectra[band][k] += v * duration;
                self.exposure[band][k] += duration;
                means[band] += v;
                counts[band] += 1;
            }
            self.means[band] += means[band] * duration;
            self.weights[band] += counts[band] as f64 * duration;
        }
        let count = counts.iter().sum::<usize>();
        if count > 0 {
            means[2] = (means[0] + means[1]) / count as f64;
        }
        for band in 0..2 {
            if counts[band] > 0 {
                means[band] /= counts[band] as f64;
            }
        }
        let relative =
            if counts.iter().all(|n| *n > 0) && means[0].norm() > 0.0 && means[1].norm() > 0.0 {
                (means[1] * means[0].conj()).arg().to_degrees()
            } else {
                f64::NAN
            };
        write!(text, "{unix:.9}\t{duration:.9}\t{elapsed:.9}")?;
        for (i, v) in means.iter().enumerate() {
            let valid = if i == 2 { count > 0 } else { counts[i] > 0 };
            let angle = if valid && v.norm() > 0.0 {
                v.arg().to_degrees()
            } else {
                f64::NAN
            };
            write!(
                text,
                "\t{:.12e}\t{:.12e}\t{:.12e}\t{angle:.9}",
                v.re,
                v.im,
                v.norm()
            )?;
            if keep && valid {
                self.amplitudes[i].push((elapsed, v.norm()));
                self.phases[i].push((elapsed, angle));
            }
        }
        writeln!(
            text,
            "\t{}\t{}\t{:.12e}\t{relative:.9}",
            counts[0],
            counts[1],
            count as f64 * self.layout.df_hz
        )?;
        if keep {
            self.relative_phase.push((elapsed, relative));
        }
        self.duration += duration;
        self.rows += 1;
        Ok(())
    }

    fn spectrum_points(&self, phase: bool) -> Vec<(&str, Points)> {
        let stride = self.layout.channels.div_ceil(MAX_PLOT_POINTS);
        (0..2)
            .map(|band| {
                let data = (0..self.layout.channels)
                    .step_by(stride)
                    .filter_map(|k| {
                        let weight = self.exposure[band][k];
                        if weight == 0.0 {
                            return None;
                        }
                        let v = self.spectra[band][k] / weight;
                        Some((
                            self.layout.frequency(band, k) / 1e6,
                            if phase { v.arg().to_degrees() } else { v.re },
                        ))
                    })
                    .collect();
                (if band == 0 { "band 1" } else { "band 2" }, data)
            })
            .collect()
    }

    fn spectral_amplitude(&self) -> Vec<(&str, Points)> {
        let mut points = self.spectrum_points(false);
        let stride = self.layout.channels.div_ceil(MAX_PLOT_POINTS);
        for (band, (_, data)) in points.iter_mut().enumerate() {
            *data = (0..self.layout.channels)
                .step_by(stride)
                .filter_map(|k| {
                    let weight = self.exposure[band][k];
                    (weight > 0.0).then(|| {
                        (
                            self.layout.frequency(band, k) / 1e6,
                            (self.spectra[band][k] / weight).norm(),
                        )
                    })
                })
                .collect();
        }
        points
    }

    fn write_spectrum(&self, path: &Path) -> Result<(), DynError> {
        let mut text = BufWriter::new(File::create(path)?);
        writeln!(
            text,
            "band\tchannel\tfrequency_hz\texposure_s\treal\timag\tamplitude\tphase_deg"
        )?;
        for band in 0..2 {
            for k in 0..self.layout.channels {
                let weight = self.exposure[band][k];
                let v = if weight > 0.0 {
                    self.spectra[band][k] / weight
                } else {
                    Complex::new(0.0, 0.0)
                };
                let angle = if weight > 0.0 && v.norm() > 0.0 {
                    v.arg().to_degrees()
                } else {
                    f64::NAN
                };
                writeln!(
                    text,
                    "{}\t{k}\t{:.12e}\t{weight:.9}\t{:.12e}\t{:.12e}\t{:.12e}\t{angle:.9}",
                    band + 1,
                    self.layout.frequency(band, k),
                    v.re,
                    v.im,
                    v.norm()
                )?;
            }
        }
        text.flush()?;
        Ok(())
    }
}

fn time_table(path: &Path) -> Result<BufWriter<File>, DynError> {
    let mut text = BufWriter::new(File::create(path)?);
    writeln!(text, "unix_s\tintegration_s\tmidpoint_from_scan_start_s\tband1_real\tband1_imag\tband1_amplitude\tband1_phase_deg\tband2_real\tband2_imag\tband2_amplitude\tband2_phase_deg\tjoint_real\tjoint_imag\tjoint_amplitude\tjoint_phase_deg\tband1_valid_channels\tband2_valid_channels\tusable_bandwidth_hz\tband2_minus_band1_phase_deg")?;
    Ok(text)
}

fn collect(
    paths: &[PathBuf; 2],
    directory: &Path,
    name: &str,
    reference_hz: Option<f64>,
    applied: Option<(&Corrections, i64)>,
) -> Result<Summary, DynError> {
    let (mut readers, layout) = reader_pair(paths)?;
    let rows = readers[0].sectors;
    let stride = plot_stride(rows);
    let mut result = Summary::new(layout);
    let mut text = time_table(&directory.join(format!("{name}-time.tsv")))?;
    let mut applied_text = applied.map(|_| {
        let mut text = BufWriter::new(File::create(directory.join("applied.tsv"))?);
        writeln!(text, "process_time_s\tmodel_reference_s\tdelay_s\trate_hz\tcommon_phase_rad\tband1_if_phase_rad\tband2_if_phase_rad\tcalibrator_coherence")?;
        Ok::<_, DynError>(text)
    }).transpose()?;
    let mut inspect = |row: usize, h: &[u8; 128], v: &[Vec<Complex<f64>>; 2]| {
        let keep = row % stride == 0 || row + 1 == rows;
        result.add(h, v, keep, &mut text)?;
        if let (Some((c, epoch)), Some(table)) = (applied, applied_text.as_mut()) {
            let (sec, nsec) = stamp(h);
            let t = (sec as i64 - epoch) as f64 + nsec as f64 * 1e-9 + f32_at(h, 112) as f64 * 0.5;
            let index = c
                .solutions
                .partition_point(|s| s.end_s <= t)
                .min(c.solutions.len() - 1);
            let s = &c.solutions[index];
            let d = s.delay_s + s.rate_hz / c.reference_hz * (t - s.reference_s);
            let p = s.phase_rad + TAU * s.rate_hz * (t - s.reference_s);
            let if_phase = c.band_phase_rad[1] - c.band_phase_rad[0];
            writeln!(
                table,
                "{t:.9}\t{:.9}\t{d:.12e}\t{:.12e}\t{p:.12e}\t{:.12e}\t{:.12e}\t{:.9}",
                s.reference_s, s.rate_hz, c.band_phase_rad[0], c.band_phase_rad[1], s.coherence
            )?;
            if keep {
                let elapsed =
                    sec as f64 + nsec as f64 * 1e-9 - result.epoch + f32_at(h, 112) as f64 * 0.5;
                for (points, value) in result.applied.iter_mut().zip([
                    d * 1e9,
                    s.rate_hz,
                    p.sin().atan2(p.cos()).to_degrees(),
                    if_phase.sin().atan2(if_phase.cos()).to_degrees(),
                ]) {
                    points.push((elapsed, value));
                }
            }
        }
        Ok(())
    };
    if let Some(reference) = reference_hz {
        // QA and joint packaging share the same COR read pass.
        // Native products already contain the full station names and scan tag.
        let stem = paths[0]
            .file_stem()
            .and_then(|s| s.to_str())
            .and_then(|s| s.strip_suffix("_multiband"))
            .ok_or("unexpected native multiband product filename")?;
        let output = directory.join(format!("{stem}_mbcx.cor"));
        drop(readers);
        pack_joint(paths, &output, reference, &mut inspect)?;
    } else {
        for row in 0..rows {
            let (ha, va) = readers[0].sector()?;
            let (hb, vb) = readers[1].sector()?;
            if ha[..16] != hb[..16] || ha[112..116] != hb[112..116] {
                return Err("unequal multiband diagnostic time grids".into());
            }
            inspect(row, &ha, &[va, vb])?;
        }
    }
    text.flush()?;
    if let Some(table) = applied_text.as_mut() {
        table.flush()?;
    }
    result.write_spectrum(&directory.join(format!("{name}-spectrum.tsv")))?;
    Ok(result)
}

fn three_series(data: &[Points; 3]) -> Vec<(&str, Points)> {
    ["band 1", "band 2", "joint"]
        .into_iter()
        .zip(data.iter().cloned())
        .collect()
}

pub(super) fn write(
    directory: &Path,
    cross: &[PathBuf; 2],
    autos: &[[PathBuf; 2]; 2],
    before: Option<&[PathBuf; 2]>,
    correction: &Corrections,
    args: &Args,
    scan_epoch: &str,
    transfer: bool,
) -> Result<(), DynError> {
    let (scan_sec, _) = crate::epoch_to_yyyydddhhmmss(scan_epoch)?;
    let visibility = collect(
        cross,
        directory,
        "visibility",
        Some(correction.reference_hz),
        Some((correction, scan_sec)),
    )?;
    let acf1 = collect(
        &autos[0],
        directory,
        "acf1",
        Some(correction.reference_hz),
        None,
    )?;
    let acf2 = collect(
        &autos[1],
        directory,
        "acf2",
        Some(correction.reference_hz),
        None,
    )?;
    // The native outputs must describe precisely the same integrations.
    if [&acf1, &acf2].iter().any(|s| {
        s.rows != visibility.rows
            || s.epoch != visibility.epoch
            || (s.duration - visibility.duration).abs() > 1e-6
    }) {
        return Err("multiband auto/cross output time grids differ".into());
    }
    let uncorrected = before
        .map(|paths| collect(paths, directory, "uncorrected", None, None))
        .transpose()?;
    let header = CorReader::open(&cross[0])?.header;
    let name = |offset| {
        String::from_utf8_lossy(&header[offset..offset + 16])
            .trim_end_matches('\0')
            .to_owned()
    };
    let station1 = name(32);
    let station2 = name(80);
    figure(
        &directory.join("visibility.png"),
        &[
            (
                "Corrected complex continuum mean",
                "Time from scan start [s]",
                "Amplitude [COR units]",
                three_series(&visibility.amplitudes),
                false,
            ),
            (
                "Corrected continuum phase",
                "Time from scan start [s]",
                "Phase [deg]",
                three_series(&visibility.phases),
                true,
            ),
        ],
    )?;
    figure(
        &directory.join("band-phase.png"),
        &[(
            "Residual band 2 - band 1 continuum phase",
            "Time from scan start [s]",
            "Phase [deg]",
            vec![("corrected", visibility.relative_phase.clone())],
            true,
        )],
    )?;
    let mut amplitudes = visibility.spectral_amplitude();
    let mut phases = visibility.spectrum_points(true);
    amplitudes[0].0 = "band 1 corrected";
    amplitudes[1].0 = "band 2 corrected";
    phases[0].0 = "band 1 corrected";
    phases[1].0 = "band 2 corrected";
    if let Some(s) = &uncorrected {
        let mut a = s.spectral_amplitude();
        let mut p = s.spectrum_points(true);
        a[0].0 = "band 1 before";
        a[1].0 = "band 2 before";
        p[0].0 = "band 1 before";
        p[1].0 = "band 2 before";
        amplitudes.extend(a);
        phases.extend(p);
    }
    figure(
        &directory.join("spectrum.png"),
        &[
            (
                "Complex time-averaged cross spectrum (RF gap retained)",
                "RF [MHz]",
                "Amplitude [COR units]",
                amplitudes,
                false,
            ),
            (
                "Time-averaged cross-spectrum phase",
                "RF [MHz]",
                "Phase [deg]",
                phases,
                true,
            ),
        ],
    )?;
    figure(
        &directory.join("autocorrelation.png"),
        &[
            (
                &format!("{station1} self spectrum"),
                "RF [MHz]",
                "Power [COR units]",
                acf1.spectrum_points(false),
                false,
            ),
            (
                &format!("{station2} self spectrum"),
                "RF [MHz]",
                "Power [COR units]",
                acf2.spectrum_points(false),
                false,
            ),
        ],
    )?;
    figure(
        &directory.join("applied.png"),
        &[
            (
                "Applied residual delay at output midpoints",
                "Time from scan start [s]",
                "Delay [ns]",
                vec![("delay", visibility.applied[0].clone())],
                false,
            ),
            (
                "Applied rate at reference RF",
                "Time from scan start [s]",
                "Rate [Hz]",
                vec![("rate", visibility.applied[1].clone())],
                false,
            ),
            (
                "Applied common phase at reference RF",
                "Time from scan start [s]",
                "Phase [deg]",
                vec![("common phase", visibility.applied[2].clone())],
                true,
            ),
            (
                "Fixed instrumental band 2 - band 1 phase",
                "Time from scan start [s]",
                "Phase [deg]",
                vec![("IF phase difference", visibility.applied[3].clone())],
                true,
            ),
        ],
    )?;
    let layout = visibility.layout;
    let bandwidth = layout.channels as f64 * layout.df_hz;
    let mut parameters = BufWriter::new(File::create(directory.join("verification.txt"))?);
    writeln!(parameters, "source={}\nxml_epoch={scan_epoch}\nfirst_output_unix_s={:.9}\nsampling_hz={}\nfft_length={}", name(128), visibility.epoch, i32_at(&header, 12), i32_at(&header, 24))?;
    writeln!(parameters, "software_version={}\nbaseline={station1},{station2}\nmode={}\ncalibrators={}\nreference_hz={:.12e}\nband1_low_hz={:.12e}\nband2_low_hz={:.12e}\nbandwidth_per_band_hz={bandwidth:.12e}\noccupied_bandwidth_hz={:.12e}\nrf_gap_hz={:.12e}\nchannel_spacing_hz={:.12e}\nchannels_per_band={}\nrows={}\nintegration_total_s={:.9}\nsolution_count={}\nband1_if_phase_rad={:.12e}\nband2_if_phase_rad={:.12e}\nsearch_delay_half_width_ns={}\nsearch_rate_half_width_hz={}\nminimum_coherence={}\nphase_mode={}\nmean_weighting=equal_nonzero_channel_weight_times_integration_duration\ncoherence_is_detection_snr=false\nplots_max_points_per_series={MAX_PLOT_POINTS}",
        env!("CARGO_PKG_VERSION"), if transfer { "calibrator_transfer" } else { "joint_fit" }, args.multiband_calibrator.join(","), correction.reference_hz, layout.low_hz[0], layout.low_hz[1], 2.0 * bandwidth, layout.low_hz[1] - layout.low_hz[0] - bandwidth, layout.df_hz, layout.channels, visibility.rows, visibility.duration, correction.solutions.len(), correction.band_phase_rad[0], correction.band_phase_rad[1], args.multiband_delay_window_ns, args.multiband_rate_window_hz, args.multiband_min_coherence, args.multiband_phase_mode)?;
    for (name, s) in [
        ("corrected", Some(&visibility)),
        ("uncorrected", uncorrected.as_ref()),
        ("acf1", Some(&acf1)),
        ("acf2", Some(&acf2)),
    ] {
        if let Some(s) = s {
            let weight = s.weights.iter().sum::<f64>();
            let mean = if weight > 0.0 {
                (s.means[0] + s.means[1]) / weight
            } else {
                Complex::new(0.0, 0.0)
            };
            writeln!(parameters, "{name}_mean_real={:.12e}\n{name}_mean_imag={:.12e}\n{name}_mean_amplitude={:.12e}\n{name}_mean_phase_deg={:.9}\n{name}_mean_usable_bandwidth_hz={:.12e}", mean.re, mean.im, mean.norm(), mean.arg().to_degrees(), weight / s.duration * layout.df_hz)?;
        }
    }
    parameters.flush()?;
    println!(
        "[multiband] ACF/XCF joint products, verification TSV and PNG: {}",
        directory.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qa_averages_complex_values_with_exposure_and_excludes_zero_bins() {
        let layout = Layout {
            low_hz: [6600e6, 8192e6],
            df_hz: 8e6,
            channels: 2,
            reference_hz: 7652e6,
        };
        let mut summary = Summary::new(layout);
        let mut text = Vec::new();
        for (index, (duration, sign)) in [(1.0_f32, 1.0), (3.0_f32, -1.0)].into_iter().enumerate() {
            let mut h = [0; 128];
            h[..4].copy_from_slice(&(100_i32 + index as i32).to_le_bytes());
            h[112..116].copy_from_slice(&duration.to_le_bytes());
            let spectra = [
                vec![Complex::new(sign, 0.0), Complex::new(0.0, 0.0)],
                vec![Complex::new(0.0, sign), Complex::new(0.0, 0.0)],
            ];
            summary.add(&h, &spectra, true, &mut text).unwrap();
        }
        assert_eq!(summary.weights, [4.0, 4.0]);
        assert_eq!(summary.exposure, [vec![4.0, 0.0], vec![4.0, 0.0]]);
        let joint = (summary.means[0] + summary.means[1]) / 8.0;
        assert_eq!(joint, Complex::new(-0.25, -0.25));
        assert_eq!(summary.spectrum_points(true)[0].1[0].0, 6600.0);
        assert_eq!(summary.spectrum_points(true)[1].1[0].0, 8192.0);
        assert_eq!(summary.relative_phase[0].1, 90.0);
        assert_eq!(summary.relative_phase[1].1, 90.0);
        let range = padded_range(std::iter::once(0.0));
        assert!(range.start < 0.0 && range.end > 0.0);
    }
}
