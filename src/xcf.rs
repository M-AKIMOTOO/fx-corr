use std::cmp::Ordering;
use std::path::Path;
use num_complex::Complex;

use crate::plot::{plot_multi_series_f64_x, plot_series_f64_x, plot_series_with_x, BLUE, GREEN};
use crate::utils::{safe_arg, unwrap_phase, DynError, FftHelper};
use std::f64::consts::PI;

pub struct CrossSpectrumResult {
    pub delay_seconds_from_phase: Option<f64>,
}

pub fn finalize_cross_spectrum(
    accumulated_cross_spec: &mut [Complex<f64>],
    helper: &FftHelper,
    fft_len: usize,
    half_spec_len: usize,
    lags: &[i32],
    sampling_rate_hz: f64,
    sky_low_mhz: f64,
    invert_lag_axis: bool,
    output_dir: &Path,
    is_lsb: bool,
    write_plots: bool,
) -> Result<CrossSpectrumResult, DynError> {
    // Zero out DC component (now at index 0) and the original DC component (shifted to Nyquist).
    accumulated_cross_spec[0] = Complex::new(0.0, 0.0);
    if fft_len / 2 < accumulated_cross_spec.len() {
        accumulated_cross_spec[fft_len / 2] = Complex::new(0.0, 0.0);
    }

    // --- Calculate Spectra for Phase/Mag Plots ---
    let mut mag_spectrum: Vec<f64> = Vec::with_capacity(accumulated_cross_spec.len());
    let mut phase_spectrum_rad: Vec<f64> = Vec::with_capacity(accumulated_cross_spec.len());
    for &bin in accumulated_cross_spec.iter() {
        mag_spectrum.push(bin.norm());
        phase_spectrum_rad.push(safe_arg(&bin));
    }
    if mag_spectrum.len() > 1 {
        let sum_of_rest: f64 = mag_spectrum.iter().skip(1).sum();
        mag_spectrum[0] = sum_of_rest / (mag_spectrum.len() - 1) as f64 / 2.0;
    }

    // --- Unwrap Phase and Fit for Delay ---
    let unwrapped_phase_spectrum_rad = if phase_spectrum_rad.len() > 1 {
        let mut unwrapped = Vec::with_capacity(phase_spectrum_rad.len());
        unwrapped.push(phase_spectrum_rad[0]);
        let tail_unwrapped = unwrap_phase(&phase_spectrum_rad[1..]);
        unwrapped.extend(tail_unwrapped);
        unwrapped
    } else {
        phase_spectrum_rad.clone()
    };

    let bw_mhz = sampling_rate_hz / 2.0 / 1e6;
    let df_hz = sampling_rate_hz / fft_len as f64;
    let freqs_mhz: Vec<f64> = (0..accumulated_cross_spec.len())
        .map(|i| {
            if is_lsb { (sky_low_mhz + bw_mhz) - (i as f64 * df_hz) / 1e6 }
            else { sky_low_mhz + (i as f64 * df_hz) / 1e6 }
        })
        .collect();

    let freqs_hz: Vec<f64> = freqs_mhz.iter().map(|f| f * 1.0e6).collect();

    let (phase_fit_line_rad, delay_seconds_from_phase) = if freqs_hz.len()
        == unwrapped_phase_spectrum_rad.len()
        && unwrapped_phase_spectrum_rad.len() > 1
    {
        let start_index = 1;
        let mut end_index = freqs_hz.len();
        if fft_len % 2 == 0 && end_index > 1 {
            end_index -= 1;
        }

        if end_index > start_index {
            let fit_freqs = &freqs_hz[start_index..end_index];
            let fit_phases = &unwrapped_phase_spectrum_rad[start_index..end_index];

            let n = fit_freqs.len() as f64;
            let sum_x: f64 = fit_freqs.iter().sum();
            let sum_y: f64 = fit_phases.iter().sum();
            let sum_xy: f64 = fit_freqs.iter().zip(fit_phases).map(|(&x, &y)| x * y).sum();
            let sum_x2: f64 = fit_freqs.iter().map(|&x| x * x).sum();
            let denominator = n * sum_x2 - sum_x * sum_x;

            if denominator.abs() > 1e-9 {
                let slope = (n * sum_xy - sum_x * sum_y) / denominator;
                let intercept = (sum_y - slope * sum_x) / n;
                let delay = -slope / (2.0 * PI);
                let fit_line: Vec<f64> = freqs_hz
                    .iter()
                    .map(|&freq| slope * freq + intercept)
                    .collect();
                (fit_line, Some(delay))
            } else {
                (vec![], None)
            }
        } else {
            (vec![], None)
        }
    } else {
        (vec![], None)
    };

    if write_plots {
        // --- Plotting Section ---
        let unwrapped_phase_spectrum_deg: Vec<f64> = unwrapped_phase_spectrum_rad
            .iter()
            .map(|r| r.to_degrees())
            .collect();
        let phase_fit_line_deg: Vec<f64> = phase_fit_line_rad.iter().map(|r| r.to_degrees()).collect();
        let phase_spectrum_deg: Vec<f64> = phase_spectrum_rad.iter().map(|r| r.to_degrees()).collect();

        let phase_plot_freqs: &[f64] = if freqs_mhz.len() > 1 {
            &freqs_mhz[1..]
        } else {
            &freqs_mhz
        };
        let phase_plot_deg: Vec<f64> = if phase_spectrum_deg.len() > 1 {
            phase_spectrum_deg[1..].to_vec()
        } else {
            phase_spectrum_deg.clone()
        };
        let unwrapped_phase_plot_deg: Vec<f64> = if unwrapped_phase_spectrum_deg.len() > 1 {
            unwrapped_phase_spectrum_deg[1..].to_vec()
        } else {
            unwrapped_phase_spectrum_deg.clone()
        };
        let phase_fit_plot_deg: Vec<f64> = if phase_fit_line_deg.len() > 1 {
            phase_fit_line_deg[1..].to_vec()
        } else {
            phase_fit_line_deg.clone()
        };

        println!("\n[plot] Writing correlation plots...");
        let mag_upper = mag_spectrum
            .iter()
            .copied()
            .fold(0.0_f64, |acc, v| if v > acc { v } else { acc });
        let mag_upper = if mag_upper > 0.0 {
            mag_upper * 1.05
        } else {
            1.0
        };

        let (peak_freq_mhz, _peak_val) = freqs_mhz
            .iter()
            .zip(mag_spectrum.iter())
            .max_by(|(_, &a), (_, &b)| a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((&0.0, &0.0));

        let mag_path = output_dir.join(format!("xcf_spectrum_magnitude_fft{}.png", fft_len));
        plot_series_f64_x(
            &freqs_mhz,
            &mag_spectrum,
            &format!("Cross-Spectrum Magnitude (peak {:.3} MHz)", peak_freq_mhz),
            &mag_path.to_string_lossy(),
            "Frequency (MHz)",
            "Magnitude",
            Some((0.0, mag_upper)),
            &format!("Cross-Spectrum (peak {:.3} MHz)", peak_freq_mhz),
        )?;

        let phs_path = output_dir.join(format!("xcf_spectrum_phase_fft{}.png", fft_len));
        plot_series_f64_x(
            phase_plot_freqs,
            &phase_plot_deg,
            "Cross-Spectrum Phase",
            &phs_path.to_string_lossy(),
            "Frequency (MHz)",
            "Phase (degrees)",
            None,
            "Cross-Spectrum Phase",
        )?;

        let unw_path = output_dir.join(format!("xcf_spectrum_unwrapped_phase_fft{}.png", fft_len));
        if !phase_fit_plot_deg.is_empty() {
            let fit_label = if let Some(delay) = delay_seconds_from_phase {
                format!("Linear fit (delay {:.3e} s)", delay)
            } else {
                "Linear fit".to_string()
            };
            plot_multi_series_f64_x(
                phase_plot_freqs,
                &[
                    (
                        &unwrapped_phase_plot_deg,
                        &BLUE,
                        "Unwrapped phase (DC removed)",
                    ),
                    (&phase_fit_plot_deg, &GREEN, &fit_label),
                ],
                "Cross-Spectrum Unwrapped Phase (fit excludes DC/Nyquist)",
                &unw_path.to_string_lossy(),
                "Frequency (MHz)",
                "Unwrapped Phase (degrees)",
                None,
            )?;
        } else {
            plot_series_f64_x(
                phase_plot_freqs,
                &unwrapped_phase_plot_deg,
                "Cross-Spectrum Unwrapped Phase (DC removed)",
                &unw_path.to_string_lossy(),
                "Frequency (MHz)",
                "Unwrapped Phase (degrees)",
                None,
                "Unwrapped Phase",
            )?;
        }

        // --- CCF Calculation and Plotting ---
        println!("\nIntegration finished. Finalizing correlation...");
        let mut full_cross_spec = vec![Complex::new(0.0_f32, 0.0_f32); fft_len];
        for i in 0..half_spec_len {
            full_cross_spec[i] = Complex::new(accumulated_cross_spec[i].re as f32, accumulated_cross_spec[i].im as f32);
        }
        if fft_len > 1 {
            let mirror_limit = if fft_len % 2 == 0 {
                half_spec_len - 1
            } else {
                half_spec_len
            };
            for mirrored in 1..mirror_limit {
                let v = accumulated_cross_spec[mirrored].conj();
                full_cross_spec[fft_len - mirrored] = Complex::new(v.re as f32, v.im as f32);
            }
        }

        helper.inverse_c2c(&mut full_cross_spec)?; // IFFT
        let cross_corr_mag: Vec<f64> = full_cross_spec.iter().map(|c| c.norm() as f64).collect();

        let shift = fft_len / 2;
        let cross_corr_mag_shifted: Vec<f64> = cross_corr_mag
            .iter()
            .cycle()
            .skip(shift)
            .take(fft_len)
            .copied()
            .collect();

        let lags_display: Vec<i32> = lags
            .iter()
            .map(|&lag| if invert_lag_axis { -lag } else { lag })
            .collect();

        let mut peaks: Vec<(i32, f64)> = lags_display
            .iter()
            .zip(cross_corr_mag_shifted.iter())
            .map(|(&l, &v)| (l, v))
            .collect();
        peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));

        let top_points: Vec<(i32, f64)> = peaks.iter().take(10).cloned().collect();
        let highlight = if top_points.is_empty() {
            None
        } else {
            Some(top_points.as_slice())
        };

        let corr_path = output_dir.join(format!("xcf_correlation_magnitude_fft{}.png", fft_len));
        let y_upper = cross_corr_mag_shifted
            .iter()
            .copied()
            .fold(0.0_f64, |acc, v| if v > acc { v } else { acc });
        let y_upper = if y_upper > 0.0 { y_upper * 1.05 } else { 1.0 };

        plot_series_with_x(
            &lags_display,
            &[(&cross_corr_mag_shifted, &BLUE)],
            "Cross-Correlation Magnitude",
            &corr_path.to_string_lossy(),
            "Lag (samples, fftshifted)",
            "Magnitude",
            Some((0.0, y_upper)),
            highlight,
        )?;

        println!("Top correlation magnitude peaks (lag, value):");
        for (lag, val) in peaks.iter().take(10) {
            println!("  lag {:>5}: {:.6}", lag, val);
        }

        // --- SNR Calculation ---
        if let Some(&(peak_lag, peak_val)) = peaks.first() {
            let exclusion_zone = 20;
            let noise_mags: Vec<f64> = lags_display
                .iter()
                .zip(cross_corr_mag_shifted.iter())
                .filter(|(&l, _)| (l - peak_lag).abs() > exclusion_zone)
                .map(|(_, &v)| v)
                .collect();
            if !noise_mags.is_empty() {
                let mut sorted_noise_mags = noise_mags.clone();
                sorted_noise_mags.sort_by(|a, b| a.total_cmp(b));
                let mid = sorted_noise_mags.len() / 2;
                let median = if sorted_noise_mags.len() % 2 == 0 {
                    (sorted_noise_mags[mid - 1] + sorted_noise_mags[mid]) / 2.0
                } else {
                    sorted_noise_mags[mid]
                };
                let mean = noise_mags.iter().sum::<f64>() / noise_mags.len() as f64;
                let variance = noise_mags
                    .iter()
                    .map(|value| {
                        let diff = mean - value;
                        diff * diff
                    })
                    .sum::<f64>()
                    / noise_mags.len() as f64;
                let std_dev = variance.sqrt();
                println!("Median correlation magnitude: {:.6}", median);
                if median > 0.0 {
                    println!("Peak/Median ratio: {:.6}", peak_val / median);
                }
                if std_dev > 0.0 {
                    println!("SNR (Peak/StdDev): {:.6}", peak_val / std_dev);
                }
            }
        }

        println!();
    }

    Ok(CrossSpectrumResult {
        delay_seconds_from_phase,
    })
}
