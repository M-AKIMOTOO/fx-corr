use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use num_complex::Complex;

use crate::plot::{plot_series_f64_x, plot_series_with_x};
use crate::utils::DynError;

fn safe_token(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' || ch == '.' {
            out.push(ch);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "fringe".to_string()
    } else {
        out
    }
}

pub fn write_ql_fringe_products(
    out_dir: &Path,
    tag: &str,
    label: &str,
    start_offset_s: f64,
    integ_s: f64,
    spectrum: &[Complex<f64>],
    obs_mhz: f64,
    df_hz: f64,
) -> Result<(), DynError> {
    if spectrum.is_empty() {
        return Ok(());
    }
    std::fs::create_dir_all(out_dir)?;

    let label = safe_token(label);
    let start_tag = format!("t{:09.3}", start_offset_s).replace('.', "p");
    let base = format!("{}_{}_{}s", tag, label, start_tag);
    let nbin = spectrum.len();
    let df_mhz = df_hz / 1.0e6;
    let freqs_mhz: Vec<f64> = (0..nbin).map(|k| obs_mhz + k as f64 * df_mhz).collect();
    let amp: Vec<f64> = spectrum.iter().map(|z| z.norm()).collect();
    let phase_deg: Vec<f64> = spectrum.iter().map(|z| z.arg().to_degrees()).collect();

    let spec_tsv = out_dir.join(format!("{}_fringe_spectrum.tsv", base));
    {
        let mut w = BufWriter::new(File::create(&spec_tsv)?);
        writeln!(w, "# freq_mhz\tre\tim\tamp\tphase_deg")?;
        for k in 0..nbin {
            writeln!(
                w,
                "{:.9}\t{:.9e}\t{:.9e}\t{:.9e}\t{:.6}",
                freqs_mhz[k], spectrum[k].re, spectrum[k].im, amp[k], phase_deg[k]
            )?;
        }
    }

    plot_series_f64_x(
        &freqs_mhz,
        &amp,
        "QL XCF Spectrum Amplitude",
        &out_dir
            .join(format!("{}_fringe_spectrum_amp.png", base))
            .to_string_lossy(),
        "Frequency (MHz)",
        "Amplitude",
        None,
        "XCF amplitude",
    )?;

    let fft_len = nbin * 2;
    let mut lag_spec = vec![Complex::new(0.0_f64, 0.0_f64); fft_len];
    lag_spec[..nbin].copy_from_slice(spectrum);
    for k in 1..nbin {
        lag_spec[fft_len - k] = spectrum[k].conj();
    }
    let mut planner = rustfft::FftPlanner::<f64>::new();
    let fft = planner.plan_fft_inverse(fft_len);
    fft.process(&mut lag_spec);
    let scale = 1.0 / fft_len as f64;
    for z in &mut lag_spec {
        *z *= scale;
    }

    let half = fft_len / 2;
    let mut lags = Vec::with_capacity(fft_len);
    let mut lag_amp = Vec::with_capacity(fft_len);
    for i in 0..fft_len {
        let idx = (i + half) % fft_len;
        lags.push(i as i32 - half as i32);
        lag_amp.push(lag_spec[idx].norm());
    }

    let lag_tsv = out_dir.join(format!("{}_fringe_lag.tsv", base));
    {
        let mut w = BufWriter::new(File::create(&lag_tsv)?);
        writeln!(w, "# start_offset_s\t{:.9}", start_offset_s)?;
        writeln!(w, "# integration_s\t{:.9}", integ_s)?;
        writeln!(w, "# lag_sample\tre\tim\tamp")?;
        for i in 0..fft_len {
            let idx = (i + half) % fft_len;
            writeln!(
                w,
                "{}\t{:.9e}\t{:.9e}\t{:.9e}",
                i as i32 - half as i32,
                lag_spec[idx].re,
                lag_spec[idx].im,
                lag_spec[idx].norm()
            )?;
        }
    }

    plot_series_with_x(
        &lags,
        &[(&lag_amp, &crate::plot::BLUE)],
        "QL XCF Lag Fringe",
        &out_dir
            .join(format!("{}_fringe_lag_amp.png", base))
            .to_string_lossy(),
        "Lag (samples)",
        "Amplitude",
        None,
        None,
    )?;

    println!(
        "[info] QL fringe: wrote {} and {} ({:.3}s integration)",
        spec_tsv.display(),
        lag_tsv.display(),
        integ_s
    );
    Ok(())
}
