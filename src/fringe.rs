use num_complex::Complex;
use crate::utils::{DynError, FftHelper};

pub fn perform_fringe_search(
    fringe_data: Vec<Vec<Complex<f64>>>,
    helper: &FftHelper,
    _sampling_rate_hz: f64,
    seconds_per_frame: f64,
    half_spec_len: usize,
) -> Result<(), DynError> {
    println!("[info] Performing 2D FFT for fringe search...");

    let num_frames = fringe_data.len();
    if num_frames == 0 {
        println!("[warn] No frames available for fringe search.");
        return Ok(());
    }
    let fft_len = helper.len();

    // Step 1: Frequency-wise IFFT (to get time-lag data)
    let mut lag_time_data: Vec<Vec<Complex<f64>>> = Vec::with_capacity(num_frames);
    for frame_cross_spec_half in fringe_data {
        let mut full_cross_spec = vec![Complex::new(0.0, 0.0); fft_len];
        full_cross_spec[..half_spec_len].copy_from_slice(&frame_cross_spec_half);
        if fft_len > 1 {
            let mirror_limit = if fft_len % 2 == 0 { half_spec_len - 1 } else { half_spec_len };
            for mirrored in 1..mirror_limit {
                full_cross_spec[fft_len - mirrored] = frame_cross_spec_half[mirrored].conj();
            }
        }
        helper.inverse_c2c(&mut full_cross_spec)?;
        lag_time_data.push(full_cross_spec);
    }

    // Step 2: Time-wise FFT (to get rate data)
    let mut lag_rate_data: Vec<Vec<Complex<f64>>> = Vec::with_capacity(fft_len);
    for lag_idx in 0..fft_len {
        let mut time_series_for_lag: Vec<Complex<f64>> = Vec::with_capacity(num_frames);
        for frame_idx in 0..num_frames {
            time_series_for_lag.push(lag_time_data[frame_idx][lag_idx]);
        }
        // Perform FFT on this time series
        // For simplicity, we can use inverse_c2c as it also performs FFT and scaling
        // A forward_c2c would be more appropriate conceptually, but the magnitude is what we care about. 
        // If rustfft had a direct forward_c2c on FftHelper, we'd use that.
        // For now, we'll use inverse_c2c and adjust interpretation if needed.
        let mut planner = rustfft::FftPlanner::new();
        let fft = planner.plan_fft_forward(num_frames);
        fft.process(&mut time_series_for_lag);
        let scale = 1.0 / num_frames as f64;
        for val in time_series_for_lag.iter_mut() {
            *val *= scale;
        }
        lag_rate_data.push(time_series_for_lag);
    }

    println!("[debug] lag_rate_data dimensions: {} (lag) x {} (rate)", lag_rate_data.len(), lag_rate_data.first().map_or(0, |v| v.len()));

    // Step X: Save lag_rate_data to file for Python plotting
    {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create("fringe_search_data.txt")?;
        writeln!(file, "# Lag_Index Rate_Index Magnitude")?;
        for (lag_idx, rate_data_for_lag) in lag_rate_data.iter().enumerate() {
            for (rate_idx, &complex_val) in rate_data_for_lag.iter().enumerate() {
                let magnitude = complex_val.norm();
                writeln!(file, "{} {} {:.6e}", lag_idx, rate_idx, magnitude)?;
            }
        }
        println!("[info] Fringe search data saved to fringe_search_data.txt");
    }

    // Step 3: Peak Detection
    let mut max_magnitude = 0.0;
    let mut peak_lag_idx = 0;
    let mut peak_rate_idx = 0;

    for (lag_idx, rate_data_for_lag) in lag_rate_data.iter().enumerate() {
        for (rate_idx, &complex_val) in rate_data_for_lag.iter().enumerate() {
            let magnitude = complex_val.norm();
            if magnitude > max_magnitude {
                max_magnitude = magnitude;
                peak_lag_idx = lag_idx;
                peak_rate_idx = rate_idx;
            }
        }
    }

    // Step 4: Convert Indices to Physical Units
    // Lag conversion (accounting for FFT shift)
    let shifted_peak_lag_idx = (peak_lag_idx + fft_len / 2) % fft_len;
    let peak_lag_samples = if shifted_peak_lag_idx >= fft_len / 2 {
        (shifted_peak_lag_idx as isize - fft_len as isize) as f64
    } else {
        shifted_peak_lag_idx as f64
    };

    // Rate conversion (accounting for FFT shift for rate axis)
    let rate_resolution = 1.0 / (num_frames as f64 * seconds_per_frame);
    let shifted_peak_rate_idx = (peak_rate_idx + num_frames / 2) % num_frames;
    let peak_rate_hz = if shifted_peak_rate_idx >= num_frames / 2 {
        (shifted_peak_rate_idx as isize - num_frames as isize) as f64 * rate_resolution
    } else {
        shifted_peak_rate_idx as f64 * rate_resolution
    };

    println!("[info] Fringe Search Peak: Lag = {:.3} samples, Rate = {:.3} Hz", peak_lag_samples, peak_rate_hz);

    Ok(())
}