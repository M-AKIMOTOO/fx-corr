//! Run with `cargo bench --offline --bench correlation` on an idle host.
#[allow(unused_imports)]
#[path = "../src/corr_kernel.rs"]
mod corr_kernel;
#[allow(dead_code, unused_imports)]
#[path = "../src/geom.rs"]
mod geom;
mod utils {
    pub type DynError = Box<dyn std::error::Error + Send + Sync>;
}
use corr_kernel::accumulate_direct_acf_xcf;
use num_complex::Complex;
use std::hint::black_box;
use std::time::Instant;

// Snapshot of the previous scalar f32 recurrence for a reproducible baseline.
fn legacy(
    s1: &[Complex<f32>],
    s2: &[Complex<f32>],
    a11: &mut [f64],
    a12: &mut [Complex<f64>],
    a22: &mut [f64],
    fr: Complex<f32>,
    mut phase: Complex<f32>,
    step: Complex<f32>,
) {
    for k in 0..s1.len() {
        a11[k] += s1[k].norm_sqr() as f64;
        a22[k] += s2[k].norm_sqr() as f64;
        let value = (s1[k] * s2[k].conj()) * fr * phase;
        a12[k] += Complex::new(value.re as f64, value.im as f64);
        phase *= step;
    }
}

fn main() {
    benchmark_delay_model();
    println!("bins legacy_ns/bin new_ns/bin speedup legacy_max_error new_max_error");
    for n in [169, 1025, 3200, 4097, 32769, 524289] {
        let s1 = vec![Complex::new(1.0_f32, 0.0); n];
        let s2 = s1.clone();
        let mut a11 = vec![0.0; n];
        let mut a12 = vec![Complex::new(0.0, 0.0); n];
        let mut a22 = a11.clone();
        let fr = Complex::new(1.0_f32, 0.0);
        let start = -0.21_f64;
        let slope = 1.3 / n as f64;
        let phase0 = Complex::from_polar(1.0_f64, start);
        let step = Complex::from_polar(1.0_f64, slope);
        let old_phase = Complex::from_polar(1.0_f32, start as f32);
        let old_step = Complex::from_polar(1.0_f32, slope as f32);
        let iterations = (8_000_000 / n).max(8);
        let mut timings = [Vec::new(), Vec::new()];
        // Alternate order to reduce thermal/frequency bias. Discard warm-up.
        for round in 0..8 {
            for choice in [round % 2, 1 - round % 2] {
                let t = Instant::now();
                for _ in 0..iterations {
                    if choice == 0 {
                        legacy(
                            black_box(&s1),
                            black_box(&s2),
                            black_box(&mut a11),
                            black_box(&mut a12),
                            black_box(&mut a22),
                            black_box(fr),
                            black_box(old_phase),
                            black_box(old_step),
                        );
                    } else {
                        accumulate_direct_acf_xcf(
                            black_box(&s1),
                            black_box(&s2),
                            black_box(&mut a11),
                            black_box(&mut a12),
                            black_box(&mut a22),
                            black_box(fr),
                            black_box(phase0),
                            black_box(step),
                        );
                    }
                }
                let ns = t.elapsed().as_secs_f64() * 1e9 / (iterations * n) as f64;
                if round > 0 {
                    timings[choice].push(ns);
                }
            }
        }
        let mut errors = [0.0_f64; 2];
        for choice in 0..2 {
            a12.fill(Complex::new(0.0, 0.0));
            if choice == 0 {
                legacy(
                    &s1, &s2, &mut a11, &mut a12, &mut a22, fr, old_phase, old_step,
                );
            } else {
                accumulate_direct_acf_xcf(&s1, &s2, &mut a11, &mut a12, &mut a22, fr, phase0, step);
            }
            for (k, value) in a12.iter().enumerate() {
                let reference = Complex::from_polar(1.0, start + slope * k as f64);
                errors[choice] = errors[choice].max((*value - reference).norm());
            }
        }
        for samples in &mut timings {
            samples.sort_by(f64::total_cmp);
        }
        let old = timings[0][3];
        let new = timings[1][3];
        println!(
            "{n} {old:.3} {new:.3} {:.2} {:.3e} {:.3e}",
            old / new,
            errors[0],
            errors[1]
        );
    }
}

fn benchmark_delay_model() {
    use geom::{EarthOrientation, GeometricDelayMode, SourceVectorMode};
    let epoch = 60977.34375;
    let a1 = geom::YAMAGU32_ECEF;
    let a2 = [-3961788.974, 3243597.492, 3790597.692];
    let eop = EarthOrientation::default();
    println!("delay model: microseconds per sample (median of 3 alternating runs)");
    for mode in [GeometricDelayMode::Anchored, GeometricDelayMode::Geocentric] {
        let mut timings = [Vec::new(), Vec::new()];
        for round in 0..4 {
            for choice in [round % 2, 1 - round % 2] {
                let begin = Instant::now();
                for i in 0..256 {
                    let mjd = black_box(epoch + i as f64 / 86400.0);
                    let delay = if choice == 0 {
                        geom::calculate_geometric_delay_and_derivatives_full_with_eop(
                            a1,
                            a2,
                            4.594776025749384,
                            -0.2282965720462692,
                            mjd,
                            epoch,
                            eop,
                            mode,
                            SourceVectorMode::MeanGast,
                        )
                        .2
                    } else {
                        geom::calculate_geometric_delay_full_with_eop(
                            a1,
                            a2,
                            4.594776025749384,
                            -0.2282965720462692,
                            mjd,
                            epoch,
                            eop,
                            mode,
                            SourceVectorMode::MeanGast,
                        )
                    };
                    black_box(delay);
                }
                if round > 0 {
                    timings[choice].push(begin.elapsed().as_secs_f64() * 1e6 / 256.0);
                }
            }
        }
        for samples in &mut timings {
            samples.sort_by(f64::total_cmp);
        }
        println!(
            "{mode:?}: full={:.3} delay_only={:.3} speedup={:.2}",
            timings[0][1],
            timings[1][1],
            timings[0][1] / timings[1][1]
        );
    }
}
