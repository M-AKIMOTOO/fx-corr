//! Contiguous FFT-bin accumulation shared by the normal correlation path.
use num_complex::Complex;

#[inline(always)]
pub(crate) fn phase_to_f32(phase: Complex<f64>) -> Complex<f32> {
    Complex::new(phase.re as f32, phase.im as f32)
}

/// Keep the oscillator in f64: an f32 unit step can accumulate appreciable
/// amplitude and phase errors over hundreds of thousands of spectral bins.
/// Eight independent recurrences remove the per-bin dependency chain and let
/// the compiler vectorize the arithmetic without architecture-specific code.
#[inline]
pub(crate) fn accumulate_direct_acf_xcf(
    spectrum1: &[Complex<f32>],
    spectrum2: &[Complex<f32>],
    acc11: &mut [f64],
    acc12: &mut [Complex<f64>],
    acc22: &mut [f64],
    fr_mix: Complex<f32>,
    phase_start: Complex<f64>,
    phase_step: Complex<f64>,
) {
    debug_assert_eq!(spectrum1.len(), spectrum2.len());
    debug_assert_eq!(spectrum1.len(), acc11.len());
    debug_assert_eq!(spectrum1.len(), acc12.len());
    debug_assert_eq!(spectrum1.len(), acc22.len());

    const LANES: usize = 8;
    let mut phases = [phase_start; LANES];
    for lane in 1..LANES {
        phases[lane] = phases[lane - 1] * phase_step;
    }
    let step2 = phase_step * phase_step;
    let step4 = step2 * step2;
    let step8 = step4 * step4;
    let full_len = spectrum1.len() / LANES * LANES;
    for base in (0..full_len).step_by(LANES) {
        let s1 = &spectrum1[base..base + LANES];
        let s2 = &spectrum2[base..base + LANES];
        let a11 = &mut acc11[base..base + LANES];
        let a12 = &mut acc12[base..base + LANES];
        let a22 = &mut acc22[base..base + LANES];
        for lane in 0..LANES {
            a11[lane] += s1[lane].norm_sqr() as f64;
            a22[lane] += s2[lane].norm_sqr() as f64;
            let value = (s1[lane] * s2[lane].conj()) * fr_mix * phase_to_f32(phases[lane]);
            a12[lane] += Complex::new(value.re as f64, value.im as f64);
        }
        for phase in &mut phases {
            *phase *= step8;
        }
    }
    for (lane, k) in (full_len..spectrum1.len()).enumerate() {
        acc11[k] += spectrum1[k].norm_sqr() as f64;
        acc22[k] += spectrum2[k].norm_sqr() as f64;
        let value = (spectrum1[k] * spectrum2[k].conj()) * fr_mix * phase_to_f32(phases[lane]);
        acc12[k] += Complex::new(value.re as f64, value.im as f64);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accumulation_matches_independent_trigonometric_reference() {
        // Includes empty input, all tail lengths and the largest operational FFT.
        for n in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 4097, 524_289] {
            for slope in [0.0, 3.0e-6, -0.00031, 0.19] {
                let s1: Vec<_> = (0..n)
                    .map(|k| {
                        Complex::new(
                            ((17 * k) % 31) as f32 / 32.0 - 0.5,
                            ((13 * k + 5) % 29) as f32 / 32.0 - 0.5,
                        )
                    })
                    .collect();
                let s2: Vec<_> = (0..n)
                    .map(|k| {
                        Complex::new(
                            ((11 * k + 7) % 23) as f32 / 32.0 - 0.5,
                            ((19 * k + 2) % 37) as f32 / 40.0 - 0.5,
                        )
                    })
                    .collect();
                let fr = Complex::from_polar(1.0_f32, 0.37);
                let start = -0.21;
                let mut a11 = vec![2.0; n];
                let mut a12 = vec![Complex::new(1.0, -2.0); n];
                let mut a22 = vec![3.0; n];
                // Verify accumulation into existing values, not just assignment.
                for _ in 0..2 {
                    accumulate_direct_acf_xcf(
                        &s1,
                        &s2,
                        &mut a11,
                        &mut a12,
                        &mut a22,
                        fr,
                        Complex::from_polar(1.0, start),
                        Complex::from_polar(1.0, slope),
                    );
                }
                for k in 0..n {
                    let expected_phase =
                        phase_to_f32(Complex::from_polar(1.0, start + k as f64 * slope));
                    let expected = (s1[k] * s2[k].conj()) * fr * expected_phase;
                    let expected = Complex::new(
                        1.0 + 2.0 * expected.re as f64,
                        -2.0 + 2.0 * expected.im as f64,
                    );
                    assert_eq!(a11[k], 2.0 + 2.0 * s1[k].norm_sqr() as f64);
                    assert_eq!(a22[k], 3.0 + 2.0 * s2[k].norm_sqr() as f64);
                    assert!(
                        (a12[k] - expected).norm() < 5e-7,
                        "n={n} k={k} slope={slope}"
                    );
                }
            }
        }
    }
}
