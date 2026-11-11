//! Native-rate FX planning: equal time intervals and equal RF bin spacing.
//! No time-sample decimation, interpolation or additional quantization.
use crate::utils::DynError;

#[derive(Clone, Copy, Debug)]
pub struct BandPlan {
    pub sampling_hz: [f64; 2],
    pub bandwidth_hz: [f64; 2],
    pub output_sampling_hz: f64,
    pub fft: [usize; 2],
    pub ratios: [usize; 2],
}

impl BandPlan {
    pub fn new(
        fft: usize,
        sampling_hz: [f64; 2],
        bandwidth_hz: [f64; 2],
    ) -> Result<Self, DynError> {
        if fft < 2 {
            return Err("real-input correlation requires an FFT length >= 2".into());
        }
        for i in 0..2 {
            if !sampling_hz[i].is_finite()
                || sampling_hz[i] <= 0.0
                || !bandwidth_hz[i].is_finite()
                || bandwidth_hz[i] <= 0.0
                || bandwidth_hz[i] > sampling_hz[i] / 2.0
            {
                return Err(format!("antenna {}: require 0 < bandwidth <= sampling/2, in Hz (sampling={}, bandwidth={})", i + 1, sampling_hz[i], bandwidth_hz[i]).into());
            }
        }
        let output_sampling_hz = 2.0 * bandwidth_hz[0].min(bandwidth_hz[1]);
        let mut ratios = [0; 2];
        let mut native_fft = [0; 2];
        for i in 0..2 {
            let ratio = sampling_hz[i] / output_sampling_hz;
            let integer = ratio.round();
            if integer < 1.0 || (ratio - integer).abs() > 1e-10 || integer > usize::MAX as f64 {
                return Err(format!("antenna {}: sampling must be an integer multiple of twice the smaller bandwidth; got ratio {ratio}. Rational-rate resampling is not implemented", i + 1).into());
            }
            ratios[i] = integer as usize;
            native_fft[i] = fft
                .checked_mul(ratios[i])
                .ok_or("native FFT length overflow")?;
            if native_fft[i] > isize::MAX as usize {
                return Err("native FFT exceeds index range".into());
            }
        }
        Ok(Self {
            sampling_hz,
            bandwidth_hz,
            output_sampling_hz,
            fft: native_fft,
            ratios,
        })
    }

    pub fn unequal(&self) -> bool {
        self.ratios != [1, 1]
    }

    /// Subtract integer sample times before conversion to f64. Large common
    /// epoch/skip offsets must not consume carrier phase precision.
    pub fn time_difference(&self, sample1: u64, sample2: u64) -> f64 {
        let ticks =
            sample2 as i128 * self.ratios[0] as i128 - sample1 as i128 * self.ratios[1] as i128;
        ticks as f64 / (self.output_sampling_hz * self.ratios[0] as f64 * self.ratios[1] as f64)
    }

    /// Native FFT phases refer to the start of each FFT window. Different
    /// physical LO origins require exp(-2*pi*i*(LO2-LO1)*t_ref), even when
    /// their bins share an RF. Reduce integer cycles before conversion to f64.
    pub fn grid_phase(
        &self,
        antenna: usize,
        reference_sample: u64,
        lo_delta_bins: isize,
    ) -> num_complex::Complex<f32> {
        let modulus = self.fft[antenna] as i128;
        let cycles = (reference_sample as i128 * lo_delta_bins as i128).rem_euclid(modulus) as f64
            / modulus as f64;
        let p = num_complex::Complex::from_polar(1.0, -std::f64::consts::TAU * cycles);
        num_complex::Complex::new(p.re as f32, p.im as f32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn contained_512_and_16_mhz_have_the_same_time_and_frequency_grid() {
        for (fs, fft) in [
            ([1024e6, 32e6], [65536, 2048]),
            ([32e6, 1024e6], [2048, 65536]),
            ([1024e6, 1024e6], [65536, 65536]),
        ] {
            let bw = if fs[0] < fs[1] {
                [16e6, 512e6]
            } else {
                [512e6, 16e6]
            };
            let plan = BandPlan::new(2048, fs, bw).unwrap();
            assert_eq!(plan.fft, fft);
            assert_eq!(plan.output_sampling_hz, 32e6);
            for i in 0..2 {
                assert_eq!(plan.sampling_hz[i] / plan.fft[i] as f64, 15625.0);
                assert_eq!(plan.fft[i] as f64 / plan.sampling_hz[i], 64e-6);
            }
        }
    }
    #[test]
    fn native_sample_separation_preserves_small_delay_after_a_large_skip() {
        let plan = BandPlan::new(2048, [1024e6, 32e6], [512e6, 16e6]).unwrap();
        let epoch_samples = 100_000_u64 * 32_000_000;
        assert_eq!(
            plan.time_difference(epoch_samples * 32 + 17, epoch_samples + 2),
            47.0 / 1024e6
        );
        assert_eq!(plan.time_difference(epoch_samples * 32, epoch_samples), 0.0);
    }

    #[test]
    fn shifted_rf_grid_has_an_exact_phase_origin_after_a_large_skip() {
        let plan = BandPlan::new(2048, [1024e6, 32e6], [512e6, 16e6]).unwrap();
        let base = 100_000_u64 * 32_000_000;
        let wide = plan.grid_phase(0, base * 32 + 32, 512);
        let narrow = plan.grid_phase(1, base + 1, 512);
        assert!((wide - num_complex::Complex::new(0.0, -1.0)).norm() < 1e-6);
        assert_eq!(wide, narrow);
        assert!((wide.conj() - plan.grid_phase(0, base * 32 + 32, -512)).norm() < 1e-6);
    }
    #[test]
    fn invalid_sample_rates_and_bandwidths_are_rejected() {
        for (fs, bw) in [
            ([f64::NAN, 32e6], [16e6, 16e6]),
            ([32e6, 32e6], [17e6, 16e6]),
            ([48e6, 32e6], [24e6, 16e6]),
            ([32e6, 32e6], [0.0, 16e6]),
        ] {
            assert!(BandPlan::new(2048, fs, bw).is_err());
        }
    }
}
