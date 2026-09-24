use std::f64::consts::PI;
use std::sync::Arc;

use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use rustfft::{num_complex::Complex, Fft, FftPlanner};
use std::error::Error;

pub type DynError = Box<dyn Error + Send + Sync>;

pub struct FftHelper {
    len: usize,
    pub forward_r2c: Arc<dyn RealToComplex<f32>>,
    pub inverse_c2r: Arc<dyn ComplexToReal<f32>>,
    pub inverse_c2c: Arc<dyn Fft<f32>>,
}

pub struct FftScratch {
    forward_r2c: Vec<Complex<f32>>,
}

impl FftHelper {
    pub fn new(len: usize) -> Self {
        let mut planner_c2c = FftPlanner::new();
        let mut planner_r2c = RealFftPlanner::new();
        let forward_r2c = planner_r2c.plan_fft_forward(len);
        let inverse_c2r = planner_r2c.plan_fft_inverse(len);
        let inverse_c2c = planner_c2c.plan_fft_inverse(len);
        Self {
            len,
            forward_r2c,
            inverse_c2r,
            inverse_c2c,
        }
    }
    pub fn inverse_c2c(&self, spectrum: &mut [Complex<f32>]) -> Result<(), DynError> {
        if spectrum.len() != self.len {
            return Err("Spectrum length mismatch".into());
        }
        self.inverse_c2c.process(spectrum);
        let scale = 1.0_f32 / self.len as f32;
        for value in spectrum.iter_mut() {
            *value *= scale;
        }
        Ok(())
    }
    pub fn forward_r2c_process(
        &self,
        input: &mut [f32],
        output: &mut [Complex<f32>],
    ) -> Result<(), DynError> {
        if input.len() != self.len || output.len() != self.len / 2 + 1 {
            return Err("Length mismatch".into());
        }
        self.forward_r2c.process(input, output)?;
        Ok(())
    }
    pub fn make_scratch(&self) -> FftScratch {
        FftScratch {
            forward_r2c: self.forward_r2c.make_scratch_vec(),
        }
    }
    pub fn forward_r2c_process_with_scratch(
        &self,
        input: &mut [f32],
        output: &mut [Complex<f32>],
        scratch: &mut FftScratch,
    ) -> Result<(), DynError> {
        if input.len() != self.len || output.len() != self.len / 2 + 1 {
            return Err("Length mismatch".into());
        }
        self.forward_r2c
            .process_with_scratch(input, output, &mut scratch.forward_r2c)?;
        Ok(())
    }
    pub fn inverse_c2r_process(
        &self,
        spectrum: &mut [Complex<f32>],
        output: &mut [f32],
    ) -> Result<(), DynError> {
        if spectrum.len() != self.len / 2 + 1 || output.len() != self.len {
            return Err("Length mismatch".into());
        }
        self.inverse_c2r.process(spectrum, output)?;
        let scale = 1.0_f32 / self.len as f32;
        for value in output.iter_mut() {
            *value *= scale;
        }
        Ok(())
    }
}

#[allow(dead_code)]
pub fn apply_integer_sample_shift_zerofill(samples: &mut [f32], shift_samples: i64) {
    if samples.is_empty() || shift_samples == 0 {
        return;
    }
    let n = samples.len();
    let s_abs = shift_samples.unsigned_abs() as usize;
    if s_abs >= n {
        samples.fill(0.0_f32);
        return;
    }
    if shift_samples > 0 {
        // Positive delay: shift right (later in time), zero-fill head.
        samples.copy_within(..(n - s_abs), s_abs);
        samples[..s_abs].fill(0.0_f32);
    } else {
        // Negative delay: shift left (earlier in time), zero-fill tail.
        samples.copy_within(s_abs.., 0);
        samples[(n - s_abs)..].fill(0.0_f32);
    }
}

pub struct DecodePlan {
    bits: usize,
    shuffle_kind: ShuffleKind,
    input_shifts: [u32; 32],
    fast_kind: DecodeFastKind,
    levels_f32: Vec<f32>,
    packed2_bytes: Option<Box<Packed2Bytes>>,
}
// Common recorder layouts permute bytes and use the same bit permutation in
// each byte. Three 4 KiB tables cover USB and both absolute LSB parities.
struct Packed2Bytes {
    source_bytes: [usize; 4],
    tables: [[[f32; 4]; 256]; 3],
}

impl Packed2Bytes {
    fn new(shifts: &[u32; 32], levels: &[f32]) -> Option<Box<Self>> {
        let source_bytes = std::array::from_fn(|i| shifts[8 * i] as usize / 8);
        for byte in 0..4 {
            for bit in 0..8 {
                if shifts[8 * byte + bit] as usize / 8 != source_bytes[byte]
                    || shifts[8 * byte + bit] % 8 != shifts[bit] % 8
                {
                    return None;
                }
            }
        }
        let mut plan = Box::new(Self {
            source_bytes,
            tables: [[[0.0; 4]; 256]; 3],
        });
        for byte in 0..256 {
            for sample in 0..4 {
                let code = ((byte >> (shifts[2 * sample] % 8)) & 1)
                    | (((byte >> (shifts[2 * sample + 1] % 8)) & 1) << 1);
                let value = levels[code];
                plan.tables[0][byte][sample] = value;
                plan.tables[1][byte][sample] = if sample % 2 == 1 { -value } else { value };
                plan.tables[2][byte][sample] = if sample % 2 == 0 { -value } else { value };
            }
        }
        Some(plan)
    }

    fn decode(
        &self,
        raw: &[u8],
        samples: usize,
        output: &mut [f32],
        lsb: bool,
        first_odd: bool,
        identity: bool,
    ) -> usize {
        let table = &self.tables[if lsb { 1 + usize::from(first_odd) } else { 0 }];
        let [b0, b1, b2, b3] = self.source_bytes;
        let words = (samples / 16).min(raw.len() / 4);
        for (src, dst) in raw[..words * 4]
            .chunks_exact(4)
            .zip(output[..words * 16].chunks_exact_mut(16))
        {
            dst[0..4].copy_from_slice(&table[src[b0] as usize]);
            dst[4..8].copy_from_slice(&table[src[b1] as usize]);
            dst[8..12].copy_from_slice(&table[src[b2] as usize]);
            dst[12..16].copy_from_slice(&table[src[b3] as usize]);
        }
        let mut written = words * 16;
        let rest = &raw[words * 4..];
        // Nonidentity shuffle must always see a complete physical 32-bit word.
        if rest.len() >= 4 || identity {
            for &source in &self.source_bytes {
                if written == samples || source >= rest.len() {
                    break;
                }
                let count = (samples - written).min(4);
                output[written..written + count]
                    .copy_from_slice(&table[rest[source] as usize][..count]);
                written += count;
            }
        }
        written
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ShuffleKind {
    Identity,
    PairSwap,
    Generic,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DecodeFastKind {
    Generic,
    Packed1,
    Packed2,
    Packed4,
    Packed8,
}

#[inline]
fn apply_shuffle_word(mut word: u32, plan: &DecodePlan) -> u32 {
    if plan.shuffle_kind == ShuffleKind::PairSwap {
        word = ((word & 0xAAAA_AAAA) >> 1) | ((word & 0x5555_5555) << 1);
    } else if plan.shuffle_kind == ShuffleKind::Generic {
        let mut shuffled = 0u32;
        for (out_bit, &in_shift) in plan.input_shifts.iter().enumerate() {
            shuffled |= ((word >> in_shift) & 1) << out_bit;
        }
        word = shuffled;
    }
    word
}

pub fn build_decode_plan(
    bits: usize,
    shuffle_in: &[usize],
    levels: &[f64],
) -> Result<DecodePlan, DynError> {
    let mut input_shifts = [0u32; 32];
    for (idx, &mapped) in shuffle_in.iter().enumerate() {
        input_shifts[idx] = mapped as u32;
    }
    let is_identity = shuffle_in.iter().enumerate().all(|(i, &v)| i == v);
    let is_pair_swap = shuffle_in.iter().enumerate().all(|(i, &v)| (i ^ 1) == v);
    let kind = if is_identity {
        ShuffleKind::Identity
    } else if is_pair_swap {
        ShuffleKind::PairSwap
    } else {
        ShuffleKind::Generic
    };
    let fast_kind = match bits {
        1 => DecodeFastKind::Packed1,
        2 => DecodeFastKind::Packed2,
        4 => DecodeFastKind::Packed4,
        8 => DecodeFastKind::Packed8,
        _ => DecodeFastKind::Generic,
    };
    let levels_f32 = levels.iter().map(|&v| v as f32).collect::<Vec<_>>();
    let packed2_bytes = if bits == 2 {
        Packed2Bytes::new(&input_shifts, &levels_f32)
    } else {
        None
    };
    Ok(DecodePlan {
        bits,
        shuffle_kind: kind,
        input_shifts,
        fast_kind,
        levels_f32,
        packed2_bytes,
    })
}

pub fn decode_block_into_with_plan(
    raw: &[u8],
    samples: usize,
    plan: &DecodePlan,
    output: &mut [f32],
    lsb_to_usb: bool,
    first_sample_odd: bool,
) -> Result<(), DynError> {
    if output.len() < samples {
        return Err(format!("output buffer too short: {} < {}", output.len(), samples).into());
    }
    let level_map = plan.levels_f32.as_slice();

    let mut out_idx = 0usize;
    let mut odd = first_sample_odd;

    match plan.fast_kind {
        DecodeFastKind::Packed1 => {
            for chunk in raw.chunks_exact(4) {
                let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                word = apply_shuffle_word(word, plan);
                for _ in 0..32 {
                    if out_idx >= samples {
                        break;
                    }
                    let code = (word & 0x1) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd {
                        val = -val;
                    }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    word >>= 1;
                }
            }
        }
        DecodeFastKind::Packed2 => {
            if let Some(bytes) = plan.packed2_bytes.as_ref() {
                out_idx = bytes.decode(
                    raw,
                    samples,
                    output,
                    lsb_to_usb,
                    first_sample_odd,
                    plan.shuffle_kind == ShuffleKind::Identity,
                );
            } else {
                for chunk in raw.chunks_exact(4) {
                    let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                    word = apply_shuffle_word(word, plan);
                    for _ in 0..16 {
                        if out_idx >= samples {
                            break;
                        }
                        let code = (word & 0x3) as usize;
                        let mut val = level_map[code];
                        if lsb_to_usb && odd {
                            val = -val;
                        }
                        output[out_idx] = val;
                        out_idx += 1;
                        odd = !odd;
                        word >>= 2;
                    }
                }
            }
        }
        DecodeFastKind::Packed4 => {
            for chunk in raw.chunks_exact(4) {
                let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                word = apply_shuffle_word(word, plan);
                for _ in 0..8 {
                    if out_idx >= samples {
                        break;
                    }
                    let code = (word & 0xF) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd {
                        val = -val;
                    }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    word >>= 4;
                }
            }
        }
        DecodeFastKind::Packed8 => {
            for chunk in raw.chunks_exact(4) {
                let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                word = apply_shuffle_word(word, plan);
                for _ in 0..4 {
                    if out_idx >= samples {
                        break;
                    }
                    let code = (word & 0xFF) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd {
                        val = -val;
                    }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    word >>= 8;
                }
            }
        }
        DecodeFastKind::Generic => {
            let bits = plan.bits;
            let code_mask = (1u64 << bits) - 1;
            let mut acc = 0u64;
            let mut acc_bits = 0;
            for chunk in raw.chunks_exact(4) {
                let word = apply_shuffle_word(
                    u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]),
                    plan,
                );
                acc |= (word as u64) << acc_bits;
                acc_bits += 32;
                while acc_bits >= bits && out_idx < samples {
                    let code = (acc & code_mask) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd {
                        val = -val;
                    }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    acc_bits -= bits;
                    acc >>= bits;
                }
            }
        }
    }

    if out_idx != samples {
        return Err(format!("decoded {} samples, expected {}", out_idx, samples).into());
    }
    Ok(())
}

pub fn quantise_frame(
    samples: &[f32],
    bits: usize,
    levels: &[f64],
    shuffle_out: &[usize],
    output: &mut Vec<u8>,
) -> Result<(), DynError> {
    output.clear();
    let code_mask = (1u64 << bits) - 1;
    let mut bit_acc = 0u64;
    let mut acc_bits = 0;
    for &val in samples.iter() {
        let mut best_idx = 0;
        let mut min_err = f32::MAX;
        for (idx, &lv) in levels.iter().enumerate() {
            let err = (val - lv as f32).abs();
            if err < min_err {
                min_err = err;
                best_idx = idx;
            }
        }
        bit_acc |= (best_idx as u64 & code_mask) << acc_bits;
        acc_bits += bits;
        while acc_bits >= 32 {
            let word = (bit_acc & 0xFFFF_FFFF) as u32;
            let mut shuffled = 0u32;
            for (pos, &target) in shuffle_out.iter().enumerate() {
                shuffled |= ((word >> pos) & 1) << target;
            }
            output.extend_from_slice(&shuffled.to_le_bytes());
            bit_acc >>= 32;
            acc_bits -= 32;
        }
    }
    if acc_bits > 0 {
        // Zero-pad the tail to the next 32-bit word to keep raw word alignment.
        let word = (bit_acc & 0xFFFF_FFFF) as u32;
        let mut shuffled = 0u32;
        for (pos, &target) in shuffle_out.iter().enumerate() {
            shuffled |= ((word >> pos) & 1) << target;
        }
        output.extend_from_slice(&shuffled.to_le_bytes());
    }
    Ok(())
}

pub fn unwrap_phase(phase: &[f64]) -> Vec<f64> {
    let mut unwrapped = Vec::with_capacity(phase.len());
    if let Some(&first) = phase.first() {
        unwrapped.push(first);
        let mut offset = 0.0;
        for i in 1..phase.len() {
            let diff = phase[i] - phase[i - 1];
            if diff > PI {
                offset -= 2.0 * PI;
            } else if diff < -PI {
                offset += 2.0 * PI;
            }
            unwrapped.push(phase[i] + offset);
        }
    }
    unwrapped
}
pub fn safe_arg(z: &Complex<f64>) -> f64 {
    if z.re == 0.0 && z.im == 0.0 {
        0.0
    } else {
        z.arg()
    }
}

#[cfg(test)]
mod decode_tests {
    use super::*;

    #[test]
    fn packed2_byte_tables_match_independent_bit_permutation() {
        let maps: Vec<Vec<usize>> = vec![
            (0..32).collect(),
            (0..32).map(|i| i ^ 7).collect(), // native VSREC
            (0..32).map(|i| i ^ 1).collect(),
            (0..32).rev().collect(),
            (0..32).map(|i| ((i / 8 + 1) % 4) * 8 + i % 8).collect(),
            (0..32).map(|i| (5 * i + 7) % 32).collect(), // generic fallback
        ];
        let levels = [-2.25, 0.0, -0.5, 1.25];
        let raw: Vec<u8> = (0..=255u8)
            .flat_map(|v| [v, v.rotate_left(1), v.wrapping_mul(37), v ^ 0xaa])
            .collect();
        for (m, map) in maps.iter().enumerate() {
            let plan = build_decode_plan(2, map, &levels).unwrap();
            assert_eq!(plan.packed2_bytes.is_some(), m != 5);
            for lsb in [false, true] {
                for first_odd in [false, true] {
                    let mut reference = Vec::new();
                    for bytes in raw.chunks_exact(4) {
                        let word = u32::from_le_bytes(bytes.try_into().unwrap());
                        for sample in 0..16 {
                            let code = ((word >> map[2 * sample]) & 1)
                                | (((word >> map[2 * sample + 1]) & 1) << 1);
                            let value = levels[code as usize] as f32;
                            reference.push(if lsb && (first_odd ^ (sample % 2 == 1)) {
                                -value
                            } else {
                                value
                            });
                        }
                    }
                    for samples in [0, 1, 3, 4, 7, 8, 15, 16, 17, 31, 4093, 4096] {
                        let mut actual = vec![123.0; samples + 3];
                        decode_block_into_with_plan(
                            &raw,
                            samples,
                            &plan,
                            &mut actual,
                            lsb,
                            first_odd,
                        )
                        .unwrap();
                        assert_eq!(
                            actual[..samples]
                                .iter()
                                .map(|v| v.to_bits())
                                .collect::<Vec<_>>(),
                            reference[..samples]
                                .iter()
                                .map(|v| v.to_bits())
                                .collect::<Vec<_>>(),
                            "map={m} samples={samples} lsb={lsb} odd={first_odd}"
                        );
                        assert_eq!(&actual[samples..], &[123.0; 3]);
                    }
                    // The identity decoder also accepts incomplete raw words.
                    for raw_len in 1..8 {
                        let samples = if m == 0 {
                            raw_len * 4
                        } else {
                            raw_len / 4 * 16
                        };
                        let mut actual = vec![0.0; samples];
                        decode_block_into_with_plan(
                            &raw[..raw_len],
                            samples,
                            &plan,
                            &mut actual,
                            lsb,
                            first_odd,
                        )
                        .unwrap();
                        assert_eq!(actual, reference[..samples]);
                        if m != 0 && raw_len % 4 != 0 {
                            let mut short = vec![0.0; raw_len * 4];
                            assert!(decode_block_into_with_plan(
                                &raw[..raw_len],
                                raw_len * 4,
                                &plan,
                                &mut short,
                                lsb,
                                first_odd
                            )
                            .is_err());
                        }
                    }
                }
            }
        }
    }
}
