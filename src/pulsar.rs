use num_complex::Complex;
use roxmltree::Node;

use crate::geom;
use crate::utils::DynError;

const DM_DELAY_S_MHZ: f64 = 4.148_808e3;

#[derive(Clone, Debug)]
pub struct PulsarConfig {
    pub name: Option<String>,
    pub epoch: Option<String>,
    pub period_s: f64,
    pub pdot: f64,
    pub pddot: f64,
    pub bins: usize,
    pub dm: f64,
    pub dedisperse: bool,
    pub time_correction: String,
}

#[derive(Clone, Debug)]
pub struct PulsarRuntime {
    cfg: PulsarConfig,
    phase_epoch_offset_s: f64,
    obs_mhz: f64,
    df_mhz: f64,
    ref_mhz: f64,
}

pub struct FoldAccum {
    bins: usize,
    nfreq: usize,
    want_phased: bool,
    phased: Vec<Vec<f64>>,
    p11: Vec<Vec<f64>>,
    p12: Vec<Vec<Complex<f64>>>,
    p22: Vec<Vec<f64>>,
    counts: Vec<Vec<u32>>,
}

fn is_tag(node: Node<'_, '_>, tag: &str) -> bool {
    node.is_element() && node.tag_name().name().eq_ignore_ascii_case(tag)
}

fn child_text<'a>(node: Node<'a, 'a>, tag: &str) -> Option<&'a str> {
    node.children()
        .find(|n| is_tag(*n, tag))
        .and_then(|n| n.text())
        .map(str::trim)
}

fn parse_f64(raw: &str, label: &str) -> Result<f64, DynError> {
    raw.trim()
        .parse::<f64>()
        .map_err(|e| format!("invalid pulsar/{label}='{raw}': {e}").into())
}

fn parse_bool(raw: &str, label: &str) -> Result<bool, DynError> {
    let v = raw.trim().to_ascii_lowercase();
    match v.as_str() {
        "" | "0" | "false" | "off" | "no" => Ok(false),
        "1" | "true" | "on" | "yes" => Ok(true),
        _ => Err(format!("invalid pulsar/{label}='{raw}' (use true/false)").into()),
    }
}

pub fn parse_xml_node(node: Option<Node<'_, '_>>) -> Result<Option<PulsarConfig>, DynError> {
    let Some(node) = node else {
        return Ok(None);
    };
    let period_s = parse_f64(
        child_text(node, "period").ok_or("pulsar/period missing")?,
        "period",
    )?;
    if !period_s.is_finite() || period_s <= 0.0 {
        return Err("pulsar/period must be positive finite seconds".into());
    }

    let bins = match child_text(node, "bins") {
        Some(v) if !v.is_empty() => v
            .parse::<usize>()
            .map_err(|e| format!("invalid pulsar/bins='{v}': {e}"))?,
        _ => 0,
    };
    if bins == 0 {
        return Err("pulsar/bins must be >= 1".into());
    }

    let cfg = PulsarConfig {
        name: child_text(node, "name").map(str::to_string),
        epoch: child_text(node, "epoch").map(str::to_string),
        period_s,
        pdot: child_text(node, "pdot")
            .map(|v| parse_f64(v, "pdot"))
            .transpose()?
            .unwrap_or(0.0),
        pddot: child_text(node, "pddot")
            .map(|v| parse_f64(v, "pddot"))
            .transpose()?
            .unwrap_or(0.0),
        bins,
        dm: child_text(node, "dm")
            .map(|v| parse_f64(v, "dm"))
            .transpose()?
            .unwrap_or(0.0),
        dedisperse: child_text(node, "dedisperse")
            .map(|v| parse_bool(v, "dedisperse"))
            .transpose()?
            .unwrap_or(false),
        time_correction: child_text(node, "time-correction")
            .or_else(|| child_text(node, "time_correction"))
            .unwrap_or("topocentric")
            .to_ascii_lowercase(),
    };
    if cfg.time_correction != "topocentric" {
        return Err(format!(
            "pulsar/time-correction='{}' is reserved for future use; only topocentric is implemented",
            cfg.time_correction
        )
        .into());
    }
    Ok(Some(cfg))
}

impl PulsarRuntime {
    pub fn new(
        cfg: PulsarConfig,
        process_epoch: &str,
        obs_mhz: f64,
        df_hz: f64,
        nfreq: usize,
    ) -> Result<Self, DynError> {
        let phase_epoch_offset_s = if let Some(epoch) = cfg.epoch.as_deref() {
            let process_mjd = geom::parse_epoch_to_mjd(process_epoch)?;
            let pulse_mjd = geom::parse_epoch_to_mjd(epoch)?;
            (process_mjd - pulse_mjd) * 86400.0
        } else {
            0.0
        };
        Ok(Self {
            cfg,
            phase_epoch_offset_s,
            obs_mhz,
            df_mhz: df_hz / 1.0e6,
            ref_mhz: obs_mhz + (nfreq.saturating_sub(1) as f64) * df_hz / 1.0e6,
        })
    }

    pub fn config(&self) -> &PulsarConfig {
        &self.cfg
    }

    pub fn bins(&self) -> usize {
        self.cfg.bins
    }

    pub fn phase_bin(&self, t_since_process_s: f64, freq_bin: usize) -> usize {
        let mut dt = self.phase_epoch_offset_s + t_since_process_s;
        if self.cfg.dedisperse && self.cfg.dm != 0.0 {
            let f_mhz = self.obs_mhz + freq_bin as f64 * self.df_mhz;
            if f_mhz > 0.0 && self.ref_mhz > 0.0 {
                let delay_s = DM_DELAY_S_MHZ
                    * self.cfg.dm
                    * (1.0 / (f_mhz * f_mhz) - 1.0 / (self.ref_mhz * self.ref_mhz));
                dt -= delay_s;
            }
        }

        let p = self.cfg.period_s;
        let f0 = 1.0 / p;
        let fdot = -self.cfg.pdot / (p * p);
        let fddot = 2.0 * self.cfg.pdot * self.cfg.pdot / (p * p * p) - self.cfg.pddot / (p * p);
        let phase_cycles = f0 * dt + 0.5 * fdot * dt * dt + (1.0 / 6.0) * fddot * dt * dt * dt;
        let frac = phase_cycles.rem_euclid(1.0);
        let bin = (frac * self.cfg.bins as f64).floor() as usize;
        bin.min(self.cfg.bins - 1)
    }
}

impl FoldAccum {
    pub fn new(bins: usize, nfreq: usize, want_phased: bool) -> Self {
        let real_grid = || vec![vec![0.0_f64; nfreq]; bins];
        let complex_grid = || vec![vec![Complex::new(0.0_f64, 0.0_f64); nfreq]; bins];
        Self {
            bins,
            nfreq,
            want_phased,
            phased: if want_phased { real_grid() } else { Vec::new() },
            p11: real_grid(),
            p12: complex_grid(),
            p22: real_grid(),
            counts: vec![vec![0_u32; nfreq]; bins],
        }
    }

    pub fn add_values(
        &mut self,
        runtime: &PulsarRuntime,
        t_since_process_s: f64,
        freq_bin: usize,
        phased: Option<f64>,
        p11: f64,
        p12: Complex<f64>,
        p22: f64,
    ) {
        if freq_bin >= self.nfreq {
            return;
        }
        let bin = runtime.phase_bin(t_since_process_s, freq_bin);
        if self.want_phased {
            if let Some(v) = phased {
                self.phased[bin][freq_bin] += v;
            }
        }
        self.p11[bin][freq_bin] += p11;
        self.p12[bin][freq_bin] += p12;
        self.p22[bin][freq_bin] += p22;
        self.counts[bin][freq_bin] = self.counts[bin][freq_bin].saturating_add(1);
    }

    pub fn merge(&mut self, other: FoldAccum) {
        for b in 0..self.bins {
            for k in 0..self.nfreq {
                if self.want_phased {
                    self.phased[b][k] += other.phased[b][k];
                }
                self.p11[b][k] += other.p11[b][k];
                self.p12[b][k] += other.p12[b][k];
                self.p22[b][k] += other.p22[b][k];
                self.counts[b][k] = self.counts[b][k].saturating_add(other.counts[b][k]);
            }
        }
    }

    pub fn count_max(&self, bin: usize, start: usize, end: usize) -> u32 {
        self.counts[bin][start..end]
            .iter()
            .copied()
            .max()
            .unwrap_or(0)
    }

    pub fn normalized_real(
        &self,
        product: FoldProduct,
        bin: usize,
        start: usize,
        end: usize,
        norm_per_frame: f64,
    ) -> Vec<Complex<f32>> {
        let src = match product {
            FoldProduct::Phased => &self.phased[bin],
            FoldProduct::P11 => &self.p11[bin],
            FoldProduct::P22 => &self.p22[bin],
        };
        (start..end)
            .map(|k| {
                let c = self.counts[bin][k] as f64;
                let v = if c > 0.0 {
                    src[k] * norm_per_frame / c
                } else {
                    0.0
                };
                Complex::new(v as f32, 0.0)
            })
            .collect()
    }

    pub fn normalized_complex(
        &self,
        bin: usize,
        start: usize,
        end: usize,
        norm_per_frame: f64,
    ) -> Vec<Complex<f32>> {
        (start..end)
            .map(|k| {
                let c = self.counts[bin][k] as f64;
                if c > 0.0 {
                    let v = self.p12[bin][k] * (norm_per_frame / c);
                    Complex::new(v.re as f32, v.im as f32)
                } else {
                    Complex::new(0.0, 0.0)
                }
            })
            .collect()
    }
}

#[derive(Clone, Copy, Debug)]
pub enum FoldProduct {
    Phased,
    P11,
    P22,
}
