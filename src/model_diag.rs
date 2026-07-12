use crate::geom::{self, EarthOrientation, GeometricDelayMode, SourceVectorMode};

const DIAGNOSTIC_RATE_STEP_S: usize = 30;
const DIAGNOSTIC_ACCEL_STEP_S: usize = 60;
const DIAGNOSTIC_JERK_STEP_S: usize = 120;
const DIAGNOSTIC_SNAP_STEP_S: usize = 300;
const DIAGNOSTIC_GUARD_S: usize = 3 * DIAGNOSTIC_SNAP_STEP_S;

pub(crate) fn diagnostic_derivative_stencil_label() -> &'static str {
    "7-point central finite differences: rate_h=30s accel_h=60s jerk_h=120s snap_h=300s"
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct DelayDerivatives {
    pub delay_s: f64,
    pub rate_sps: f64,
    pub accel_sps2: f64,
    pub jerk_sps3: f64,
    pub snap_sps4: f64,
}

impl DelayDerivatives {
    pub(crate) fn from_epoch_polynomial(
        delay_s: f64,
        rate_sps: f64,
        accel_sps2: f64,
        jerk_sps3: f64,
        snap_sps4: f64,
        t_s: f64,
    ) -> Self {
        let t2 = t_s * t_s;
        let t3 = t2 * t_s;
        let t4 = t2 * t2;
        Self {
            delay_s: delay_s
                + rate_sps * t_s
                + 0.5 * accel_sps2 * t2
                + jerk_sps3 * t3 / 6.0
                + snap_sps4 * t4 / 24.0,
            rate_sps: rate_sps + accel_sps2 * t_s + 0.5 * jerk_sps3 * t2 + snap_sps4 * t3 / 6.0,
            accel_sps2: accel_sps2 + jerk_sps3 * t_s + 0.5 * snap_sps4 * t2,
            jerk_sps3: jerk_sps3 + snap_sps4 * t_s,
            snap_sps4,
        }
    }

    pub(crate) fn add(self, rhs: Self) -> Self {
        Self {
            delay_s: self.delay_s + rhs.delay_s,
            rate_sps: self.rate_sps + rhs.rate_sps,
            accel_sps2: self.accel_sps2 + rhs.accel_sps2,
            jerk_sps3: self.jerk_sps3 + rhs.jerk_sps3,
            snap_sps4: self.snap_sps4 + rhs.snap_sps4,
        }
    }

    pub(crate) fn sub(self, rhs: Self) -> Self {
        Self {
            delay_s: self.delay_s - rhs.delay_s,
            rate_sps: self.rate_sps - rhs.rate_sps,
            accel_sps2: self.accel_sps2 - rhs.accel_sps2,
            jerk_sps3: self.jerk_sps3 - rhs.jerk_sps3,
            snap_sps4: self.snap_sps4 - rhs.snap_sps4,
        }
    }

    pub(crate) fn phase_derivatives(self, frequency_hz: f64) -> Self {
        Self {
            delay_s: -frequency_hz * self.delay_s,
            rate_sps: -frequency_hz * self.rate_sps,
            accel_sps2: -frequency_hz * self.accel_sps2,
            jerk_sps3: -frequency_hz * self.jerk_sps3,
            snap_sps4: -frequency_hz * self.snap_sps4,
        }
    }

    pub(crate) fn phase_polynomial_deg(self, frequency_hz: f64) -> [f64; 5] {
        [
            -360.0 * frequency_hz * self.delay_s,
            -360.0 * frequency_hz * self.rate_sps,
            -180.0 * frequency_hz * self.accel_sps2,
            -60.0 * frequency_hz * self.jerk_sps3,
            -15.0 * frequency_hz * self.snap_sps4,
        ]
    }
}

#[derive(Clone, Copy)]
pub(crate) struct GeometryContext {
    pub ant1_ecef_m: [f64; 3],
    pub ant2_ecef_m: [f64; 3],
    pub ra_j2000_rad: f64,
    pub dec_j2000_rad: f64,
    pub reference_mjd_utc: f64,
}

#[derive(Clone, Copy)]
pub(crate) struct ModelVariant {
    pub source_mode: SourceVectorMode,
    pub delay_mode: GeometricDelayMode,
    pub eop: EarthOrientation,
    pub eop_zero: bool,
}

impl ModelVariant {
    pub(crate) fn key(self) -> String {
        if self.eop_zero {
            return format!(
                "{}__{}__eop-zero",
                source_mode_key(self.source_mode),
                delay_mode_key(self.delay_mode)
            );
        }
        format!(
            "{}__{}",
            source_mode_key(self.source_mode),
            delay_mode_key(self.delay_mode)
        )
    }
}

pub(crate) fn source_mode_key(mode: SourceVectorMode) -> &'static str {
    match mode {
        SourceVectorMode::MeanGast => "mean-gast",
        SourceVectorMode::PnmGast => "pnm-gast",
        SourceVectorMode::PnmEra => "pnm-era-engineering",
    }
}

pub(crate) fn delay_mode_key(mode: GeometricDelayMode) -> &'static str {
    match mode {
        GeometricDelayMode::Anchored => "anchored",
        GeometricDelayMode::Barycentric => "barycentric",
        GeometricDelayMode::VlbiMinus => "vlbi-minus",
        GeometricDelayMode::VlbiPlus => "vlbi-plus",
        GeometricDelayMode::Geocentric => "geocentric",
    }
}

pub(crate) fn all_model_variants(
    active_source_mode: SourceVectorMode,
    active_delay_mode: GeometricDelayMode,
    active_eop: EarthOrientation,
) -> Vec<ModelVariant> {
    let active = ModelVariant {
        source_mode: active_source_mode,
        delay_mode: active_delay_mode,
        eop: active_eop,
        eop_zero: false,
    };
    let mut out = vec![active];
    for source_mode in [
        SourceVectorMode::MeanGast,
        SourceVectorMode::PnmGast,
        SourceVectorMode::PnmEra,
    ] {
        for delay_mode in [
            GeometricDelayMode::Anchored,
            GeometricDelayMode::Barycentric,
            GeometricDelayMode::VlbiMinus,
            GeometricDelayMode::VlbiPlus,
            GeometricDelayMode::Geocentric,
        ] {
            if source_mode == active_source_mode && delay_mode == active_delay_mode {
                continue;
            }
            out.push(ModelVariant {
                source_mode,
                delay_mode,
                eop: active_eop,
                eop_zero: false,
            });
        }
    }
    out.push(ModelVariant {
        source_mode: active_source_mode,
        delay_mode: active_delay_mode,
        eop: EarthOrientation {
            dut1_s: 0.0,
            tt_minus_utc_s: active_eop.tt_minus_utc_s,
            xp_arcsec: 0.0,
            yp_arcsec: 0.0,
        },
        eop_zero: true,
    });
    out
}

fn geometric_delay_at(ctx: GeometryContext, variant: ModelVariant, elapsed_process_s: f64) -> f64 {
    let mjd_utc = ctx.reference_mjd_utc + elapsed_process_s / 86400.0;
    let (ra, dec) = match variant.source_mode {
        SourceVectorMode::MeanGast => geom::precess_j2000_to_mean_of_date(
            ctx.ra_j2000_rad,
            ctx.dec_j2000_rad,
            mjd_utc + variant.eop.tt_minus_utc_s / 86400.0,
        ),
        SourceVectorMode::PnmGast | SourceVectorMode::PnmEra => {
            (ctx.ra_j2000_rad, ctx.dec_j2000_rad)
        }
    };
    geom::calculate_geometric_delay_and_derivatives_full_with_eop(
        ctx.ant1_ecef_m,
        ctx.ant2_ecef_m,
        ra,
        dec,
        mjd_utc,
        ctx.reference_mjd_utc,
        variant.eop,
        variant.delay_mode,
        variant.source_mode,
    )
    .2
}

fn solve_5x5(mut a: [[f64; 5]; 5], mut b: [f64; 5]) -> Option<[f64; 5]> {
    for col in 0..5 {
        let mut pivot = col;
        let mut pivot_abs = a[col][col].abs();
        for (row, values) in a.iter().enumerate().skip(col + 1) {
            let value = values[col].abs();
            if value > pivot_abs {
                pivot = row;
                pivot_abs = value;
            }
        }
        if pivot_abs < 1.0e-24 {
            return None;
        }
        if pivot != col {
            a.swap(col, pivot);
            b.swap(col, pivot);
        }
        let inv = 1.0 / a[col][col];
        for value in &mut a[col][col..] {
            *value *= inv;
        }
        b[col] *= inv;
        for row in 0..5 {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            if factor == 0.0 {
                continue;
            }
            for j in col..5 {
                a[row][j] -= factor * a[col][j];
            }
            b[row] -= factor * b[col];
        }
    }
    Some(b)
}

pub(crate) fn local_quartic_derivatives(
    delay_grid: &[f64],
    center_idx: usize,
    radius: usize,
) -> (f64, f64, f64, f64) {
    let center = delay_grid[center_idx];
    let start = center_idx.saturating_sub(radius);
    let end = (center_idx + radius).min(delay_grid.len().saturating_sub(1));
    let mut normal = [[0.0_f64; 5]; 5];
    let mut rhs = [0.0_f64; 5];
    for (idx, value) in delay_grid.iter().enumerate().take(end + 1).skip(start) {
        let x = idx as f64 - center_idx as f64;
        let y = *value - center;
        let mut powers = [1.0_f64; 9];
        for k in 1..powers.len() {
            powers[k] = powers[k - 1] * x;
        }
        for row in 0..5 {
            rhs[row] += powers[row] * y;
            for col in 0..5 {
                normal[row][col] += powers[row + col];
            }
        }
    }
    if let Some(coefficients) = solve_5x5(normal, rhs) {
        (
            coefficients[1],
            2.0 * coefficients[2],
            6.0 * coefficients[3],
            24.0 * coefficients[4],
        )
    } else {
        (0.0, 0.0, 0.0, 0.0)
    }
}

fn stable_diagnostic_derivatives(
    center_delay_s: f64,
    mut relative_sample: impl FnMut(isize) -> f64,
) -> DelayDerivatives {
    let centered = |value: f64| value - center_delay_s;

    let h1 = DIAGNOSTIC_RATE_STEP_S as isize;
    let fm3 = centered(relative_sample(-3 * h1));
    let fm2 = centered(relative_sample(-2 * h1));
    let fm1 = centered(relative_sample(-h1));
    let fp1 = centered(relative_sample(h1));
    let fp2 = centered(relative_sample(2 * h1));
    let fp3 = centered(relative_sample(3 * h1));
    let rate_sps =
        (-fm3 + 9.0 * fm2 - 45.0 * fm1 + 45.0 * fp1 - 9.0 * fp2 + fp3) / (60.0 * h1 as f64);

    let h2 = DIAGNOSTIC_ACCEL_STEP_S as isize;
    let fm3 = centered(relative_sample(-3 * h2));
    let fm2 = centered(relative_sample(-2 * h2));
    let fm1 = centered(relative_sample(-h2));
    let fp1 = centered(relative_sample(h2));
    let fp2 = centered(relative_sample(2 * h2));
    let fp3 = centered(relative_sample(3 * h2));
    let accel_sps2 = (2.0 * fm3 - 27.0 * fm2 + 270.0 * fm1 + 270.0 * fp1 - 27.0 * fp2 + 2.0 * fp3)
        / (180.0 * (h2 as f64).powi(2));

    let h3 = DIAGNOSTIC_JERK_STEP_S as isize;
    let fm3 = centered(relative_sample(-3 * h3));
    let fm2 = centered(relative_sample(-2 * h3));
    let fm1 = centered(relative_sample(-h3));
    let fp1 = centered(relative_sample(h3));
    let fp2 = centered(relative_sample(2 * h3));
    let fp3 = centered(relative_sample(3 * h3));
    let jerk_sps3 =
        (fm3 - 8.0 * fm2 + 13.0 * fm1 - 13.0 * fp1 + 8.0 * fp2 - fp3) / (8.0 * (h3 as f64).powi(3));

    let h4 = DIAGNOSTIC_SNAP_STEP_S as isize;
    let fm3 = centered(relative_sample(-3 * h4));
    let fm2 = centered(relative_sample(-2 * h4));
    let fm1 = centered(relative_sample(-h4));
    let fp1 = centered(relative_sample(h4));
    let fp2 = centered(relative_sample(2 * h4));
    let fp3 = centered(relative_sample(3 * h4));
    let snap_sps4 = (-fm3 + 12.0 * fm2 - 39.0 * fm1 - 39.0 * fp1 + 12.0 * fp2 - fp3)
        / (6.0 * (h4 as f64).powi(4));

    DelayDerivatives {
        delay_s: center_delay_s,
        rate_sps,
        accel_sps2,
        jerk_sps3,
        snap_sps4,
    }
}

pub(crate) fn evaluate_at(
    ctx: GeometryContext,
    variant: ModelVariant,
    elapsed_process_s: f64,
) -> DelayDerivatives {
    let center = geometric_delay_at(ctx, variant, elapsed_process_s);
    stable_diagnostic_derivatives(center, |offset| {
        geometric_delay_at(ctx, variant, elapsed_process_s + offset as f64)
    })
}

pub(crate) fn evaluate_integer_series(
    ctx: GeometryContext,
    variant: ModelVariant,
    start_elapsed_process_s: f64,
    rows: usize,
) -> Vec<DelayDerivatives> {
    let guard = DIAGNOSTIC_GUARD_S as isize;
    let end = rows.saturating_sub(1) as isize + guard;
    let delay_grid: Vec<f64> = (-guard..=end)
        .map(|offset| geometric_delay_at(ctx, variant, start_elapsed_process_s + offset as f64))
        .collect();
    (0..rows)
        .map(|row| {
            let center_idx = row + DIAGNOSTIC_GUARD_S;
            let center = delay_grid[center_idx];
            stable_diagnostic_derivatives(center, |offset| {
                delay_grid[(center_idx as isize + offset) as usize]
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{local_quartic_derivatives, stable_diagnostic_derivatives, DelayDerivatives};

    #[test]
    fn quartic_fit_recovers_delay_derivatives() {
        let delay = 2.5;
        let rate = -3.0e-4;
        let accel = 4.0e-7;
        let jerk = -5.0e-10;
        let snap = 6.0e-13;
        let grid: Vec<f64> = (-10..=10)
            .map(|t| {
                DelayDerivatives::from_epoch_polynomial(delay, rate, accel, jerk, snap, t as f64)
                    .delay_s
            })
            .collect();
        let (fit_rate, fit_accel, fit_jerk, fit_snap) = local_quartic_derivatives(&grid, 10, 10);
        assert!((fit_rate - rate).abs() < 1.0e-14);
        assert!((fit_accel - accel).abs() < 1.0e-14);
        assert!((fit_jerk - jerk).abs() < 1.0e-14);
        assert!((fit_snap - snap).abs() < 1.0e-14);
    }

    #[test]
    fn phase_polynomial_uses_factorials_and_correction_sign() {
        let d = DelayDerivatives {
            delay_s: 1.0,
            rate_sps: 2.0,
            accel_sps2: 3.0,
            jerk_sps3: 4.0,
            snap_sps4: 5.0,
        };
        let c = d.phase_polynomial_deg(10.0);
        assert_eq!(c, [-3600.0, -7200.0, -5400.0, -2400.0, -750.0]);
    }

    #[test]
    fn wide_diagnostic_stencil_recovers_quartic_derivatives() {
        let delay = 2.5;
        let rate = -3.0e-4;
        let accel = 4.0e-7;
        let jerk = -5.0e-10;
        let snap = 6.0e-13;
        let value = |t: isize| {
            DelayDerivatives::from_epoch_polynomial(delay, rate, accel, jerk, snap, t as f64)
                .delay_s
        };
        let fit = stable_diagnostic_derivatives(value(0), value);
        assert!((fit.delay_s - delay).abs() < 1.0e-15);
        assert!((fit.rate_sps - rate).abs() < 1.0e-14);
        assert!((fit.accel_sps2 - accel).abs() < 1.0e-14);
        assert!((fit.jerk_sps3 - jerk).abs() < 1.0e-14);
        assert!((fit.snap_sps4 - snap).abs() < 1.0e-14);
    }
}
