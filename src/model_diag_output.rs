use crate::geom::{self, EarthOrientation, GeometricDelayMode, SourceVectorMode};
use crate::model_diag::{
    all_model_variants, diagnostic_derivative_stencil_label, evaluate_at, evaluate_integer_series,
    DelayDerivatives, GeometryContext, ModelVariant,
};
use crate::utils::DynError;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Clone)]
pub(crate) struct DiagnosticInput {
    pub ant1_name: String,
    pub ant2_name: String,
    pub epoch: String,
    pub geometry: GeometryContext,
    pub active_source_mode: SourceVectorMode,
    pub active_delay_mode: GeometricDelayMode,
    pub active_eop: EarthOrientation,
    pub start_elapsed_process_s: f64,
    pub total_skip_s: f64,
    pub model_time_offset_s: f64,
    pub duration_s: f64,
    pub frame_dt_s: f64,
    pub total_frames: usize,
    pub fs_hz: f64,
    pub obs_hz: f64,
    pub carrier_hz: f64,
    pub clock_epoch_coefficients: DelayDerivatives,
    pub other_epoch_coefficients: DelayDerivatives,
    pub read_align_reference_s: f64,
    pub d_seek_s: f64,
    pub residual_target: String,
}

#[derive(Clone)]
pub(crate) struct VariantCheckpoint {
    pub key: String,
    pub geom: DelayDerivatives,
    pub delta: DelayDerivatives,
    pub delta_phase_from_start_deg: f64,
    pub delta_phase_after_linear_deg: f64,
    pub delta_phase_after_quadratic_deg: f64,
}

#[derive(Clone)]
pub(crate) struct CheckpointReport {
    pub label: &'static str,
    pub frame_idx: usize,
    pub scan_elapsed_s: f64,
    pub data_elapsed_process_s: f64,
    pub data_mjd_utc: f64,
    pub model_elapsed_process_s: f64,
    pub model_mjd_utc: f64,
    pub geom: DelayDerivatives,
    pub clock: DelayDerivatives,
    pub other: DelayDerivatives,
    pub total: DelayDerivatives,
    pub ant1_az_deg: f64,
    pub ant1_el_deg: f64,
    pub ant1_csc_el: f64,
    pub ant1_csc_el_valid: bool,
    pub ant2_az_deg: f64,
    pub ant2_el_deg: f64,
    pub ant2_csc_el: f64,
    pub ant2_csc_el_valid: bool,
    pub u_lambda: f64,
    pub v_lambda: f64,
    pub w_lambda: f64,
    pub variants: Vec<VariantCheckpoint>,
}

pub(crate) struct DiagnosticReport {
    pub path: PathBuf,
    pub rows: usize,
    pub variants: usize,
    pub checkpoints: Vec<CheckpointReport>,
}

fn require_finite(value: f64, field: &str, model: &str, row: &str) -> Result<(), DynError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(format!(
            "non-finite diagnostic value: field={field} model={model} row={row} value={value:?}"
        )
        .into())
    }
}

fn validate_derivatives(
    value: DelayDerivatives,
    field: &str,
    model: &str,
    row: &str,
) -> Result<(), DynError> {
    require_finite(value.delay_s, &format!("{field}.delay_s"), model, row)?;
    require_finite(value.rate_sps, &format!("{field}.rate_sps"), model, row)?;
    require_finite(value.accel_sps2, &format!("{field}.accel_sps2"), model, row)?;
    require_finite(value.jerk_sps3, &format!("{field}.jerk_sps3"), model, row)?;
    require_finite(value.snap_sps4, &format!("{field}.snap_sps4"), model, row)?;
    Ok(())
}

fn source_coordinates(
    geometry: GeometryContext,
    variant: ModelVariant,
    elapsed_process_s: f64,
) -> (f64, f64, f64) {
    let mjd_utc = geometry.reference_mjd_utc + elapsed_process_s / 86400.0;
    let (ra, dec) = match variant.source_mode {
        SourceVectorMode::MeanGast => geom::precess_j2000_to_mean_of_date(
            geometry.ra_j2000_rad,
            geometry.dec_j2000_rad,
            mjd_utc + variant.eop.tt_minus_utc_s / 86400.0,
        ),
        SourceVectorMode::PnmGast | SourceVectorMode::PnmEra => {
            (geometry.ra_j2000_rad, geometry.dec_j2000_rad)
        }
    };
    (ra, dec, mjd_utc)
}

#[derive(Clone, Copy)]
struct AzElBasisValues {
    az1_deg: f64,
    el1_deg: f64,
    csc1: f64,
    csc1_valid: bool,
    az2_deg: f64,
    el2_deg: f64,
    csc2: f64,
    csc2_valid: bool,
    u_lambda: f64,
    v_lambda: f64,
    w_lambda: f64,
}

fn csc_with_validity(elevation_rad: f64) -> (f64, bool) {
    let sin_el = elevation_rad.sin();
    if sin_el > 1.0e-12 {
        (1.0 / sin_el, true)
    } else {
        // Zero is a finite placeholder. Consumers must check the adjacent
        // validity column before interpreting it as a mapping factor.
        (0.0, false)
    }
}

fn az_el_and_basis(
    input: &DiagnosticInput,
    active: ModelVariant,
    elapsed_process_s: f64,
) -> AzElBasisValues {
    let (ra, dec, mjd_utc) = source_coordinates(input.geometry, active, elapsed_process_s);
    let (az1, el1) = geom::source_az_el_with_eop_mode(
        input.geometry.ant1_ecef_m,
        ra,
        dec,
        mjd_utc,
        active.eop,
        active.source_mode,
    );
    let (az2, el2) = geom::source_az_el_with_eop_mode(
        input.geometry.ant2_ecef_m,
        ra,
        dec,
        mjd_utc,
        active.eop,
        active.source_mode,
    );
    let basis = geom::baseline_phase_basis(
        input.geometry.ant1_ecef_m,
        input.geometry.ant2_ecef_m,
        ra,
        dec,
        mjd_utc,
        active.eop,
        active.source_mode,
        input.obs_hz,
    );
    let (csc1, csc1_valid) = csc_with_validity(el1);
    let (csc2, csc2_valid) = csc_with_validity(el2);
    AzElBasisValues {
        az1_deg: az1.to_degrees(),
        el1_deg: el1.to_degrees(),
        csc1,
        csc1_valid,
        az2_deg: az2.to_degrees(),
        el2_deg: el2.to_degrees(),
        csc2,
        csc2_valid,
        u_lambda: basis.u_lambda,
        v_lambda: basis.v_lambda,
        w_lambda: basis.w_lambda,
    }
}

fn validate_az_el_basis(value: AzElBasisValues, model: &str, row: &str) -> Result<(), DynError> {
    for (field, number) in [
        ("ant1_az_deg", value.az1_deg),
        ("ant1_el_deg", value.el1_deg),
        ("ant1_csc_el", value.csc1),
        ("ant2_az_deg", value.az2_deg),
        ("ant2_el_deg", value.el2_deg),
        ("ant2_csc_el", value.csc2),
        ("u_lambda", value.u_lambda),
        ("v_lambda", value.v_lambda),
        ("w_lambda", value.w_lambda),
    ] {
        require_finite(number, field, model, row)?;
    }
    Ok(())
}

fn sampling_times(duration_s: f64) -> Result<Vec<f64>, DynError> {
    require_finite(duration_s, "duration_s", "input", "metadata")?;
    if duration_s < 0.0 {
        return Err(format!("invalid diagnostic duration_s={duration_s}: expected >= 0").into());
    }
    let whole = duration_s.floor();
    if whole > (usize::MAX - 1) as f64 {
        return Err(format!("diagnostic duration_s={duration_s} is too large").into());
    }
    let whole_seconds = whole as usize;
    let mut times = Vec::with_capacity(whole_seconds.saturating_add(2));
    times.extend((0..=whole_seconds).map(|second| second as f64));
    if duration_s - whole > 1.0e-9 {
        times.push(duration_s);
    }
    Ok(times)
}

pub(crate) fn diagnostic_checkpoint_indices(total_frames: usize) -> Vec<usize> {
    if total_frames == 0 {
        return Vec::new();
    }
    let last = total_frames - 1;
    let mut indices = vec![0, last / 4, last / 2, 3 * last / 4, last];
    indices.sort_unstable();
    indices.dedup();
    indices
}

fn checkpoint_label(position: usize, count: usize) -> &'static str {
    match (position, count) {
        (0, _) => "start",
        (p, c) if p + 1 == c => "end",
        (1, 5) => "quarter",
        (2, 5) => "mid",
        (3, 5) => "three-quarter",
        _ => "intermediate",
    }
}

fn push_derivative_columns(
    columns: &mut Vec<String>,
    value: DelayDerivatives,
    fs_hz: f64,
    field: &str,
    model: &str,
    row: &str,
) -> Result<(), DynError> {
    let scaled = [
        ("delay_sample", value.delay_s * fs_hz),
        ("rate_sample_s", value.rate_sps * fs_hz),
        ("accel_sample_s2", value.accel_sps2 * fs_hz),
        ("jerk_sample_s3", value.jerk_sps3 * fs_hz),
        ("snap_sample_s4", value.snap_sps4 * fs_hz),
    ];
    for (suffix, number) in scaled {
        require_finite(number, &format!("{field}.{suffix}"), model, row)?;
        columns.push(format!("{number:+.17e}"));
    }
    Ok(())
}

fn tsv_header_columns(variants: &[ModelVariant]) -> Vec<String> {
    let mut columns = [
        "scan_elapsed_s",
        "data_elapsed_process_s",
        "data_mjd_utc",
        "model_elapsed_process_s",
        "model_mjd_utc",
        "ant1_az_deg",
        "ant1_el_deg",
        "ant1_csc_el",
        "ant1_csc_el_valid",
        "ant2_az_deg",
        "ant2_el_deg",
        "ant2_csc_el",
        "ant2_csc_el_valid",
        "u_lambda",
        "v_lambda",
        "w_lambda",
    ]
    .into_iter()
    .map(str::to_string)
    .collect::<Vec<_>>();
    for component in ["active_geom", "clock", "other", "active_total"] {
        for suffix in [
            "delay_sample",
            "rate_sample_s",
            "accel_sample_s2",
            "jerk_sample_s3",
            "snap_sample_s4",
        ] {
            columns.push(format!("{component}_{suffix}"));
        }
    }
    columns.extend([
        "active_phase_cycles".to_string(),
        "active_phase_rate_hz".to_string(),
        "active_phase_accel_hz_s".to_string(),
        "active_phase_jerk_hz_s2".to_string(),
        "active_phase_snap_hz_s3".to_string(),
    ]);
    for variant in variants {
        let key = variant.key();
        for suffix in [
            "geom_delay_sample",
            "geom_rate_sample_s",
            "geom_accel_sample_s2",
            "geom_jerk_sample_s3",
            "geom_snap_sample_s4",
            "delta_delay_sample",
            "delta_phase_rate_hz",
            "delta_phase_accel_hz_s",
            "delta_phase_jerk_hz_s2",
            "delta_phase_snap_hz_s3",
            "delta_phase_from_start_deg",
            "delta_phase_after_linear_deg",
            "delta_phase_after_quadratic_deg",
        ] {
            columns.push(format!("{key}_{suffix}"));
        }
    }
    columns
}

fn validate_input(input: &DiagnosticInput) -> Result<(), DynError> {
    let metadata = "metadata";
    for (index, value) in input.geometry.ant1_ecef_m.iter().copied().enumerate() {
        require_finite(value, &format!("ant1_ecef_m[{index}]"), "input", metadata)?;
    }
    for (index, value) in input.geometry.ant2_ecef_m.iter().copied().enumerate() {
        require_finite(value, &format!("ant2_ecef_m[{index}]"), "input", metadata)?;
    }
    for (field, value) in [
        ("ra_j2000_rad", input.geometry.ra_j2000_rad),
        ("dec_j2000_rad", input.geometry.dec_j2000_rad),
        ("reference_mjd_utc", input.geometry.reference_mjd_utc),
        ("start_elapsed_process_s", input.start_elapsed_process_s),
        ("total_skip_s", input.total_skip_s),
        ("model_time_offset_s", input.model_time_offset_s),
        ("duration_s", input.duration_s),
        ("frame_dt_s", input.frame_dt_s),
        ("fs_hz", input.fs_hz),
        ("obs_hz", input.obs_hz),
        ("carrier_hz", input.carrier_hz),
        ("read_align_reference_s", input.read_align_reference_s),
        ("d_seek_s", input.d_seek_s),
        ("eop.dut1_s", input.active_eop.dut1_s),
        ("eop.tt_minus_utc_s", input.active_eop.tt_minus_utc_s),
        ("eop.xp_arcsec", input.active_eop.xp_arcsec),
        ("eop.yp_arcsec", input.active_eop.yp_arcsec),
    ] {
        require_finite(value, field, "input", metadata)?;
    }
    for (field, value) in [
        ("frame_dt_s", input.frame_dt_s),
        ("fs_hz", input.fs_hz),
        ("obs_hz", input.obs_hz),
    ] {
        if value <= 0.0 {
            return Err(format!(
                "invalid diagnostic value: field={field} model=input row={metadata} value={value:?}; expected > 0"
            )
            .into());
        }
    }
    if input.duration_s < 0.0 {
        return Err(format!(
            "invalid diagnostic value: field=duration_s model=input row={metadata} value={:?}; expected >= 0",
            input.duration_s
        )
        .into());
    }
    validate_derivatives(
        input.clock_epoch_coefficients,
        "clock_epoch_coefficients",
        "input",
        metadata,
    )?;
    validate_derivatives(
        input.other_epoch_coefficients,
        "other_epoch_coefficients",
        "input",
        metadata,
    )?;
    Ok(())
}

fn evaluate_sampling_series(
    geometry: GeometryContext,
    variant: ModelVariant,
    start_elapsed_process_s: f64,
    times: &[f64],
) -> Vec<DelayDerivatives> {
    let integer_rows = times.iter().take_while(|time| time.fract() == 0.0).count();
    let mut values =
        evaluate_integer_series(geometry, variant, start_elapsed_process_s, integer_rows);
    values.extend(
        times[integer_rows..]
            .iter()
            .map(|time| evaluate_at(geometry, variant, start_elapsed_process_s + *time)),
    );
    values
}

#[derive(Clone, Copy)]
struct PolynomialFit {
    t0: f64,
    scale: f64,
    coefficients: [f64; 3],
}

impl PolynomialFit {
    fn evaluate(self, time_s: f64) -> f64 {
        let x = (time_s - self.t0) / self.scale;
        self.coefficients[0] + self.coefficients[1] * x + self.coefficients[2] * x * x
    }
}

#[derive(Clone, Copy)]
struct VariantDetrend {
    linear: PolynomialFit,
    quadratic: PolynomialFit,
}

fn polynomial_fit(
    times: &[f64],
    values: &[f64],
    requested_degree: usize,
    model: &str,
) -> Result<PolynomialFit, DynError> {
    if times.is_empty() || times.len() != values.len() {
        return Err(format!(
            "invalid detrend input: field=delta_delay model={model} row=all times={} values={}",
            times.len(),
            values.len()
        )
        .into());
    }
    let degree = requested_degree.min(times.len() - 1).min(2);
    let n = degree + 1;
    let t0 = times[0];
    let scale = (times[times.len() - 1] - t0).abs().max(1.0);
    let mut normal = [[0.0_f64; 3]; 3];
    let mut rhs = [0.0_f64; 3];
    for (row, (&time, &value)) in times.iter().zip(values).enumerate() {
        require_finite(time, "detrend_time_s", model, &row.to_string())?;
        require_finite(value, "delta_delay_s", model, &row.to_string())?;
        let x = (time - t0) / scale;
        let powers = [1.0, x, x * x, x * x * x, x * x * x * x];
        for i in 0..n {
            rhs[i] += powers[i] * value;
            for j in 0..n {
                normal[i][j] += powers[i + j];
            }
        }
    }
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if normal[row][col].abs() > normal[pivot][col].abs() {
                pivot = row;
            }
        }
        if normal[pivot][col].abs() < 1.0e-24 {
            return Err(format!(
                "singular detrend fit: field=delta_delay model={model} row=all degree={degree}"
            )
            .into());
        }
        if pivot != col {
            normal.swap(pivot, col);
            rhs.swap(pivot, col);
        }
        let inverse = 1.0 / normal[col][col];
        for value in &mut normal[col][col..n] {
            *value *= inverse;
        }
        rhs[col] *= inverse;
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = normal[row][col];
            for j in col..n {
                normal[row][j] -= factor * normal[col][j];
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    for (index, coefficient) in rhs.iter().copied().enumerate().take(n) {
        require_finite(
            coefficient,
            &format!("detrend_coefficient[{index}]"),
            model,
            "all",
        )?;
    }
    Ok(PolynomialFit {
        t0,
        scale,
        coefficients: rhs,
    })
}

pub(crate) fn generate_model_diagnostics(
    input: &DiagnosticInput,
    path: &Path,
) -> Result<DiagnosticReport, DynError> {
    validate_input(input)?;
    let sample_times = sampling_times(input.duration_s)?;
    let rows = sample_times.len();
    let variants = all_model_variants(
        input.active_source_mode,
        input.active_delay_mode,
        input.active_eop,
    );
    let series: Vec<Vec<DelayDerivatives>> = variants
        .iter()
        .map(|variant| {
            evaluate_sampling_series(
                input.geometry,
                *variant,
                input.start_elapsed_process_s,
                &sample_times,
            )
        })
        .collect();
    let active = variants[0];
    let active_series = &series[0];
    for (variant, values) in variants.iter().zip(&series) {
        let model = variant.key();
        for (row, value) in values.iter().copied().enumerate() {
            validate_derivatives(value, "geom", &model, &row.to_string())?;
        }
    }
    let detrends = variants
        .iter()
        .enumerate()
        .map(|(variant_index, variant)| {
            let delta_delays = series[variant_index]
                .iter()
                .zip(active_series)
                .map(|(value, active_value)| value.sub(*active_value).delay_s)
                .collect::<Vec<_>>();
            let model = variant.key();
            Ok(VariantDetrend {
                linear: polynomial_fit(&sample_times, &delta_delays, 1, &model)?,
                quadratic: polynomial_fit(&sample_times, &delta_delays, 2, &model)?,
            })
        })
        .collect::<Result<Vec<_>, DynError>>()?;

    let data_start_mjd_utc = input.geometry.reference_mjd_utc + input.total_skip_s / 86400.0;
    let model_start_mjd_utc =
        input.geometry.reference_mjd_utc + input.start_elapsed_process_s / 86400.0;
    require_finite(
        data_start_mjd_utc,
        "data_start_mjd_utc",
        "input",
        "metadata",
    )?;
    require_finite(
        model_start_mjd_utc,
        "model_start_mjd_utc",
        "input",
        "metadata",
    )?;

    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "# yi-corr one-pass geometric model diagnostics")?;
    writeln!(writer, "# epoch={}", input.epoch)?;
    writeln!(
        writer,
        "# antenna1={} ecef_m={:+.9},{:+.9},{:+.9}",
        input.ant1_name,
        input.geometry.ant1_ecef_m[0],
        input.geometry.ant1_ecef_m[1],
        input.geometry.ant1_ecef_m[2]
    )?;
    writeln!(
        writer,
        "# antenna2={} ecef_m={:+.9},{:+.9},{:+.9}",
        input.ant2_name,
        input.geometry.ant2_ecef_m[0],
        input.geometry.ant2_ecef_m[1],
        input.geometry.ant2_ecef_m[2]
    )?;
    writeln!(
        writer,
        "# reference_mjd_utc={:.12} ra_j2000_rad={:.17e} dec_j2000_rad={:.17e}",
        input.geometry.reference_mjd_utc, input.geometry.ra_j2000_rad, input.geometry.dec_j2000_rad
    )?;
    writeln!(
        writer,
        "# total_skip_s={:+.17e} model_time_offset_s={:+.17e} duration_s={:.9} frame_dt_s={:.17e} frames={}",
        input.total_skip_s,
        input.model_time_offset_s,
        input.duration_s,
        input.frame_dt_s,
        input.total_frames
    )?;
    writeln!(
        writer,
        "# data_time_axis=start_elapsed_process_s={:+.17e} start_mjd_utc={:.12}; scan_elapsed_s=0 is the leading data timestamp",
        input.total_skip_s, data_start_mjd_utc
    )?;
    writeln!(
        writer,
        "# model_time_axis=start_elapsed_process_s={:+.17e} start_mjd_utc={:.12}; includes model_time_offset_s={:+.17e}",
        input.start_elapsed_process_s, model_start_mjd_utc, input.model_time_offset_s
    )?;
    writeln!(
        writer,
        "# time_origin=scan_elapsed_s is measured from the first data timestamp; FFT frame midpoint is only a delay-model evaluation point and is not the phase origin"
    )?;
    writeln!(
        writer,
        "# phase_diagnostics=delta_phase_from_start_deg is the primary exact model difference; after_linear and after_quadratic remove full-scan least-squares delay polynomials"
    )?;
    writeln!(
        writer,
        "# diagnostic_derivatives={}",
        diagnostic_derivative_stencil_label()
    )?;
    writeln!(
        writer,
        "# fs_hz={:.9} obs_hz={:.9} carrier_hz={:.9}",
        input.fs_hz, input.obs_hz, input.carrier_hz
    )?;
    writeln!(
        writer,
        "# active_source={} active_delay={} eop_dut1_s={:+.17e} eop_tt_utc_s={:+.17e} eop_xp_arcsec={:+.17e} eop_yp_arcsec={:+.17e}",
        active.source_mode.label(),
        active.delay_mode.label(),
        active.eop.dut1_s,
        active.eop.tt_minus_utc_s,
        active.eop.xp_arcsec,
        active.eop.yp_arcsec
    )?;
    writeln!(
        writer,
        "# phase_convention=correction_cycles=-frequency_hz*delay_s; rate=-f*tau1; accel=-f*tau2; jerk=-f*tau3"
    )?;
    writeln!(
        writer,
        "# pnm-era-engineering is a sign/control path, not the standard CIO+ERA transform"
    )?;
    writeln!(
        writer,
        "# not_modeled=troposphere,ionosphere,solid-earth-tide,ocean-loading,pole-tide,station-velocity,axis-offset,gravitational-delay"
    )?;
    writeln!(
        writer,
        "# read_align_reference_s={:+.17e} d_seek_s={:+.17e} residual_target={}",
        input.read_align_reference_s, input.d_seek_s, input.residual_target
    )?;

    let header_columns = tsv_header_columns(&variants);
    let expected_columns = header_columns.len();
    writeln!(writer, "{}", header_columns.join("\t"))?;

    for (row, &scan_elapsed_s) in sample_times.iter().enumerate() {
        let row_label = format!("{row}@{scan_elapsed_s:.9}s");
        let data_elapsed_process_s = input.total_skip_s + scan_elapsed_s;
        let data_mjd_utc = input.geometry.reference_mjd_utc + data_elapsed_process_s / 86400.0;
        let model_elapsed_process_s = input.start_elapsed_process_s + scan_elapsed_s;
        let model_mjd_utc = input.geometry.reference_mjd_utc + model_elapsed_process_s / 86400.0;
        for (field, number) in [
            ("scan_elapsed_s", scan_elapsed_s),
            ("data_elapsed_process_s", data_elapsed_process_s),
            ("data_mjd_utc", data_mjd_utc),
            ("model_elapsed_process_s", model_elapsed_process_s),
            ("model_mjd_utc", model_mjd_utc),
        ] {
            require_finite(number, field, "active", &row_label)?;
        }
        let geom_active = active_series[row];
        let clock = DelayDerivatives::from_epoch_polynomial(
            input.clock_epoch_coefficients.delay_s,
            input.clock_epoch_coefficients.rate_sps,
            input.clock_epoch_coefficients.accel_sps2,
            input.clock_epoch_coefficients.jerk_sps3,
            input.clock_epoch_coefficients.snap_sps4,
            model_elapsed_process_s,
        );
        let other = DelayDerivatives::from_epoch_polynomial(
            input.other_epoch_coefficients.delay_s,
            input.other_epoch_coefficients.rate_sps,
            input.other_epoch_coefficients.accel_sps2,
            input.other_epoch_coefficients.jerk_sps3,
            input.other_epoch_coefficients.snap_sps4,
            model_elapsed_process_s,
        );
        let total = geom_active.add(clock).add(other);
        let phase = total.phase_derivatives(input.obs_hz);
        validate_derivatives(clock, "clock", "active", &row_label)?;
        validate_derivatives(other, "other", "active", &row_label)?;
        validate_derivatives(total, "total", "active", &row_label)?;
        validate_derivatives(phase, "phase", "active", &row_label)?;
        let sky = az_el_and_basis(input, active, model_elapsed_process_s);
        validate_az_el_basis(sky, "active", &row_label)?;
        let mut columns = vec![
            format!("{scan_elapsed_s:.6}"),
            format!("{data_elapsed_process_s:.9}"),
            format!("{data_mjd_utc:.12}"),
            format!("{model_elapsed_process_s:.9}"),
            format!("{model_mjd_utc:.12}"),
            format!("{:.9}", sky.az1_deg),
            format!("{:.9}", sky.el1_deg),
            format!("{:+.17e}", sky.csc1),
            u8::from(sky.csc1_valid).to_string(),
            format!("{:.9}", sky.az2_deg),
            format!("{:.9}", sky.el2_deg),
            format!("{:+.17e}", sky.csc2),
            u8::from(sky.csc2_valid).to_string(),
            format!("{:+.17e}", sky.u_lambda),
            format!("{:+.17e}", sky.v_lambda),
            format!("{:+.17e}", sky.w_lambda),
        ];
        push_derivative_columns(
            &mut columns,
            geom_active,
            input.fs_hz,
            "active_geom",
            "active",
            &row_label,
        )?;
        push_derivative_columns(
            &mut columns,
            clock,
            input.fs_hz,
            "clock",
            "active",
            &row_label,
        )?;
        push_derivative_columns(
            &mut columns,
            other,
            input.fs_hz,
            "other",
            "active",
            &row_label,
        )?;
        push_derivative_columns(
            &mut columns,
            total,
            input.fs_hz,
            "active_total",
            "active",
            &row_label,
        )?;
        for number in [
            phase.delay_s,
            phase.rate_sps,
            phase.accel_sps2,
            phase.jerk_sps3,
            phase.snap_sps4,
        ] {
            columns.push(format!("{number:+.17e}"));
        }
        for (variant_idx, values) in series.iter().enumerate() {
            let model = variants[variant_idx].key();
            let value = values[row];
            let delta = value.sub(geom_active);
            let delta0 = values[0].sub(active_series[0]);
            let phase_from_start_deg = -360.0 * input.obs_hz * (delta.delay_s - delta0.delay_s);
            let phase_after_linear_deg = -360.0
                * input.obs_hz
                * (delta.delay_s - detrends[variant_idx].linear.evaluate(scan_elapsed_s));
            let phase_after_quadratic_deg = -360.0
                * input.obs_hz
                * (delta.delay_s - detrends[variant_idx].quadratic.evaluate(scan_elapsed_s));
            let delta_phase = delta.phase_derivatives(input.obs_hz);
            validate_derivatives(delta, "delta", &model, &row_label)?;
            validate_derivatives(delta_phase, "delta_phase", &model, &row_label)?;
            for (field, number) in [
                ("delta_phase_from_start_deg", phase_from_start_deg),
                ("delta_phase_after_linear_deg", phase_after_linear_deg),
                ("delta_phase_after_quadratic_deg", phase_after_quadratic_deg),
            ] {
                require_finite(number, field, &model, &row_label)?;
            }
            push_derivative_columns(
                &mut columns,
                value,
                input.fs_hz,
                "variant_geom",
                &model,
                &row_label,
            )?;
            let delta_delay_sample = delta.delay_s * input.fs_hz;
            require_finite(delta_delay_sample, "delta_delay_sample", &model, &row_label)?;
            columns.extend([
                format!("{delta_delay_sample:+.17e}"),
                format!("{:+.17e}", delta_phase.rate_sps),
                format!("{:+.17e}", delta_phase.accel_sps2),
                format!("{:+.17e}", delta_phase.jerk_sps3),
                format!("{:+.17e}", delta_phase.snap_sps4),
                format!("{phase_from_start_deg:+.17e}"),
                format!("{phase_after_linear_deg:+.17e}"),
                format!("{phase_after_quadratic_deg:+.17e}"),
            ]);
        }
        if columns.len() != expected_columns {
            return Err(format!(
                "TSV column mismatch: field=column_count model=all row={row_label} expected={expected_columns} actual={}",
                columns.len()
            )
            .into());
        }
        writeln!(writer, "{}", columns.join("\t"))?;
    }
    writer.flush()?;

    let indices = diagnostic_checkpoint_indices(input.total_frames);
    let mut checkpoints = Vec::with_capacity(indices.len());
    let first_active = active_series[0];
    let first_deltas: Vec<DelayDerivatives> = series
        .iter()
        .enumerate()
        .map(|(variant_index, values)| {
            let delta = values[0].sub(first_active);
            validate_derivatives(delta, "delta", &variants[variant_index].key(), "data-start")?;
            Ok(delta)
        })
        .collect::<Result<Vec<_>, DynError>>()?;
    for (position, frame_idx) in indices.iter().copied().enumerate() {
        let scan_elapsed_s = (frame_idx as f64 + 0.5) * input.frame_dt_s;
        let row_label = format!("checkpoint:{position}@{scan_elapsed_s:.9}s");
        let data_elapsed_process_s = input.total_skip_s + scan_elapsed_s;
        let data_mjd_utc = input.geometry.reference_mjd_utc + data_elapsed_process_s / 86400.0;
        let model_elapsed_process_s = input.start_elapsed_process_s + scan_elapsed_s;
        let model_mjd_utc = input.geometry.reference_mjd_utc + model_elapsed_process_s / 86400.0;
        for (field, number) in [
            ("scan_elapsed_s", scan_elapsed_s),
            ("data_elapsed_process_s", data_elapsed_process_s),
            ("data_mjd_utc", data_mjd_utc),
            ("model_elapsed_process_s", model_elapsed_process_s),
            ("model_mjd_utc", model_mjd_utc),
        ] {
            require_finite(number, field, "active", &row_label)?;
        }
        let geom_active = evaluate_at(input.geometry, active, model_elapsed_process_s);
        let clock = DelayDerivatives::from_epoch_polynomial(
            input.clock_epoch_coefficients.delay_s,
            input.clock_epoch_coefficients.rate_sps,
            input.clock_epoch_coefficients.accel_sps2,
            input.clock_epoch_coefficients.jerk_sps3,
            input.clock_epoch_coefficients.snap_sps4,
            model_elapsed_process_s,
        );
        let other = DelayDerivatives::from_epoch_polynomial(
            input.other_epoch_coefficients.delay_s,
            input.other_epoch_coefficients.rate_sps,
            input.other_epoch_coefficients.accel_sps2,
            input.other_epoch_coefficients.jerk_sps3,
            input.other_epoch_coefficients.snap_sps4,
            model_elapsed_process_s,
        );
        let total = geom_active.add(clock).add(other);
        validate_derivatives(geom_active, "geom", "active", &row_label)?;
        validate_derivatives(clock, "clock", "active", &row_label)?;
        validate_derivatives(other, "other", "active", &row_label)?;
        validate_derivatives(total, "total", "active", &row_label)?;
        let sky = az_el_and_basis(input, active, model_elapsed_process_s);
        validate_az_el_basis(sky, "active", &row_label)?;
        let variant_points = variants
            .iter()
            .enumerate()
            .map(
                |(variant_idx, variant)| -> Result<VariantCheckpoint, DynError> {
                    let model = variant.key();
                    let geom_value = evaluate_at(input.geometry, *variant, model_elapsed_process_s);
                    let delta = geom_value.sub(geom_active);
                    let delta0 = first_deltas[variant_idx];
                    let delta_phase_from_start_deg =
                        -360.0 * input.obs_hz * (delta.delay_s - delta0.delay_s);
                    let delta_phase_after_linear_deg = -360.0
                        * input.obs_hz
                        * (delta.delay_s - detrends[variant_idx].linear.evaluate(scan_elapsed_s));
                    let delta_phase_after_quadratic_deg = -360.0
                        * input.obs_hz
                        * (delta.delay_s
                            - detrends[variant_idx].quadratic.evaluate(scan_elapsed_s));
                    validate_derivatives(geom_value, "geom", &model, &row_label)?;
                    validate_derivatives(delta, "delta", &model, &row_label)?;
                    for (field, number) in [
                        ("delta_phase_from_start_deg", delta_phase_from_start_deg),
                        ("delta_phase_after_linear_deg", delta_phase_after_linear_deg),
                        (
                            "delta_phase_after_quadratic_deg",
                            delta_phase_after_quadratic_deg,
                        ),
                    ] {
                        require_finite(number, field, &model, &row_label)?;
                    }
                    Ok(VariantCheckpoint {
                        key: model,
                        geom: geom_value,
                        delta,
                        delta_phase_from_start_deg,
                        delta_phase_after_linear_deg,
                        delta_phase_after_quadratic_deg,
                    })
                },
            )
            .collect::<Result<Vec<_>, DynError>>()?;
        checkpoints.push(CheckpointReport {
            label: checkpoint_label(position, indices.len()),
            frame_idx,
            scan_elapsed_s,
            data_elapsed_process_s,
            data_mjd_utc,
            model_elapsed_process_s,
            model_mjd_utc,
            geom: geom_active,
            clock,
            other,
            total,
            ant1_az_deg: sky.az1_deg,
            ant1_el_deg: sky.el1_deg,
            ant1_csc_el: sky.csc1,
            ant1_csc_el_valid: sky.csc1_valid,
            ant2_az_deg: sky.az2_deg,
            ant2_el_deg: sky.el2_deg,
            ant2_csc_el: sky.csc2,
            ant2_csc_el_valid: sky.csc2_valid,
            u_lambda: sky.u_lambda,
            v_lambda: sky.v_lambda,
            w_lambda: sky.w_lambda,
            variants: variant_points,
        });
    }

    Ok(DiagnosticReport {
        path: path.to_path_buf(),
        rows,
        variants: variants.len(),
        checkpoints,
    })
}

#[cfg(test)]
mod tests {
    use super::{diagnostic_checkpoint_indices, sampling_times};

    #[test]
    fn checkpoints_cover_scan_without_duplicates() {
        assert!(diagnostic_checkpoint_indices(0).is_empty());
        assert_eq!(diagnostic_checkpoint_indices(1), vec![0]);
        assert_eq!(diagnostic_checkpoint_indices(2), vec![0, 1]);
        assert_eq!(
            diagnostic_checkpoint_indices(3600),
            vec![0, 899, 1799, 2699, 3599]
        );
    }

    #[test]
    fn sampling_times_never_extend_past_duration() {
        assert_eq!(sampling_times(0.0).unwrap(), vec![0.0]);
        assert_eq!(sampling_times(0.2).unwrap(), vec![0.0, 0.2]);
        assert_eq!(sampling_times(2.0).unwrap(), vec![0.0, 1.0, 2.0]);
        assert_eq!(sampling_times(2.1).unwrap(), vec![0.0, 1.0, 2.0, 2.1]);
        assert_eq!(sampling_times(2.0 + 5.0e-13).unwrap(), vec![0.0, 1.0, 2.0]);
    }
}
