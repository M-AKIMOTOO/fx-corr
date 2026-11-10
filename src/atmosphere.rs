//! Station atmosphere inputs and first-order propagation models.
//!
//! Providers return zenith neutral paths (m) and slant TEC (TECU). The FX
//! pipeline keeps the positive ionospheric group delay separate from its
//! negative phase delay. No frequency-bin remapping belongs in this module.

use crate::utils::DynError;
use clap::ValueEnum;
use num_complex::Complex;
use roxmltree::Node;
use std::collections::{HashMap, HashSet};
use std::f64::consts::{PI, TAU};
use std::path::{Path, PathBuf};

pub const C_MPS: f64 = 299_792_458.0;
pub const TEC_COEFFICIENT: f64 = 40.3e16 / C_MPS;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum IonosphereMode {
    Off,
    Input,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum TroposphereMode {
    Off,
    Nominal,
    Input,
}

impl IonosphereMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Input => "TEC-input",
        }
    }
}

impl TroposphereMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Nominal => "nominal-troposphere",
            Self::Input => "input-troposphere/Niell",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StationValues {
    pub stec_tecu: Option<f64>,
    pub vtec_tecu: Option<f64>,
    pub zhd_m: Option<f64>,
    pub zwd_m: Option<f64>,
    pub pressure_hpa: Option<f64>,
    pub temperature_k: Option<f64>,
    pub vapour_pressure_hpa: Option<f64>,
}

impl StationValues {
    fn has_tec(self) -> bool {
        self.stec_tecu.is_some() || self.vtec_tecu.is_some()
    }

    fn has_trop(self) -> bool {
        self.zhd_m.is_some() || self.pressure_hpa.is_some()
    }

    fn shape(self) -> [bool; 7] {
        [
            self.stec_tecu.is_some(),
            self.vtec_tecu.is_some(),
            self.zhd_m.is_some(),
            self.zwd_m.is_some(),
            self.pressure_hpa.is_some(),
            self.temperature_k.is_some(),
            self.vapour_pressure_hpa.is_some(),
        ]
    }

    fn validate(self) -> Result<(), DynError> {
        for (name, value, low, high) in [
            ("stec_tecu", self.stec_tecu, 0.0, f64::MAX),
            ("vtec_tecu", self.vtec_tecu, 0.0, f64::MAX),
            ("zhd_m", self.zhd_m, 0.0, 10.0),
            ("zwd_m", self.zwd_m, 0.0, 10.0),
            ("pressure_hpa", self.pressure_hpa, 0.0, 1200.0),
            ("temperature_k", self.temperature_k, 150.0, 350.0),
            ("vapour_pressure_hpa", self.vapour_pressure_hpa, 0.0, 120.0),
        ] {
            if let Some(v) = value {
                if !v.is_finite() || v < low || v > high {
                    return Err(
                        format!("atmosphere {name}={v}: expected finite {low}..{high}").into(),
                    );
                }
            }
        }
        if self.stec_tecu.is_some() && self.vtec_tecu.is_some() {
            return Err("specify only one of stec_tecu and vtec_tecu per station".into());
        }
        let zenith = self.zhd_m.is_some() || self.zwd_m.is_some();
        let weather = self.pressure_hpa.is_some()
            || self.temperature_k.is_some()
            || self.vapour_pressure_hpa.is_some();
        if zenith && weather {
            return Err("atmosphere: zenith paths and weather values cannot be mixed".into());
        }
        if zenith && (self.zhd_m.is_none() || self.zwd_m.is_none()) {
            return Err(
                "atmosphere: both zhd_m and zwd_m are required (use 0 for zero wet path)".into(),
            );
        }
        if weather
            && (self.pressure_hpa.is_none()
                || self.temperature_k.is_none()
                || self.vapour_pressure_hpa.is_none())
        {
            return Err(
                "atmosphere: pressure_hpa, temperature_k and vapour_pressure_hpa are all required"
                    .into(),
            );
        }
        if let (Some(p), Some(e)) = (self.pressure_hpa, self.vapour_pressure_hpa) {
            if e > p {
                return Err("atmosphere: water vapour pressure exceeds total pressure".into());
            }
        }
        if !self.has_tec() && !self.has_trop() {
            return Err("atmosphere station has no TEC, zenith paths or weather values".into());
        }
        Ok(())
    }

    fn interpolate(self, other: Self, fraction: f64) -> Self {
        let lerp = |a: Option<f64>, b: Option<f64>| a.zip(b).map(|(x, y)| x + (y - x) * fraction);
        Self {
            stec_tecu: lerp(self.stec_tecu, other.stec_tecu),
            vtec_tecu: lerp(self.vtec_tecu, other.vtec_tecu),
            zhd_m: lerp(self.zhd_m, other.zhd_m),
            zwd_m: lerp(self.zwd_m, other.zwd_m),
            pressure_hpa: lerp(self.pressure_hpa, other.pressure_hpa),
            temperature_k: lerp(self.temperature_k, other.temperature_k),
            vapour_pressure_hpa: lerp(self.vapour_pressure_hpa, other.vapour_pressure_hpa),
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AtmosphereConfig {
    pub file: Option<PathBuf>,
    pub ionosphere: Option<IonosphereMode>,
    pub troposphere: Option<TroposphereMode>,
    pub shell_height_km: Option<f64>,
    pub min_elevation_deg: Option<f64>,
    pub stations: [Option<StationValues>; 2],
}

pub fn parse_xml(
    root: Node<'_, '_>,
    station1: Node<'_, '_>,
    station2: Node<'_, '_>,
) -> Result<AtmosphereConfig, DynError> {
    fn atmosphere<'a, 'input>(parent: Node<'a, 'input>) -> Option<Node<'a, 'input>> {
        parent.children().find(|n| n.has_tag_name("atmosphere"))
    }
    for parent in [root, station1, station2] {
        if parent
            .children()
            .filter(|n| n.has_tag_name("atmosphere"))
            .count()
            > 1
        {
            return Err("only one atmosphere element is allowed in each schedule/station".into());
        }
        if let Some(node) = atmosphere(parent) {
            if node.children().any(|n| n.is_element()) {
                return Err("atmosphere settings must use attributes, not child elements".into());
            }
        }
    }
    let mut config = AtmosphereConfig::default();
    if let Some(node) = atmosphere(root) {
        for attr in node.attributes() {
            match attr.name() {
                "file" => config.file = Some(PathBuf::from(attr.value())),
                "ionosphere" => {
                    config.ionosphere = Some(
                        IonosphereMode::from_str(attr.value(), true).map_err(|e| e.to_string())?,
                    );
                }
                "troposphere" => {
                    config.troposphere = Some(
                        TroposphereMode::from_str(attr.value(), true).map_err(|e| e.to_string())?,
                    );
                }
                "shell-height-km" => config.shell_height_km = Some(attr.value().parse()?),
                "min-elevation-deg" => config.min_elevation_deg = Some(attr.value().parse()?),
                name => {
                    return Err(format!("unknown schedule/atmosphere attribute '{name}'").into())
                }
            }
        }
    }
    for (i, node) in [station1, station2].into_iter().enumerate() {
        if let Some(n) = atmosphere(node) {
            let mut values = StationValues::default();
            for attr in n.attributes() {
                let field = match attr.name() {
                    "stec-tecu" => &mut values.stec_tecu,
                    "vtec-tecu" => &mut values.vtec_tecu,
                    "zhd-m" => &mut values.zhd_m,
                    "zwd-m" => &mut values.zwd_m,
                    "pressure-hpa" => &mut values.pressure_hpa,
                    "temperature-k" => &mut values.temperature_k,
                    "vapour-pressure-hpa" => &mut values.vapour_pressure_hpa,
                    name => {
                        return Err(format!("unknown station/atmosphere attribute '{name}'").into())
                    }
                };
                *field = Some(attr.value().parse()?);
            }
            values.validate()?;
            config.stations[i] = Some(values);
        }
    }
    Ok(config)
}

#[derive(Clone, Debug)]
struct TimedValues {
    mjd_utc: f64,
    values: StationValues,
}

#[derive(Clone, Debug)]
enum StationSeries {
    Constant(StationValues),
    TimeSeries(Vec<TimedValues>),
}

impl StationSeries {
    fn values_at(&self, mjd_utc: f64) -> StationValues {
        match self {
            Self::Constant(v) => *v,
            Self::TimeSeries(rows) => {
                let right = rows.partition_point(|r| r.mjd_utc < mjd_utc);
                if right == 0 {
                    return rows[0].values;
                }
                if right == rows.len() {
                    return rows[right - 1].values;
                }
                let (a, b) = (&rows[right - 1], &rows[right]);
                a.values
                    .interpolate(b.values, (mjd_utc - a.mjd_utc) / (b.mjd_utc - a.mjd_utc))
            }
        }
    }

    fn validate_span(&self, station: &str, first: f64, last: f64) -> Result<(), DynError> {
        if let Self::TimeSeries(rows) = self {
            // Sub-microsecond UTC/MJD conversion roundoff only; never extrapolate observations.
            let tolerance = 1e-11;
            if first < rows[0].mjd_utc - tolerance
                || last > rows[rows.len() - 1].mjd_utc + tolerance
            {
                return Err(format!(
                    "atmosphere station {station}: UTC MJD {first:.12}..{last:.12} outside input coverage {:.12}..{:.12}; extrapolation is disabled",
                    rows[0].mjd_utc, rows[rows.len() - 1].mjd_utc
                ).into());
            }
        }
        Ok(())
    }
}

fn read_csv(path: &Path) -> Result<HashMap<String, StationSeries>, DynError> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty() && !line.trim_start().starts_with('#'));
    let (_, header) = lines.next().ok_or("empty atmosphere CSV")?;
    let columns: Vec<_> = header
        .trim_start_matches('\u{feff}')
        .split(',')
        .map(str::trim)
        .collect();
    let allowed = [
        "station",
        "mjd_utc",
        "stec_tecu",
        "vtec_tecu",
        "zhd_m",
        "zwd_m",
        "pressure_hpa",
        "temperature_k",
        "vapour_pressure_hpa",
    ];
    let mut seen = HashSet::new();
    for column in &columns {
        if !allowed.contains(column) || !seen.insert(*column) {
            return Err(format!("unknown or duplicate atmosphere CSV column '{column}'").into());
        }
    }
    if !seen.contains("station") || !seen.contains("mjd_utc") {
        return Err("atmosphere CSV requires station,mjd_utc columns; blank mjd_utc means an explicit constant".into());
    }
    let mut grouped: HashMap<String, Vec<(Option<f64>, StationValues)>> = HashMap::new();
    for (line_index, line) in lines {
        let parse = || -> Result<(String, Option<f64>, StationValues), DynError> {
            let cells: Vec<_> = line.split(',').map(str::trim).collect();
            if cells.len() != columns.len() {
                return Err(
                    "column count differs from header (quoted CSV fields are not supported)".into(),
                );
            }
            let get = |key| {
                columns
                    .iter()
                    .position(|c| *c == key)
                    .map(|i| cells[i])
                    .filter(|s| !s.is_empty())
            };
            let value = |key| get(key).map(str::parse::<f64>).transpose();
            let station = get("station").ok_or("missing station key")?.to_string();
            let mjd = value("mjd_utc")?;
            if let Some(m) = mjd {
                if !m.is_finite() || !(40_000.0..100_000.0).contains(&m) {
                    return Err(
                        "mjd_utc must be a finite UTC Modified Julian Date, 40000..100000".into(),
                    );
                }
            }
            let values = StationValues {
                stec_tecu: value("stec_tecu")?,
                vtec_tecu: value("vtec_tecu")?,
                zhd_m: value("zhd_m")?,
                zwd_m: value("zwd_m")?,
                pressure_hpa: value("pressure_hpa")?,
                temperature_k: value("temperature_k")?,
                vapour_pressure_hpa: value("vapour_pressure_hpa")?,
            };
            values.validate()?;
            Ok((station, mjd, values))
        };
        let (key, mjd, values) =
            parse().map_err(|e| format!("{}:{}: {e}", path.display(), line_index + 1))?;
        grouped.entry(key).or_default().push((mjd, values));
    }
    let mut result = HashMap::new();
    for (station, mut rows) in grouped {
        if rows.len() == 1 && rows[0].0.is_none() {
            result.insert(station, StationSeries::Constant(rows[0].1));
            continue;
        }
        if rows.len() < 2 || rows.iter().any(|r| r.0.is_none()) {
            return Err(format!(
                "atmosphere {station}: use one row with blank mjd_utc, or at least two dated rows"
            )
            .into());
        }
        rows.sort_by(|a, b| a.0.unwrap().total_cmp(&b.0.unwrap()));
        for pair in rows.windows(2) {
            if pair[0].0 == pair[1].0 || pair[0].1.shape() != pair[1].1.shape() {
                return Err(format!(
                    "atmosphere {station}: duplicate epochs or changing input fields"
                )
                .into());
            }
        }
        result.insert(
            station,
            StationSeries::TimeSeries(
                rows.into_iter()
                    .map(|(m, values)| TimedValues {
                        mjd_utc: m.unwrap(),
                        values,
                    })
                    .collect(),
            ),
        );
    }
    Ok(result)
}

#[derive(Clone, Copy, Debug)]
pub struct Site {
    pub latitude_rad: f64,
    pub height_m: f64,
}

impl Site {
    pub fn from_ecef(xyz: [f64; 3]) -> Self {
        let (latitude_rad, longitude) = crate::geom::ecef_geodetic_lat_lon(xyz);
        let e2 = (1.0 / 298.257_223_563) * (2.0 - 1.0 / 298.257_223_563);
        let height_m = xyz[0] * latitude_rad.cos() * longitude.cos()
            + xyz[1] * latitude_rad.cos() * longitude.sin()
            + xyz[2] * latitude_rad.sin()
            - 6_378_137.0 * (1.0 - e2 * latitude_rad.sin().powi(2)).sqrt();
        Self {
            latitude_rad,
            height_m,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct AtmosphereSample {
    /// Non-dispersive input troposphere only; legacy nominal mode stays in the geometric table.
    pub troposphere_delay_s: f64,
    /// 40.3 * (STEC2-STEC1) * 1e16 / c, in seconds * Hz^2.
    pub ionosphere_s_hz2: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct AtmosphereKnot {
    pub elapsed_s: f64,
    pub sample: AtmosphereSample,
}

/// The table includes both one-second geometry updates and every input epoch.
/// Preserve fractional final intervals and short-cadence TEC/weather changes.
pub fn interpolate_table(table: &[AtmosphereKnot], elapsed_s: f64) -> AtmosphereSample {
    let right = table.partition_point(|k| k.elapsed_s < elapsed_s);
    if right == 0 {
        return table[0].sample;
    }
    if right == table.len() {
        return table[right - 1].sample;
    }
    let (a, b) = (&table[right - 1], &table[right]);
    let fraction = (elapsed_s - a.elapsed_s) / (b.elapsed_s - a.elapsed_s);
    AtmosphereSample {
        troposphere_delay_s: a.sample.troposphere_delay_s
            + (b.sample.troposphere_delay_s - a.sample.troposphere_delay_s) * fraction,
        ionosphere_s_hz2: a.sample.ionosphere_s_hz2
            + (b.sample.ionosphere_s_hz2 - a.sample.ionosphere_s_hz2) * fraction,
    }
}

#[derive(Clone, Debug)]
pub struct AtmosphereRuntime {
    pub ionosphere: IonosphereMode,
    pub troposphere: TroposphereMode,
    pub shell_height_km: f64,
    pub min_elevation_deg: f64,
    pub input_label: String,
    keys: [String; 2],
    sites: [Site; 2],
    stations: [Option<StationSeries>; 2],
}

impl AtmosphereRuntime {
    pub fn new(
        config: &AtmosphereConfig,
        keys: [String; 2],
        xyz: [[f64; 3]; 2],
        vlbi: bool,
    ) -> Result<Self, DynError> {
        let (stations, input_label) = if let Some(path) = config.file.as_ref() {
            let all = read_csv(path)?;
            let mut selected = [None, None];
            for (i, key) in keys.iter().enumerate() {
                // Clone supports autocorrelation: both legs may have the same station key.
                selected[i] = Some(all.get(key).cloned().ok_or_else(|| {
                    format!("{}: atmosphere station key '{key}' missing", path.display())
                })?);
            }
            (selected, path.display().to_string())
        } else {
            (
                config.stations.map(|v| v.map(StationSeries::Constant)),
                "station XML constants".to_string(),
            )
        };
        let first = |i: usize| {
            stations[i]
                .as_ref()
                .map(|s| s.values_at(0.0))
                .unwrap_or_default()
        };
        let ionosphere = config.ionosphere.unwrap_or_else(|| {
            if (0..2).any(|i| first(i).has_tec()) {
                IonosphereMode::Input
            } else {
                IonosphereMode::Off
            }
        });
        let troposphere = config.troposphere.unwrap_or_else(|| {
            if (0..2).any(|i| first(i).has_trop()) {
                TroposphereMode::Input
            } else if vlbi {
                TroposphereMode::Nominal
            } else {
                TroposphereMode::Off
            }
        });
        for (i, key) in keys.iter().enumerate() {
            if ionosphere == IonosphereMode::Input && !first(i).has_tec() {
                return Err(format!("ionosphere=input: station {key} requires stec_tecu or vtec_tecu (use explicit 0 for zero TEC)").into());
            }
            if troposphere == TroposphereMode::Input && !first(i).has_trop() {
                return Err(format!("troposphere=input: station {key} requires zhd_m+zwd_m or pressure_hpa+temperature_k+vapour_pressure_hpa").into());
            }
        }
        let shell_height_km = config.shell_height_km.unwrap_or(450.0);
        let min_elevation_deg = config.min_elevation_deg.unwrap_or(5.0);
        if !shell_height_km.is_finite() || !(100.0..=2000.0).contains(&shell_height_km) {
            return Err("ionosphere shell height must be 100..2000 km".into());
        }
        if !min_elevation_deg.is_finite() || !(3.0..=90.0).contains(&min_elevation_deg) {
            return Err("atmosphere minimum elevation must be 3..90 degrees".into());
        }
        let runtime = Self {
            ionosphere,
            troposphere,
            shell_height_km,
            min_elevation_deg,
            input_label,
            keys,
            sites: xyz.map(Site::from_ecef),
            stations,
        };
        if runtime.needs_elevation()
            && xyz.iter().any(|p| {
                let radius = p[0].hypot(p[1]).hypot(p[2]);
                !radius.is_finite() || !(6.0e6..7.0e6).contains(&radius)
            })
        {
            return Err("VTEC/troposphere mapping requires finite Earth-surface station ECEF coordinates in metres".into());
        }
        Ok(runtime)
    }

    pub fn has_input(&self) -> bool {
        self.ionosphere == IonosphereMode::Input || self.troposphere == TroposphereMode::Input
    }

    pub fn needs_elevation(&self) -> bool {
        self.troposphere == TroposphereMode::Input
            || (self.ionosphere == IonosphereMode::Input
                && self
                    .stations
                    .iter()
                    .flatten()
                    .any(|s| s.values_at(0.0).vtec_tecu.is_some()))
    }

    pub fn label(&self) -> String {
        if self.ionosphere == IonosphereMode::Off {
            self.troposphere.label().into()
        } else if self.troposphere == TroposphereMode::Off {
            "ionosphere/TEC-input".into()
        } else {
            format!(
                "{}+ionosphere/{}",
                self.troposphere.label(),
                self.ionosphere.label()
            )
        }
    }

    pub fn validate_span(&self, first_mjd: f64, last_mjd: f64) -> Result<(), DynError> {
        if self.has_input() {
            for (key, series) in self.keys.iter().zip(&self.stations) {
                if let Some(s) = series {
                    s.validate_span(key, first_mjd, last_mjd)?;
                }
            }
        }
        Ok(())
    }

    pub fn input_epochs(&self) -> Vec<f64> {
        self.stations
            .iter()
            .flatten()
            .flat_map(|s| match s {
                StationSeries::Constant(_) => Vec::new(),
                StationSeries::TimeSeries(rows) => rows.iter().map(|r| r.mjd_utc).collect(),
            })
            .collect()
    }

    pub fn log_inputs(&self, mjd: f64) {
        println!("[atmosphere] ionosphere={} troposphere={} provider={} shell_height_km={} min_elevation_deg={}", self.ionosphere.label(), self.troposphere.label(), self.input_label, self.shell_height_km, self.min_elevation_deg);
        for (key, series) in self.keys.iter().zip(&self.stations) {
            if let Some(s) = series {
                println!(
                    "[atmosphere:input] station={} mjd_utc={:.12} values={:?} provider={:?}",
                    key,
                    mjd,
                    s.values_at(mjd),
                    s
                );
            }
        }
    }

    pub fn evaluate(&self, mjd: f64, elevation: [f64; 2]) -> Result<AtmosphereSample, DynError> {
        let mut tec = [0.0; 2];
        let mut path = [0.0; 2];
        for i in 0..2 {
            let values = self.stations[i]
                .as_ref()
                .map(|s| s.values_at(mjd))
                .unwrap_or_default();
            let needs_mapping = self.troposphere == TroposphereMode::Input
                || (self.ionosphere == IonosphereMode::Input && values.vtec_tecu.is_some());
            if needs_mapping
                && (!elevation[i].is_finite() || elevation[i].to_degrees() < self.min_elevation_deg)
            {
                return Err(format!("atmosphere station {} at UTC MJD {mjd:.12}: elevation {:.6} deg below {:.6} deg model limit", self.keys[i], elevation[i].to_degrees(), self.min_elevation_deg).into());
            }
            if self.ionosphere == IonosphereMode::Input {
                tec[i] = if let Some(stec) = values.stec_tecu {
                    stec
                } else {
                    values.vtec_tecu.unwrap()
                        * thin_shell_mapping(
                            elevation[i],
                            self.sites[i].height_m,
                            self.shell_height_km * 1000.0,
                        )?
                };
            }
            if self.troposphere == TroposphereMode::Input {
                let (zhd, zwd) = if let (Some(dry), Some(wet)) = (values.zhd_m, values.zwd_m) {
                    (dry, wet)
                } else {
                    saastamoinen_zenith(
                        self.sites[i],
                        values.pressure_hpa.unwrap(),
                        values.temperature_k.unwrap(),
                        values.vapour_pressure_hpa.unwrap(),
                    )
                };
                let (dry_map, wet_map) = niell_mapping(self.sites[i], elevation[i], mjd)?;
                path[i] = zhd * dry_map + zwd * wet_map;
            }
        }
        let sample = AtmosphereSample {
            troposphere_delay_s: (path[1] - path[0]) / C_MPS,
            ionosphere_s_hz2: TEC_COEFFICIENT * (tec[1] - tec[0]),
        };
        if !sample.ionosphere_s_hz2.is_finite() || !sample.troposphere_delay_s.is_finite() {
            return Err("non-finite atmosphere delay; check TEC and station coordinates".into());
        }
        Ok(sample)
    }
}

pub fn thin_shell_mapping(
    elevation: f64,
    station_height_m: f64,
    shell_height_m: f64,
) -> Result<f64, DynError> {
    let ratio = (6_378_137.0 + station_height_m) / (6_378_137.0 + shell_height_m);
    if ratio >= 1.0 || ratio <= 0.0 || !ratio.is_finite() {
        return Err("ionosphere thin shell must be above the station".into());
    }
    Ok(1.0 / (1.0 - (ratio * elevation.cos()).powi(2)).sqrt())
}

/// Saastamoinen hydrostatic/wet zenith paths, m. Water vapour pressure is
/// supplied explicitly in hPa rather than inferring humidity silently.
pub fn saastamoinen_zenith(site: Site, p_hpa: f64, t_k: f64, e_hpa: f64) -> (f64, f64) {
    let dry = 0.0022768 * p_hpa
        / (1.0 - 0.00266 * (2.0 * site.latitude_rad).cos() - 0.00028 * site.height_m / 1000.0);
    let wet = 0.002277 * (1255.0 / t_k + 0.05) * e_hpa;
    (dry, wet)
}

/// Niell (1996), corrected continued fraction, latitude interpolation and
/// southern-hemisphere season. Numerical coefficients are from Table 3.
pub fn niell_mapping(site: Site, elevation: f64, mjd: f64) -> Result<(f64, f64), DynError> {
    const COEF: [[f64; 5]; 9] = [
        [
            1.2769934e-3,
            1.2683230e-3,
            1.2465397e-3,
            1.2196049e-3,
            1.2045996e-3,
        ],
        [
            2.9153695e-3,
            2.9152299e-3,
            2.9288445e-3,
            2.9022565e-3,
            2.9024912e-3,
        ],
        [
            62.610505e-3,
            62.837393e-3,
            63.721774e-3,
            63.824265e-3,
            64.258455e-3,
        ],
        [0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5],
        [0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5],
        [0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5],
        [
            5.8021897e-4,
            5.6794847e-4,
            5.8118019e-4,
            5.9727542e-4,
            6.1641693e-4,
        ],
        [
            1.4275268e-3,
            1.5138625e-3,
            1.4572752e-3,
            1.5007428e-3,
            1.7599082e-3,
        ],
        [
            4.3472961e-2,
            4.6729510e-2,
            4.3908931e-2,
            4.4626982e-2,
            5.4736038e-2,
        ],
    ];
    let mut year = 2000 + ((mjd - 51544.0) / 365.2425).floor() as i32;
    let jan1 = |year| crate::geom::parse_epoch_to_mjd(&format!("{year}/001 00:00:00"));
    if mjd < jan1(year)? {
        year -= 1;
    }
    if mjd >= jan1(year + 1)? {
        year += 1;
    }
    let doy = mjd - jan1(year)? + 1.0;
    let season =
        (TAU * ((doy - 28.0) / 365.25 + if site.latitude_rad < 0.0 { 0.5 } else { 0.0 })).cos();
    let latitude_index = ((site.latitude_rad.abs().to_degrees() - 15.0) / 15.0).clamp(0.0, 4.0);
    let lo = (latitude_index.floor() as usize).min(3);
    let frac = latitude_index - lo as f64;
    let interp = |row: usize| COEF[row][lo] * (1.0 - frac) + COEF[row][lo + 1] * frac;
    let map = |a, b, c| {
        let s = elevation.sin();
        (1.0 + a / (1.0 + b / (1.0 + c))) / (s + a / (s + b / (s + c)))
    };
    let dry = map(
        interp(0) - interp(3) * season,
        interp(1) - interp(4) * season,
        interp(2) - interp(5) * season,
    ) + (1.0 / elevation.sin() - map(2.53e-5, 5.49e-3, 1.14e-3)) * site.height_m / 1000.0;
    let wet = map(interp(6), interp(7), interp(8));
    Ok((dry, wet))
}

/// Positive group delay at the reference RF. The existing non-dispersive FX
/// correction contributes -2*pi*f*group_delay. Add the spectral phase below
/// to replace that phase with the required +2*pi*D/f for X1 * conj(X2).
#[inline]
pub fn ionosphere_group_delay(coefficient: f64, reference_hz: f64) -> f64 {
    if coefficient == 0.0 {
        0.0
    } else {
        coefficient / reference_hz.powi(2)
    }
}

pub fn ionosphere_phase_basis(frequency_hz: f64, reference_hz: f64) -> f64 {
    TAU * (1.0 / frequency_hz + frequency_hz / reference_hz.powi(2))
}

#[inline]
pub fn ionosphere_phase(coefficient: f64, phase_basis: f64) -> Complex<f64> {
    if coefficient == 0.0 {
        Complex::new(1.0, 0.0)
    } else {
        Complex::from_polar(1.0, (coefficient * phase_basis).rem_euclid(2.0 * PI))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(values: [StationValues; 2]) -> AtmosphereRuntime {
        AtmosphereRuntime::new(
            &AtmosphereConfig {
                stations: values.map(Some),
                ..Default::default()
            },
            ["A".into(), "B".into()],
            [[6_378_137.0, 0.0, 0.0]; 2],
            false,
        )
        .unwrap()
    }

    #[test]
    fn ionosphere_group_and_phase_have_opposite_signs_and_physical_scaling() {
        let d = TEC_COEFFICIENT; // One TECU: 40.3 m group path at 100 MHz.
        let f0 = 100e6;
        assert!((ionosphere_group_delay(d, f0) * C_MPS - 40.3).abs() < 1e-13);
        assert!((ionosphere_group_delay(d, f0 * 2.0) * C_MPS - 40.3 / 4.0).abs() < 1e-13);
        for f in [50e6, 80e6, 100e6, 150e6, 300e6] {
            // Incoming X1 conj(X2) has the plasma's negative phase delay.
            let raw_phase = -TAU * d / f;
            let group_alignment_phase = -TAU * f * ionosphere_group_delay(d, f0);
            let spectral = ionosphere_phase(d, ionosphere_phase_basis(f, f0));
            let corrected = Complex::from_polar(1.0, raw_phase + group_alignment_phase) * spectral;
            assert!((corrected - Complex::new(1.0, 0.0)).norm() < 1e-12);
            // The total correction is +2*pi*D/f, whose negative phase slope
            // removes a positive group delay. This catches double counting.
            let total = |freq| {
                -TAU * freq * ionosphere_group_delay(d, f0) + d * ionosphere_phase_basis(freq, f0)
            };
            let derivative = (total(f + 100.0) - total(f - 100.0)) / 200.0;
            assert!((derivative / -TAU - d / f.powi(2)).abs() < 1e-15);
        }
    }

    #[test]
    fn equal_tec_cancels_and_reversing_baseline_conjugates_phase() {
        let a = StationValues {
            stec_tecu: Some(10.0),
            ..Default::default()
        };
        let b = StationValues {
            stec_tecu: Some(11.0),
            ..Default::default()
        };
        let forward = input([a, b]).evaluate(60337.0, [0.0; 2]).unwrap();
        let reverse = input([b, a]).evaluate(60337.0, [0.0; 2]).unwrap();
        assert_eq!(forward.ionosphere_s_hz2, -reverse.ionosphere_s_hz2);
        assert_eq!(
            input([a, a])
                .evaluate(60337.0, [0.0; 2])
                .unwrap()
                .ionosphere_s_hz2,
            0.0
        );
        let p = ionosphere_phase_basis(150e6, 160e6);
        assert!(
            (ionosphere_phase(forward.ionosphere_s_hz2, p).conj()
                - ionosphere_phase(reverse.ionosphere_s_hz2, p))
            .norm()
                < 1e-13
        );
    }

    #[test]
    fn thin_shell_is_unity_at_zenith_and_maps_vertical_to_slant() {
        assert_eq!(thin_shell_mapping(PI / 2.0, 0.0, 450e3).unwrap(), 1.0);
        let mapping = thin_shell_mapping(30f64.to_radians(), 0.0, 450e3).unwrap();
        assert!(mapping > 1.7 && mapping < 1.71);
        let a = StationValues {
            vtec_tecu: Some(10.0),
            ..Default::default()
        };
        let b = StationValues {
            stec_tecu: Some(10.0 * mapping),
            ..Default::default()
        };
        assert!(
            input([a, b])
                .evaluate(60337.0, [30f64.to_radians(); 2])
                .unwrap()
                .ionosphere_s_hz2
                .abs()
                < 1e-5
        );
        assert!(input([a, b]).evaluate(60337.0, [0.0; 2]).is_err());
    }

    #[test]
    fn neutral_paths_are_non_dispersive_and_zenith_mapping_is_normalized() {
        let site = Site {
            latitude_rad: 45f64.to_radians(),
            height_m: 0.0,
        };
        let (dry, wet) = saastamoinen_zenith(site, 1013.25, 293.15, 10.0);
        assert!((dry - 2.3069676).abs() < 1e-12);
        assert!((wet - 0.09861880018761725).abs() < 1e-12);
        let (md, mw) = niell_mapping(site, PI / 2.0, 60337.0).unwrap();
        assert_eq!((md, mw), (1.0, 1.0));
        let a = StationValues {
            zhd_m: Some(2.3),
            zwd_m: Some(0.1),
            ..Default::default()
        };
        let b = StationValues {
            zhd_m: Some(2.4),
            zwd_m: Some(0.2),
            ..Default::default()
        };
        let sample = input([a, b]).evaluate(60337.0, [PI / 2.0; 2]).unwrap();
        assert!((sample.troposphere_delay_s * C_MPS - 0.2).abs() < 1e-14);
        assert_eq!(sample.ionosphere_s_hz2, 0.0);
    }

    #[test]
    fn niell_matches_independent_author_c_vectors_across_seasons_and_hemispheres() {
        let epochs = [
            "2024/028 00:00:00",
            "2024/210 12:00:00",
            "2024/366 23:59:59",
        ]
        .map(|s| crate::geom::parse_epoch_to_mjd(s).unwrap());
        let reference = include_str!("../reference/atmosphere/validation/rtklib-niell.csv");
        for line in reference.lines().skip(1) {
            let row: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
            let site = Site {
                latitude_rad: row[1].to_radians(),
                height_m: row[2],
            };
            let (dry, wet) =
                niell_mapping(site, row[3].to_radians(), epochs[row[0] as usize]).unwrap();
            assert!((dry - row[4]).abs() < 5e-12, "dry {line}: {dry}");
            assert!((wet - row[5]).abs() < 5e-12, "wet {line}: {wet}");
        }
    }

    #[test]
    fn table_preserves_subsecond_input_knots_and_fractional_scan_end() {
        let knot = |t, d| AtmosphereKnot {
            elapsed_s: t,
            sample: AtmosphereSample {
                troposphere_delay_s: d,
                ionosphere_s_hz2: d * 10.0,
            },
        };
        let table = [
            knot(-0.2, 2.0),
            knot(0.0, 4.0),
            knot(0.125, 5.0),
            knot(0.5, 2.0),
        ];
        assert_eq!(interpolate_table(&table, -0.1).troposphere_delay_s, 3.0);
        assert_eq!(interpolate_table(&table, 0.0625).ionosphere_s_hz2, 45.0);
        assert_eq!(interpolate_table(&table, 0.3125).troposphere_delay_s, 3.5);
    }

    #[test]
    fn csv_interpolates_utc_and_rejects_extrapolation_and_incomplete_models() {
        let path = std::env::temp_dir().join(format!(
            "fx-corr-atmosphere-{}-{}.csv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(
            &path,
            "station,mjd_utc,stec_tecu\nA,60337,10\nA,60338,12\nB,,11\n",
        )
        .unwrap();
        let series = read_csv(&path).unwrap();
        assert_eq!(series["A"].values_at(60337.5).stec_tecu, Some(11.0));
        assert!(series["A"].validate_span("A", 60337.0, 60338.0).is_ok());
        assert!(series["A"].validate_span("A", 60336.9, 60338.0).is_err());
        assert!(series["B"].validate_span("B", 60336.9, 60339.0).is_ok());
        for invalid in [
            "station,mjd_utc,stec_tecu\nA,60337,10\nA,60337,12\n",
            "station,mjd_utc,pressure_hpa\nA,,1013\n",
            "station,mjd_utc,stec_tecu,vtec_tecu\nA,,10,10\n",
            "station,mjd_utc,stec_tecu\nA,60337,10\nA,,12\n",
            "station,mjd_utc,stec_tecu\nA,,NaN\n",
            "station,mjd_utc,tec\nA,,10\n",
        ] {
            std::fs::write(&path, invalid).unwrap();
            assert!(read_csv(&path).is_err(), "{invalid}");
        }
        std::fs::remove_file(path).unwrap();
    }
}
