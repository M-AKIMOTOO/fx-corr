use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::utils::DynError;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ProcessEntry {
    pub epoch: String,
    pub skip_sec: f64,
    pub length_sec: Option<f64>,
    pub object: Option<String>,
    pub stations: Option<String>,
    pub ra: Option<String>,
    pub dec: Option<String>,
}

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct IFileData {
    pub ra: String,
    pub dec: String,
    pub epoch: Option<String>,
    pub source: Option<String>,
    pub stream_label: Option<String>,
    pub fft: Option<usize>,
    pub sampling_hz: Option<f64>,
    pub ant1_bit: Option<usize>,
    pub ant2_bit: Option<usize>,
    pub ant1_level: Option<String>,
    pub ant2_level: Option<String>,
    pub ant1_shuffle: Option<String>,
    pub ant2_shuffle: Option<String>,
    pub obsfreq_mhz: Option<f64>,
    pub clock_delay_s: Option<f64>,
    pub clock_rate_sps: Option<f64>,
    pub ant1_clock_delay_s: Option<f64>,
    pub ant2_clock_delay_s: Option<f64>,
    pub ant1_clock_rate_sps: Option<f64>,
    pub ant2_clock_rate_sps: Option<f64>,
    pub ant1_sideband: Option<String>,
    pub ant2_sideband: Option<String>,
    pub ant1_rotation_hz: Option<f64>,
    pub ant2_rotation_hz: Option<f64>,
    pub ant1_center_mhz: Option<f64>,
    pub ant2_center_mhz: Option<f64>,
    pub ant1_bw_mhz: Option<f64>,
    pub ant2_bw_mhz: Option<f64>,
    pub ant1_station_name: Option<String>,
    pub ant2_station_name: Option<String>,
    pub ant1_station_key: Option<String>,
    pub ant2_station_key: Option<String>,
    pub ant1_ecef_m: Option<[f64; 3]>,
    pub ant2_ecef_m: Option<[f64; 3]>,
    pub process_epochs: Vec<String>,
    pub process_skip_sec: Option<f64>,
    pub process_length_sec: Option<f64>,
    pub processes: Vec<ProcessEntry>,
}

fn parse_optional_f64(params: &HashMap<String, String>, keys: &[&str]) -> Result<Option<f64>, DynError> {
    for key in keys {
        if let Some(value) = params.get(*key) {
            return Ok(Some(value.trim().parse::<f64>()?));
        }
    }
    Ok(None)
}

fn parse_clock_pair(clock_raw: &str) -> Result<(Option<f64>, Option<f64>), DynError> {
    let parts: Vec<&str> = clock_raw
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .collect();
    if parts.is_empty() { return Ok((None, None)); }
    let delay = Some(parts[0].parse::<f64>()?);
    let rate = if parts.len() >= 2 { Some(parts[1].parse::<f64>()?) } else { None };
    Ok((delay, rate))
}

fn parse_xyz_triplet(raw: &str) -> Result<[f64; 3], DynError> {
    let values: Vec<f64> = raw
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>())
        .collect::<Result<Vec<_>, _>>()?;
    if values.len() != 3 {
        return Err("ECEF xyz must contain exactly 3 values".into());
    }
    Ok([values[0], values[1], values[2]])
}

fn parse_optional_xyz(
    params: &HashMap<String, String>,
    triplet_keys: &[&str],
    x_keys: &[&str],
    y_keys: &[&str],
    z_keys: &[&str],
) -> Result<Option<[f64; 3]>, DynError> {
    for key in triplet_keys {
        if let Some(v) = params.get(*key) {
            return Ok(Some(parse_xyz_triplet(v)?));
        }
    }
    let x = parse_optional_f64(params, x_keys)?;
    let y = parse_optional_f64(params, y_keys)?;
    let z = parse_optional_f64(params, z_keys)?;
    Ok(match (x, y, z) {
        (Some(xv), Some(yv), Some(zv)) => Some([xv, yv, zv]),
        _ => None,
    })
}

fn parse_ifile_kv(path: &PathBuf) -> Result<IFileData, DynError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut params = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.splitn(2, '#').next().unwrap_or("").trim();
        if line.is_empty() || line.starts_with(';') { continue; }
        if let Some(index) = line.find('=') {
            let (key, value) = line.split_at(index);
            let key = key.trim().to_ascii_lowercase().replace('_', "");
            let value = value.trim_start_matches('=').trim().trim_matches('"').trim_matches('\'').to_string();
            params.insert(key, value);
        }
    }
    let mut clock_delay_s = parse_optional_f64(&params, &["clockdelay", "clockoffset"])?;
    let mut clock_rate_sps = parse_optional_f64(&params, &["clockrate", "clockdrift"])?;
    let ant1_clock_delay_s = parse_optional_f64(&params, &["ant1clockdelay", "clock1delay", "station1clockdelay"])?;
    let ant2_clock_delay_s = parse_optional_f64(&params, &["ant2clockdelay", "clock2delay", "station2clockdelay"])?;
    let ant1_clock_rate_sps = parse_optional_f64(&params, &["ant1clockrate", "clock1rate", "station1clockrate"])?;
    let ant2_clock_rate_sps = parse_optional_f64(&params, &["ant2clockrate", "clock2rate", "station2clockrate"])?;
    let obsfreq_mhz = parse_optional_f64(&params, &["obsfreq", "skyfreq", "freq"])?;
    let ant1_ecef_m = parse_optional_xyz(
        &params,
        &["ant1xyz", "ant1ecef", "station1xyz"],
        &["ant1x", "station1x", "x1"],
        &["ant1y", "station1y", "y1"],
        &["ant1z", "station1z", "z1"],
    )?;
    let ant2_ecef_m = parse_optional_xyz(
        &params,
        &["ant2xyz", "ant2ecef", "station2xyz"],
        &["ant2x", "station2x", "x2"],
        &["ant2y", "station2y", "y2"],
        &["ant2z", "station2z", "z2"],
    )?;
    if let Some(clock_raw) = params.get("clock") {
        let (d, r) = parse_clock_pair(clock_raw)?;
        if clock_delay_s.is_none() { clock_delay_s = d; }
        if clock_rate_sps.is_none() { clock_rate_sps = r; }
    }
    let ra = params.get("ra").or_else(|| params.get("srcra")).ok_or("ra missing")?.clone();
    let dec = params.get("dec").or_else(|| params.get("srcdec")).ok_or("dec missing")?.clone();
    let epoch = params.get("epoch").or_else(|| params.get("scanepoch")).cloned();
    let process_skip_sec = parse_optional_f64(&params, &["skip", "processskip"])?;
    let process_length_sec = parse_optional_f64(&params, &["length", "processlength"])?;
    Ok(IFileData {
        ra: ra.clone(),
        dec: dec.clone(),
        epoch: epoch.clone(),
        source: params.get("source").or_else(|| params.get("object")).cloned(),
        stream_label: params
            .get("streamlabel")
            .or_else(|| params.get("label"))
            .cloned(),
        fft: params.get("fft").map(|v| v.parse()).transpose()?,
        sampling_hz: params.get("samplinghz").or_else(|| params.get("fs")).map(|v| v.parse()).transpose()?,
        ant1_bit: params.get("bit").map(|v| v.parse()).transpose()?,
        ant2_bit: params.get("bit").map(|v| v.parse()).transpose()?,
        ant1_level: params.get("level").cloned(),
        ant2_level: params.get("level").cloned(),
        ant1_shuffle: params.get("shuffle").cloned(),
        ant2_shuffle: params.get("shuffle").cloned(),
        obsfreq_mhz, clock_delay_s, clock_rate_sps,
        ant1_clock_delay_s,
        ant2_clock_delay_s,
        ant1_clock_rate_sps,
        ant2_clock_rate_sps,
        ant1_sideband: params.get("sideband").cloned(),
        ant2_sideband: params.get("sideband").cloned(),
        ant1_rotation_hz: parse_optional_f64(&params, &["ant1rotationhz", "rotation1hz", "rotation"])?,
        ant2_rotation_hz: parse_optional_f64(&params, &["ant2rotationhz", "rotation2hz", "rotation"])?,
        ant1_center_mhz: None, ant2_center_mhz: None, ant1_bw_mhz: None, ant2_bw_mhz: None,
        ant1_station_name: params.get("ant1").or_else(|| params.get("station1")).cloned(),
        ant2_station_name: params.get("ant2").or_else(|| params.get("station2")).cloned(),
        ant1_station_key: params.get("ant1key").or_else(|| params.get("station1key")).cloned(),
        ant2_station_key: params.get("ant2key").or_else(|| params.get("station2key")).cloned(),
        ant1_ecef_m,
        ant2_ecef_m,
        process_epochs: epoch.clone().into_iter().collect(),
        process_skip_sec,
        process_length_sec,
        processes: vec![ProcessEntry {
            epoch: epoch.clone().unwrap_or_else(|| "2000".to_string()),
            skip_sec: process_skip_sec.unwrap_or(0.0),
            length_sec: process_length_sec,
            object: params.get("source").or_else(|| params.get("object")).cloned(),
            stations: params.get("stations").cloned(),
            ra: Some(ra.clone()),
            dec: Some(dec.clone()),
        }],
    })
}

pub fn parse_ifile(path: &PathBuf) -> Result<IFileData, DynError> {
    if path.extension().and_then(|s| s.to_str()).map(|s| s.eq_ignore_ascii_case("xml")).unwrap_or(false) {
        return crate::xml::parse_xml_schedule(path);
    }
    parse_ifile_kv(path)
}
