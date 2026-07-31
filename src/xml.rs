use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use roxmltree::{Document, Node};

use crate::ifile::{IFileData, ProcessEntry};
use crate::pulsar;
use crate::utils::DynError;

#[derive(Clone, Debug)]
struct ClockEntry {
    delay_s: f64,
    rate_sps: f64,
    accel_sps2: f64,
    jerk_sps3: f64,
    snap_sps4: f64,
    epoch: Option<String>,
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

fn parse_f64_text(raw: &str) -> Result<f64, DynError> {
    let normalized = raw.trim().replace(
        ['\u{2212}', '\u{ff0d}', '\u{fe63}', '\u{2013}', '\u{2014}'],
        "-",
    );
    normalized.parse::<f64>().map_err(|e| {
        format!("invalid numeric value '{raw}' (normalized '{normalized}'): {e}").into()
    })
}

fn parse_optional_child_f64(node: Node<'_, '_>, tag: &str, default: f64) -> Result<f64, DynError> {
    match child_text(node, tag) {
        Some(v) if !v.is_empty() => parse_f64_text(v),
        _ => Ok(default),
    }
}

fn parse_optional_child_f64_any(
    node: Node<'_, '_>,
    tags: &[&str],
    default: f64,
) -> Result<f64, DynError> {
    for tag in tags {
        if let Some(v) = child_text(node, tag) {
            if !v.is_empty() {
                return parse_f64_text(v);
            }
        }
    }
    Ok(default)
}

fn parse_opt_f64(node: Node<'_, '_>, tag: &str) -> Result<Option<f64>, DynError> {
    Ok(match child_text(node, tag) {
        Some(v) if !v.is_empty() => Some(parse_f64_text(v)?),
        _ => None,
    })
}

fn parse_opt_f64_any(node: Node<'_, '_>, tags: &[&str]) -> Result<Option<f64>, DynError> {
    for tag in tags {
        if let Some(value) = parse_opt_f64(node, tag)? {
            return Ok(Some(value));
        }
    }
    Ok(None)
}

fn parse_station_ecef(node: Node<'_, '_>) -> Result<Option<[f64; 3]>, DynError> {
    let x = parse_opt_f64(node, "pos-x")?;
    let y = parse_opt_f64(node, "pos-y")?;
    let z = parse_opt_f64(node, "pos-z")?;
    Ok(match (x, y, z) {
        (Some(xv), Some(yv), Some(zv)) => Some([xv, yv, zv]),
        _ => None,
    })
}

fn process_station_keys<'a>(
    process: Node<'a, 'a>,
    station_by_key: &HashMap<String, Node<'a, 'a>>,
) -> Result<Vec<String>, DynError> {
    let keys = if let Some(stations) = child_text(process, "stations") {
        stations
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
    } else {
        let mut keys = station_by_key.keys().cloned().collect::<Vec<_>>();
        keys.sort();
        keys
    };
    if keys.len() < 2 {
        return Err("process/stations must contain at least two station keys".into());
    }
    for key in &keys {
        if !station_by_key.contains_key(key) {
            return Err(format!("station key '{}' not found", key).into());
        }
    }
    Ok(keys)
}

fn process_baseline_pairs(
    process: Node<'_, '_>,
    keys: &[String],
    station_by_key: &HashMap<String, Node<'_, '_>>,
) -> Result<Vec<(String, String)>, DynError> {
    let baseline_text = process
        .children()
        .filter(|n| is_tag(*n, "baseline"))
        .filter_map(|n| n.text())
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join(",");

    if baseline_text.is_empty() {
        let mut pairs = Vec::new();
        for i in 0..keys.len() {
            for j in (i + 1)..keys.len() {
                pairs.push((keys[i].clone(), keys[j].clone()));
            }
        }
        return Ok(pairs);
    }

    let raw_tokens = baseline_text
        .split(|c: char| c == ',' || c == ';' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>();
    let key_set = keys.iter().cloned().collect::<HashSet<_>>();
    let mut required = HashSet::new();
    for i in 0..keys.len() {
        for j in (i + 1)..keys.len() {
            required.insert(canonical_pair(&keys[i], &keys[j]));
        }
    }
    let mut seen = HashSet::new();
    let mut pairs = Vec::new();
    for raw in raw_tokens {
        let compact = raw
            .chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .collect::<String>();
        if compact.is_empty() {
            continue;
        }
        if compact.len() % 2 != 0 {
            return Err(format!(
                "process/baseline entry '{}' must contain 2 station keys",
                raw
            )
            .into());
        }
        let chars = compact.chars().collect::<Vec<_>>();
        for pair in chars.chunks(2) {
            let k1 = pair[0].to_string();
            let k2 = pair[1].to_string();
            if k1 == k2 {
                return Err(
                    format!("process/baseline '{}' uses the same station twice", raw).into(),
                );
            }
            if !station_by_key.contains_key(&k1) {
                return Err(format!("station key '{}' in process/baseline not found", k1).into());
            }
            if !station_by_key.contains_key(&k2) {
                return Err(format!("station key '{}' in process/baseline not found", k2).into());
            }
            if !key_set.contains(&k1) || !key_set.contains(&k2) {
                return Err(format!(
                    "process/baseline '{}' is not contained in process/stations",
                    raw
                )
                .into());
            }
            let closure_pair = canonical_pair(&k1, &k2);
            if !seen.insert(closure_pair.clone()) {
                return Err(format!(
                    "process/baseline contains duplicate closure pair {}{}",
                    closure_pair.0, closure_pair.1
                )
                .into());
            }
            pairs.push((k1, k2));
        }
    }
    if pairs.is_empty() {
        return Err("process/baseline did not contain any station-key pairs".into());
    }
    if seen != required {
        let mut missing = required
            .difference(&seen)
            .map(|(a, b)| format!("{a}{b}"))
            .collect::<Vec<_>>();
        missing.sort();
        let mut extra = seen
            .difference(&required)
            .map(|(a, b)| format!("{a}{b}"))
            .collect::<Vec<_>>();
        extra.sort();
        let mut detail = Vec::new();
        if !missing.is_empty() {
            detail.push(format!("missing {}", missing.join(",")));
        }
        if !extra.is_empty() {
            detail.push(format!("extra {}", extra.join(",")));
        }
        return Err(format!(
            "process/baseline must form closure over process/stations ({})",
            detail.join("; ")
        )
        .into());
    }
    Ok(pairs)
}

fn canonical_pair(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_string(), b.to_string())
    } else {
        (b.to_string(), a.to_string())
    }
}

fn push_process_entries(
    out: &mut Vec<ProcessEntry>,
    process: Node<'_, '_>,
    keys: &[String],
    station_by_key: &HashMap<String, Node<'_, '_>>,
    source_radec_by_name: &HashMap<String, (String, String)>,
) -> Result<(), DynError> {
    let epoch = child_text(process, "epoch")
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "2000".to_string());
    let skip_sec = child_text(process, "skip")
        .map(parse_f64_text)
        .transpose()?
        .unwrap_or(0.0);
    let length_sec = child_text(process, "length")
        .map(parse_f64_text)
        .transpose()?;
    let object = child_text(process, "object").map(|s| s.trim().to_string());
    let (ra, dec) = object
        .as_ref()
        .and_then(|obj| source_radec_by_name.get(obj))
        .map(|(ra, dec)| (Some(ra.clone()), Some(dec.clone())))
        .unwrap_or((None, None));

    for (ant1_key, ant2_key) in process_baseline_pairs(process, keys, station_by_key)? {
        let stations = format!("{}{}", ant1_key, ant2_key);
        out.push(ProcessEntry {
            epoch: epoch.clone(),
            skip_sec,
            length_sec,
            object: object.clone(),
            stations: Some(stations),
            ra: ra.clone(),
            dec: dec.clone(),
            ant1_station_key: Some(ant1_key),
            ant2_station_key: Some(ant2_key),
        });
    }
    Ok(())
}

pub fn parse_xml_schedule(path: &PathBuf) -> Result<IFileData, DynError> {
    parse_xml_schedule_for_process(path, None)
}

pub fn parse_xml_schedule_for_process(
    path: &PathBuf,
    selected_process_index: Option<usize>,
) -> Result<IFileData, DynError> {
    let xml = std::fs::read_to_string(path)?;
    let doc = Document::parse(&xml)?;

    let eop_node = doc
        .descendants()
        .find(|n| is_tag(*n, "eop") || is_tag(*n, "earth-orientation"));
    let delay_model_node = doc
        .descendants()
        .find(|n| is_tag(*n, "delay-model") || is_tag(*n, "delay_model"));
    let eop_file = eop_node.and_then(|n| {
        n.attribute("file")
            .or_else(|| n.attribute("path"))
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(PathBuf::from)
            .or_else(|| {
                child_text(n, "file")
                    .or_else(|| child_text(n, "path"))
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(PathBuf::from)
            })
    });
    let dut1_s = eop_node
        .map(|n| parse_opt_f64_any(n, &["dut1"]))
        .transpose()?
        .flatten();
    let tt_utc_s = eop_node
        .map(|n| parse_opt_f64_any(n, &["tt-utc", "tt_utc", "ttutc"]))
        .transpose()?
        .flatten();
    let xp_arcsec = eop_node
        .map(|n| parse_opt_f64_any(n, &["xp", "xpolar", "xp_arcsec"]))
        .transpose()?
        .flatten();
    let yp_arcsec = eop_node
        .map(|n| parse_opt_f64_any(n, &["yp", "ypolar", "yp_arcsec"]))
        .transpose()?
        .flatten();
    let model_time_offset_s = delay_model_node
        .map(|n| {
            parse_opt_f64_any(
                n,
                &[
                    "time-offset",
                    "time_offset",
                    "model-time-offset",
                    "model_time_offset",
                ],
            )
        })
        .transpose()?
        .flatten();

    let process_nodes: Vec<_> = doc
        .descendants()
        .filter(|n| is_tag(*n, "process"))
        .collect();
    if process_nodes.is_empty() {
        return Err("process node not found in XML schedule".into());
    }

    let station_by_key: HashMap<String, Node<'_, '_>> = doc
        .descendants()
        .filter(|n| is_tag(*n, "station"))
        .filter_map(|n| n.attribute("key").map(|k| (k.trim().to_string(), n)))
        .collect();

    let terminal_by_name: HashMap<String, Node<'_, '_>> = doc
        .descendants()
        .filter(|n| is_tag(*n, "terminal"))
        .filter_map(|n| n.attribute("name").map(|name| (name.trim().to_string(), n)))
        .collect();

    let source_radec_by_name: HashMap<String, (String, String)> = doc
        .descendants()
        .filter(|n| is_tag(*n, "source"))
        .filter_map(|n| {
            let name = n.attribute("name")?.trim().to_string();
            let ra = child_text(n, "ra")?.trim().to_string();
            let dec = child_text(n, "dec")?.trim().to_string();
            Some((name, (ra, dec)))
        })
        .collect();

    let mut processes = Vec::new();
    for process in &process_nodes {
        let keys = process_station_keys(*process, &station_by_key)?;
        push_process_entries(
            &mut processes,
            *process,
            &keys,
            &station_by_key,
            &source_radec_by_name,
        )?;
    }
    if processes.is_empty() {
        return Err("no station baselines found in XML process entries".into());
    }

    let selected_process = if let Some(idx) = selected_process_index {
        processes.get(idx).cloned().ok_or_else(|| {
            format!(
                "--process-index {} out of range (0..{})",
                idx,
                processes.len().saturating_sub(1)
            )
        })?
    } else {
        processes[0].clone()
    };
    let ant1_key = selected_process
        .ant1_station_key
        .clone()
        .ok_or("selected process missing ant1 station key")?;
    let ant2_key = selected_process
        .ant2_station_key
        .clone()
        .ok_or("selected process missing ant2 station key")?;

    let station1 = *station_by_key
        .get(&ant1_key)
        .ok_or_else(|| format!("station key '{}' not found", ant1_key))?;
    let station2 = *station_by_key
        .get(&ant2_key)
        .ok_or_else(|| format!("station key '{}' not found", ant2_key))?;

    let ant1_station_name = child_text(station1, "name").map(|s| s.to_string());
    let ant2_station_name = child_text(station2, "name").map(|s| s.to_string());
    let ant1_ecef_m = parse_station_ecef(station1)?;
    let ant2_ecef_m = parse_station_ecef(station2)?;

    let term1_name = child_text(station1, "terminal").map(|s| s.to_string());
    let term2_name = child_text(station2, "terminal").map(|s| s.to_string());
    let term1 = term1_name
        .as_deref()
        .and_then(|name| terminal_by_name.get(name).copied());
    let term2 = term2_name
        .as_deref()
        .and_then(|name| terminal_by_name.get(name).copied());

    let sampling_hz = term1
        .and_then(|t| child_text(t, "speed"))
        .or_else(|| term2.and_then(|t| child_text(t, "speed")))
        .map(|v| v.parse::<f64>())
        .transpose()?;
    let ant1_bit = term1
        .and_then(|t| child_text(t, "bit"))
        .map(|v| v.parse::<usize>())
        .transpose()?;
    let ant2_bit = term2
        .and_then(|t| child_text(t, "bit"))
        .map(|v| v.parse::<usize>())
        .transpose()?;
    let ant1_level = term1
        .and_then(|t| child_text(t, "level"))
        .map(|s| s.to_string());
    let ant2_level = term2
        .and_then(|t| child_text(t, "level"))
        .map(|s| s.to_string());

    let shuffle_by_key: HashMap<String, String> = doc
        .descendants()
        .filter(|n| is_tag(*n, "shuffle"))
        .filter_map(|n| {
            let key = n.attribute("key")?.trim().to_string();
            let value = n.text()?.trim().to_string();
            Some((key, value))
        })
        .collect();
    let ant1_shuffle = shuffle_by_key.get(&ant1_key).cloned();
    let ant2_shuffle = shuffle_by_key.get(&ant2_key).cloned();

    let stream = doc
        .descendants()
        .find(|n| is_tag(*n, "stream"))
        .ok_or("stream node not found in XML schedule")?;
    let stream_label = child_text(stream, "label").map(|s| s.to_string());
    let frequency_hz = child_text(stream, "frequency")
        .ok_or("stream/frequency missing")?
        .parse::<f64>()?;
    let obsfreq_mhz = Some(frequency_hz.abs() / 1e6);
    let fft = child_text(stream, "fft")
        .map(|v| v.parse::<usize>())
        .transpose()?;
    let output_hz = child_text(stream, "output")
        .map(parse_f64_text)
        .transpose()?;
    let inband = child_text(stream, "inband")
        .map(|v| v.parse::<usize>())
        .transpose()?;
    let pulsar = pulsar::parse_xml_node(stream.children().find(|n| is_tag(*n, "pulsar")))?;

    let default_sideband = "USB";
    let special_by_key: HashMap<String, Node<'_, '_>> = stream
        .children()
        .filter(|n| is_tag(*n, "special"))
        .filter_map(|n| n.attribute("key").map(|k| (k.trim().to_string(), n)))
        .collect();
    let sp1 = special_by_key.get(&ant1_key).copied();
    let sp2 = special_by_key.get(&ant2_key).copied();
    let ant1_sideband = Some(
        sp1.and_then(|n| child_text(n, "sideband"))
            .unwrap_or(default_sideband)
            .to_uppercase(),
    );
    let ant2_sideband = Some(
        sp2.and_then(|n| child_text(n, "sideband"))
            .unwrap_or(default_sideband)
            .to_uppercase(),
    );
    let ant1_rotation_hz = sp1
        .map(|n| parse_opt_f64(n, "rotation"))
        .transpose()?
        .flatten();
    let ant2_rotation_hz = sp2
        .map(|n| parse_opt_f64(n, "rotation"))
        .transpose()?
        .flatten();
    let ant1_rotation2_hz = sp1
        .map(|n| parse_opt_f64(n, "rotation2"))
        .transpose()?
        .flatten();
    let ant2_rotation2_hz = sp2
        .map(|n| parse_opt_f64(n, "rotation2"))
        .transpose()?
        .flatten();

    let clock_by_key: HashMap<String, ClockEntry> = doc
        .descendants()
        .filter(|n| is_tag(*n, "clock"))
        .map(|n| -> Result<Option<(String, ClockEntry)>, DynError> {
            let Some(key_attr) = n.attribute("key") else {
                return Ok(None);
            };
            let key = key_attr.trim().to_string();
            let delay = match parse_optional_child_f64(n, "delay", 0.0) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };
            let rate = match parse_optional_child_f64(n, "rate", 0.0) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };
            let accel = match parse_optional_child_f64_any(
                n,
                &["accel", "acel", "acceleration", "accel_sps2"],
                0.0,
            ) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };
            let jerk = match parse_optional_child_f64_any(n, &["jerk", "jerk_sps3"], 0.0) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };
            let snap = match parse_optional_child_f64_any(n, &["snap", "snap_sps4"], 0.0) {
                Ok(v) => v,
                Err(e) => return Err(e),
            };
            let epoch = child_text(n, "epoch").map(|s| s.to_string());
            Ok(Some((
                key,
                ClockEntry {
                    delay_s: delay,
                    rate_sps: rate,
                    accel_sps2: accel,
                    jerk_sps3: jerk,
                    snap_sps4: snap,
                    epoch,
                },
            )))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
        .collect();
    let clock1 = clock_by_key.get(&ant1_key).cloned().unwrap_or(ClockEntry {
        delay_s: 0.0,
        rate_sps: 0.0,
        accel_sps2: 0.0,
        jerk_sps3: 0.0,
        snap_sps4: 0.0,
        epoch: None,
    });
    let clock2 = clock_by_key.get(&ant2_key).cloned().unwrap_or(ClockEntry {
        delay_s: 0.0,
        rate_sps: 0.0,
        accel_sps2: 0.0,
        jerk_sps3: 0.0,
        snap_sps4: 0.0,
        epoch: None,
    });
    let delay1 = clock1.delay_s;
    let rate1 = clock1.rate_sps;
    let accel1 = clock1.accel_sps2;
    let jerk1 = clock1.jerk_sps3;
    let snap1 = clock1.snap_sps4;
    let delay2 = clock2.delay_s;
    let rate2 = clock2.rate_sps;
    let accel2 = clock2.accel_sps2;
    let jerk2 = clock2.jerk_sps3;
    let snap2 = clock2.snap_sps4;
    let clock_delay_s = Some(delay2 - delay1);
    let clock_rate_sps = Some(rate2 - rate1);
    let clock_accel_sps2 = Some(accel2 - accel1);
    let clock_jerk_sps3 = Some(jerk2 - jerk1);
    let clock_snap_sps4 = Some(snap2 - snap1);

    let source_name = selected_process.object.clone();
    let source_node = if let Some(name) = source_name.as_deref() {
        doc.descendants()
            .find(|n| is_tag(*n, "source") && n.attribute("name").map(str::trim) == Some(name))
            .or_else(|| doc.descendants().find(|n| is_tag(*n, "source")))
    } else {
        doc.descendants().find(|n| is_tag(*n, "source"))
    }
    .ok_or("source node not found in XML schedule")?;
    let ra = selected_process
        .ra
        .clone()
        .or_else(|| child_text(source_node, "ra").map(|s| s.to_string()))
        .ok_or("source/ra missing")?;
    let dec = selected_process
        .dec
        .clone()
        .or_else(|| child_text(source_node, "dec").map(|s| s.to_string()))
        .ok_or("source/dec missing")?;

    let process_epochs = processes
        .iter()
        .map(|p| p.epoch.clone())
        .collect::<Vec<_>>();
    let process_skip_sec = Some(selected_process.skip_sec);
    let process_length_sec = selected_process.length_sec;

    let bw_mhz = sampling_hz.map(|fs| fs / 2e6);
    let ant1_bw_mhz = bw_mhz;
    let ant2_bw_mhz = bw_mhz;
    let ant1_center_mhz = if let (Some(obs), Some(bw)) = (obsfreq_mhz, bw_mhz) {
        Some(obs + ant1_rotation_hz.unwrap_or(0.0) / 1e6 + 0.5 * bw)
    } else {
        None
    };
    let ant2_center_mhz = if let (Some(obs), Some(bw)) = (obsfreq_mhz, bw_mhz) {
        Some(obs + ant2_rotation_hz.unwrap_or(0.0) / 1e6 + 0.5 * bw)
    } else {
        None
    };

    Ok(IFileData {
        ra,
        dec,
        epoch: Some(selected_process.epoch.clone()),
        source: source_name,
        stream_label,
        fft,
        output_sec: output_hz,
        inband,
        pulsar,
        sampling_hz,
        ant1_bit,
        ant2_bit,
        ant1_level,
        ant2_level,
        ant1_shuffle,
        ant2_shuffle,
        obsfreq_mhz,
        clock_delay_s,
        clock_rate_sps,
        clock_accel_sps2,
        clock_jerk_sps3,
        clock_snap_sps4,
        ant1_clock_delay_s: Some(delay1),
        ant2_clock_delay_s: Some(delay2),
        ant1_clock_rate_sps: Some(rate1),
        ant2_clock_rate_sps: Some(rate2),
        ant1_clock_accel_sps2: Some(accel1),
        ant2_clock_accel_sps2: Some(accel2),
        ant1_clock_jerk_sps3: Some(jerk1),
        ant2_clock_jerk_sps3: Some(jerk2),
        ant1_clock_snap_sps4: Some(snap1),
        ant2_clock_snap_sps4: Some(snap2),
        ant1_clock_epoch: clock1.epoch,
        ant2_clock_epoch: clock2.epoch,
        ant1_sideband,
        ant2_sideband,
        ant1_rotation_hz,
        ant2_rotation_hz,
        ant1_rotation2_hz,
        ant2_rotation2_hz,
        ant1_center_mhz,
        ant2_center_mhz,
        ant1_bw_mhz,
        ant2_bw_mhz,
        ant1_station_name,
        ant2_station_name,
        ant1_station_key: Some(ant1_key),
        ant2_station_key: Some(ant2_key),
        ant1_ecef_m,
        ant2_ecef_m,
        process_epochs,
        process_skip_sec,
        process_length_sec,
        processes,
        model_time_offset_s,
        dut1_s,
        tt_utc_s,
        xp_arcsec,
        yp_arcsec,
        eop_file,
    })
}

#[derive(Clone, Debug, PartialEq)]
pub struct ClockPolynomial {
    pub epoch: String,
    pub delay_s: f64,
    pub rate_sps: f64,
    pub accel_sps2: f64,
    pub jerk_sps3: f64,
    pub snap_sps4: f64,
}

pub fn write_schedule_with_clock_polynomial(
    input: &Path,
    output: &Path,
    station_key: &str,
    clock: &ClockPolynomial,
) -> Result<(), DynError> {
    if input == output
        || (output.exists() && std::fs::canonicalize(input)? == std::fs::canonicalize(output)?)
    {
        return Err("phase-calibrated schedule output must differ from input XML".into());
    }
    if station_key.is_empty()
        || !station_key
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "_-.".contains(c))
    {
        return Err(format!("invalid XML station key for clock update: {}", station_key).into());
    }
    for (name, value) in [
        ("delay", clock.delay_s),
        ("rate", clock.rate_sps),
        ("accel", clock.accel_sps2),
        ("jerk", clock.jerk_sps3),
        ("snap", clock.snap_sps4),
    ] {
        if !value.is_finite() {
            return Err(format!("non-finite phase-calibrated clock {}", name).into());
        }
    }

    let text = std::fs::read_to_string(input)?;
    let range = {
        let doc = Document::parse(&text)?;
        let node = doc
            .descendants()
            .find(|n| is_tag(*n, "clock") && n.attribute("key").map(str::trim) == Some(station_key))
            .ok_or_else(|| format!("clock key {} not found in schedule", station_key))?;
        node.range()
    };
    let line_start = text[..range.start].rfind("\n").map_or(0, |i| i + 1);
    let indent = &text[line_start..range.start];
    if !indent.chars().all(char::is_whitespace) {
        return Err("could not determine XML indentation for clock update".into());
    }
    let child_indent = format!("{}  ", indent);
    let replacement = format!(
        "<clock key=\"{}\">\n{}<epoch>{}</epoch>\n{}<delay>{:+.17e}</delay>\n{}<rate>{:+.17e}</rate>\n{}<accel>{:+.17e}</accel>\n{}<jerk>{:+.17e}</jerk>\n{}<snap>{:+.17e}</snap>\n{}</clock>",
        station_key,
        child_indent,
        clock.epoch,
        child_indent,
        clock.delay_s,
        child_indent,
        clock.rate_sps,
        child_indent,
        clock.accel_sps2,
        child_indent,
        clock.jerk_sps3,
        child_indent,
        clock.snap_sps4,
        indent
    );
    let mut updated = String::with_capacity(text.len() + replacement.len());
    updated.push_str(&text[..range.start]);
    updated.push_str(&replacement);
    updated.push_str(&text[range.end..]);
    Document::parse(&updated)
        .map_err(|e| format!("generated phase-calibrated XML is invalid: {}", e))?;

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let partial = PathBuf::from(format!("{}.part", output.display()));
    std::fs::write(&partial, updated)?;
    std::fs::rename(&partial, output)?;
    Ok(())
}

pub fn write_example_xml(path: &PathBuf) -> Result<(), DynError> {
    let template = r#"<?xml version="1.0" encoding="UTF-8"?>
<!-- phased_array example schedule -->
<!-- This file has two parts: -->
<!-- 1) Copy/Paste examples (commented). -->
<!-- 2) Minimal editable template (active XML). -->
<schedule>
  <!-- ==================== COPY/PASTE EXAMPLES (COMMENTED) ==================== -->
  <!-- Station examples -->
  <!--
  <!-- Any single-char key is allowed -->
  <station key='O'><name>KASHIM34</name><pos-x>-3997650.05799</pos-x><pos-y>+3276690.07124</pos-y><pos-z>+3724278.43114</pos-z><terminal>ADS3000</terminal></station>
  <station key='H'><name>HITACH32</name><pos-x>-3961788.9740</pos-x><pos-y>+3243597.4920</pos-y><pos-z>+3790597.6920</pos-z><terminal>ADS3000</terminal></station>
  <station key='T'><name>TAKAHA32</name><pos-x>-3961881.8250</pos-x><pos-y>+3243372.4800</pos-y><pos-z>+3790687.4490</pos-z><terminal>ADS3000</terminal></station>
  <station key='K'><name>YAMAGU32</name><pos-x>-3502544.587</pos-x><pos-y>+3950966.235</pos-y><pos-z>+3566381.192</pos-z><terminal>VSREC</terminal></station>
  <station key='L'><name>YAMAGU34</name><pos-x>-3502567.576</pos-x><pos-y>+3950885.734</pos-y><pos-z>+3566449.115</pos-z><terminal>VSREC</terminal></station>
  -->
  <!-- Terminal examples -->
  <!--
  <!-- ADS3000/ADS1000 use this level order -->
  <terminal name='ADS1000'><speed>1024000000</speed><channel>1</channel><bit>2</bit><level>-1.5,+0.5,-0.5,+1.5</level></terminal>
  <terminal name='ADS3000'><speed>1024000000</speed><channel>1</channel><bit>2</bit><level>-1.5,+0.5,-0.5,+1.5</level></terminal>
  <!-- VSREC uses this level order -->
  <terminal name='VSREC'><speed>1024000000</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,+0.5,+1.5</level></terminal>
  -->
  <!-- VSREC bit shuffle example -->
  <!-- <shuffle key='K'>24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle> -->

  <!-- ==================== MINIMAL EDITABLE TEMPLATE (ACTIVE) ==================== -->
  <!-- EDIT: station key/name/xyz/terminal -->
  <station key="A">
    <name>ANT1</name>
    <pos-x>-3502544.587</pos-x>
    <pos-y>+3950966.235</pos-y>
    <pos-z>+3566381.192</pos-z>
    <terminal>term1</terminal>
  </station>
  <station key="B">
    <name>ANT2</name>
    <pos-x>-3502567.576</pos-x>
    <pos-y>+3950885.734</pos-y>
    <pos-z>+3566449.115</pos-z>
    <terminal>term2</terminal>
  </station>

  <!-- EDIT: clock delay/rate [s], [s/s] -->
  <clock key="A"><delay>0</delay><rate>0</rate></clock>
  <clock key="B"><delay>0</delay><rate>0</rate></clock>

  <!-- EDIT: Earth orientation and delay-model reference -->
  <!-- Use either explicit values or a local IERS finals2000A.data file. -->
  <eop file="data/eop/finals2000A.data"/>
  <delay-model><time-offset>0.0</time-offset></delay-model>

  <!-- EDIT: terminal settings -->
  <terminal name="term1"><speed>1024000000</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,+0.5,+1.5</level></terminal>
  <terminal name="term2"><speed>1024000000</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,+0.5,+1.5</level></terminal>

  <!-- EDIT: bit shuffle -->
  <shuffle key="A">31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0</shuffle>
  <shuffle key="B">31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0</shuffle>

  <!-- EDIT: source -->
  <source name="TARGET"><ra>00h00m00.00000</ra><dec>+00d00'00.0000</dec></source>

  <!-- EDIT: stream -->
  <stream>
    <frequency>0</frequency>
    <fft>16384</fft>
    <output>1.0</output>
    <inband>1</inband>
    <special key="A"><rotation>0</rotation></special>
    <special key="B"><rotation>0</rotation></special>
  </stream>

  <!-- EDIT: process -->
  <process><epoch>2000/001 00:00:00</epoch><length>1</length><object>TARGET</object><stations>AB</stations><baseline>AB</baseline></process>
</schedule>
"#;
    std::fs::write(path, template)?;
    Ok(())
}
