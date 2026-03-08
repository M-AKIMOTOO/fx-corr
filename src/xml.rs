use std::collections::HashMap;
use std::path::PathBuf;

use roxmltree::{Document, Node};

use crate::ifile::{IFileData, ProcessEntry};
use crate::utils::DynError;

#[derive(Clone, Debug)]
struct ClockEntry {
    delay_s: f64,
    rate_sps: f64,
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

fn parse_opt_f64(node: Node<'_, '_>, tag: &str) -> Result<Option<f64>, DynError> {
    Ok(match child_text(node, tag) {
        Some(v) if !v.is_empty() => Some(v.parse::<f64>()?),
        _ => None,
    })
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

fn pick_ant_keys<'a>(
    process: Node<'a, 'a>,
    station_by_key: &HashMap<String, Node<'a, 'a>>,
) -> Result<(String, String), DynError> {
    if let Some(stations) = child_text(process, "stations") {
        let keys: Vec<String> = stations
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_string())
            .collect();
        if keys.len() >= 2 {
            return Ok((keys[0].clone(), keys[1].clone()));
        }
    }
    let mut keys: Vec<String> = station_by_key.keys().cloned().collect();
    keys.sort();
    if keys.len() >= 2 {
        return Ok((keys[0].clone(), keys[1].clone()));
    }
    Err("failed to determine ant1/ant2 station keys from XML".into())
}

pub fn parse_xml_schedule(path: &PathBuf) -> Result<IFileData, DynError> {
    let xml = std::fs::read_to_string(path)?;
    let doc = Document::parse(&xml)?;

    let process_nodes: Vec<_> = doc.descendants().filter(|n| is_tag(*n, "process")).collect();
    let process = *process_nodes
        .first()
        .ok_or("process node not found in XML schedule")?;

    let process_epochs = process_nodes
        .iter()
        .filter_map(|p| child_text(*p, "epoch"))
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

    let source_name = child_text(process, "object").map(|s| s.trim().to_string());
    let epoch = child_text(process, "epoch").map(|s| s.trim().to_string());
    let process_skip_sec = child_text(process, "skip")
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>())
        .transpose()?;
    let process_length_sec = child_text(process, "length")
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<f64>())
        .transpose()?;

    let station_by_key: HashMap<String, Node<'_, '_>> = doc
        .descendants()
        .filter(|n| is_tag(*n, "station"))
        .filter_map(|n| n.attribute("key").map(|k| (k.trim().to_string(), n)))
        .collect();

    let (ant1_key, ant2_key) = pick_ant_keys(process, &station_by_key)?;
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

    let terminal_by_name: HashMap<String, Node<'_, '_>> = doc
        .descendants()
        .filter(|n| is_tag(*n, "terminal"))
        .filter_map(|n| n.attribute("name").map(|name| (name.trim().to_string(), n)))
        .collect();

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

    let ant1_level = term1.and_then(|t| child_text(t, "level")).map(|s| s.to_string());
    let ant2_level = term2.and_then(|t| child_text(t, "level")).map(|s| s.to_string());

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

    let default_sideband = if frequency_hz < 0.0 { "LSB" } else { "USB" };
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
    let ant1_rotation_hz = sp1.map(|n| parse_opt_f64(n, "rotation")).transpose()?.flatten();
    let ant2_rotation_hz = sp2.map(|n| parse_opt_f64(n, "rotation")).transpose()?.flatten();

    let clock_by_key: HashMap<String, ClockEntry> = doc
        .descendants()
        .filter(|n| is_tag(*n, "clock"))
        .filter_map(|n| {
            let key = n.attribute("key")?.trim().to_string();
            let delay = child_text(n, "delay")
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(0.0);
            let rate = child_text(n, "rate")
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(0.0);
            let epoch = child_text(n, "epoch").map(|s| s.to_string());
            Some((
                key,
                ClockEntry {
                    delay_s: delay,
                    rate_sps: rate,
                    epoch,
                },
            ))
        })
        .collect();
    let clock1 = clock_by_key.get(&ant1_key).cloned().unwrap_or(ClockEntry {
        delay_s: 0.0,
        rate_sps: 0.0,
        epoch: None,
    });
    let clock2 = clock_by_key.get(&ant2_key).cloned().unwrap_or(ClockEntry {
        delay_s: 0.0,
        rate_sps: 0.0,
        epoch: None,
    });
    let delay1 = clock1.delay_s;
    let rate1 = clock1.rate_sps;
    let delay2 = clock2.delay_s;
    let rate2 = clock2.rate_sps;
    let clock_delay_s = Some(delay2 - delay1);
    let clock_rate_sps = Some(rate2 - rate1);

    let source_node = if let Some(name) = source_name.as_deref() {
        doc.descendants()
            .find(|n| is_tag(*n, "source") && n.attribute("name").map(str::trim) == Some(name))
            .or_else(|| doc.descendants().find(|n| is_tag(*n, "source")))
    } else {
        doc.descendants().find(|n| is_tag(*n, "source"))
    }
    .ok_or("source node not found in XML schedule")?;
    let ra = child_text(source_node, "ra").ok_or("source/ra missing")?.to_string();
    let dec = child_text(source_node, "dec").ok_or("source/dec missing")?.to_string();

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

    let processes: Vec<ProcessEntry> = process_nodes
        .iter()
        .filter_map(|p| {
            let epoch = child_text(*p, "epoch")
                .map(|s| s.trim().to_string())
                .unwrap_or_else(|| "2000".to_string());
            let skip_sec = child_text(*p, "skip")
                .and_then(|s| s.trim().parse::<f64>().ok())
                .unwrap_or(0.0);
            let length_sec = child_text(*p, "length")
                .and_then(|s| s.trim().parse::<f64>().ok());
            let object = child_text(*p, "object").map(|s| s.trim().to_string());
            let stations = child_text(*p, "stations").map(|s| s.trim().to_string());
            if let Some(st) = stations.as_ref() {
                let stn = st.replace(char::is_whitespace, "");
                if !stn.contains(&ant1_key) || !stn.contains(&ant2_key) {
                    return None;
                }
            }
            let (ra_opt, dec_opt) = if let Some(obj) = object.as_ref() {
                if let Some((ra, dec)) = source_radec_by_name.get(obj) {
                    (Some(ra.clone()), Some(dec.clone()))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };
            Some(ProcessEntry {
                epoch,
                skip_sec,
                length_sec,
                object,
                stations,
                ra: ra_opt,
                dec: dec_opt,
            })
        })
        .collect();

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
        epoch,
        source: source_name,
        stream_label,
        fft,
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
        ant1_clock_delay_s: Some(delay1),
        ant2_clock_delay_s: Some(delay2),
        ant1_clock_rate_sps: Some(rate1),
        ant2_clock_rate_sps: Some(rate2),
        ant1_clock_epoch: clock1.epoch,
        ant2_clock_epoch: clock2.epoch,
        ant1_sideband,
        ant2_sideband,
        ant1_rotation_hz,
        ant2_rotation_hz,
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
    })
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
    <special key="A"><rotation>0</rotation><sideband>LSB</sideband></special>
    <special key="B"><rotation>0</rotation><sideband>LSB</sideband></special>
  </stream>

  <!-- EDIT: process -->
  <process><epoch>2000/001 00:00:00</epoch><length>1</length><object>TARGET</object><stations>AB</stations></process>
</schedule>
"#;
    std::fs::write(path, template)?;
    Ok(())
}
