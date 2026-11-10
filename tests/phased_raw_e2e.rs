use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_temp_dir() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("yi-phasedarray-e2e-{}-{nonce}", std::process::id()))
}

fn write_schedule_with_clock_delay(
    path: &Path,
    station1: &str,
    station2: &str,
    station2_clock_delay_s: f64,
) {
    let xml = format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<schedule>
  <station key="A"><name>{station1}</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z><terminal>term</terminal></station>
  <station key="B"><name>{station2}</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z><terminal>term</terminal></station>
  <clock key="A"><delay>0</delay><rate>0</rate></clock>
  <clock key="B"><delay>{station2_clock_delay_s}</delay><rate>0</rate></clock>
  <terminal name="term"><speed>8192</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,0.5,1.5</level></terminal>
  <shuffle key="A">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <shuffle key="B">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <source name="TARGET"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>
  <stream>
    <frequency>100000000</frequency><fft>128</fft><output>8</output>
    <special key="A"><rotation>0</rotation><sideband>LSB</sideband></special>
    <special key="B"><rotation>0</rotation><sideband>LSB</sideband></special>
  </stream>
  <process><epoch>2000/001 00:00:00</epoch><length>1</length><object>TARGET</object><stations>AB</stations><baseline>AB</baseline></process>
</schedule>
"#
    );
    fs::write(path, xml).unwrap();
}

fn write_schedule(path: &Path, station1: &str, station2: &str) {
    write_schedule_with_clock_delay(path, station1, station2, 0.0);
}

#[test]
fn ionosphere_corrects_injected_dispersive_voltage_without_changing_frequency_bins() {
    use std::f64::consts::TAU;
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    fs::create_dir_all(&raw_dir).unwrap();
    let schedule = root.join("test.xml");
    let nfft = 1024;
    let fs_hz = 32_768_000.0;
    let rf_hz = 40_000_000.0;
    let tones = [20usize, 73, 147, 212, 289, 371, 443];
    // XML shuffle lists output bits MSB first; 31..0 is byte/sample identity.
    let identity = (0..32)
        .rev()
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let levels = (0..256)
        .map(|i| (i as f64 - 127.5).to_string())
        .collect::<Vec<_>>()
        .join(",");
    for (case, sideband, delta_tec, rot, tec_rate) in [
        (0, "LSB", 1.0, 0.0, 0.0),
        (1, "USB", -1.0, 0.0, 0.0),
        (2, "USB", 1.0, 1_024_000.0, 0.0),
        (3, "LSB", 1.0, 0.0, 20.0),
    ] {
        write_schedule(&schedule, "ANT1", "ANT2");
        let xml = fs::read_to_string(&schedule).unwrap()
            .replace("<speed>8192</speed>", "<speed>32768000</speed>")
            .replace("<bit>2</bit>", "<bit>8</bit>")
            .replace("<level>-1.5,-0.5,0.5,1.5</level>", &format!("<level>{levels}</level>"))
            .replace("24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7", &identity)
            .replace("<frequency>100000000</frequency>", "<frequency>40000000</frequency>")
            .replace("<fft>128</fft>", "<fft>1024</fft>")
            .replace("<output>8</output>", "<output>250</output>")
            .replace("<length>1</length>", "<length>0.004</length>")
            .replace("<sideband>LSB</sideband>", &format!("<sideband>{sideband}</sideband>"))
            .replace("<special key=\"B\"><rotation>0</rotation>", &format!("<special key=\"B\"><rotation>{rot}</rotation>"))
            .replace("<name>ANT1</name>", &format!("<name>ANT1</name><atmosphere stec-tecu=\"{}\"/>", if delta_tec > 0.0 { 0 } else { 1 }))
            .replace("<name>ANT2</name>", &format!("<name>ANT2</name><atmosphere stec-tecu=\"{}\"/>", if delta_tec > 0.0 { 1 } else { 0 }));
        let xml = if tec_rate != 0.0 {
            let csv = root.join("tec.csv");
            let last = 51544.0 + 0.004 / 86400.0;
            fs::write(csv, format!("station,mjd_utc,stec_tecu\nA,51544,0\nB,51544,1\nA,{last:.12},0\nB,{last:.12},1.08\n")).unwrap();
            xml.replace("<schedule>", "<schedule><atmosphere file=\"tec.csv\"/>")
        } else {
            xml
        };
        fs::write(&schedule, &xml).unwrap();
        // Independent sky/voltage model. A positive TEC causes phase advance
        // +2*pi*40.3*TEC/(c*f) in station 2, with frequency in physical RF.
        let d = 40.3e16 * delta_tec / 299_792_458.0;
        for (station, ionized, rotation) in [("ANT1", false, 0.0), ("ANT2", true, rot)] {
            let sample_count = if tec_rate == 0.0 { nfft } else { nfft * 128 };
            let frame: Vec<u8> = (0..sample_count)
                .map(|sample| {
                    let frame_mid_s =
                        (sample / nfft) as f64 * nfft as f64 / fs_hz + 0.5 * nfft as f64 / fs_hz;
                    let d = d + 40.3e16 / 299_792_458.0 * tec_rate * frame_mid_s;
                    let mut value = 0.0;
                    for (j, &bin) in tones.iter().enumerate() {
                        let f = rf_hz + bin as f64 * fs_hz / nfft as f64;
                        let native_bin = bin as f64 - rotation * nfft as f64 / fs_hz;
                        if native_bin <= 0.0 {
                            continue;
                        }
                        let phase = TAU * native_bin * sample as f64 / nfft as f64
                            + j as f64 * 0.271
                            + if ionized { TAU * d / f } else { 0.0 };
                        value += 8.0 * phase.cos();
                    }
                    if sideband == "LSB" && sample % 2 == 1 {
                        value = -value;
                    }
                    (value + 127.5).round().clamp(0.0, 255.0) as u8
                })
                .collect();
            let raw = frame.repeat(if tec_rate == 0.0 { 128 } else { 1 });
            fs::write(raw_dir.join(format!("{station}_2000001000000.raw")), raw).unwrap();
        }
        let mut corrected = Vec::new();
        for mode in ["off", "on", "debug", "fold", "phased"] {
            fs::write(
                &schedule,
                if mode == "fold" {
                    xml.replace(
                        "<schedule>",
                        "<schedule><pulsar><period>1</period><bins>1</bins><dm>0</dm></pulsar>",
                    )
                } else {
                    xml.clone()
                },
            )
            .unwrap();
            let phased = mode == "phased";
            let out = root.join(format!("{case}-{mode}"));
            let mut cmd = Command::new(if phased {
                env!("CARGO_BIN_EXE_yi-phasedarray")
            } else {
                env!("CARGO_BIN_EXE_yi-corr")
            });
            cmd.args([
                "--sc",
                schedule.to_str().unwrap(),
                "--raw",
                raw_dir.to_str().unwrap(),
                if phased { "--output" } else { "--cor" },
                out.to_str().unwrap(),
                "--cpu",
                "1",
                "--no-affinity",
                "--ionosphere",
                if mode == "off" { "off" } else { "input" },
            ])
            .env("YI_EOP_MODE", "none")
            .env("YI_ANT2_GRID_OFFSET", "0")
            .env("YI_FRINGE_STOP_FREQ_MODE", "low")
            .env("YI_GEOM_DELAY_MODE", "anchored");
            if mode == "debug" {
                cmd.arg("--debug");
            }
            if phased {
                cmd.args(["--phased-name", "ARRAY", "--phased-diagnostics"]);
            }
            let result = cmd.output().unwrap();
            assert!(
                result.status.success(),
                "case={case} {mode}: {}\n{}",
                String::from_utf8_lossy(&result.stdout),
                String::from_utf8_lossy(&result.stderr)
            );
            let product = fs::read_dir(&out)
                .unwrap()
                .map(|e| e.unwrap().path())
                .find(|p| {
                    p.file_name()
                        .unwrap()
                        .to_string_lossy()
                        .starts_with("ANT1_ANT2_")
                        && p.extension().and_then(|s| s.to_str()) == Some("cor")
                })
                .unwrap();
            let bytes = fs::read(product).unwrap();
            let mut phases = Vec::new();
            for bin in tones
                .into_iter()
                .filter(|bin| *bin as f64 * fs_hz / nfft as f64 > rot)
            {
                let start = 256 + 128 + bin * 8;
                let re = f32::from_le_bytes(bytes[start..start + 4].try_into().unwrap()) as f64;
                let im = f32::from_le_bytes(bytes[start + 4..start + 8].try_into().unwrap()) as f64;
                assert!(re.hypot(im) > 1e-8, "tone disappeared from bin {bin}");
                phases.push(im.atan2(re));
            }
            if mode == "off" {
                assert!(
                    phases.iter().any(|p| p.abs() > 0.5),
                    "injection must require correction"
                );
            } else {
                assert!(
                    phases.iter().all(|p| p.abs() < 0.015),
                    "case={case} {mode} residual phases {phases:?}"
                );
                if mode == "on" {
                    corrected = phases;
                } else {
                    assert!(
                        phases
                            .iter()
                            .zip(&corrected)
                            .all(|(a, b)| (a - b).abs() < 1e-6),
                        "phased/mapped correction diverged"
                    );
                }
            }
        }
    }
    fs::remove_dir_all(root).unwrap();
}

fn corr_info_paths(directory: &Path) -> Vec<PathBuf> {
    fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().and_then(|ext| ext.to_str()) == Some("info"))
        .collect()
}

#[test]
fn atmosphere_xml_and_relative_csv_inputs_match_and_invalid_coverage_fails() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    fs::create_dir_all(&raw_dir).unwrap();
    let schedule = root.join("test.xml");
    write_schedule(&schedule, "ANT1", "ANT2");
    let base_xml = fs::read_to_string(&schedule)
        .unwrap()
        .replace("00h00m00.0", "15h30m00.0")
        .replace("+00d00'00.0", "+35d00'00.0");
    let raw: Vec<u8> = (0..4096).map(|i| ((i * 53 + i / 7) % 256) as u8).collect();
    for station in ["ANT1", "ANT2"] {
        fs::write(raw_dir.join(format!("{station}_2000001000000.raw")), &raw).unwrap();
    }
    let mut reference = std::collections::BTreeMap::new();
    for provider in [
        "xml",
        "csv-constant",
        "csv-dated",
        "cli-override",
        "off",
        "bare",
        "bad-coverage",
    ] {
        let out = root.join(provider);
        let csv = root.join("weather.csv");
        let header = "station,mjd_utc,vtec_tecu,pressure_hpa,temperature_k,vapour_pressure_hpa\n";
        let mut csv_text = header.to_string();
        if provider == "csv-constant" || provider == "cli-override" {
            csv_text += "A,,10,1000,288.15,10\nB,,10.1,990,289.15,11\n";
        } else {
            let end = 51544.0
                + if provider == "bad-coverage" {
                    0.25
                } else {
                    2.0
                } / 86400.0;
            for t in [51544.0, end] {
                csv_text += &format!("A,{t:.12},10,1000,288.15,10\nB,{t:.12},10.1,990,289.15,11\n");
            }
        }
        fs::write(&csv, csv_text).unwrap();
        let xml = match provider {
            "xml" | "off" => base_xml
                .replace("<name>ANT1</name>", "<name>ANT1</name><atmosphere vtec-tecu=\"10\" pressure-hpa=\"1000\" temperature-k=\"288.15\" vapour-pressure-hpa=\"10\"/>")
                .replace("<name>ANT2</name>", "<name>ANT2</name><atmosphere vtec-tecu=\"10.1\" pressure-hpa=\"990\" temperature-k=\"289.15\" vapour-pressure-hpa=\"11\"/>"),
            "bare" => base_xml.clone(),
            _ => base_xml.replace("<schedule>", &format!("<schedule><atmosphere file=\"{}\"/>", if provider == "cli-override" {"absent.csv"} else {"weather.csv"})),
        };
        fs::write(&schedule, xml).unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_yi-corr"));
        // Run outside the schedule directory: XML file references must remain
        // relative to the original XML, rather than the process cwd.
        command
            .current_dir(&raw_dir)
            .args([
                "--sc",
                schedule.to_str().unwrap(),
                "--raw",
                raw_dir.to_str().unwrap(),
                "--cor",
                out.to_str().unwrap(),
                "--cpu",
                "1",
                "--no-affinity",
                "--skip",
                "0.125",
            ])
            .env("YI_EOP_MODE", "none")
            .env("YI_GEOM_DELAY_MODE", "anchored")
            .env("YI_FRINGE_STOP_FREQ_MODE", "low");
        if provider == "cli-override" {
            command.args(["--atmosphere-file", csv.to_str().unwrap()]);
        }
        if provider == "off" {
            command.args(["--ionosphere", "off", "--troposphere", "off"]);
        }
        let result = command.output().unwrap();
        if provider == "bad-coverage" {
            assert!(!result.status.success());
            assert!(String::from_utf8_lossy(&result.stderr).contains("outside input coverage"));
            assert!(read_corr_info(&corr_info_paths(&out)[0]).ends_with("# status=error\n"));
            continue;
        }
        assert!(
            result.status.success(),
            "{provider}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let products: std::collections::BTreeMap<_, _> = fs::read_dir(&out)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("cor"))
            .map(|p| (p.file_name().unwrap().to_owned(), fs::read(p).unwrap()))
            .collect();
        if provider == "xml" || provider == "off" {
            reference = products;
        } else {
            assert_eq!(products, reference, "{provider}");
        }
    }
    fs::remove_dir_all(root).unwrap();
}

fn read_corr_info(path: &Path) -> String {
    let bytes = fs::read(path).unwrap();
    assert_eq!(&bytes[..4], b"PK\x03\x04", "not ZIP: {}", path.display());
    let mut archive = zip::ZipArchive::new(std::io::Cursor::new(bytes))
        .expect("complete ZIP including the central directory");
    assert_eq!(
        archive.len(),
        1,
        "zipcat must be able to omit the member name"
    );
    let mut member = archive.by_name("info.txt").unwrap();
    assert_eq!(member.compression(), zip::CompressionMethod::Deflated);
    assert!(!member.encrypted());
    let mut text = String::new();
    member
        .read_to_string(&mut text)
        .expect("UTF-8 text and valid ZIP CRC");
    text
}

#[test]
fn corr_saves_compressed_details_and_shows_each_scan_concisely() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let out = root.join("cor");
    fs::create_dir_all(&raw_dir).unwrap();
    let schedule = root.join("test.xml");
    write_schedule(&schedule, "ANT1", "ANT2");
    let xml = fs::read_to_string(&schedule).unwrap()
        .replace("<source name=\"TARGET\">", "<source name=\"CAL\"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>\n<source name=\"TARGET\">")
        .replace("</schedule>", "<process><epoch>2000/001 00:00:10</epoch><length>1</length><object>CAL</object><stations>AB</stations></process>\n</schedule>");
    fs::write(&schedule, xml).unwrap();
    let raw: Vec<u8> = (0..2048).map(|i| ((i * 53 + i / 7) % 256) as u8).collect();
    for station in ["ANT1", "ANT2"] {
        for tag in ["2000001000000", "2000001000010"] {
            fs::write(raw_dir.join(format!("{station}_{tag}.raw")), &raw).unwrap();
        }
    }
    let mut reference_products = std::collections::BTreeMap::new();
    for save_stdout in [false, true] {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_yi-corr"));
        cmd.current_dir(&root)
            .args([
                "--sc",
                schedule.to_str().unwrap(),
                "--raw",
                raw_dir.to_str().unwrap(),
                "--cor",
                out.to_str().unwrap(),
                "--cpu",
                "1",
                "--no-affinity",
                "--fft",
                "256",
                "--model-diagnostics",
            ])
            .env("YI_EOP_MODE", "none");
        if save_stdout {
            cmd.arg("--stdout");
        }
        let result = cmd.output().unwrap();
        let stdout = String::from_utf8_lossy(&result.stdout);
        assert!(
            result.status.success(),
            "{stdout}\n{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(stdout.matches("[scan] source=").count(), 2, "{stdout}");
        assert_eq!(stdout.matches("[done] source=").count(), 2, "{stdout}");
        for expected in [
            "source=TARGET baseline=ANT1-ANT2 epoch=2000/001 00:00:00",
            "source=CAL baseline=ANT1-ANT2 epoch=2000/001 00:00:10",
            "FFT=256 sampling=0.008192MS/s bits=2/2 sideband=LSB/LSB",
            "RF=100.000000MHz BW=0.004096MHz df=32.000000Hz",
            "integration=0.125000000s (requested 8.000000Hz; XML) sectors=8",
            "compute-threads=1 affinity=disabled",
            "100.0%",
        ] {
            assert!(stdout.contains(expected), "missing {expected:?}:\n{stdout}");
        }
        for hidden in ["[info]", "Antenna Parameters", "geom-diag", "[mk_xml-arg]"] {
            assert!(
                !stdout.contains(hidden),
                "detail leaked: {hidden}\n{stdout}"
            );
        }
        assert!(
            !stdout.contains('\r'),
            "redirected output must use complete lines"
        );
        assert!(stdout.lines().count() < 30, "{stdout}");
        let paths = corr_info_paths(&out);
        assert_eq!(
            paths.len(),
            if save_stdout { 2 } else { 1 },
            "logs must not overwrite"
        );
        for path in paths {
            let info = read_corr_info(&path);
            for expected in [
                "# schedule XML:",
                "# resolved CLI arguments:",
                "[info] CPU allocation:",
                "[info] Parsed schedule XML parameters:",
                "[mk_xml-arg] --scan=1",
                "[mk_xml-arg] --scan=2",
                "[geom-diag:source]",
                "# status=ok",
            ] {
                assert!(
                    info.contains(expected),
                    "missing {expected:?} in {}",
                    path.display()
                );
            }
            assert_eq!(info.matches("[info] Synth timing summary:").count(), 2);
        }
        let products: std::collections::BTreeMap<_, _> = fs::read_dir(&out)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().and_then(|ext| ext.to_str()) == Some("cor"))
            .map(|path| {
                (
                    path.file_name().unwrap().to_owned(),
                    fs::read(&path).unwrap(),
                )
            })
            .collect();
        assert_eq!(products.len(), 6);
        if save_stdout {
            assert_eq!(
                products, reference_products,
                "logging must preserve correlation products"
            );
            let logs: Vec<_> = fs::read_dir(root.join("stdout")).unwrap().collect();
            assert_eq!(logs.len(), 1);
            assert_eq!(
                fs::read_to_string(logs[0].as_ref().unwrap().path()).unwrap(),
                stdout
            );
        } else {
            reference_products = products;
            assert!(!root.join("stdout").exists());
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn corr_finalizes_info_even_after_a_processing_error() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let out = root.join("cor");
    fs::create_dir_all(&raw_dir).unwrap();
    let schedule = root.join("test.xml");
    write_schedule(&schedule, "ANT1", "ANT2");
    let xml = fs::read_to_string(&schedule)
        .unwrap()
        .replace("<output>8</output>", "<output>0</output>");
    fs::write(&schedule, xml).unwrap();
    for station in ["ANT1", "ANT2"] {
        fs::write(
            raw_dir.join(format!("{station}_2000001000000.raw")),
            vec![0x39; 2048],
        )
        .unwrap();
    }
    let result = Command::new(env!("CARGO_BIN_EXE_yi-corr"))
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--cor",
            out.to_str().unwrap(),
            "--cpu",
            "1",
            "--no-affinity",
        ])
        .env("YI_EOP_MODE", "none")
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("stream/output must be positive finite")
    );
    let paths = corr_info_paths(&out);
    assert_eq!(paths.len(), 1);
    let info = read_corr_info(&paths[0]);
    assert!(info.contains("[error] stream/output must be positive finite"));
    assert!(info.ends_with("# status=error\n"));
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn multiband_calibrator_corrects_targets_without_fitting_them() {
    let root = unique_temp_dir();
    fs::create_dir_all(&root).unwrap();
    let schedules = [root.join("low.xml"), root.join("high.xml")];
    let raw_dirs = [root.join("c"), root.join("x")];
    for band in 0..2 {
        fs::create_dir_all(&raw_dirs[band]).unwrap();
        write_schedule_with_clock_delay(&schedules[band], "ANT1", "ANT2", 2e-6);
        let xml = fs::read_to_string(&schedules[band]).unwrap()
            .replace("<frequency>100000000</frequency>", &format!("<frequency>{}</frequency>", 100000000 + band * 16384))
            .replace("<output>8</output>", "<output>2</output>")
            .replace("<source name=\"TARGET\">", "<source name=\"CAL\"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>\n<source name=\"TARGET\">")
            .replace("<object>TARGET</object>", "<object>CAL</object>")
            .replace("</schedule>", "<process><epoch>2000/001 00:00:10</epoch><length>1</length><object>TARGET</object><stations>AB</stations><baseline>AB</baseline></process>\n</schedule>");
        fs::write(&schedules[band], xml).unwrap();
        let mut state = 0x1234_5678_u32 + band as u32;
        let mut raw = vec![0_u8; 2048];
        for byte in &mut raw {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            *byte = (state >> 24) as u8;
        }
        for tag in ["2000001000000", "2000001000010"] {
            for station in ["ANT1", "ANT2"] {
                fs::write(raw_dirs[band].join(format!("{station}_{tag}.raw")), &raw).unwrap();
            }
        }
    }
    let low_xml = fs::read_to_string(&schedules[0]).unwrap();
    let high_xml = fs::read_to_string(&schedules[1]).unwrap();
    let low_doc = roxmltree::Document::parse(&low_xml).unwrap();
    let high_doc = roxmltree::Document::parse(&high_xml).unwrap();
    let stream = |doc: &roxmltree::Document<'_>| {
        doc.root_element()
            .children()
            .find(|n| n.has_tag_name("stream"))
            .unwrap()
            .range()
    };
    let low_range = stream(&low_doc);
    let cx_xml = root.join("cx.xml");
    let mut cx = low_xml.clone();
    cx.replace_range(low_range.clone(), &format!(
        "<multiband calibrator=\"CAL\"><band name=\"x\">{}</band><band name=\"c\">{}</band></multiband>",
        &high_xml[stream(&high_doc)], &low_xml[low_range]
    ));
    fs::write(&cx_xml, cx).unwrap();
    let mut reference_joint = Vec::new();
    for mode in ["direct", "mapped", "single-xml"] {
        let out = root.join(mode);
        let mut command = Command::new(env!("CARGO_BIN_EXE_yi-corr"));
        command.args([
            "--sc",
            if mode == "single-xml" {
                &cx_xml
            } else {
                &schedules[0]
            }
            .to_str()
            .unwrap(),
            "--raw",
            if mode == "single-xml" {
                &root
            } else {
                &raw_dirs[0]
            }
            .to_str()
            .unwrap(),
            "--cor",
            out.to_str().unwrap(),
            "--multiband-delay-window-ns",
            "10000",
            "--multiband-solve-integration",
            "0.0625",
            "--multiband-window",
            "1",
            "--cpu",
            "1",
        ]);
        if mode != "single-xml" {
            command.args([
                "--multiband-schedule",
                schedules[1].to_str().unwrap(),
                "--multiband-raw-directory",
                raw_dirs[1].to_str().unwrap(),
                "--multiband-calibrator",
                "CAL",
            ]);
        }
        if mode == "mapped" {
            command.arg("--debug");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "stdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&result.stdout),
            String::from_utf8_lossy(&result.stderr)
        );
        let scan = out.join("multiband/scan0001");
        assert!(
            !scan.join("solve-band1").exists(),
            "target must not be fringe-fitted"
        );
        let joint_path = scan.join("ANT1_ANT2_2000001000010_mbcx.cor");
        let bytes = fs::read(&joint_path).unwrap();
        let mut product_names = fs::read_dir(&scan)
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("cor"))
            .map(|p| p.file_name().unwrap().to_str().unwrap().to_owned())
            .collect::<Vec<_>>();
        product_names.sort();
        assert_eq!(
            product_names,
            [
                "ANT1_ANT1_2000001000010_mbcx.cor",
                "ANT1_ANT2_2000001000010_mbcx.cor",
                "ANT2_ANT2_2000001000010_mbcx.cor",
            ]
        );
        assert!(out
            .join("multiband/scan0000/ANT1_ANT2_2000001000000_mbcx.cor")
            .is_file());
        if mode == "direct" {
            reference_joint = bytes.clone();
        } else {
            assert_eq!(bytes, reference_joint, "{mode}");
        }
        assert_eq!(&bytes[..8], b"YIMBCOR\0");
        assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()), 1);
        assert_eq!(
            f64::from_le_bytes(bytes[16..24].try_into().unwrap()),
            8192.0
        );
        assert_eq!(
            f64::from_le_bytes(bytes[48..56].try_into().unwrap()),
            100000000.0
        );
        assert_eq!(
            f64::from_le_bytes(bytes[304..312].try_into().unwrap()),
            100016384.0
        );
        assert_eq!(bytes.len(), 544 + 2 * 2 * (128 + 64 * 8));
        for station in [1, 2] {
            let joint =
                fs::read(scan.join(format!("ANT{station}_ANT{station}_2000001000010_mbcx.cor")))
                    .unwrap();
            assert_eq!(joint.len(), bytes.len());
            for band in 0..2 {
                let h = &joint[32 + band * 256..288 + band * 256];
                assert_eq!(&h[32..48], &h[80..96]);
                let native = fs::read_dir(scan.join(format!("band{}", band + 1)))
                    .unwrap()
                    .map(|e| fs::read(e.unwrap().path()).unwrap())
                    .find(|b| b.len() >= 256 && &b[..256] == h)
                    .unwrap();
                for row in 0..2 {
                    let begin = 544 + (row * 2 + band) * (128 + 64 * 8);
                    assert_eq!(
                        &joint[begin..begin + 128 + 64 * 8],
                        &native[256 + row * (128 + 64 * 8)..256 + (row + 1) * (128 + 64 * 8)]
                    );
                    for v in joint[begin + 128..begin + 128 + 64 * 8].chunks_exact(8) {
                        assert!(f32::from_le_bytes(v[..4].try_into().unwrap()) >= 0.0);
                        assert_eq!(f32::from_le_bytes(v[4..].try_into().unwrap()), 0.0);
                    }
                }
            }
        }
        for graph in [
            "solutions.png",
            "applied.png",
            "visibility.png",
            "band-phase.png",
            "spectrum.png",
            "autocorrelation.png",
        ] {
            let png = fs::read(scan.join(graph)).unwrap();
            assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        }
        assert!(fs::read_to_string(scan.join("verification.txt"))
            .unwrap()
            .contains("mode=calibrator_transfer"));
        let applied = fs::read_to_string(scan.join("applied.tsv")).unwrap();
        let times = applied
            .lines()
            .skip(1)
            .map(|s| s.split('\t').next().unwrap().parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(
            times,
            vec![0.25, 0.75],
            "target output times must use the target epoch"
        );
        assert!(out
            .join("multiband/scan0000/uncorrected-spectrum.tsv")
            .exists());
        assert!(!scan.join("uncorrected-spectrum.tsv").exists());
        for row in 0..2 {
            for band in 0..2 {
                let begin = 544 + (row * 2 + band) * (128 + 64 * 8) + 128;
                let sum = bytes[begin..begin + 64 * 8].chunks_exact(8).fold(
                    (0.0_f64, 0.0_f64),
                    |(r, i), z| {
                        (
                            r + f32::from_le_bytes(z[..4].try_into().unwrap()) as f64,
                            i + f32::from_le_bytes(z[4..].try_into().unwrap()) as f64,
                        )
                    },
                );
                assert!(sum.0 > 0.0 && sum.1.abs() < sum.0 * 0.001, "{sum:?}");
            }
        }
        let exported = Command::new("python3")
            .args([
                concat!(env!("CARGO_MANIFEST_DIR"), "/tools/multiband_visibility.py"),
                joint_path.to_str().unwrap(),
                "--average",
                scan.join("average.tsv").to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            exported.status.success(),
            "{}",
            String::from_utf8_lossy(&exported.stderr)
        );
        assert_eq!(
            fs::read_to_string(scan.join("average.tsv"))
                .unwrap()
                .lines()
                .count(),
            3
        );
        let time_mean = Command::new("python3")
            .args([
                concat!(env!("CARGO_MANIFEST_DIR"), "/tools/multiband_visibility.py"),
                joint_path.to_str().unwrap(),
                "--average",
                scan.join("time-mean.tsv").to_str().unwrap(),
                "--time-average",
                "--band-scales",
                "2",
                "1",
                "--band-weights",
                "1",
                "3",
            ])
            .output()
            .unwrap();
        assert!(
            time_mean.status.success(),
            "{}",
            String::from_utf8_lossy(&time_mean.stderr)
        );
        let means = fs::read_to_string(scan.join("average.tsv")).unwrap();
        let qa = fs::read_to_string(scan.join("visibility-time.tsv")).unwrap();
        for (export, diagnostic) in means.lines().skip(1).zip(qa.lines().skip(1)) {
            let a: Vec<f64> = export.split('\t').map(|v| v.parse().unwrap()).collect();
            let b: Vec<f64> = diagnostic.split('\t').map(|v| v.parse().unwrap()).collect();
            assert!((a[3] - b[11]).abs() < 1e-9);
            assert!((a[4] - b[12]).abs() < 1e-9);
            assert_eq!(a[10], b[17]);
        }
        let expected = means
            .lines()
            .skip(1)
            .map(|line| {
                let values: Vec<f64> = line.split('\t').map(|v| v.parse().unwrap()).collect();
                (2.0 * values[6] + 3.0 * values[8]) / 4.0
            })
            .sum::<f64>()
            / 2.0;
        let time_means = fs::read_to_string(scan.join("time-mean.tsv")).unwrap();
        let values: Vec<f64> = time_means
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .map(|v| v.parse().unwrap())
            .collect();
        assert_eq!(time_means.lines().count(), 2);
        assert!((values[1] - 1.0).abs() < 1e-6);
        assert!((values[3] - expected).abs() < 1e-9);
    }
    fs::remove_dir_all(root).unwrap();
}

fn write_three_station_schedule(path: &Path) {
    let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<schedule>
  <station key="A"><name>YAMAGU32</name><pos-x>0</pos-x><pos-y>0</pos-y><pos-z>0</pos-z><terminal>term</terminal></station>
  <station key="B"><name>YAMAGU34</name><pos-x>110</pos-x><pos-y>0</pos-y><pos-z>0</pos-z><terminal>term</terminal></station>
  <station key="C"><name>HITACH32</name><pos-x>873000</pos-x><pos-y>0</pos-y><pos-z>0</pos-z><terminal>term</terminal></station>
  <clock key="A"><delay>0</delay><rate>0</rate></clock>
  <clock key="B"><delay>0</delay><rate>0</rate></clock>
  <clock key="C"><delay>0</delay><rate>0</rate></clock>
  <terminal name="term"><speed>8192</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,0.5,1.5</level></terminal>
  <shuffle key="A">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <shuffle key="B">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <shuffle key="C">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <source name="TARGET"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>
  <stream>
    <label>validation</label><frequency>100000000</frequency><fft>128</fft><output>8</output>
    <special key="A"><rotation>0</rotation><sideband>LSB</sideband></special>
    <special key="B"><rotation>0</rotation><sideband>LSB</sideband></special>
    <special key="C"><rotation>0</rotation><sideband>LSB</sideband></special>
  </stream>
  <process><epoch>2000/001 00:00:00</epoch><length>1</length><object>TARGET</object><stations>ABC</stations><baseline>AB AC BC</baseline></process>
</schedule>
"#;
    fs::write(path, xml).unwrap();
}

fn write_phasecal_schedule(path: &Path) {
    let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<schedule>
  <station key="A"><name>YAMAGU32</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z><terminal>term</terminal></station>
  <station key="B"><name>YAMAGU34</name><pos-x>-3502567.576</pos-x><pos-y>3950885.734</pos-y><pos-z>3566449.115</pos-z><terminal>term</terminal></station>
  <clock key="A"><epoch>2000/001 00:00:00</epoch><delay>0</delay><rate>0</rate><accel>0</accel></clock>
  <clock key="B"><epoch>2000/001 00:00:00</epoch><delay>1e-6</delay><rate>2e-12</rate><accel>3e-18</accel></clock>
  <terminal name="term"><speed>8192</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,0.5,1.5</level></terminal>
  <shuffle key="A">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <shuffle key="B">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
  <source name="GAIN"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>
  <source name="TARGET"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>
  <stream>
    <label>phasecal</label><frequency>100000000</frequency><fft>128</fft><output>8</output>
    <special key="A"><rotation>0</rotation><sideband>LSB</sideband></special>
    <special key="B"><rotation>0</rotation><sideband>LSB</sideband></special>
  </stream>
  <process><epoch>2000/001 00:00:00</epoch><length>1</length><object>GAIN</object><stations>AB</stations><baseline>AB</baseline></process>
  <process><epoch>2000/001 00:00:10</epoch><length>1</length><object>TARGET</object><stations>AB</stations><baseline>AB</baseline></process>
  <process><epoch>2000/001 00:00:20</epoch><length>1</length><object>GAIN</object><stations>AB</stations><baseline>AB</baseline></process>
  <process><epoch>2000/001 00:00:30</epoch><length>1</length><object>TARGET</object><stations>AB</stations><baseline>AB</baseline></process>
  <process><epoch>2000/001 00:00:40</epoch><length>1</length><object>GAIN</object><stations>AB</stations><baseline>AB</baseline></process>
</schedule>
"#;
    fs::write(path, xml).unwrap();
}

#[cfg(unix)]
fn write_fake_frinz(path: &Path) {
    use std::os::unix::fs::PermissionsExt;

    let script = r##"#!/bin/sh
input=
mode=
loops=1
while [ "$#" -gt 0 ]; do
  case "$1" in
    --input) input="$2"; shift 2 ;;
    --loop) loops="$2"; shift 2 ;;
    --search) mode="$2"; shift 2 ;;
    *) shift ;;
  esac
done
stem=$(basename "$input" .cor)
if [ "$mode" = peak ]; then
  i=1
  while [ "$i" -le "$loops" ]; do
    delay=+0.000000
    phase=+0.000000
    rate=+0.000000
    case "$stem" in
      *gain_phasecal_merged*)
        delay=+0.250000
        if [ "$i" -eq 7 ]; then delay=+99.000000; fi
        ;;
      *gain_rate_candidate_merged*)
        phase=+0.100000
        rate=+0.000010
        ;;
      *gain_acel_candidate_merged*)
        phase=+0.100000
        rate=+0.000010
        ;;
    esac
    printf "%s\n" "2000/001 00:00:00 L GAIN 1.00 0.100000 20.0 $phase 0.010000 $delay $rate 1 2 3 4 5 6 51544.0 - False False"
    i=$((i + 1))
  done
  exit 0
fi
parent=$(dirname "$input")
if [ "$mode" = rate ]; then
  out="$parent/frinZ/acel_search"
  mkdir -p "$out"
  printf "%s\n" "# Fitted: y = 1.256637061e-3 * x + 2.000000000e-1" > "$out/"$stem"_step1_linear.txt"
  printf "%s\n" "# Corrected Rate: 2.000000000e-4 (from x / (2 * PI))" >> "$out/"$stem"_step1_linear.txt"
  exit 0
fi
if [ "$mode" = acel ]; then
  out="$parent/frinZ/acel_search"
  mkdir -p "$out"
  printf "%s\n" "# Fitted: y = 3.141592654e-6 * x^2 + 1.256637061e-3 * x + 2.000000000e-1" > "$out/"$stem"_step1_quadric.txt"
  printf "%s\n" "# Corrected Acel (Hz/s): 1.000000000e-6 (from x^2 / PI)" >> "$out/"$stem"_step1_quadric.txt"
  printf "%s\n" "# Corrected Rate (Hz): 2.000000000e-4 (from x / (2 * PI))" >> "$out/"$stem"_step1_quadric.txt"
  exit 0
fi
exit 1
"##;
    fs::write(path, script).unwrap();
    let mut permissions = fs::metadata(path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(path, permissions).unwrap();
}

#[test]
fn phased_raw_is_complete_and_can_be_read_by_yi_corr() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let phased_dir = root.join("phased");
    let cor_dir = root.join("cor");
    fs::create_dir_all(&raw_dir).unwrap();
    fs::create_dir_all(&phased_dir).unwrap();
    fs::create_dir_all(&cor_dir).unwrap();

    let schedule = root.join("phased.xml");
    write_schedule(&schedule, "ANT1", "ANT2");
    let tag = "2000001000000";

    // One second at 8192 sample/s and two bits/sample is 2048 bytes.
    let mut state = 0x1234_5678_u32;
    let mut raw = vec![0_u8; 2048];
    for byte in &mut raw {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        *byte = (state >> 24) as u8;
    }
    fs::write(raw_dir.join(format!("ANT1_{tag}.raw")), &raw).unwrap();
    fs::write(raw_dir.join(format!("ANT2_{tag}.raw")), &raw).unwrap();

    let phased = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            phased_dir.to_str().unwrap(),
            "--phased-name",
            "ARRAY",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        phased.status.success(),
        "yi-phasedarray failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&phased.stdout),
        String::from_utf8_lossy(&phased.stderr)
    );

    let phased_raw = phased_dir.join(format!("ARRAY_{tag}.raw"));
    let phased_meta = phased_dir.join(format!("ARRAY_{tag}.raw.meta"));
    assert_eq!(fs::metadata(&phased_raw).unwrap().len(), 2048);
    let phased_bytes = fs::read(&phased_raw).unwrap();
    let mismatches = phased_bytes
        .iter()
        .zip(raw.iter())
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        mismatches <= 32,
        "native VSREC/LSB packed order changed across the phased round trip: {mismatches} mismatched bytes"
    );
    assert!(phased_meta.is_file());
    assert!(!phased_dir.join(format!("ARRAY_{tag}.raw.part")).exists());
    assert!(!phased_dir
        .join(format!("ARRAY_ARRAY_{tag}_phasedarray.cor"))
        .exists());
    let meta = fs::read_to_string(&phased_meta).unwrap();
    assert!(meta.contains("format=yi-phasedarray-raw-v1"));
    assert!(meta.contains(&format!("software_version={}", env!("CARGO_PKG_VERSION"))));
    assert!(meta.contains("virtual_station=ARRAY"));
    assert!(meta.contains("native_format_station=ANT1"));
    assert!(meta.contains("bit=2"));
    assert!(meta.contains("sideband=LSB"));
    assert!(meta.contains("rotation_hz=0.000000000"));
    assert!(meta.contains(
        "shuffle_external=24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7"
    ));

    let corr_schedule = root.join("corr.xml");
    write_schedule(&corr_schedule, "ARRAY", "ANT1");
    let original_raw = raw_dir.join(format!("ANT1_{tag}.raw"));
    let corr = Command::new(env!("CARGO_BIN_EXE_yi-corr"))
        .args([
            "--sc",
            corr_schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--ant1",
            phased_raw.to_str().unwrap(),
            "--ant2",
            original_raw.to_str().unwrap(),
            "--cor",
            cor_dir.to_str().unwrap(),
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        corr.status.success(),
        "yi-corr could not read phased raw:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&corr.stdout),
        String::from_utf8_lossy(&corr.stderr)
    );
    let xcf = cor_dir.join(format!("ARRAY_ANT1_{tag}_phasedarray.cor"));
    assert!(fs::metadata(xcf).unwrap().len() > 0);

    let _ = fs::remove_dir_all(root);
}

#[test]
fn phased_array_writes_gain_phasecal_off_and_on_xcf() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let output_dir = root.join("phased");
    fs::create_dir_all(&raw_dir).unwrap();
    fs::create_dir_all(&output_dir).unwrap();

    let uncalibrated = root.join("gain-uncalibrated.xml");
    let calibrated = root.join("gain-calibrated.xml");
    write_schedule_with_clock_delay(&uncalibrated, "YAMAGU32", "YAMAGU34", 0.0);
    write_schedule_with_clock_delay(&calibrated, "YAMAGU32", "YAMAGU34", 8.0 / 8192.0);
    let tag = "2000001000000";
    let mut state = 0x2468_ace1_u32;
    let mut raw = vec![0_u8; 2048];
    for byte in &mut raw {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        *byte = (state >> 24) as u8;
    }
    fs::write(raw_dir.join(format!("YAMAGU32_{tag}.raw")), &raw).unwrap();
    fs::write(raw_dir.join(format!("YAMAGU34_{tag}.raw")), &raw).unwrap();

    let phased = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .args([
            "--sc",
            calibrated.to_str().unwrap(),
            "--gain-uncalibrated-sc",
            uncalibrated.to_str().unwrap(),
            "--gain-reference-key",
            "A",
            "--gain-corrected-key",
            "B",
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--phased-name",
            "YAMAGU66",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        phased.status.success(),
        "gain comparison failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&phased.stdout),
        String::from_utf8_lossy(&phased.stderr)
    );
    let stdout = String::from_utf8_lossy(&phased.stdout);
    assert!(stdout.contains("[gain-corr] phase-reference gain comparison enabled"));
    assert!(stdout.contains("reference=A fixed; corrected=B"));
    assert!(stdout.contains("[gain-corr] complete:"));

    let gain_dir = output_dir.join("gain_correlation");
    let off = gain_dir.join(format!("YAMAGU32_YAMAGU34_{tag}_gain_phasecal_off.cor"));
    let on = gain_dir.join(format!("YAMAGU32_YAMAGU34_{tag}_gain_phasecal_on.cor"));
    assert!(fs::metadata(&off).unwrap().len() > 0);
    assert!(fs::metadata(&on).unwrap().len() > 0);
    assert_ne!(fs::read(off).unwrap(), fs::read(on).unwrap());
    for unwanted in [
        format!("YAMAGU32_YAMAGU32_{tag}_gain_phasecal_off.cor"),
        format!("YAMAGU34_YAMAGU34_{tag}_gain_phasecal_off.cor"),
        format!("YAMAGU32_YAMAGU32_{tag}_gain_phasecal_on.cor"),
        format!("YAMAGU34_YAMAGU34_{tag}_gain_phasecal_on.cor"),
    ] {
        assert!(!gain_dir.join(unwanted).exists());
    }

    let phased_raw = output_dir.join(format!("YAMAGU66_{tag}.raw"));
    assert_eq!(fs::metadata(phased_raw).unwrap().len(), raw.len() as u64);
    let metadata = fs::read_to_string(output_dir.join(format!("YAMAGU66_{tag}.raw.meta"))).unwrap();
    assert!(metadata.contains("gain_uncalibrated_schedule="));
    assert!(metadata.contains("gain_correlation_directory="));

    let _ = fs::remove_dir_all(root);
}

#[test]
fn phased_raw_preserves_scan_length_across_nonzero_delay_seek() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let phased_dir = root.join("phased");
    fs::create_dir_all(&raw_dir).unwrap();
    fs::create_dir_all(&phased_dir).unwrap();

    let schedule = root.join("phased-delay.xml");
    // Eight samples of clock delay. Before scan-length preservation this
    // reduced the common input overlap below 64 complete 128-sample frames and
    // produced a raw file one FFT frame shorter than the one-second inputs.
    write_schedule_with_clock_delay(&schedule, "ANT1", "ANT2", 8.0 / 8192.0);
    let tag = "2000001000000";
    let raw = vec![0x1b_u8; 2048];
    fs::write(raw_dir.join(format!("ANT1_{tag}.raw")), &raw).unwrap();
    fs::write(raw_dir.join(format!("ANT2_{tag}.raw")), &raw).unwrap();

    let phased = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            phased_dir.to_str().unwrap(),
            "--phased-name",
            "ARRAY",
            "--phased-diagnostics",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        phased.status.success(),
        "yi-phasedarray failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&phased.stdout),
        String::from_utf8_lossy(&phased.stderr)
    );

    let stdout = String::from_utf8_lossy(&phased.stdout);
    assert!(
        stdout.contains("Requested scan-length preservation"),
        "missing boundary-padding diagnostic:\n{stdout}"
    );
    let phased_raw = phased_dir.join(format!("ARRAY_{tag}.raw"));
    assert_eq!(
        fs::metadata(phased_raw).unwrap().len(),
        raw.len() as u64,
        "virtual-station raw must retain the nominal one-second scan length"
    );
    for product in [
        format!("ANT1_ANT1_{tag}_phasedarray.cor"),
        format!("ANT1_ANT2_{tag}_phasedarray.cor"),
        format!("ANT2_ANT2_{tag}_phasedarray.cor"),
        format!("ARRAY_ARRAY_{tag}_phasedarray.cor"),
    ] {
        assert!(
            fs::metadata(phased_dir.join(&product)).unwrap().len() > 0,
            "missing phased diagnostic product {product}"
        );
    }

    let _ = fs::remove_dir_all(root);
}

#[test]
fn phased_validation_runs_complete_three_station_workflow() {
    if !Command::new("python3")
        .args(["-c", "import numpy"])
        .status()
        .is_ok_and(|status| status.success())
    {
        eprintln!("skipping phased-validation E2E because python3/numpy is unavailable");
        return;
    }

    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let output_dir = root.join("output");
    fs::create_dir_all(&raw_dir).unwrap();
    fs::create_dir_all(&output_dir).unwrap();
    let schedule = root.join("three-station.xml");
    write_three_station_schedule(&schedule);
    let tag = "2000001000000";
    let mut state = 0x9e37_79b9_u32;
    let mut raw = vec![0_u8; 2048];
    for byte in &mut raw {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        *byte = (state >> 24) as u8;
    }
    for station in ["YAMAGU32", "YAMAGU34", "HITACH32"] {
        fs::write(raw_dir.join(format!("{station}_{tag}.raw")), &raw).unwrap();
    }

    let phased = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--phased-name",
            "YAMAGU66",
            "--phased-validation",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        phased.status.success(),
        "three-station phased validation failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&phased.stdout),
        String::from_utf8_lossy(&phased.stderr)
    );
    let stdout = String::from_utf8_lossy(&phased.stdout);
    assert!(stdout.contains("phased pair: YAMAGU32-YAMAGU34"));
    assert!(stdout.contains("stage 4/4: build complete visibility NPZ"));
    assert!(stdout.contains("[validation] complete:"));

    assert_eq!(
        fs::metadata(output_dir.join(format!("YAMAGU66_{tag}.raw")))
            .unwrap()
            .len(),
        raw.len() as u64
    );
    assert!(
        fs::metadata(output_dir.join(format!("YAMAGU66_{tag}_visibility_validation.npz")))
            .unwrap()
            .len()
            > 0
    );
    let validation_dir = output_dir.join("phased_validation");
    for product in [
        format!("YAMAGU32_HITACH32_{tag}_validation.cor"),
        format!("YAMAGU34_HITACH32_{tag}_validation.cor"),
        format!("YAMAGU66_HITACH32_{tag}_validation.cor"),
        format!("YAMAGU66_YAMAGU66_{tag}_validation.cor"),
    ] {
        assert!(
            fs::metadata(validation_dir.join(&product)).unwrap().len() > 0,
            "missing validation product {product}"
        );
    }

    let _ = fs::remove_dir_all(root);
}

#[cfg(unix)]
#[test]
fn automatic_gain_phasecal_runs_frinz_updates_l_and_synthesizes_every_scan() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let output_dir = root.join("output");
    fs::create_dir_all(&raw_dir).unwrap();
    fs::create_dir_all(&output_dir).unwrap();
    let schedule = root.join("auto-phasecal.xml");
    write_phasecal_schedule(&schedule);
    let fake_frinz = root.join("fake-frinz");
    write_fake_frinz(&fake_frinz);

    let mut state = 0x1357_9bdf_u32;
    let mut raw = vec![0_u8; 2048];
    for byte in &mut raw {
        state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        *byte = (state >> 24) as u8;
    }
    let tags = [
        "2000001000000",
        "2000001000010",
        "2000001000020",
        "2000001000030",
        "2000001000040",
    ];
    for tag in tags {
        fs::write(raw_dir.join(format!("YAMAGU32_{tag}.raw")), &raw).unwrap();
        fs::write(raw_dir.join(format!("YAMAGU34_{tag}.raw")), &raw).unwrap();
    }

    let phased = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .env("YI_FRINZ", &fake_frinz)
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--phased-name",
            "YAMAGU66",
            "--gain-phasecal",
            "--gain-source",
            "GAIN",
            "--gain-reference-key",
            "A",
            "--gain-corrected-key",
            "B",
            "--gain-fringe-length",
            "1",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        phased.status.success(),
        "automatic gain phasecal failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&phased.stdout),
        String::from_utf8_lossy(&phased.stderr)
    );
    let stdout = String::from_utf8_lossy(&phased.stdout);
    assert!(stdout.contains("[phasecal] stage 1/7: phasecal-off gain correlations"));
    assert!(stdout.contains("[phasecal] stage 2/7: merge gain scans and run frinZ --search peak"));
    assert!(stdout.contains("[phasecal] stage 3/7: correlate gain scans with median group delay"));
    assert!(stdout.contains("[phasecal] stage 4/7: solve rate and acceleration candidates"));
    assert!(stdout.contains("[phasecal] frinZ candidate 1/2: --search rate"));
    assert!(stdout.contains("[phasecal] frinZ candidate 2/2: --search acel"));
    assert!(stdout.contains("[phasecal] selected rate model by minimum BIC"));
    assert!(stdout.contains("[phasecal] stage 5/7: final phasecal-on gain correlations"));
    assert!(stdout.contains("[phasecal] stage 6/7: synthesize corrected target and gain scans"));
    assert!(stdout.contains("[phasecal] stage 7/7: complete"));

    let corrected = output_dir.join("auto-phasecal_B_phasecal.xml");
    let corrected_xml = fs::read_to_string(&corrected).unwrap();
    assert!(corrected_xml.contains("<clock key=\"A\"><epoch>2000/001 00:00:00</epoch><delay>0</delay><rate>0</rate><accel>0</accel></clock>"));
    assert!(corrected_xml.contains("<clock key=\"B\">"));
    assert!(corrected_xml.contains("<epoch>2000/001 00:00:00</epoch>"));
    assert!(!corrected_xml
        .contains("<clock key=\"B\"><epoch>2000/001 00:00:00</epoch><delay>1e-6</delay>"));

    let gain_dir = output_dir.join("gain_correlation");
    for tag in [tags[0], tags[2], tags[4]] {
        for label in ["gain_phasecal_off", "gain_delay_on", "gain_phasecal_on"] {
            assert!(
                fs::metadata(gain_dir.join(format!("YAMAGU32_YAMAGU34_{tag}_{label}.cor")))
                    .unwrap()
                    .len()
                    > 256
            );
        }
    }
    let solution = fs::read_to_string(gain_dir.join("gain_phasecal_solution.txt")).unwrap();
    assert!(solution.contains("gain_scans=3"));
    assert!(solution.contains("fringe_windows=24"));
    assert!(solution.contains("format=yi-phasedarray-gain-phasecal-v4"));
    assert!(solution.contains("peak_delay_count=24"));
    assert!(solution.contains("peak_delay_samples_median=+2.50000000000000000e-1"));
    assert!(solution.contains("peak_delay_samples_max=+9.90000000000000000e1"));
    assert!(solution.contains("group_delay_s=+3.05175781250000000e-5"));
    assert!(solution.contains("phase_unit=radian"));
    assert!(solution.contains("selected_model=rate"));
    assert!(solution.contains("rate_phase_std_rad="));
    assert!(solution.contains("acel_phase_std_rad="));
    assert!(solution.contains("rate_bic="));
    assert!(solution.contains("acel_bic="));
    assert!(solution.contains("corrected_schedule="));
    assert!(gain_dir
        .join("YAMAGU32_YAMAGU34_gain_phasecal_merged.cor")
        .is_file());
    assert!(gain_dir
        .join("YAMAGU32_YAMAGU34_gain_delay_on_merged.cor")
        .is_file());
    assert!(gain_dir.join("auto-phasecal_B_group_delay.xml").is_file());

    for tag in tags {
        assert_eq!(
            fs::metadata(output_dir.join(format!("YAMAGU66_{tag}.raw")))
                .unwrap()
                .len(),
            raw.len() as u64
        );
    }

    let resume_path = gain_dir.join("gain_phasecal.resume");
    let resume = fs::read_to_string(&resume_path).unwrap();
    assert!(resume.contains("format=yi-phasedarray-gain-phasecal-resume-v1"));
    assert!(resume.contains("completed=off:2000001000000"));
    assert!(resume.contains("completed=workflow:complete"));

    // Simulate upgrading from a pre-resume release: complete phasecal-off .cor
    // files exist but no state file does. They must be adopted, not recomputed.
    fs::remove_file(&resume_path).unwrap();
    let legacy_resume = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .env("YI_FRINZ", &fake_frinz)
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--phased-name",
            "YAMAGU66",
            "--gain-phasecal",
            "--gain-source",
            "GAIN",
            "--gain-reference-key",
            "A",
            "--gain-corrected-key",
            "B",
            "--gain-fringe-length",
            "1",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(legacy_resume.status.success());
    let legacy_stdout = String::from_utf8_lossy(&legacy_resume.stdout);
    assert!(legacy_stdout.contains("[resume] adopted complete legacy phasecal-off scan"));

    // A subsequent identical run must reuse every expensive scan stage.
    let resumed = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .env("YI_FRINZ", &fake_frinz)
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--phased-name",
            "YAMAGU66",
            "--gain-phasecal",
            "--gain-source",
            "GAIN",
            "--gain-reference-key",
            "A",
            "--gain-corrected-key",
            "B",
            "--gain-fringe-length",
            "1",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(
        resumed.status.success(),
        "resume run failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&resumed.stdout),
        String::from_utf8_lossy(&resumed.stderr)
    );
    let resumed_stdout = String::from_utf8_lossy(&resumed.stdout);
    assert!(resumed_stdout.contains("[resume] reuse phasecal-off scan"));
    assert!(resumed_stdout.contains("[resume] reuse group-delay scan"));
    assert!(resumed_stdout.contains("[resume] reuse rate candidate scan"));
    assert!(resumed_stdout.contains("[resume] reuse acel candidate scan"));
    assert!(resumed_stdout.contains("[resume] reuse phasecal-on scan"));
    assert!(resumed_stdout.contains("[resume] reuse phased scan"));

    let changed_conditions = Command::new(env!("CARGO_BIN_EXE_yi-phasedarray"))
        .env("YI_FRINZ", &fake_frinz)
        .args([
            "--sc",
            schedule.to_str().unwrap(),
            "--raw",
            raw_dir.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--phased-name",
            "YAMAGU66",
            "--gain-phasecal",
            "--gain-source",
            "GAIN",
            "--gain-reference-key",
            "A",
            "--gain-corrected-key",
            "B",
            "--gain-fringe-length",
            "2",
            "--cpu",
            "1",
        ])
        .output()
        .unwrap();
    assert!(!changed_conditions.status.success());
    assert!(String::from_utf8_lossy(&changed_conditions.stderr)
        .contains("gain resume conditions changed"));

    let _ = fs::remove_dir_all(root);
}

#[test]
fn no_affinity_preserves_products_and_ignores_pinning_configuration() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    let home = root.join("home");
    let cargo_home = root.join("cargo");
    fs::create_dir_all(&raw_dir).unwrap();
    fs::create_dir_all(&home).unwrap();
    fs::create_dir_all(cargo_home.join("tmp")).unwrap();
    let schedule = root.join("test.xml");
    write_schedule(&schedule, "ANT1", "ANT2");
    let raw: Vec<u8> = (0..2048).map(|i| ((i * 53 + i / 7) % 256) as u8).collect();
    for station in ["ANT1", "ANT2"] {
        fs::write(raw_dir.join(format!("{station}_2000001000000.raw")), &raw).unwrap();
    }
    let config = cargo_home.join("tmp/yi-corr-affinity.txt");
    let registry = home.join(".yi-corr");
    for phased in [false, true] {
        for cpu in ["1", "3"] {
            let mut reference = std::collections::BTreeMap::new();
            for unpinned in [false, true] {
                let out = root.join(format!("out-{phased}-{cpu}-{unpinned}"));
                let mut cmd = Command::new(if phased {
                    env!("CARGO_BIN_EXE_yi-phasedarray")
                } else {
                    env!("CARGO_BIN_EXE_yi-corr")
                });
                cmd.args([
                    "--sc",
                    schedule.to_str().unwrap(),
                    "--raw",
                    raw_dir.to_str().unwrap(),
                    "--cor",
                    out.to_str().unwrap(),
                    "--cpu",
                    cpu,
                    "--chunk-frames",
                    "3",
                    "--pipeline-depth",
                    "1",
                    "--usb",
                ])
                .env("HOME", &home)
                .env("CARGO_HOME", &cargo_home)
                .env_remove("YI_READER_CORE");
                if unpinned {
                    // Invalid inputs must not even be read in this mode.
                    fs::write(&config, "invalid affinity specification").unwrap();
                    fs::write(&registry, "registry must remain untouched").unwrap();
                    cmd.arg("--no-affinity").env("YI_READER_CORE", "invalid");
                }
                let result = cmd.output().unwrap();
                let stdout = String::from_utf8_lossy(&result.stdout);
                assert!(
                    result.status.success(),
                    "stdout:\n{stdout}\nstderr:\n{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                assert!(stdout.contains(if cpu == "1" {
                    "compute-threads=1"
                } else {
                    "compute-threads=2"
                }));
                if unpinned {
                    assert!(stdout.contains("affinity=disabled"));
                    assert!(!stdout
                        .lines()
                        .any(|line| line.starts_with("[info] CPU allocation:")
                            && line.contains("I/O CPU=")));
                    assert!(!stdout.contains("Cross-process CPU reservation:"));
                    assert_eq!(
                        fs::read_to_string(&registry).unwrap(),
                        "registry must remain untouched"
                    );
                    fs::remove_file(&config).unwrap();
                    fs::remove_file(&registry).unwrap();
                }
                let products: std::collections::BTreeMap<_, _> = fs::read_dir(&out)
                    .unwrap()
                    .map(|e| e.unwrap().path())
                    .filter(|p| {
                        matches!(p.extension().and_then(|s| s.to_str()), Some("cor" | "raw"))
                    })
                    .map(|p| (p.file_name().unwrap().to_owned(), fs::read(&p).unwrap()))
                    .collect();
                assert!(!products.is_empty());
                if unpinned {
                    assert_eq!(products, reference, "phased={phased} cpu={cpu}");
                } else {
                    reference = products;
                }
            }
        }
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn small_io_chunks_preserve_integration_and_phased_raw_with_delay_drift() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    fs::create_dir_all(&raw_dir).unwrap();
    let raw: Vec<u8> = (0..4096)
        .map(|i| ((i * 53 + i / 7 + 79) % 256) as u8)
        .collect();
    for station in ["ANT1", "ANT2"] {
        fs::write(raw_dir.join(format!("{station}_2000001000000.raw")), &raw).unwrap();
    }
    for rate in [-0.04, 0.0, 0.04] {
        let schedule = root.join("test.xml");
        write_schedule_with_clock_delay(&schedule, "ANT1", "ANT2", 13.25 / 8192.);
        let xml = fs::read_to_string(&schedule).unwrap();
        fs::write(
            &schedule,
            xml.replace(
                "</delay><rate>0</rate></clock>\n  <terminal",
                &format!("</delay><rate>{rate}</rate></clock>\n  <terminal"),
            ),
        )
        .unwrap();
        for phased in [false, true] {
            let mut reference = std::collections::BTreeMap::new();
            for (run, chunk, usb) in [(0, "1000", false), (1, "3", false), (2, "3", true)] {
                let out = root.join(format!("out-{rate}-{phased}-{run}"));
                let mut cmd = Command::new(if phased {
                    env!("CARGO_BIN_EXE_yi-phasedarray")
                } else {
                    env!("CARGO_BIN_EXE_yi-corr")
                });
                cmd.args([
                    "--sc",
                    schedule.to_str().unwrap(),
                    "--raw",
                    raw_dir.to_str().unwrap(),
                    if phased { "--output" } else { "--cor" },
                    out.to_str().unwrap(),
                    "--cpu",
                    "2",
                    "--chunk-frames",
                    chunk,
                    "--pipeline-depth",
                    "1",
                ]);
                if phased {
                    cmd.args(["--phased-name", "ARRAY"]);
                }
                if usb {
                    cmd.arg("--usb");
                }
                let result = cmd.output().unwrap();
                assert!(
                    result.status.success(),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let products: std::collections::BTreeMap<_, _> = fs::read_dir(&out)
                    .unwrap()
                    .map(|e| e.unwrap().path())
                    .filter(|p| {
                        matches!(p.extension().and_then(|s| s.to_str()), Some("cor" | "raw"))
                    })
                    .map(|p| (p.file_name().unwrap().to_owned(), fs::read(&p).unwrap()))
                    .collect();
                assert!(!products.is_empty());
                if run == 0 {
                    reference = products;
                } else {
                    assert_eq!(products, reference, "rate={rate} phased={phased} usb={usb}");
                }
            }
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn reused_workers_preserve_rotated_mixed_bit_products() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    fs::create_dir_all(&raw_dir).unwrap();
    for (station, multiplier) in [("ANT1", 53), ("ANT2", 71)] {
        let raw: Vec<u8> = (0..8192)
            .map(|i| ((i * multiplier + i / 7 + 79) % 256) as u8)
            .collect();
        fs::write(raw_dir.join(format!("{station}_2000001000000.raw")), raw).unwrap();
    }
    for (bits1, bits2, rot1, rot2, sideband) in [
        (2, 2, "512", "256", "LSB"),
        (2, 4, "0", "-512", "USB"),
        (4, 2, "-256", "0", "LSB"),
    ] {
        let schedule = root.join("test.xml");
        write_schedule_with_clock_delay(&schedule, "ANT1", "ANT2", 13.25 / 8192.0);
        let xml = fs::read_to_string(&schedule)
            .unwrap()
            .replacen("<terminal>term</terminal>", "<terminal>term1</terminal>", 1)
            .replacen("<terminal>term</terminal>", "<terminal>term2</terminal>", 1);
        let terminal = |name: &str, bits: usize| {
            let count = 1usize << bits;
            let levels = (0..count)
                .map(|i| (i as f64 - (count - 1) as f64 / 2.0).to_string())
                .collect::<Vec<_>>()
                .join(",");
            format!("<terminal name=\"{name}\"><speed>8192</speed><channel>1</channel><bit>{bits}</bit><level>{levels}</level></terminal>")
        };
        let xml = xml.replace(
            "<terminal name=\"term\"><speed>8192</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,0.5,1.5</level></terminal>",
            &(terminal("term1", bits1) + &terminal("term2", bits2)),
        ).replace(
            "<special key=\"A\"><rotation>0</rotation><sideband>LSB</sideband></special>",
            &format!("<special key=\"A\"><rotation>{rot1}</rotation><sideband>{sideband}</sideband></special>"),
        ).replace(
            "<special key=\"B\"><rotation>0</rotation><sideband>LSB</sideband></special>",
            &format!("<special key=\"B\"><rotation>{rot2}</rotation><sideband>{sideband}</sideband></special>"),
        );
        fs::write(&schedule, xml).unwrap();
        for phased in [false, true] {
            let binary = if phased { "yi-phasedarray" } else { "yi-corr" };
            let current = if phased {
                env!("CARGO_BIN_EXE_yi-phasedarray")
            } else {
                env!("CARGO_BIN_EXE_yi-corr")
            };
            let mut runs = vec![
                (PathBuf::from(current), "1", "1000"),
                (PathBuf::from(current), "4", "5"),
            ];
            // Optional comparison with a saved pre-change release build.
            if let Some(directory) = std::env::var_os("FX_CORR_BASELINE_DIR") {
                runs.push((PathBuf::from(directory).join(binary), "1", "1000"));
            }
            let mut reference = std::collections::BTreeMap::new();
            for (run, (executable, cpu, chunk)) in runs.into_iter().enumerate() {
                let out = root.join(format!("out-{bits1}-{bits2}-{phased}-{run}"));
                let result = Command::new(executable)
                    .args([
                        "--sc",
                        schedule.to_str().unwrap(),
                        "--raw",
                        raw_dir.to_str().unwrap(),
                        "--output",
                        out.to_str().unwrap(),
                        "--cpu",
                        cpu,
                        "--no-affinity",
                        "--chunk-frames",
                        chunk,
                        "--pipeline-depth",
                        "1",
                        "--phased-name",
                        "ARRAY",
                    ])
                    .args(phased.then_some("--phased-diagnostics"))
                    .output()
                    .unwrap();
                assert!(
                    result.status.success(),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let products: std::collections::BTreeMap<_, _> = fs::read_dir(&out)
                    .unwrap()
                    .map(|e| e.unwrap().path())
                    .filter(|p| {
                        matches!(p.extension().and_then(|s| s.to_str()), Some("cor" | "raw"))
                    })
                    .map(|p| (p.file_name().unwrap().to_owned(), fs::read(p).unwrap()))
                    .collect();
                assert!(!products.is_empty());
                if run == 0 {
                    reference = products;
                } else {
                    assert_eq!(
                        products, reference,
                        "bits={bits1}/{bits2} phased={phased} run={run}"
                    );
                }
            }
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn vlbi_flag_overrides_model_env_and_applies_troposphere() {
    let root = unique_temp_dir();
    let raw_dir = root.join("raw");
    fs::create_dir_all(&raw_dir).unwrap();
    let raw: Vec<u8> = (0..4096)
        .map(|i| ((i * 53 + i / 7 + 79) % 256) as u8)
        .collect();
    for station in ["ANT1", "ANT2"] {
        fs::write(raw_dir.join(format!("{station}_2000001000000.raw")), &raw).unwrap();
    }
    for (baseline, xyz) in [
        (
            "short",
            "<pos-x>-3502567.576</pos-x><pos-y>3950885.734</pos-y><pos-z>3566449.115</pos-z>",
        ),
        (
            "long",
            "<pos-x>-3961788.974</pos-x><pos-y>3243597.492</pos-y><pos-z>3790597.692</pos-z>",
        ),
    ] {
        let schedule = root.join("test.xml");
        write_schedule(&schedule, "ANT1", "ANT2");
        let xml = fs::read_to_string(&schedule).unwrap();
        let original = "<station key=\"B\"><name>ANT2</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z>";
        assert!(xml.contains(original));
        fs::write(
            &schedule,
            xml.replace(
                original,
                &format!("<station key=\"B\"><name>ANT2</name>{xyz}"),
            ),
        )
        .unwrap();
        let mut products = Vec::new();
        for run in 0..5 {
            let out = root.join(format!("{baseline}-{run}"));
            let mut cmd = Command::new(env!("CARGO_BIN_EXE_yi-corr"));
            cmd.args([
                "--sc",
                schedule.to_str().unwrap(),
                "--raw",
                raw_dir.to_str().unwrap(),
                "--cor",
                out.to_str().unwrap(),
                "--cpu",
                "2",
            ])
            .env_remove("YI_SOURCE_VECTOR_MODE")
            .env_remove("YI_GEOM_DELAY_MODE")
            .env("YI_EOP_MODE", "none");
            if run == 1 {
                cmd.env("YI_SOURCE_VECTOR_MODE", "mean-gast")
                    .env("YI_GEOM_DELAY_MODE", "anchored");
            }
            if run == 2 {
                cmd.env("YI_SOURCE_VECTOR_MODE", "pnm-gast")
                    .env("YI_GEOM_DELAY_MODE", "vlbi-minus");
            }
            if run == 3 {
                cmd.arg("--vlbi")
                    .env("YI_SOURCE_VECTOR_MODE", "mean-gast")
                    .env("YI_GEOM_DELAY_MODE", "vlbi-plus");
            }
            if run == 4 {
                cmd.arg("--vlbi");
            }
            let result = cmd.output().unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            if run == 3 {
                let paths = corr_info_paths(&out);
                assert_eq!(paths.len(), 1);
                let info = read_corr_info(&paths[0]);
                assert!(info.contains("VLBI preset:"));
                assert!(info.contains("VLBI atmosphere: nominal troposphere enabled"));
                assert!(String::from_utf8_lossy(&result.stdout)
                    .contains("atmosphere=nominal-troposphere"));
            }
            let files: std::collections::BTreeMap<_, _> = fs::read_dir(&out)
                .unwrap()
                .map(|e| e.unwrap().path())
                .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("cor"))
                .map(|p| (p.file_name().unwrap().to_owned(), fs::read(&p).unwrap()))
                .collect();
            assert!(!files.is_empty());
            products.push(files);
        }
        assert_eq!(
            products[0], products[1],
            "legacy defaults changed: {baseline}"
        );
        assert_ne!(
            products[2], products[3],
            "VLBI troposphere did not affect the model: {baseline}"
        );
        assert_eq!(
            products[3], products[4],
            "VLBI preset did not override conflicting environment models: {baseline}"
        );
        assert_ne!(
            products[0], products[3],
            "model switch had no effect: {baseline}"
        );
    }
    fs::remove_dir_all(root).unwrap();
}
