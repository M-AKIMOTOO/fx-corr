use std::{
    f64::consts::TAU,
    fs,
    io::Read,
    path::{Path, PathBuf},
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

fn temp_root() -> PathBuf {
    std::env::temp_dir().join(format!(
        "fx-corr-unequal-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ))
}

fn info(dir: &Path) -> String {
    let path = fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .find(|p| p.extension().and_then(|s| s.to_str()) == Some("info"))
        .unwrap();
    let mut zip = zip::ZipArchive::new(fs::File::open(path).unwrap()).unwrap();
    let mut text = String::new();
    zip.by_name("info.txt")
        .unwrap()
        .read_to_string(&mut text)
        .unwrap();
    text
}

fn spectrum(path: &Path, nfft: usize) -> Vec<num_complex::Complex<f64>> {
    let bytes = fs::read(path).unwrap();
    assert_eq!(bytes.len(), 256 + 128 + nfft / 2 * 8);
    bytes[384..]
        .chunks_exact(8)
        .map(|b| {
            num_complex::Complex::new(
                f32::from_le_bytes(b[..4].try_into().unwrap()) as f64,
                f32::from_le_bytes(b[4..].try_into().unwrap()) as f64,
            )
        })
        .collect()
}

#[test]
fn native_1024_and_32_msps_correlate_only_the_contained_16_mhz() {
    let root = temp_root();
    let raw = root.join("raw");
    fs::create_dir_all(&raw).unwrap();
    let nfft = 256usize;
    let narrow_fs = 32_000_000.0;
    let wide_fs = 1_024_000_000.0;
    let df = narrow_fs / nfft as f64;
    let epoch = "2000/001 00:00:00";
    // Lower edge, middle and upper edge containment; reverse station order;
    // mixed sidebands and both ordinary 8-bit and packed 2-bit data.
    for (case, wide_first, offset_mhz, sidebands, bits) in [
        (0, true, 0.0, ["USB", "USB"], 8),
        (1, true, 8.0, ["LSB", "USB"], 8),
        (2, false, 496.0, ["LSB", "LSB"], 8),
        (3, false, 8.0, ["USB", "LSB"], 2),
        (4, true, 8.0, ["USB", "LSB"], 8),
    ] {
        let rates = if wide_first {
            [wide_fs, narrow_fs]
        } else {
            [narrow_fs, wide_fs]
        };
        let rotations = if wide_first {
            [0.0, offset_mhz * 1e6]
        } else {
            [offset_mhz * 1e6, 0.0]
        };
        let delay = if wide_first { 3.25e-9 } else { -2.75e-9 };
        let tones: &[usize] = if bits == 2 { &[23] } else { &[5, 23, 71, 117] };
        let output_low = 6600e6 + offset_mhz * 1e6;
        let levels = if bits == 8 {
            (0..256)
                .map(|i| (i as f64 - 127.5).to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            "-1.5,-0.5,0.5,1.5".into()
        };
        let shuffle = (0..32)
            .rev()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let station = |key: &str, name: &str, terminal: &str| {
            format!("<station key=\"{key}\"><name>{name}</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z><terminal>{terminal}</terminal></station>")
        };
        let xml = format!(
            r#"<schedule>
{}
{}
<clock key="A"><delay>0</delay><rate>0</rate></clock>
<clock key="B"><delay>{delay:.17e}</delay><rate>0</rate></clock>
<terminal name="t1"><speed>{}</speed><bit>{bits}</bit><level>{levels}</level></terminal>
<terminal name="t2"><speed>{}</speed><bit>{bits}</bit><level>{levels}</level></terminal>
<shuffle key="A">{shuffle}</shuffle><shuffle key="B">{shuffle}</shuffle>
<source name="SYNTH"><ra>00h00m00</ra><dec>+00d00m00</dec></source>
<stream><frequency>6600000000</frequency><fft>{nfft}</fft><output>7812.5</output>
<special key="A"><rotation>{}</rotation><sideband>{}</sideband></special>
<special key="B"><rotation>{}</rotation><sideband>{}</sideband></special>
</stream>
<process><epoch>{epoch}</epoch><length>{}</length><object>SYNTH</object><stations>AB</stations></process>
</schedule>"#,
            station("A", "ANT1", "t1"),
            station("B", "ANT2", "t2"),
            rates[0],
            rates[1],
            rotations[0],
            sidebands[0],
            rotations[1],
            sidebands[1],
            if case == 0 { 0.000128 } else { 0.000112 }
        );
        let xml = if case == 4 {
            xml.replace(
                "<name>ANT1</name>",
                "<name>ANT1</name><atmosphere stec-tecu=\"0\"/>",
            )
            .replace(
                "<name>ANT2</name>",
                "<name>ANT2</name><atmosphere stec-tecu=\"1\"/>",
            )
        } else {
            xml
        };
        for ant in 0..2 {
            let native_n = nfft * (rates[ant] / narrow_fs) as usize;
            let mut bytes = vec![0u8; native_n * 16 * bits / 8];
            for sample in 0..native_n * 16 {
                let t = sample as f64 / rates[ant];
                let mut v = 0.0;
                for (j, &bin) in tones.iter().enumerate() {
                    let rf = output_low + bin as f64 * df;
                    let baseband = rf - (6600e6 + rotations[ant]);
                    let phase = TAU * baseband * t + j as f64 * 0.27
                        - if ant == 1 { TAU * rf * delay } else { 0.0 }
                        + if ant == 1 && case == 4 {
                            TAU * (40.3e16 / 299_792_458.0) / rf
                        } else {
                            0.0
                        };
                    v += if bits == 8 { 12.0 } else { 1.3 } * phase.cos();
                }
                // A strong wide-band line outside the common RF must never
                // fold into the selected 16 MHz or change its amplitudes.
                if rates[ant] == wide_fs && bits == 8 {
                    v += 25.0 * (TAU * 250e6 * t + 0.83).cos();
                }
                if sidebands[ant] == "LSB" && sample % 2 != 0 {
                    v = -v;
                }
                if bits == 8 {
                    bytes[sample] = (v + 127.5).round().clamp(0.0, 255.0) as u8;
                } else {
                    let code = if v < -1.0 {
                        0
                    } else if v < 0.0 {
                        1
                    } else if v < 1.0 {
                        2
                    } else {
                        3
                    };
                    bytes[sample / 4] |= code << (2 * (sample % 4));
                }
            }
            fs::write(raw.join(format!("ANT{}_2000001000000.raw", ant + 1)), bytes).unwrap();
        }
        let mut reference = Vec::new();
        for mode in ["direct", "debug", "fold", "short-chunks"] {
            let schedule = root.join(format!("{case}-{mode}.xml"));
            fs::write(
                &schedule,
                if mode == "fold" {
                    xml.replace(
                        "<stream>",
                        "<stream><pulsar><period>1</period><bins>1</bins><dm>0</dm></pulsar>",
                    )
                } else {
                    xml.clone()
                },
            )
            .unwrap();
            let out = root.join(format!("{case}-{mode}"));
            let mut cmd = Command::new(env!("CARGO_BIN_EXE_yi-corr"));
            for (key, _) in
                std::env::vars_os().filter(|(k, _)| k.to_string_lossy().starts_with("YI_"))
            {
                cmd.env_remove(key);
            }
            cmd.args([
                "--sc",
                schedule.to_str().unwrap(),
                "--raw",
                raw.to_str().unwrap(),
                "--cor",
                out.to_str().unwrap(),
                "--cpu",
                "1",
                "--no-affinity",
            ])
            .env("YI_EOP_MODE", "none");
            if case != 0 {
                cmd.args(["--skip", "0.00001778125"]);
            }
            if mode == "debug" {
                cmd.arg("--debug");
            }
            if mode == "short-chunks" {
                cmd.args(["--chunk-frames", "3", "--read-batch-mib", "1", "--usb"]);
            }
            let result = cmd.output().unwrap();
            let stdout = String::from_utf8_lossy(&result.stdout);
            assert!(
                result.status.success(),
                "case={case} {mode}: {stdout}\n{}",
                String::from_utf8_lossy(&result.stderr)
            );
            assert!(stdout.contains("BW=16.000000MHz"), "{stdout}");
            assert!(
                stdout.contains("sampling=1024.000000MS/s FFT=8192"),
                "{stdout}"
            );
            assert!(
                stdout.contains("sampling=32.000000MS/s FFT=256"),
                "{stdout}"
            );
            assert!(stdout.contains("FFT-window=0.000008000s"), "{stdout}");
            let details = info(&out);
            assert!(!details.contains("skipped "), "{details}");
            if mode == "fold" {
                assert!(details.contains("Pulsar folding:"));
            }
            let path = fs::read_dir(&out)
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
            let bytes = fs::read(&path).unwrap();
            // Header RF and sampling describe the actual narrow output grid.
            assert_eq!(
                f64::from_le_bytes(bytes[16..24].try_into().unwrap()),
                output_low
            );
            assert_eq!(
                i32::from_le_bytes(bytes[12..16].try_into().unwrap()),
                narrow_fs as i32
            );
            assert_eq!(
                i32::from_le_bytes(bytes[24..28].try_into().unwrap()),
                nfft as i32
            );
            let cross = spectrum(&path, nfft);
            for &bin in tones {
                assert!(
                    cross[bin].norm() > 1e-4,
                    "lost line at {bin}: {:?}",
                    cross[bin]
                );
                assert!(
                    cross[bin].arg().abs() < 0.03,
                    "case={case} {mode} bin={bin} phase={} amplitude={}",
                    cross[bin].arg(),
                    cross[bin].norm()
                );
                if bits == 8 {
                    let expected = 12.0_f64.powi(2) / (2.0 * (256.0_f64.powi(2) - 1.0) / 12.0);
                    assert!((cross[bin].norm()/expected - 1.0).abs() < 0.025, "incorrect native FFT normalization: case={case} bin={bin} got={} expected={expected}", cross[bin].norm());
                }
            }
            if mode == "direct" {
                reference = cross;
            } else {
                for &bin in tones {
                    assert!(
                        (cross[bin] - reference[bin]).norm() < 2e-7,
                        "chunk/grid dependent phase or power at {bin}"
                    );
                }
            }
        }
    }
    fs::remove_dir_all(root).unwrap();
}
