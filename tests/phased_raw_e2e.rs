use std::fs;
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

    let script = r#"#!/bin/sh
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
if [ "$mode" = peak ]; then
  i=1
  while [ "$i" -le "$loops" ]; do
    delay=+0.250000
    if [ "$i" -eq 7 ]; then delay=+99.000000; fi
    printf '%s\n' "2000/001 00:00:00 L GAIN 1.00 0.100000 20.0 +0.0 0.010000 $delay +0.000000 1 2 3 4 5 6 51544.0 - False False"
    i=$((i + 1))
  done
  exit 0
fi
parent=$(dirname "$input")
stem=$(basename "$input" .cor)
out="$parent/frinZ/acel_search"
mkdir -p "$out"
printf '%s\n' '# Fitted: y = 1.000000e-6 * x^2 + 2.000000e-3 * x + 2.000000e-1' > "$out/${stem}_step1_quadric.txt"
printf '%s\n' '# Corrected Acel (Hz/s): 3.183099e-7 (from x^2 / PI)' >> "$out/${stem}_step1_quadric.txt"
printf '%s\n' '# Corrected Rate (Hz): 3.183099e-4 (from x / (2 * PI))' >> "$out/${stem}_step1_quadric.txt"
"#;
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
    assert!(meta.contains("software_version=3.5.3"));
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
    assert!(stdout.contains("[phasecal] stage 4/7: run frinZ --search acel"));
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
    assert!(solution.contains("format=yi-phasedarray-gain-phasecal-v2"));
    assert!(solution.contains("peak_delay_count=24"));
    assert!(solution.contains("peak_delay_samples_median=+2.50000000000000000e-1"));
    assert!(solution.contains("peak_delay_samples_max=+9.90000000000000000e1"));
    assert!(solution.contains("group_delay_s=+3.05175781250000000e-5"));
    assert!(solution.contains("phase_c0_rad=+2.00000000000000011e-1"));
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
