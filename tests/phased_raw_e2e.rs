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
    assert!(meta.contains("software_version=3.2.4"));
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
        stdout.contains("Phased scan-length preservation"),
        "missing boundary-padding diagnostic:\n{stdout}"
    );
    let phased_raw = phased_dir.join(format!("ARRAY_{tag}.raw"));
    assert_eq!(
        fs::metadata(phased_raw).unwrap().len(),
        raw.len() as u64,
        "virtual-station raw must retain the nominal one-second scan length"
    );

    let _ = fs::remove_dir_all(root);
}
