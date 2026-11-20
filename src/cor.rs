use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use num_complex::Complex;

use crate::geom;
use crate::utils::DynError;

const FILE_HEADER_SIZE: usize = 256;
const SECTOR_HEADER_SIZE: usize = 128;
const EFFECTIVE_INTEG_TIME_OFFSET: usize = 112;
const NUM_SECTOR_OFFSET: u64 = 28;
const ST_TIME_OFFSET: usize = 0;
const ST_NSEC_OFFSET: usize = 4;
const ET_TIME_OFFSET: usize = 8;
const ET_NSEC_OFFSET: usize = 12;
#[derive(Clone, Copy, Debug)]
pub struct CorStation<'a> {
    pub name: &'a str,
    pub code: u8,
    pub ecef_m: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct CorHeaderConfig {
    pub sampling_speed_hz: i32,
    pub observing_frequency_hz: f64,
    pub fft_point: i32,
    pub number_of_sector_hint: i32,
    pub clock_reference_unix_sec: i64,
    pub source_name: String,
    pub source_ra_rad: f64,
    pub source_dec_rad: f64,
}

pub struct CorWriter {
    path: PathBuf,
    writer: BufWriter<File>,
    fft_half: usize,
    sectors_written: i32,
}

fn write_i32_le(buf: &mut [u8], offset: usize, value: i32) {
    if offset + 4 <= buf.len() {
        buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    }
}

fn write_u32_le(buf: &mut [u8], offset: usize, value: u32) {
    if offset + 4 <= buf.len() {
        buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    }
}

fn write_f32_le(buf: &mut [u8], offset: usize, value: f32) {
    if offset + 4 <= buf.len() {
        buf[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    }
}

fn write_f64_le(buf: &mut [u8], offset: usize, value: f64) {
    if offset + 8 <= buf.len() {
        buf[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
    }
}

fn write_fixed_ascii(buf: &mut [u8], offset: usize, len: usize, value: &str) {
    if offset + len > buf.len() {
        return;
    }
    let dst = &mut buf[offset..offset + len];
    dst.fill(0);
    let src = value.as_bytes();
    let n = src.len().min(len);
    dst[..n].copy_from_slice(&src[..n]);
}

fn split_sec_nsec(unix_sec: i64, extra_sec: f64) -> Result<(u32, u32), DynError> {
    if !extra_sec.is_finite() {
        return Err("non-finite sector time offset".into());
    }

    // Adding a small offset to a present-day Unix timestamp loses hundreds of ns.
    let whole = extra_sec.floor();
    let ns = ((extra_sec - whole) * 1_000_000_000.0).round() as u32;
    let sec = (whole as i128)
        .checked_add(i128::from(unix_sec))
        .and_then(|s| s.checked_add(i128::from(ns / 1_000_000_000)))
        .and_then(|s| u32::try_from(s).ok())
        .ok_or("sector sec out of u32 range for .cor sector header")?;
    Ok((sec, ns % 1_000_000_000))
}

fn build_file_header(
    cfg: &CorHeaderConfig,
    st1: CorStation<'_>,
    st2: CorStation<'_>,
) -> Result<[u8; FILE_HEADER_SIZE], DynError> {
    if cfg.fft_point <= 0 || (cfg.fft_point % 2) != 0 {
        return Err("cor header fft_point must be positive even integer".into());
    }
    if cfg.number_of_sector_hint < 0 {
        return Err("cor header number_of_sector_hint must be non-negative".into());
    }
    let clock_ref_sec_u32 = u32::try_from(cfg.clock_reference_unix_sec)
        .map_err(|_| "clock_reference_unix_sec out of u32 range for .cor header")?;

    let mut out = [0u8; FILE_HEADER_SIZE];

    // Legacy gico-compatible marker often seen in .cor files.
    out[0..4].copy_from_slice(&[0x83, 0xF9, 0xA2, 0x3E]);
    write_i32_le(&mut out, 4, 0x0103_0000);
    write_i32_le(&mut out, 8, 1);
    write_i32_le(&mut out, 12, cfg.sampling_speed_hz);
    write_f64_le(&mut out, 16, cfg.observing_frequency_hz);
    write_i32_le(&mut out, 24, cfg.fft_point);
    write_i32_le(&mut out, 28, cfg.number_of_sector_hint);

    write_fixed_ascii(&mut out, 32, 16, st1.name);
    write_f64_le(&mut out, 48, st1.ecef_m[0]);
    write_f64_le(&mut out, 56, st1.ecef_m[1]);
    write_f64_le(&mut out, 64, st1.ecef_m[2]);
    write_fixed_ascii(&mut out, 72, 1, &(st1.code as char).to_string());

    write_fixed_ascii(&mut out, 80, 16, st2.name);
    write_f64_le(&mut out, 96, st2.ecef_m[0]);
    write_f64_le(&mut out, 104, st2.ecef_m[1]);
    write_f64_le(&mut out, 112, st2.ecef_m[2]);
    write_fixed_ascii(&mut out, 120, 1, &(st2.code as char).to_string());

    write_fixed_ascii(&mut out, 128, 16, &cfg.source_name);
    write_f64_le(&mut out, 144, cfg.source_ra_rad);
    write_f64_le(&mut out, 152, cfg.source_dec_rad);
    // Keep legacy clock sec/nsec placeholders for both stations.
    // Offsets follow gico Header/Clock layout:
    //   st1 sec/nsec: [160..168], st1 poly: [168..208]
    //   st2 sec/nsec: [208..216], st2 poly: [216..256]
    write_u32_le(&mut out, 160, clock_ref_sec_u32);
    write_u32_le(&mut out, 164, 0);
    write_u32_le(&mut out, 208, clock_ref_sec_u32);
    write_u32_le(&mut out, 212, 0);
    // clock polynomial terms are left as zero.

    Ok(out)
}

impl CorWriter {
    pub fn create(
        path: &Path,
        cfg: &CorHeaderConfig,
        st1: CorStation<'_>,
        st2: CorStation<'_>,
    ) -> Result<Self, DynError> {
        let header = build_file_header(cfg, st1, st2)?;
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writer.write_all(&header)?;
        Ok(Self {
            path: path.to_path_buf(),
            writer,
            fft_half: (cfg.fft_point / 2) as usize,
            sectors_written: 0,
        })
    }

    pub fn write_sector_at(
        &mut self,
        timestamp_unix_sec: i64,
        start_offset_s: f64,
        duration_s: f64,
        effective_integ_time_s: f64,
        spectrum: &[Complex<f32>],
    ) -> Result<(), DynError> {
        if spectrum.len() != self.fft_half {
            return Err(format!(
                "sector spectrum length mismatch: expected {}, got {}",
                self.fft_half,
                spectrum.len()
            )
            .into());
        }

        let mut sector_header = [0u8; SECTOR_HEADER_SIZE];
        // Keep signed start time at offset 0 for compatibility with frinZ readers.
        let (st_sec, st_nsec) = split_sec_nsec(timestamp_unix_sec, start_offset_s)?;
        let st_i32 = i32::try_from(st_sec)
            .map_err(|_| "sector start sec out of i32 range for .cor sector header")?;
        write_i32_le(&mut sector_header, ST_TIME_OFFSET, st_i32);
        write_u32_le(&mut sector_header, ST_NSEC_OFFSET, st_nsec);

        let (et_sec, et_nsec) = split_sec_nsec(timestamp_unix_sec, start_offset_s + duration_s)?;
        write_u32_le(&mut sector_header, ET_TIME_OFFSET, et_sec);
        write_u32_le(&mut sector_header, ET_NSEC_OFFSET, et_nsec);
        write_f32_le(
            &mut sector_header,
            EFFECTIVE_INTEG_TIME_OFFSET,
            effective_integ_time_s as f32,
        );
        self.writer.write_all(&sector_header)?;
        for z in spectrum {
            self.writer.write_all(&z.re.to_le_bytes())?;
            self.writer.write_all(&z.im.to_le_bytes())?;
        }
        self.sectors_written += 1;
        Ok(())
    }

    pub fn finalize(mut self) -> Result<PathBuf, DynError> {
        self.writer.flush()?;
        {
            let file = self.writer.get_mut();
            file.seek(SeekFrom::Start(NUM_SECTOR_OFFSET))?;
            file.write_all(&self.sectors_written.to_le_bytes())?;
            file.flush()?;
        }
        Ok(self.path)
    }
}

pub fn epoch_seconds_parts(epoch: &str) -> Result<(i64, f64), DynError> {
    // Convert the minute boundary separately so a fractional second never
    // passes through a full Julian date or a floating-point Unix timestamp.
    let epoch = epoch.trim();
    if let Some((minute, seconds)) = epoch.rsplit_once(':') {
        let seconds: f64 = seconds.trim_end_matches(['Z', 'z']).parse()?;
        if !(0.0..61.0).contains(&seconds) {
            return Err("invalid epoch seconds".into());
        }
        let mjd = geom::parse_epoch_to_mjd(&format!("{minute}:00"))?;
        let base = ((mjd - 40587.0) * 86400.0).round() as i64;
        return Ok((base + seconds.trunc() as i64, seconds.fract()));
    }
    let mjd = geom::parse_epoch_to_mjd(epoch)?;
    let seconds = (mjd - mjd.floor()) * 86400.0;
    Ok((
        (mjd.floor() as i64 - 40587) * 86400 + seconds.trunc() as i64,
        seconds.fract(),
    ))
}

pub fn epoch_to_unix_seconds(epoch: &str) -> Result<i64, DynError> {
    let (sec, fraction) = epoch_seconds_parts(epoch)?;
    Ok(sec + fraction.round() as i64)
}

pub fn unix_seconds_to_yyyydddhhmmss(unix_seconds: i64) -> Result<String, DynError> {
    let ts = i64::clamp(unix_seconds, i64::from(i32::MIN), i64::from(i32::MAX)) as libc::time_t;
    let mut tm_out = std::mem::MaybeUninit::<libc::tm>::uninit();
    let ptr = unsafe { libc::gmtime_r(&ts, tm_out.as_mut_ptr()) };
    if ptr.is_null() {
        return Err("gmtime_r failed while formatting epoch".into());
    }
    let tm = unsafe { tm_out.assume_init() };
    let year = tm.tm_year + 1900;
    let doy = tm.tm_yday + 1;
    let month_day_time = format!("{:02}{:02}{:02}", tm.tm_hour, tm.tm_min, tm.tm_sec);
    Ok(format!("{:04}{:03}{}", year, doy, month_day_time))
}

pub fn epoch_to_yyyydddhhmmss(epoch: &str) -> Result<(i64, String), DynError> {
    let unix = epoch_to_unix_seconds(epoch)?;
    let tag = unix_seconds_to_yyyydddhhmmss(unix)?;
    Ok((unix, tag))
}
