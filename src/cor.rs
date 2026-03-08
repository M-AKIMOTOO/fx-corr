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
const GEO1_SEC_OFFSET: usize = 16;
const GEO2_SEC_OFFSET: usize = 64;
const AMP0_OFFSET: usize = 116;
const AMP1_OFFSET: usize = 120;
const PHS0_OFFSET: usize = 124;
const PHS1_OFFSET: usize = 126;

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

#[derive(Clone, Copy, Debug, Default)]
pub struct CorClockModel {
    pub sec: u32,
    pub nsec: u32,
    pub delay: f64,
    pub rate: f64,
    pub acel: f64,
    pub jerk: f64,
    pub snap: f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct CorSectorModel {
    pub station1: CorClockModel,
    pub station2: CorClockModel,
    pub amp: [f32; 2],
    pub phs: [i16; 2],
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

fn write_i16_le(buf: &mut [u8], offset: usize, value: i16) {
    if offset + 2 <= buf.len() {
        buf[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
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

    let mut sec = unix_sec as f64 + extra_sec;
    let mut nsec = sec.fract() * 1_000_000_000.0;
    sec = sec.floor();

    if nsec < 0.0 {
        nsec += 1_000_000_000.0;
        sec -= 1.0;
    } else if nsec >= 1_000_000_000.0 {
        nsec -= 1_000_000_000.0;
        sec += 1.0;
    }

    let sec_i64 = sec as i64;
    let nsec_u32 = nsec.round() as i64;
    let (sec_i64, nsec_u32) = if nsec_u32 >= 1_000_000_000 {
        (sec_i64 + 1, 0)
    } else {
        (sec_i64, nsec_u32)
    };
    let sec_u32 = u32::try_from(sec_i64)
        .map_err(|_| "sector sec out of u32 range for .cor sector header")?;
    let nsec_u32 = u32::try_from(nsec_u32)
        .map_err(|_| "sector nsec out of u32 range for .cor sector header")?;
    Ok((sec_u32, nsec_u32))
}

fn write_clock_model(buf: &mut [u8], base_offset: usize, model: &CorClockModel) {
    write_u32_le(buf, base_offset, model.sec);
    write_u32_le(buf, base_offset + 4, model.nsec);
    write_f64_le(buf, base_offset + 8, model.delay);
    write_f64_le(buf, base_offset + 16, model.rate);
    write_f64_le(buf, base_offset + 24, model.acel);
    write_f64_le(buf, base_offset + 32, model.jerk);
    write_f64_le(buf, base_offset + 40, model.snap);
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
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        let header = build_file_header(cfg, st1, st2)?;
        writer.write_all(&header)?;
        Ok(Self {
            path: path.to_path_buf(),
            writer,
            fft_half: (cfg.fft_point / 2) as usize,
            sectors_written: 0,
        })
    }

    #[allow(dead_code)]
    pub fn write_sector(
        &mut self,
        timestamp_unix_sec: i64,
        effective_integ_time_s: f32,
        spectrum: &[Complex<f32>],
    ) -> Result<(), DynError> {
        self.write_sector_with_model(timestamp_unix_sec, effective_integ_time_s, spectrum, None)
    }

    pub fn write_sector_with_model(
        &mut self,
        timestamp_unix_sec: i64,
        effective_integ_time_s: f32,
        spectrum: &[Complex<f32>],
        model: Option<CorSectorModel>,
    ) -> Result<(), DynError> {
        if spectrum.len() != self.fft_half {
            return Err(format!(
                "sector spectrum length mismatch: expected {}, got {}",
                self.fft_half,
                spectrum.len()
            )
            .into());
        }
        let ts_i32 = i32::try_from(timestamp_unix_sec)
            .map_err(|_| "timestamp out of i32 range for .cor sector header")?;

        let mut sector_header = [0u8; SECTOR_HEADER_SIZE];
        // Keep signed start time at offset 0 for compatibility with frinZ readers.
        write_i32_le(&mut sector_header, ST_TIME_OFFSET, ts_i32);
        let (_st_sec, st_nsec) = split_sec_nsec(timestamp_unix_sec, 0.0)?;
        let integ_sec = if effective_integ_time_s.is_finite() && effective_integ_time_s > 0.0 {
            effective_integ_time_s as f64
        } else {
            1.0
        };
        let (et_sec, et_nsec) = split_sec_nsec(timestamp_unix_sec, integ_sec)?;
        write_u32_le(&mut sector_header, ST_NSEC_OFFSET, st_nsec);
        write_u32_le(&mut sector_header, ET_TIME_OFFSET, et_sec);
        write_u32_le(&mut sector_header, ET_NSEC_OFFSET, et_nsec);
        if let Some(model) = model {
            write_clock_model(&mut sector_header, GEO1_SEC_OFFSET, &model.station1);
            write_clock_model(&mut sector_header, GEO2_SEC_OFFSET, &model.station2);
            write_f32_le(&mut sector_header, AMP0_OFFSET, model.amp[0]);
            write_f32_le(&mut sector_header, AMP1_OFFSET, model.amp[1]);
            write_i16_le(&mut sector_header, PHS0_OFFSET, model.phs[0]);
            write_i16_le(&mut sector_header, PHS1_OFFSET, model.phs[1]);
        }
        write_f32_le(
            &mut sector_header,
            EFFECTIVE_INTEG_TIME_OFFSET,
            effective_integ_time_s,
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

fn mjd_to_unix_seconds(mjd: f64) -> i64 {
    ((mjd - 40587.0) * 86400.0).round() as i64
}

pub fn epoch_to_unix_seconds(epoch: &str) -> Result<i64, DynError> {
    let mjd = geom::parse_epoch_to_mjd(epoch)?;
    Ok(mjd_to_unix_seconds(mjd))
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
