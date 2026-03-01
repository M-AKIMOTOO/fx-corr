use crate::utils::DynError;
use erfa::aliases::{eraAnp, eraC2s, eraGmst06, eraGst06a, eraPmat06, eraRxp, eraS2c};

const C: f64 = 299792458.0; // Speed of light in m/s

// Fixed ECEF coordinates for Yamaguchi interferometer (meters)
pub const YAMAGU32_ECEF: [f64; 3] = [-3502544.587, 3950966.235, 3566381.192];
pub const YAMAGU34_ECEF: [f64; 3] = [-3502567.576, 3950885.734, 3566449.115];

// --- VLBI Geometric Definitions ---
// Celestial direction vector (s): Unit vector pointing from the observer towards the celestial source.
// Antenna position vector (b): Vector from the reference point (e.g., Earth's center) to the antenna.
// Geometric delay (tau): Time difference in signal arrival.
//   tau_i = -(b_i . s) / C
//   Positive tau_i means signal arrives late. Negative tau_i means signal arrives early.
// Baseline vector (B): Vector from antenna 1 to antenna 2 (B = b2 - b1).
// Baseline geometric delay (tau_baseline): Delay of antenna 2 relative to antenna 1.
//   tau_baseline = tau_2 - tau_1
// ----------------------------------

fn parse_packed_hhmmss(value: f64) -> Result<(f64, f64, f64), DynError> { 
    let abs = value.abs();
    let h = (abs / 10000.0).floor();
    let rem = abs - h * 10000.0;
    let m = (rem / 100.0).floor();
    let s = rem - m * 100.0;
    if !(0.0..60.0).contains(&m) || !(0.0..60.0).contains(&s) {
        return Err("Invalid hhmmss/ddmmss value".into());
    }
    Ok((h, m, s))
}

fn parse_hms_like(input: &str) -> Result<(f64, f64, f64), DynError> {
    let cleaned = input
        .replace('h', " ")
        .replace('m', " ")
        .replace('s', " ")
        .replace(':', " ");
    let parts: Vec<&str> = cleaned.split_whitespace().collect();
    match parts.len() {
        1 => {
            let packed = parts[0].parse::<f64>()?;
            parse_packed_hhmmss(packed)
        }
        3 => {
            let h = parts[0].parse::<f64>()?;
            let m = parts[1].parse::<f64>()?;
            let s = parts[2].parse::<f64>()?;
            if !(0.0..60.0).contains(&m) || !(0.0..60.0).contains(&s) {
                return Err("Invalid hms value".into());
            }
            Ok((h, m, s))
        }
        _ => Err("Invalid hms format".into()),
    }
}

fn parse_dms_like(input: &str) -> Result<(f64, f64, f64), DynError> {
    let cleaned = input
        .replace('d', " ")
        .replace('m', " ")
        .replace('s', " ")
        .replace('\'', " ")
        .replace('\"', " ")
        .replace(':', " ");
    let parts: Vec<&str> = cleaned.split_whitespace().collect();
    match parts.len() {
        1 => {
            let packed = parts[0].parse::<f64>()?;
            parse_packed_hhmmss(packed)
        }
        3 => {
            let d = parts[0].parse::<f64>()?;
            let m = parts[1].parse::<f64>()?;
            let s = parts[2].parse::<f64>()?;
            if !(0.0..60.0).contains(&m) || !(0.0..60.0).contains(&s) {
                return Err("Invalid dms value".into());
            }
            Ok((d, m, s))
        }
        _ => Err("Invalid dms format".into()),
    }
}

// Function to parse RA string to radians.
// Supports: hms markers, hh:mm:ss, hhmmss, or decimal degrees.
pub fn parse_ra(ra_str: &str) -> Result<f64, DynError> {
    let raw = ra_str.trim().to_lowercase();
    if raw.is_empty() {
        return Err("Empty RA".into());
    }

    let has_hms_marker = raw.contains('h')
        || raw.contains('m')
        || raw.contains('s')
        || raw.contains(':')
        || raw.contains(' ');

    let hours = if has_hms_marker {
        let (h, m, s) = parse_hms_like(&raw)?;
        h + m / 60.0 + s / 3600.0
    } else {
        let v = raw.parse::<f64>()?;
        if v.abs() >= 10000.0 {
            let (h, m, s) = parse_packed_hhmmss(v)?;
            h + m / 60.0 + s / 3600.0
        } else {
            // Decimal degrees
            return Ok(v.to_radians());
        }
    };

    Ok((hours * 15.0).to_radians())
}

// Function to parse Dec string to radians.
// Supports: dms markers, dd:mm:ss, ddmmss, or decimal degrees.
pub fn parse_dec(dec_str: &str) -> Result<f64, DynError> {
    let raw = dec_str.trim().to_lowercase();
    if raw.is_empty() {
        return Err("Empty Dec".into());
    }

    let sign = if raw.starts_with('-') { -1.0 } else { 1.0 };
    let stripped = raw.trim_start_matches(|c: char| c == '+' || c == '-');

    let has_dms_marker = stripped.contains('d')
        || stripped.contains('m')
        || stripped.contains('s')
        || stripped.contains('\'')
        || stripped.contains('\"')
        || stripped.contains(':')
        || stripped.contains(' ');

    let degrees = if has_dms_marker {
        let (d, m, s) = parse_dms_like(stripped)?;
        d + m / 60.0 + s / 3600.0
    } else {
        let v = stripped.parse::<f64>()?;
        if v.abs() >= 10000.0 {
            let (d, m, s) = parse_packed_hhmmss(v)?;
            d + m / 60.0 + s / 3600.0
        } else {
            return Ok((sign * v).to_radians());
        }
    };

    Ok((sign * degrees).to_radians())
}

// Parse an ISO 8601-like string and calculate MJD.
// Example: 2024-02-12T15:52:00Z
fn parse_iso_epoch_to_mjd(epoch_str: &str) -> Result<f64, DynError> {
    let normalized = epoch_str.replace('z', "Z");
    let parts: Vec<&str> = normalized.split(&['-', 'T', ':', 'Z'][..]).collect();
    if parts.len() < 6 {
        return Err("Invalid epoch format".into());
    }
    let year = parts[0].parse::<i32>()?;
    let month = parts[1].parse::<u32>()?;
    let day = parts[2].parse::<u32>()?;
    let hour = parts[3].parse::<u32>()?;
    let minute = parts[4].parse::<u32>()?;
    let second = parts[5].parse::<f64>()?;

    // Formula from https://en.wikipedia.org/wiki/Julian_day#Julian_day_number_calculation
    let (y, m) = if month <= 2 {
        (year - 1, month + 12)
    } else {
        (year, month)
    };

    let a = y / 100;
    let b = 2 - a + (a / 4);

    let jd_int = (365.25 * (y + 4716) as f64).floor() as i32
        + (30.6001 * (m + 1) as f64).floor() as i32
        + day as i32
        + b
        - 1524;

    let frac_day = (hour as f64 / 24.0) + (minute as f64 / 1440.0) + (second / 86400.0);

    // Julian day starts at noon; shift by -0.5 so 00:00 UTC maps correctly.
    let julian_day = jd_int as f64 - 0.5 + frac_day;

    Ok(julian_day - 2400000.5) // Convert to MJD
}

fn is_leap_year(year: i32) -> bool {
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
}

fn doy_to_month_day(year: i32, doy: u32) -> Result<(u32, u32), DynError> {
    let month_days = if is_leap_year(year) {
        [31u32, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        [31u32, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };
    let max_doy: u32 = month_days.iter().sum();
    if doy == 0 || doy > max_doy {
        return Err("DOY out of range".into());
    }
    let mut rem = doy;
    for (idx, md) in month_days.iter().enumerate() {
        if rem <= *md {
            return Ok(((idx + 1) as u32, rem));
        }
        rem -= *md;
    }
    Err("Failed to convert DOY to month/day".into())
}

fn parse_year_doy_epoch_to_mjd(epoch_str: &str) -> Result<f64, DynError> {
    // Supports formats like:
    //   YYYY/DDD
    //   YYYY/DDD HH:MM:SS
    //   YYYY/DDD HH:MM:SS.sss
    let mut it = epoch_str.split_whitespace();
    let date_part = it.next().ok_or("Invalid YYYY/DDD epoch: missing date")?;
    let time_part = it.next().unwrap_or("00:00:00");
    let (y_str, ddd_str) = date_part
        .split_once('/')
        .ok_or("Invalid YYYY/DDD epoch: expected '/'")?;
    let year = y_str.parse::<i32>()?;
    let doy = ddd_str.parse::<u32>()?;
    let (month, day) = doy_to_month_day(year, doy)?;
    let iso = format!("{year:04}-{month:02}-{day:02}T{time_part}Z");
    parse_iso_epoch_to_mjd(&iso)
}

// Parse epoch to MJD.
// Supports:
// - ISO datetime (e.g. 2024-02-12T15:52:00Z)
// - MJD numeric (e.g. 60350.0)
// - Year / fractional year (e.g. 2000, 2024.5)
pub fn parse_epoch_to_mjd(epoch_str: &str) -> Result<f64, DynError> {
    let raw = epoch_str.trim();
    if raw.is_empty() {
        return Err("Empty epoch".into());
    }
    if raw.contains('/') {
        return parse_year_doy_epoch_to_mjd(raw);
    }
    if raw.contains('-')
        || raw.contains('T')
        || raw.contains(':')
        || raw.contains('Z')
        || raw.contains('z')
    {
        return parse_iso_epoch_to_mjd(raw);
    }
    let numeric = if raw.starts_with('J') || raw.starts_with('j') {
        raw[1..].parse::<f64>()?
    } else {
        raw.parse::<f64>()?
    };

    if (40000.0..100000.0).contains(&numeric) {
        // Interpreted as MJD
        return Ok(numeric);
    }
    if (1800.0..3000.0).contains(&numeric) {
        // Interpreted as Julian year referenced to J2000.0
        return Ok(51544.5 + (numeric - 2000.0) * 365.25);
    }

    Err("Unsupported epoch format. Use ISO datetime, MJD, or year (e.g. 2000)".into())
}

fn mjd_to_jd_split(mjd: f64) -> (f64, f64) {
    // Keep MJD as the second part to preserve precision in downstream ERFA calls.
    (2400000.5, mjd)
}

/// Precess J2000 mean RA/Dec to mean equator/equinox of date.
/// Uses ERFA IAU 2006 precession matrix.
pub fn precess_j2000_to_mean_of_date(ra_j2000: f64, dec_j2000: f64, mjd: f64) -> (f64, f64) {
    let (date1, date2) = mjd_to_jd_split(mjd);
    let p_j2000 = eraS2c(ra_j2000, dec_j2000);
    let rbp = eraPmat06(date1, date2);
    let p_date = eraRxp(rbp, p_j2000);
    let (ra_date, dec_date) = eraC2s(p_date);
    (eraAnp(ra_date), dec_date)
}

// Calculate Greenwich Mean Sidereal Time (GMST) in radians from MJD
#[allow(dead_code)]
pub fn mjd_to_gmst(mjd: f64) -> f64 {
    let (ut1a, ut1b) = mjd_to_jd_split(mjd);
    // We do not have external UT1-UTC and TT-UTC here.
    // For delay tracking, use the same split for UT1 and TT in ERFA.
    eraAnp(eraGmst06(ut1a, ut1b, ut1a, ut1b))
}

fn source_vector_itrs(ra: f64, dec: f64, mjd: f64) -> [f64; 3] {
    // Mean-of-date equinox path:
    // source coordinates are precessed to date, then rotated by apparent sidereal time.
    let (ut1a, ut1b) = mjd_to_jd_split(mjd);
    let gst = eraAnp(eraGst06a(ut1a, ut1b, ut1a, ut1b));
    let ha = gst - ra;
    let (sd, cd) = dec.sin_cos();
    let (sh, ch) = ha.sin_cos();
    [cd * ch, -cd * sh, sd]
}

// Helper function to calculate single antenna delay
fn calculate_single_antenna_delay(
    ant_xyz: [f64; 3],
    ra: f64,  // radians
    dec: f64, // radians
    mjd: f64,
) -> f64 {
    let s = source_vector_itrs(ra, dec, mjd);
    -1.0 * (ant_xyz[0] * s[0] + ant_xyz[1] * s[1] + ant_xyz[2] * s[2]) / C
}

// Helper function to calculate baseline geometric delay directly.
// This avoids subtracting two large antenna delays for short baselines.
fn calculate_baseline_delay(
    ant1_xyz: [f64; 3],
    ant2_xyz: [f64; 3],
    ra: f64,
    dec: f64,
    mjd: f64,
) -> f64 {
    let s = source_vector_itrs(ra, dec, mjd);

    let bx = ant2_xyz[0] - ant1_xyz[0];
    let by = ant2_xyz[1] - ant1_xyz[1];
    let bz = ant2_xyz[2] - ant1_xyz[2];

    -1.0 * (bx * s[0] + by * s[1] + bz * s[2]) / C
}

#[allow(dead_code)]
pub fn calculate_antenna_delay_and_derivatives(
    ant_xyz: [f64; 3],
    ra: f64,  // radians
    dec: f64, // radians
    mjd: f64,
) -> (f64, f64, f64) {
    let dt_s = 1.0;
    let dt_mjd = dt_s / 86400.0;
    let delay_minus_dt = calculate_single_antenna_delay(ant_xyz, ra, dec, mjd - dt_mjd);
    let delay_t = calculate_single_antenna_delay(ant_xyz, ra, dec, mjd);
    let delay_plus_dt = calculate_single_antenna_delay(ant_xyz, ra, dec, mjd + dt_mjd);
    let rate = (delay_plus_dt - delay_minus_dt) / (2.0 * dt_s);
    let accel = (delay_plus_dt - 2.0 * delay_t + delay_minus_dt) / (dt_s * dt_s);
    (delay_t, rate, accel)
}

pub fn calculate_geometric_delay_and_derivatives(
    ant1_xyz: [f64; 3],
    ant2_xyz: [f64; 3],
    ra: f64,  // radians
    dec: f64, // radians
    mjd: f64,
) -> (f64, f64, f64, f64, f64) {
    // (delay_ant1, delay_ant2, geometric_delay, geometric_rate, geometric_accel)
    let dt_s = 1.0; // Time step for numerical differentiation (1 second)
    let dt_mjd = dt_s / 86400.0; // Convert dt to MJD units

    // Calculate delays at t-dt, t, t+dt
    let geom_delay_minus_dt = calculate_baseline_delay(ant1_xyz, ant2_xyz, ra, dec, mjd - dt_mjd);

    let delay1_t = calculate_single_antenna_delay(ant1_xyz, ra, dec, mjd);
    let delay2_t = calculate_single_antenna_delay(ant2_xyz, ra, dec, mjd);
    let geom_delay_t = calculate_baseline_delay(ant1_xyz, ant2_xyz, ra, dec, mjd);

    let geom_delay_plus_dt = calculate_baseline_delay(ant1_xyz, ant2_xyz, ra, dec, mjd + dt_mjd);

    // Calculate derivatives
    let geometric_rate = (geom_delay_plus_dt - geom_delay_minus_dt) / (2.0 * dt_s);
    let geometric_accel =
        (geom_delay_plus_dt - 2.0 * geom_delay_t + geom_delay_minus_dt) / (dt_s * dt_s);

    (
        delay1_t,
        delay2_t,
        geom_delay_t,
        geometric_rate,
        geometric_accel,
    )
}

#[cfg(test)]
mod tests {
    use super::{mjd_to_gmst, parse_dec, parse_epoch_to_mjd, parse_ra};

    #[test]
    fn parse_ra_supports_hhmmss_and_hms() {
        let ra_hms = parse_ra("16h42m58.8s").unwrap();
        let ra_packed = parse_ra("164258.8").unwrap();
        assert!((ra_hms - ra_packed).abs() < 1e-12);
    }

    #[test]
    fn parse_dec_supports_ddmmss_and_dms() {
        let dec_dms = parse_dec("+39d48m36.0s").unwrap();
        let dec_packed = parse_dec("+394836.0").unwrap();
        assert!((dec_dms - dec_packed).abs() < 1e-12);
    }

    #[test]
    fn parse_epoch_supports_default_year_and_mjd() {
        let j2000_mjd = parse_epoch_to_mjd("2000").unwrap();
        assert!((j2000_mjd - 51544.5).abs() < 1e-9);
        let mjd = parse_epoch_to_mjd("60350.25").unwrap();
        assert!((mjd - 60350.25).abs() < 1e-12);
    }

    #[test]
    fn parse_epoch_iso_is_utc_without_12h_offset() {
        let mjd = parse_epoch_to_mjd("2025-09-29T08:38:00Z").unwrap();
        assert!((mjd - 60947.35972222222).abs() < 1e-9);
    }

    #[test]
    fn parse_epoch_supports_year_doy_format() {
        let mjd = parse_epoch_to_mjd("2024/043 15:52:00").unwrap();
        // 2024-02-12 15:52:00 UTC
        assert!((mjd - 60352.66111111111).abs() < 1e-9);
    }

    #[test]
    fn gmst_progresses_at_earth_rotation_rate() {
        let mjd = 60352.66111;
        let gmst_0 = mjd_to_gmst(mjd);
        let gmst_1 = mjd_to_gmst(mjd + 1.0 / 86400.0);
        let mut d = gmst_1 - gmst_0;
        if d < 0.0 {
            d += 2.0 * std::f64::consts::PI;
        }
        // Sidereal angular speed ~ 7.2921159e-5 rad/s
        assert!(d > 7.2e-5 && d < 7.4e-5, "unexpected gmst delta rad/s: {d}");
    }
}
