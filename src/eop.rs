use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::utils::DynError;

#[derive(Clone, Debug)]
pub struct EopValues {
    pub dut1_s: f64,
    pub tt_utc_s: f64,
    pub xp_arcsec: f64,
    pub yp_arcsec: f64,
    pub source: String,
    pub mjd0: f64,
    pub mjd1: f64,
    pub kind0: EopKind,
    pub kind1: EopKind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EopKind {
    BulletinA,
    BulletinB,
}

impl EopKind {
    pub fn as_str(self) -> &'static str {
        match self {
            EopKind::BulletinA => "Bulletin A",
            EopKind::BulletinB => "Bulletin B",
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct EopRow {
    mjd: f64,
    xp_a_arcsec: Option<f64>,
    yp_a_arcsec: Option<f64>,
    dut1_a_s: Option<f64>,
    xp_b_arcsec: Option<f64>,
    yp_b_arcsec: Option<f64>,
    dut1_b_s: Option<f64>,
}

impl EopRow {
    fn selected(self) -> Option<(f64, f64, f64, EopKind)> {
        match (self.xp_b_arcsec, self.yp_b_arcsec, self.dut1_b_s) {
            (Some(xp), Some(yp), Some(dut1)) => Some((xp, yp, dut1, EopKind::BulletinB)),
            _ => match (self.xp_a_arcsec, self.yp_a_arcsec, self.dut1_a_s) {
                (Some(xp), Some(yp), Some(dut1)) => Some((xp, yp, dut1, EopKind::BulletinA)),
                _ => None,
            },
        }
    }
}

fn field(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end).unwrap_or("").trim()
}

fn opt_f64(line: &str, start: usize, end: usize) -> Option<f64> {
    let raw = field(line, start, end);
    if raw.is_empty() {
        None
    } else {
        raw.parse::<f64>().ok()
    }
}

fn parse_finals2000a_line(line: &str) -> Option<EopRow> {
    let mjd = opt_f64(line, 7, 15)?;
    Some(EopRow {
        mjd,
        xp_a_arcsec: opt_f64(line, 18, 27),
        yp_a_arcsec: opt_f64(line, 37, 46),
        dut1_a_s: opt_f64(line, 58, 68),
        xp_b_arcsec: opt_f64(line, 134, 144),
        yp_b_arcsec: opt_f64(line, 144, 154),
        dut1_b_s: opt_f64(line, 154, 165),
    })
}

fn tt_minus_utc_s(mjd_utc: f64) -> f64 {
    const LEAPS: &[(f64, f64)] = &[
        (41317.0, 10.0),
        (41499.0, 11.0),
        (41683.0, 12.0),
        (42048.0, 13.0),
        (42413.0, 14.0),
        (42778.0, 15.0),
        (43144.0, 16.0),
        (43509.0, 17.0),
        (43874.0, 18.0),
        (44239.0, 19.0),
        (44786.0, 20.0),
        (45151.0, 21.0),
        (45516.0, 22.0),
        (46247.0, 23.0),
        (47161.0, 24.0),
        (47892.0, 25.0),
        (48257.0, 26.0),
        (48804.0, 27.0),
        (49169.0, 28.0),
        (49534.0, 29.0),
        (50083.0, 30.0),
        (50630.0, 31.0),
        (51179.0, 32.0),
        (53736.0, 33.0),
        (54832.0, 34.0),
        (56109.0, 35.0),
        (57204.0, 36.0),
        (57754.0, 37.0),
    ];
    let mut tai_utc = 10.0;
    for (mjd_eff, value) in LEAPS {
        if mjd_utc >= *mjd_eff {
            tai_utc = *value;
        } else {
            break;
        }
    }
    tai_utc + 32.184
}

pub fn interpolate_finals2000a(path: &Path, mjd_utc: f64) -> Result<EopValues, DynError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(row) = parse_finals2000a_line(&line) {
            if row.selected().is_some() {
                rows.push(row);
            }
        }
    }
    if rows.len() < 2 {
        return Err(format!(
            "EOP file '{}' does not contain enough finals2000A rows",
            path.display()
        )
        .into());
    }
    rows.sort_by(|a, b| {
        a.mjd
            .partial_cmp(&b.mjd)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let first = rows.first().unwrap().mjd;
    let last = rows.last().unwrap().mjd;
    if mjd_utc < first || mjd_utc > last {
        return Err(format!(
            "EOP MJD {:.6} is outside '{}' range {:.1}..{:.1}",
            mjd_utc,
            path.display(),
            first,
            last
        )
        .into());
    }

    for pair in rows.windows(2) {
        let a = pair[0];
        let b = pair[1];
        if mjd_utc < a.mjd || mjd_utc > b.mjd {
            continue;
        }
        let (xp0, yp0, dut10, kind0) = a.selected().unwrap();
        let (xp1, yp1, dut11, kind1) = b.selected().unwrap();
        let w = if b.mjd == a.mjd {
            0.0
        } else {
            (mjd_utc - a.mjd) / (b.mjd - a.mjd)
        };
        let lerp = |x0: f64, x1: f64| x0 + (x1 - x0) * w;
        return Ok(EopValues {
            dut1_s: lerp(dut10, dut11),
            tt_utc_s: tt_minus_utc_s(mjd_utc),
            xp_arcsec: lerp(xp0, xp1),
            yp_arcsec: lerp(yp0, yp1),
            source: path.display().to_string(),
            mjd0: a.mjd,
            mjd1: b.mjd,
            kind0,
            kind1,
        });
    }

    Err(format!(
        "EOP MJD {:.6} did not bracket any rows in '{}'",
        mjd_utc,
        path.display()
    )
    .into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn put(buf: &mut [u8], start: usize, end: usize, value: &str) {
        let width = end - start;
        let text = format!("{value:>width$}");
        buf[start..end].copy_from_slice(text.as_bytes());
    }

    fn finals_line(mjd: f64, xp_a: f64, yp_a: f64, dut1_a: f64) -> String {
        let mut buf = vec![b' '; 170];
        put(&mut buf, 7, 15, &format!("{mjd:.0}"));
        put(&mut buf, 18, 27, &format!("{xp_a:.6}"));
        put(&mut buf, 37, 46, &format!("{yp_a:.6}"));
        put(&mut buf, 58, 68, &format!("{dut1_a:.7}"));
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn tt_utc_after_2017_is_69_184() {
        assert!((tt_minus_utc_s(60977.0) - 69.184).abs() < 1e-12);
    }

    #[test]
    fn interpolates_bulletin_a_rows() {
        let path = std::env::temp_dir().join(format!(
            "fx_corr_eop_test_{}_{}.data",
            std::process::id(),
            60977
        ));
        std::fs::write(
            &path,
            format!(
                "{}
{}
",
                finals_line(60977.0, 0.10, 0.30, 0.11),
                finals_line(60978.0, 0.20, 0.40, 0.13)
            ),
        )
        .unwrap();
        let eop = interpolate_finals2000a(&path, 60977.5).unwrap();
        let _ = std::fs::remove_file(&path);
        assert!((eop.xp_arcsec - 0.15).abs() < 1e-12);
        assert!((eop.yp_arcsec - 0.35).abs() < 1e-12);
        assert!((eop.dut1_s - 0.12).abs() < 1e-12);
        assert_eq!(eop.kind0, EopKind::BulletinA);
        assert_eq!(eop.kind1, EopKind::BulletinA);
    }
}
