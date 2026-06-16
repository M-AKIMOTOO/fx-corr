#!/usr/bin/env python3
"""Compare frinZ residual delay with an approximate solid Earth tide signal.

This is a diagnostic tool.  It uses a simple degree-2 elastic Earth tide model
for Sun and Moon displacement and projects the differential station motion onto
the source vector already dumped by yi-corr as phase_db{xyz}_cycles_per_m.
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import re
import sys
from pathlib import Path

import numpy as np

C_MPS = 299_792_458.0
RE_M = 6_378_136.6
GM_SUN_OVER_EARTH = 332_946.0487
GM_MOON_OVER_EARTH = 0.0123000371
H2 = 0.6078
L2 = 0.0847

DEFAULT_ANT1 = np.array([-3502544.587, 3950966.235, 3566381.192], dtype=float)
DEFAULT_ANT2 = np.array([-3961788.974, 3243597.492, 3790597.692], dtype=float)


def parse_epoch(value: str) -> dt.datetime:
    value = value.strip()
    m = re.fullmatch(r"(\d{4})/(\d{1,3})[ T](\d{1,2}):(\d{2}):(\d{2}(?:\.\d*)?)", value)
    if m:
        year = int(m.group(1))
        doy = int(m.group(2))
        hour = int(m.group(3))
        minute = int(m.group(4))
        sec = float(m.group(5))
        base = dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc) + dt.timedelta(days=doy - 1)
        return base + dt.timedelta(hours=hour, minutes=minute, seconds=sec)
    m = re.fullmatch(r"(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})", value)
    if m:
        year = int(m.group(1))
        doy = int(m.group(2))
        hour = int(m.group(3))
        minute = int(m.group(4))
        second = int(m.group(5))
        return (
            dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc)
            + dt.timedelta(days=doy - 1, hours=hour, minutes=minute, seconds=second)
        )
    iso = value.replace("Z", "+00:00")
    parsed = dt.datetime.fromisoformat(iso)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def epoch_from_delay_model_name(path: Path) -> dt.datetime | None:
    m = re.search(r"_(\d{4}\d{3}\d{6})\.tsv$", path.name)
    if not m:
        return None
    return parse_epoch(m.group(1))


def read_delay_model(path: Path) -> tuple[list[str], np.ndarray]:
    header: list[str] | None = None
    rows: list[list[float]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").strip().split()
                if "t_sec" in cols and "phase_dbx_cycles_per_m" in cols:
                    header = cols
                continue
            rows.append([float(x) for x in line.split()])
    if header is None:
        raise SystemExit(f"{path}: delay-model header lacks phase_db*_cycles_per_m columns")
    return header, np.array(rows, dtype=float)


def try_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def read_frinz_resdelay(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows: list[tuple[float, float]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 7:
                vals = [try_float(p) for p in parts]
                if all(v is not None for v in vals):
                    # frinZ --add plot TSV:
                    # t_sec amp_pct snr phase_deg noise_pct res_delay_sample res_rate_hz
                    rows.append((float(vals[0]), float(vals[5])))
                    continue
            if len(parts) < 10:
                continue
            # frinZ text rows: epoch label source len amp snr phase noise res-delay res-rate ...
            phase = try_float(parts[6])
            res_delay = try_float(parts[8])
            if phase is not None and res_delay is not None:
                rows.append((len(rows) * 10.0, res_delay))
                continue
            # add_plot TSV fallback: use likely numeric columns by name-independent layout.
            nums = [try_float(p) for p in parts]
            nums = [x for x in nums if x is not None]
            if len(nums) >= 9:
                rows.append((nums[0], nums[8]))
    if not rows:
        raise SystemExit(f"{path}: could not parse frinZ residual delay rows")
    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        raise ValueError("zero vector")
    return v / n


def solid_tide_displacement(station_ecef: np.ndarray, body_itrs_m: np.ndarray, gm_ratio: float) -> np.ndarray:
    rhat = unit(station_ecef)
    bhat = unit(body_itrs_m)
    dist = float(np.linalg.norm(body_itrs_m))
    dot = float(np.dot(rhat, bhat))
    scale = gm_ratio * (RE_M**4) / (dist**3)
    radial = 0.5 * H2 * (3.0 * dot * dot - 1.0) * rhat
    transverse = 3.0 * L2 * dot * (bhat - dot * rhat)
    return scale * (radial + transverse)


def body_itrs_vectors(times_utc: list[dt.datetime]) -> tuple[np.ndarray, np.ndarray]:
    from astropy import units as u
    from astropy.coordinates import ITRS, solar_system_ephemeris, get_body
    from astropy.time import Time

    with solar_system_ephemeris.set("builtin"):
        t = Time(times_utc)
        sun = get_body("sun", t).transform_to(ITRS(obstime=t))
        moon = get_body("moon", t).transform_to(ITRS(obstime=t))
    sun_xyz = sun.cartesian.xyz.to_value(u.m).T
    moon_xyz = moon.cartesian.xyz.to_value(u.m).T
    return np.array(sun_xyz, dtype=float), np.array(moon_xyz, dtype=float)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("delay_model", type=Path)
    ap.add_argument("--frinz", type=Path, help="frinZ text output or add_plot TSV for comparison")
    ap.add_argument("--epoch", help="UTC start epoch, e.g. 2025/302 08:15:00")
    ap.add_argument("--fs-hz", type=float, default=1.024e9)
    ap.add_argument("--obsfreq-hz", type=float, default=6.6e9)
    args = ap.parse_args()

    start = parse_epoch(args.epoch) if args.epoch else epoch_from_delay_model_name(args.delay_model)
    if start is None:
        raise SystemExit("could not infer epoch from filename; pass --epoch 'YYYY/DDD HH:MM:SS'")

    header, dm = read_delay_model(args.delay_model)
    idx = {name: i for i, name in enumerate(header)}
    t_sec = dm[:, idx["t_sec"]]
    phase_db = np.column_stack(
        [
            dm[:, idx["phase_dbx_cycles_per_m"]],
            dm[:, idx["phase_dby_cycles_per_m"]],
            dm[:, idx["phase_dbz_cycles_per_m"]],
        ]
    )
    source_itrs = phase_db * (C_MPS / args.obsfreq_hz)
    times = [start + dt.timedelta(seconds=float(t)) for t in t_sec]
    sun_itrs, moon_itrs = body_itrs_vectors(times)

    sun_delay = []
    moon_delay = []
    total_delay = []
    for i in range(len(times)):
        ds1_sun = solid_tide_displacement(DEFAULT_ANT1, sun_itrs[i], GM_SUN_OVER_EARTH)
        ds2_sun = solid_tide_displacement(DEFAULT_ANT2, sun_itrs[i], GM_SUN_OVER_EARTH)
        ds1_moon = solid_tide_displacement(DEFAULT_ANT1, moon_itrs[i], GM_MOON_OVER_EARTH)
        ds2_moon = solid_tide_displacement(DEFAULT_ANT2, moon_itrs[i], GM_MOON_OVER_EARTH)
        db_sun = ds2_sun - ds1_sun
        db_moon = ds2_moon - ds1_moon
        s = source_itrs[i]
        sun_samp = -float(np.dot(db_sun, s)) / C_MPS * args.fs_hz
        moon_samp = -float(np.dot(db_moon, s)) / C_MPS * args.fs_hz
        sun_delay.append(sun_samp)
        moon_delay.append(moon_samp)
        total_delay.append(sun_samp + moon_samp)

    sun_delay = np.array(sun_delay)
    moon_delay = np.array(moon_delay)
    total_delay = np.array(total_delay)
    total_rel = total_delay - total_delay[0]

    print(f"# epoch_utc {start.isoformat()}")
    print("# model approximate_degree2_solid_earth_tide_sun_moon h2=0.6078 l2=0.0847")
    print(f"# rows {len(t_sec)}")
    print(f"# solid_tide_total_span_sample {total_rel.min():+.9e} {total_rel.max():+.9e}")

    if args.frinz:
        ft, fres = read_frinz_resdelay(args.frinz)
        tide_i = np.interp(ft, t_sec, total_rel)
        fres_rel = fres - fres[0]
        if len(ft) >= 3:
            corr = float(np.corrcoef(fres_rel, tide_i)[0, 1])
        else:
            corr = float("nan")
        scale = float(np.dot(tide_i, fres_rel) / np.dot(tide_i, tide_i)) if np.dot(tide_i, tide_i) else float("nan")
        rms = float(np.sqrt(np.mean((fres_rel - scale * tide_i) ** 2)))
        print(f"# frinz_resdelay_rel_span_sample {fres_rel.min():+.9e} {fres_rel.max():+.9e}")
        print(f"# corr_resdelay_tide {corr:+.6f}")
        print(f"# fit_resdelay_rel ~= {scale:+.6f} * solid_tide_rel ; rms_sample {rms:.9e}")

    print("# t_sec\tsolid_tide_sample\tsun_sample\tmoon_sample\tsolid_tide_rel_sample")
    for t, tot, sun, moon, rel in zip(t_sec, total_delay, sun_delay, moon_delay, total_rel):
        print(f"{t:.6f}\t{tot:+.12e}\t{sun:+.12e}\t{moon:+.12e}\t{rel:+.12e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
