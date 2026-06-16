#!/usr/bin/env python3
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description="Print IERS EOP values for yi-corr.")
    ap.add_argument("epoch", help="UTC epoch, e.g. 2025-10-29T08:15:00")
    ap.add_argument(
        "--allow-stale",
        action="store_true",
        help="allow stale/predictive local IERS_Auto values without downloading",
    )
    ap.add_argument(
        "--download",
        action="store_true",
        help="allow astropy to download fresh IERS data",
    )
    args = ap.parse_args()

    try:
        from astropy import units as u
        from astropy.time import Time
        from astropy.utils import iers
    except Exception as e:
        print(f"astropy is required: {e}", file=sys.stderr)
        return 1

    iers.conf.auto_download = bool(args.download)
    if args.allow_stale:
        iers.conf.auto_max_age = None

    t = Time(args.epoch, scale="utc")
    tab = iers.IERS_Auto.open()
    dut1 = t.get_delta_ut1_utc(iers_table=tab)
    xp, yp = tab.pm_xy(t)
    dut1_s = dut1.to_value(u.s) if hasattr(dut1, "to_value") else float(dut1)
    xp_arcsec = xp.to_value(u.arcsec)
    yp_arcsec = yp.to_value(u.arcsec)
    tt_utc_s = (t.tai.jd - t.utc.jd) * 86400.0 + 32.184

    print(f"dut1_s {dut1_s:.12f}")
    print(f"tt_utc_s {tt_utc_s:.12f}")
    print(f"xp_arcsec {xp_arcsec:.12f}")
    print(f"yp_arcsec {yp_arcsec:.12f}")
    print(
        f"<eop><dut1>{dut1_s:.12f}</dut1><tt-utc>{tt_utc_s:.12f}</tt-utc>"
        f"<xp>{xp_arcsec:.12f}</xp><yp>{yp_arcsec:.12f}</yp></eop>"
    )
    print(
        "YI_DUT1_S={:.12f} YI_TT_UTC_S={:.12f} "
        "YI_XP_ARCSEC={:.12f} YI_YP_ARCSEC={:.12f}".format(
            dut1_s, tt_utc_s, xp_arcsec, yp_arcsec
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
