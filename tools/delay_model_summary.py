#!/usr/bin/env python3
import argparse
import csv
import math
import sys
from pathlib import Path


def linfit(xs, ys):
    n = len(xs)
    if n < 2:
        return float('nan'), float('nan')
    x0 = xs[0]
    xx = [x - x0 for x in xs]
    sx = sum(xx)
    sy = sum(ys)
    sxx = sum(x * x for x in xx)
    sxy = sum(x * y for x, y in zip(xx, ys))
    den = n * sxx - sx * sx
    if den == 0.0:
        return float('nan'), float('nan')
    slope = (n * sxy - sx * sy) / den
    intercept = (sy - slope * sx) / n
    return intercept, slope


def unwrap_deg(phases):
    out = []
    prev = None
    off = 0.0
    for p in phases:
        q = p + off
        if prev is not None:
            while q - prev > 180.0:
                off -= 360.0
                q = p + off
            while q - prev < -180.0:
                off += 360.0
                q = p + off
        out.append(q)
        prev = q
    return out


def read_rows(path):
    with open(path, newline='') as f:
        header = None
        rows = []
        for line in f:
            if line.startswith('#'):
                if line.startswith('# t_sec'):
                    header = line[2:].strip().split('\t')
                continue
            if not line.strip():
                continue
            if header is None:
                raise ValueError(f'{path}: header not found')
            vals = line.rstrip('\n').split('\t')
            rows.append({k: float(v) for k, v in zip(header, vals)})
        return rows


def main():
    ap = argparse.ArgumentParser(description='Summarize yi-corr delay_model_*.tsv files.')
    ap.add_argument('files', nargs='+')
    ap.add_argument('--fs', type=float, default=1024.0e6, help='sample rate Hz')
    ap.add_argument('--obsfreq', type=float, default=6.6e9, help='reference observing frequency Hz')
    ap.add_argument('--duration', type=float, default=120.0, help='seconds to summarize from each file')
    args = ap.parse_args()

    print('# file n model_rate_mean_Hz model_rate_start_Hz model_rate_end_Hz model_rate_slope_Hz_s accel_mean_sample_s2 residual0_sample residualN_sample residual_slope_sample_s dump_phase_rate_Hz carrier_phase_rate_Hz carrier_hz')
    for name in args.files:
        path = Path(name)
        rows = [r for r in read_rows(path) if r['t_sec'] <= args.duration]
        if len(rows) < 2:
            continue
        ts = [r['t_sec'] for r in rows]
        rate_hz = [-r['model_rate_sample_s'] * args.obsfreq / args.fs for r in rows]
        accel = [r['model_accel_sample_s2'] for r in rows]
        residual = [r['residual_sample'] for r in rows]
        dump_phase = unwrap_deg([r['fringe_phase_deg_wrapped'] for r in rows])
        carrier_phase = unwrap_deg([r.get('carrier_phase_deg_wrapped', r['fringe_phase_deg_wrapped']) for r in rows])
        carrier_hz = rows[0].get('carrier_hz', args.obsfreq)
        _, rate_slope = linfit(ts, rate_hz)
        _, residual_slope = linfit(ts, residual)
        _, dump_phase_slope = linfit(ts, dump_phase)
        _, carrier_phase_slope = linfit(ts, carrier_phase)
        print(
            f"{path.name} {len(rows):3d} "
            f"{sum(rate_hz)/len(rate_hz):+.9f} {rate_hz[0]:+.9f} {rate_hz[-1]:+.9f} {rate_slope:+.9e} "
            f"{sum(accel)/len(accel):+.9f} "
            f"{residual[0]:+.6f} {residual[-1]:+.6f} {residual_slope:+.9f} "
            f"{dump_phase_slope/360.0:+.9f} "
            f"{carrier_phase_slope/360.0:+.9f} {carrier_hz:+.9f}"
        )
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
