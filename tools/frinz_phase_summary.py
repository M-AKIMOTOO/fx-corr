#!/usr/bin/env python3
import argparse
import math
import re
import sys
from collections import defaultdict

ROW_RE = re.compile(
    r"^\s*(\d{4}/\d{3}\s+\d{2}:\d{2}:\d{2})\s+"  # epoch
    r"(\S+)\s+(\S+)\s+"                              # label source
    r"([+-]?\d+(?:\.\d+)?)\s+"                       # length
    r"([+-]?\d+(?:\.\d+)?)\s+"                       # amp
    r"([+-]?\d+(?:\.\d+)?)\s+"                       # snr
    r"([+-]?\d+(?:\.\d+)?)\s+"                       # phase deg
    r"([+-]?\d+(?:\.\d+)?)\s+"                       # noise
    r"([+-]?\d+(?:\.\d+)?)\s+"                       # res delay sample
    r"([+-]?\d+(?:\.\d+)?)"                           # res rate hz
)


def parse_hms(epoch):
    _, hms = epoch.split()
    h, m, s = [int(x) for x in hms.split(':')]
    return h * 3600 + m * 60 + s


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


def main():
    ap = argparse.ArgumentParser(description='Summarize frinZ phase/delay/rate rows by scan.')
    ap.add_argument('--fs', type=float, default=1024.0e6, help='sample rate Hz for delay drift conversion')
    ap.add_argument('--scan-gap', type=float, default=300.0, help='seconds separating scans')
    args = ap.parse_args()

    rows = []
    for line in sys.stdin:
        m = ROW_RE.match(line)
        if not m:
            continue
        epoch, label, source = m.group(1), m.group(2), m.group(3)
        rows.append({
            'epoch': epoch,
            't': parse_hms(epoch),
            'label': label,
            'source': source,
            'length': float(m.group(4)),
            'amp': float(m.group(5)),
            'snr': float(m.group(6)),
            'phase_deg': float(m.group(7)),
            'noise': float(m.group(8)),
            'res_delay_sample': float(m.group(9)),
            'res_rate_hz': float(m.group(10)),
        })

    if not rows:
        print('no frinZ rows parsed', file=sys.stderr)
        return 1

    scans = []
    cur = []
    last_t = None
    day_offset = 0
    for r in rows:
        t = r['t'] + day_offset
        if last_t is not None and t < last_t - 12 * 3600:
            day_offset += 24 * 3600
            t = r['t'] + day_offset
        if cur and t - last_t > args.scan_gap:
            scans.append(cur)
            cur = []
        rr = dict(r)
        rr['tabs'] = t
        cur.append(rr)
        last_t = t
    if cur:
        scans.append(cur)

    print('# scan start n amp_mean snr_mean phase0_deg phaseN_deg phase_slope_deg_s phase_rate_hz res_rate_mean_hz res_rate_slope_hz_s delay0_sample delayN_sample delay_slope_sample_s delay_slope_equiv_hz')
    for sc in scans:
        ts = [r['tabs'] for r in sc]
        phase = unwrap_deg([r['phase_deg'] for r in sc])
        delay = [r['res_delay_sample'] for r in sc]
        rates = [r['res_rate_hz'] for r in sc]
        amps = [r['amp'] for r in sc]
        snrs = [r['snr'] for r in sc]
        _, ph_slope = linfit(ts, phase)
        _, delay_slope = linfit(ts, delay)
        _, rate_slope = linfit(ts, rates)
        phase_rate_hz = ph_slope / 360.0
        delay_equiv_hz = delay_slope * 6.6e9 / args.fs
        print(
            f"{sc[0]['epoch']} {len(sc):2d} "
            f"{sum(amps)/len(amps):.6f} {sum(snrs)/len(snrs):.2f} "
            f"{phase[0]:+.3f} {phase[-1]:+.3f} "
            f"{ph_slope:+.6f} {phase_rate_hz:+.9f} "
            f"{sum(rates)/len(rates):+.9f} {rate_slope:+.9e} "
            f"{delay[0]:+.6f} {delay[-1]:+.6f} "
            f"{delay_slope:+.9f} {delay_equiv_hz:+.9f}"
        )
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
