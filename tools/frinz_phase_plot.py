#!/usr/bin/env python3
import argparse
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
    h, m, s = [int(x) for x in hms.split(":")]
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
        return float("nan"), float("nan")
    x0 = xs[0]
    xx = [x - x0 for x in xs]
    sx = sum(xx)
    sy = sum(ys)
    sxx = sum(x * x for x in xx)
    sxy = sum(x * y for x, y in zip(xx, ys))
    den = n * sxx - sx * sx
    if den == 0.0:
        return float("nan"), float("nan")
    slope = (n * sxy - sx * sy) / den
    intercept = (sy - slope * sx) / n
    return intercept, slope


def read_lines(paths):
    if not paths:
        yield from sys.stdin
        return
    for path in paths:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            yield from f


def parse_rows(lines):
    rows = []
    for line in lines:
        m = ROW_RE.match(line)
        if not m:
            continue
        epoch = m.group(1)
        rows.append(
            {
                "epoch": epoch,
                "t": parse_hms(epoch),
                "label": m.group(2),
                "source": m.group(3),
                "length": float(m.group(4)),
                "amp": float(m.group(5)),
                "snr": float(m.group(6)),
                "phase_deg": float(m.group(7)),
                "noise": float(m.group(8)),
                "res_delay_sample": float(m.group(9)),
                "res_rate_hz": float(m.group(10)),
            }
        )
    return rows


def split_scans(rows, scan_gap):
    scans = []
    cur = []
    last_t = None
    day_offset = 0
    for r in rows:
        t = r["t"] + day_offset
        if last_t is not None and t < last_t - 12 * 3600:
            day_offset += 24 * 3600
            t = r["t"] + day_offset
        if cur and t - last_t > scan_gap:
            scans.append(cur)
            cur = []
        rr = dict(r)
        rr["tabs"] = t
        cur.append(rr)
        last_t = t
    if cur:
        scans.append(cur)
    return scans


def main():
    ap = argparse.ArgumentParser(
        description="Plot wrapped and unwrapped frinZ phase rows."
    )
    ap.add_argument("files", nargs="*", help="frinZ output text files; stdin if omitted")
    ap.add_argument("-o", "--output", default="frinz_phase.png", help="output PNG/PDF/SVG")
    ap.add_argument("--scan-gap", type=float, default=300.0, help="seconds separating scans")
    ap.add_argument("--title", default=None, help="figure title")
    ap.add_argument("--no-lines", action="store_true", help="plot points only")
    ap.add_argument("--dpi", type=int, default=150, help="raster output DPI")
    args = ap.parse_args()

    rows = parse_rows(read_lines(args.files))
    if not rows:
        print("no frinZ rows parsed", file=sys.stderr)
        return 1

    scans = split_scans(rows, args.scan_gap)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, constrained_layout=True)
    ax_wrap, ax_unwrap = axes

    for idx, sc in enumerate(scans, 1):
        t0 = sc[0]["tabs"]
        xs = [r["tabs"] - t0 for r in sc]
        wrapped = [r["phase_deg"] for r in sc]
        unwrapped = unwrap_deg(wrapped)
        _, slope = linfit([r["tabs"] for r in sc], unwrapped)
        rate_hz = slope / 360.0 if math.isfinite(slope) else float("nan")
        label = f"{sc[0]['epoch']}  ({rate_hz:+.6f} Hz)"
        style = "o" if args.no_lines else "o-"
        ax_wrap.plot(xs, wrapped, style, ms=3, lw=1, label=label)
        ax_unwrap.plot(xs, unwrapped, style, ms=3, lw=1, label=label)

    title = args.title or "frinZ phase"
    ax_wrap.set_title(title)
    ax_wrap.set_ylabel("wrapped phase [deg]")
    ax_unwrap.set_ylabel("unwrapped phase [deg]")
    ax_unwrap.set_xlabel("time from scan start [s]")
    ax_wrap.grid(True, alpha=0.3)
    ax_unwrap.grid(True, alpha=0.3)
    ax_wrap.legend(fontsize="small", loc="best")
    ax_unwrap.legend(fontsize="small", loc="best")

    out = Path(args.output)
    fig.savefig(out, dpi=args.dpi)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
