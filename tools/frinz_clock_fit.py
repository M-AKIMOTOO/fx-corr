#!/usr/bin/env python3
import argparse
import math
import re
import sys

ROW_RE = re.compile(
    r"^\s*(\d{4}/\d{3}\s+\d{2}:\d{2}:\d{2})\s+"
    r"(\S+)\s+(\S+)\s+"
    r"([+-]?\d+(?:\.\d+)?)\s+"
    r"([+-]?\d+(?:\.\d+)?)\s+"
    r"([+-]?\d+(?:\.\d+)?)\s+"
    r"([+-]?\d+(?:\.\d+)?)\s+"
    r"([+-]?\d+(?:\.\d+)?)\s+"
    r"([+-]?\d+(?:\.\d+)?)\s+"
    r"([+-]?\d+(?:\.\d+)?)"
)


def hms_seconds(epoch):
    _, hms = epoch.split()
    h, m, s = [int(x) for x in hms.split(":")]
    return h * 3600 + m * 60 + s


def parse_rows(paths):
    rows = []
    streams = [sys.stdin] if not paths else [open(p) for p in paths]
    try:
        for stream in streams:
            for line in stream:
                m = ROW_RE.match(line)
                if not m:
                    continue
                epoch = m.group(1)
                rows.append(
                    {
                        "epoch": epoch,
                        "tod": hms_seconds(epoch),
                        "label": m.group(2),
                        "source": m.group(3),
                        "length": float(m.group(4)),
                        "amp": float(m.group(5)),
                        "snr": float(m.group(6)),
                        "phase_deg": float(m.group(7)),
                        "noise": float(m.group(8)),
                        "delay_sample": float(m.group(9)),
                        "rate_hz": float(m.group(10)),
                    }
                )
    finally:
        for stream in streams:
            if stream is not sys.stdin:
                stream.close()
    rows.sort(key=lambda r: r["tod"])
    day_offset = 0
    prev = None
    for r in rows:
        t = r["tod"] + day_offset
        if prev is not None and t < prev - 12 * 3600:
            day_offset += 24 * 3600
            t = r["tod"] + day_offset
        r["tabs"] = t
        prev = t
    return rows


def solve_linear(a, b):
    n = len(b)
    a = [row[:] for row in a]
    b = b[:]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-30:
            raise ValueError("singular fit matrix")
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            b[col], b[pivot] = b[pivot], b[col]
        inv = 1.0 / a[col][col]
        for j in range(col, n):
            a[col][j] *= inv
        b[col] *= inv
        for r in range(n):
            if r == col:
                continue
            f = a[r][col]
            if f == 0.0:
                continue
            for j in range(col, n):
                a[r][j] -= f * a[col][j]
            b[r] -= f * b[col]
    return b


def least_squares(rows):
    normal = [[0.0] * 3 for _ in range(3)]
    rhs = [0.0] * 3
    for x, y, w in rows:
        for i in range(3):
            rhs[i] += w * x[i] * y
            for j in range(3):
                normal[i][j] += w * x[i] * x[j]
    return solve_linear(normal, rhs)


def main():
    ap = argparse.ArgumentParser(
        description="Fit XML clock delay/rate/acel correction from frinZ --search rows."
    )
    ap.add_argument("files", nargs="*", help="frinZ text output files; stdin if omitted")
    ap.add_argument("--fs", type=float, default=1024.0e6, help="sample rate [Hz]")
    ap.add_argument("--obsfreq", type=float, default=6.6e9, help="fringe carrier frequency [Hz]")
    ap.add_argument(
        "--t0",
        type=float,
        default=None,
        help="fit epoch in seconds of day; default is first row epoch",
    )
    ap.add_argument("--current-delay", type=float, default=0.0, help="current XML delay [s]")
    ap.add_argument("--current-rate", type=float, default=0.0, help="current XML rate [s/s]")
    ap.add_argument("--current-acel", type=float, default=0.0, help="current XML acel [s/s^2]")
    ap.add_argument(
        "--rate-weight",
        type=float,
        default=1.0,
        help="relative weight for Res-Rate equations",
    )
    ap.add_argument(
        "--delay-weight",
        type=float,
        default=1.0,
        help="relative weight for Res-Delay equations",
    )
    args = ap.parse_args()

    rows = parse_rows(args.files)
    if not rows:
        print("no frinZ rows parsed", file=sys.stderr)
        return 1

    t0 = rows[0]["tabs"] if args.t0 is None else args.t0
    eqs_delay = []
    eqs_both = []
    for r in rows:
        t = r["tabs"] - t0 + 0.5 * r["length"]
        delay_s = r["delay_sample"] / args.fs
        rate_sps = r["rate_hz"] / args.obsfreq
        xd = [1.0, t, 0.5 * t * t]
        xr = [0.0, 1.0, t]
        eqs_delay.append((xd, delay_s, args.delay_weight))
        eqs_both.append((xd, delay_s, args.delay_weight))
        eqs_both.append((xr, rate_sps, args.rate_weight))

    d_only = least_squares(eqs_delay)
    d_both = least_squares(eqs_both)

    def emit(name, delta):
        delay = args.current_delay + delta[0]
        rate = args.current_rate + delta[1]
        acel = args.current_acel + delta[2]
        print(f"# {name}")
        print(f"delta_delay_s {delta[0]:+.12e}")
        print(f"delta_rate_sps {delta[1]:+.12e}")
        print(f"delta_acel_sps2 {delta[2]:+.12e}")
        print(f"xml_delay_s {delay:.12e}")
        print(f"xml_rate_sps {rate:.12e}")
        print(f"xml_acel_sps2 {acel:.12e}")
        print(
            f"<delay>{delay:.12e}</delay><rate>{rate:.12e}</rate><acel>{acel:.12e}</acel>"
        )

    print(f"# rows {len(rows)}")
    print(f"# t0_epoch {rows[0]['epoch'] if args.t0 is None else args.t0}")
    print(f"# fs {args.fs:.9e} obsfreq {args.obsfreq:.9e}")
    emit("fit_from_delay_only", d_only)
    emit("fit_from_delay_and_rate", d_both)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
