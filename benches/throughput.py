#!/usr/bin/env python3
"""Compare release binaries on identical synthetic packed input (stdlib only).

python3 benches/throughput.py --before /tmp/before --after target/release
Each directory must contain yi-corr and yi-phasedarray. Runs alternate order,
discard warm-up, report median wall/compute times, and compare RAW/.cor bytes.
"""
import argparse
import hashlib
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--fft", type=int, nargs="+", default=[128, 4096, 65536])
    parser.add_argument("--cpu", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--samples", type=int, default=2**25)
    parser.add_argument("--rounds", type=int, default=3)
    args = parser.parse_args()
    binaries = [args.before.resolve(), args.after.resolve()]
    if args.samples <= 0 or args.samples % max(args.fft) or args.rounds < 1:
        parser.error("samples must be positive and divisible by every FFT; rounds must be positive")
    for n in args.fft:
        if n < 32 or n & (n - 1) or args.samples % n:
            parser.error("FFT lengths must be powers of two >= 32 and divide samples")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("FX_", "YI_"))}
    with tempfile.TemporaryDirectory(prefix="fx-corr-throughput-") as tmp:
        root = Path(tmp)
        raw = root / "raw"
        raw.mkdir()
        # Repeatable pseudorandom bytes without numpy or a costly sample generator.
        pattern = b"".join(hashlib.sha256(str(i).encode()).digest() for i in range(2048))
        payload = (pattern * ((args.samples // 4 + len(pattern) - 1) // len(pattern)))[:args.samples // 4]
        for station in ["ANT1", "ANT2"]:
            (raw / f"{station}_2000001000000.raw").write_bytes(payload)
        del payload
        print("fft cpu mode before_wall_s after_wall_s speedup before_compute_s after_compute_s products_equal", flush=True)
        for n in args.fft:
            schedule = root / "schedule.xml"
            schedule.write_text(f'''<schedule>
<station key="A"><name>ANT1</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z><terminal>term</terminal></station>
<station key="B"><name>ANT2</name><pos-x>-3502544.587</pos-x><pos-y>3950966.235</pos-y><pos-z>3566381.192</pos-z><terminal>term</terminal></station>
<clock key="A"><delay>0</delay><rate>0</rate></clock><clock key="B"><delay>0</delay><rate>0</rate></clock>
<terminal name="term"><speed>{args.samples}</speed><channel>1</channel><bit>2</bit><level>-1.5,-0.5,0.5,1.5</level></terminal>
<shuffle key="A">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
<shuffle key="B">24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7</shuffle>
<source name="TARGET"><ra>00h00m00.0</ra><dec>+00d00'00.0</dec></source>
<stream><frequency>100000000</frequency><fft>{n}</fft><output>8</output>
<special key="A"><rotation>0</rotation><sideband>LSB</sideband></special>
<special key="B"><rotation>0</rotation><sideband>LSB</sideband></special></stream>
<process><epoch>2000/001 00:00:00</epoch><length>1</length><object>TARGET</object><stations>AB</stations><baseline>AB</baseline></process>
</schedule>''')
            for cpu in args.cpu:
                for mode in ["corr", "phasedarray"]:
                    wall, compute, products = [[], []], [[], []], [None, None]
                    for round_no in range(args.rounds + 1):
                        for variant in [round_no % 2, 1 - round_no % 2]:
                            out = root / f"out-{variant}"
                            out.mkdir(exist_ok=True)
                            for old in out.iterdir():
                                old.unlink()
                            phased = mode == "phasedarray"
                            command = [str(binaries[variant] / ("yi-phasedarray" if phased else "yi-corr")),
                                       "--sc", str(schedule), "--raw", str(raw),
                                       "--output" if phased else "--cor", str(out),
                                       "--cpu", str(cpu), "--no-affinity"]
                            if phased:
                                command.extend(["--phased-name", "ARRAY"])
                            begin = time.perf_counter()
                            result = subprocess.run(command, env=env, capture_output=True, text=True)
                            elapsed = time.perf_counter() - begin
                            if result.returncode:
                                raise RuntimeError(result.stdout + result.stderr)
                            match = re.search(r"Synth timing summary:.*?compute=([0-9.]+)s", result.stdout)
                            if round_no:
                                wall[variant].append(elapsed)
                                if match:
                                    compute[variant].append(float(match[1]))
                            products[variant] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                 for p in out.iterdir() if p.suffix in [".raw", ".cor"]}
                    before, after = map(statistics.median, wall)
                    c = [f"{statistics.median(v):.6f}" if v else "n/a" for v in compute]
                    equal = bool(products[0]) and products[0] == products[1]
                    print(f"{n} {cpu} {mode} {before:.6f} {after:.6f} {before / after:.2f} {c[0]} {c[1]} {equal}", flush=True)
                    if not equal:
                        raise AssertionError(f"products differ: fft={n} cpu={cpu} mode={mode}")


if __name__ == "__main__":
    main()
