#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

MAS_PER_RAD = 180.0 / math.pi * 3600.0 * 1000.0


def parse_delay_model(path):
    header = None
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        if line.startswith('#'):
            if line.startswith('# t_sec'):
                header = line[2:].split('\t')
            continue
        parts = line.split('\t')
        if header and len(parts) == len(header):
            row = {k: float(v) for k, v in zip(header, parts)}
        else:
            raise ValueError(
                f'{path}: delay model is missing the new baseline phase-basis columns; '
                'rerun with the rebuilt target/release/yi-corr and YI_DUMP_DELAY_MODEL=1'
            )
        rows.append(row)
    if not rows:
        raise ValueError(f'no delay-model rows in {path}')
    return rows


def parse_add_plot(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        p = line.split('\t')
        rows.append({
            't_sec': float(p[0]),
            'amp_pct': float(p[1]),
            'snr': float(p[2]),
            'phase_deg': float(p[3]),
            'noise_pct': float(p[4]),
            'res_delay_sample': float(p[5]),
            'res_rate_hz': float(p[6]),
        })
    if not rows:
        raise ValueError(f'no frinZ add-plot rows in {path}')
    return rows


def unwrap_deg(rows):
    prev = None
    offset = 0.0
    out = []
    for r in rows:
        x = r['phase_deg'] + offset
        if prev is not None:
            while x - prev > 180.0:
                offset -= 360.0
                x = r['phase_deg'] + offset
            while x - prev < -180.0:
                offset += 360.0
                x = r['phase_deg'] + offset
        prev = x
        rr = dict(r)
        rr['phase_unwrapped_deg'] = x
        rr['phase_cycles'] = x / 360.0
        out.append(rr)
    return out


def solve_linear(a, b):
    n = len(b)
    a = [row[:] for row in a]
    b = b[:]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1.0e-30:
            raise ValueError('singular fit matrix')
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


def least_squares(design, y, ridge=0.0):
    m = len(design[0])
    normal = [[0.0] * m for _ in range(m)]
    rhs = [0.0] * m
    for row, val in zip(design, y):
        for i in range(m):
            rhs[i] += row[i] * val
            for j in range(m):
                normal[i][j] += row[i] * row[j]
    if ridge != 0.0:
        for i in range(m):
            normal[i][i] += ridge
    coef = solve_linear(normal, rhs)
    residual = [val - sum(c * x for c, x in zip(coef, row)) for row, val in zip(design, y)]
    return coef, residual


def rms_deg(res_cycles):
    return math.sqrt(sum((r * 360.0) ** 2 for r in res_cycles) / len(res_cycles))


def ptp_deg(res_cycles):
    return (max(res_cycles) - min(res_cycles)) * 360.0


def pair_arg(raw):
    if ':' not in raw:
        raise argparse.ArgumentTypeError('pair must be ADD_PLOT_TSV:DELAY_MODEL_TSV')
    a, b = raw.split(':', 1)
    return a, b


def main():
    ap = argparse.ArgumentParser(
        description='Fit residual frinZ phase to RA/Dec or ECEF baseline-vector error using delay_model basis columns.'
    )
    ap.add_argument('pairs', nargs='+', type=pair_arg, help='ADD_PLOT_TSV:DELAY_MODEL_TSV')
    ap.add_argument('--sample-time', choices=['start', 'center'], default='start')
    ap.add_argument('--db-ridge', type=float, default=1.0e-10, help='ridge term for rank-deficient ECEF dB fit')
    args = ap.parse_args()

    all_rows = []
    for scan_idx, (add_path, dm_path) in enumerate(args.pairs):
        dm = parse_delay_model(dm_path)
        by_sec = {round(r['t_sec']): r for r in dm}
        add = unwrap_deg(parse_add_plot(add_path))
        for r in add:
            t = r['t_sec'] + (5.0 if args.sample_time == 'center' else 0.0)
            key = round(t)
            if key not in by_sec:
                raise ValueError(f'{dm_path}: no delay-model row near t={t}')
            b = by_sec[key]
            all_rows.append({**r, **b, 'scan_idx': scan_idx})

    scans = sorted(set(r['scan_idx'] for r in all_rows))
    y = [r['phase_cycles'] for r in all_rows]

    def scan_offsets(r):
        return [1.0 if r['scan_idx'] == s else 0.0 for s in scans]

    design = [scan_offsets(r) + [r['phase_dra_cycles_per_rad'], r['phase_ddec_cycles_per_rad']] for r in all_rows]
    coef, res = least_squares(design, y)
    dra = coef[len(scans)]
    ddec = coef[len(scans) + 1]
    print('# phase = per_scan_offset + dRA*phase_dra + dDec*phase_ddec')
    print(f'dRA_rad {dra:+.12e}')
    print(f'dDec_rad {ddec:+.12e}')
    print(f'dRA_mas {dra * MAS_PER_RAD:+.6f}')
    print(f'dDec_mas {ddec * MAS_PER_RAD:+.6f}')
    print(f'rms_deg {rms_deg(res):.6f}')
    print(f'ptp_deg {ptp_deg(res):.6f}')

    design = []
    y_db = []
    for scan in scans:
        sub = [r for r in all_rows if r['scan_idx'] == scan]
        mean_phase = sum(r['phase_cycles'] for r in sub) / len(sub)
        mean_dbx = sum(r['phase_dbx_cycles_per_m'] for r in sub) / len(sub)
        mean_dby = sum(r['phase_dby_cycles_per_m'] for r in sub) / len(sub)
        mean_dbz = sum(r['phase_dbz_cycles_per_m'] for r in sub) / len(sub)
        for r in sub:
            design.append([
                r['phase_dbx_cycles_per_m'] - mean_dbx,
                r['phase_dby_cycles_per_m'] - mean_dby,
                r['phase_dbz_cycles_per_m'] - mean_dbz,
            ])
            y_db.append(r['phase_cycles'] - mean_phase)
    coef, res = least_squares(design, y_db, args.db_ridge)
    dbx, dby, dbz = coef[:3]
    norm = math.sqrt(dbx * dbx + dby * dby + dbz * dbz)
    print('# scan-mean removed phase = dB_ecef dot scan-mean removed phase_db')
    print(f'dBx_m {dbx:+.12e}')
    print(f'dBy_m {dby:+.12e}')
    print(f'dBz_m {dbz:+.12e}')
    print(f'dB_norm_m {norm:.12e}')
    print(f'dBx_cm {dbx * 100.0:+.6f}')
    print(f'dBy_cm {dby * 100.0:+.6f}')
    print(f'dBz_cm {dbz * 100.0:+.6f}')
    print(f'dB_norm_cm {norm * 100.0:.6f}')
    print(f'rms_deg {rms_deg(res):.6f}')
    print(f'ptp_deg {ptp_deg(res):.6f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
