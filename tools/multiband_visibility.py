#!/usr/bin/env python3
"""Stream yi-corr joint.mbcor v1 into frequency-aware spectra or complex means.

No third-party dependencies. Band weights are inverse per-channel noise
variances after applying the band scales. Exact zero bins are excluded from
the mean because yi-corr uses them for channels outside the common band.
"""
import argparse
import contextlib
import csv
import json
import math
import struct


def read_exact(stream, count):
    data = stream.read(count)
    if len(data) != count:
        raise ValueError("truncated multiband file")
    return data


def export(args):
    with contextlib.ExitStack() as stack:
        stream = stack.enter_context(open(args.input, "rb"))
        magic, version, bands, bandwidth, reference = struct.unpack("<8sIIdd", read_exact(stream, 32))
        if magic != b"YIMBCOR\0" or version != 1 or bands != 2:
            raise ValueError("expected yi-corr .mbcor version 1, two bands")
        layouts = []
        for _ in range(bands):
            header = read_exact(stream, 256)
            fs, low, fft, rows = struct.unpack_from("<idii", header, 12)
            if header[:4] != b"\x83\xf9\xa2\x3e" or fs <= 0 or fft < 2 or fft % 2 or rows < 1 or not math.isfinite(low):
                raise ValueError("invalid embedded .cor header")
            layouts.append((low, fs / fft, fft // 2, rows))
        if layouts[0][1:] != layouts[1][1:] or not 0 < layouts[0][0] + layouts[0][1] * layouts[0][2] <= layouts[1][0]:
            raise ValueError("inconsistent multiband frequency grids")
        if not math.isfinite(reference) or reference <= 0 or bandwidth != sum(df * n for _, df, n, _ in layouts):
            raise ValueError("inconsistent joint frequency metadata")
        averages = spectra = None
        if args.average:
            averages = csv.writer(stack.enter_context(open(args.average, "w", newline="")), delimiter="\t")
            averages.writerow(["unix_s", "integration_s", "reference_hz", "real", "imag", "amplitude",
                               "band1_real", "band1_imag", "band2_real", "band2_imag", "usable_bandwidth_hz"])
        if args.spectra:
            spectra = csv.writer(stack.enter_context(open(args.spectra, "w", newline="")), delimiter="\t")
            spectra.writerow(["unix_s", "integration_s", "band", "frequency_hz", "real", "imag", "weight"])
        time_sum, time_weight, total_duration = 0j, 0.0, 0.0
        band_time_sum, band_time_weight = [0j, 0j], [0.0, 0.0]
        first_unix, bandwidth_duration = None, 0.0
        for _ in range(layouts[0][3]):
            means, counts, boundaries = [], [], []
            usable_bandwidth = 0.0
            for band, (low, df, channels, _) in enumerate(layouts):
                sector = read_exact(stream, 128)
                sec, nsec = struct.unpack_from("<iI", sector)
                duration, = struct.unpack_from("<f", sector, 112)
                if nsec >= 1_000_000_000 or not math.isfinite(duration) or duration <= 0:
                    raise ValueError("invalid sector time/integration")
                boundaries.append((sector[:16], sector[112:116]))
                unix_s = sec + nsec * 1e-9
                total, count = 0j, 0
                for k, (real, imag) in enumerate(struct.iter_unpack("<ff", read_exact(stream, channels * 8))):
                    if not math.isfinite(real) or not math.isfinite(imag):
                        raise ValueError("nonfinite visibility")
                    valid = real != 0 or imag != 0
                    value = complex(real, imag) * args.band_scales[band]
                    if valid:
                        total += value
                        count += 1
                    if spectra:
                        spectra.writerow([f"{unix_s:.9f}", duration, band + 1, low + k * df,
                                          value.real, value.imag, args.band_weights[band] if valid else 0])
                means.append(total / count if count else 0j)
                counts.append(count)
                usable_bandwidth += count * df
            if boundaries[0] != boundaries[1]:
                raise ValueError("unequal multiband time boundaries")
            weights = [n * w for n, w in zip(counts, args.band_weights)]
            denominator = sum(weights)
            if not denominator:
                raise ValueError("no usable weighted channels")
            joint = sum(z * w for z, w in zip(means, weights)) / denominator
            if averages and not args.time_average:
                averages.writerow([f"{unix_s:.9f}", duration, reference, joint.real, joint.imag, abs(joint),
                                   means[0].real, means[0].imag, means[1].real, means[1].imag, usable_bandwidth])
            if args.time_average:
                if first_unix is None:
                    first_unix = unix_s
                time_sum += joint * denominator * duration
                time_weight += denominator * duration
                total_duration += duration
                bandwidth_duration += usable_bandwidth * duration
                for band in range(bands):
                    band_time_sum[band] += means[band] * weights[band] * duration
                    band_time_weight[band] += weights[band] * duration
        if args.time_average:
            joint = time_sum / time_weight
            band_mean = [z / w if w else 0j for z, w in zip(band_time_sum, band_time_weight)]
            averages.writerow([f"{first_unix:.9f}", total_duration, reference, joint.real, joint.imag, abs(joint),
                               band_mean[0].real, band_mean[0].imag, band_mean[1].real, band_mean[1].imag,
                               bandwidth_duration / total_duration])
        if stream.read(1):
            raise ValueError("unexpected trailing multiband data")
        print(json.dumps({"format": "yi-multiband-v1", "reference_hz": reference, "occupied_bandwidth_hz": bandwidth,
                          "bands": [{"low_hz": low, "channel_spacing_hz": df, "channels": n, "rows": rows}
                                    for low, df, n, rows in layouts]}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input")
    parser.add_argument("--average", help="write complex weighted band/channel means as TSV")
    parser.add_argument("--spectra", help="write every channel with its actual RF frequency as TSV")
    parser.add_argument("--time-average", action="store_true", help="write one coherent mean over all integrations; requires --average")
    parser.add_argument("--band-weights", type=float, nargs=2, default=[1.0, 1.0], help="inverse noise variance per scaled channel")
    parser.add_argument("--band-scales", type=float, nargs=2, default=[1.0, 1.0], help="amplitude scale for each band before averaging")
    args = parser.parse_args()
    if args.time_average and not args.average:
        parser.error("--time-average requires --average")
    if any(not math.isfinite(v) or v <= 0 for v in args.band_weights + args.band_scales):
        parser.error("band weights and scales must be finite and positive")
    try:
        export(args)
    except (OSError, ValueError, struct.error) as exc:
        parser.exit(1, f"multiband_visibility: {exc}\n")


if __name__ == "__main__":
    main()
