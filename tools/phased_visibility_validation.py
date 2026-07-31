#!/usr/bin/env python3
"""Build a complete yi-phasedarray visibility-validation NPZ from .cor files."""

from __future__ import annotations

import argparse
import json
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

FILE_HEADER_BYTES = 256
SECTOR_HEADER_BYTES = 128


@dataclass(frozen=True)
class CorData:
    path: Path
    sampling_hz: float
    observing_hz: float
    fft: int
    station1: str
    station2: str
    ecef1_m: np.ndarray
    ecef2_m: np.ndarray
    source: str
    time_unix_s: np.ndarray
    integration_s: np.ndarray
    visibility: np.ndarray

    @property
    def frequency_if_mhz(self) -> np.ndarray:
        return np.arange(self.fft // 2, dtype=np.float64) * (
            self.sampling_hz / self.fft / 1.0e6
        )

    @property
    def baseline_length_m(self) -> float:
        return float(np.linalg.norm(self.ecef2_m - self.ecef1_m))


def _ascii_field(buf: bytes, start: int, length: int) -> str:
    return buf[start : start + length].split(b"\0", 1)[0].decode("ascii").strip()


def read_cor(path: Path) -> CorData:
    raw = path.read_bytes()
    if len(raw) < FILE_HEADER_BYTES:
        raise ValueError(f"{path}: shorter than the 256-byte .cor header")
    sampling_hz = float(struct.unpack_from("<i", raw, 12)[0])
    observing_hz = struct.unpack_from("<d", raw, 16)[0]
    fft = struct.unpack_from("<i", raw, 24)[0]
    sectors = struct.unpack_from("<i", raw, 28)[0]
    if sampling_hz <= 0 or fft <= 0 or fft % 2 or sectors < 0:
        raise ValueError(
            f"{path}: invalid header fs={sampling_hz} fft={fft} sectors={sectors}"
        )
    nfreq = fft // 2
    sector_bytes = SECTOR_HEADER_BYTES + nfreq * np.dtype("<c8").itemsize
    expected = FILE_HEADER_BYTES + sectors * sector_bytes
    if len(raw) != expected:
        raise ValueError(f"{path}: size={len(raw)} expected={expected}")

    times = np.empty(sectors, dtype=np.float64)
    integrations = np.empty(sectors, dtype=np.float64)
    vis = np.empty((sectors, nfreq), dtype=np.complex64)
    offset = FILE_HEADER_BYTES
    for index in range(sectors):
        sec = struct.unpack_from("<i", raw, offset)[0]
        nsec = struct.unpack_from("<I", raw, offset + 4)[0]
        times[index] = sec + nsec * 1.0e-9
        integrations[index] = struct.unpack_from("<f", raw, offset + 112)[0]
        data_offset = offset + SECTOR_HEADER_BYTES
        vis[index] = np.frombuffer(raw, dtype="<c8", count=nfreq, offset=data_offset)
        offset += sector_bytes

    return CorData(
        path=path,
        sampling_hz=sampling_hz,
        observing_hz=observing_hz,
        fft=fft,
        station1=_ascii_field(raw, 32, 16),
        station2=_ascii_field(raw, 80, 16),
        ecef1_m=np.array(struct.unpack_from("<3d", raw, 48), dtype=np.float64),
        ecef2_m=np.array(struct.unpack_from("<3d", raw, 96), dtype=np.float64),
        source=_ascii_field(raw, 128, 16),
        time_unix_s=times,
        integration_s=integrations,
        visibility=vis,
    )


def orient_visibility(data: CorData, left: str, right: str) -> np.ndarray:
    if data.station1 == left and data.station2 == right:
        return data.visibility.astype(np.complex128)
    if data.station1 == right and data.station2 == left:
        return data.visibility.conj().astype(np.complex128)
    raise ValueError(
        f"{data.path}: header baseline {data.station1}-{data.station2}, "
        f"expected {left}-{right} (either order)"
    )


def orient_acf(data: CorData, station: str) -> np.ndarray:
    if data.station1 != station or data.station2 != station:
        raise ValueError(
            f"{data.path}: header baseline {data.station1}-{data.station2}, "
            f"expected {station}-{station}"
        )
    return data.visibility.real.astype(np.float64)


def station_ecef(data: CorData, station: str) -> np.ndarray:
    if data.station1 == station:
        return data.ecef1_m
    if data.station2 == station:
        return data.ecef2_m
    raise ValueError(f"{data.path}: station {station} is not present in the header")


def validate_common_grid(reference: CorData, datasets: list[CorData]) -> None:
    for data in datasets:
        if data.fft != reference.fft:
            raise ValueError(f"{data.path}: FFT differs from {reference.path}")
        if not np.isclose(data.sampling_hz, reference.sampling_hz, rtol=0, atol=0.5):
            raise ValueError(f"{data.path}: sampling frequency differs from reference")
        if not np.isclose(data.observing_hz, reference.observing_hz, rtol=0, atol=0.5):
            raise ValueError(f"{data.path}: observing frequency differs from reference")
        if data.visibility.shape != reference.visibility.shape:
            raise ValueError(
                f"{data.path}: time/frequency shape differs from reference"
            )
        if not np.allclose(
            data.time_unix_s, reference.time_unix_s, rtol=0, atol=1.0e-6
        ):
            raise ValueError(f"{data.path}: sector timestamps differ from reference")
        if not np.allclose(
            data.integration_s, reference.integration_s, rtol=0, atol=1.0e-6
        ):
            raise ValueError(f"{data.path}: integration times differ from reference")


def normalized_visibility(
    cross: np.ndarray, auto_left: np.ndarray, auto_right: np.ndarray
) -> np.ndarray:
    power = auto_left * auto_right
    valid = np.isfinite(power) & (power > 0.0)
    result = np.full(cross.shape, np.nan + 1j * np.nan, dtype=np.complex128)
    result[valid] = cross[valid] / np.sqrt(power[valid])
    return result


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(valid) < 2:
        return float("nan")
    xv = x[valid].astype(np.float64)
    yv = y[valid].astype(np.float64)
    xv -= xv.mean()
    yv -= yv.mean()
    denom = np.linalg.norm(xv) * np.linalg.norm(yv)
    return float(np.dot(xv, yv) / denom) if denom > 0.0 else float("nan")


def pearson_along(x: np.ndarray, y: np.ndarray, axis: int) -> np.ndarray:
    out_shape = x.shape[1 - axis]
    out = np.empty(out_shape, dtype=np.float64)
    for index in range(out_shape):
        if axis == 1:
            out[index] = pearson(x[index, :], y[index, :])
        else:
            out[index] = pearson(x[:, index], y[:, index])
    return out


def complex_coherence(x: np.ndarray, y: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(x.real) & np.isfinite(y.real)
    if not np.any(valid):
        return float("nan")
    xv = x[valid]
    yv = y[valid]
    denom = np.linalg.norm(xv) * np.linalg.norm(yv)
    return float(abs(np.vdot(yv, xv)) / denom) if denom > 0.0 else float("nan")


def fit_complex_scale(
    predicted: np.ndarray, measured: np.ndarray, mask: np.ndarray
) -> tuple[complex, float]:
    valid = mask & np.isfinite(predicted.real) & np.isfinite(measured.real)
    if not np.any(valid):
        return complex(np.nan, np.nan), float("nan")
    pv = predicted[valid]
    mv = measured[valid]
    denom = np.vdot(pv, pv).real
    if denom <= 0.0:
        return complex(np.nan, np.nan), float("nan")
    scale = np.vdot(pv, mv) / denom
    measured_norm = np.linalg.norm(mv)
    residual = np.linalg.norm(mv - scale * pv)
    return complex(scale), float(
        residual / measured_norm
    ) if measured_norm > 0 else float("nan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine the complete three-station and phased-station .cor set into "
            "one compressed NumPy validation archive."
        )
    )
    parser.add_argument("--y32-y34", type=Path, required=True)
    parser.add_argument("--y32-h", type=Path, required=True)
    parser.add_argument("--y34-h", type=Path, required=True)
    parser.add_argument("--y66-h", type=Path, required=True)
    parser.add_argument("--y32-acf", type=Path, required=True)
    parser.add_argument("--y34-acf", type=Path, required=True)
    parser.add_argument("--h-acf", type=Path, required=True)
    parser.add_argument("--y66-acf", type=Path, required=True)
    parser.add_argument(
        "--y66-prequant-acf",
        type=Path,
        help="optional yi-phasedarray ACF before raw-output requantization",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--y32-name", default="YAMAGU32")
    parser.add_argument("--y34-name", default="YAMAGU34")
    parser.add_argument("--h-name", default="HITACH32")
    parser.add_argument("--y66-name", default="YAMAGU66")
    parser.add_argument("--weight-y32", type=float, default=2.0**-0.5)
    parser.add_argument("--weight-y34", type=float, default=2.0**-0.5)
    parser.add_argument("--freq-min", type=float, default=None, help="minimum IF MHz")
    parser.add_argument("--freq-max", type=float, default=None, help="maximum IF MHz")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    files = {
        "y32_y34": read_cor(args.y32_y34),
        "y32_h": read_cor(args.y32_h),
        "y34_h": read_cor(args.y34_h),
        "y66_h": read_cor(args.y66_h),
        "y32_acf": read_cor(args.y32_acf),
        "y34_acf": read_cor(args.y34_acf),
        "h_acf": read_cor(args.h_acf),
        "y66_acf": read_cor(args.y66_acf),
    }
    if args.y66_prequant_acf is not None:
        files["y66_prequant_acf"] = read_cor(args.y66_prequant_acf)
    reference = files["y32_h"]
    validate_common_grid(reference, list(files.values()))

    y32_y34 = orient_visibility(files["y32_y34"], args.y32_name, args.y34_name)
    y32_h = orient_visibility(files["y32_h"], args.y32_name, args.h_name)
    y34_h = orient_visibility(files["y34_h"], args.y34_name, args.h_name)
    y66_h = orient_visibility(files["y66_h"], args.y66_name, args.h_name)
    y32_acf = orient_acf(files["y32_acf"], args.y32_name)
    y34_acf = orient_acf(files["y34_acf"], args.y34_name)
    h_acf = orient_acf(files["h_acf"], args.h_name)
    y66_acf = orient_acf(files["y66_acf"], args.y66_name)
    y66_prequant_acf = (
        orient_acf(files["y66_prequant_acf"], args.y66_name)
        if "y66_prequant_acf" in files
        else None
    )

    coefficient_y32_y34 = normalized_visibility(y32_y34, y32_acf, y34_acf)
    coefficient_y32_h = normalized_visibility(y32_h, y32_acf, h_acf)
    coefficient_y34_h = normalized_visibility(y34_h, y34_acf, h_acf)
    coefficient_y66_h = normalized_visibility(y66_h, y66_acf, h_acf)

    frequency_if_mhz = reference.frequency_if_mhz
    frequency_mask = np.ones(frequency_if_mhz.shape, dtype=bool)
    if args.freq_min is not None:
        frequency_mask &= frequency_if_mhz >= args.freq_min
    if args.freq_max is not None:
        frequency_mask &= frequency_if_mhz <= args.freq_max
    analysis_mask = np.broadcast_to(frequency_mask, y32_h.shape).copy()

    amp_y32_h = np.abs(coefficient_y32_h)
    amp_y34_h = np.abs(coefficient_y34_h)
    amp_y66_h = np.abs(coefficient_y66_h)
    predicted_raw = args.weight_y32 * y32_h + args.weight_y34 * y34_h
    predicted_coefficient = (
        args.weight_y32 * coefficient_y32_h + args.weight_y34 * coefficient_y34_h
    )
    fitted_scale_raw, residual_raw = fit_complex_scale(
        predicted_raw, y66_h, analysis_mask
    )
    fitted_scale_coefficient, residual_coefficient = fit_complex_scale(
        predicted_coefficient, coefficient_y66_h, analysis_mask
    )

    closure = y32_y34 * y34_h * np.conj(y32_h)
    closure_phase_rad = np.angle(closure)
    long_amp_ratio = np.divide(
        amp_y32_h,
        amp_y34_h,
        out=np.full_like(amp_y32_h, np.nan),
        where=np.isfinite(amp_y34_h) & (amp_y34_h > 0.0),
    )
    amp_pearson_all = pearson(amp_y32_h[analysis_mask], amp_y34_h[analysis_mask])
    amp_pearson_time = pearson_along(
        np.where(analysis_mask, amp_y32_h, np.nan),
        np.where(analysis_mask, amp_y34_h, np.nan),
        axis=1,
    )
    amp_pearson_frequency = pearson_along(
        np.where(analysis_mask, amp_y32_h, np.nan),
        np.where(analysis_mask, amp_y34_h, np.nan),
        axis=0,
    )
    long_complex_coherence = complex_coherence(
        coefficient_y32_h, coefficient_y34_h, analysis_mask
    )
    phased_prediction_coherence = complex_coherence(
        coefficient_y66_h, predicted_coefficient, analysis_mask
    )
    selected_ratio = long_amp_ratio[analysis_mask]
    long_amp_ratio_median = (
        float(np.nanmedian(selected_ratio))
        if np.any(np.isfinite(selected_ratio))
        else float("nan")
    )
    selected_closure_phase = closure_phase_rad[analysis_mask]
    selected_closure_phase = selected_closure_phase[np.isfinite(selected_closure_phase)]
    closure_resultant = (
        np.mean(np.exp(1j * selected_closure_phase))
        if selected_closure_phase.size
        else complex(np.nan, np.nan)
    )

    metadata = {
        "schema": "yi-phasedarray-visibility-validation-v1",
        "source": reference.source,
        "stations": {
            "y32": args.y32_name,
            "y34": args.y34_name,
            "h": args.h_name,
            "y66": args.y66_name,
        },
        "input_files": {key: str(data.path) for key, data in files.items()},
        "weights": {"y32": args.weight_y32, "y34": args.weight_y34},
        "analysis_if_mhz": [args.freq_min, args.freq_max],
        "visibility_convention": "V_left_right = E_left * conj(E_right)",
    }
    payload = {
        "schema_version": np.array(1, dtype=np.int32),
        "metadata_json": np.array(json.dumps(metadata, sort_keys=True)),
        "time_unix_s": reference.time_unix_s,
        "integration_s": reference.integration_s,
        "sampling_hz": np.array(reference.sampling_hz),
        "observing_frequency_hz": np.array(reference.observing_hz),
        "fft": np.array(reference.fft, dtype=np.int32),
        "frequency_if_mhz": frequency_if_mhz,
        "analysis_frequency_mask": frequency_mask,
        "ecef_y32_m": station_ecef(files["y32_h"], args.y32_name),
        "ecef_y34_m": station_ecef(files["y34_h"], args.y34_name),
        "ecef_h_m": station_ecef(files["y32_h"], args.h_name),
        "baseline_y32_y34_m": np.array(files["y32_y34"].baseline_length_m),
        "baseline_y32_h_m": np.array(files["y32_h"].baseline_length_m),
        "baseline_y34_h_m": np.array(files["y34_h"].baseline_length_m),
        "baseline_y66_h_m": np.array(files["y66_h"].baseline_length_m),
        "long_baseline_length_difference_m": np.array(
            files["y32_h"].baseline_length_m - files["y34_h"].baseline_length_m
        ),
        "acf_y32": y32_acf,
        "acf_y34": y34_acf,
        "acf_h": h_acf,
        "acf_y66": y66_acf,
        "visibility_y32_y34": y32_y34,
        "visibility_y32_h": y32_h,
        "visibility_y34_h": y34_h,
        "visibility_y66_h": y66_h,
        "visibility_amplitude_y32_y34": np.abs(y32_y34),
        "visibility_amplitude_y32_h": np.abs(y32_h),
        "visibility_amplitude_y34_h": np.abs(y34_h),
        "visibility_amplitude_y66_h": np.abs(y66_h),
        "visibility_phase_y32_y34_rad": np.angle(y32_y34),
        "visibility_phase_y32_h_rad": np.angle(y32_h),
        "visibility_phase_y34_h_rad": np.angle(y34_h),
        "visibility_phase_y66_h_rad": np.angle(y66_h),
        "visibility_coefficient_y32_y34": coefficient_y32_y34,
        "visibility_coefficient_y32_h": coefficient_y32_h,
        "visibility_coefficient_y34_h": coefficient_y34_h,
        "visibility_coefficient_y66_h": coefficient_y66_h,
        "visibility_coefficient_amplitude_y32_h": amp_y32_h,
        "visibility_coefficient_amplitude_y34_h": amp_y34_h,
        "visibility_coefficient_amplitude_y66_h": amp_y66_h,
        "long_baseline_amplitude_ratio_y32_over_y34": long_amp_ratio,
        "long_baseline_amplitude_pearson_per_time": amp_pearson_time,
        "long_baseline_amplitude_pearson_per_frequency": amp_pearson_frequency,
        "long_baseline_amplitude_pearson_all": np.array(amp_pearson_all),
        "long_baseline_amplitude_ratio_median": np.array(long_amp_ratio_median),
        "long_baseline_complex_coherence": np.array(long_complex_coherence),
        "closure_complex_y32_y34_h": closure,
        "closure_phase_y32_y34_h_rad": closure_phase_rad,
        "closure_phase_circular_mean_rad": np.array(np.angle(closure_resultant)),
        "closure_phase_resultant_length": np.array(abs(closure_resultant)),
        "predicted_visibility_y66_h": predicted_raw,
        "predicted_visibility_coefficient_y66_h": predicted_coefficient,
        "prediction_fit_scale_raw": np.array(fitted_scale_raw),
        "prediction_fit_scale_coefficient": np.array(fitted_scale_coefficient),
        "prediction_normalized_residual_raw": np.array(residual_raw),
        "prediction_normalized_residual_coefficient": np.array(residual_coefficient),
        "prediction_complex_coherence": np.array(phased_prediction_coherence),
    }
    if y66_prequant_acf is not None:
        payload["acf_y66_prequant"] = y66_prequant_acf
        payload["acf_y66_postquant_over_prequant"] = np.divide(
            y66_acf,
            y66_prequant_acf,
            out=np.full_like(y66_acf, np.nan),
            where=np.isfinite(y66_prequant_acf) & (y66_prequant_acf != 0.0),
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **payload)

    print(f"[info] wrote {args.output}")
    print(
        "[result] baselines [m]: "
        f"Y32-Y34={files['y32_y34'].baseline_length_m:.3f} "
        f"Y32-H={files['y32_h'].baseline_length_m:.3f} "
        f"Y34-H={files['y34_h'].baseline_length_m:.3f}"
    )
    print(f"[result] long-baseline amplitude Pearson r={amp_pearson_all:.9f}")
    print(f"[result] long-baseline complex coherence={long_complex_coherence:.9f}")
    print(
        "[result] phased prediction: "
        f"coherence={phased_prediction_coherence:.9f} "
        f"normalized-residual={residual_coefficient:.9e}"
    )


if __name__ == "__main__":
    main()
