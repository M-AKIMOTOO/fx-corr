#!/usr/bin/env python3

import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "phased_visibility_validation.py"


def write_cor(path, left, right, xyz_left, xyz_right, visibility):
    sectors, nfreq = visibility.shape
    fft = 2 * nfreq
    header = bytearray(256)
    header[0:4] = bytes([0x83, 0xF9, 0xA2, 0x3E])
    struct.pack_into("<i", header, 12, 1024)
    struct.pack_into("<d", header, 16, 6.6e9)
    struct.pack_into("<i", header, 24, fft)
    struct.pack_into("<i", header, 28, sectors)
    header[32 : 32 + len(left)] = left.encode()
    struct.pack_into("<3d", header, 48, *xyz_left)
    header[80 : 80 + len(right)] = right.encode()
    struct.pack_into("<3d", header, 96, *xyz_right)
    header[128:134] = b"TARGET"

    with path.open("wb") as stream:
        stream.write(header)
        for index, spectrum in enumerate(visibility):
            sector = bytearray(128)
            struct.pack_into("<i", sector, 0, 1_700_000_000 + index)
            struct.pack_into("<I", sector, 4, 0)
            struct.pack_into("<f", sector, 112, 1.0)
            stream.write(sector)
            stream.write(np.asarray(spectrum, dtype="<c8").tobytes())


class PhasedVisibilityValidationTest(unittest.TestCase):
    def test_exact_three_station_and_phased_solution(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shape = (3, 8)
            amplitude = np.arange(1, 1 + np.prod(shape), dtype=np.float64).reshape(
                shape
            )
            g32 = np.exp(0.2j)
            g34 = np.exp(-0.35j)
            gh = np.exp(0.05j)
            v32h = amplitude * g32 * np.conj(gh)
            v34h = amplitude * g34 * np.conj(gh)
            v32v34 = amplitude * g32 * np.conj(g34)
            weight = 2.0**-0.5
            predicted = weight * v32h + weight * v34h
            phased_scale = 1.3 * np.exp(0.17j)
            v66h = phased_scale * predicted
            ones = np.ones(shape, dtype=np.complex128)

            xyz32 = (0.0, 0.0, 0.0)
            xyz34 = (110.0, 0.0, 0.0)
            xyzh = (873_000.0, 0.0, 0.0)
            paths = {
                name: root / f"{name}.cor"
                for name in (
                    "y32_y34",
                    "y32_h",
                    "y34_h",
                    "y66_h",
                    "y32_acf",
                    "y34_acf",
                    "h_acf",
                    "y66_acf",
                )
            }
            write_cor(paths["y32_y34"], "YAMAGU32", "YAMAGU34", xyz32, xyz34, v32v34)
            write_cor(paths["y32_h"], "YAMAGU32", "HITACH32", xyz32, xyzh, v32h)
            # Exercise automatic reversal/conjugation.
            write_cor(
                paths["y34_h"],
                "HITACH32",
                "YAMAGU34",
                xyzh,
                xyz34,
                np.conj(v34h),
            )
            write_cor(paths["y66_h"], "YAMAGU66", "HITACH32", xyz32, xyzh, v66h)
            write_cor(paths["y32_acf"], "YAMAGU32", "YAMAGU32", xyz32, xyz32, ones)
            write_cor(paths["y34_acf"], "YAMAGU34", "YAMAGU34", xyz34, xyz34, ones)
            write_cor(paths["h_acf"], "HITACH32", "HITACH32", xyzh, xyzh, ones)
            write_cor(paths["y66_acf"], "YAMAGU66", "YAMAGU66", xyz32, xyz32, ones)

            output = root / "validation.npz"
            command = [sys.executable, str(TOOL)]
            for option, key in (
                ("--y32-y34", "y32_y34"),
                ("--y32-h", "y32_h"),
                ("--y34-h", "y34_h"),
                ("--y66-h", "y66_h"),
                ("--y32-acf", "y32_acf"),
                ("--y34-acf", "y34_acf"),
                ("--h-acf", "h_acf"),
                ("--y66-acf", "y66_acf"),
            ):
                command.extend((option, str(paths[key])))
            command.extend(("--output", str(output)))
            subprocess.run(command, check=True, capture_output=True, text=True)

            with np.load(output) as archive:
                self.assertEqual(int(archive["schema_version"]), 1)
                self.assertAlmostEqual(float(archive["baseline_y32_y34_m"]), 110.0)
                self.assertAlmostEqual(
                    float(archive["long_baseline_amplitude_pearson_all"]), 1.0, places=7
                )
                self.assertAlmostEqual(
                    float(archive["long_baseline_complex_coherence"]), 1.0, places=7
                )
                self.assertAlmostEqual(
                    float(archive["prediction_complex_coherence"]), 1.0, places=7
                )
                self.assertLess(
                    float(archive["prediction_normalized_residual_raw"]), 1.0e-7
                )
                np.testing.assert_allclose(
                    archive["closure_phase_y32_y34_h_rad"], 0.0, atol=1.0e-6
                )


if __name__ == "__main__":
    unittest.main()
