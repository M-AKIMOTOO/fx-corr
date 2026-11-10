#!/usr/bin/env python3
"""Check a ZIP .info with Python, the user's zipcat, and zedit's actual ZIP engine.

The zedit sources are read only; the probe and edited round-trip copies are
created under --output. No observation or archived reference is edited.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("info", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--zipcat", type=Path, default=Path("/home/akimoto/program/c/zipcat/zipcat"))
    parser.add_argument("--zedit", type=Path, default=Path("/home/akimoto/program/c/zedit"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.info) as archive:
        assert archive.namelist() == ["info.txt"]
        member = archive.getinfo("info.txt")
        assert member.compress_type == zipfile.ZIP_DEFLATED
        assert not member.flag_bits & 1
        assert archive.testzip() is None
        expected = archive.read("info.txt")
        expected.decode("utf-8")
        compressed_bytes = member.compress_size
    zipcat = subprocess.run([str(args.zipcat), str(args.info)], capture_output=True, check=True)
    assert zipcat.stdout == expected, "zipcat text differs"
    probe = args.output.resolve() / "zedit-probe"
    subprocess.run(["gcc", "-O2", "-std=c11", "-D_POSIX_C_SOURCE=200809L", "-Wall", "-Wextra",
                    "-I", str(args.zedit / "src"), str(args.zedit / "tests/probe.c"),
                    str(args.zedit / "src/zip.c"), "-o", str(probe)], check=True)
    zedit = subprocess.run([str(probe), "dump", str(args.info), "0"], capture_output=True, check=True)
    assert zedit.stdout == expected, "zedit text differs"
    replacement = args.output / "edited.txt"
    edited = args.output / "edited.info"
    replacement.write_bytes(expected + "\n# zedit UTF-8 確認\n".encode("utf-8"))
    subprocess.run([str(probe), "replace", str(args.info), "0", str(replacement), str(edited)], check=True)
    with zipfile.ZipFile(edited) as archive:
        assert archive.testzip() is None
        assert archive.read("info.txt") == replacement.read_bytes()
    result = dict(archive=str(args.info.resolve()), sha256=hashlib.sha256(args.info.read_bytes()).hexdigest(),
                  member="info.txt", uncompressed_bytes=len(expected), compressed_bytes=compressed_bytes,
                  zipcat_read_matches=True, zedit_engine_read_matches=True, zedit_edit_roundtrip=True)
    (args.output / "compatibility.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
