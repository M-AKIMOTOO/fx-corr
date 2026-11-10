#!/usr/bin/env python3
"""Copy a legacy gzip .info log into a ZIP .info without changing the source."""
import argparse
import codecs
import gzip
from pathlib import Path
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, help="default: SOURCE-STEM.zip.info; must not exist")
    args = parser.parse_args()
    output = args.output or args.source.with_name(args.source.stem + ".zip.info")
    with args.source.open("rb") as stream:
        if stream.read(2) != b"\x1f\x8b":
            parser.error("source must be a legacy gzip .info")
    decoder = codecs.getincrementaldecoder("utf-8")()
    with output.open("xb") as destination:
        with gzip.open(args.source, "rb") as source, zipfile.ZipFile(
                destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=False) as archive:
            with archive.open("info.txt", "w") as member:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    decoder.decode(block)
                    member.write(block)
                decoder.decode(b"", final=True)
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise ValueError("ZIP checksum verification failed")
    print(output)


if __name__ == "__main__":
    main()
