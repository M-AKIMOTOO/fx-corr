#!/usr/bin/env python3
"""Generate independent Niell/Saastamoinen vectors with the author's RTKLIB C.

Reference code is cached under /tmp, never copied into the Rust implementation.
Only the input/output vectors and source hashes are saved in reference/.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import urllib.request

PROJECT = Path(__file__).resolve().parents[1]
CACHE = Path("/tmp/fx-corr-atmosphere-rtklib")
DEST = PROJECT / "reference/atmosphere/validation"


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    sources = []
    for name in ("rtkcmn.c", "rtklib.h"):
        url = "https://raw.githubusercontent.com/tomojitakasu/RTKLIB/master/src/" + name
        path = CACHE / name
        if not path.exists():
            with urllib.request.urlopen(url, timeout=30) as response:
                path.write_bytes(response.read())
        sources.append(dict(url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    probe = CACHE / "probe.c"
    probe.write_text(r'''
#include "rtkcmn.c"
int main(void) {
    const double epochs[][6]={{2024,1,28,0,0,0},{2024,7,28,12,0,0},{2024,12,31,23,59,59}};
    const double lats[]={15.0,37.0,-37.0,75.0};
    const double heights[]={0.0,1200.0};
    const double elevations[]={5.0,30.0,90.0};
    for(int t=0;t<3;t++) for(int l=0;l<4;l++) for(int h=0;h<2;h++) for(int e=0;e<3;e++) {
        double pos[]={lats[l]*D2R,0,heights[h]}, azel[]={0,elevations[e]*D2R}, wet;
        gtime_t time=epoch2time(epochs[t]);
        double dry=nmf(time,pos,azel,&wet);
        printf("%d,%.1f,%.1f,%.1f,%.17g,%.17g\n",t,lats[l],heights[h],elevations[e],dry,wet);
    }
    return 0;
}
''')
    binary = CACHE / "probe"
    subprocess.run(["gcc", "-w", "-O2", "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections", str(probe), "-o", str(binary), "-lm", "-lpthread"], check=True)
    output = subprocess.check_output([str(binary)], text=True)
    rows = []
    for line in output.splitlines():
        t, lat, height, elevation, dry, wet = line.split(",")
        rows.append(dict(epoch_index=int(t), latitude_deg=float(lat), height_m=float(height), elevation_deg=float(elevation), dry_mapping=float(dry), wet_mapping=float(wet)))
    DEST.mkdir(parents=True, exist_ok=True)
    payload = dict(reference="T. Takasu RTKLIB, nmf() original C implementation", sources=sources, epochs_utc=["2024/028 00:00:00", "2024/210 12:00:00", "2024/366 23:59:59"], rows=rows)
    (DEST / "rtklib-niell.json").write_text(json.dumps(payload, indent=2) + "\n")
    (DEST / "rtklib-niell.csv").write_text("epoch_index,latitude_deg,height_m,elevation_deg,dry_mapping,wet_mapping\n" + output)
    print(f"saved {len(rows)} independent C reference vectors")


if __name__ == "__main__":
    main()
