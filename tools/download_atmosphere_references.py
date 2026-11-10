#!/usr/bin/env python3
"""Download published atmosphere papers; record the actual edition and source.

Publisher endpoints are tried first. Institutional copies of the published
article are acceptable fallbacks; failures stay visible in manifest.json.
Observation raw data are never copied by this tool.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

DEST = Path(__file__).resolve().parents[1] / "reference/atmosphere/literature"
PAPERS = [
    dict(id="degasperin-2018-ionosphere", title="The effect of the ionosphere on ultra-low-frequency radio-interferometric observations", authors="F. de Gasperin et al.", year=2018,
         journal="Astronomy & Astrophysics 615, A179", publisher="EDP Sciences", doi="10.1051/0004-6361/201833012", role="First-order dispersive phase; limitations of neglecting higher orders and Faraday rotation",
         urls=["https://www.aanda.org/articles/aa/pdf/2018/07/aa33012-18.pdf", "https://www.aanda.org/10.1051/0004-6361/201833012/pdf", "https://arxiv.org/pdf/1804.07947"]),
    dict(id="niell-1996-mapping", title="Global mapping functions for the atmosphere delay at radio wavelengths", authors="A. E. Niell", year=1996,
         journal="Journal of Geophysical Research: Solid Earth 101(B2), 3227–3246", publisher="American Geophysical Union / Wiley", doi="10.1029/95JB03048", role="Niell hydrostatic/wet mapping coefficients and season/height dependence; use corrected equations (4)/(5)",
         urls=["https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/95JB03048", "https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/95JB03048", "https://articles.adsabs.harvard.edu/pdf/1996JGR...101.3227N", "https://www.haystack.mit.edu/geo/pubs/NMF_JGR.pdf", "https://safe.nrao.edu/wiki/pub/Main/RefBendDelayCalc/Niell_Global.pdf"]),
    dict(id="saastamoinen-1972-atmosphere", title="Atmospheric Correction for the Troposphere and Stratosphere in Radio Ranging Satellites", authors="J. Saastamoinen", year=1972,
         journal="The Use of Artificial Satellites for Geodesy, Geophysical Monograph Series 15, 247–251", publisher="American Geophysical Union / Wiley", doi="10.1029/GM015p0247", role="Pressure/temperature/water-vapour neutral atmosphere delay",
         urls=["https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/GM015p0247", "https://onlinelibrary.wiley.com/wol1/doi/10.1029/GM015p0247/fullpdf", "https://articles.adsabs.harvard.edu/pdf/1972GMS....15..247S"]),
    dict(id="landskron-boehm-2018-vmf3-gpt3", title="VMF3/GPT3: refined discrete and empirical troposphere mapping functions", authors="D. Landskron and J. Böhm", year=2018,
         journal="Journal of Geodesy 92, 349–360", publisher="Springer Nature", doi="10.1007/s00190-017-1066-2", role="Modern neutral-atmosphere reference, revised hydrostatic zenith formula and future VMF/GPT provider extension (VMF3/GPT3 not implemented)",
         urls=["https://link.springer.com/content/pdf/10.1007/s00190-017-1066-2.pdf"]),
    dict(id="yan-2021-ionosphere-propagation", title="Revisiting the light time correction in gravimetric missions like GRACE and GRACE follow-on", authors="Y. Yan, V. Müller, G. Heinzel and M. Zhong", year=2021,
         journal="Journal of Geodesy 95, 48", publisher="Springer Nature", doi="10.1007/s00190-021-01498-5", role="Published first-order plasma phase advancement and code/group delay sign, section 3.3; satellite geometry model is not used by yi-corr",
         urls=["https://link.springer.com/content/pdf/10.1007/s00190-021-01498-5.pdf"]),
    dict(id="hofmeister-boehm-2017-neutral-atmosphere", title="Application of ray-traced tropospheric slant delays to geodetic VLBI analysis", authors="A. Hofmeister and J. Böhm", year=2017,
         journal="Journal of Geodesy 91, 945–964", publisher="Springer Nature", doi="10.1007/s00190-017-1000-7", role="Neutral-atmosphere slant/zenith paths and limits of mapping approximations; future ray-tracing provider reference (ray tracing not implemented)",
         urls=["https://link.springer.com/content/pdf/10.1007/s00190-017-1000-7.pdf"]),
    dict(id="soja-2015-vlbi-wet-delay", title="Tropospheric delay determination by Kalman filtering VLBI data", authors="B. Soja et al.", year=2015,
         journal="Earth, Planets and Space 67, 144", publisher="Springer Nature / SpringerOpen", doi="10.1186/s40623-015-0293-0", role="Temporal variability of neutral-atmosphere wet delay and observational input; Kalman estimation is not implemented",
         urls=["https://link.springer.com/content/pdf/10.1186/s40623-015-0293-0.pdf"]),
]


def fetch(paper):
    paper = dict(paper)
    target = DEST / (paper["id"] + ".pdf")
    metadata = target.with_suffix(".json")
    if target.exists() and metadata.exists():
        saved = json.loads(metadata.read_text())
        if hashlib.sha256(target.read_bytes()).hexdigest() != saved["sha256"]:
            raise ValueError(f"reference changed: {target}")
        saved.setdefault("edition", "author manuscript of a journal-published paper" if "arxiv.org" in saved["resolved_url"] else "published article PDF (publisher or institutional archive)")
        metadata.write_text(json.dumps(saved, ensure_ascii=False, indent=2) + "\n")
        print(f"verified {target.name}", flush=True)
        return saved
    failures = []
    for url in paper["urls"]:
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (fx-corr scientific reference archive)"})
            with urllib.request.urlopen(request, timeout=25) as response:
                body = response.read()
                resolved = response.url
            if not body.startswith(b"%PDF-"):
                raise ValueError("response is not a PDF")
            target.write_bytes(body)
            edition = "author manuscript of a journal-published paper" if "arxiv.org" in resolved else "published article PDF (publisher or institutional archive)"
            paper.update(file=target.name, edition=edition, bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), download_url=url, resolved_url=resolved, downloaded_utc=datetime.now(timezone.utc).isoformat(), failures=failures, status="downloaded")
            metadata.write_text(json.dumps(paper, ensure_ascii=False, indent=2) + "\n")
            print(f"saved {target.name}: {len(body):,} bytes from {resolved}", flush=True)
            return paper
        except Exception as exc:
            failures.append(dict(url=url, error=str(exc)))
    paper.update(status="unavailable", failures=failures)
    print(f"unavailable {paper['id']}: {failures}", flush=True)
    return paper


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        papers = list(pool.map(fetch, PAPERS))
    (DEST.parent / "manifest.json").write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(), papers=papers), ensure_ascii=False, indent=2) + "\n")
    (DEST.parent / "SHA256SUMS").write_text("".join(f"{p['sha256']}  literature/{p['file']}\n" for p in papers if p["status"] == "downloaded"))
    print(f"downloaded {sum(p['status'] == 'downloaded' for p in papers)}/{len(papers)} published papers")


if __name__ == "__main__":
    main()
