"""Stage 7: download the public X-ray data of CGCG 312-012 from HEASARC.

    python scripts/fetch_xray.py [--list]

No login is needed. Files go to data/cache/s7/ (not in git); URL, size, time
and SHA-256 of every file are recorded in data/manifest.json.

  nustar/60160365002/{auxil,hk,event_uf,event_cl}   NuSTAR, 2022-05-01
  swift/00040946001/{auxil,xrt}                     Swift XRT, 2010-12-29
  swift/00081051001/{auxil,xrt}                     Swift XRT, 2022-05-01
  bat157/                                           BAT 157-month spectrum, response, catalogue
"""

from __future__ import annotations

import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, irsa, xray  # noqa: E402

HEADERS = {"User-Agent": "irecho/0.1 (python-requests)"}
HREF = re.compile(r'href="([^"?/][^"]*)"')


def listing(session, url: str) -> list[str]:
    resp = session.get(url, timeout=120, headers=HEADERS)
    resp.raise_for_status()
    return [h for h in HREF.findall(resp.text) if not h.startswith(("http", "mailto"))]


def crawl(session, url: str, rel: str = "") -> list[tuple[str, str]]:
    """(relative path, url) of every file below a directory listing."""
    out = []
    for name in listing(session, url):
        if name.endswith("/"):
            out += crawl(session, url + name, rel + name)
        else:
            out.append((rel + name, url + name))
    return out


def download(session, url: str, path: Path, key: str) -> dict:
    meta_path = path.with_suffix(path.suffix + ".meta.json")
    if path.exists() and meta_path.exists():
        return json.loads(meta_path.read_text())
    for attempt in range(4):
        try:
            resp = session.get(url, timeout=600, headers=HEADERS)
            resp.raise_for_status()
            break
        except Exception as exc:  # noqa: BLE001
            if attempt == 3:
                raise
            print("retry", url, exc, flush=True)
            time.sleep(5 * (attempt + 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(resp.content)
    meta = {"url": url, "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "n_rows": 0, "sha256": irsa.sha256_bytes(resp.content), "bytes": len(resp.content)}
    meta_path.write_text(json.dumps(meta, indent=1) + "\n")
    irsa._update_manifest(key, meta, config.MANIFEST_PATH)
    return meta


def main() -> None:
    import requests

    only_list = "--list" in sys.argv
    session = requests.Session()
    jobs = []
    for group, base, subdirs in xray.ARCHIVE_DIRS:
        for sub in subdirs:
            for rel, url in crawl(session, f"{base}{sub}/", f"{sub}/"):
                jobs.append((f"{group}/{rel}", url))
    for name, url in xray.BAT157_FILES.items():
        jobs.append((f"bat157/{name}", url))
    print(len(jobs), "files", flush=True)
    total = 0
    for i, (rel, url) in enumerate(jobs):
        if only_list:
            print(rel)
            continue
        meta = download(session, url, xray.S7_DIR / rel, f"s7/{rel}")
        total += meta["bytes"]
        if i % 20 == 0 or i == len(jobs) - 1:
            print(i + 1, len(jobs), rel, f"{total / 1e6:.1f} MB", flush=True)
    print("done", f"{total / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
