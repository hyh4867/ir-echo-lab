"""Stage 7 (step 5): what the literature and public catalogues say about the
X-ray emission of CGCG 312-012. Run after docs/XRAY_PLAN.md was committed.

    python scripts/xray_literature.py

1. SIMBAD reference list of the object (titles), marked for X-ray keywords.
2. VizieR: rows of the object in the BASS X-ray catalogue (Ricci et al. 2017,
   J/ApJS/233/17) and in Swift-XRT source catalogues (2SXPS, IX/58).
3. The object's line in the Swift-BAT 157-month catalogue (already downloaded).
Output: results/s7_literature.json. No login is used.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, xray  # noqa: E402

SIMBAD = "https://simbad.harvard.edu/simbad/sim-tap/sync"
VIZIER_ASU = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"  # the CDS TAP service timed out on 2026-11-17
KEYWORDS = ("x-ray", "xray", "nustar", "swift", "bat ", "bass", "hard x", "xmm", "chandra", "rosat", "absor",
            "obscur", "compton", "70 month", "70-month", "105-month", "105 month", "flare", "transient", "tidal",
            "changing", "echo", "outburst", "variab")
HEAD = {"User-Agent": "irecho/0.1 (python-requests)"}


def tap(url: str, query: str, timeout: float = 120.0) -> list[dict]:
    resp = requests.post(url, data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "tsv", "QUERY": query},
                         timeout=timeout, headers=HEAD)
    resp.raise_for_status()
    lines = [ln for ln in resp.text.splitlines() if ln.strip()]
    cols = [c.strip('"') for c in lines[0].split("\t")]
    return [dict(zip(cols, (v.strip().strip('"') for v in ln.split("\t")))) for ln in lines[1:]]


def asu(params: dict) -> dict:
    """VizieR ASU query; returns {table name: [rows]} with empty columns dropped."""
    resp = requests.get(VIZIER_ASU, params=dict(params, **{"-out.all": "", "-out.max": "20"}), timeout=120, headers=HEAD)
    resp.raise_for_status()
    tables, name, block = {}, None, []

    def flush():
        if name and len(block) >= 4:
            cols = block[0].split("\t")
            rows = [dict(zip(cols, (v.strip() for v in ln.split("\t")))) for ln in block[3:]]
            tables[name] = [{k: v for k, v in r.items() if v} for r in rows]

    for ln in resp.text.splitlines():
        if ln.startswith("#Name:") and "/" in ln and name is None or ln.startswith("#Table"):
            flush()
            block = []
            name = ln.split()[-1].rstrip(":") if ln.startswith("#Table") else name
        elif ln.strip() and not ln.startswith("#"):
            block.append(ln)
    flush()
    return tables


def attempt(label: str, fn):
    try:
        return {"status": "ok", "rows": fn()}
    except Exception as exc:  # noqa: BLE001
        print(label, "failed:", type(exc).__name__, str(exc)[:200])
        return {"status": f"failed: {type(exc).__name__}", "rows": []}


def main() -> None:
    out = {"date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "object": xray.TARGET["counterpart"],
           "keywords": list(KEYWORDS)}
    refs = attempt("simbad", lambda: tap(
        SIMBAD, "SELECT ref.bibcode, ref.\"year\", ref.title FROM ref JOIN has_ref ON ref.oidbib = has_ref.oidbibref "
                "JOIN ident ON has_ref.oidref = ident.oidref WHERE ident.id = 'Z 312-12'"))
    rows = sorted(refs["rows"], key=lambda r: (r.get("year", ""), r.get("bibcode", "")))
    for r in rows:
        r["keyword_hit"] = [k.strip() for k in KEYWORDS if k in r.get("title", "").lower()]
    out["simbad_references"] = {"service": SIMBAD, "identifier": "Z 312-12", "status": refs["status"],
                                "n": len(rows), "n_with_keywords": sum(1 for r in rows if r["keyword_hit"]),
                                "rows": rows}
    ra, dec = xray.TARGET["ra_deg"], xray.TARGET["dec_deg"]
    out["vizier"] = {"service": VIZIER_ASU}
    for label, params in (
            ("bass_xray_ricci2017_by_position", {"-source": "J/ApJS/233/17", "-c": f"{ra:.4f}+{dec:.4f}", "-c.rs": "60"}),
            ("bass_xray_ricci2017_by_name", {"-source": "J/ApJS/233/17", "SWIFT": "J0929.7+6232"}),
            ("2sxps", {"-source": "IX/58/2sxps", "-c": f"{ra:.4f}+{dec:.4f}", "-c.rs": "30"}),
            ("lsxps", {"-source": "IX/69", "-c": f"{ra:.4f}+{dec:.4f}", "-c.rs": "30"})):
        res = attempt(label, lambda q=params: asu(q))
        out["vizier"][label] = {"params": params, "status": res["status"], "tables": res["rows"]}
    cat = xray.S7_DIR / "bat157" / "BAT_157m_catalog_20250525.txt"
    if cat.exists():
        lines = cat.read_text(errors="replace").splitlines()
        out["bat157_catalogue"] = {"file": "data/cache/s7/bat157/" + cat.name,
                                   "header": [ln for ln in lines[:40] if ln.startswith(("#", "|")) or "Name" in ln][:6],
                                   "row": [ln for ln in lines if "J0929.7+6232" in ln]}
    (config.RESULTS_DIR / "s7_literature.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print("simbad", out["simbad_references"]["n"], "refs,", out["simbad_references"]["n_with_keywords"], "with keywords")
    for r in rows:
        if r["keyword_hit"]:
            print(" ", r["bibcode"], r["title"][:110], r["keyword_hit"])
    for label, rec in out["vizier"].items():
        if isinstance(rec, dict):
            for t, rows_ in rec["tables"].items():
                print(label, t, len(rows_))
                for r in rows_:
                    print("   ", json.dumps(r, ensure_ascii=False)[:1200])
    print(out.get("bat157_catalogue", {}).get("row"))


if __name__ == "__main__":
    main()
