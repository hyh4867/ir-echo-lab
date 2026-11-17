"""Rewrite the hash table in docs/DATA.md from data/manifest.json.

    python scripts/update_data_doc.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config  # noqa: E402

START, END = "<!-- manifest:start -->", "<!-- manifest:end -->"
GROUPED_PREFIXES = ("s5/neowise/", "s5/allwise_mep/", "s5/untimely/", "s5/ztf/", "s6/asassn/", "s6/bat157/",
                    "s6b/neowise/", "s6b/allwise_mep/", "s6b/untimely/",
                    "allwise_cat/s5_targets_", "allwise_cat/s5_fields_",
                    "allwise_cat/s6b_targets_", "allwise_cat/s6b_fields_", "catalogs/")


def main() -> None:
    manifest = json.loads(config.MANIFEST_PATH.read_text())
    lines = ["| 캐시 파일 (data/cache/ 아래) | 행 수 | 받은 시각 (UTC) | SHA-256 |", "|---|---|---|---|"]
    # Stage 5a files (one per field and survey, plus batched catalogue queries)
    # are too many to list; they are summarised per group below the table and
    # every one of them is in data/manifest.json.
    groups: dict[str, list[dict]] = {}
    for key in sorted(manifest):
        m = manifest[key]
        group = next((g for g in GROUPED_PREFIXES if key.startswith(g)), None)
        if group is not None:
            groups.setdefault(group, []).append(m)
            continue
        lines.append(f"| `{key}` | {m['n_rows']} | {m['retrieved_utc']} | `{m['sha256']}` |")
    total_rows = sum(m["n_rows"] for m in manifest.values())
    total_bytes = sum(m.get("bytes", 0) for m in manifest.values())
    lines.append("")
    lines.append(f"파일 {len(manifest)}개, 모두 {total_rows:,}행, {total_bytes / 1e6:.2f} MB.")
    if groups:
        lines += ["", "5a·5b단계의 파일은 수가 많아 묶음별 합계만 적는다. 파일별 조회 조건, 시각, 행 수, "
                  "SHA-256은 `data/manifest.json`에 있다.", "",
                  "| 묶음 (data/cache/ 아래) | 파일 수 | 행 수 | 크기 (MB) | 받은 시각 (UTC) |",
                  "|---|---|---|---|---|"]
        for group in sorted(groups):
            ms = groups[group]
            times = sorted(m["retrieved_utc"] for m in ms)
            lines.append(f"| `{group}*` | {len(ms)} | {sum(m['n_rows'] for m in ms):,} | "
                         f"{sum(m.get('bytes', 0) for m in ms) / 1e6:.1f} | {times[0]} – {times[-1]} |")
    doc_path = config.PROJECT_ROOT / "docs" / "DATA.md"
    doc = doc_path.read_text()
    head, rest = doc.split(START)
    _, tail = rest.split(END)
    doc_path.write_text(head + START + "\n" + "\n".join(lines) + "\n" + END + tail)
    print(f"updated {doc_path.name}: {len(manifest)} files")


if __name__ == "__main__":
    main()
