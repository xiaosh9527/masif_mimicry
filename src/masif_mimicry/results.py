from __future__ import annotations

import csv
import html
from pathlib import Path


def find_result_csvs(run_dir: Path) -> list[Path]:
    return sorted(run_dir.glob("**/*.csv"))


def build_report(run_dir: Path) -> Path:
    csv_paths = find_result_csvs(run_dir)
    rows: list[dict[str, str]] = []
    for path in csv_paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                row["source_file"] = str(path.relative_to(run_dir))
                rows.append(row)
    preferred = ["P1_id", "P2_id", "MaSIF-score", "ca_clash", "heavy_atom_clash"]
    columns = preferred + sorted({key for row in rows for key in row} - set(preferred))
    table = "<p>No result CSV files were found yet.</p>"
    if rows:
        head = "".join(f"<th>{html.escape(col)}</th>" for col in columns)
        body = "".join(
            "<tr>" + "".join(f"<td>{html.escape(row.get(col, ''))}</td>" for col in columns) + "</tr>"
            for row in rows
        )
        table = f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>MaSIF-mimicry results</title>
<style>body{{font:15px system-ui;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#18212b}}
table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:.5rem;border:1px solid #d8dee4;text-align:left}}
th{{background:#f4f7fa;position:sticky;top:0}}tbody tr:nth-child(even){{background:#fafbfc}}</style></head>
<body><h1>MaSIF-mimicry results</h1><p>{len(rows)} candidates from {len(csv_paths)} files.</p>{table}</body></html>
"""
    output = run_dir / "report.html"
    output.write_text(page, encoding="utf-8")
    return output
