"""Export every sheet of the spreadsheets under data/original-research to CSV.

Cell values are the values cached in the workbook (formulas are not re-evaluated).
The CSV files are written next to each workbook in a ``csv/`` sub-folder, named
``<workbook>__<sheet>.csv``.  The workbooks remain the authoritative copies.

Usage:  python scripts/export_xlsx_to_csv.py
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1] / "data" / "original-research"


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9=.]+", "_", s).strip("_")


def main() -> None:
    for wb_path in sorted(ROOT.rglob("*.xlsx")):
        out_dir = wb_path.parent / "csv"
        out_dir.mkdir(exist_ok=True)
        wb = openpyxl.load_workbook(wb_path, data_only=True)
        for ws in wb.worksheets:
            rows = [["" if v is None else v for v in r] for r in ws.iter_rows(values_only=True)]
            while rows and all(v == "" for v in rows[-1]):
                rows.pop()
            target = out_dir / f"{slug(wb_path.stem)}__{slug(ws.title)}.csv"
            with open(target, "w", newline="") as fh:
                csv.writer(fh).writerows(rows)
            print(f"{target.relative_to(ROOT.parent.parent)}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
