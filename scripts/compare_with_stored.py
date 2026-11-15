"""Compare two result directories written by ``rsi_partial.experiments.run``.

Numeric columns must agree exactly (same seed protocol => bit-identical); columns with
non-numeric content (None versus NaN encodings, strings) are compared after mapping
None -> NaN.  Manifests are compared on their configuration and seed only (host,
platform, library versions and timestamps legitimately differ between machines).

Usage:  python scripts/compare_with_stored.py STORED_DIR NEW_DIR
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def compare_parquet(a: Path, b: Path) -> list[str]:
    da, db = pd.read_parquet(a), pd.read_parquet(b)
    problems: list[str] = []
    if list(da.columns) != list(db.columns) or len(da) != len(db):
        return [f"{a.name}: shape/columns differ"]
    for c in da.columns:
        x, y = da[c], db[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            if not np.array_equal(x.to_numpy(float), y.to_numpy(float), equal_nan=True):
                problems.append(f"{a.name}:{c} numeric difference (max abs "
                                f"{np.nanmax(np.abs(x.to_numpy(float) - y.to_numpy(float))):.3g})")
        else:
            xs = x.astype(object).where(x.notna(), None)
            ys = y.astype(object).where(y.notna(), None)
            if not xs.equals(ys):
                problems.append(f"{a.name}:{c} differs")
    return problems


def compare_manifest(a: Path, b: Path) -> list[str]:
    ma, mb = json.loads(a.read_text()), json.loads(b.read_text())
    return [f"{a.name}:{k} differs" for k in ("experiment_id", "root_seed", "config_grid")
            if ma.get(k) != mb.get(k)]


def main(stored: str, new: str) -> int:
    stored, new = Path(stored), Path(new)
    problems: list[str] = []
    for p in sorted(stored.glob("*.parquet")):
        q = new / p.name
        problems += compare_parquet(p, q) if q.exists() else [f"{p.name}: missing in {new}"]
    for p in sorted(stored.glob("*manifest.json")):
        q = new / p.name
        problems += compare_manifest(p, q) if q.exists() else [f"{p.name}: missing in {new}"]
    for line in problems or ["identical (numeric columns bit-for-bit; manifests: config and seed)"]:
        print(line)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
