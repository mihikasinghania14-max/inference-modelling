"""Result I/O: Parquet (CSV fallback), deterministic filenames, JSON manifests.

Guarantees (design brief sec.7):
  * deterministic filenames keyed by ``experiment_id`` + config hash;
  * nothing overwrites silently (``write_records`` refuses an existing path unless
    ``overwrite=True``);
  * every run writes a JSON manifest sidecar (root seed, code version, library versions,
    host, UTC timestamp, full config grid).
"""

from __future__ import annotations

import json
import platform
import socket
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .. import __version__
from .schema import validate_record


def _records_to_frame(records: Iterable[dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame([validate_record(r) for r in records])


def result_path(out_dir: str | Path, experiment_id: str, config_hash: str,
                *, ext: str = "parquet") -> Path:
    return Path(out_dir) / f"{experiment_id}__{config_hash}.{ext}"


def write_records(
    records: Iterable[dict[str, Any]],
    out_dir: str | Path,
    experiment_id: str,
    config_hash: str,
    *,
    prefer_parquet: bool = True,
    overwrite: bool = False,
) -> Path:
    """Write per-(config,replicate) records; return the path actually written.

    Falls back to CSV if pyarrow/parquet is unavailable.  Refuses to clobber an existing
    file unless ``overwrite=True``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = _records_to_frame(records)

    ext = "parquet" if prefer_parquet else "csv"
    path = result_path(out_dir, experiment_id, config_hash, ext=ext)
    if path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite {path} (pass overwrite=True)")
    if prefer_parquet:
        try:
            frame.to_parquet(path, index=False)
            return path
        except Exception:
            path = result_path(out_dir, experiment_id, config_hash, ext="csv")
    if path.exists() and not overwrite:
        raise FileExistsError(f"refusing to overwrite {path} (pass overwrite=True)")
    frame.to_csv(path, index=False)
    return path


def read_records(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def library_versions() -> dict[str, str]:
    import scipy
    vers = {"numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__}
    try:
        import pyarrow
        vers["pyarrow"] = pyarrow.__version__
    except Exception:
        vers["pyarrow"] = "unavailable"
    return vers


def write_manifest(
    out_dir: str | Path,
    experiment_id: str,
    *,
    root_seed: int,
    config_grid: Any,
    extra: dict[str, Any] | None = None,
) -> Path:
    """Write a JSON manifest sidecar for a run (fully reconstructible provenance)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "experiment_id": experiment_id,
        "root_seed": int(root_seed),
        "code_version": __version__,
        "library_versions": library_versions(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_grid": config_grid,
    }
    if extra:
        manifest["extra"] = extra
    path = out_dir / f"{experiment_id}__manifest.json"
    path.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    return path
