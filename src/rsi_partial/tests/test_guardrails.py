"""Structural epistemic-guardrail tests + storage/scoring/figure smoke.

The most important is the IMPORT BOUNDARY: no measuring layer may import the theory
module ``analysis.predictions`` -- estimators measure, predictions live only in analysis
(design brief; architecture design-invariant 2).  A violation here means a theoretical value
could leak into an estimator.
"""

from __future__ import annotations

import ast
from pathlib import Path

import numpy as np

import rsi_partial

PKG_ROOT = Path(rsi_partial.__file__).parent
MEASURING_LAYERS = ("systems", "observables", "observe", "inference", "measure")


def _imports_predictions(py_path: Path) -> bool:
    tree = ast.parse(py_path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if "analysis.predictions" in node.module or node.module.endswith("predictions"):
                if "analysis" in node.module:
                    return True
        if isinstance(node, ast.Import):
            for a in node.names:
                if "analysis.predictions" in a.name:
                    return True
    return False


def test_no_measuring_layer_imports_predictions():
    """Grep-able separation: nothing in the measuring layers imports analysis.predictions."""
    offenders = []
    for layer in MEASURING_LAYERS:
        for py in (PKG_ROOT / layer).rglob("*.py"):
            if _imports_predictions(py):
                offenders.append(str(py.relative_to(PKG_ROOT)))
    assert not offenders, f"theory leaked into measuring layers: {offenders}"


def test_metadata_returns_none_where_unknown():
    """Ising exposes no faked equilibrium metadata (measured, not analytic)."""
    from rsi_partial.systems.ising2d import Ising2D
    ise = Ising2D()
    assert ise.O_infinity(1024) is None
    assert ise.t_mix(1024) is None


def test_cards_offequilibrium_rho2_not_faked():
    """Cards.rho2 returns the exact equilibrium 2 but None off-equilibrium (empirical form
    lives in analysis, not the system)."""
    from rsi_partial.systems.cards import Cards
    assert Cards().rho2(1.0) == 2.0
    assert Cards().rho2(0.6) is None


def test_scoring_logic():
    from rsi_partial.analysis.bootstrap import CI
    from rsi_partial.analysis.score import FALSIFY, INCONCLUSIVE, SUPPORT, score_ci
    assert score_ci(CI(0.0, -0.05, 0.05), support=(-0.1, 0.1), falsify=(0.25, np.inf)) == SUPPORT
    assert score_ci(CI(0.4, 0.35, 0.45), support=(-0.1, 0.1), falsify=(0.25, np.inf)) == FALSIFY
    assert score_ci(CI(0.2, 0.05, 0.35), support=(-0.1, 0.1), falsify=(0.25, np.inf)) == INCONCLUSIVE


def test_bootstrap_ci_covers_truth():
    rng = np.random.default_rng(0)
    data = rng.normal(5.0, 1.0, size=2000)
    from rsi_partial.analysis.bootstrap import bootstrap_ci
    ci = bootstrap_ci(data, np.mean, n_boot=500, rng=rng)
    assert ci.lo < 5.0 < ci.hi


def test_storage_roundtrip_and_no_silent_overwrite(tmp_path):
    from rsi_partial.storage.io import read_records, write_records
    recs = [{"experiment_id": "t", "config_hash": "abcd", "system": "cards",
             "observable": "AD", "ground_truth": "exact_marginal", "n": 512,
             "regime": "equilibrium", "phi_target": 1.0, "m": 32, "replicate": r,
             "seed_128": "00", "code_version": "0.1.0", "Ohat_m": float(r)} for r in range(5)]
    p = write_records(recs, tmp_path, "t", "abcd")
    df = read_records(p)
    assert len(df) == 5 and "Ohat_m" in df.columns
    import pytest
    with pytest.raises(FileExistsError):
        write_records(recs, tmp_path, "t", "abcd")          # no silent overwrite


def test_figures_render(tmp_path):
    from rsi_partial.analysis import figures as F
    fig = F.fig_precision_vs_n([128, 512, 2048], [2.01, 1.99, 2.00], target=2.0)
    out = F.save(fig, str(tmp_path / "prec.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0
