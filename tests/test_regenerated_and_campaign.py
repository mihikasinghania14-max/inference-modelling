"""Tests for the code and data added in the regeneration step.

The Monte Carlo tests use small sizes; the data tests check that the stored regenerated files
are internally consistent and agree with exact results, not that they match any paper.
"""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rsi_shuffle import grid, theory

ROOT = Path(__file__).resolve().parents[1]
REGEN = ROOT / "data" / "extension" / "regenerated"


def _load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- grid with shared decks ---------------------------------------------------------------

def test_share_decks_same_layout_and_distribution():
    kw = dict(sizes=[16, 32], taus=[0.0, 1.0], R=400, m_values={16: [4, 16], 32: [8, 32]})
    a = grid.run_grid(seed=1, share_decks=False, **kw)
    b = grid.run_grid(seed=1, share_decks=True, **kw)
    assert [(r["N"], r["M"], r["tau"]) for r in a] == [(r["N"], r["M"], r["tau"]) for r in b]
    for r in b:
        if r["tau"] == 1.0 and r["M"] < r["N"]:
            exact = theory.estimator_variance(r["N"], r["M"])
            assert abs(r["variance"] / exact - 1) < 0.35          # R = 400: sd of ratio ~ 7%
            assert abs(r["mean_AD"] - float(theory.ad_inf(r["N"]))) < 5 * r["standard_error"] + 0.05


def test_grid_generator_arguments_are_not_exhausted():
    rows = grid.run_grid(sizes=iter([8]), taus=iter([0.0, 1.0]), R=50, seed=0,
                         m_values={8: [2, 8]}, share_decks=False)
    assert len(rows) == 4


# ---- exact Curie-Weiss moments --------------------------------------------------------------

@pytest.fixture(scope="module")
def cw():
    return _load(ROOT / "experiments" / "extension" / "campaign" / "curie_weiss_exact.py")


def test_cw_high_temperature_limits(cw):
    r = cw.moments(20000, 3.0)
    assert r["R_absM"] == pytest.approx(math.pi / 2 - 1, rel=0.01)    # |Gaussian| ratio
    assert r["chi_M"] == pytest.approx(1 / (1 - 1 / 3.0), rel=0.01)   # chi = 1/(1 - 1/T)


def test_cw_low_temperature_floor_vanishes(cw):
    assert cw.moments(5000, 0.7)["R_absM"] < 1e-3


def test_cw_critical_scaling(cw):
    a, b = cw.moments(2000, 1.0)["var_M"], cw.moments(32000, 1.0)["var_M"]
    assert math.log(b / a) / math.log(16) == pytest.approx(-0.5, abs=0.03)


# ---- stored regenerated Ising results ------------------------------------------------------

@pytest.mark.skipif(not (REGEN / "ising_critical_sizes.csv").exists(), reason="regenerated data absent")
def test_ising_regenerated_internal_consistency():
    d = pd.read_csv(REGEN / "ising_critical_sizes.csv")
    assert d["gate_pass"].all()
    assert (d["L"] ** 2 == d["n"]).all()
    assert d["R_E"].is_monotonic_decreasing
    # energy floor: effective exponent close to the exact Kaufman slope (-0.8637 for 16 <= L <= 128)
    s = json.loads((REGEN / "ising_critical_slopes.json").read_text())
    lo, hi = s["bootstrap_95ci"]["R_E"]
    assert lo < -0.8637 < hi or abs(s["slopes_vs_log_n"]["R_E"] + 0.8637) < 0.01
    # magnetisation floor is flat
    assert abs(s["slopes_vs_log_n"]["R_M"]) < 0.02
    assert d.loc[d["L"] >= 16, "R_M"].between(0.068, 0.076).all()
    # Binder cumulant stays near the critical value 0.8559 for L >= 16
    assert np.allclose(d.loc[d["L"] >= 16, "U4"], 0.8559, atol=0.005)
