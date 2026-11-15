"""Auto-scoring: turn a measured value + bootstrap CI into a support/falsify/inconclusive
verdict against the pre-registered bands in :mod:`rsi_partial.analysis.predictions`.

Verdict rule (Phase-4 sec.4.3): a prediction "passes" only if its CI lies inside the
support band; "fails" only if its CI lies inside the falsification band; otherwise
"inconclusive -> increase R".  This is the ONLY place measurement meets theory.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .bootstrap import CI
from .predictions import PREDICTIONS, Prediction

SUPPORT = "support"
FALSIFY = "falsify"
INCONCLUSIVE = "inconclusive"


def _ci_within(ci: CI, band: tuple[float, float] | None) -> bool:
    if band is None:
        return False
    lo, hi = band
    return ci.lo >= lo and ci.hi <= hi


@dataclass(frozen=True)
class Verdict:
    prediction_id: str
    status: str                 # epistemic tag of the prediction (EST/EMP/CONJ/PCV)
    verdict: str                # support | falsify | inconclusive
    quantity: str
    measured: float
    ci_lo: float
    ci_hi: float
    support_band: tuple[float, float] | None
    falsify_band: tuple[float, float] | None
    note: str = ""

    def as_row(self) -> dict[str, Any]:
        return {
            "prediction": self.prediction_id,
            "status": self.status,
            "verdict": self.verdict,
            "quantity": self.quantity,
            "measured": self.measured,
            "ci_lo": self.ci_lo,
            "ci_hi": self.ci_hi,
            "support_band": self.support_band,
            "falsify_band": self.falsify_band,
            "note": self.note,
        }


def score_ci(
    ci: CI,
    *,
    support: tuple[float, float] | None,
    falsify: tuple[float, float] | None,
) -> str:
    """Core verdict logic for a CI against support / falsify bands."""
    in_support = _ci_within(ci, support)
    in_falsify = _ci_within(ci, falsify)
    if in_support and not in_falsify:
        return SUPPORT
    if in_falsify and not in_support:
        return FALSIFY
    return INCONCLUSIVE


def score_prediction(
    prediction_id: str,
    ci: CI,
    *,
    support: tuple[float, float] | None = None,
    falsify: tuple[float, float] | None = None,
    note: str = "",
) -> Verdict:
    """Score one prediction.  ``support``/``falsify`` override the registered bands (needed
    for P2, whose band depends on the observable's universality-class ``a``)."""
    pred: Prediction = PREDICTIONS[prediction_id]
    sup = support if support is not None else pred.support
    fal = falsify if falsify is not None else pred.falsify
    verdict = score_ci(ci, support=sup, falsify=fal)
    return Verdict(
        prediction_id=prediction_id, status=pred.status, verdict=verdict,
        quantity=pred.quantity, measured=ci.value, ci_lo=ci.lo, ci_hi=ci.hi,
        support_band=sup, falsify_band=fal, note=note or pred.notes,
    )


def verdict_table(verdicts: list[Verdict]) -> list[dict[str, Any]]:
    """Machine-readable verdict table (list of rows); ready for DataFrame / JSON."""
    return [v.as_row() for v in verdicts]
