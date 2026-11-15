"""Canonical record schema for stored results (architecture sec.3).

One row per ``(config, replicate)``.  Aggregated summaries (one row per config, with
bootstrap CIs) are produced in the analysis layer.  ``RECORD_FIELDS`` is the source of
truth; :func:`validate_record` coerces + checks a dict against it.
"""

from __future__ import annotations

from typing import Any

# field -> (python type, nullable)
RECORD_FIELDS: dict[str, tuple[type, bool]] = {
    "experiment_id": (str, False),
    "config_hash": (str, False),
    "system": (str, False),
    "observable": (str, False),
    "ground_truth": (str, False),
    "n": (int, False),
    "L": (int, True),
    "T": (float, True),
    "T_over_Tc": (float, True),
    "regime": (str, False),
    "phi_target": (float, False),
    "theta": (float, True),
    "m": (int, False),
    "delta": (float, True),
    "epsilon": (float, True),
    "replicate": (int, False),
    "seed_128": (str, False),
    "code_version": (str, False),
    "O_full": (float, True),
    "Ohat_m": (float, False),
    "phi_hat": (float, True),
    "Vx": (float, True),
    "samp_var_raw": (float, True),
    "samp_var_fpc": (float, True),
    "skew": (float, True),
    "exkurt": (float, True),
    "post_mean": (float, True),
    "post_sd": (float, True),
    "ci_lo": (float, True),
    "ci_hi": (float, True),
    "reject_H0": (bool, True),
    "tau_int": (float, True),
    "R_eff": (float, True),
    "equilibrated_flag": (bool, True),
    "access_reads": (int, True),
    "rho2_hat": (float, True),
    "Rn_hat": (float, True),
}

#: Aggregated-summary fields (one row per config; each headline quantity + bootstrap CI).
SUMMARY_FIELDS: tuple[str, ...] = (
    "experiment_id", "config_hash", "system", "observable", "ground_truth",
    "n", "L", "T", "T_over_Tc", "regime", "phi_target", "m", "delta", "epsilon", "R",
    "var_Ohat", "sampling_term", "floor", "Var_nu_O", "chi_n", "R_n",
    "rho2_phi", "I_eff", "mstar", "mstar_kind",
    "value", "ci_lo", "ci_hi",           # generic headline value + bootstrap CI
)


def validate_record(rec: dict[str, Any], *, strict: bool = False) -> dict[str, Any]:
    """Coerce a record dict to the schema; fill missing nullable fields with ``None``.

    ``strict=True`` raises on unknown keys or non-nullable missing keys.
    """
    out: dict[str, Any] = {}
    for field, (typ, nullable) in RECORD_FIELDS.items():
        if field not in rec or rec[field] is None:
            if not nullable and strict:
                raise ValueError(f"missing non-nullable field {field!r}")
            out[field] = None
        else:
            v = rec[field]
            try:
                out[field] = None if v is None else typ(v)
            except (TypeError, ValueError):
                if strict:
                    raise
                out[field] = None
    if strict:
        unknown = set(rec) - set(RECORD_FIELDS)
        if unknown:
            raise ValueError(f"unknown fields: {sorted(unknown)}")
    return out
