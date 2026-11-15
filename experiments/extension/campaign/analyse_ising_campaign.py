"""Summarise the critical-Ising campaign written by run_ising_campaign.sh.

For every size L the 100 block means of each of the two independent starts are read.
Reported per size (blocks of both starts pooled):

  R[|M|]  = Var|M| / E|M|^2        relative floor of the magnetisation magnitude
  R[E]    = Var E / E^2            relative floor of the energy density
  chi[M]  = n Var|M|               susceptibility of |M|   (n = L^2)
  U4      = <M^2>^2 / <M^4>        Binder ratio

with delete-one-block jackknife errors, and an equilibration gate: the two starts must
agree on R[|M|] and R[E] within 4 combined standard errors.

Slopes d log R / d log n are fitted over sizes L >= 16 with a block-bootstrap 95% interval
(resampling the blocks of every size and every start, refitting the weighted-free
least-squares line).

    python analyse_ising_campaign.py BLOCK_DIR OUT_PREFIX
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

CRITICAL_FIT_MIN_L = 16


def stats(v):
    m1, m2, m4, e1, e2 = v
    return {"R_M": (m2 - m1**2) / m1**2, "R_E": (e2 - e1**2) / e1**2, "var_M": m2 - m1**2,
            "U4": m2**2 / m4, "mean_absM": m1, "mean_E": e1}


def jackknife(blocks: np.ndarray, key: str):
    nb = len(blocks)
    full = stats(blocks.mean(0))[key]
    loo = np.array([stats(np.delete(blocks, i, 0).mean(0))[key] for i in range(nb)])
    return full, float(np.sqrt((nb - 1) / nb * np.sum((loo - loo.mean()) ** 2)))


def load(block_dir: Path):
    data: dict[int, dict[int, np.ndarray]] = {}
    for f in sorted(block_dir.glob("L*_start*.dat")):
        L, st = map(int, re.match(r"L(\d+)_start(\d)\.dat", f.name).groups())
        data.setdefault(L, {})[st] = np.loadtxt(f)
    return data


def slope(ns, ys):
    return float(np.polyfit(np.log(ns), np.log(ys), 1)[0])


def main(block_dir: str, prefix: str) -> None:
    data = load(Path(block_dir))
    sizes = sorted(data)
    rows, pooled = [], {}
    rng = np.random.default_rng(20261115)
    for L in sizes:
        starts = data[L]
        allb = np.vstack(list(starts.values()))
        pooled[L] = starts
        row = {"L": L, "n": L * L, "blocks": len(allb)}
        for key in ("R_M", "R_E", "U4"):
            row[key], row[key + "_se"] = jackknife(allb, key)
        row["chi_M"] = L * L * stats(allb.mean(0))["var_M"]
        row["mean_absM"] = float(allb[:, 0].mean()); row["mean_E"] = float(allb[:, 3].mean())
        gate = {}
        if len(starts) == 2:
            a, b = starts[0], starts[1]
            for key in ("R_M", "R_E"):
                (ra, sa), (rb, sb) = jackknife(a, key), jackknife(b, key)
                gate[key] = abs(ra - rb) / np.hypot(sa, sb)
        row["gate_z_R_M"] = gate.get("R_M", float("nan")); row["gate_z_R_E"] = gate.get("R_E", float("nan"))
        row["gate_pass"] = bool(all(z < 4 for z in gate.values())) if gate else False
        rows.append(row)

    fit = [r for r in rows if r["L"] >= CRITICAL_FIT_MIN_L]
    ns = np.array([r["n"] for r in fit], float)
    point = {"R_M": slope(ns, [r["R_M"] for r in fit]), "R_E": slope(ns, [r["R_E"] for r in fit]),
             "chi_M": slope(ns, [r["chi_M"] for r in fit])}
    boot = {k: [] for k in point}
    for _ in range(2000):
        vals = {k: [] for k in point}
        for r in fit:
            blocks = np.vstack(list(pooled[r["L"]].values()))
            res = blocks[rng.integers(0, len(blocks), len(blocks))].mean(0)
            s = stats(res)
            vals["R_M"].append(s["R_M"]); vals["R_E"].append(s["R_E"]); vals["chi_M"].append(r["n"] * s["var_M"])
        for k in point:
            boot[k].append(slope(ns, vals[k]))
    ci = {k: [float(np.percentile(boot[k], 2.5)), float(np.percentile(boot[k], 97.5))] for k in point}

    out = Path(prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["L", "n", "blocks", "R_M", "R_M_se", "R_E", "R_E_se", "chi_M", "U4", "U4_se", "mean_absM",
            "mean_E", "gate_z_R_M", "gate_z_R_E", "gate_pass"]
    with open(f"{out}_sizes.csv", "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    summary = {"fit_sizes_L_ge": CRITICAL_FIT_MIN_L, "slopes_vs_log_n": point,
               "bootstrap_95ci": ci, "a_M_estimate": point["chi_M"], "bootstrap_resamples": 2000,
               "note": "regenerated with experiments/extension/campaign; not the original runs"}
    Path(f"{out}_slopes.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(open(f"{out}_sizes.csv").read())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
