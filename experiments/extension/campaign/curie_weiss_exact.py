"""Exact equilibrium moments of the Curie-Weiss model, no sampling involved.

Spins s_i = +-1, i = 1..N, weight  exp( N m^2 / (2T) )  with  m = (1/N) sum s_i  (J = 1,
so T_c = 1).  For each N and T the distribution of m is summed exactly over the N + 1
magnetisation values, using log-binomials for stability.

Output columns:  T, N, var_M (= E M^2, the mean is 0), E_absM, var_absM,
                 R_absM = Var|M| / E|M|^2   (relative floor of |M|),
                 chi_M = N var_M.

The slope of log var_M against log N is a - 1 (a = 0 above T_c, 1/2 at T_c, 1 below).
R_absM tends to pi/2 - 1 above T_c, to sqrt(2) - 1 at T_c (see verification/cw.py) and to 0
below T_c, where |M| concentrates at the spontaneous magnetisation.

    python curie_weiss_exact.py OUT.csv
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.special import gammaln

SIZES = [100, 150, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 10000, 20000, 50000]
TEMPS = {"0.7 (below T_c)": 0.7, "0.9": 0.9, "1.0 (T_c)": 1.0, "1.5": 1.5, "2.0 (above T_c)": 2.0}


def moments(N: int, T: float) -> dict[str, float]:
    k = np.arange(N + 1)
    m = (2 * k - N) / N
    lw = gammaln(N + 1) - gammaln(k + 1) - gammaln(N - k + 1) + N * m**2 / (2 * T)
    w = np.exp(lw - lw.max()); w /= w.sum()
    a1, a2 = w @ np.abs(m), w @ m**2
    return {"var_M": float(a2), "E_absM": float(a1), "var_absM": float(a2 - a1**2),
            "R_absM": float((a2 - a1**2) / a1**2), "chi_M": float(N * a2)}


def main(out: str) -> None:
    rows = []
    for label, T in TEMPS.items():
        for N in SIZES:
            rows.append({"T": T, "regime": label, "N": N, **moments(N, T)})
    with open(out, "w") as fh:
        cols = ["T", "regime", "N", "var_M", "E_absM", "var_absM", "R_absM", "chi_M"]
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(f'"{r[c]}"' if c == "regime" else f"{r[c]:.10g}" for c in cols) + "\n")
    print(f"{'T':>5} {'slope log Var(M) vs log N (= a-1)':>36}  fit range")
    for label, T in TEMPS.items():
        s = [r for r in rows if r["T"] == T and 100 <= r["N"] <= 3000]
        sl = np.polyfit(np.log([r["N"] for r in s]), np.log([r["var_M"] for r in s]), 1)[0]
        print(f"{T:>5} {sl:>36.4f}  N = 100..3000")


if __name__ == "__main__":
    main(sys.argv[1])
