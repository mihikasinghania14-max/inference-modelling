"""Summarise the block means written by ``wolff.c``.

``wolff`` prints 100 lines ``mean|M|  mean M^2  mean M^4`` (one per block of equal
length).  This script reports the relative floor  R[|M|] = Var(|M|)/E[|M|]^2  and the
Binder ratio  <M^2>^2/<M^4>, with delete-one-block jackknife standard errors.

Usage:  python wolff_analyze.py out.dat [out2.dat ...]
"""
import sys

import numpy as np


def jackknife(blocks: np.ndarray, stat):
    nb = len(blocks)
    full = stat(blocks.mean(axis=0))
    loo = np.array([stat(np.delete(blocks, i, axis=0).mean(axis=0)) for i in range(nb)])
    se = np.sqrt((nb - 1) / nb * np.sum((loo - loo.mean()) ** 2))
    return full, se


def r_abs_m(v):      # v = (E|M|, E M^2, E M^4)
    return (v[1] - v[0] ** 2) / v[0] ** 2


def binder(v):
    return v[1] ** 2 / v[2]


if __name__ == "__main__":
    for path in sys.argv[1:]:
        b = np.loadtxt(path)
        r, rse = jackknife(b, r_abs_m)
        u, use = jackknife(b, binder)
        print(f"{path}: blocks={len(b)}  R[|M|]={r:.4f} +/- {rse:.4f}   <M^2>^2/<M^4>={u:.4f} +/- {use:.4f}")
