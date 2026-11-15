"""Declarative, FROZEN Phase-5 config grids (Phase-4 sec.4.2) -- staged.

This is a data module: the exact experimental matrix, transcribed from the frozen design.
Runners consume it; nothing here computes.  Grids are exposed as plain structures so a
run is fully specified by (stage, smoke|full) + a root seed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# ---- shared grids (Phase-4 sec.4.2) ---------------------------------------

SIZES_EXACT = [2 ** 7, 2 ** 9, 2 ** 11, 2 ** 13, 2 ** 15]      # 128 .. 32768
SIZES_EXACT_VALIDATORS = [3000, 12000]                         # non-power-of-two controls
ISING_L = [16, 24, 32, 48, 64, 96]                            # + 128 for critical M-vs-E
ISING_L_CRIT_EXTRA = 128
MEANFIELD_N = [10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5]
S5_L = [256, 1024, 4096]
S5_XI = [2, 4, 8, 16, 32, 64]
S6_N = [512, 4096, 32768]

M_GEOMETRIC = [8, 16, 32, 64, 128, 256, 512, 1024]            # truncated to m <= n/2, +census
DELTAS = [0.20, 0.15, 0.10, 0.075, 0.05]
EPSILONS = [0.20, 0.15, 0.10, 0.075, 0.05]
ISING_T_OVER_TC = [2.0, 1.4, 1.15, 1.05, 1.00]               # + 0.85 broken-symmetry control
ISING_T_CONTROL = 0.85
THETA_GRID = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]        # off-equilibrium mixing levels

# ---- replication targets (Phase-4 sec.4.3) --------------------------------

R_VARIANCE = 2000        # variance / rho^2 configs
R_DETECTION = 5000       # m* / detection-threshold configs
R_ISING_LARGE = 400      # expensive Ising L>=96 (wider CIs accepted)
R_SMOKE = 300            # reduced-R smoke pass before full R


def m_grid_for(n: int) -> list[int]:
    """Geometric m grid truncated to ``m <= n/2`` plus the census ``m=n``."""
    grid = [m for m in M_GEOMETRIC if m <= n // 2]
    if n not in grid:
        grid.append(n)
    return sorted(set(grid))


@dataclass(frozen=True)
class SystemSpec:
    system: str
    observables: tuple[str, ...]
    ground_truth: str
    predicted_rho2_a: str                     # human-readable (ρ², a) prediction
    sizes: tuple[int, ...]
    params: dict[str, Any] = field(default_factory=dict)


# ---- the systems / observables table (Phase-4 sec.4.1) --------------------

SYSTEMS: dict[str, SystemSpec] = {
    "S1": SystemSpec("product_chain", ("magnetization",), "exact_marginal",
                     "(1, 0)", tuple(SIZES_EXACT + SIZES_EXACT_VALIDATORS)),
    "S2": SystemSpec("cards", ("AD",), "exact_marginal",
                     "(2, 0)", tuple(SIZES_EXACT + SIZES_EXACT_VALIDATORS)),
    "S2p": SystemSpec("cards", ("fixed_point_fraction",), "exact_marginal",
                      "(~1/n, 0) SNR-degenerate control", tuple(SIZES_EXACT)),
    "S3": SystemSpec("ising2d", ("abs_magnetization", "energy_density", "block_mean"),
                     "mcmc_cluster", "(Θ(1),0) off-crit; M:a=7/8, E:a~0 at T_c",
                     tuple(L * L for L in ISING_L)),
    "S4": SystemSpec("curie_weiss", ("magnetization",), "metropolis",
                     "T<T_c: a=1; T=T_c: a=1/2", tuple(MEANFIELD_N)),
    "S5": SystemSpec("ising1d", ("magnetization",), "transfer_matrix",
                     "a interpolates 0-> with xi/L", tuple(L for L in S5_L)),
    "S6": SystemSpec("ssep", ("window_density",), "exact_marginal",
                     "conserved -> slow modes; test a", tuple(S6_N)),
}


# ---- staged plan (design brief staging; architecture sec.7) --------------------

STAGES: dict[str, dict[str, Any]] = {
    "stage0": {
        "title": "Validation: core + S1,S2,S3 + S2' observable + analytical & Phase-3 gates",
        "systems": ["S1", "S2", "S2p", "S3"],
        "gate": "analytical + Phase-3 regression suites green",
    },
    "stage1": {
        "title": "Primary-1 & Primary-5 on a=0 systems (size-free m*)",
        "systems": ["S1", "S2", "S3"],           # S3 at high T
        "predictions": ["P3", "P5"],
    },
    "stage2": {
        "title": "Primary-2 & Primary-3: Ising T_c magnetization; energy-vs-magnetization",
        "systems": ["S3"],
        "predictions": ["P2", "P4", "P7"],
    },
    "stage3": {
        "title": "Primary-4: the 2x2 factorial (adds S2', S4)",
        "systems": ["S2", "S2p", "S4"],
        "predictions": ["P1", "P6"],
    },
    "stage4": {
        "title": "Secondary: S5, S6, information saturation, self-calibration",
        "systems": ["S5", "S6"],
        "predictions": ["P8"],
    },
}

# Systems implemented so far.  The others register onto the same interface.
#   S1 chain, S2 cards, S2' fixed-point, S3 Ising (Wolff GT + checkerboard Glauber),
#   S4 Curie-Weiss (exact mean-field sampler; a=1 hard floor measured, CX1).
IMPLEMENTED_SYSTEMS = {"S1", "S2", "S2p", "S3", "S4"}
