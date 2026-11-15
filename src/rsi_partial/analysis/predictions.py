"""The theory layer -- the ONLY place theoretical predictions live.

Epistemic guardrail (design brief; Phase-4 sec.8.2): estimators MEASURE raw quantities;
comparison to any predicted value happens ONLY here and in :mod:`rsi_partial.analysis.score`.
Nothing in ``systems/``, ``observables/``, ``observe/``, ``inference/``, ``measure/``
may import this module.  (A test enforces the import boundary.)

Every value carries an epistemic status tag:
  * ``EST``  -- established mathematically (proven / exactly derivable)
  * ``EMP``  -- established empirically in Phase 3 (measured, robust across tested range)
  * ``CONJ`` -- conjectured (plausible, unproven)
  * ``PCV``  -- pending computational validation (Phase 5 targets)

References: RSI paper (``singhania_mihika_keerthana_finalpaper.pdf``),
``Phase4_conjecture_and_design.md`` (frozen spec), ``Phase3_findings.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.special import ellipk
from scipy.stats import norm

# ---------------------------------------------------------------------------
# Detection / estimation constants
# ---------------------------------------------------------------------------

#: One-sided detection design used throughout Phase 5 (Phase-4 sec.3):
#: alpha = 0.05 (size, one-sided), beta = 0.20 (Type-II), power = 0.80.
ALPHA_DEFAULT = 0.05
BETA_DEFAULT = 0.20


def z(p: float) -> float:
    """Standard-normal quantile ``z_p`` (``P(Z <= z_p) = p``)."""
    return float(norm.ppf(p))


def kappa(alpha: float = ALPHA_DEFAULT, beta: float = BETA_DEFAULT) -> float:
    """Detection constant ``kappa = z_{1-alpha} + z_{1-beta}`` (one-sided).

    For ``alpha=0.05, beta=0.20``: ``kappa ~ 2.487``, ``kappa^2 ~ 6.18``.  [EST]
    """
    return z(1.0 - alpha) + z(1.0 - beta)


def kappa_sq(alpha: float = ALPHA_DEFAULT, beta: float = BETA_DEFAULT) -> float:
    k = kappa(alpha, beta)
    return k * k


def detection_mstar(delta: float, rho2: float, *, alpha=ALPHA_DEFAULT, beta=BETA_DEFAULT,
                    R_n: float = 0.0) -> float:
    """Predicted detection sample complexity.  [EST leading-order / EMP constant]

    Floor-corrected form (Phase-4 sec.1.2, sec.3 P4):
        ``m*(delta) = kappa^2 / [ rho^2 (delta^2 - kappa^2 R_n) ]``.
    With ``R_n=0`` this is the size-free ``kappa^2/(rho^2 delta^2)``.  Returns ``inf``
    when ``delta <= kappa*sqrt(R_n)`` (infeasible: no ``m`` certifies).
    """
    k2 = kappa_sq(alpha, beta)
    denom = rho2 * (delta * delta - k2 * R_n)
    if denom <= 0.0:
        return math.inf
    return k2 / denom


def estimation_mstar(epsilon: float, rho2: float, *, alpha=ALPHA_DEFAULT) -> float:
    """Predicted estimation (credible-interval) sample complexity.  [EST, RSI Eq. 5.7]

    ``m*(epsilon) = z_{1-alpha/2}^2 / (rho^2 epsilon^2)``.  For cards (``rho^2=2``)
    this is ``z_{0.975}^2/(2 eps^2) ~ 1.92/eps^2`` -- distinct from the detection
    constant (modelling note M2).
    """
    zc = z(1.0 - alpha / 2.0)
    return (zc * zc) / (rho2 * epsilon * epsilon)


# ---------------------------------------------------------------------------
# S1 product chain -- EXACT [EST]
# ---------------------------------------------------------------------------

class Chain:
    O_infinity = 0.5                    # [EST]
    rho2_eq = 1.0                       # [EST]
    a = 0                               # [EST] (independent coordinates)

    @staticmethod
    def rho2(phi: float) -> float:      # [EST] exact
        return phi / (2.0 - phi)

    @staticmethod
    def var_phihat(phi: float, m: int) -> float:  # [EST] exact
        return phi * (2.0 - phi) / m

    @staticmethod
    def detection_const() -> float:     # kappa^2 / rho^2_eq = kappa^2  (~6.18)
        return kappa_sq() / Chain.rho2_eq


# ---------------------------------------------------------------------------
# S2 cards (random transposition) -- mixed EST / EMP
# ---------------------------------------------------------------------------

class Cards:
    rho2_eq = 2.0                       # [EST] RSI sec.6  (E[D]^2/Var(D) -> 2)
    a = 0                               # [EMP->CONJ] via R_n ~ 1/n

    @staticmethod
    def AD_infinity(n: int) -> float:   # [EST] Prop. 2.6 / App. A
        return (n * n - 1) / (3.0 * n)

    @staticmethod
    def phi_of_t(t: int, n: int) -> float:  # [EST] exact mean law / single-card kernel
        return 1.0 - (1.0 - 2.0 / n) ** t

    @staticmethod
    def rho2(phi: float) -> float:      # [EMP] leading-order off-equilibrium SNR law
        """``rho^2(phi) = phi / (3/2 - phi)`` -- Phase-3 empirical, leading order in n.

        NOT exact (only ``rho^2(1)=2`` is [EST]).  Lives here, never returned as truth
        by the system/observable layer.
        """
        return phi / (1.5 - phi)

    @staticmethod
    def var_phihat(phi: float, m: int) -> float:  # [EMP] leading order
        return phi * (3.0 - 2.0 * phi) / (2.0 * m)

    @staticmethod
    def detection_const() -> float:     # [EMP] kappa^2 / rho^2_eq  (~3.09)
        return kappa_sq() / Cards.rho2_eq

    @staticmethod
    def R_n(n: int) -> float:           # [EMP->CONJ] cards floor ~ 1/n (order only)
        """Strongly-suggested scaling ``R_n ~ c/n``; the constant is NOT pinned.

        Provided for scoring the *slope* (a-1 = -1), never as a point target.
        """
        return 1.0 / n


# Exact per-card equilibrium moments (used by the exact-variance analytical test).
def cards_var_D(n: int) -> float:
    """Equilibrium ``Var(D)`` for a single card, ``D=|sigma(i)-i|`` with i uniform. [EST]

    ``E[D]=AD_inf``, ``E[D^2] = (n^2-1)/6`` (RSI App. E.2), so
    ``Var(D) = (n^2-1)/6 - AD_inf^2``.
    """
    ad = Cards.AD_infinity(n)
    return (n * n - 1) / 6.0 - ad * ad


def cards_exact_var_Ohat(n: int, m: int, EVx: float, floor: float) -> float:
    """Exact ``Var(Ohat_m)`` per RSI App. E.1/E.3 (identity star), given the *measured*
    ``E_sigma[V_sigma]`` and floor ``Var_sigma(AD)``:

        ``Var(Ohat_m) = (n-m)/(m(n-1)) E_sigma[V_sigma] + Var_sigma(AD)``.  [EST]
    """
    fpc_term = (n - m) / (m * (n - 1)) * EVx
    return fpc_term + floor


# ---------------------------------------------------------------------------
# S3 2-D Ising -- Onsager exact equilibrium references [EST]
# ---------------------------------------------------------------------------

T_C_ISING = 2.0 / math.log(1.0 + math.sqrt(2.0))   # ~2.269185


def ising_onsager_energy(T: float) -> float:
    """Exact Onsager internal energy per site ``u(T)`` for the 2-D square Ising
    model (``J=1``, zero field), periodic infinite lattice.  [EST]

    ``u = -coth(2b)[1 + (2/pi)(2 tanh^2(2b) - 1) K(kappa^2)]``,
    ``b=1/T``, modulus ``kappa = 2 sinh(2b)/cosh^2(2b)``, ``K`` complete elliptic
    integral of the first kind (scipy parameter convention: ``ellipk(kappa^2)``).
    """
    b = 1.0 / T
    s2 = math.sinh(2.0 * b)
    c2 = math.cosh(2.0 * b)
    t2 = math.tanh(2.0 * b)
    kappa_m = 2.0 * s2 / (c2 * c2)          # modulus k
    coeff = 2.0 * t2 * t2 - 1.0             # vanishes exactly at T_c (kappa_m -> 1)
    K = float(ellipk(kappa_m * kappa_m))    # scipy takes m = k^2; K -> inf at T_c
    coth2 = c2 / s2
    # Removable 0*inf singularity at T_c: coeff -> 0 linearly while K -> inf only
    # logarithmically, so coeff*K -> 0.  Guard against the numerical nan.
    term = 0.0 if not math.isfinite(K) else coeff * K
    return -coth2 * (1.0 + (2.0 / math.pi) * term)


def ising_onsager_magnetization(T: float) -> float:
    """Exact spontaneous magnetization ``m(T)`` (Yang) for ``T < T_c``; ``0`` for
    ``T >= T_c``.  [EST]

    ``m = (1 - sinh(2b)^{-4})^{1/8}``, ``b=1/T``.
    """
    if T >= T_C_ISING:
        return 0.0
    b = 1.0 / T
    s = math.sinh(2.0 * b)
    val = 1.0 - s ** (-4)
    return val ** 0.125 if val > 0.0 else 0.0


def curie_weiss_m0(T: float, J: float = 1.0) -> float:
    """Spontaneous magnetization ``m0(T)`` of the Curie-Weiss model: the positive root of
    ``m = tanh(J m / T)`` for ``T < T_c = J``; ``0`` for ``T >= J``.  [EST]"""
    if T >= J:
        return 0.0
    m = 0.9
    for _ in range(200):                       # fixed-point iteration converges below T_c
        m = math.tanh(J * m / T)
    return m


class CurieWeiss:
    T_c = 1.0                                  # J = 1
    # Magnetization-susceptibility exponent a (chi_N[M] ~ N^a). [CONJ/PCV -- CX1 targets]
    a_below_Tc = 1.0                           # extensive: Var_nu(M) -> m0^2 = Theta(1)
    a_at_Tc = 0.5                              # critical mean-field fluctuations
    a_above_Tc = 0.0                           # self-averaging: Var_nu(M) ~ 1/N

    @staticmethod
    def m0(T: float) -> float:
        return curie_weiss_m0(T)


class Ising:
    T_c = T_C_ISING
    # Observable-susceptibility exponents at T_c (2-D Ising universality). [CONJ/PCV]
    a_magnetization_Tc = 7.0 / 8.0     # gamma/nu = 7/4 over d=2 -> R_n ~ n^{-1/16}
    a_energy_Tc = 0.0                  # log-divergent specific heat -> a ~ 0 (+log)
    a_offcritical = 0.0                # finite xi -> summable


# ---------------------------------------------------------------------------
# Falsifiable predictions P1-P8 (Phase-4 sec.3) -- support / falsify bands
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Prediction:
    """A single falsifiable prediction with pre-registered support/falsify bands.

    ``support``/``falsify`` are closed intervals on the *measured quantity* named by
    ``quantity``.  ``scorer`` (optional) overrides the default CI-in-band logic for
    predictions that are not simple interval membership.
    """

    pid: str
    status: str
    statement: str
    quantity: str
    support: tuple[float, float] | None = None
    falsify: tuple[float, float] | None = None
    notes: str = ""
    scorer: Callable | None = field(default=None, repr=False)


PREDICTIONS: dict[str, Prediction] = {
    "P1": Prediction(
        "P1", "CONJ",
        "Large-m variance floor equals whole-system fluctuation Var_nu(O).",
        quantity="plateau_over_Var_nu_O_ratio",
        support=(0.90, 1.10),           # plateau within +/-10% of Var_nu(O)
        falsify=(-math.inf, 0.90),      # plateau below Var_nu(O) breaks identity
        notes="Measure Var(Ohat_m) plateau and independent Var_nu(O); ratio -> 1.",
    ),
    "P2": Prediction(
        "P2", "PCV",
        "Susceptibility exponent controls the floor: slope of log R_n vs log n = a-1.",
        quantity="floor_slope",
        # support band is per-observable (a-1 +/- 0.1); default is set by the scorer
        notes="Ising |M| at T_c: a=7/8 -> slope -1/8; energy: a~0 -> slope ~ -1.",
    ),
    "P3": Prediction(
        "P3", "EMP",
        "m* is size-free for a=0 systems: m* ~ n^b with |b|<0.1; m* delta^2 -> kappa^2/rho^2.",
        quantity="mstar_exponent_b",
        support=(-0.10, 0.10),
        falsify=(0.25, math.inf),
        notes="CI must exclude 0.25 to support; b>0.25 falsifies for an a=0 system.",
    ),
    "P4": Prediction(
        "P4", "CONJ",
        "m* breaks when a>0: infeasibility boundary tracks delta = kappa*sqrt(R_n).",
        quantity="infeasibility_boundary_ratio",
        support=(0.85, 1.15),           # observed boundary / (kappa sqrt(R_n)) within +/-15%
        falsify=None,
        notes="Falsify if m* at T_c is flat in n (no floor effect).",
    ),
    "P5": Prediction(
        "P5", "CONJ",
        "Size-independence is a property of (rho^2,a), not the system.",
        quantity="mstar_delta2_agreement",
        notes="m*delta^2 agrees across a=0 systems; M vs E at T_c differ. Composite score.",
    ),
    "P6": Prediction(
        "P6", "PCV",
        "Correlation range, not exchangeability, is operative.",
        quantity="composite",
        notes="Ising high-T (a=0) scale-free; Curie-Weiss T<T_c (a=1) size-limited.",
    ),
    "P7": Prediction(
        "P7", "PCV",
        "Ising temperature sweep traces a continuous crossover of size-dependence.",
        quantity="composite",
        notes="m*(delta,n) flat for T>=1.4 T_c, bending up as T->T_c, tracking xi/L.",
    ),
    "P8": Prediction(
        "P8", "CONJ",
        "Effective information saturates at I_max = 1/R_n; crossover m_x = 1/(rho^2 R_n).",
        quantity="Imax_over_inv_Rn_ratio",
        support=(0.85, 1.15),           # fitted I_max matches 1/R_n within +/-15%
        notes="I_eff(m) = rho^2 m/(1+rho^2 m R_n); plateau = 1/R_n (from P1).",
    ),
}


def floor_slope_band(a: float, tol: float = 0.10) -> tuple[float, float]:
    """Support band for P2: slope = a-1 within +/-tol."""
    return (a - 1.0 - tol, a - 1.0 + tol)
