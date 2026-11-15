"""Observables.  Each implements :class:`rsi_partial.core.interfaces.Observable`.

Decoupled from systems: the same estimator/detector code runs unchanged by swapping
the per-component readout.  cards -> {AD, fixed_point_fraction}; chain -> {magnetization};
ising -> {abs_magnetization, energy_density, block_mean}.
"""
