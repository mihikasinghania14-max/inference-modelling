"""
Average Displacement at Theoretical Mixing Time: Simulation vs Theory
Reads data from plot.xlsx and produces a publication-quality figure.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Load data ─────────────────────────────────────────────────────────────────
df_raw = pd.read_excel('/mnt/user-data/uploads/plot.xlsx', header=0)

# Split into powers-of-2 block (rows 0–12, skip row 0 which is n=0 placeholder)
# and intermediate block (rows 15–22)
pow2 = df_raw.iloc[1:13].copy()
inter = df_raw.iloc[15:23].copy()

# Rename columns for convenience
cols = ['n', 'theory', 'mean_sim', 'sd']
pow2.columns = cols
inter.columns = cols

# Cast to numeric
for df in (pow2, inter):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')

pow2  = pow2.dropna()
inter = inter.dropna()

# Combined sorted set for the theoretical curve (all n values)
all_n = np.sort(np.concatenate([pow2['n'].values, inter['n'].values]))

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))

# 1. Theoretical curve — thin black line through all n values
all_theory = (all_n**2 - 1) / (3 * all_n)
ax.plot(
    all_n, all_theory,
    color='black', linewidth=1.6, linestyle='-',
    zorder=2, label=r'Theoretical: $AD(\infty) = \frac{n^2-1}{3n}$',
)

# 2. Powers-of-2 simulation — blue points connected by a line
ax.errorbar(
    pow2['n'], pow2['mean_sim'],
    yerr=pow2['sd'],
    fmt='o-',
    color='#2563EB',
    ecolor='#93C5FD',
    elinewidth=1.4, capsize=4, capthick=1.4,
    markersize=6, linewidth=1.4,
    zorder=4, label='Simulation — powers of 2',
)

# 3. Intermediate n values — orange points on the same line
ax.errorbar(
    inter['n'], inter['mean_sim'],
    yerr=inter['sd'],
    fmt='s',                        # square markers, no connecting line
    color='#F97316',
    ecolor='#FED7AA',
    elinewidth=1.4, capsize=4, capthick=1.4,
    markersize=6,
    zorder=5, label='Simulation — intermediate $n$ values',
)

# ── Axes formatting ───────────────────────────────────────────────────────────
ax.set_xscale('log', base=2)

# Force the exact n values we care about to appear as x-ticks
pow2_ticks = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
ax.set_xticks(pow2_ticks)
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.tick_params(axis='x', rotation=45)

ax.set_xlabel('Number of Cards ($n$)', fontsize=13)
ax.set_ylabel('Average Displacement (AD)', fontsize=13)
ax.set_title(
    'Average Displacement at Theoretical Mixing Time: Simulation vs Theory',
    fontsize=13, pad=12,
)

ax.legend(fontsize=10, loc='upper left')
ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
ax.set_xlim(1.5, 5500)

fig.tight_layout()
fig.savefig('/mnt/user-data/outputs/ad_simulation_vs_theory.png', dpi=200, bbox_inches='tight')
print("Saved.")
