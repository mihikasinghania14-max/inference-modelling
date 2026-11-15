#!/usr/bin/env bash
# Regenerate the stored Stage-1 "cards" smoke results and compare them with the copy in
# results/outputs/extension-cards-smoke.  Takes a few minutes on a laptop.
#   bash scripts/reproduce_extension_cards_smoke.sh [outdir]
set -euo pipefail
cd "$(dirname "$0")/.."
out="${1:-$(mktemp -d)}"
python -m rsi_partial.experiments.run snr_floor --system cards --sizes 512 2048 8192 --smoke --out "$out"
python -m rsi_partial.experiments.run detection --system cards --sizes 2048 8192 --deltas 0.15 0.10 --smoke --out "$out"
python scripts/compare_with_stored.py results/outputs/extension-cards-smoke "$out"
