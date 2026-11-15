#!/usr/bin/env bash
# Critical 2D Ising campaign with the Wolff sampler (ising_wolff.c).
# Two independent starts (random, ordered) per size; block means are written to $OUT.
#   bash run_ising_campaign.sh [outdir] [jobs]
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here/../../../data/extension/regenerated/ising_critical_blocks}"
jobs="${2:-2}"
mkdir -p "$out"
cc -O2 -o "$out/ising_wolff" "$here/ising_wolff.c" -lm
# size  measurements per run
cat > "$out/plan.txt" <<PLAN
8 400000
12 400000
16 400000
24 400000
32 400000
48 200000
64 150000
96 80000
128 60000
PLAN
run_one() { # L nmeas start seed
  "$1/ising_wolff" "$2" "$3" "$5" "$4" > "$1/L$2_start$4.dat" 2> "$1/L$2_start$4.err"
}
export -f run_one
while read -r L N; do
  echo "$out $L $N 0 $((1000 + L))"
  echo "$out $L $N 1 $((2000 + L))"
done < "$out/plan.txt" | xargs -P "$jobs" -L 1 bash -c 'run_one "$@"' _
rm -f "$out/ising_wolff"
echo "done: $(ls "$out"/*.dat | wc -l) block files in $out"
