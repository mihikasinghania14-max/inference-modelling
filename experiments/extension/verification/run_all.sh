#!/usr/bin/env bash
# Re-run every verification script and compare with the stored outputs in
# results/outputs/extension-verification/.   Usage:  bash run_all.sh [outdir]
# Requires: python3 with numpy, scipy, sympy, mpmath; a C compiler (cc/gcc) for wolff.c.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
ref="$here/../../../results/outputs/extension-verification"
out="${1:-$(mktemp -d)}"
mkdir -p "$out"
cd "$here"
status=0
for s in shuffle_algebra cw ising_exact shuffle_sim field; do
  python3 -I "$s.py" > "$out/$s.txt"
  if diff -q "$out/$s.txt" "$ref/$s.txt" >/dev/null; then echo "OK   $s.py identical to stored output"
  else echo "DIFF $s.py differs from stored output (see $out/$s.txt)"; status=1; fi
done
cc -O2 -o "$out/wolff" wolff.c -lm
for L in 16 32; do
  "$out/wolff" "$L" 200000 1 > "$out/wolff_L$L.dat" 2> "$out/wolff_L$L.err"
  if cmp -s "$out/wolff_L$L.dat" "$ref/wolff_L$L.dat"; then echo "OK   wolff L=$L seed 1 identical to stored output"
  else echo "DIFF wolff L=$L differs from stored output"; status=1; fi
done
python3 wolff_analyze.py "$out"/wolff_L16.dat "$out"/wolff_L32.dat | tee "$out/wolff_summary.txt"
exit $status
