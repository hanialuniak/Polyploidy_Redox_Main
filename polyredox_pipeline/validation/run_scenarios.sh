#!/usr/bin/env bash
# Simulate each scenario, run the full pipeline, and print the verdicts.
# Usage: bash validation/run_scenarios.sh OUTPUT_DIR
set -euo pipefail
out=${1:?output directory}
export PYTHONPATH="$(cd "$(dirname "$0")/.." && pwd)"
for scenario in effect null general; do
  python3 -m polyredox simulate --outdir "$out/$scenario" --scenario "$scenario" --seed 11
  python3 -m polyredox run "$out/$scenario/config.yaml" 2> "$out/$scenario/run.log"
  echo "=== $scenario"
  python3 - "$out/$scenario/results" <<'P'
import json, sys, pandas as pd
r = sys.argv[1]
v = pd.read_csv(f"{r}/redox_verdicts.tsv", sep="\t")
print(v[["programme", "beta", "ci_low", "ci_high", "p_adj", "verdict"]].round(3).to_string(index=False))
m = json.load(open(f"{r}/ml_summary.json"))
print({k: round(m[k], 3) for k in ("ap_technical", "ap_technical_plus_redox", "ap_gain", "null_gain_q95", "p_empirical")})
print(pd.read_csv(f"{r}/doublet_excess.tsv", sep="\t")[["cell_class", "observed_polyploid_rate", "expected_doublet_rate", "excess_ci_low"]].round(3).to_string(index=False))
P
done
