#!/usr/bin/env bash
# Generate QE powder-Raman inputs for all bundled example structures.
# Usage: bash run_examples.sh [output_dir]
set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$SKILL_DIR/examples/out}"
PY="${PYTHON:-python3}"

mkdir -p "$OUT"

run() {
  local cif="$1" prefix="$2"
  echo "==> $prefix ($cif)"
  "$PY" "$SKILL_DIR/scripts/qe_raman.py" prepare "$SKILL_DIR/examples/$cif" \
      --outdir "$OUT/$prefix" --prefix "$prefix"
}

run "LiF_mp-1138.cif"     "LiF"
run "MgH2_mp-23710.cif"   "MgH2"
run "Al2O3_mp-1143.cif"   "Al2O3"
run "AlAs_mp-2172.cif"    "AlAs"
run "CdS_mp-672.cif"      "CdS"
run "MoS2_mp-2815.cif"    "MoS2"

echo "Done. Inputs written under $OUT"
