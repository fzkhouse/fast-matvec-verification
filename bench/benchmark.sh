#!/usr/bin/env bash
# Compatibility wrapper for the unified benchmark entry point.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

[[ -n "${SCHEME:-}" && -n "${ROWS:-}" && -n "${DATA_ROOT:-}" ]] || {
  echo "set SCHEME, ROWS, and DATA_ROOT; e.g. SCHEME=kyber512 ROWS=10000 DATA_ROOT=data/kyber $0" >&2
  exit 2
}

"$ROOT/run.sh" --scheme "$SCHEME" --rows "$ROWS" --data-root "$DATA_ROOT" \
  --runs "${RUNS:-9}" --impl "${IMPL:-both}" --out "${OUT:-$ROOT/results/${SCHEME}-benchmark.csv}"
