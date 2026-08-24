#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SCHEME=""
IMPL="both"
ROWS=""
RUNS=9
DATA_ROOT=""
OUT=""

usage() {
  cat <<'EOF'
usage: ./run.sh --scheme NAME --rows LIST --data-root DIRECTORY [options]

options:
  --impl reference|baseline|both   default: both
  --runs N                         default: 9
  --out FILE                       default: results/NAME.csv

schemes:
  falcon512    paper q=12289,n=512,K=L=1
  falcon1024   paper flattened dimension 1024: q=12289,n=512,K=L=2
  kyber512     paper q=3329,n=256,K=5,L=9
  kyber1024    q=3329,n=256,K=L=4
  dilithium2   q=8380417,n=256,K=L=4
  dilithium5   paper q=8380417,n=256,K=8,L=7
  rpt          paper dimensions: 512,2048,4096,8192

Large fixtures are intentionally not committed. --data-root identifies an
explicit generated/downloaded fixture directory; no private absolute path is
embedded in this script.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scheme) SCHEME="$2"; shift 2 ;;
    --impl) IMPL="$2"; shift 2 ;;
    --rows) ROWS="$2"; shift 2 ;;
    --runs) RUNS="$2"; shift 2 ;;
    --data-root) DATA_ROOT="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$SCHEME" && -n "$ROWS" && -n "$DATA_ROOT" ]] || { usage >&2; exit 2; }
[[ "$IMPL" == reference || "$IMPL" == baseline || "$IMPL" == both ]] || { echo "invalid --impl" >&2; exit 2; }
[[ "$RUNS" =~ ^[1-9][0-9]*$ ]] || { echo "--runs must be a positive integer" >&2; exit 2; }
[[ -n "$OUT" ]] || OUT="$ROOT/results/${SCHEME}.csv"

make -C "$ROOT" all >/dev/null
python3 "$ROOT/bench/paper_bench.py" \
  --scheme "$SCHEME" --impl "$IMPL" --rows "$ROWS" --runs "$RUNS" \
  --data-root "$DATA_ROOT" --out "$OUT"
echo "results: $OUT"
