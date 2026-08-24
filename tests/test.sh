#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
make -C "$ROOT" all >/dev/null

# Falcon and Kyber paper fixtures are intentionally external to the release;
# make sure every published executable is present, then use generated,
# self-contained Dilithium fixtures to check both supported parameter modes.
for bin in falcon-reference falcon-baseline kyber-reference kyber-baseline \
           dilithium-reference dilithium-baseline dilithium5-gen \
           dilithium2-reference dilithium2-baseline dilithium2-gen; do
  [[ -x "$ROOT/build/$bin" ]]
done

TMPDIR_TEST="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_TEST"' EXIT
for mode in 2 5; do
  scheme="dilithium${mode}"
  if [[ "$mode" == 5 ]]; then
    ref="$ROOT/build/dilithium-reference"
    base="$ROOT/build/dilithium-baseline"
    gen="$ROOT/build/dilithium5-gen"
  else
    ref="$ROOT/build/dilithium2-reference"
    base="$ROOT/build/dilithium2-baseline"
    gen="$ROOT/build/dilithium2-gen"
  fi
  mkdir -p "$TMPDIR_TEST/$scheme"
  "$gen" 32 19 256 "$TMPDIR_TEST/$scheme" >/dev/null
  suffix=0032
  "$ref" "$TMPDIR_TEST/$scheme/mat_z_${suffix}.txt" \
      "$TMPDIR_TEST/$scheme/mat_u_${suffix}.txt" "$TMPDIR_TEST/$scheme/table_${suffix}.txt" | grep -q 'result=PASS'
  "$base" "$TMPDIR_TEST/$scheme/mat_z_${suffix}.txt" \
      "$TMPDIR_TEST/$scheme/mat_u_${suffix}.txt" "$TMPDIR_TEST/$scheme/mat_A_${suffix}.txt" \
      "$TMPDIR_TEST/$scheme/table_${suffix}.txt" | grep -q 'result=PASS'
done

"$ROOT/run.sh" --help >/dev/null
echo "published binaries and Dilithium-2/-5 fixtures passed"
