# batch-vrfy-open

Portable, reproducible reference experiments for batched verification of
lattice-style matrix-vector relations. This is an experimental research
artifact, not a production cryptographic library.

## Quick start

```sh
make
./run.sh --scheme falcon512 --impl both --rows 1000,3000,5000,7000 --data-root /path/to/falcon-fixtures
./run.sh --scheme kyber512 --impl both --rows 10000,30000,50000,100000 --data-root /path/to/kyber-fixtures
./run.sh --scheme rpt --impl both --rows 512,2048,4096,8192 --data-root /path/to/rpt-fixtures
make test
```

Results are written below `results/`; generated fixtures are ignored by git.
Large paper fixtures are intentionally not committed, so every benchmark
command takes an explicit `--data-root` rather than embedding a private path.

## Implementations

`reference` performs the RLT or RPT batch check; `baseline` is the direct native
row-by-row implementation actually used by the original experiment scripts.
The source mapping is recorded in `PROVENANCE.md`.  The public reference path
is intentionally limited by `OPTIMIZATION_POLICY.md`; performance crossover
points are reported rather than hidden.

The retained `paper-*-9runs.csv` files report the complete nine-run results.
The Kyber reference retains its final NTT-domain equality check and is faster
than the retained native baseline on the released paper grid.

The challenge is a research benchmark mechanism, not a complete security
proof or a constant-time cryptographic implementation. Do not use this code as
a replacement for an audited Falcon, Kyber, or Dilithium implementation.

## Parameters and data

The paper comparison grid is in `configs/paper.toml`. The imported sources
cover Falcon, Kyber, Dilithium, and RPT. `build/dilithium5-gen` and
`build/dilithium2-gen` create coherent `A/z/u/table` fixtures; do not combine
historical files produced in separate generator runs. RPT uses the original
four files `mat_B_N.txt`, `vec_zu_N.txt`, `mat_C_N.txt`, and `mat_ro_N.txt`.

## Benchmark output

```sh
SCHEME=rpt ROWS=512,2048,4096,8192 DATA_ROOT=/path/to/rpt-fixtures RUNS=9 bench/benchmark.sh
```

The benchmark records pass counts and median/mean online verification time.
Preprocessing/data-generation time is kept outside the online timing window.
Some portable paths report zero cycles until a platform-specific counter is
added; time statistics remain available on every supported platform.

`results/*-smoke.csv` contains one-run developer smoke results only.  These
files demonstrate matching input and successful verification; they are not the
paper's final nine-run statistics and must not be used to redraw paper tables.

## Licensing

Original project code is Apache-2.0. See `NOTICE` before adding any third-party
source. Third-party Falcon, Kyber, and Dilithium implementations must retain
their own license and attribution files and must not be represented as this
project's code.
