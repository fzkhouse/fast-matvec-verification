# batch-vrfy-open

Source artifact for the ASIACRYPT 2026 paper on batched verification of
lattice-style matrix-vector relations. This repository contains experimental
research software, not a production cryptographic library.

## Artifact badge scope

This release is submitted for the **Artifacts Available** badge only. It makes
the source code associated with the paper publicly retrievable. It is not a
submission for the Artifacts Functional or Artifacts Reproduced badges, and it
does not claim that a reviewer can reproduce every paper result from this
source archive alone.

Large experiment fixtures and locally generated benchmark CSV files are not
part of the source release. They are ignored by Git. Commands that use those
fixtures therefore require an explicit `--data-root` supplied by the user.

## Contents and paper correspondence

- `reference/` contains the released RLT and RPT reference implementations.
- `baseline/` contains experiment-specific direct verification baselines.
- `third_party/` contains the upstream Falcon, Kyber, and Dilithium components
  linked by the experiment programs.
- `bench/` and `run.sh` contain benchmark orchestration and CSV summarization.
- `configs/paper.toml` records the grid associated with Table 4 of the paper;
  its dimension sweep also records the 1024-dimensional settings associated
  with Table 2.
- `tests/` contains the release smoke test.

The detailed mapping from released sources to the original experiment tree is
in `PROVENANCE.md`. `OPTIMIZATION_POLICY.md` describes which optimizations are
included in the public reference implementation.

## Platform and dependencies

The paper experiments were developed and run on the following platform:

- CPU: Intel Xeon Platinum 8163 at 2.50 GHz;
- memory: 64 GB;
- operating system: Ubuntu 18.04 LTS, x86-64;
- compiler: GCC 7.3;
- build tool: GNU Make;
- scripts: Bash and Python 3.

The released build requires Linux on x86-64. AVX2 is required by all published
experiment binaries. The Kyber targets additionally require FMA, BMI2, POPCNT,
and AES-NI CPU support. A C11 compiler, GNU Make, Bash, Python 3, the standard C
library, and `libm` are required. The source release was also clean-built and
smoke-tested with GCC 13.3, GNU Make 4.3, and Python 3.12 on Ubuntu 24.04.

## Build and smoke test

```sh
make
make test
```

A successful smoke test ends with:

```text
published binaries and Dilithium-2/-5 fixtures passed
```

The smoke test is included as a basic integrity check for the source release;
it is not a reproduction of the paper's complete experiments.

## Running with external fixtures

The unified command has the following form:

```sh
./run.sh --scheme NAME --impl reference|baseline|both \
  --rows LIST --runs 9 --data-root /path/to/fixtures --out results/output.csv
```

Supported scheme names and expected file naming are documented by:

```sh
./run.sh --help
```

For example, with separately obtained RPT fixtures:

```sh
./run.sh --scheme rpt --impl both --rows 512,2048,4096,8192 \
  --runs 9 --data-root /path/to/rpt-fixtures
```

The CSV output records the scheme, implementation, dimensions, run count,
successful verification count, and timing summary. File loading and table
decoding occur outside the reported online verification window.

## Auxiliary scripts for numerical evaluation of probability bounds

The `probability_calculation_python/` directory contains auxiliary Python
scripts for numerical evaluations used in the paper:

- `Cw distribution simulation based on Gaussian.py` estimates the distribution
  of $\|C\omega\|$ using empirical Monte Carlo simulations.
- `Abort Probability Analysis Based on Chi-Distribution over Gaussian.py`
  computes the numerical bound $B_2$ for $\mathcal{B}_{\Omega}$.

The directory also contains `norm_distribution_analysis.pdf`, which documents
the associated numerical analysis.

## Security notice

The code is a research benchmark mechanism. It is not an audited,
constant-time implementation and must not replace production Falcon, Kyber,
Dilithium, or other cryptographic libraries.

## Licensing and provenance

Original project code is released under Apache License 2.0; see `LICENSE`.
Third-party source remains under its upstream license. See `NOTICE`,
`third_party/README.md`, and the license file beside each third-party source
tree. `PROVENANCE.md` records the experiment-source mapping.
