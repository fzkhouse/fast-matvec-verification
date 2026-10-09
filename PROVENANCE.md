# Source and comparison provenance

The public artifact is derived only from the AC26 experiment tree supplied with
this repository.  The following mapping is authoritative for benchmark
comparisons.

| Scheme | Public RLT source | Baseline retained in this release | Original benchmark role |
| --- | --- | --- | --- |
| Falcon-512 | `reference/falcon` | `third_party/falcon512avx2/bench_matvec_falcon_native.c` | native matrix-vector equality baseline |
| Kyber-512 (K=5,L=9) | `reference/kyber` | `baseline/kyber/bench_kyber_native_mv.c` | native matrix-vector baseline |
| Dilithium-5 | `reference/dilithium` | `baseline/dilithium/main_baseline.c` | native row-by-row NTT baseline |
| RPT | `reference/rpt/rpt_verify.c` (`RPT_REFERENCE=1`) | same source with `RPT_REFERENCE=0` | original direct `B*z-u` comparison branch |

The public Kyber reference uses the last complete-C verification implementation
from AC26 history (`05dd808`) as the starting point. It restores the final
NTT-domain equality check; the private binary fast path that omitted that check
is deliberately not published. The release includes two small,
parameter-independent x86-64 assembly helpers for vector accumulation, with C
fallbacks selected by `USE_KYBER_ASM=0`. It excludes AC26's private
parameter-specialized assembly kernels, special K=5/L=9 dispatch, packed
binary input path, prefetch policies, and unpublished tuning constants.

The public Falcon and Dilithium reference trees were copied from the original
experiment paths and retain their ordinary C/AVX2 computation paths. Their
baseline sources are not simplified or substituted. The exact third-party
source-package locations and licenses are recorded in `third_party/README.md`.

The historical Dilithium files under the original experiment data directory
must not be used as a single benchmark fixture: their `A`, `z`, `u`, and table
artifacts were produced at different times and are not a coherent input set.
The public `dilithium5-gen` target generates all four files in one run.

The public RPT implementation retains the original experiment's two
arithmetic paths: the RPT statistic `u*C_i-z*rho_i` and the direct baseline
`B*z-u`. It replaces only the original verbose file-reading/printing code with
a checked contiguous loader; it does not use sparse, packed, or
parameter-specialized kernels.
