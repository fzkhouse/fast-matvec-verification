# Public optimization policy

The benchmark starts from the complete, readable RLT algorithm and compares it
with the exact baseline sources listed in `PROVENANCE.md`.  A baseline is never
weakened, replaced, or configured differently to improve a comparison.

If a public reference path is not competitive on a paper parameter point, the
following additions may be enabled in this order and must be documented in the
benchmark CSV metadata:

1. compiler optimisation and ordinary AVX2 vector operations;
2. preallocated per-invocation workspaces;
3. offline precomputation tables already required by the algorithm;
4. straightforward cache-friendly loop order.

The Kyber reference currently also uses two small, parameter-independent x86-64
assembly helpers for vector accumulation.  They have the same public input and
output contract as their C fallbacks, are selected at build time, and do not
skip the final NTT equality check.

The public artifact excludes private parameter-specialized assembly kernels,
fixed-parameter dispatch, packed binary layouts, fusion across verification
stages, prefetch policies, and unpublished tuning constants.  A result that
uses an excluded technique must not be reported as a public-reference result.

Small batches may legitimately be slower because online fixed costs dominate.
Every released result reports the full crossover.  A scheme that has no
public-reference advantage at a paper-scale point is labelled
`reference-readable` and is not presented as a performance claim.
