# Third-party sources

The directories below contain source copied from the NIST post-quantum
cryptography submission packages used by the original experiments. They are
not relicensed by this project. Their license files and all notices in source
headers must be retained when the artifact is redistributed.

| Directory | Upstream package | Upstream project | License |
| --- | --- | --- | --- |
| `falcon512avx2/` | Falcon Round 3, `Optimized_Implementation/falcon512/falcon512avx2` | https://falcon-sign.info/ | MIT |
| `kyber512/` | `NIST-PQ-Submission-Kyber-20201001`, `Additional_Implementations/avx2/crypto_kem/kyber512` | https://pq-crystals.org/kyber/ | Public Domain; embedded components retain file-level notices |
| `dilithium5/` | CRYSTALS-Dilithium Round 3, `Additional_Implementations/avx2/crypto_sign/dilithium5` | https://pq-crystals.org/dilithium/ | Public Domain/CC0, Apache-2.0, or GPL-2.0; embedded components retain file-level notices |

The Falcon core is copied without algorithmic changes; the artifact adds
`bench_matvec_falcon_native.c` as its experiment-specific matrix-vector
baseline. The Kyber and Dilithium trees contain the upstream components linked
by the top-level Makefile. Project-specific reference and baseline drivers are
kept outside `third_party/`.

See the repository-level `NOTICE` and the `LICENSE` file in each directory for
the applicable terms.
