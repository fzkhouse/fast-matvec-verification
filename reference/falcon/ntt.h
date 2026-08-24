#ifndef AC26_NTT_H
#define AC26_NTT_H

#include <stdint.h>

typedef struct {
    uint32_t n;
    uint32_t q;
    uint32_t primitive_root;
    uint32_t psi;
    uint32_t psi_inv;
    uint32_t omega;
    uint32_t omega_inv;
    uint32_t n_inv;
    uint64_t qinv;            // NEW: floor(2^64 / q)
    uint32_t *psi_pows;
    uint32_t *psi_inv_pows;

    uint32_t *tw_fwd;
    uint32_t *tw_inv;
} ntt_ctx_t;

int ntt_ctx_init(ntt_ctx_t *ctx, uint32_t n, uint32_t q);
void ntt_ctx_free(ntt_ctx_t *ctx);
int ntt_forward(uint32_t *a, const ntt_ctx_t *ctx);
int ntt_forward_bitrev(uint32_t *a, const ntt_ctx_t *ctx); // new: DIF forward, bit-reversed output
int ntt_inverse(uint32_t *a, const ntt_ctx_t *ctx);
// negacyclic NTT for R_q = Z_q[x]/(x^n + 1)
int ntt_forward(uint32_t *a, const ntt_ctx_t *ctx);
int ntt_inverse(uint32_t *a, const ntt_ctx_t *ctx);

#endif