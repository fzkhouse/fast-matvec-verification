#define _POSIX_C_SOURCE 200112L
#include "table.h"
#include <stdlib.h>
#include <string.h>

static void* aligned_malloc(size_t size) {
    void *p = NULL;
    if (posix_memalign(&p, 32, size) != 0) return NULL;
    return p;
}

static uint32_t rng32_next(uint32_t *s) {
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}

static void poly_rand(poly_t *p, uint32_t *s) {
    for (int i = 0; i < DILITHIUM_N; i++) p->coeffs[i] = (int32_t)(rng32_next(s) % DILITHIUM_Q);
}

int table_generate(const char *path, uint32_t R, uint32_t T, const poly_t *A) {
    table_t tb = {0};
    tb.q = DILITHIUM_Q; tb.n = DILITHIUM_N;
    tb.R = R; tb.T = T; tb.K_ = DILITHIUM_K; tb.L_ = DILITHIUM_L;

    size_t Acount = (size_t)R * tb.K_ * tb.L_;
    size_t Gcount = (size_t)R * tb.K_;
    size_t ATcount = (size_t)R * tb.L_;

    tb.randoms = (uint32_t*)malloc((size_t)R * T * sizeof(uint32_t));
    tb.A = (poly_t*)aligned_malloc(Acount * sizeof(poly_t));
    tb.gamma = (poly_t*)aligned_malloc(Gcount * sizeof(poly_t));     /* NTT 域 */
    tb.ATgamma = (poly_t*)aligned_malloc(ATcount * sizeof(poly_t));  /* NTT 域 */
    if (!tb.randoms || !tb.A || !tb.gamma || !tb.ATgamma) return -1;

    for (uint32_t r = 0; r < R; r++) {
        memcpy(&tb.A[r * tb.K_ * tb.L_], A, tb.K_ * tb.L_ * sizeof(poly_t));
    }

    /* 预计算 A 的 NTT（一次） */
    poly_t A_ntt[DILITHIUM_K * DILITHIUM_L];
    for (uint32_t i = 0; i < tb.K_ * tb.L_; i++) {
        A_ntt[i] = A[i];
        poly_ntt(&A_ntt[i]);
    }

    uint32_t s = 0xC0FFEEu;
    for (uint32_t i = 0; i < R * T; i++) tb.randoms[i] = rng32_next(&s) % DILITHIUM_Q;

    for (uint32_t i = 0; i < Gcount; i++) {
        poly_rand(&tb.gamma[i], &s);
        poly_ntt(&tb.gamma[i]);
    }

    for (uint32_t r = 0; r < R; r++) {
        for (uint32_t l = 0; l < tb.L_; l++) {
            poly_zero(&tb.ATgamma[r * tb.L_ + l]);
            for (uint32_t k = 0; k < tb.K_; k++) {
                const poly_t *Akl = &A_ntt[k * tb.L_ + l];
                const poly_t *gk  = &tb.gamma[r * tb.K_ + k];
                poly_pointwise_accum(&tb.ATgamma[r * tb.L_ + l], Akl, gk);
            }
        }
    }

    int rc = table_save(path, &tb);
    table_free(&tb);
    return rc;
}
