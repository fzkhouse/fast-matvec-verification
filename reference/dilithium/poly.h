#ifndef DILITHIUM_BRLT_POLY_H
#define DILITHIUM_BRLT_POLY_H

#include "params.h"
#include <stdint.h>

typedef struct {
    __attribute__((aligned(32)))
    int32_t coeffs[DILITHIUM_N];
} poly_t;

void poly_zero(poly_t *a);
void poly_add_to(poly_t *acc, const poly_t *a);
void poly_mul_scalar_add(poly_t *acc, const poly_t *a, uint32_t scalar);

void poly_negacyclic_mul(const poly_t *a, const poly_t *b, poly_t *out);
int64_t poly_dot_int64(const poly_t *a, const poly_t *b);

void poly_adjoint_mul_accum(const poly_t *a, const poly_t *g, poly_t *acc);

/* NTT helpers (Dilithium AVX2) */
void poly_ntt(poly_t *a);
void poly_invntt_tomont(poly_t *a);
void poly_pointwise_accum(poly_t *acc, const poly_t *a_ntt, const poly_t *b_ntt);
int  poly_equal_modq(const poly_t *a, const poly_t *b);

/* Montgomery -> standard coeffs */
void poly_from_mont(poly_t *a);

#endif