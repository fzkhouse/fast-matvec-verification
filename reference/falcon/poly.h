#ifndef AC26_POLY_H
#define AC26_POLY_H

#include <stdint.h>
#include "ntt.h"

typedef struct {
    uint32_t n;
    uint32_t *c;
} poly_t;

int poly_init(poly_t *p, uint32_t n);
void poly_free(poly_t *p);
void poly_zero(poly_t *p);
void poly_copy(poly_t *dst, const poly_t *src);

void poly_rand(poly_t *p, uint32_t q);
void poly_add_inplace(poly_t *a, const poly_t *b, uint32_t q);
void poly_sub_inplace(poly_t *a, const poly_t *b, uint32_t q);
void poly_add_scaled_inplace(poly_t *acc, const poly_t *x, uint32_t s, uint32_t q);

int poly_ntt_inplace(poly_t *p, const ntt_ctx_t *ctx);
int poly_intt_inplace(poly_t *p, const ntt_ctx_t *ctx);

// acc += a (*) b  (a,b 在 NTT 域；Hadamard)
void poly_hadamard_accumulate(poly_t *acc, const poly_t *a, const poly_t *b, uint32_t q);

int poly_equal(const poly_t *a, const poly_t *b);

#endif