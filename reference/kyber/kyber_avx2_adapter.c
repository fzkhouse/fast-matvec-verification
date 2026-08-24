#include "kyber_avx2_adapter.h"
#include <string.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"

#ifndef KYB_N
#define KYB_N KYBER_N
#endif

#ifndef KYB_Q
#define KYB_Q KYBER_Q
#endif

#if defined(__GNUC__) || defined(__clang__)
#define ALIGN32 __attribute__((aligned(32)))
#else
#define ALIGN32
#endif

static inline int16_t norm_q_i16(int32_t x) {
    int32_t v = x % KYB_Q;
    if (v < 0) v += KYB_Q;
    return (int16_t)v;
}

void kyb_poly_zero(poly *p) {
    memset(p, 0, sizeof(*p));
}

void kyb_poly_from_u32(poly *dst, const uint32_t *src) {
    for (int i = 0; i < KYB_N; i++) {
        dst->coeffs[i] = norm_q_i16((int32_t)(src[i] % KYB_Q));
    }
}

void kyb_poly_from_u64_modq(poly *dst, const uint64_t *src) {
    for (int i = 0; i < KYB_N; i++) {
        dst->coeffs[i] = norm_q_i16((int32_t)(src[i] % KYB_Q));
    }
}

void kyb_poly_to_u32(uint32_t *dst, const poly *src) {
    for (int i = 0; i < KYB_N; i++) {
        int32_t v = src->coeffs[i] % KYB_Q;
        if (v < 0) v += KYB_Q;
        dst[i] = (uint32_t)v;
    }
}

void kyb_poly_ntt(poly *p) {
    poly_ntt(p);
    poly_reduce(p);
}

void kyb_poly_reduce(poly *p) {
    poly_reduce(p);
}

void kyb_poly_basemul_acc(poly *acc, const poly *a_ntt, const poly *b_ntt) {
    ALIGN32 poly t;
    poly_basemul_montgomery(&t, a_ntt, b_ntt);
    poly_add(acc, acc, &t);
    poly_reduce(acc); // 关键：恢复每次累加后约减，保证正确性
}

int kyb_poly_equal_modq(const poly *a, const poly *b) {
    for (int i = 0; i < KYB_N; i++) {
        int32_t x = a->coeffs[i] % KYB_Q;
        int32_t y = b->coeffs[i] % KYB_Q;
        if (x < 0) x += KYB_Q;
        if (y < 0) y += KYB_Q;
        if (x != y) return 0;
    }
    return 1;
}

int kyb_poly_equal_reduced(const poly *a, const poly *b) {
    return memcmp(a->coeffs, b->coeffs, sizeof(a->coeffs)) == 0;
}