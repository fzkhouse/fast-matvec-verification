#include "poly.h"
#include "../../third_party/dilithium5/ntt.h"
#include "../../third_party/dilithium5/consts.h"

/* Dilithium q=8380417; QINV = q^(-1) mod 2^32 */
#define QINV 58728449u

static inline int32_t fqred(int64_t x) {
    x %= (int64_t)DILITHIUM_Q;
    if (x < 0) x += DILITHIUM_Q;
    return (int32_t)x;
}

static inline int32_t montgomery_reduce(int64_t a) {
    int32_t t = (int32_t)a * (int32_t)QINV;
    t = (int32_t)((a - (int64_t)t * DILITHIUM_Q) >> 32);
    return t;
}

void poly_zero(poly_t *a) {
    for (int i = 0; i < DILITHIUM_N; i++) a->coeffs[i] = 0;
}

void poly_add_to(poly_t *acc, const poly_t *a) {
    for (int i = 0; i < DILITHIUM_N; i++) {
        int64_t t = (int64_t)acc->coeffs[i] + a->coeffs[i];
        acc->coeffs[i] = fqred(t);
    }
}

void poly_mul_scalar_add(poly_t *acc, const poly_t *a, uint32_t scalar) {
    for (int i = 0; i < DILITHIUM_N; i++) {
        int64_t t = (int64_t)a->coeffs[i] * scalar;
        int64_t s = (int64_t)acc->coeffs[i] + t;
        acc->coeffs[i] = fqred(s);
    }
}

void poly_negacyclic_mul(const poly_t *a, const poly_t *b, poly_t *out) {
    poly_t ta = *a;
    poly_t tb = *b;
    poly_t tc;

    ntt_avx(ta.coeffs, qdata);
    ntt_avx(tb.coeffs, qdata);
    pointwise_avx(tc.coeffs, ta.coeffs, tb.coeffs, qdata);
    invntt_avx(tc.coeffs, qdata);

    /* invntt 输出为 Montgomery 域，转回标准系数域 */
    for (int i = 0; i < DILITHIUM_N; i++) {
        int32_t v = montgomery_reduce((int64_t)tc.coeffs[i]);
        if (v < 0) v += DILITHIUM_Q;
        out->coeffs[i] = v;
    }
}

int64_t poly_dot_int64(const poly_t *a, const poly_t *b) {
    int64_t sum = 0;
    for (int i = 0; i < DILITHIUM_N; i++) {
        sum += (int64_t)a->coeffs[i] * b->coeffs[i];
        sum %= DILITHIUM_Q;
    }
    if (sum < 0) sum += DILITHIUM_Q;
    return sum;
}

void poly_adjoint_mul_accum(const poly_t *a, const poly_t *g, poly_t *acc) {
    for (int j = 0; j < DILITHIUM_N; j++) {
        int64_t s = acc->coeffs[j];
        for (int i = 0; i < DILITHIUM_N; i++) {
            int idx = j + i;
            int sign = 1;
            if (idx >= DILITHIUM_N) { idx -= DILITHIUM_N; sign = -1; }
            s += (int64_t)a->coeffs[i] * g->coeffs[idx] * (int64_t)sign;
        }
        acc->coeffs[j] = fqred(s);
    }
}

/* ===== NTT helpers ===== */
void poly_ntt(poly_t *a) {
    ntt_avx(a->coeffs, qdata);
}

void poly_invntt_tomont(poly_t *a) {
    invntt_avx(a->coeffs, qdata);
}

void poly_pointwise_accum(poly_t *acc, const poly_t *a_ntt, const poly_t *b_ntt) {
    poly_t tmp;
    pointwise_avx(tmp.coeffs, a_ntt->coeffs, b_ntt->coeffs, qdata);
    for (int i = 0; i < DILITHIUM_N; i++) {
        int64_t t = (int64_t)acc->coeffs[i] + tmp.coeffs[i];
        acc->coeffs[i] = fqred(t);
    }
}

int poly_equal_modq(const poly_t *a, const poly_t *b) {
    for (int i = 0; i < DILITHIUM_N; i++) {
        if (fqred(a->coeffs[i]) != fqred(b->coeffs[i])) return 0;
    }
    return 1;
}

void poly_from_mont(poly_t *a) {
    for (int i = 0; i < DILITHIUM_N; i++) {
        int32_t v = montgomery_reduce((int64_t)a->coeffs[i]);
        if (v < 0) v += DILITHIUM_Q;
        a->coeffs[i] = v;
    }
}
