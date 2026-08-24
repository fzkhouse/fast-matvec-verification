#include "poly.h"
#include <stdlib.h>
#include <string.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

static inline uint32_t mod_add(uint32_t a, uint32_t b, uint32_t q) {
    uint64_t t = (uint64_t)a + b;
    if (t >= q) t -= q;
    return (uint32_t)t;
}

static inline uint32_t mod_sub(uint32_t a, uint32_t b, uint32_t q) {
    return (a >= b) ? (a - b) : (a + q - b);
}

static inline uint32_t mod_mul(uint32_t a, uint32_t b, uint32_t q) {
    return (uint32_t)(((uint64_t)a * (uint64_t)b) % (uint64_t)q);
}

int poly_init(poly_t *p, uint32_t n) {
    if (!p || n == 0) return -1;
    p->n = n;
    p->c = (uint32_t *)calloc(n, sizeof(uint32_t));
    return p->c ? 0 : -1;
}

void poly_free(poly_t *p) {
    if (!p) return;
    free(p->c);
    p->c = NULL;
    p->n = 0;
}

void poly_zero(poly_t *p) {
    memset(p->c, 0, (size_t)p->n * sizeof(uint32_t));
}

void poly_copy(poly_t *dst, const poly_t *src) {
    memcpy(dst->c, src->c, (size_t)src->n * sizeof(uint32_t));
}

void poly_rand(poly_t *p, uint32_t q) {
    for (uint32_t i = 0; i < p->n; ++i) {
        p->c[i] = (uint32_t)(rand() % (int)q);
    }
}

void poly_add_inplace(poly_t *a, const poly_t *b, uint32_t q) {
#ifdef __AVX2__
    uint32_t i = 0;
    const __m256i vq = _mm256_set1_epi32((int)q);
    const __m256i vq_minus_1 = _mm256_set1_epi32((int)(q - 1u));

    for (; i + 8 <= a->n; i += 8) {
        __m256i va = _mm256_loadu_si256((const __m256i *)&a->c[i]);
        __m256i vb = _mm256_loadu_si256((const __m256i *)&b->c[i]);

        __m256i vsum = _mm256_add_epi32(va, vb);      // < 2q
        __m256i vsub = _mm256_sub_epi32(vsum, vq);    // vsum - q
        __m256i mask = _mm256_cmpgt_epi32(vsum, vq_minus_1); // vsum >= q

        __m256i vres = _mm256_blendv_epi8(vsum, vsub, mask);
        _mm256_storeu_si256((__m256i *)&a->c[i], vres);
    }
    for (; i < a->n; ++i) {
        a->c[i] = mod_add(a->c[i], b->c[i], q);
    }
#else
    for (uint32_t i = 0; i < a->n; ++i) {
        a->c[i] = mod_add(a->c[i], b->c[i], q);
    }
#endif
}

void poly_sub_inplace(poly_t *a, const poly_t *b, uint32_t q) {
#ifdef __AVX2__
    uint32_t i = 0;
    const __m256i vq = _mm256_set1_epi32((int)q);

    for (; i + 8 <= a->n; i += 8) {
        __m256i va = _mm256_loadu_si256((const __m256i *)&a->c[i]);
        __m256i vb = _mm256_loadu_si256((const __m256i *)&b->c[i]);

        __m256i vdiff = _mm256_sub_epi32(va, vb);
        __m256i mask = _mm256_cmpgt_epi32(vb, va); // va < vb
        __m256i vaddq = _mm256_and_si256(mask, vq);
        __m256i vres = _mm256_add_epi32(vdiff, vaddq);

        _mm256_storeu_si256((__m256i *)&a->c[i], vres);
    }
    for (; i < a->n; ++i) {
        a->c[i] = mod_sub(a->c[i], b->c[i], q);
    }
#else
    for (uint32_t i = 0; i < a->n; ++i) {
        a->c[i] = mod_sub(a->c[i], b->c[i], q);
    }
#endif
}

void poly_add_scaled_inplace(poly_t *acc, const poly_t *x, uint32_t s, uint32_t q) {
    for (uint32_t i = 0; i < acc->n; ++i) {
        uint32_t t = mod_mul(x->c[i], s, q);
        acc->c[i] = mod_add(acc->c[i], t, q);
    }
}

int poly_ntt_inplace(poly_t *p, const ntt_ctx_t *ctx) {
    return ntt_forward(p->c, ctx);
}

int poly_intt_inplace(poly_t *p, const ntt_ctx_t *ctx) {
    return ntt_inverse(p->c, ctx);
}

void poly_hadamard_accumulate(poly_t *acc, const poly_t *a, const poly_t *b, uint32_t q) {
#ifdef __AVX2__
    uint32_t i = 0;
    for (; i + 8 <= acc->n; i += 8) {
        __m256i va = _mm256_loadu_si256((const __m256i *)&a->c[i]);
        __m256i vb = _mm256_loadu_si256((const __m256i *)&b->c[i]);

        __m256i prod_even = _mm256_mul_epu32(va, vb);
        __m256i va_odd = _mm256_srli_epi64(va, 32);
        __m256i vb_odd = _mm256_srli_epi64(vb, 32);
        __m256i prod_odd = _mm256_mul_epu32(va_odd, vb_odd);

        uint64_t pe[4], po[4];
        _mm256_storeu_si256((__m256i *)pe, prod_even);
        _mm256_storeu_si256((__m256i *)po, prod_odd);

        acc->c[i + 0] = mod_add(acc->c[i + 0], (uint32_t)(pe[0] % q), q);
        acc->c[i + 1] = mod_add(acc->c[i + 1], (uint32_t)(po[0] % q), q);
        acc->c[i + 2] = mod_add(acc->c[i + 2], (uint32_t)(pe[1] % q), q);
        acc->c[i + 3] = mod_add(acc->c[i + 3], (uint32_t)(po[1] % q), q);
        acc->c[i + 4] = mod_add(acc->c[i + 4], (uint32_t)(pe[2] % q), q);
        acc->c[i + 5] = mod_add(acc->c[i + 5], (uint32_t)(po[2] % q), q);
        acc->c[i + 6] = mod_add(acc->c[i + 6], (uint32_t)(pe[3] % q), q);
        acc->c[i + 7] = mod_add(acc->c[i + 7], (uint32_t)(po[3] % q), q);
    }
    for (; i < acc->n; ++i) {
        acc->c[i] = mod_add(acc->c[i], mod_mul(a->c[i], b->c[i], q), q);
    }
#else
    for (uint32_t i = 0; i < acc->n; ++i) {
        acc->c[i] = mod_add(acc->c[i], mod_mul(a->c[i], b->c[i], q), q);
    }
#endif
}

int poly_equal(const poly_t *a, const poly_t *b) {
    if (a->n != b->n) return 0;
    return memcmp(a->c, b->c, (size_t)a->n * sizeof(uint32_t)) == 0;
}