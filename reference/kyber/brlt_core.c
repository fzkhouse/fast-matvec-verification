#include "brlt_core.h"
#include "kyber_avx2_adapter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

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

#define GAT(g,r,c) ((g)->data + (((size_t)(r) * (g)->cols + (c)) * (g)->n))
#define GAT16(g,r,c) ((g)->data + (((size_t)(r) * (g)->cols + (c)) * (g)->n))

#if defined(KYBER_REFERENCE_USE_ASM)
extern void brlt_accum_u32_avx2(uint32_t *acc, const uint32_t *src, uint32_t n);
extern void brlt_accum_r_u64_u32_asm(uint64_t *acc, const uint32_t *src,
                                     uint32_t n, uint64_t r);
#endif

static inline uint32_t mod_q(uint64_t x) { return (uint32_t)(x % KYB_Q); }

static inline uint32_t xorshift32(uint32_t *s) {
    uint32_t x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}

static int read_u32(FILE *fp, uint32_t *x) {
    return (fscanf(fp, "%u", x) == 1) ? 0 : -1;
}

static int grid_init(poly_grid_u32_t *g, uint32_t rows, uint32_t cols, uint32_t n) {
    g->rows = rows;
    g->cols = cols;
    g->n = n;
    g->data = (uint32_t *)calloc((size_t)rows * cols * n, sizeof(uint32_t));
    return g->data ? 0 : -1;
}

void grid_free_u32(poly_grid_u32_t *g) {
    if (!g) return;
    free(g->data);
    memset(g, 0, sizeof(*g));
}

void grid_free_u16(poly_grid_u16_t *g) {
    if (!g) return;
    free(g->data);
    memset(g, 0, sizeof(*g));
}

int grid_load_u32(const char *path, poly_grid_u32_t *g, uint32_t expect_q, uint32_t expect_n) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    uint32_t q, n, rows, cols;
    if (read_u32(fp, &q) || read_u32(fp, &n) || read_u32(fp, &rows) || read_u32(fp, &cols)) {
        fclose(fp);
        return -1;
    }
    if (q != expect_q || n != expect_n) {
        fclose(fp);
        return -1;
    }
    if (grid_init(g, rows, cols, n) != 0) {
        fclose(fp);
        return -1;
    }

    for (uint32_t r = 0; r < rows; r++) {
        for (uint32_t c = 0; c < cols; c++) {
            uint32_t *dst = GAT(g, r, c);
            for (uint32_t i = 0; i < n; i++) {
                uint32_t v;
                if (read_u32(fp, &v)) {
                    fclose(fp);
                    return -1;
                }
                dst[i] = v % q;
            }
        }
    }

    fclose(fp);
    return 0;
}

int grid_load_u16(const char *path, poly_grid_u16_t *g, uint32_t expect_q, uint32_t expect_n) {
    memset(g, 0, sizeof(*g));
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    uint32_t q;
    if (read_u32(fp, &q) || read_u32(fp, &g->n) ||
        read_u32(fp, &g->rows) || read_u32(fp, &g->cols) ||
        q != expect_q || g->n != expect_n) {
        fclose(fp);
        return -1;
    }
    const size_t count = (size_t)g->rows * g->cols * g->n;
    g->data = calloc(count, sizeof(*g->data));
    if (!g->data) { fclose(fp); return -1; }
    for (size_t i = 0; i < count; i++) {
        uint32_t value;
        if (read_u32(fp, &value)) { grid_free_u16(g); fclose(fp); return -1; }
        g->data[i] = (uint16_t)(value % q);
    }
    fclose(fp);
    return 0;
}

void table_free_u32(table_u32_t *tb) {
    if (!tb) return;
    grid_free_u32(&tb->gamma_ntt);
    grid_free_u32(&tb->ATg_ntt);
    memset(tb, 0, sizeof(*tb));
}

int table_load_u32(const char *path, table_u32_t *tb, uint32_t expect_q, uint32_t expect_n) {
    memset(tb, 0, sizeof(*tb));

    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    uint32_t q, n, T_in;
    if (read_u32(fp, &q) || read_u32(fp, &n) || read_u32(fp, &tb->R) ||
        read_u32(fp, &T_in) || read_u32(fp, &tb->K) || read_u32(fp, &tb->L)) {
        fclose(fp);
        return -1;
    }

    tb->file_T = T_in;
    if (q != expect_q || n != expect_n) {
        fclose(fp);
        return -1;
    }

    if (grid_init(&tb->gamma_ntt, tb->R, tb->K, n) != 0) {
        fclose(fp);
        return -1;
    }
    if (grid_init(&tb->ATg_ntt, tb->R, tb->L, n) != 0) {
        fclose(fp);
        return -1;
    }

    for (uint32_t r = 0; r < tb->R; r++) {
        for (uint32_t t = 0; t < T_in; t++) {
            uint32_t dummy;
            if (read_u32(fp, &dummy)) {
                fclose(fp);
                return -1;
            }
        }

        for (uint32_t k = 0; k < tb->K; k++) {
            uint32_t *dst = GAT(&tb->gamma_ntt, r, k);
            for (uint32_t i = 0; i < n; i++) {
                uint32_t v;
                if (read_u32(fp, &v)) {
                    fclose(fp);
                    return -1;
                }
                dst[i] = v % q;
            }
        }

        for (uint32_t l = 0; l < tb->L; l++) {
            uint32_t *dst = GAT(&tb->ATg_ntt, r, l);
            for (uint32_t i = 0; i < n; i++) {
                uint32_t v;
                if (read_u32(fp, &v)) {
                    fclose(fp);
                    return -1;
                }
                dst[i] = v % q;
            }
        }
    }

    fclose(fp);
    return 0;
}

static int aligned_zalloc(void **ptr, size_t align, size_t bytes) {
    *ptr = NULL;
    if (posix_memalign(ptr, align, bytes) != 0) return -1;
    memset(*ptr, 0, bytes);
    return 0;
}

static inline void accum_u32(uint32_t *acc, const uint32_t *src, uint32_t n) {
#if defined(KYBER_REFERENCE_USE_ASM)
    brlt_accum_u32_avx2(acc, src, n);
#else
#ifdef __AVX2__
    uint32_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256i a = _mm256_load_si256((const __m256i *)(acc + i));
        __m256i b = _mm256_loadu_si256((const __m256i *)(src + i));
        _mm256_store_si256((__m256i *)(acc + i), _mm256_add_epi32(a, b));
    }
    for (; i < n; i++) acc[i] += src[i];
#else
    for (uint32_t i = 0; i < n; i++) acc[i] += src[i];
#endif
#endif
}

/* Generic four-way AVX2 update used by the cache-bounded round group below.
 * `count` is runtime data, not a Kyber parameter specialization. */
static inline void accum_u32_group(uint32_t *const acc[4], uint32_t count,
                                   const uint32_t *src, uint32_t n) {
#ifdef __AVX2__
    uint32_t i = 0;
    for (; i + 8 <= n; i += 8) {
        const __m256i vsrc = _mm256_loadu_si256((const __m256i *)(src + i));
        for (uint32_t slot = 0; slot < count; slot++) {
            const __m256i va = _mm256_loadu_si256((const __m256i *)(acc[slot] + i));
            _mm256_storeu_si256((__m256i *)(acc[slot] + i), _mm256_add_epi32(va, vsrc));
        }
    }
    for (; i < n; i++) {
        for (uint32_t slot = 0; slot < count; slot++) acc[slot][i] += src[i];
    }
#else
    for (uint32_t slot = 0; slot < count; slot++) accum_u32(acc[slot], src, n);
#endif
}

static inline void accum_u32_four(uint32_t *a0, uint32_t *a1, uint32_t *a2,
                                  uint32_t *a3, const uint32_t *src, uint32_t n) {
#ifdef __AVX2__
    uint32_t i = 0;
    for (; i + 8 <= n; i += 8) {
        const __m256i s = _mm256_loadu_si256((const __m256i *)(src + i));
        __m256i x = _mm256_loadu_si256((const __m256i *)(a0 + i));
        _mm256_storeu_si256((__m256i *)(a0 + i), _mm256_add_epi32(x, s));
        x = _mm256_loadu_si256((const __m256i *)(a1 + i));
        _mm256_storeu_si256((__m256i *)(a1 + i), _mm256_add_epi32(x, s));
        x = _mm256_loadu_si256((const __m256i *)(a2 + i));
        _mm256_storeu_si256((__m256i *)(a2 + i), _mm256_add_epi32(x, s));
        x = _mm256_loadu_si256((const __m256i *)(a3 + i));
        _mm256_storeu_si256((__m256i *)(a3 + i), _mm256_add_epi32(x, s));
    }
    for (; i < n; i++) { a0[i] += src[i]; a1[i] += src[i]; a2[i] += src[i]; a3[i] += src[i]; }
#else
    accum_u32(a0, src, n); accum_u32(a1, src, n);
    accum_u32(a2, src, n); accum_u32(a3, src, n);
#endif
}

/* A safe compact bucket representation: every update is reduced modulo q,
 * unlike a wrapping uint16 accumulator. */
static inline void accum_u16_mod_four(uint16_t *a0, uint16_t *a1, uint16_t *a2,
                                      uint16_t *a3, const uint16_t *src, uint32_t n) {
#ifdef __AVX2__
    const __m256i vq = _mm256_set1_epi16((short)KYB_Q);
    const __m256i vq1 = _mm256_set1_epi16((short)(KYB_Q - 1));
    uint32_t i = 0;
    for (; i + 16 <= n; i += 16) {
        const __m256i s = _mm256_loadu_si256((const __m256i *)(src + i));
        uint16_t *dsts[4] = { a0, a1, a2, a3 };
        for (uint32_t slot = 0; slot < 4; slot++) {
            const __m256i old = _mm256_loadu_si256((const __m256i *)(dsts[slot] + i));
            const __m256i sum = _mm256_add_epi16(old, s);
            const __m256i reduced = _mm256_sub_epi16(sum, vq);
            const __m256i geq = _mm256_cmpgt_epi16(sum, vq1);
            _mm256_storeu_si256((__m256i *)(dsts[slot] + i),
                                _mm256_blendv_epi8(sum, reduced, geq));
        }
    }
    for (; i < n; i++) {
        const uint16_t v = (uint16_t)src[i];
        uint16_t *dsts[4] = { a0, a1, a2, a3 };
        for (uint32_t slot = 0; slot < 4; slot++) {
            uint16_t sum = (uint16_t)(dsts[slot][i] + v);
            dsts[slot][i] = (sum >= KYB_Q) ? (uint16_t)(sum - KYB_Q) : sum;
        }
    }
#else
    for (uint32_t i = 0; i < n; i++) {
        const uint16_t v = (uint16_t)src[i];
        uint16_t *dsts[4] = { a0, a1, a2, a3 };
        for (uint32_t slot = 0; slot < 4; slot++) {
            uint16_t sum = (uint16_t)(dsts[slot][i] + v);
            dsts[slot][i] = (sum >= KYB_Q) ? (uint16_t)(sum - KYB_Q) : sum;
        }
    }
#endif
}

static inline void accum_u16_mod_group(uint16_t *const dst[4], uint32_t count,
                                       const uint16_t *src, uint32_t n) {
#ifdef __AVX2__
    const __m256i vq = _mm256_set1_epi16((short)KYB_Q);
    const __m256i vq1 = _mm256_set1_epi16((short)(KYB_Q - 1));
    uint32_t i = 0;
    for (; i + 16 <= n; i += 16) {
        const __m256i s = _mm256_loadu_si256((const __m256i *)(src + i));
        for (uint32_t slot = 0; slot < count; slot++) {
            const __m256i old = _mm256_loadu_si256((const __m256i *)(dst[slot] + i));
            const __m256i sum = _mm256_add_epi16(old, s);
            const __m256i reduced = _mm256_sub_epi16(sum, vq);
            const __m256i geq = _mm256_cmpgt_epi16(sum, vq1);
            _mm256_storeu_si256((__m256i *)(dst[slot] + i),
                                _mm256_blendv_epi8(sum, reduced, geq));
        }
    }
    for (; i < n; i++) {
        const uint16_t v = (uint16_t)src[i];
        for (uint32_t slot = 0; slot < count; slot++) {
            uint16_t sum = (uint16_t)(dst[slot][i] + v);
            dst[slot][i] = (sum >= KYB_Q) ? (uint16_t)(sum - KYB_Q) : sum;
        }
    }
#else
    for (uint32_t slot = 0; slot < count; slot++) {
        for (uint32_t i = 0; i < n; i++) {
            uint16_t sum = (uint16_t)(dst[slot][i] + src[i]);
            dst[slot][i] = (sum >= KYB_Q) ? (uint16_t)(sum - KYB_Q) : sum;
        }
    }
#endif
}

static inline void accum_r_u64(uint64_t *acc, const uint32_t *src, uint32_t n, uint64_t r) {
#if defined(KYBER_REFERENCE_USE_ASM)
    brlt_accum_r_u64_u32_asm(acc, src, n, r);
#else
#ifdef __AVX2__
    const __m256i vr = _mm256_set1_epi64x((long long)r);
    uint32_t i = 0;
    for (; i + 8 <= n; i += 8) {
        const __m256i vsrc = _mm256_loadu_si256((const __m256i *)(src + i));
        const __m128i lo = _mm256_castsi256_si128(vsrc);
        const __m128i hi = _mm256_extracti128_si256(vsrc, 1);
        const __m256i src_lo = _mm256_cvtepu32_epi64(lo);
        const __m256i src_hi = _mm256_cvtepu32_epi64(hi);
        const __m256i acc_lo = _mm256_loadu_si256((const __m256i *)(acc + i));
        const __m256i acc_hi = _mm256_loadu_si256((const __m256i *)(acc + i + 4));
        _mm256_storeu_si256((__m256i *)(acc + i),
                             _mm256_add_epi64(acc_lo, _mm256_mul_epu32(src_lo, vr)));
        _mm256_storeu_si256((__m256i *)(acc + i + 4),
                             _mm256_add_epi64(acc_hi, _mm256_mul_epu32(src_hi, vr)));
    }
    for (; i < n; i++) acc[i] += r * (uint64_t)src[i];
#else
    for (uint32_t i = 0; i < n; i++) acc[i] += r * (uint64_t)src[i];
#endif
#endif
}

static void derive_random_vec(uint32_t *rv, uint32_t T, uint32_t rid) {
    uint32_t s = 0x9E3779B9u ^ (rid + 1u);
    for (uint32_t i = 0; i < T; i++) rv[i] = xorshift32(&s) % KYB_Q;
}

int brlt_ctx_init_kyber512(brlt_ctx_t *ctx, const table_u32_t *tb, uint32_t T) {
    memset(ctx, 0, sizeof(*ctx));

    if (!tb) return -1;
    if (tb->K == 0 || tb->L == 0) return -1;
    if (tb->gamma_ntt.n != KYB_N || tb->ATg_ntt.n != KYB_N) return -1;
    if (T == 0 || (T & (T - 1u)) != 0) return -1; // power-of-two

    ctx->T = T;
    ctx->n = KYB_N;
    ctx->K = tb->K;
    ctx->L = tb->L;
    ctx->R = tb->R;

    size_t sz_acc_c  = (size_t)ctx->T * ctx->L * ctx->n * sizeof(uint32_t);
    size_t sz_acc_d  = (size_t)ctx->T * ctx->K * ctx->n * sizeof(uint32_t);
    size_t sz_lazy_z = (size_t)ctx->L * ctx->n * sizeof(uint64_t);
    size_t sz_lazy_u = (size_t)ctx->K * ctx->n * sizeof(uint64_t);
    size_t sz_rv     = (size_t)ctx->T * sizeof(uint32_t);
    ctx->group_words = (size_t)ctx->T * (ctx->L + ctx->K) * ctx->n;
    size_t sz_group  = 4u * ctx->group_words * sizeof(uint16_t);
    size_t sz_gamma  = (size_t)ctx->R * ctx->K * sizeof(poly);
    size_t sz_atg    = (size_t)ctx->R * ctx->L * sizeof(poly);

    if (aligned_zalloc((void **)&ctx->acc_c, 32, sz_acc_c) != 0 ||
        aligned_zalloc((void **)&ctx->acc_d, 32, sz_acc_d) != 0 ||
        aligned_zalloc((void **)&ctx->lazy_z, 32, sz_lazy_z) != 0 ||
        aligned_zalloc((void **)&ctx->lazy_u, 32, sz_lazy_u) != 0 ||
        aligned_zalloc((void **)&ctx->rv, 32, sz_rv) != 0 ||
        aligned_zalloc((void **)&ctx->group, 32, sz_group) != 0 ||
        aligned_zalloc((void **)&ctx->gamma_poly, 32, sz_gamma) != 0 ||
        aligned_zalloc((void **)&ctx->ATg_poly, 32, sz_atg) != 0) {
        brlt_ctx_free(ctx);
        return -1;
    }

    for (uint32_t r = 0; r < ctx->R; r++) {
        for (uint32_t k = 0; k < ctx->K; k++) {
            kyb_poly_from_u32(&ctx->gamma_poly[(size_t)r * ctx->K + k], GAT(&tb->gamma_ntt, r, k));
        }
        for (uint32_t l = 0; l < ctx->L; l++) {
            kyb_poly_from_u32(&ctx->ATg_poly[(size_t)r * ctx->L + l], GAT(&tb->ATg_ntt, r, l));
        }
    }

    return 0;
}

void brlt_ctx_free(brlt_ctx_t *ctx) {
    if (!ctx) return;
    free(ctx->acc_c);      ctx->acc_c = NULL;
    free(ctx->acc_d);      ctx->acc_d = NULL;
    free(ctx->lazy_z);     ctx->lazy_z = NULL;
    free(ctx->lazy_u);     ctx->lazy_u = NULL;
    free(ctx->rv);         ctx->rv = NULL;
    free(ctx->group);      ctx->group = NULL;
    free(ctx->gamma_poly); ctx->gamma_poly = NULL;
    free(ctx->ATg_poly);   ctx->ATg_poly = NULL;
    memset(ctx, 0, sizeof(*ctx));
}

int brlt_verify_kyber512_ctx(
    brlt_ctx_t *ctx,
    const poly_grid_u16_t *z,
    const poly_grid_u16_t *u,
    uint32_t bucket_times,
    uint32_t brlt_times
) {
    if (!ctx || !z || !u) return -1;
    if (ctx->n != z->n || ctx->n != u->n) return -1;
    if (ctx->L != z->cols || ctx->K != u->cols) return -1;
    if (z->rows != u->rows) return -1;
    if (ctx->R < bucket_times * brlt_times) return -1;

    const uint32_t n = ctx->n, K = ctx->K, L = ctx->L, T = ctx->T;
    const size_t c_words = (size_t)T * L * n;
    const size_t d_words = (size_t)T * K * n;
    const size_t round_words = c_words + d_words;
    uint32_t tmp[KYB_N];
    uint32_t rng = 1u;
    int ok = 1;

    /*
     * Cache-bounded generic fusion.  Four independently checked rounds share
     * one pass over z/u, but retain 32-bit accumulators so no wraparound or
     * scheme-specific packing is involved.
     */
    for (uint32_t first = 0; first < bucket_times && ok; first += 4) {
        const uint32_t count = (bucket_times - first < 4) ? bucket_times - first : 4;
        /* Allocated before timing.  Only the live round slots are cleared. */
        uint16_t *group = ctx->group;
        memset(group, 0, (size_t)count * round_words * sizeof(*group));

        for (uint32_t i = 0; i < z->rows; i++) {
            uint32_t bucket[4];
            uint16_t *round[4];
            for (uint32_t slot = 0; slot < count; slot++) {
                round[slot] = group + (size_t)slot * round_words;
                bucket[slot] = xorshift32(&rng) & (T - 1u);
            }
            for (uint32_t l = 0; l < L; l++) {
                uint16_t *dst[4];
                for (uint32_t slot = 0; slot < count; slot++) {
                    dst[slot] = round[slot] + ((size_t)bucket[slot] * L + l) * n;
                }
                if (count == 4) accum_u16_mod_four(dst[0], dst[1], dst[2], dst[3], GAT16(z, i, l), n);
                else accum_u16_mod_group(dst, count, GAT16(z, i, l), n);
            }
            for (uint32_t k = 0; k < K; k++) {
                uint16_t *dst[4];
                for (uint32_t slot = 0; slot < count; slot++) {
                    dst[slot] = round[slot] + c_words + ((size_t)bucket[slot] * K + k) * n;
                }
                if (count == 4) accum_u16_mod_four(dst[0], dst[1], dst[2], dst[3], GAT16(u, i, k), n);
                else accum_u16_mod_group(dst, count, GAT16(u, i, k), n);
            }
        }

        for (uint32_t slot = 0; slot < count && ok; slot++) {
            const uint32_t bi = first + slot;
            const uint16_t *round = group + (size_t)slot * round_words;
            for (size_t j = 0; j < c_words; j++) ctx->acc_c[j] = round[j];
            for (size_t j = 0; j < d_words; j++) ctx->acc_d[j] = round[c_words + j];

            for (uint32_t r = 0; r < brlt_times && ok; r++) {
                const uint32_t rid = bi * brlt_times + r;
                derive_random_vec(ctx->rv, T, rid);
                memset(ctx->lazy_z, 0, (size_t)L * n * sizeof(uint64_t));
                memset(ctx->lazy_u, 0, (size_t)K * n * sizeof(uint64_t));
                for (uint32_t b = 0; b < T; b++) {
                    const uint64_t rr = ctx->rv[b];
                    for (uint32_t l = 0; l < L; l++) {
                        accum_r_u64(&ctx->lazy_z[(size_t)l * n],
                                    &ctx->acc_c[((size_t)b * L + l) * n], n, rr);
                    }
                    for (uint32_t k = 0; k < K; k++) {
                        accum_r_u64(&ctx->lazy_u[(size_t)k * n],
                                    &ctx->acc_d[((size_t)b * K + k) * n], n, rr);
                    }
                }

                ALIGN32 poly sumz_ntt[L], sumu_ntt[K], lhs, rhs;
                kyb_poly_zero(&lhs);
                kyb_poly_zero(&rhs);
                for (uint32_t l = 0; l < L; l++) {
                    for (uint32_t j = 0; j < n; j++) {
                        tmp[j] = mod_q(ctx->lazy_z[(size_t)l * n + j]);
                    }
                    kyb_poly_from_u32(&sumz_ntt[l], tmp);
                    kyb_poly_ntt(&sumz_ntt[l]);
                }
                for (uint32_t k = 0; k < K; k++) {
                    for (uint32_t j = 0; j < n; j++) {
                        tmp[j] = mod_q(ctx->lazy_u[(size_t)k * n + j]);
                    }
                    kyb_poly_from_u32(&sumu_ntt[k], tmp);
                    kyb_poly_ntt(&sumu_ntt[k]);
                }
                for (uint32_t l = 0; l < L; l++) {
                    const poly *atl = &ctx->ATg_poly[(size_t)rid * L + l];
                    kyb_poly_basemul_acc(&lhs, &sumz_ntt[l], atl);
                }
                for (uint32_t k = 0; k < K; k++) {
                    const poly *gk = &ctx->gamma_poly[(size_t)rid * K + k];
                    kyb_poly_basemul_acc(&rhs, &sumu_ntt[k], gk);
                }
                if (!kyb_poly_equal_modq(&lhs, &rhs)) ok = 0;
            }
        }
    }

    return ok ? 1 : 0;
}
