#include "brlt_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

static void* aligned_malloc(size_t size) {
    void *p = NULL;
    if (posix_memalign(&p, 32, size) != 0) return NULL;
    return p;
}

static inline int32_t fqred64_pos(uint64_t x) {
    const uint64_t q = (uint64_t)DILITHIUM_Q;
    const unsigned __int128 mu = ((unsigned __int128)1 << 64) / q;

    uint64_t t = (uint64_t)(((unsigned __int128)x * mu) >> 64);
    uint64_t r = x - t * q;
    if (r >= q) r -= q;
    return (int32_t)r;
}

static inline int32_t fqred64(int64_t x) {
    const uint64_t q = (uint64_t)DILITHIUM_Q;
    const unsigned __int128 mu = ((unsigned __int128)1 << 64) / q;

    int64_t r;
    if (x >= 0) {
        uint64_t t = (uint64_t)(((unsigned __int128)x * mu) >> 64);
        r = (int64_t)x - (int64_t)t * (int64_t)q;
    } else {
        uint64_t t = (uint64_t)(((unsigned __int128)(-x) * mu) >> 64);
        r = (int64_t)(-x) - (int64_t)t * (int64_t)q;
        r = -r;
    }

    if (r >= (int64_t)q) r -= (int64_t)q;
    if (r < 0) r += (int64_t)q;
    return (int32_t)r;
}

#if defined(__AVX2__)
static inline void add_i32_to_i64_avx2(const int32_t *src, int64_t *dst) {
    for (int j = 0; j < DILITHIUM_N; j += 8) {
        __m256i v = _mm256_loadu_si256((const __m256i*)(src + j));

        __m128i v_lo = _mm256_castsi256_si128(v);
        __m128i v_hi = _mm256_extracti128_si256(v, 1);

        __m256i v64_lo = _mm256_cvtepi32_epi64(v_lo);
        __m256i v64_hi = _mm256_cvtepi32_epi64(v_hi);

        __m256i d0 = _mm256_loadu_si256((const __m256i*)(dst + j + 0));
        __m256i d1 = _mm256_loadu_si256((const __m256i*)(dst + j + 4));

        d0 = _mm256_add_epi64(d0, v64_lo);
        d1 = _mm256_add_epi64(d1, v64_hi);

        _mm256_storeu_si256((__m256i*)(dst + j + 0), d0);
        _mm256_storeu_si256((__m256i*)(dst + j + 4), d1);
    }
}
#endif

#if defined(__AVX2__)
/* AVX2: 64x64 -> high64 (unsigned) using 32-bit partial products */
static inline __m256i mulhi_epu64_avx2(__m256i x, __m256i y) {
    const __m256i mask32 = _mm256_set1_epi64x(0xFFFFFFFFULL);

    __m256i x_lo = _mm256_and_si256(x, mask32);
    __m256i x_hi = _mm256_srli_epi64(x, 32);
    __m256i y_lo = _mm256_and_si256(y, mask32);
    __m256i y_hi = _mm256_srli_epi64(y, 32);

    __m256i p0 = _mm256_mul_epu32(x_lo, y_lo);
    __m256i p1 = _mm256_mul_epu32(x_lo, y_hi);
    __m256i p2 = _mm256_mul_epu32(x_hi, y_lo);
    __m256i p3 = _mm256_mul_epu32(x_hi, y_hi);

    __m256i p0_hi = _mm256_srli_epi64(p0, 32);
    __m256i p1_lo = _mm256_and_si256(p1, mask32);
    __m256i p2_lo = _mm256_and_si256(p2, mask32);

    __m256i s = _mm256_add_epi64(p0_hi, _mm256_add_epi64(p1_lo, p2_lo));
    __m256i carry = _mm256_srli_epi64(s, 32);

    __m256i p1_hi = _mm256_srli_epi64(p1, 32);
    __m256i p2_hi = _mm256_srli_epi64(p2, 32);

    __m256i hi = _mm256_add_epi64(p3, _mm256_add_epi64(p1_hi, _mm256_add_epi64(p2_hi, carry)));
    return hi;
}
#endif

/* AVX2 归约（src 非负）；无 AVX2 时走标量 */
static inline void reduce_acc_to_poly(const int64_t *restrict src, int32_t *restrict dst) {
#if defined(__AVX2__)
    const __m256i vq  = _mm256_set1_epi64x((int64_t)DILITHIUM_Q);
    const __m256i vmu = _mm256_set1_epi64x((int64_t)(((unsigned __int128)1 << 64) / (uint64_t)DILITHIUM_Q));
    const __m256i vq_minus1 = _mm256_set1_epi64x((int64_t)DILITHIUM_Q - 1);

    for (int j = 0; j < DILITHIUM_N; j += 4) {
        __m256i x  = _mm256_loadu_si256((const __m256i*)(src + j));
        __m256i t  = mulhi_epu64_avx2(x, vmu);
        __m256i tq = _mm256_mul_epu32(t, vq);
        __m256i r  = _mm256_sub_epi64(x, tq);

        __m256i gt = _mm256_cmpgt_epi64(r, vq_minus1);
        __m256i r2 = _mm256_sub_epi64(r, vq);
        r = _mm256_blendv_epi8(r, r2, gt);

        dst[j + 0] = (int32_t)_mm256_extract_epi64(r, 0);
        dst[j + 1] = (int32_t)_mm256_extract_epi64(r, 1);
        dst[j + 2] = (int32_t)_mm256_extract_epi64(r, 2);
        dst[j + 3] = (int32_t)_mm256_extract_epi64(r, 3);
    }
#else
    for (int j = 0; j < DILITHIUM_N; j++) dst[j] = fqred64_pos((uint64_t)src[j]);
#endif
}

static inline void madd_i64(int64_t *restrict dst, const int64_t *restrict src, int64_t coeff) {
    for (int j = 0; j < DILITHIUM_N; j++) dst[j] += coeff * src[j];
}

static inline void mmul_i64(int64_t *restrict dst, const int64_t *restrict src, int64_t coeff) {
    for (int j = 0; j < DILITHIUM_N; j++) dst[j] = coeff * src[j];
}

int brlt_ctx_init(brlt_ctx_t *ctx, table_t *tb) {
    if (!ctx || !tb) return -1;
    ctx->tb = tb;
    ctx->T = tb->T;
    return 0;
}

int brlt_verify(const brlt_ctx_t *ctx, const poly_grid_t *z, const poly_grid_t *u) {
    if (!ctx || !ctx->tb || !z || !u) return -1;
    const table_t *tb = ctx->tb;
    if (z->rows != u->rows || z->cols != tb->L_ || u->cols != tb->K_) return -1;

    size_t z_bucket_bytes = (size_t)tb->L_ * sizeof(poly_t);
    size_t u_bucket_bytes = (size_t)tb->K_ * sizeof(poly_t);
    size_t zsum_size = (size_t)tb->L_ * tb->T * DILITHIUM_N;  // [L][T][N]
    size_t usum_size = (size_t)tb->K_ * tb->T * DILITHIUM_N;  // [K][T][N]
    size_t zsum_bytes = zsum_size * sizeof(int64_t);
    size_t usum_bytes = usum_size * sizeof(int64_t);
    size_t zacc_size = (size_t)tb->L_ * DILITHIUM_N;
    size_t uacc_size = (size_t)tb->K_ * DILITHIUM_N;
    size_t zacc_bytes = zacc_size * sizeof(int64_t);
    size_t uacc_bytes = uacc_size * sizeof(int64_t);

    poly_t *z_bucket = aligned_malloc(z_bucket_bytes);
    poly_t *u_bucket = aligned_malloc(u_bucket_bytes);
    int64_t *z_sum = aligned_malloc(zsum_bytes);
    int64_t *u_sum = aligned_malloc(usum_bytes);
    int64_t *z_acc = aligned_malloc(zacc_bytes);
    int64_t *u_acc = aligned_malloc(uacc_bytes);
    if (!z_bucket || !u_bucket || !z_sum || !u_sum || !z_acc || !u_acc) {
        free(z_bucket); free(u_bucket); free(z_sum);
        free(u_sum); free(z_acc); free(u_acc);
        return -1;
    }

    memset(z_sum, 0, zsum_bytes);
    memset(u_sum, 0, usum_bytes);

    // accumulate into [L][T][N] and [K][T][N]
    for (uint32_t t = 0; t < tb->T; t++) {
        for (uint32_t i = t; i < z->rows; i += tb->T) {
            const poly_t *zrow = &z->data[i * tb->L_];
            const poly_t *urow = &u->data[i * tb->K_];

            for (uint32_t l = 0; l < tb->L_; l++) {
                const int32_t *src = zrow[l].coeffs;
                int64_t *dst = z_sum + ((size_t)l * tb->T + t) * DILITHIUM_N;
#if defined(__AVX2__)
                add_i32_to_i64_avx2(src, dst);
#else
                for (int j = 0; j < DILITHIUM_N; j++) dst[j] += src[j];
#endif
            }
            for (uint32_t k = 0; k < tb->K_; k++) {
                const int32_t *src = urow[k].coeffs;
                int64_t *dst = u_sum + ((size_t)k * tb->T + t) * DILITHIUM_N;
#if defined(__AVX2__)
                add_i32_to_i64_avx2(src, dst);
#else
                for (int j = 0; j < DILITHIUM_N; j++) dst[j] += src[j];
#endif
            }
        }
    }

    int ok = 1;
    for (uint32_t rid = 0; rid < tb->R; rid++) {
        if (tb->T == 0) {
            memset(z_acc, 0, zacc_bytes);
            memset(u_acc, 0, uacc_bytes);
        } else {
            // z_acc: first t writes, rest t adds
            for (uint32_t l = 0; l < tb->L_; l++) {
                int64_t *dst = z_acc + (size_t)l * DILITHIUM_N;
                int64_t coeff0 = (int64_t)tb->randoms[rid * tb->T + 0];
                const int64_t *src0 = z_sum + ((size_t)l * tb->T + 0) * DILITHIUM_N;
                mmul_i64(dst, src0, coeff0);
                for (uint32_t t = 1; t < tb->T; t++) {
                    int64_t coeff = (int64_t)tb->randoms[rid * tb->T + t];
                    const int64_t *src = z_sum + ((size_t)l * tb->T + t) * DILITHIUM_N;
                    madd_i64(dst, src, coeff);
                }
            }

            // u_acc: first t writes, rest t adds
            for (uint32_t k = 0; k < tb->K_; k++) {
                int64_t *dst = u_acc + (size_t)k * DILITHIUM_N;
                int64_t coeff0 = (int64_t)tb->randoms[rid * tb->T + 0];
                const int64_t *src0 = u_sum + ((size_t)k * tb->T + 0) * DILITHIUM_N;
                mmul_i64(dst, src0, coeff0);
                for (uint32_t t = 1; t < tb->T; t++) {
                    int64_t coeff = (int64_t)tb->randoms[rid * tb->T + t];
                    const int64_t *src = u_sum + ((size_t)k * tb->T + t) * DILITHIUM_N;
                    madd_i64(dst, src, coeff);
                }
            }
        }

        for (uint32_t l = 0; l < tb->L_; l++) {
            int64_t *src = z_acc + l * DILITHIUM_N;
            reduce_acc_to_poly(src, z_bucket[l].coeffs);
        }
        for (uint32_t k = 0; k < tb->K_; k++) {
            int64_t *src = u_acc + k * DILITHIUM_N;
            reduce_acc_to_poly(src, u_bucket[k].coeffs);
        }

        for (uint32_t l = 0; l < tb->L_; l++) poly_ntt(&z_bucket[l]);
        for (uint32_t k = 0; k < tb->K_; k++) poly_ntt(&u_bucket[k]);

        poly_t LHS, RHS;
        poly_zero(&LHS);
        poly_zero(&RHS);

        for (uint32_t l = 0; l < tb->L_; l++) {
            poly_pointwise_accum(&LHS, &z_bucket[l], &tb->ATgamma[rid * tb->L_ + l]);
        }
        for (uint32_t k = 0; k < tb->K_; k++) {
            poly_pointwise_accum(&RHS, &u_bucket[k], &tb->gamma[rid * tb->K_ + k]);
        }

        poly_invntt_tomont(&LHS);
        poly_invntt_tomont(&RHS);
        poly_from_mont(&LHS);
        poly_from_mont(&RHS);

        if (!poly_equal_modq(&LHS, &RHS)) {
            fprintf(stderr, "rid=%u FAIL\n", rid);
            ok = 0;
            break;
        }
    }

    free(z_bucket);
    free(u_bucket);
    free(z_sum);
    free(u_sum);
    free(z_acc);
    free(u_acc);
    return ok;
}
