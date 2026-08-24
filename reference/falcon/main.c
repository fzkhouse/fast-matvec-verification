#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>

#include "params.h"
#include "poly.h"
#include "ntt.h"

#ifdef __AVX2__
#include <immintrin.h>
#endif

#if defined(__x86_64__) || defined(__i386__)
#define HAVE_RDTSC 1
#else
#define HAVE_RDTSC 0
#endif

#define CELL(g, r, c) ((g)->cell[(size_t)(r) * (g)->cols + (c)])

typedef struct {
    uint32_t rows, cols;
    poly_t *cell;
} poly_grid_t;

typedef struct {
    uint32_t R, file_T, K, L;  // file_T only for parsing/info; runtime does not depend on it
    poly_grid_t gamma_ntt;     // R x K
    poly_grid_t ATg_ntt;       // R x L
} table_data_t;

typedef struct {
    uint32_t T, L, K, n;
    uint32_t q;
    uint64_t qinv;
    poly_t *sum_z;
    poly_t *sum_u;
    uint64_t *lazy_z; // L * n
    uint64_t *lazy_u; // K * n
} brlt_ws_t;

typedef struct {
    uint32_t T, L, K, n;
    uint32_t *acc_c; // T * L * n
    uint32_t *acc_d; // T * K * n
} bucket_ws_t;

static inline void accum_src_u32_inplace(uint32_t *restrict acc,
                                         const uint32_t *restrict src,
                                         uint32_t n) {
#ifdef __AVX2__
    uint32_t j = 0;
    for (; j + 8 <= n; j += 8) {
        __m256i va = _mm256_loadu_si256((const __m256i *)(acc + j));
        __m256i vb = _mm256_loadu_si256((const __m256i *)(src + j));
        va = _mm256_add_epi32(va, vb);
        _mm256_storeu_si256((__m256i *)(acc + j), va);
    }
    for (; j < n; ++j) acc[j] += src[j];
#else
    for (uint32_t j = 0; j < n; ++j) acc[j] += src[j];
#endif
}

static inline uint32_t add_mod_q(uint32_t a, uint32_t b, uint32_t q) {
    uint32_t s = a + b;
    if (s >= q) s -= q;
    return s;
}

static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

static inline uint64_t now_cycles(void) {
#if HAVE_RDTSC
    uint32_t hi, lo;
    __asm__ __volatile__("lfence\nrdtsc" : "=a"(lo), "=d"(hi) :: "memory");
    return ((uint64_t)hi << 32) | lo;
#else
    return 0;
#endif
}

static inline uint32_t mod_u64_barrett(uint64_t x, uint32_t q, uint64_t qinv) {
    __uint128_t z = (__uint128_t)x * qinv;
    uint64_t t = (uint64_t)(z >> 64);
    uint64_t r = x - t * (uint64_t)q;
    if (r >= q) r -= q;
    if (r >= q) r -= q;
    return (uint32_t)r;
}

static inline uint32_t xorshift32(uint32_t *s) {
    uint32_t x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}

/* Falcon vrfy.c constants for q=12289 */
static inline uint32_t falcon_monty_mul_12289(uint32_t x, uint32_t y) {
    const uint32_t Q = 12289u;
    const uint32_t Q0I = 12287u; /* -1/q mod 2^16 */
    uint32_t z = x * y;
    uint32_t w = ((z * Q0I) & 0xFFFFu) * Q;
    z = (z + w) >> 16;
    z -= Q;
    z += Q & -(z >> 31);
    return z;
}

static inline uint32_t point_mul_ntt(uint32_t a, uint32_t b, uint32_t q, uint64_t qinv, int table_is_falcon_monty) {
    if (table_is_falcon_monty && q == 12289u) {
        return falcon_monty_mul_12289(a, b);
    }
    return mod_u64_barrett((uint64_t)a * (uint64_t)b, q, qinv);
}

static inline void accum_r_times_src_u64(uint64_t *restrict acc,
                                         const uint32_t *restrict src,
                                         uint32_t n, uint64_t r) {
#ifdef __AVX2__
    uint32_t j = 0;
    const __m256i vr = _mm256_set1_epi64x((long long)r);

    for (; j + 8 <= n; j += 8) {
        const __m256i vsrc32 = _mm256_loadu_si256((const __m256i *)(src + j));
        const __m128i lo128 = _mm256_castsi256_si128(vsrc32);
        const __m128i hi128 = _mm256_extracti128_si256(vsrc32, 1);

        const __m256i vsrc_lo64 = _mm256_cvtepu32_epi64(lo128);
        const __m256i vsrc_hi64 = _mm256_cvtepu32_epi64(hi128);

        __m256i vacc_lo = _mm256_loadu_si256((const __m256i *)(acc + j));
        __m256i vacc_hi = _mm256_loadu_si256((const __m256i *)(acc + j + 4));

        const __m256i vmul_lo = _mm256_mul_epu32(vsrc_lo64, vr);
        const __m256i vmul_hi = _mm256_mul_epu32(vsrc_hi64, vr);

        vacc_lo = _mm256_add_epi64(vacc_lo, vmul_lo);
        vacc_hi = _mm256_add_epi64(vacc_hi, vmul_hi);

        _mm256_storeu_si256((__m256i *)(acc + j), vacc_lo);
        _mm256_storeu_si256((__m256i *)(acc + j + 4), vacc_hi);
    }
    for (; j < n; ++j) acc[j] += r * (uint64_t)src[j];
#else
    for (uint32_t j = 0; j < n; ++j) acc[j] += r * (uint64_t)src[j];
#endif
}

static int grid_init(poly_grid_t *g, uint32_t rows, uint32_t cols, uint32_t n) {
    g->rows = rows;
    g->cols = cols;
    g->cell = (poly_t *)calloc((size_t)rows * cols, sizeof(poly_t));
    if (!g->cell) return -1;
    for (uint32_t i = 0; i < rows; ++i) {
        for (uint32_t j = 0; j < cols; ++j) {
            if (poly_init(&CELL(g, i, j), n) != 0) return -1;
        }
    }
    return 0;
}

static void grid_free(poly_grid_t *g) {
    if (!g || !g->cell) return;
    for (uint32_t i = 0; i < g->rows; ++i) {
        for (uint32_t j = 0; j < g->cols; ++j) {
            poly_free(&CELL(g, i, j));
        }
    }
    free(g->cell);
    memset(g, 0, sizeof(*g));
}

static int read_u32(FILE *fp, uint32_t *x) {
    return fscanf(fp, "%u", x) == 1 ? 0 : -1;
}

static int load_grid_file(const char *path, poly_grid_t *g, uint32_t expect_q, uint32_t expect_n) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    uint32_t q, n, rows, cols;
    if (read_u32(fp, &q) || read_u32(fp, &n) || read_u32(fp, &rows) || read_u32(fp, &cols)) { fclose(fp); return -1; }
    if (q != expect_q || n != expect_n) { fclose(fp); return -1; }

    if (grid_init(g, rows, cols, n) != 0) { fclose(fp); return -1; }

    for (uint32_t i = 0; i < rows; ++i) {
        for (uint32_t j = 0; j < cols; ++j) {
            for (uint32_t t = 0; t < n; ++t) {
                uint32_t v;
                if (read_u32(fp, &v)) { fclose(fp); return -1; }
                CELL(g, i, j).c[t] = v % q;
            }
        }
    }

    fclose(fp);
    return 0;
}

static void table_free(table_data_t *tb) {
    if (!tb) return;
    grid_free(&tb->gamma_ntt);
    grid_free(&tb->ATg_ntt);
    memset(tb, 0, sizeof(*tb));
}

static int load_table_file(const char *path, table_data_t *tb, uint32_t expect_q, uint32_t expect_n) {
    memset(tb, 0, sizeof(*tb));
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    uint32_t q, n, T_in;
    if (read_u32(fp, &q) || read_u32(fp, &n) ||
        read_u32(fp, &tb->R) || read_u32(fp, &T_in) || read_u32(fp, &tb->K) || read_u32(fp, &tb->L)) {
        fclose(fp); return -1;
    }
    tb->file_T = T_in;
    if (q != expect_q || n != expect_n) { fclose(fp); return -1; }

    if (grid_init(&tb->gamma_ntt, tb->R, tb->K, n) != 0) { fclose(fp); return -1; }
    if (grid_init(&tb->ATg_ntt, tb->R, tb->L, n) != 0) { fclose(fp); return -1; }

    for (uint32_t r = 0; r < tb->R; ++r) {
        for (uint32_t t = 0; t < T_in; ++t) {
            uint32_t dummy;
            if (read_u32(fp, &dummy)) { fclose(fp); return -1; } // ignore randv from file
        }
        for (uint32_t k = 0; k < tb->K; ++k) {
            for (uint32_t i = 0; i < n; ++i) {
                uint32_t v;
                if (read_u32(fp, &v)) { fclose(fp); return -1; }
                CELL(&tb->gamma_ntt, r, k).c[i] = v % q;
            }
        }
        for (uint32_t l = 0; l < tb->L; ++l) {
            for (uint32_t i = 0; i < n; ++i) {
                uint32_t v;
                if (read_u32(fp, &v)) { fclose(fp); return -1; }
                CELL(&tb->ATg_ntt, r, l).c[i] = v % q;
            }
        }
    }

    fclose(fp);
    return 0;
}

static void derive_random_vec(uint32_t *rv, uint32_t T, uint32_t q, uint32_t rid) {
    uint32_t s = 0x9E3779B9u ^ (rid + 1u);
    for (uint32_t i = 0; i < T; ++i) {
        rv[i] = xorshift32(&s) % q;
    }
}

static int brlt_ws_init(brlt_ws_t *ws, uint32_t T, uint32_t L, uint32_t K, uint32_t n, uint32_t q) {
    memset(ws, 0, sizeof(*ws));
    ws->T = T; ws->L = L; ws->K = K; ws->n = n; ws->q = q;
    ws->qinv = (uint64_t)((((__uint128_t)1) << 64) / q);

    ws->sum_z = (poly_t *)calloc(L, sizeof(poly_t));
    ws->sum_u = (poly_t *)calloc(K, sizeof(poly_t));
    ws->lazy_z = (uint64_t *)calloc((size_t)L * n, sizeof(uint64_t));
    ws->lazy_u = (uint64_t *)calloc((size_t)K * n, sizeof(uint64_t));

    if (!ws->sum_z || !ws->sum_u || !ws->lazy_z || !ws->lazy_u) return -1;
    for (uint32_t i = 0; i < L; ++i) if (poly_init(&ws->sum_z[i], n) != 0) return -1;
    for (uint32_t i = 0; i < K; ++i) if (poly_init(&ws->sum_u[i], n) != 0) return -1;
    return 0;
}

static void brlt_ws_free(brlt_ws_t *ws) {
    if (!ws) return;
    for (uint32_t i = 0; i < ws->L; ++i) poly_free(&ws->sum_z[i]);
    for (uint32_t i = 0; i < ws->K; ++i) poly_free(&ws->sum_u[i]);
    free(ws->sum_z); free(ws->sum_u);
    free(ws->lazy_z); free(ws->lazy_u);
    memset(ws, 0, sizeof(*ws));
}

static int bucket_ws_init(bucket_ws_t *bw, uint32_t T, uint32_t L, uint32_t K, uint32_t n) {
    memset(bw, 0, sizeof(*bw));
    bw->T = T; bw->L = L; bw->K = K; bw->n = n;
    bw->acc_c = (uint32_t *)calloc((size_t)T * L * n, sizeof(uint32_t));
    bw->acc_d = (uint32_t *)calloc((size_t)T * K * n, sizeof(uint32_t));
    return (bw->acc_c && bw->acc_d) ? 0 : -1;
}

static void bucket_ws_free(bucket_ws_t *bw) {
    if (!bw) return;
    free(bw->acc_c);
    free(bw->acc_d);
    memset(bw, 0, sizeof(*bw));
}

static void bucket_poly_lazy(const poly_grid_t *z, const poly_grid_t *u, poly_grid_t *c, poly_grid_t *d,
                             uint32_t T, uint32_t q, uint64_t qinv, uint32_t *rng_state, bucket_ws_t *bw) {
    const uint32_t L = c->cols, K = d->cols, n = bw->n;

    memset(bw->acc_c, 0, (size_t)T * L * n * sizeof(uint32_t));
    memset(bw->acc_d, 0, (size_t)T * K * n * sizeof(uint32_t));

    for (uint32_t i = 0; i < z->rows; ++i) {
        uint32_t b = xorshift32(rng_state) & (T - 1u);
        for (uint32_t l = 0; l < L; ++l) {
            const uint32_t *src = CELL(z, i, l).c;
            uint32_t *acc = &bw->acc_c[((size_t)b * L + l) * n];
            accum_src_u32_inplace(acc, src, n);
        }
        for (uint32_t k = 0; k < K; ++k) {
            const uint32_t *src = CELL(u, i, k).c;
            uint32_t *acc = &bw->acc_d[((size_t)b * K + k) * n];
            accum_src_u32_inplace(acc, src, n);
        }
    }

    for (uint32_t b = 0; b < T; ++b) {
        for (uint32_t l = 0; l < L; ++l) {
            uint32_t *dst = CELL(c, b, l).c;
            const uint32_t *acc = &bw->acc_c[((size_t)b * L + l) * n];
            for (uint32_t j = 0; j < n; ++j) dst[j] = mod_u64_barrett((uint64_t)acc[j], q, qinv);
        }
        for (uint32_t k = 0; k < K; ++k) {
            uint32_t *dst = CELL(d, b, k).c;
            const uint32_t *acc = &bw->acc_d[((size_t)b * K + k) * n];
            for (uint32_t j = 0; j < n; ++j) dst[j] = mod_u64_barrett((uint64_t)acc[j], q, qinv);
        }
    }
}

static int brlt_common_poly_ws(const poly_grid_t *c, const poly_grid_t *d,
                               const poly_t *gamma_ntt, const poly_t *AT_gamma_ntt,
                               const uint32_t *random_vec,
                               const ntt_ctx_t *ntt, uint32_t q, brlt_ws_t *ws,
                               int table_is_falcon_monty) {
    const uint32_t n = ws->n, T = ws->T, L = ws->L, K = ws->K;

    memset(ws->lazy_z, 0, (size_t)L * n * sizeof(uint64_t));
    memset(ws->lazy_u, 0, (size_t)K * n * sizeof(uint64_t));
    for (uint32_t i = 0; i < T; ++i) {
        uint64_t r = random_vec[i];
        for (uint32_t l = 0; l < L; ++l) {
            accum_r_times_src_u64(&ws->lazy_z[(size_t)l * n], CELL(c, i, l).c, n, r);
        }
        for (uint32_t k = 0; k < K; ++k) {
            accum_r_times_src_u64(&ws->lazy_u[(size_t)k * n], CELL(d, i, k).c, n, r);
        }
    }

    for (uint32_t l = 0; l < L; ++l) {
        for (uint32_t j = 0; j < n; ++j) {
            ws->sum_z[l].c[j] = mod_u64_barrett(ws->lazy_z[(size_t)l * n + j], q, ws->qinv);
        }
    }
    for (uint32_t k = 0; k < K; ++k) {
        for (uint32_t j = 0; j < n; ++j) {
            ws->sum_u[k].c[j] = mod_u64_barrett(ws->lazy_u[(size_t)k * n + j], q, ws->qinv);
        }
    }

    for (uint32_t l = 0; l < L; ++l) {
        if (ntt_forward_bitrev(ws->sum_z[l].c, ntt) != 0) return -1;
    }
    for (uint32_t k = 0; k < K; ++k) {
        if (ntt_forward_bitrev(ws->sum_u[k].c, ntt) != 0) return -1;
    }

    uint32_t diff = 0;
    for (uint32_t j = 0; j < n; ++j) {
        uint32_t az = 0, au = 0;
        for (uint32_t l = 0; l < L; ++l) {
            uint32_t p = point_mul_ntt(ws->sum_z[l].c[j], AT_gamma_ntt[l].c[j], q, ws->qinv, table_is_falcon_monty);
            az = add_mod_q(az, p, q);
        }
        for (uint32_t k = 0; k < K; ++k) {
            uint32_t p = point_mul_ntt(ws->sum_u[k].c[j], gamma_ntt[k].c[j], q, ws->qinv, table_is_falcon_monty);
            au = add_mod_q(au, p, q);
        }
        diff |= (az ^ au);
    }
    return diff == 0 ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <scheme> <mat_A_ntt.txt> <mat_z.txt> <mat_u.txt> <table.txt> [--falcon-monty-table]\n", argv[0]);
        return 1;
    }

    int table_is_falcon_monty = 0;
    if (argc >= 7 && strcmp(argv[6], "--falcon-monty-table") == 0) {
        table_is_falcon_monty = 1;
    }

    rlt_params_t P;
    if (params_from_name(argv[1], &P) != 0) return 1;

    ntt_ctx_t ntt;
    if (ntt_ctx_init(&ntt, P.n, P.q) != 0) return 1;

    poly_grid_t A_ntt = {0}, z = {0}, u = {0};
    if (load_grid_file(argv[2], &A_ntt, P.q, P.n) != 0) return 1;
    if (load_grid_file(argv[3], &z, P.q, P.n) != 0) return 1;
    if (load_grid_file(argv[4], &u, P.q, P.n) != 0) return 1;

    if (A_ntt.rows == 0 || A_ntt.cols == 0 || A_ntt.rows > 10 || A_ntt.cols > 10) return 1;
    if (z.cols != A_ntt.cols || u.cols != A_ntt.rows || z.rows != u.rows) return 1;

    uint32_t K = A_ntt.rows, L = A_ntt.cols, N_rows = z.rows;

    const uint32_t Lambda = 128u, m = 7u, T = 1u << m;
    const uint32_t bucket_times = (Lambda + m - 2u) / (m - 1u);
    const uint32_t brlt_times = (uint32_t)ceil((double)m / log2((double)P.q / 2.0));
    const uint32_t total_needed = bucket_times * brlt_times;

    table_data_t tb;
    if (load_table_file(argv[5], &tb, P.q, P.n) != 0) return 1;
    if (tb.R < total_needed || tb.K != K || tb.L != L) return 1;

    poly_grid_t c = {0}, d = {0};
    if (grid_init(&c, T, L, P.n) != 0) return 1;
    if (grid_init(&d, T, K, P.n) != 0) return 1;

    brlt_ws_t ws;
    if (brlt_ws_init(&ws, T, L, K, P.n, P.q) != 0) return 1;

    bucket_ws_t bws;
    if (bucket_ws_init(&bws, T, L, K, P.n) != 0) return 1;

    uint32_t rng_state = 1u;
    int ok_all = 1;

    uint32_t *rv = (uint32_t *)malloc((size_t)T * sizeof(uint32_t));
    if (!rv) return 1;

    /* compute window begin marker (outside timed region) */
    printf("[compute_begin]\n");
    fflush(stdout);

    uint64_t c0 = now_cycles();
    uint64_t t0 = now_ns();

    for (uint32_t bi = 0; bi < bucket_times && ok_all; ++bi) {
        bucket_poly_lazy(&z, &u, &c, &d, T, P.q, ws.qinv, &rng_state, &bws);

        for (uint32_t r = 0; r < brlt_times; ++r) {
            uint32_t rid = bi * brlt_times + r;
            derive_random_vec(rv, T, P.q, rid);
            int pass = brlt_common_poly_ws(
                &c, &d,
                &CELL(&tb.gamma_ntt, rid, 0),
                &CELL(&tb.ATg_ntt, rid, 0),
                rv, &ntt, P.q, &ws,
                table_is_falcon_monty
            );
            if (pass != 1) { ok_all = 0; break; }
        }
    }

    uint64_t t1 = now_ns();
    uint64_t c1 = now_cycles();

    /* compute window end marker (outside timed region) */
    printf("[compute_end]\n");
    fflush(stdout);

    double sec_total = (double)(t1 - t0) / 1e9;
    uint64_t cyc_total = (c1 >= c0) ? (c1 - c0) : 0;

    printf("scheme=%s q=%u n=%u N=%u K=%u L=%u T=%u\n", P.name, P.q, P.n, N_rows, K, L, T);
    printf("result=%s compute_time=%.6f s compute_cycles=%llu\n",
           ok_all ? "PASS" : "FAIL",
           sec_total,
           (unsigned long long)cyc_total);

    if (table_is_falcon_monty) {
        printf("[mode] falcon_monty_pointmul=ON (q must be 12289 and table coefficients must be NTT+Monty)\n");
    } else {
        printf("[mode] falcon_monty_pointmul=OFF\n");
    }

#if !HAVE_RDTSC
    printf("[warn] rdtsc not available on this architecture; compute_cycles=0\n");
#endif

    free(rv);
    bucket_ws_free(&bws);
    brlt_ws_free(&ws);
    grid_free(&c); grid_free(&d);
    table_free(&tb);
    grid_free(&A_ntt); grid_free(&z); grid_free(&u);
    ntt_ctx_free(&ntt);

    return ok_all ? 0 : 2;
}
