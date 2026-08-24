#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "inner.h"

#define Q 12289u

#if defined(__x86_64__) || defined(__i386__)
#define HAVE_RDTSC 1
#else
#define HAVE_RDTSC 0
#endif

typedef struct {
    uint32_t q, n, rows, cols;
    uint16_t *data;
} poly_grid_t;

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

static int read_u32(FILE *fp, uint32_t *x) {
    return (fscanf(fp, "%u", x) == 1) ? 0 : -1;
}

static int load_grid_file(const char *path, poly_grid_t *g) {
    size_t i, total;
    FILE *fp;

    memset(g, 0, sizeof(*g));
    fp = fopen(path, "r");
    if (!fp) return -1;

    if (read_u32(fp, &g->q) || read_u32(fp, &g->n)
        || read_u32(fp, &g->rows) || read_u32(fp, &g->cols))
    {
        fclose(fp);
        return -1;
    }

    total = (size_t)g->rows * g->cols * g->n;
    g->data = (uint16_t *)malloc(total * sizeof(uint16_t));
    if (!g->data) {
        fclose(fp);
        return -1;
    }

    for (i = 0; i < total; i++) {
        uint32_t v;
        if (read_u32(fp, &v)) {
            fclose(fp);
            return -1;
        }
        g->data[i] = (uint16_t)(v % g->q);
    }

    fclose(fp);
    return 0;
}

static void free_grid(poly_grid_t *g) {
    if (!g) return;
    free(g->data);
    memset(g, 0, sizeof(*g));
}

static int ilog2_pow2(uint32_t n) {
    int lg = 0;
    if (n == 0 || (n & (n - 1u)) != 0) return -1;
    while ((1u << lg) != n) lg++;
    return lg;
}

int main(int argc, char **argv) {
    poly_grid_t Acoef = {0}, Z = {0}, U = {0};
    uint16_t *A_nttm = NULL;
    size_t apoly_cnt, p, an;
    int logn, ok;
    uint64_t t0, t1, c0, c1, cyc_total;

    if (argc != 4) {
        fprintf(stderr, "usage: %s <mat_A_coeff.txt> <mat_z.txt> <mat_u.txt>\n", argv[0]);
        return 1;
    }

    if (load_grid_file(argv[1], &Acoef) != 0) return 1;
    if (load_grid_file(argv[2], &Z) != 0) return 1;
    if (load_grid_file(argv[3], &U) != 0) return 1;

    if (Acoef.q != Q || Z.q != Q || U.q != Q) return 1;
    if (Acoef.n != Z.n || Acoef.n != U.n) return 1;
    if (Z.cols != Acoef.cols) return 1;
    if (U.cols != Acoef.rows) return 1;
    if (Z.rows != U.rows) return 1;

    logn = ilog2_pow2(Acoef.n);
    if (logn < 0) return 1;

    /* 计时外：A_coeff -> Falcon NTT+Monty */
    apoly_cnt = (size_t)Acoef.rows * Acoef.cols;
    an = (size_t)Acoef.n;
    A_nttm = (uint16_t *)malloc(apoly_cnt * an * sizeof(uint16_t));
    if (!A_nttm) return 1;

    memcpy(A_nttm, Acoef.data, apoly_cnt * an * sizeof(uint16_t));
    for (p = 0; p < apoly_cnt; p++) {
        Zf(to_ntt_monty)(A_nttm + p * an, (unsigned)logn);
    }

    /* 只计总验算区间 */
    c0 = now_cycles();
    t0 = now_ns();
    ok = Zf(matvec_verify_eq_same_input)(
        A_nttm, Z.data, U.data,
        Acoef.rows, Acoef.cols, Z.rows,
        (unsigned)logn);
    t1 = now_ns();
    c1 = now_cycles();

    cyc_total = (c1 >= c0) ? (c1 - c0) : 0;

    printf("Falcon-native matvec-eq: %s\n", ok ? "PASS" : "FAIL");
    printf("N=%u K=%u L=%u n=%u\n", Z.rows, Acoef.rows, Acoef.cols, Acoef.n);
    printf("time_matvec_eq=%.6f s cycles_matvec_eq=%llu\n",
           (double)(t1 - t0) / 1e9,
           (unsigned long long)cyc_total);

#if !HAVE_RDTSC
    printf("[warn] rdtsc not available on this architecture; cycles_matvec_eq=0\n");
#endif

    free(A_nttm);
    free_grid(&Acoef);
    free_grid(&Z);
    free_grid(&U);
    return ok ? 0 : 2;
}