#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "params.h"
#include "poly.h"

#define KYB_RINV 169

typedef struct {
    int q, n, rows, cols;
    poly *data; // rows * cols
} poly_grid_t;

#if defined(__x86_64__) || defined(__i386__)
static inline uint64_t now_cycles(void) {
    uint32_t hi, lo;
    __asm__ __volatile__("lfence\nrdtsc" : "=a"(lo), "=d"(hi) :: "memory");
    return ((uint64_t)hi << 32) | lo;
}
#else
static inline uint64_t now_cycles(void) { return 0; }
#endif

static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

static int mod_q(int x) {
    int v = x % KYBER_Q;
    return (v < 0) ? (v + KYBER_Q) : v;
}

static int16_t mul_mod_q_i16(int16_t a, int16_t b) {
    int32_t x = (int32_t)a * (int32_t)b;
    int32_t v = x % KYBER_Q;
    if (v < 0) v += KYBER_Q;
    return (int16_t)v;
}

static void poly_scale_rinv(poly *p) {
    for (int i = 0; i < KYBER_N; i++) {
        p->coeffs[i] = mul_mod_q_i16(p->coeffs[i], (int16_t)KYB_RINV);
    }
}

static void poly_zero_local(poly *p) {
    memset(p, 0, sizeof(*p));
}

static int read_poly_line(FILE *fp, poly *p) {
    for (int i = 0; i < KYBER_N; i++) {
        int v;
        if (fscanf(fp, "%d", &v) != 1) return -1;
        p->coeffs[i] = (int16_t)mod_q(v);
    }
    return 0;
}

static int load_grid(const char *path, poly_grid_t *g) {
    memset(g, 0, sizeof(*g));
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    if (fscanf(fp, "%d %d %d %d", &g->q, &g->n, &g->rows, &g->cols) != 4) {
        fclose(fp);
        return -1;
    }
    if (g->q != KYBER_Q || g->n != KYBER_N || g->rows <= 0 || g->cols <= 0) {
        fclose(fp);
        return -1;
    }

    size_t cnt = (size_t)g->rows * (size_t)g->cols;

    if (posix_memalign((void **)&g->data, 32, cnt * sizeof(poly)) != 0 || !g->data) {
        fclose(fp);
        g->data = NULL;
        return -1;
    }
    memset(g->data, 0, cnt * sizeof(poly));

    for (size_t i = 0; i < cnt; i++) {
        if (read_poly_line(fp, &g->data[i]) != 0) {
            free(g->data);
            g->data = NULL;
            fclose(fp);
            return -1;
        }
    }

    fclose(fp);
    return 0;
}

static void free_grid(poly_grid_t *g) {
    free(g->data);
    memset(g, 0, sizeof(*g));
}

static int poly_equal_mod_q(const poly *a, const poly *b) {
    for (int i = 0; i < KYBER_N; i++) {
        if (mod_q(a->coeffs[i]) != mod_q(b->coeffs[i])) return 0;
    }
    return 1;
}

static int verify_once(const poly *A_ntt, const poly *z, const poly *u, int rows, int K, int L) {
    poly *z_ntt = NULL;
    if (posix_memalign((void **)&z_ntt, 32, (size_t)L * sizeof(poly)) != 0) {
        return 0;
    }

    for (int r = 0; r < rows; r++) {
        for (int l = 0; l < L; l++) {
            z_ntt[l] = z[r * L + l];
            poly_ntt(&z_ntt[l]);
            poly_reduce(&z_ntt[l]);
        }

        for (int k = 0; k < K; k++) {
            poly acc;
            poly_zero_local(&acc);

            for (int l = 0; l < L; l++) {
                poly tmp;
                const poly *Akl = &A_ntt[k * L + l];
                poly_basemul_montgomery(&tmp, Akl, &z_ntt[l]);
                poly_add(&acc, &acc, &tmp);
                poly_reduce(&acc);
            }

            poly_invntt_tomont(&acc);
            poly_scale_rinv(&acc);
            poly_reduce(&acc);

            if (!poly_equal_mod_q(&acc, &u[r * K + k])) {
                free(z_ntt);
                return 0;
            }
        }
    }

    free(z_ntt);
    return 1;
}

static void usage(const char *p) {
    fprintf(stderr,
        "Usage:\n"
        "  %s <mat_A_ntt.txt> <mat_z.txt> <mat_u.txt>\n"
        "  %s <mat_A_ntt.txt> <mat_z.txt> <mat_u.txt> <warmup> <repeat>\n"
        "  %s <mat_A_ntt.txt> <mat_z.txt> <mat_u.txt> <table_ignored> <warmup> <repeat>\n",
        p, p, p);
}

int main(int argc, char **argv) {
    if (!(argc == 4 || argc == 6 || argc == 7)) {
        usage(argv[0]);
        return 1;
    }

    const char *A_path = argv[1];
    const char *z_path = argv[2];
    const char *u_path = argv[3];

    int warmup = 2;
    int repeat = 9;

    if (argc == 6) {
        warmup = atoi(argv[4]);
        repeat = atoi(argv[5]);
    } else if (argc == 7) {
        warmup = atoi(argv[5]);
        repeat = atoi(argv[6]);
    }

    if (warmup < 0 || repeat <= 0) {
        fprintf(stderr, "invalid warmup/repeat\n");
        return 1;
    }

    poly_grid_t A = {0}, Z = {0}, U = {0};
    if (load_grid(A_path, &A) != 0 || load_grid(z_path, &Z) != 0 || load_grid(u_path, &U) != 0) {
        fprintf(stderr, "load failed\n");
        free_grid(&A); free_grid(&Z); free_grid(&U);
        return 1;
    }

    int K = A.rows;
    int L = A.cols;

    if (!(Z.cols == L && U.cols == K && Z.rows == U.rows)) {
        fprintf(stderr, "z/u shape mismatch\n");
        free_grid(&A); free_grid(&Z); free_grid(&U);
        return 1;
    }

    for (int i = 0; i < warmup; i++) {
        if (!verify_once(A.data, Z.data, U.data, Z.rows, K, L)) {
            fprintf(stderr, "warmup verify failed\n");
            free_grid(&A); free_grid(&Z); free_grid(&U);
            return 2;
        }
    }

    printf("[compute_begin]\n");
    fflush(stdout);

    uint64_t t0 = now_ns();
    uint64_t c0 = now_cycles();
    int ok = 1;
    for (int i = 0; i < repeat; i++) {
        if (!verify_once(A.data, Z.data, U.data, Z.rows, K, L)) {
            ok = 0;
            break;
        }
    }
    uint64_t t1 = now_ns();
    uint64_t c1 = now_cycles();

    double total_s = (double)(t1 - t0) / 1e9;
    double avg_s = total_s / (double)repeat;
    uint64_t total_cycles = (c1 >= c0) ? (c1 - c0) : 0ull;
    uint64_t avg_cycles = repeat ? (total_cycles / repeat) : 0ull;

    printf("[compute_end]\n");
    printf("scheme=kyber_native_mv rows=%d K=%d L=%d warmup=%d repeat=%d\n",
           Z.rows, K, L, warmup, repeat);
    printf("result=%s compute_time=%.6f s compute_total_time=%.6f s compute_cycles=%llu avg_cycles=%llu\n",
           ok ? "PASS" : "FAIL", avg_s, total_s,
           (unsigned long long)total_cycles, (unsigned long long)avg_cycles);

    free_grid(&A); free_grid(&Z); free_grid(&U);
    return ok ? 0 : 2;
}