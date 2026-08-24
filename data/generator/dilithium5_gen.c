#define _POSIX_C_SOURCE 200112L
#include "table.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

static size_t total_allocated = 0;
static size_t peak_allocated = 0;

static void* aligned_malloc(size_t size) {
    void *p = NULL;
    if (posix_memalign(&p, 32, size) != 0) return NULL;
    total_allocated += size;
    if (total_allocated > peak_allocated) peak_allocated = total_allocated;
    return p;
}

static void aligned_free(void *p, size_t size) {
    free(p);
    total_allocated -= size;
}

static uint32_t rng32_next(uint32_t *s) {
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}

static void poly_rand(poly_t *p, uint32_t *s) {
    for (int i = 0; i < DILITHIUM_N; i++) p->coeffs[i] = (int32_t)(rng32_next(s) % DILITHIUM_Q);
}

static uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

static uint64_t current_rss_kb(void) {
    FILE *fp = fopen("/proc/self/status", "r");
    if (!fp) return 0;

    char line[256];
    uint64_t rss = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            unsigned long long tmp = 0;
            if (sscanf(line + 6, "%llu", &tmp) == 1) {
                rss = (uint64_t)tmp;
            }
            break;
        }
    }

    fclose(fp);
    return rss;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s <rows> <R> <T> <out_dir>\n", argv[0]);
        return 1;
    }
    uint32_t rows = (uint32_t)atoi(argv[1]);
    uint32_t R = (uint32_t)atoi(argv[2]);
    uint32_t T = (uint32_t)atoi(argv[3]);
    const char *out = argv[4];

    poly_grid_t A = { .rows = DILITHIUM_K, .cols = DILITHIUM_L };
    poly_grid_t z = { .rows = rows, .cols = DILITHIUM_L };
    poly_grid_t u = { .rows = rows, .cols = DILITHIUM_K };

    size_t A_size = (size_t)A.rows * A.cols * sizeof(poly_t);
    size_t z_size = (size_t)z.rows * z.cols * sizeof(poly_t);
    size_t u_size = (size_t)u.rows * u.cols * sizeof(poly_t);

    A.data = (poly_t*)aligned_malloc(A_size);
    if (!A.data) return 1;

    uint32_t s = 12345;
    for (uint32_t i = 0; i < A.rows * A.cols; i++) poly_rand(&A.data[i], &s);

    /* 预计算 A 的 NTT */
    size_t A_ntt_size = (size_t)A.rows * A.cols * sizeof(poly_t);
    poly_t *A_ntt = (poly_t*)aligned_malloc(A_ntt_size);
    for (uint32_t i = 0; i < A.rows * A.cols; i++) {
        A_ntt[i] = A.data[i];
        poly_ntt(&A_ntt[i]);
    }

    char path_A[512], path_z[512], path_u[512], path_tb[512];
    snprintf(path_A, sizeof(path_A), "%s/mat_A_%04u.txt", out, rows);
    snprintf(path_z, sizeof(path_z), "%s/mat_z_%04u.txt", out, rows);
    snprintf(path_u, sizeof(path_u), "%s/mat_u_%04u.txt", out, rows);
    snprintf(path_tb, sizeof(path_tb), "%s/table_%04u.txt", out, rows);

    if (grid_save(path_A, &A) != 0) return 1;

    printf("[table_begin]\n");
    fflush(stdout);

    uint64_t table_t0 = now_ns();
    uint64_t table_rss_sum = 0;
    uint64_t table_rss_cnt = 0;
    uint64_t table_rss_peak = 0;

    uint64_t rss = current_rss_kb();
    table_rss_sum += rss;
    table_rss_cnt += 1;
    if (rss > table_rss_peak) table_rss_peak = rss;

    if (table_generate(path_tb, R, T, A.data) != 0) return 1;

    rss = current_rss_kb();
    table_rss_sum += rss;
    table_rss_cnt += 1;
    if (rss > table_rss_peak) table_rss_peak = rss;

    uint64_t table_elapsed_ns = now_ns() - table_t0;
    double table_avg_rss_kb = (table_rss_cnt == 0) ? 0.0 : (double)table_rss_sum / (double)table_rss_cnt;

    printf("table_time=%.6f s table_avg_rss_kb=%.0f table_peak_rss_kb=%llu\n",
           (double)table_elapsed_ns / 1e9,
           table_avg_rss_kb,
           (unsigned long long)table_rss_peak);
    printf("[table_end]\n");
    fflush(stdout);

    // 现在生成 z 和 u
    z.data = (poly_t*)aligned_malloc(z_size);
    u.data = (poly_t*)aligned_malloc(u_size);
    if (!z.data || !u.data) return 1;

    for (uint32_t i = 0; i < z.rows * z.cols; i++) poly_rand(&z.data[i], &s);

    /* 用 NTT 累加生成 u（与验证路径一致） */
    size_t z_ntt_size = (size_t)z.cols * sizeof(poly_t);
    poly_t *z_ntt = (poly_t*)aligned_malloc(z_ntt_size);
    for (uint32_t i = 0; i < u.rows; i++) {
        for (uint32_t l = 0; l < DILITHIUM_L; l++) {
            z_ntt[l] = z.data[i * DILITHIUM_L + l];
            poly_ntt(&z_ntt[l]);
        }

        for (uint32_t k = 0; k < DILITHIUM_K; k++) {
            poly_t acc_ntt; poly_zero(&acc_ntt);
            for (uint32_t l = 0; l < DILITHIUM_L; l++) {
                const poly_t *Akl = &A_ntt[k * DILITHIUM_L + l];
                poly_pointwise_accum(&acc_ntt, Akl, &z_ntt[l]);
            }
            poly_invntt_tomont(&acc_ntt);
            poly_from_mont(&acc_ntt);
            u.data[i * DILITHIUM_K + k] = acc_ntt;
        }
    }
    aligned_free(z_ntt, z_ntt_size);

    if (grid_save(path_z, &z) != 0) return 1;
    if (grid_save(path_u, &u) != 0) return 1;

    printf("generated:\n%s\n%s\n%s\n%s\n", path_A, path_z, path_u, path_tb);

    aligned_free(A.data, A_size);
    aligned_free(z.data, z_size);
    aligned_free(u.data, u_size);
    aligned_free(A_ntt, A_ntt_size);
    return 0;
}
