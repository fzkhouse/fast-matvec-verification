// main_baseline.c - Baseline逐行验证：使用Dilithium NTT计算A*z=u
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "io.h"
#include <x86intrin.h>

/* 获取当前时间（纳秒） */
static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

/* 将系数域多项式复制到输出（z和u本身就在系数域） */
static inline void copy_poly(poly_t *dst, const poly_t *src) {
    for (int i = 0; i < DILITHIUM_N; i++) {
        dst->coeffs[i] = src->coeffs[i];
    }
}

static void* aligned_malloc(size_t size) {
    void *p = NULL;
    if (posix_memalign(&p, 32, size) != 0) return NULL;
    return p;
}

/* 使用Dilithium NTT逐行验证 A*z = u
 *
 * z: 系数域向量 [rows][L]
 * u: 系数域向量 [rows][K]
 * A: 系数域矩阵 [K][L]
 *
 * 验证方式：对每行i，计算 A*z_i，并与u_i比较
 */
int baseline_verify_ntt(const poly_t *A_ntt, const poly_grid_t *z, const poly_grid_t *u) {
    poly_t z_ntt[DILITHIUM_L];   /* z的NTT形式 */
    poly_t Az_ntt[DILITHIUM_K];  /* A*z的NTT形式 */
    poly_t Az[DILITHIUM_K];      /* 转回系数域后的A*z */

    /* 逐行验证每个(z_i, u_i)对 */
    for (uint32_t i = 0; i < z->rows; i++) {
        /* 1. 对z_i的L个多项式做NTT */
        for (uint32_t l = 0; l < DILITHIUM_L; l++) {
            copy_poly(&z_ntt[l], &z->data[i * DILITHIUM_L + l]);
            poly_ntt(&z_ntt[l]);
        }

        /* 2. 在NTT域计算 Az = A * z_i */
        for (uint32_t k = 0; k < DILITHIUM_K; k++) {
            poly_zero(&Az_ntt[k]);
            /* Az[k] = sum_l A_ntt[k][l] * z_ntt[l] */
            for (uint32_t l = 0; l < DILITHIUM_L; l++) {
                poly_pointwise_accum(&Az_ntt[k], &A_ntt[k * DILITHIUM_L + l], &z_ntt[l]);
            }
        }

        /* 3. 逆NTT转回系数域 */
        for (uint32_t k = 0; k < DILITHIUM_K; k++) {
            poly_invntt_tomont(&Az_ntt[k]);
            poly_from_mont(&Az_ntt[k]);
            Az[k] = Az_ntt[k];
        }

        /* 4. 与u_i比较 */
        for (uint32_t k = 0; k < DILITHIUM_K; k++) {
            poly_t *uk = &u->data[i * DILITHIUM_K + k];
            if (!poly_equal_modq(&Az[k], uk)) {
                fprintf(stderr, "Row %u, poly %u: verification failed\n", i, k);
                return 0;
            }
        }
    }

    return 1;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s <mat_z.txt> <mat_u.txt> <mat_A.txt> <table.txt>\n", argv[0]);
        fprintf(stderr, "\n  Input files:\n");
        fprintf(stderr, "    mat_z.txt   - z vector (coefficient domain), rows x L\n");
        fprintf(stderr, "    mat_u.txt   - u vector (coefficient domain), rows x K\n");
        fprintf(stderr, "    mat_A.txt   - A matrix (coefficient domain), K x L\n");
        fprintf(stderr, "    table.txt   - BRLT table (used for reference, not for baseline)\n");
        return 1;
    }

    const char *z_path = argv[1];
    const char *u_path = argv[2];
    const char *A_path = argv[3];
    const char *tb_path = argv[4];

    poly_grid_t z = {0}, u = {0}, A = {0};
    table_t tb = {0};

    /* 加载数据（不计入计算时间） */
    printf("Loading data...\n");
    if (grid_load(z_path, &z, DILITHIUM_Q, DILITHIUM_N) != 0) {
        fprintf(stderr, "Failed to load %s\n", z_path);
        return 1;
    }
    if (grid_load(u_path, &u, DILITHIUM_Q, DILITHIUM_N) != 0) {
        fprintf(stderr, "Failed to load %s\n", u_path);
        return 1;
    }
    if (grid_load(A_path, &A, DILITHIUM_Q, DILITHIUM_N) != 0) {
        fprintf(stderr, "Failed to load %s\n", A_path);
        return 1;
    }
    if (table_load(tb_path, &tb, DILITHIUM_Q, DILITHIUM_N) != 0) {
        fprintf(stderr, "Failed to load %s\n", tb_path);
        return 1;
    }

    /* 打印参数信息 */
    printf("scheme=dilithium%d q=%u n=%u Nrows=%u K=%u L=%u\n",
           DILITHIUM_MODE, DILITHIUM_Q, DILITHIUM_N, z.rows, DILITHIUM_K, DILITHIUM_L);

    /* 预计算A的NTT形式（一次性） */
    printf("Precomputing A NTT form...\n");
    poly_t *A_ntt = (poly_t*)aligned_malloc(DILITHIUM_K * DILITHIUM_L * sizeof(poly_t));
    if (!A_ntt) {
        fprintf(stderr, "Failed to allocate A_ntt\n");
        return 1;
    }
    for (uint32_t k = 0; k < DILITHIUM_K; k++) {
        for (uint32_t l = 0; l < DILITHIUM_L; l++) {
            copy_poly(&A_ntt[k * DILITHIUM_L + l], &A.data[k * DILITHIUM_L + l]);
            poly_ntt(&A_ntt[k * DILITHIUM_L + l]);
        }
    }

    /* 运行逐行NTT验证（记录计算时间与CPU周期） */
    printf("\n[Baseline NTT verification - One-by-one A*z=u]\n");
    printf("[compute_begin]\n");
    struct timespec tstart, tend;
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    uint64_t cycles_start = __rdtsc();

    int ok = baseline_verify_ntt(A_ntt, &z, &u);

    uint64_t cycles_end = __rdtsc();
    clock_gettime(CLOCK_MONOTONIC, &tend);
    printf("[compute_end]\n");

    double elapsed_s = (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9;
    uint64_t cycles = cycles_end - cycles_start;

    printf("result=%s\n", ok ? "PASS" : "FAIL");
    printf("time=%.6f s\n", elapsed_s);
    printf("cycles=%llu\n", (unsigned long long)cycles);

    /* 释放内存 */
    free(A_ntt);
    table_free(&tb);
    grid_free(&A);
    grid_free(&z);
    grid_free(&u);

    return ok ? 0 : 1;
}
