#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "brlt_core.h"

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

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s kyber512 <mat_A_ntt.txt> <mat_z.txt> <mat_u.txt> <table.txt>\n", argv[0]);
        return 1;
    }

    const uint32_t Q = 3329, N = 256;
    const uint32_t Lambda = 128u, m = 8u, T = 1u << m;
    const uint32_t bucket_times = (Lambda + m - 2u) / (m - 1u);
    const uint32_t brlt_times = (uint32_t)ceil((double)m / log2((double)Q / 2.0));

    poly_grid_u32_t A_ntt = {0};
    poly_grid_u16_t z = {0}, u = {0};
    table_u32_t tb = {0};
    brlt_ctx_t ctx = {0};

    if (grid_load_u32(argv[2], &A_ntt, Q, N) != 0) return 1;
    /* z/u are canonical q-sized coefficients.  Keep this compact form so
       the timed verifier does not repeatedly narrow 32-bit text inputs. */
    if (grid_load_u16(argv[3], &z, Q, N) != 0) return 1;
    if (grid_load_u16(argv[4], &u, Q, N) != 0) return 1;
    if (table_load_u32(argv[5], &tb, Q, N) != 0) return 1;

    if (z.cols != A_ntt.cols || u.cols != A_ntt.rows || z.rows != u.rows) return 1;
    if (tb.K != A_ntt.rows || tb.L != A_ntt.cols || tb.R < bucket_times * brlt_times) return 1;

    // 预分配与预解码放到计时外
    if (brlt_ctx_init_kyber512(&ctx, &tb, T) != 0) return 1;

    printf("[compute_begin]\n");
    fflush(stdout);

    uint64_t c0 = now_cycles();
    uint64_t t0 = now_ns();

    int ok = brlt_verify_kyber512_ctx(&ctx, &z, &u, bucket_times, brlt_times);

    uint64_t t1 = now_ns();
    uint64_t c1 = now_cycles();

    printf("[compute_end]\n");
    printf("scheme=%s q=%u n=%u Nrows=%u K=%u L=%u T=%u\n", argv[1], Q, N, z.rows, tb.K, tb.L, T);
    printf("result=%s compute_time=%.6f s compute_cycles=%llu\n",
           (ok == 1 ? "PASS" : (ok == 0 ? "FAIL" : "ERROR")),
           (double)(t1 - t0) / 1e9,
           (unsigned long long)((c1 >= c0) ? (c1 - c0) : 0ull));

    brlt_ctx_free(&ctx);
    table_free_u32(&tb);
    grid_free_u32(&A_ntt);
    grid_free_u16(&z);
    grid_free_u16(&u);

    return (ok == 1) ? 0 : 2;
}
