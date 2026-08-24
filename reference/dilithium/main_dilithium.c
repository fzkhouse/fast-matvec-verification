#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <x86intrin.h>  // for __rdtsc
#include "brlt_core.h"

static inline uint64_t now_ns(void) {
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
    if (argc != 4) {
        fprintf(stderr, "usage: %s <mat_z.txt> <mat_u.txt> <table.txt>\n", argv[0]);
        return 1;
    }

    poly_grid_t z = {0}, u = {0};
    table_t tb = {0};
    brlt_ctx_t ctx = {0};

    if (grid_load(argv[1], &z, DILITHIUM_Q, DILITHIUM_N) != 0) return 1;
    if (grid_load(argv[2], &u, DILITHIUM_Q, DILITHIUM_N) != 0) return 1;
    if (table_load(argv[3], &tb, DILITHIUM_Q, DILITHIUM_N) != 0) return 1;

    if (brlt_ctx_init(&ctx, &tb) != 0) return 1;

    printf("[compute_begin]\n");
    uint64_t t0 = now_ns();
    uint64_t cycles_start = __rdtsc();
    uint64_t rss_start = current_rss_kb();
    int ok = brlt_verify(&ctx, &z, &u);
    uint64_t rss_end = current_rss_kb();
    uint64_t cycles_end = __rdtsc();
    uint64_t t1 = now_ns();
    printf("[compute_end]\n");

    uint64_t cycles = cycles_end - cycles_start;
    double avg_rss_kb = (double)(rss_start + rss_end) / 2.0;
    uint64_t peak_rss_kb = (rss_start > rss_end) ? rss_start : rss_end;

    printf("scheme=dilithium%d q=%u n=%u Nrows=%u K=%u L=%u\n",
           DILITHIUM_MODE, DILITHIUM_Q, DILITHIUM_N, z.rows, DILITHIUM_K, DILITHIUM_L);
    printf("result=%s time=%.6f s cycles=%llu avg_rss_kb=%.0f peak_rss_kb=%llu\n",
           ok ? "PASS" : "FAIL", (double)(t1 - t0) / 1e9,
           (unsigned long long)cycles, avg_rss_kb,
           (unsigned long long)peak_rss_kb);

    table_free(&tb);
    grid_free(&z);
    grid_free(&u);
    return ok ? 0 : 2;
}
