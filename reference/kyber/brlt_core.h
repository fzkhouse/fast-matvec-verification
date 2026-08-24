#ifndef BRLT_CORE_H
#define BRLT_CORE_H

#include <stddef.h>
#include <stdint.h>
#include "kyber_avx2_adapter.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t rows, cols, n;
    uint32_t *data; // [rows][cols][n]
} poly_grid_u32_t;

typedef struct {
    uint32_t rows, cols, n;
    uint16_t *data; // [rows][cols][n], canonical coefficients
} poly_grid_u16_t;

typedef struct {
    uint32_t R, K, L;
    uint32_t file_T;
    poly_grid_u32_t gamma_ntt; // [R][K][n]
    poly_grid_u32_t ATg_ntt;   // [R][L][n]
} table_u32_t;

typedef struct {
    uint32_t T, n, K, L, R;

    // workspace
    uint32_t *acc_c;  // [T][L][n]
    uint32_t *acc_d;  // [T][K][n]
    uint64_t *lazy_z; // [L][n]
    uint64_t *lazy_u; // [K][n]
    uint32_t *rv;     // [T]
    uint16_t *group;  // [4][T][L+K][n], zeroed once per fused round group
    size_t group_words;

    // decoded table polys
    poly *gamma_poly; // [R][K]
    poly *ATg_poly;   // [R][L]
} brlt_ctx_t;

void grid_free_u32(poly_grid_u32_t *g);
int  grid_load_u32(const char *path, poly_grid_u32_t *g, uint32_t expect_q, uint32_t expect_n);
void grid_free_u16(poly_grid_u16_t *g);
int  grid_load_u16(const char *path, poly_grid_u16_t *g, uint32_t expect_q, uint32_t expect_n);

void table_free_u32(table_u32_t *tb);
int  table_load_u32(const char *path, table_u32_t *tb, uint32_t expect_q, uint32_t expect_n);

int  brlt_ctx_init_kyber512(brlt_ctx_t *ctx, const table_u32_t *tb, uint32_t T);
void brlt_ctx_free(brlt_ctx_t *ctx);

int  brlt_verify_kyber512_ctx(
    brlt_ctx_t *ctx,
    const poly_grid_u16_t *z,
    const poly_grid_u16_t *u,
    uint32_t bucket_times,
    uint32_t brlt_times
);

#ifdef __cplusplus
}
#endif

#endif
