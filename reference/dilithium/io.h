#ifndef DILITHIUM_BRLT_IO_H
#define DILITHIUM_BRLT_IO_H

#include "poly.h"
#include <stdint.h>

typedef struct {
    uint32_t rows, cols;
    poly_t *data;
} poly_grid_t;

typedef struct {
    uint32_t q, n, R, T, K_, L_;
    poly_t *A;
    poly_t *gamma;
    poly_t *ATgamma;
    uint32_t *randoms;
} table_t;

int grid_save(const char *path, const poly_grid_t *g);
int grid_load(const char *path, poly_grid_t *g, uint32_t expect_q, uint32_t expect_n);
void grid_free(poly_grid_t *g);

int table_save(const char *path, const table_t *tb);
int table_load(const char *path, table_t *tb, uint32_t expect_q, uint32_t expect_n);
void table_free(table_t *tb);

#endif