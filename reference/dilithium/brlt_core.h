#ifndef DILITHIUM_BRLT_CORE_H
#define DILITHIUM_BRLT_CORE_H

#include "io.h"

typedef struct {
    table_t *tb;
    uint32_t T;
} brlt_ctx_t;

int brlt_ctx_init(brlt_ctx_t *ctx, table_t *tb);
int brlt_verify(const brlt_ctx_t *ctx, const poly_grid_t *z, const poly_grid_t *u);

#endif