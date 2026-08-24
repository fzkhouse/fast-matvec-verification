#ifndef DILITHIUM_BRLT_TABLE_H
#define DILITHIUM_BRLT_TABLE_H

#include "io.h"

int table_generate(const char *path, uint32_t R, uint32_t T, const poly_t *A);

#endif