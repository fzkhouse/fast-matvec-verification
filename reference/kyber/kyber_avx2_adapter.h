#ifndef KYBER_AVX2_ADAPTER_H
#define KYBER_AVX2_ADAPTER_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

#ifdef __cplusplus
extern "C" {
#endif

void kyb_poly_zero(poly *p);
void kyb_poly_from_u32(poly *dst, const uint32_t *src);
void kyb_poly_from_u64_modq(poly *dst, const uint64_t *src);
void kyb_poly_to_u32(uint32_t *dst, const poly *src);

void kyb_poly_ntt(poly *p);
void kyb_poly_reduce(poly *p);

void kyb_poly_basemul_acc(poly *acc, const poly *a_ntt, const poly *b_ntt);
int  kyb_poly_equal_modq(const poly *a, const poly *b);
int  kyb_poly_equal_reduced(const poly *a, const poly *b);

#ifdef __cplusplus
}
#endif

#endif