#ifndef AC26_PARAMS_H
#define AC26_PARAMS_H

#include <stdint.h>
#include <string.h>

typedef enum {
    SCHEME_FALCON512 = 0,
    SCHEME_FALCON1024,
    SCHEME_KYBER,
    SCHEME_DILITHIUM
} scheme_t;

typedef struct {
    scheme_t scheme;
    const char *name;
    uint32_t q;
    uint32_t n;        // NTT length: 256 or 512 (this project scope)
    uint32_t K_default;
    uint32_t L_default;
} rlt_params_t;

static inline int params_from_name(const char *s, rlt_params_t *p) {
    if (!s || !p) return -1;
    if (strcmp(s, "falcon512") == 0) {
        *p = (rlt_params_t){SCHEME_FALCON512, "falcon512", 12289u, 512u, 1u, 1u};
        return 0;
    }
    if (strcmp(s, "falcon1024") == 0) {
        // 按你的当前目标，仍用 512 长度核；后续可扩到 1024
        *p = (rlt_params_t){SCHEME_FALCON1024, "falcon1024", 12289u, 512u, 1u, 1u};
        return 0;
    }
    if (strcmp(s, "kyber") == 0) {
        *p = (rlt_params_t){SCHEME_KYBER, "kyber", 3329u, 256u, 2u, 2u};
        return 0;
    }
    if (strcmp(s, "dilithium") == 0) {
        *p = (rlt_params_t){SCHEME_DILITHIUM, "dilithium", 8380417u, 256u, 4u, 4u};
        return 0;
    }
    return -1;
}

#endif