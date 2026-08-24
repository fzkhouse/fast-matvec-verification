#include "ntt.h"
#include <stdlib.h>
#include <string.h>

/* ---------- forward declarations ---------- */
static uint32_t mod_add(uint32_t a, uint32_t b, uint32_t q);
static uint32_t mod_sub(uint32_t a, uint32_t b, uint32_t q);
static uint32_t mod_pow(uint32_t a, uint32_t e, uint32_t q);
static uint32_t mod_inv(uint32_t a, uint32_t q);
static uint32_t gcd_u32(uint32_t a, uint32_t b);
static int factor_u32(uint32_t x, uint32_t *f, int maxf);
static int find_primitive_root(uint32_t q, uint32_t *g_out);
static void bit_reverse(uint32_t *a, uint32_t n);
static inline uint32_t mod_u64_barrett(uint64_t x, uint32_t q, uint64_t qinv);
static inline uint32_t mod_mul_fast(uint32_t a, uint32_t b, uint32_t q, uint64_t qinv);
static int build_twiddle_tables(ntt_ctx_t *ctx);
static void ntt_cyclic(uint32_t *a, const ntt_ctx_t *ctx, int inverse);
static void ntt_cyclic_fwd_dif_bitrev(uint32_t *a, const ntt_ctx_t *ctx);

/* ---------- basic modular ops ---------- */
static uint32_t mod_add(uint32_t a, uint32_t b, uint32_t q) {
    uint64_t t = (uint64_t)a + b;
    if (t >= q) t -= q;
    return (uint32_t)t;
}

static uint32_t mod_sub(uint32_t a, uint32_t b, uint32_t q) {
    return (a >= b) ? (a - b) : (a + q - b);
}

static uint32_t mod_pow(uint32_t a, uint32_t e, uint32_t q) {
    uint64_t r = 1, x = a % q;
    while (e) {
        if (e & 1u) r = (r * x) % q;
        x = (x * x) % q;
        e >>= 1u;
    }
    return (uint32_t)r;
}

static uint32_t mod_inv(uint32_t a, uint32_t q) {
    return mod_pow(a, q - 2u, q);
}

static uint32_t gcd_u32(uint32_t a, uint32_t b) {
    while (b) {
        uint32_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static int factor_u32(uint32_t x, uint32_t *f, int maxf) {
    int n = 0;
    for (uint32_t p = 2; p * p <= x; ++p) {
        if (x % p == 0) {
            if (n >= maxf) return -1;
            f[n++] = p;
            while (x % p == 0) x /= p;
        }
    }
    if (x > 1) {
        if (n >= maxf) return -1;
        f[n++] = x;
    }
    return n;
}

static int find_primitive_root(uint32_t q, uint32_t *g_out) {
    uint32_t phi = q - 1u;
    uint32_t fac[32];
    int nf = factor_u32(phi, fac, 32);
    if (nf <= 0) return -1;

    for (uint32_t g = 2; g < q; ++g) {
        int ok = 1;
        for (int i = 0; i < nf; ++i) {
            if (mod_pow(g, phi / fac[i], q) == 1u) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            *g_out = g;
            return 0;
        }
    }
    return -1;
}

static void bit_reverse(uint32_t *a, uint32_t n) {
    uint32_t j = 0;
    for (uint32_t i = 1; i < n; ++i) {
        uint32_t bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            uint32_t t = a[i];
            a[i] = a[j];
            a[j] = t;
        }
    }
}

static inline uint32_t mod_u64_barrett(uint64_t x, uint32_t q, uint64_t qinv) {
    __uint128_t z = (__uint128_t)x * qinv;
    uint64_t t = (uint64_t)(z >> 64);
    uint64_t r = x - t * (uint64_t)q;
    if (r >= q) r -= q;
    if (r >= q) r -= q;
    return (uint32_t)r;
}

static inline uint32_t mod_mul_fast(uint32_t a, uint32_t b, uint32_t q, uint64_t qinv) {
    return mod_u64_barrett((uint64_t)a * (uint64_t)b, q, qinv);
}

/* ---------- twiddle tables ---------- */
static int build_twiddle_tables(ntt_ctx_t *ctx) {
    const uint32_t n = ctx->n, q = ctx->q;
    const uint64_t qinv = ctx->qinv;
    ctx->tw_fwd = (uint32_t *)malloc((size_t)(n - 1) * sizeof(uint32_t));
    ctx->tw_inv = (uint32_t *)malloc((size_t)(n - 1) * sizeof(uint32_t));
    if (!ctx->tw_fwd || !ctx->tw_inv) return -1;

    uint32_t off = 0;
    for (uint32_t len = 2; len <= n; len <<= 1u) {
        uint32_t half = len >> 1u;
        uint32_t step = n / len;
        uint32_t wlen_f = mod_pow(ctx->omega, step, q);
        uint32_t wlen_i = mod_pow(ctx->omega_inv, step, q);

        uint32_t wf = 1u, wi = 1u;
        for (uint32_t j = 0; j < half; ++j) {
            ctx->tw_fwd[off + j] = wf;
            ctx->tw_inv[off + j] = wi;
            wf = mod_mul_fast(wf, wlen_f, q, qinv);
            wi = mod_mul_fast(wi, wlen_i, q, qinv);
        }
        off += half;
    }
    return 0;
}

/* ---------- original CT forward/inverse path ---------- */
static void ntt_cyclic(uint32_t *a, const ntt_ctx_t *ctx, int inverse) {
    uint32_t n = ctx->n, q = ctx->q;
    uint64_t qinv = ctx->qinv;
    const uint32_t *tw = inverse ? ctx->tw_inv : ctx->tw_fwd;

    bit_reverse(a, n);

    uint32_t off = 0;
    for (uint32_t len = 2; len <= n; len <<= 1u) {
        uint32_t half = len >> 1u;
        const uint32_t *tw_stage = &tw[off];

        for (uint32_t i = 0; i < n; i += len) {
            for (uint32_t j = 0; j < half; ++j) {
                uint32_t u = a[i + j];
                uint32_t v = mod_mul_fast(a[i + j + half], tw_stage[j], q, qinv);
                a[i + j] = mod_add(u, v, q);
                a[i + j + half] = mod_sub(u, v, q);
            }
        }
        off += half;
    }

    if (inverse) {
        for (uint32_t i = 0; i < n; ++i) {
            a[i] = mod_mul_fast(a[i], ctx->n_inv, q, qinv);
        }
    }
}

/* ---------- new DIF forward: bit-reversed output ---------- */
static void ntt_cyclic_fwd_dif_bitrev(uint32_t *a, const ntt_ctx_t *ctx) {
    const uint32_t n = ctx->n, q = ctx->q;
    const uint64_t qinv = ctx->qinv;

    for (uint32_t len = n; len >= 2; len >>= 1u) {
        uint32_t half = len >> 1u;
        const uint32_t *tw_stage = &ctx->tw_fwd[half - 1u];

        for (uint32_t i = 0; i < n; i += len) {
            for (uint32_t j = 0; j < half; ++j) {
                uint32_t u = a[i + j];
                uint32_t v = a[i + j + half];
                a[i + j] = mod_add(u, v, q);
                uint32_t t = mod_sub(u, v, q);
                a[i + j + half] = mod_mul_fast(t, tw_stage[j], q, qinv);
            }
        }
    }
}

/* ---------- public API ---------- */
int ntt_ctx_init(ntt_ctx_t *ctx, uint32_t n, uint32_t q) {
    if (!ctx || n == 0 || (n & (n - 1u)) != 0) return -1;
    if (((q - 1u) % (2u * n)) != 0) return -1;
    if (gcd_u32(n, q) != 1u) return -1;

    memset(ctx, 0, sizeof(*ctx));
    ctx->n = n;
    ctx->q = q;
    ctx->qinv = (uint64_t)((((__uint128_t)1) << 64) / q);

    if (find_primitive_root(q, &ctx->primitive_root) != 0) return -1;

    uint32_t exp = (q - 1u) / (2u * n);
    ctx->psi = mod_pow(ctx->primitive_root, exp, q);
    ctx->omega = mod_mul_fast(ctx->psi, ctx->psi, q, ctx->qinv);

    if (mod_pow(ctx->psi, n, q) != (q - 1u)) return -1;
    if (mod_pow(ctx->omega, n, q) != 1u) return -1;

    ctx->psi_inv = mod_inv(ctx->psi, q);
    ctx->omega_inv = mod_inv(ctx->omega, q);
    ctx->n_inv = mod_inv(n % q, q);

    ctx->psi_pows = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    ctx->psi_inv_pows = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
    if (!ctx->psi_pows || !ctx->psi_inv_pows) {
        ntt_ctx_free(ctx);
        return -1;
    }

    ctx->psi_pows[0] = 1u;
    ctx->psi_inv_pows[0] = 1u;
    for (uint32_t i = 1; i < n; ++i) {
        ctx->psi_pows[i] = mod_mul_fast(ctx->psi_pows[i - 1], ctx->psi, q, ctx->qinv);
        ctx->psi_inv_pows[i] = mod_mul_fast(ctx->psi_inv_pows[i - 1], ctx->psi_inv, q, ctx->qinv);
    }

    if (build_twiddle_tables(ctx) != 0) {
        ntt_ctx_free(ctx);
        return -1;
    }
    return 0;
}

void ntt_ctx_free(ntt_ctx_t *ctx) {
    if (!ctx) return;
    free(ctx->psi_pows);
    free(ctx->psi_inv_pows);
    free(ctx->tw_fwd);
    free(ctx->tw_inv);
    memset(ctx, 0, sizeof(*ctx));
}

int ntt_forward(uint32_t *a, const ntt_ctx_t *ctx) {
    if (!a || !ctx) return -1;
    for (uint32_t i = 0; i < ctx->n; ++i) {
        a[i] = mod_mul_fast(a[i], ctx->psi_pows[i], ctx->q, ctx->qinv);
    }
    ntt_cyclic(a, ctx, 0);
    return 0;
}

int ntt_forward_bitrev(uint32_t *a, const ntt_ctx_t *ctx) {
    if (!a || !ctx) return -1;
    for (uint32_t i = 0; i < ctx->n; ++i) {
        a[i] = mod_mul_fast(a[i], ctx->psi_pows[i], ctx->q, ctx->qinv);
    }
    ntt_cyclic_fwd_dif_bitrev(a, ctx);
    return 0;
}

int ntt_inverse(uint32_t *a, const ntt_ctx_t *ctx) {
    if (!a || !ctx) return -1;
    ntt_cyclic(a, ctx, 1);
    for (uint32_t i = 0; i < ctx->n; ++i) {
        a[i] = mod_mul_fast(a[i], ctx->psi_inv_pows[i], ctx->q, ctx->qinv);
    }
    return 0;
}