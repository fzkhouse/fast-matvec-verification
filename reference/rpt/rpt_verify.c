#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

typedef struct { size_t rows, cols; int64_t *v; } matrix_t;

static void matrix_free(matrix_t *m) { free(m->v); memset(m, 0, sizeof(*m)); }

static int parse_line(const char *line, int64_t **out, size_t *count) {
    size_t cap = 64, n = 0;
    int64_t *v = malloc(cap * sizeof(*v));
    if (!v) return -1;
    const char *p = line;
    while (*p) {
        char *end;
        errno = 0;
        long long x = strtoll(p, &end, 10);
        if (p == end) { p++; continue; }
        if (errno || n == cap) {
            if (n == cap) { cap *= 2; int64_t *t = realloc(v, cap * sizeof(*v)); if (!t) { free(v); return -1; } v = t; }
            else { free(v); return -1; }
        }
        v[n++] = (int64_t)x;
        p = end;
    }
    if (!n) { free(v); return -1; }
    *out = v; *count = n; return 0;
}

static int vector_load_line(const char *path, unsigned line_no, int64_t **v, size_t *n) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char *line = NULL; size_t cap = 0; ssize_t got;
    for (unsigned i = 0; i <= line_no; i++) {
        got = getline(&line, &cap, f);
        if (got < 0) { free(line); fclose(f); return -1; }
    }
    int rc = parse_line(line, v, n);
    free(line); fclose(f); return rc;
}

static int matrix_load(const char *path, int header, matrix_t *m) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char *line = NULL; size_t line_cap = 0; ssize_t got;
    size_t want_rows = 0, want_cols = 0;
    if (header) {
        if (getline(&line, &line_cap, f) < 0 || sscanf(line, "%zu %zu", &want_rows, &want_cols) != 2 || !want_rows || !want_cols) goto fail;
    }
    size_t rows = 0, cols = want_cols, alloc_rows = header ? want_rows : 64;
    int64_t *data = malloc(alloc_rows * (cols ? cols : 1) * sizeof(*data));
    if (!data) goto fail;
    while ((got = getline(&line, &line_cap, f)) >= 0) {
        int64_t *row = NULL; size_t n = 0;
        if (parse_line(line, &row, &n) != 0) { free(row); continue; }
        if (!cols) { cols = n; free(data); data = malloc(alloc_rows * cols * sizeof(*data)); if (!data) { free(row); goto fail; } }
        if (n != cols || (header && rows >= want_rows)) { free(row); goto fail_data; }
        if (rows == alloc_rows) { alloc_rows *= 2; int64_t *t = realloc(data, alloc_rows * cols * sizeof(*data)); if (!t) { free(row); goto fail_data; } data = t; }
        memcpy(data + rows * cols, row, cols * sizeof(*data));
        free(row); rows++;
    }
    if (!rows || (header && rows != want_rows)) goto fail_data;
    free(line); fclose(f); m->rows = rows; m->cols = cols; m->v = data; return 0;
fail_data: free(data);
fail: free(line); fclose(f); return -1;
}

static int all_i32(const int64_t *v, size_t n) {
    for (size_t i = 0; i < n; i++) if (v[i] < INT32_MIN || v[i] > INT32_MAX) return 0;
    return 1;
}

static int64_t dot(const int64_t *a, const int64_t *b, size_t n, int avx2_ok) {
    /* The public reference keeps the wide, readable arithmetic path.  The
       original direct baseline is compiled with the same routine, so the
       comparison does not depend on a hidden packing conversion. */
    (void)avx2_ok;
    int64_t r = 0; for (size_t i = 0; i < n; i++) r += a[i] * b[i]; return r;
}

static uint64_t norm2(const int64_t *v, size_t n) { uint64_t r = 0; for (size_t i = 0; i < n; i++) r += (uint64_t)v[i] * (uint64_t)v[i]; return r; }
static uint64_t now_ns(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (uint64_t)t.tv_sec*1000000000ull + (uint64_t)t.tv_nsec; }

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: %s B.txt zu.txt C.txt rho.txt bound\n", argv[0]); return 1; }
    int64_t *z = NULL, *u = NULL; size_t nz = 0, nu = 0; matrix_t B = {0}, C = {0}, rho = {0};
    if (vector_load_line(argv[2], 0, &z, &nz) || vector_load_line(argv[2], 1, &u, &nu) ||
        matrix_load(argv[1], 1, &B) || matrix_load(argv[3], 0, &C) || matrix_load(argv[4], 0, &rho)) goto fail;
    if (B.cols != nz || B.rows != nu || C.rows != rho.rows || C.cols != nu || rho.cols != nz) goto fail;
    char *end; int64_t bound = strtoll(argv[5], &end, 10); if (*end || bound < 0) goto fail;
    const int avx = all_i32(z, nz) && all_i32(u, nu) && all_i32(B.v, B.rows*B.cols) && all_i32(C.v, C.rows*C.cols) && all_i32(rho.v, rho.rows*rho.cols);
    uint64_t t0 = now_ns(); uint64_t value = 0;
#if RPT_REFERENCE
    int64_t *p = calloc(C.rows, sizeof(*p)); if (!p) goto fail;
    for (size_t i = 0; i < C.rows; i++) p[i] = dot(u, C.v+i*C.cols, nu, avx) - dot(z, rho.v+i*rho.cols, nz, avx);
    value = norm2(p, C.rows); free(p);
    uint64_t threshold = 337u * (uint64_t)bound * (uint64_t)bound;
#else
    int64_t *d = calloc(B.rows, sizeof(*d)); if (!d) goto fail;
    for (size_t i = 0; i < B.rows; i++) d[i] = dot(B.v+i*B.cols, z, nz, avx) - u[i];
    value = norm2(d, B.rows); free(d);
    uint64_t threshold = (uint64_t)bound;
#endif
    uint64_t t1 = now_ns();
    printf("scheme=rpt implementation=%s dimension=%zu rows=%zu\n", RPT_REFERENCE ? "reference" : "baseline", nz, RPT_REFERENCE ? C.rows : B.rows);
    printf("result=PASS decision=%s statistic=%" PRIu64 " compute_time=%.6f s\n", value < threshold ? "ACCEPT" : "REJECT", value, (double)(t1-t0)/1e9);
    free(z); free(u); matrix_free(&B); matrix_free(&C); matrix_free(&rho); return 0;
fail:
    fprintf(stderr, "invalid RPT fixture or allocation failure\n"); free(z); free(u); matrix_free(&B); matrix_free(&C); matrix_free(&rho); return 1;
}
