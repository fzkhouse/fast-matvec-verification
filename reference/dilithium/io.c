#define _POSIX_C_SOURCE 200112L
#include "io.h"
#include <stdio.h>
#include <stdlib.h>

static void* aligned_malloc(size_t size) {
    void *p = NULL;
    if (posix_memalign(&p, 32, size) != 0) return NULL;
    return p;
}

int grid_save(const char *path, const poly_grid_t *g) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "%u %u %u %u\n", DILITHIUM_Q, DILITHIUM_N, g->rows, g->cols);
    for (uint32_t i = 0; i < g->rows * g->cols; i++) {
        for (int j = 0; j < DILITHIUM_N; j++) {
            fprintf(f, "%d ", g->data[i].coeffs[j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}

int grid_load(const char *path, poly_grid_t *g, uint32_t expect_q, uint32_t expect_n) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    uint32_t q,n,r,c;
    if (fscanf(f, "%u %u %u %u", &q,&n,&r,&c) != 4) { fclose(f); return -1; }
    if (q != expect_q || n != expect_n) { fclose(f); return -1; }
    g->rows = r; g->cols = c;
    g->data = (poly_t*)aligned_malloc((size_t)r * c * sizeof(poly_t));
    if (!g->data) { fclose(f); return -1; }
    for (uint32_t i = 0; i < r * c; i++) {
        for (int j = 0; j < DILITHIUM_N; j++) {
            if (fscanf(f, "%d", &g->data[i].coeffs[j]) != 1) { fclose(f); return -1; }
        }
    }
    fclose(f);
    return 0;
}

void grid_free(poly_grid_t *g) {
    if (g && g->data) free(g->data);
    g->data = NULL;
}

int table_save(const char *path, const table_t *tb) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "%u %u %u %u %u %u\n", tb->q, tb->n, tb->R, tb->T, tb->K_, tb->L_);
    for (uint32_t i = 0; i < tb->R * tb->T; i++) fprintf(f, "%u ", tb->randoms[i]);
    fprintf(f, "\n");
    for (uint32_t i = 0; i < tb->R * tb->K_ * tb->L_; i++) {
        for (int j = 0; j < DILITHIUM_N; j++) fprintf(f, "%d ", tb->A[i].coeffs[j]);
        fprintf(f, "\n");
    }
    for (uint32_t i = 0; i < tb->R * tb->K_; i++) {
        for (int j = 0; j < DILITHIUM_N; j++) fprintf(f, "%d ", tb->gamma[i].coeffs[j]);
        fprintf(f, "\n");
    }
    for (uint32_t i = 0; i < tb->R * tb->L_; i++) {
        for (int j = 0; j < DILITHIUM_N; j++) fprintf(f, "%d ", tb->ATgamma[i].coeffs[j]);
        fprintf(f, "\n");
    }
    fclose(f);
    return 0;
}

int table_load(const char *path, table_t *tb, uint32_t expect_q, uint32_t expect_n) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    if (fscanf(f, "%u %u %u %u %u %u", &tb->q, &tb->n, &tb->R, &tb->T, &tb->K_, &tb->L_) != 6) { fclose(f); return -1; }
    if (tb->q != expect_q || tb->n != expect_n) { fclose(f); return -1; }
    tb->randoms = (uint32_t*)malloc((size_t)tb->R * tb->T * sizeof(uint32_t));
    tb->A = (poly_t*)aligned_malloc((size_t)tb->R * tb->K_ * tb->L_ * sizeof(poly_t));
    tb->gamma = (poly_t*)aligned_malloc((size_t)tb->R * tb->K_ * sizeof(poly_t));
    tb->ATgamma = (poly_t*)aligned_malloc((size_t)tb->R * tb->L_ * sizeof(poly_t));
    if (!tb->randoms || !tb->A || !tb->gamma || !tb->ATgamma) { fclose(f); return -1; }
    for (uint32_t i = 0; i < tb->R * tb->T; i++) if (fscanf(f, "%u", &tb->randoms[i]) != 1) { fclose(f); return -1; }
    for (uint32_t i = 0; i < tb->R * tb->K_ * tb->L_; i++)
        for (int j = 0; j < DILITHIUM_N; j++) if (fscanf(f, "%d", &tb->A[i].coeffs[j]) != 1) { fclose(f); return -1; }
    for (uint32_t i = 0; i < tb->R * tb->K_; i++)
        for (int j = 0; j < DILITHIUM_N; j++) if (fscanf(f, "%d", &tb->gamma[i].coeffs[j]) != 1) { fclose(f); return -1; }
    for (uint32_t i = 0; i < tb->R * tb->L_; i++)
        for (int j = 0; j < DILITHIUM_N; j++) if (fscanf(f, "%d", &tb->ATgamma[i].coeffs[j]) != 1) { fclose(f); return -1; }
    fclose(f);
    return 0;
}

void table_free(table_t *tb) {
    if (!tb) return;
    if (tb->A) free(tb->A);
    if (tb->gamma) free(tb->gamma);
    if (tb->ATgamma) free(tb->ATgamma);
    if (tb->randoms) free(tb->randoms);
    tb->A = tb->gamma = tb->ATgamma = NULL;
    tb->randoms = NULL;
}
