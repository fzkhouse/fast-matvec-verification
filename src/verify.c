#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct { const char *name; uint32_t q, n, k, l; } params_t;
typedef struct { uint32_t q, n, rows, cols; uint32_t *v; } grid_t;

static const params_t PARAMS[] = {
    {"falcon512",12289,512,1,1}, {"falcon1024",12289,1024,2,2},
    {"kyber512",3329,256,2,2}, {"kyber1024",3329,256,4,4},
    {"dilithium2",8380417,256,4,4}, {"dilithium5",8380417,256,8,7}
};

static const params_t *get_params(const char *name) {
    for (size_t i=0; i<sizeof(PARAMS)/sizeof(PARAMS[0]); ++i)
        if (strcmp(PARAMS[i].name, name) == 0) return &PARAMS[i];
    return NULL;
}

static int load_grid(const char *path, grid_t *g) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno)); return 0; }
    if (fscanf(f, "%u%u%u%u", &g->q, &g->n, &g->rows, &g->cols) != 4) { fclose(f); return 0; }
    size_t count = (size_t)g->rows * g->cols * g->n;
    g->v = malloc(count * sizeof(*g->v));
    if (!g->v) { fclose(f); return 0; }
    for (size_t i=0; i<count; ++i) {
        unsigned long long x;
        if (fscanf(f, "%llu", &x) != 1) { free(g->v); g->v=NULL; fclose(f); return 0; }
        g->v[i] = (uint32_t)(x % g->q);
    }
    fclose(f); return 1;
}

static void free_grid(grid_t *g) { free(g->v); memset(g, 0, sizeof(*g)); }

static inline uint64_t tick_ns(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t)t.tv_sec * 1000000000ull + (uint64_t)t.tv_nsec;
}

static inline uint64_t rng_next(uint64_t *s) {
    *s ^= *s << 13; *s ^= *s >> 7; *s ^= *s << 17; return *s;
}

/* Negacyclic multiplication in a portable, deliberately readable form. */
static void poly_mul(const uint32_t *a, const uint32_t *b, uint32_t *out, uint32_t n, uint32_t q) {
    memset(out, 0, (size_t)n * sizeof(*out));
    for (uint32_t i=0; i<n; ++i) {
        uint64_t ai=a[i];
        for (uint32_t j=0; j<n; ++j) {
            uint32_t p=i+j;
            int64_t term=(int64_t)(ai*b[j]) * (p<n ? 1 : -1);
            int64_t x=(int64_t)out[p%n]+term;
            x%=q; if (x<0) x+=q; out[p%n]=(uint32_t)x;
        }
    }
}

static int compatible(const params_t *p, const grid_t *a, const grid_t *z, const grid_t *u) {
    return a->q==p->q && a->n==p->n && a->rows==p->k && a->cols==p->l &&
           z->q==p->q && z->n==p->n && z->cols==p->l &&
           u->q==p->q && u->n==p->n && u->cols==p->k && z->rows==u->rows;
}

static int baseline(const params_t *p, const grid_t *a, const grid_t *z, const grid_t *u) {
    uint32_t *acc=malloc((size_t)p->n*sizeof(*acc));
    if (!acc) return 0;
    for (uint32_t r=0; r<z->rows; ++r) {
        for (uint32_t i=0; i<p->k; ++i) {
            memset(acc,0,(size_t)p->n*sizeof(*acc));
            for (uint32_t j=0; j<p->l; ++j) {
                uint32_t *tmp=malloc((size_t)p->n*sizeof(*tmp));
                if (!tmp) { free(acc); return 0; }
                poly_mul(a->v+((size_t)(i*p->l+j)*p->n), z->v+((size_t)(r*p->l+j)*p->n), tmp,p->n,p->q);
                for (uint32_t t=0;t<p->n;++t) acc[t]=(acc[t]+tmp[t])%p->q;
                free(tmp);
            }
            if (memcmp(acc,u->v+((size_t)(r*p->k+i)*p->n),(size_t)p->n*sizeof(*acc))!=0) { free(acc); return 0; }
        }
    }
    free(acc); return 1;
}

/* Batch random linearization: aggregate all rows first, then multiply A once. */
static int reference(const params_t *p, const grid_t *a, const grid_t *z, const grid_t *u) {
    size_t polybytes=(size_t)p->n*sizeof(uint32_t);
    uint32_t *zz=calloc((size_t)p->l*p->n,sizeof(*zz));
    uint32_t *uu=calloc((size_t)p->k*p->n,sizeof(*uu));
    uint32_t *tmp=malloc(polybytes);
    if (!zz || !uu || !tmp) { free(zz);free(uu);free(tmp);return 0; }
    uint64_t seed=0x9e3779b97f4a7c15ull;
    for (uint32_t r=0;r<z->rows;++r) {
        uint32_t c=(uint32_t)(rng_next(&seed)%p->q);
        for (uint32_t j=0;j<p->l;++j) for (uint32_t t=0;t<p->n;++t)
            zz[j*p->n+t]=(uint32_t)((zz[j*p->n+t]+(uint64_t)c*z->v[((size_t)(r*p->l+j)*p->n)+t])%p->q);
        for (uint32_t i=0;i<p->k;++i) for (uint32_t t=0;t<p->n;++t)
            uu[i*p->n+t]=(uint32_t)((uu[i*p->n+t]+(uint64_t)c*u->v[((size_t)(r*p->k+i)*p->n)+t])%p->q);
    }
    for (uint32_t i=0;i<p->k;++i) {
        uint32_t *acc=calloc(p->n,sizeof(*acc)); if(!acc){free(zz);free(uu);free(tmp);return 0;}
        for (uint32_t j=0;j<p->l;++j) {
            poly_mul(a->v+((size_t)(i*p->l+j)*p->n),zz+j*p->n,tmp,p->n,p->q);
            for(uint32_t t=0;t<p->n;++t) acc[t]=(acc[t]+tmp[t])%p->q;
        }
        if(memcmp(acc,uu+i*p->n,polybytes)!=0){free(acc);free(zz);free(uu);free(tmp);return 0;}
        free(acc);
    }
    free(zz);free(uu);free(tmp);return 1;
}

int main(int argc, char **argv) {
    if (argc!=5) { fprintf(stderr,"usage: %s <scheme> <reference|baseline> <data-dir> <seed>\n",argv[0]); return 2; }
    const params_t *p=get_params(argv[1]); if(!p || (strcmp(argv[2],"reference")&&strcmp(argv[2],"baseline"))){fprintf(stderr,"invalid scheme or implementation\n");return 2;}
    char path[4096]; grid_t a={0},z={0},u={0};
    snprintf(path,sizeof(path),"%s/A.txt",argv[3]); if(!load_grid(path,&a))return 2;
    snprintf(path,sizeof(path),"%s/z.txt",argv[3]); if(!load_grid(path,&z)){free_grid(&a);return 2;}
    snprintf(path,sizeof(path),"%s/u.txt",argv[3]); if(!load_grid(path,&u)){free_grid(&a);free_grid(&z);return 2;}
    if(!compatible(p,&a,&z,&u)){fprintf(stderr,"data header does not match %s\n",p->name);free_grid(&a);free_grid(&z);free_grid(&u);return 2;}
    uint64_t begin=tick_ns();
    int ok=!strcmp(argv[2],"reference") ? reference(p,&a,&z,&u) : baseline(p,&a,&z,&u);
    double sec=(double)(tick_ns()-begin)/1e9;
    printf("scheme=%s implementation=%s rows=%u K=%u L=%u q=%u n=%u\n",p->name,argv[2],z.rows,p->k,p->l,p->q,p->n);
    printf("result=%s compute_time=%.6f s compute_cycles=0\n",ok?"PASS":"FAIL",sec);
    free_grid(&a);free_grid(&z);free_grid(&u); return ok?0:1;
}
