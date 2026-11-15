/*
 * Wolff single-cluster sampler for the 2D Ising model at the critical point on an
 * L x L torus (J = 1, T_c = 2 / ln(1 + sqrt 2)).
 *
 * usage: ising_wolff L nmeas seed [start]       start = 0 random (default), 1 ordered
 *
 * Between two measurements a fixed number k of clusters is flipped, k = ceil(n / <C>),
 * with <C> the mean cluster size estimated in a warm-up phase.  The schedule does not
 * depend on the cluster sizes realised after the warm-up, so the measurements are not
 * size-biased.  Each measurement records |M|, M^2, M^4, E, E^2 (per spin).  The output is
 * NB = 100 lines of block means:   |M|  M^2  M^4  E  E^2
 * Run parameters (L, k, <C>) go to stderr.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

static uint64_t s[2];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next(void) {            /* xoroshiro128+ */
    uint64_t s0 = s[0], s1 = s[1], r = s0 + s1;
    s1 ^= s0; s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); s[1] = rotl(s1, 36);
    return r;
}
static inline double U(void) { return (next() >> 11) * (1.0 / 9007199254740992.0); }

static int L, n, *sp, *stack;
static double padd;

static int flip_cluster(void) {
    int seed = (int)(U() * n), sg = sp[seed], top = 0, cs = 0;
    sp[seed] = -sg; stack[top++] = seed;
    while (top) {
        int x = stack[--top]; cs++;
        int xi = x % L, yi = x / L;
        int nb[4] = { yi * L + (xi + 1) % L, yi * L + (xi + L - 1) % L,
                      ((yi + 1) % L) * L + xi, ((yi + L - 1) % L) * L + xi };
        for (int j = 0; j < 4; j++) {
            int y = nb[j];
            if (sp[y] == sg && U() < padd) { sp[y] = -sg; stack[top++] = y; }
        }
    }
    return cs;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s L nmeas seed [start]\n", argv[0]); return 1; }
    L = atoi(argv[1]); long nmeas = atol(argv[2]);
    uint64_t seed = (uint64_t)atol(argv[3]); int start = argc > 4 ? atoi(argv[4]) : 0;
    s[0] = seed * 0x9E3779B97F4A7C15ULL + 1; s[1] = 0xDEADBEEF12345ULL ^ s[0];
    for (int i = 0; i < 20; i++) next();
    n = L * L; double K = 0.5 * log(1 + sqrt(2.0)); padd = 1 - exp(-2 * K);
    sp = malloc(n * sizeof(int)); stack = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) sp[i] = start ? 1 : ((U() < 0.5) ? 1 : -1);

    long tot = 0;
    for (int w = 0; w < 20000; w++) { int cs = flip_cluster(); if (w >= 10000) tot += cs; }
    double meanC = tot / 10000.0; int k = (int)ceil(n / meanC); if (k < 1) k = 1;

    enum { NB = 100 };
    long per = nmeas / NB; double b[NB][5] = {{0}};
    for (long t = 0; t < nmeas; t++) {
        for (int f = 0; f < k; f++) flip_cluster();
        long M = 0, B = 0;                                /* B = sum of aligned bonds */
        for (int i = 0; i < n; i++) {
            int xi = i % L, yi = i / L;
            M += sp[i];
            B += sp[i] * (sp[yi * L + (xi + 1) % L] + sp[((yi + 1) % L) * L + xi]);
        }
        double m = fabs((double)M / n), e = -(double)B / n;
        int q = (int)(t / per); if (q >= NB) q = NB - 1;
        b[q][0] += m; b[q][1] += m * m; b[q][2] += m * m * m * m; b[q][3] += e; b[q][4] += e * e;
    }
    for (int q = 0; q < NB; q++)
        printf("%.15g %.15g %.15g %.15g %.15g\n", b[q][0] / per, b[q][1] / per, b[q][2] / per,
               b[q][3] / per, b[q][4] / per);
    fprintf(stderr, "L=%d k=%d meanC=%.1f nmeas=%ld seed=%lu start=%d\n", L, k, meanC, nmeas,
            (unsigned long)seed, start);
    return 0;
}
