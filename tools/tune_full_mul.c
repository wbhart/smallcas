#define _POSIX_C_SOURCE 200809L
#include "smallcas.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef sc_value *(*mul_fn)(sc_context *, const sc_value *, const sc_value *);

typedef struct {
    size_t loss, win, cut;
    double loss_ratio, win_ratio;
} tune_result;

static double now(void)
{
    struct timespec t;

    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;

    return (x > y) - (x < y);
}

static double median(double *a, size_t n)
{
    qsort(a, n, sizeof(double), cmp_double);
    return a[n / 2];
}

static sc_value *random_poly(sc_context *ctx, sc_parent *r, size_t n, size_t bits,
                             gmp_randstate_t state)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);

    if (f == NULL)
        return NULL;
    for (size_t i = 0; i < n; i++) {
        mpz_urandomb(f->data.zz_poly.coeff[i], state, (mp_bitcnt_t)bits);
        if (bits != 0)
            mpz_setbit(f->data.zz_poly.coeff[i], (mp_bitcnt_t)(bits - 1));
        if (i & 1)
            mpz_neg(f->data.zz_poly.coeff[i], f->data.zz_poly.coeff[i]);
    }
    return f;
}

static size_t log2ceil(size_t n)
{
    size_t k = 0;

    if (n != 0)
        n--;
    while (n != 0)
        k++, n >>= 1;
    return k;
}

static sc_value *mul_ks(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t small = a->data.zz_poly.length < b->data.zz_poly.length ?
                   a->data.zz_poly.length : b->data.zz_poly.length;
    size_t ba = sc_zz_poly_max_abs_bits_raw(a), bb = sc_zz_poly_max_abs_bits_raw(b);
    mp_bitcnt_t bits = (mp_bitcnt_t)(ba + bb + log2ceil(small) + 1);

    return sc_zz_poly_mul_ks(ctx, a, b, bits);
}

static sc_value *mul_toom3(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t n = a->data.zz_poly.length > b->data.zz_poly.length ?
               a->data.zz_poly.length : b->data.zz_poly.length;
    sc_zz_poly_toom3_ws ws;

    if (!sc_zz_poly_toom3_ws_init(ctx, &ws, a, b, (n + 2) / 3))
        return NULL;
    return sc_zz_poly_mul_toom3(ctx, &ws);
}

static double run(sc_context *ctx, mul_fn fn, const sc_value *a, const sc_value *b,
                  size_t reps)
{
    double t = now();

    for (size_t i = 0; i < reps; i++) {
        sc_value *r = fn(ctx, a, b);

        if (r == NULL) {
            fprintf(stderr, "tuning multiplication failed: %s\n", ctx->error);
            exit(1);
        }
        sc_value_free(r);
    }
    return now() - t;
}

static void timed_one(sc_context *ctx, mul_fn fn, const sc_value *a,
                      const sc_value *b, double *sum)
{
    double t = now();
    sc_value *r = fn(ctx, a, b);

    *sum += now() - t;
    if (r == NULL) {
        fprintf(stderr, "tuning multiplication failed: %s\n", ctx->error);
        exit(1);
    }
    sc_value_free(r);
}

static double pair_sample(sc_context *ctx, mul_fn old, mul_fn new,
                          const sc_value *a, const sc_value *b, size_t reps,
                          int reverse)
{
    double ta = 0.0, tb = 0.0;

    for (size_t i = 0; i < reps; i++) {
        if ((i ^ (size_t)reverse) & 1) {
            timed_one(ctx, new, a, b, &tb);
            timed_one(ctx, old, a, b, &ta);
        } else {
            timed_one(ctx, old, a, b, &ta);
            timed_one(ctx, new, a, b, &tb);
        }
    }
    return tb / ta;
}

static double ratio(sc_context *ctx, mul_fn old, mul_fn new, const sc_value *a,
                    const sc_value *b, double *relmad)
{
    double v[15], d[15], ta, tb, med;
    size_t reps = 1, ns = 0;

    while (reps < 4096) {
        ta = run(ctx, old, a, b, reps);
        tb = run(ctx, new, a, b, reps);
        if (ta + tb >= 0.010)
            break;
        reps <<= 1;
    }
    for (ns = 0; ns < 15; ns++) {
        v[ns] = pair_sample(ctx, old, new, a, b, reps, (int)(ns & 1));
        if (ns >= 4) {
            double q[15];

            for (size_t i = 0; i <= ns; i++)
                q[i] = v[i];
            med = median(q, ns + 1);
            for (size_t i = 0; i <= ns; i++)
                d[i] = fabs(v[i] - med);
            *relmad = median(d, ns + 1) / med;
            if (*relmad <= 0.01)
                return med;
        }
    }
    med = median(v, ns);
    for (size_t i = 0; i < ns; i++)
        d[i] = fabs(v[i] - med);
    *relmad = median(d, ns) / med;
    return med;
}

static size_t next_size(size_t n)
{
    size_t d = n / 3;

    return n + (d != 0 ? d : 1);
}

static tune_result tune_pair(sc_context *ctx, sc_parent *r, gmp_randstate_t state,
                             const char *name, mul_fn old, mul_fn new, size_t bits,
                             int bits_equal_n, size_t lo, size_t hi, size_t fallback)
{
    tune_result tr = { 0, 0, fallback, 0.0, 0.0 };

    printf("\n%s\n", name);
    for (size_t n = lo; n <= hi;) {
        size_t b = bits_equal_n ? n : bits;
        sc_value *a = random_poly(ctx, r, n, b, state);
        sc_value *c = random_poly(ctx, r, n, b, state);
        double mad = 0.0, q;

        if (a == NULL || c == NULL) {
            fprintf(stderr, "out of memory while tuning\n");
            exit(1);
        }
        q = ratio(ctx, old, new, a, c, &mad);
        printf("  n=%-6zu bits=%-6zu new/old=%7.4f  MAD=%5.2f%%\n",
               n, b, q, 100.0 * mad);
        sc_value_free_many(2, a, c);
        if (q >= 1.05) {
            tr.loss = n;
            tr.loss_ratio = q;
        } else if (q <= 0.95) {
            tr.win = n;
            tr.win_ratio = q;
            if (tr.loss != 0) {
                tr.cut = tr.loss + (tr.win - tr.loss) / 2;
                break;
            }
            tr.cut = n;
            break;
        }
        if (n == hi)
            break;
        {
            size_t next = next_size(n);

            n = next > hi ? hi : next;
        }
    }
    if (tr.win == 0) {
        if (fallback == (size_t)-1)
            puts("  no 5% win found; keeping disabled");
        else
            printf("  no 5%% win found; keeping %zu\n", fallback);
    }
    else if (tr.loss == 0)
        printf("  already >5%% faster at first point; cutoff <= %zu\n", tr.cut);
    else
        printf("  5%% bracket [%zu, %zu], midpoint %zu\n", tr.loss, tr.win, tr.cut);
    return tr;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "PolynomialRing(ZZ)", SC_PARENT_POLY, &SC_ZZ, "x" };
    gmp_randstate_t state;
    tune_result ks, kar, low, toom, ntt, ssa;
    size_t ntt_cut, toom_fallback;

    sc_context_init(&ctx);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20260923);
    sc_tune_mul_ntt_cutoff = (size_t)-1;
    sc_tune_mul_ssa_cutoff = (size_t)-1;

    puts("Full multiplication tuning: 5% loss/win bracket, midpoint cutoff, 1% MAD.");
    ks = tune_pair(&ctx, &r, state, "classical -> Kronecker (8-bit coefficients)",
                   sc_zz_poly_mul_classical, mul_ks, 8, 0, 8, 256,
                   sc_tune_mul_ks_cutoff);
    sc_tune_mul_ks_cutoff = ks.cut;
    kar = tune_pair(&ctx, &r, state, "classical -> Karatsuba (256-bit coefficients)",
                    sc_zz_poly_mul_classical, sc_zz_poly_mul_karatsuba,
                    256, 0, 4, 256, sc_tune_mul_karatsuba_cutoff);
    sc_tune_mul_karatsuba_cutoff = kar.cut;
    low = tune_pair(&ctx, &r, state, "classical -> Karatsuba (64-bit coefficients)",
                    sc_zz_poly_mul_classical, sc_zz_poly_mul_karatsuba,
                    64, 0, 8, 384, sc_tune_mul_karatsuba_low_bits_cutoff);
    sc_tune_mul_karatsuba_low_bits_cutoff = low.cut;
    toom_fallback = sc_tune_mul_toom3_cutoff;
    sc_tune_mul_toom3_cutoff = (size_t)-1;
    toom = tune_pair(&ctx, &r, state, "lower dispatcher -> Toom-3 (256-bit coefficients)",
                     sc_zz_poly_mul, mul_toom3, 256, 0, 24, 768, toom_fallback);
    sc_tune_mul_toom3_cutoff = toom.cut;
    ntt = tune_pair(&ctx, &r, state, "lower dispatcher -> CRT-NTT (48-bit coefficients)",
                    sc_zz_poly_mul, sc_zz_poly_mul_ntt, 48, 0, 96, 16384, (size_t)-1);
    ntt_cut = ntt.cut;
    if (ntt.win != 0) {
        sc_value *a = random_poly(&ctx, &r, ntt.cut, 48, state);
        sc_value *b = random_poly(&ctx, &r, ntt.cut, 48, state);
        size_t np = sc_zz_poly_ntt_nprimes(a, b);

        if (np != 0)
            ntt_cut = ntt.cut / np;
        sc_value_free_many(2, a, b);
    }
    ssa = tune_pair(&ctx, &r, state, "lower dispatcher -> SSA (bits = length)",
                    sc_zz_poly_mul, sc_zz_poly_mul_ssa, 0, 1, 32, 2048, (size_t)-1);

    puts("\nSuggested full-product entries for include/tuning.h:");
    printf("#define SC_MUL_KS_CUTOFF ((size_t)%zu)\n", ks.cut);
    printf("#define SC_MUL_TOOM3_CUTOFF ((size_t)%zu)\n", toom.cut);
    printf("#define SC_MUL_KARATSUBA_CUTOFF ((size_t)%zu)\n", kar.cut);
    printf("#define SC_MUL_KARATSUBA_LOW_BITS ((size_t)128)\n");
    printf("#define SC_MUL_KARATSUBA_LOW_BITS_CUTOFF ((size_t)%zu)\n", low.cut);
    if (ntt.win == 0)
        puts("#define SC_MUL_NTT_CUTOFF ((size_t)-1)");
    else
        printf("#define SC_MUL_NTT_CUTOFF ((size_t)%zu) /* length / CRT primes */\n",
               ntt_cut);
    if (ssa.win == 0)
        puts("#define SC_MUL_SSA_CUTOFF ((size_t)-1)");
    else
        printf("#define SC_MUL_SSA_CUTOFF ((size_t)%zu)\n", ssa.cut);

    gmp_randclear(state);
    sc_context_clear(&ctx);
    return 0;
}
