#define _POSIX_C_SOURCE 200809L
#include "smallcas.h"
#include "smallcas_fft.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
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

typedef struct {
    sc_fft_mod m;
    sc_fft_plan p;
    mp_ptr a, scratch;
} fft_case;

static int fft_case_init(fft_case *c, unsigned logn)
{
    mp_bitcnt_t bits = (mp_bitcnt_t)1 << (logn - 1);
    size_t words;

    memset(c, 0, sizeof(*c));
    if (!sc_fft_mod_init(&c->m, 1, bits, 2, logn) ||
        !sc_fft_plan_init(&c->p, logn, &c->m))
        goto fail;
    if (c->p.len > SIZE_MAX / (size_t)c->m.n)
        goto fail;
    words = c->p.len * (size_t)c->m.n;
    c->a = calloc(words, sizeof(mp_limb_t));
    c->scratch = calloc(4 * (size_t)c->m.n + 1, sizeof(mp_limb_t));
    if (c->a == NULL || c->scratch == NULL)
        goto fail;
    for (size_t i = 0; i < c->p.len; i++)
        sc_fft_set_ui(sc_fft_entry(c->a, i, &c->m),
                      (mp_limb_t)(1103515245u * (unsigned)i + 12345u), &c->m);
    return 1;
fail:
    free(c->scratch);
    free(c->a);
    sc_fft_plan_clear(&c->p);
    sc_fft_mod_clear(&c->m);
    memset(c, 0, sizeof(*c));
    return 0;
}

static void fft_case_clear(fft_case *c)
{
    free(c->scratch);
    free(c->a);
    sc_fft_plan_clear(&c->p);
    sc_fft_mod_clear(&c->m);
}

static void fft_roundtrip(fft_case *c, int mfa)
{
    if (mfa) {
        sc_fft_forward_mfa(c->a, &c->p, &c->m, c->scratch);
        sc_fft_inverse_mfa(c->a, &c->p, &c->m, c->scratch);
    } else {
        sc_fft_forward(c->a, &c->p, &c->m, c->scratch);
        sc_fft_inverse(c->a, &c->p, &c->m, c->scratch);
    }
}

static double fft_run(fft_case *c, int mfa, size_t reps)
{
    double t = now();

    for (size_t i = 0; i < reps; i++)
        fft_roundtrip(c, mfa);
    return now() - t;
}

static double fft_pair_sample(fft_case *c, size_t reps, int reverse)
{
    double ta = 0.0, tb = 0.0;

    for (size_t i = 0; i < reps; i++) {
        double t;

        if ((i ^ (size_t)reverse) & 1) {
            t = now(), fft_roundtrip(c, 1), tb += now() - t;
            t = now(), fft_roundtrip(c, 0), ta += now() - t;
        } else {
            t = now(), fft_roundtrip(c, 0), ta += now() - t;
            t = now(), fft_roundtrip(c, 1), tb += now() - t;
        }
    }
    return tb / ta;
}

static double fft_ratio(fft_case *c, double *relmad)
{
    double v[15], d[15], ta, tb, med;
    size_t reps = 1, ns = 0;

    while (reps < 4096) {
        ta = fft_run(c, 0, reps);
        tb = fft_run(c, 1, reps);
        if (ta + tb >= 0.010)
            break;
        reps <<= 1;
    }
    for (ns = 0; ns < 15; ns++) {
        v[ns] = fft_pair_sample(c, reps, (int)(ns & 1));
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

static unsigned tune_mfa(void)
{
    unsigned loss = 0, win = 0, cut = UINT_MAX;

    puts("\nSSA radix-2 -> MFA: minimum SSA Fermat ring, forward+inverse.");
    for (unsigned logn = SC_FFT_MFA_BASE_LOG + 1; logn <= 15; logn++) {
        fft_case c;
        double mad = 0.0, q, mib;

        if (!fft_case_init(&c, logn)) {
            fprintf(stderr, "out of memory while tuning MFA at logN=%u\n", logn);
            exit(1);
        }
        mib = (double)c.p.len * (double)c.m.n * sizeof(mp_limb_t) / 1048576.0;
        q = fft_ratio(&c, &mad);
        printf("  logN=%-2u N=%-6zu vector=%7.1f MiB MFA/radix=%7.4f MAD=%5.2f%%\n",
               logn, c.p.len, mib, q, 100.0 * mad);
        fft_case_clear(&c);
        if (q >= 1.05)
            loss = logn;
        else if (q <= 0.95) {
            win = logn;
            cut = loss != 0 ? (loss + win + 1) / 2 : win;
            break;
        }
    }
    if (win == 0)
        puts("  no 5% MFA win through logN=15; disabling MFA");
    else if (loss == 0)
        printf("  already >5%% faster at first winning depth; cutoff <= %u\n", cut);
    else
        printf("  5%% bracket [%u, %u], cutoff depth %u\n", loss, win, cut);
    return cut;
}

static int write_tuning(const tune_result *ks, const tune_result *toom,
                        const tune_result *kar, const tune_result *low,
                        const tune_result *ntt, size_t ntt_cut,
                        const tune_result *ssa, unsigned mfa_cut)
{
    const char *tmp = "include/tuning.h.tmp";
    const char *dst = "include/tuning.h";
    FILE *f = fopen(tmp, "w");
    int bad;

    if (f == NULL) {
        perror(tmp);
        return 0;
    }
    fputs("/* Machine-local; generated by make tune. */\n", f);
    fputs("#ifndef SMALLCAS_TUNING_LOCAL_H\n#define SMALLCAS_TUNING_LOCAL_H\n\n", f);
    fputs("#ifndef SC_TUNE\n", f);
    fprintf(f, "#define SC_MUL_KS_CUTOFF ((size_t)%zu)\n", ks->cut);
    fprintf(f, "#define SC_MUL_TOOM3_CUTOFF ((size_t)%zu)\n", toom->cut);
    fprintf(f, "#define SC_MUL_KARATSUBA_CUTOFF ((size_t)%zu)\n", kar->cut);
    fprintf(f, "#define SC_MUL_KARATSUBA_LOW_BITS ((size_t)128)\n");
    fprintf(f, "#define SC_MUL_KARATSUBA_LOW_BITS_CUTOFF ((size_t)%zu)\n", low->cut);
    if (ntt->win == 0)
        fputs("#define SC_MUL_NTT_CUTOFF ((size_t)-1)\n", f);
    else
        fprintf(f, "#define SC_MUL_NTT_CUTOFF ((size_t)%zu)\n", ntt_cut);
    if (ssa->win == 0)
        fputs("#define SC_MUL_SSA_CUTOFF ((size_t)-1)\n", f);
    else
        fprintf(f, "#define SC_MUL_SSA_CUTOFF ((size_t)%zu)\n", ssa->cut);
    if (mfa_cut == UINT_MAX)
        fputs("#define SC_SSA_MFA_CUTOFF_LOG ((unsigned)-1)\n", f);
    else
        fprintf(f, "#define SC_SSA_MFA_CUTOFF_LOG ((unsigned)%u)\n", mfa_cut);
    fputs("#endif\n\n#include \"tuning_defaults.h\"\n\n#endif\n", f);
    bad = ferror(f);
    if (fclose(f) != 0)
        bad = 1;
    if (bad || rename(tmp, dst) != 0) {
        if (!bad)
            perror(dst);
        remove(tmp);
        return 0;
    }
    printf("\nWrote %s.  It is machine-local and ignored by Git.\n", dst);
    return 1;
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
    unsigned mfa_cut;

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
    mfa_cut = tune_mfa();
    sc_tune_ssa_mfa_cutoff_log = mfa_cut;
    ssa = tune_pair(&ctx, &r, state, "lower dispatcher -> SSA (bits = length)",
                    sc_zz_poly_mul, sc_zz_poly_mul_ssa, 0, 1, 32, 2048, (size_t)-1);
    if (!write_tuning(&ks, &toom, &kar, &low, &ntt, ntt_cut, &ssa, mfa_cut)) {
        gmp_randclear(state);
        sc_context_clear(&ctx);
        return 1;
    }

    gmp_randclear(state);
    sc_context_clear(&ctx);
    return 0;
}
