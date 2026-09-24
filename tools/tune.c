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
typedef sc_value *(*short_fn)(sc_context *, const sc_value *, const sc_value *, size_t);

typedef struct {
    size_t loss, win, cut;
    double loss_ratio, win_ratio;
} tune_result;

typedef struct {
    tune_result series_dc, inverse_newton, series_newton;
    tune_result bidir_base, mulders_base;
    tune_result quo_dc, quo_mulders, quo_newton;
    tune_result divrem_dc, divrem_mulders, divrem_newton;
    tune_result divexact_bidir, pseudodiv_fast, pseudorem_fast;
} division_tuning;

typedef struct {
    tune_result subresultant;
} gcd_tuning;

/* Division crossovers are deliberately bounded: some may never occur in a
   useful range for ZZ[x], and tuning must terminate even in that case. */
#define DIV_TUNE_MAX_POINTS 20u
#define DIV_TUNE_POINT_BUDGET 1.5

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
            fprintf(stderr, "tuned operation failed: %s\n", ctx->error);
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
        fprintf(stderr, "tuned operation failed: %s\n", ctx->error);
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

static double ratio_bounded(sc_context *ctx, mul_fn old, mul_fn new,
                            const sc_value *a, const sc_value *b,
                            double *relmad)
{
    double v[15], d[15], ta, tb, med, start;
    size_t reps = 1, ns = 0;

    while (reps < 4096) {
        ta = run(ctx, old, a, b, reps);
        tb = run(ctx, new, a, b, reps);
        if (ta + tb >= 0.010 || ta + tb >= DIV_TUNE_POINT_BUDGET / 4.0)
            break;
        reps <<= 1;
    }
    start = now();
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
            if (*relmad <= 0.01 || now() - start >= DIV_TUNE_POINT_BUDGET)
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
                        const tune_result *ssa, unsigned mfa_cut,
                        const tune_result *lowdc, size_t low_ntt, size_t low_ssa,
                        size_t high_ntt, size_t high_ssa,
                        const tune_result *midclass, const tune_result *mid63,
                        size_t mid_ntt, size_t mid_ssa,
                        const division_tuning *div, const gcd_tuning *gcd)
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
    fprintf(f, "#define SC_MULLOW_DC_CUTOFF ((size_t)%zu)\n", lowdc->cut);
#define WRITE_CUTOFF(name, value) do { \
        if ((value) == (size_t)-1) \
            fprintf(f, "#define %s ((size_t)-1)\n", (name)); \
        else \
            fprintf(f, "#define %s ((size_t)%zu)\n", (name), (value)); \
    } while (0)
    WRITE_CUTOFF("SC_MULLOW_NTT_CUTOFF", low_ntt);
    WRITE_CUTOFF("SC_MULLOW_SSA_CUTOFF", low_ssa);
    WRITE_CUTOFF("SC_MULHIGH_NTT_CUTOFF", high_ntt);
    WRITE_CUTOFF("SC_MULHIGH_SSA_CUTOFF", high_ssa);
    fprintf(f, "#define SC_MULMID_CLASSICAL_CUTOFF ((size_t)%zu)\n", midclass->cut);
    fprintf(f, "#define SC_MULMID_TOOM63_CUTOFF ((size_t)%zu)\n", mid63->cut);
    WRITE_CUTOFF("SC_MULMID_NTT_CUTOFF", mid_ntt);
    WRITE_CUTOFF("SC_MULMID_SSA_CUTOFF", mid_ssa);
    WRITE_CUTOFF("SC_SERIES_QUO_DC_CUTOFF", div->series_dc.cut);
    WRITE_CUTOFF("SC_INV_SERIES_NEWTON_CUTOFF", div->inverse_newton.cut);
    WRITE_CUTOFF("SC_SERIES_QUO_NEWTON_CUTOFF", div->series_newton.cut);
    WRITE_CUTOFF("SC_BIDIR_QUO_CUTOFF", div->bidir_base.cut);
    WRITE_CUTOFF("SC_MULDERS_QUO_CUTOFF", div->mulders_base.cut);
    WRITE_CUTOFF("SC_QUO_MULDERS_CUTOFF", div->quo_mulders.cut);
    WRITE_CUTOFF("SC_DIVREM_MULDERS_CUTOFF", div->divrem_mulders.cut);
    WRITE_CUTOFF("SC_DIVREM_DC_CUTOFF", div->divrem_dc.cut);
    WRITE_CUTOFF("SC_QUO_NEWTON_CUTOFF", div->quo_newton.cut);
    WRITE_CUTOFF("SC_DIVREM_NEWTON_CUTOFF", div->divrem_newton.cut);
    WRITE_CUTOFF("SC_QUO_DC_CUTOFF", div->quo_dc.cut);
    WRITE_CUTOFF("SC_DIVEXACT_BIDIR_CUTOFF", div->divexact_bidir.cut);
    WRITE_CUTOFF("SC_PSEUDODIV_FAST_CUTOFF", div->pseudodiv_fast.cut);
    WRITE_CUTOFF("SC_PSEUDOREM_FAST_CUTOFF", div->pseudorem_fast.cut);
    WRITE_CUTOFF("SC_GCD_SUBRESULTANT_CUTOFF", gcd->subresultant.cut);
#undef WRITE_CUTOFF
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


static sc_value *mullow_dispatch(sc_context *ctx, const sc_value *a,
                                 const sc_value *b, size_t n)
{
    return sc_zz_poly_mullow(ctx, a, b, n);
}

static sc_value *mulhigh_dispatch(sc_context *ctx, const sc_value *a,
                                  const sc_value *b, size_t n)
{
    return sc_zz_poly_mulhigh(ctx, a, b, n);
}

static sc_value *fft_window(sc_context *ctx, const sc_value *a, const sc_value *b,
                            size_t n, int high, int ntt)
{
    sc_value *p = ntt ? sc_zz_poly_mul_ntt(ctx, a, b) : sc_zz_poly_mul_ssa(ctx, a, b);
    size_t total = a->data.zz_poly.length + b->data.zz_poly.length - 1;
    sc_value view, *r;

    if (p == NULL)
        return NULL;
    view = sc_zz_poly_view(p, high ? total - n : 0, n);
    r = sc_value_copy_checked(ctx, &view);
    sc_value_free(p);
    return r;
}

static sc_value *mullow_ntt(sc_context *ctx, const sc_value *a,
                            const sc_value *b, size_t n)
{
    return fft_window(ctx, a, b, n, 0, 1);
}

static sc_value *mullow_ssa(sc_context *ctx, const sc_value *a,
                            const sc_value *b, size_t n)
{
    return fft_window(ctx, a, b, n, 0, 0);
}

static sc_value *mulhigh_ntt(sc_context *ctx, const sc_value *a,
                             const sc_value *b, size_t n)
{
    return fft_window(ctx, a, b, n, 1, 1);
}

static sc_value *mulhigh_ssa(sc_context *ctx, const sc_value *a,
                             const sc_value *b, size_t n)
{
    return fft_window(ctx, a, b, n, 1, 0);
}

static sc_value *mulmid_classical(sc_context *ctx, const sc_value *a,
                                  const sc_value *b, size_t n)
{
    return sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
}

static sc_value *mulmid_toom42(sc_context *ctx, const sc_value *a,
                               const sc_value *b, size_t n)
{
    return n & 1 ? sc_zz_poly_mulmid_toom42_odd(ctx, a, b, n) :
                   sc_zz_poly_mulmid_toom42(ctx, a, b, n);
}

static sc_value *mulmid_toom63(sc_context *ctx, const sc_value *a,
                               const sc_value *b, size_t n)
{
    sc_zz_poly_toom63_ws ws;

    if (n % 3 != 0)
        return sc_zz_poly_mulmid_toom63_tail(ctx, a, b, n);
    if (!sc_zz_poly_toom63_ws_init(ctx, &ws, a, b, n / 3))
        return NULL;
    return sc_zz_poly_mulmid_toom63(ctx, &ws);
}

static sc_value *mulmid_dispatch(sc_context *ctx, const sc_value *a,
                                 const sc_value *b, size_t n)
{
    return sc_zz_poly_mulmid_balanced(ctx, a, b, n);
}

static void timed_short(sc_context *ctx, short_fn fn, const sc_value *a,
                        const sc_value *b, size_t n, double *sum)
{
    double t = now();
    sc_value *r = fn(ctx, a, b, n);

    *sum += now() - t;
    if (r == NULL) {
        fprintf(stderr, "tuning short product failed: %s\n", ctx->error);
        exit(1);
    }
    sc_value_free(r);
}

static double short_run(sc_context *ctx, short_fn fn, const sc_value *a,
                        const sc_value *b, size_t n, size_t reps)
{
    double t = now();

    for (size_t i = 0; i < reps; i++) {
        sc_value *r = fn(ctx, a, b, n);

        if (r == NULL) {
            fprintf(stderr, "tuning short product failed: %s\n", ctx->error);
            exit(1);
        }
        sc_value_free(r);
    }
    return now() - t;
}

static double short_pair_sample(sc_context *ctx, short_fn old, short_fn new,
                                const sc_value *a, const sc_value *b, size_t n,
                                size_t reps, int reverse)
{
    double ta = 0.0, tb = 0.0;

    for (size_t i = 0; i < reps; i++) {
        if ((i ^ (size_t)reverse) & 1) {
            timed_short(ctx, new, a, b, n, &tb);
            timed_short(ctx, old, a, b, n, &ta);
        } else {
            timed_short(ctx, old, a, b, n, &ta);
            timed_short(ctx, new, a, b, n, &tb);
        }
    }
    return tb / ta;
}

static double short_ratio(sc_context *ctx, short_fn old, short_fn new,
                          const sc_value *a, const sc_value *b, size_t n,
                          double *relmad)
{
    double v[15], d[15], ta, tb, med;
    size_t reps = 1, ns = 0;

    while (reps < 4096) {
        ta = short_run(ctx, old, a, b, n, reps);
        tb = short_run(ctx, new, a, b, n, reps);
        if (ta + tb >= 0.010)
            break;
        reps <<= 1;
    }
    for (ns = 0; ns < 15; ns++) {
        v[ns] = short_pair_sample(ctx, old, new, a, b, n, reps, (int)(ns & 1));
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

static double short_ratio_bounded(sc_context *ctx, short_fn old, short_fn new,
                                  const sc_value *a, const sc_value *b, size_t n,
                                  double *relmad)
{
    double v[15], d[15], ta, tb, med, start;
    size_t reps = 1, ns = 0;

    while (reps < 4096) {
        ta = short_run(ctx, old, a, b, n, reps);
        tb = short_run(ctx, new, a, b, n, reps);
        if (ta + tb >= 0.010 || ta + tb >= DIV_TUNE_POINT_BUDGET / 4.0)
            break;
        reps <<= 1;
    }
    start = now();
    for (ns = 0; ns < 15; ns++) {
        v[ns] = short_pair_sample(ctx, old, new, a, b, n, reps, (int)(ns & 1));
        if (ns >= 4) {
            double q[15];

            for (size_t i = 0; i <= ns; i++)
                q[i] = v[i];
            med = median(q, ns + 1);
            for (size_t i = 0; i <= ns; i++)
                d[i] = fabs(v[i] - med);
            *relmad = median(d, ns + 1) / med;
            if (*relmad <= 0.01 || now() - start >= DIV_TUNE_POINT_BUDGET)
                return med;
        }
    }
    med = median(v, ns);
    for (size_t i = 0; i < ns; i++)
        d[i] = fabs(v[i] - med);
    *relmad = median(d, ns) / med;
    return med;
}

static tune_result tune_short_pair(sc_context *ctx, sc_parent *r,
                                   gmp_randstate_t state, const char *name,
                                   short_fn old, short_fn new, size_t bits,
                                   int bits_equal_n, int middle,
                                   size_t lo, size_t hi, size_t fallback)
{
    tune_result tr = { 0, 0, fallback, 0.0, 0.0 };

    printf("\n%s\n", name);
    for (size_t n = lo; n <= hi;) {
        size_t b = bits_equal_n ? n : bits;
        sc_value *a = random_poly(ctx, r, middle ? 2 * n - 1 : n, b, state);
        sc_value *c = random_poly(ctx, r, n, b, state);
        double mad = 0.0, q;

        if (a == NULL || c == NULL) {
            fprintf(stderr, "out of memory while tuning short product\n");
            exit(1);
        }
        q = short_ratio(ctx, old, new, a, c, n, &mad);
        printf("  n=%-6zu bits=%-6zu new/old=%7.4f  MAD=%5.2f%%\n",
               n, b, q, 100.0 * mad);
        sc_value_free_many(2, a, c);
        if (q >= 1.05) {
            tr.loss = n;
            tr.loss_ratio = q;
        } else if (q <= 0.95) {
            tr.win = n;
            tr.win_ratio = q;
            tr.cut = tr.loss != 0 ? tr.loss + (n - tr.loss) / 2 : n;
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
    } else if (tr.loss == 0)
        printf("  already >5%% faster at first point; cutoff <= %zu\n", tr.cut);
    else
        printf("  5%% bracket [%zu, %zu], midpoint %zu\n", tr.loss, tr.win, tr.cut);
    return tr;
}

static int make_division_case(sc_context *ctx, sc_parent *r,
                              gmp_randstate_t state, size_t n, size_t bits,
                              int balanced, int unit_lead, int exact,
                              sc_value **ap, sc_value **bp)
{
    size_t bn = balanced ? n : (n + 1) / 2, rn = exact || bn < 2 ? 0 : bn - 1;
    sc_value *b = random_poly(ctx, r, bn, bits, state);
    sc_value *q = random_poly(ctx, r, n, bits, state);
    sc_value *rem = rn ? random_poly(ctx, r, rn, bits, state) : NULL;
    sc_value *a = b && q ? sc_zz_poly_mul_classical(ctx, b, q) : NULL;
    sc_value *sum = a && rem ? sc_zz_poly_add(ctx, a, rem) : NULL;

    if (unit_lead && b != NULL)
        mpz_set_ui(b->data.zz_poly.coeff[bn - 1], 1);
    if (unit_lead) {
        sc_value_free_many(2, a, sum);
        a = b && q ? sc_zz_poly_mul_classical(ctx, b, q) : NULL;
        sum = a && rem ? sc_zz_poly_add(ctx, a, rem) : NULL;
    }
    if (rem != NULL) {
        sc_value_free(a);
        a = sum;
        sum = NULL;
    }
    sc_value_free_many(3, q, rem, sum);
    if (a == NULL || b == NULL) {
        sc_value_free_many(2, a, b);
        return 0;
    }
    *ap = a;
    *bp = b;
    return 1;
}

static int make_pseudodiv_case(sc_context *ctx, sc_parent *r,
                               gmp_randstate_t state, size_t n, size_t bits,
                               sc_value **ap, sc_value **bp)
{
    *ap = random_poly(ctx, r, 2 * n - 1, bits, state);
    *bp = random_poly(ctx, r, n, bits, state);
    if (*ap != NULL && *bp != NULL)
        return 1;
    sc_value_free_many(2, *ap, *bp);
    return 0;
}

static int make_series_case(sc_context *ctx, sc_parent *r,
                            gmp_randstate_t state, size_t n, size_t bits,
                            int unit, sc_value **ap, sc_value **bp)
{
    sc_value *b = random_poly(ctx, r, n, bits, state);
    sc_value *q = random_poly(ctx, r, n, bits, state), *a;

    if (b == NULL || q == NULL)
        return sc_value_free_many_null(2, b, q), 0;
    if (unit)
        mpz_set_ui(b->data.zz_poly.coeff[0], 1);
    else
        mpz_set_ui(b->data.zz_poly.coeff[0], 3);
    a = sc_zz_poly_mullow_classical(ctx, b, q, n);
    sc_value_free(q);
    if (a == NULL)
        return sc_value_free_many_null(1, b), 0;
    *ap = a;
    *bp = b;
    return 1;
}

static tune_result tune_div_pair(sc_context *ctx, sc_parent *r,
                                 gmp_randstate_t state, const char *name,
                                 mul_fn old, mul_fn new, size_t bits,
                                 int balanced, int unit_lead, int exact,
                                 size_t lo, size_t hi, size_t fallback)
{
    tune_result tr = { 0, 0, fallback, 0.0, 0.0 };
    size_t streak = 0, first_win = 0, points = 0;

    printf("\n%s\n", name);
    printf("  bounded search: n <= %zu, at most %u sampled sizes\n", hi, DIV_TUNE_MAX_POINTS);
    for (size_t n = lo; n <= hi && points < DIV_TUNE_MAX_POINTS; points++) {
        sc_value *a = NULL, *b = NULL;
        double mad = 0.0, q;

        if (!make_division_case(ctx, r, state, n, bits, balanced,
                                unit_lead, exact, &a, &b)) {
            fprintf(stderr, "out of memory while tuning division\n");
            exit(1);
        }
        q = ratio_bounded(ctx, old, new, a, b, &mad);
        printf("  qn=%-6zu bits=%-4zu new/old=%7.4f  MAD=%5.2f%%\n",
               n, bits, q, 100.0 * mad);
        sc_value_free_many(2, a, b);
        if (q >= 1.05) {
            tr.loss = n;
            tr.loss_ratio = q;
            streak = 0;
        } else if (q <= 0.95 && tr.loss != 0) {
            if (streak++ == 0)
                first_win = n;
            if (streak >= 2) {
                tr.win = first_win;
                tr.win_ratio = q;
                tr.cut = tr.loss ? tr.loss + (first_win - tr.loss) / 2 : first_win;
                break;
            }
        } else
            streak = 0;
        if (n == hi)
            break;
        n = next_size(n) > hi ? hi : next_size(n);
    }
    if (tr.win == 0 && points >= DIV_TUNE_MAX_POINTS)
        puts("  sample-count bound reached before a sustained crossover");
    if (tr.win == 0) {
        if (fallback == (size_t)-1)
            puts("  no sustained 5% win found; keeping disabled");
        else
            printf("  no sustained 5%% win found; keeping %zu\n", fallback);
    }
    else if (tr.loss == 0)
        printf("  already >5%% faster at first point; cutoff <= %zu\n", tr.cut);
    else
        printf("  5%% bracket [%zu, %zu], midpoint %zu\n", tr.loss, tr.win, tr.cut);
    return tr;
}

static tune_result tune_pseudo_pair(sc_context *ctx, sc_parent *r,
                                    gmp_randstate_t state, const char *name,
                                    mul_fn old, mul_fn new, size_t bits,
                                    size_t lo, size_t hi, size_t fallback)
{
    tune_result tr = { 0, 0, fallback, 0.0, 0.0 };
    size_t streak = 0, first_win = 0, points = 0;

    printf("\n%s\n", name);
    printf("  bounded search: n <= %zu, at most %u sampled sizes\n", hi, DIV_TUNE_MAX_POINTS);
    for (size_t n = lo; n <= hi && points < DIV_TUNE_MAX_POINTS; points++) {
        sc_value *a = NULL, *b = NULL;
        double mad = 0.0, q;

        if (!make_pseudodiv_case(ctx, r, state, n, bits, &a, &b)) {
            fprintf(stderr, "out of memory while tuning pseudo-division\n");
            exit(1);
        }
        q = ratio_bounded(ctx, old, new, a, b, &mad);
        printf("  qn=%-6zu bits=%-4zu new/old=%7.4f  MAD=%5.2f%%\n",
               n, bits, q, 100.0 * mad);
        sc_value_free_many(2, a, b);
        if (q >= 1.05)
            tr.loss = n, tr.loss_ratio = q, streak = 0;
        else if (q <= 0.95 && tr.loss != 0) {
            if (streak++ == 0)
                first_win = n;
            if (streak >= 2) {
                tr.win = first_win, tr.win_ratio = q;
                tr.cut = tr.loss ? tr.loss + (first_win - tr.loss) / 2 : first_win;
                break;
            }
        } else
            streak = 0;
        if (n == hi)
            break;
        n = next_size(n) > hi ? hi : next_size(n);
    }
    if (tr.win == 0 && points >= DIV_TUNE_MAX_POINTS)
        puts("  sample-count bound reached before a sustained crossover");
    if (tr.win == 0) {
        if (fallback == (size_t)-1)
            puts("  no sustained 5% win found; keeping disabled");
        else
            printf("  no sustained 5%% win found; keeping %zu\n", fallback);
    }
    else if (tr.loss == 0)
        printf("  already >5%% faster at first point; cutoff <= %zu\n", tr.cut);
    else
        printf("  5%% bracket [%zu, %zu], midpoint %zu\n", tr.loss, tr.win, tr.cut);
    return tr;
}

static sc_value *inv_series_classical_tune(sc_context *ctx, const sc_value *a,
                                            const sc_value *b, size_t n)
{
    (void)a;
    return sc_zz_poly_inv_series_classical(ctx, b, n);
}

static sc_value *inv_series_newton_tune(sc_context *ctx, const sc_value *a,
                                         const sc_value *b, size_t n)
{
    (void)a;
    return sc_zz_poly_inv_series_newton(ctx, b, n);
}

static tune_result tune_series_pair(sc_context *ctx, sc_parent *r,
                                    gmp_randstate_t state, const char *name,
                                    short_fn old, short_fn new, size_t bits,
                                    int unit, size_t lo, size_t hi,
                                    size_t fallback)
{
    tune_result tr = { 0, 0, fallback, 0.0, 0.0 };
    size_t streak = 0, first_win = 0, points = 0;

    printf("\n%s\n", name);
    printf("  bounded search: n <= %zu, at most %u sampled sizes\n", hi, DIV_TUNE_MAX_POINTS);
    for (size_t n = lo; n <= hi && points < DIV_TUNE_MAX_POINTS; points++) {
        sc_value *a = NULL, *b = NULL;
        double mad = 0.0, q;

        if (!make_series_case(ctx, r, state, n, bits, unit, &a, &b)) {
            fprintf(stderr, "out of memory while tuning series division\n");
            exit(1);
        }
        q = short_ratio_bounded(ctx, old, new, a, b, n, &mad);
        printf("  n=%-6zu bits=%-4zu new/old=%7.4f  MAD=%5.2f%%\n",
               n, bits, q, 100.0 * mad);
        sc_value_free_many(2, a, b);
        if (q >= 1.05)
            tr.loss = n, tr.loss_ratio = q, streak = 0;
        else if (q <= 0.95 && tr.loss != 0) {
            if (streak++ == 0)
                first_win = n;
            if (streak >= 2) {
                tr.win = first_win, tr.win_ratio = q;
                tr.cut = tr.loss ? tr.loss + (first_win - tr.loss) / 2 : first_win;
                break;
            }
        } else
            streak = 0;
        if (n == hi)
            break;
        n = next_size(n) > hi ? hi : next_size(n);
    }
    if (tr.win == 0 && points >= DIV_TUNE_MAX_POINTS)
        puts("  sample-count bound reached before a sustained crossover");
    if (tr.win == 0) {
        if (fallback == (size_t)-1)
            puts("  no sustained 5% win found; keeping disabled");
        else
            printf("  no sustained 5%% win found; keeping %zu\n", fallback);
    }
    else if (tr.loss == 0)
        printf("  already >5%% faster at first point; cutoff <= %zu\n", tr.cut);
    else
        printf("  5%% bracket [%zu, %zu], midpoint %zu\n", tr.loss, tr.win, tr.cut);
    return tr;
}

static sc_value *divexact_bidir_tune(sc_context *ctx, const sc_value *a,
                                     const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t qn = an - bn + 1, ln = (qn + 1) / 2, midn = bn - 1;
    sc_value *q = sc_zz_poly_quo_bidirectional(ctx, a, b);
    sc_value *mid = q ? sc_zz_poly_mulmid(ctx, b, q, ln, midn) : NULL;

    if (q == NULL || mid == NULL)
        return sc_value_free_many_null(2, q, mid);
    sc_value_free(mid);
    return q;
}

static division_tuning tune_division(sc_context *ctx, sc_parent *r,
                                     gmp_randstate_t state)
{
    division_tuning d;
    size_t series_dc0 = sc_tune_series_quo_dc_cutoff;
    size_t bidir0 = sc_tune_bidir_quo_cutoff;
    size_t mulders0 = sc_tune_mulders_quo_cutoff;
    size_t quo_dc0 = sc_tune_quo_dc_cutoff;
    size_t divrem_dc0 = sc_tune_divrem_dc_cutoff;
    size_t divexact0 = sc_tune_divexact_bidir_cutoff;
    size_t lo;

    memset(&d, 0, sizeof(d));
    puts("\nDivision tuning: bounded 5% loss -> two-win crossover searches.");
    sc_tune_series_quo_dc_cutoff = 1;
    d.series_dc = tune_series_pair(ctx, r, state,
        "series quotient classical -> divide-and-conquer (256-bit coefficients)",
        sc_zz_poly_series_quo_classical, sc_zz_poly_series_quo_dc,
        256, 0, 8, 768, series_dc0);
    sc_tune_series_quo_dc_cutoff = d.series_dc.cut;

    sc_tune_inv_series_newton_cutoff = 1;
    d.inverse_newton = tune_series_pair(ctx, r, state,
        "series inverse classical -> Newton (8-bit coefficients)",
        inv_series_classical_tune, inv_series_newton_tune,
        8, 1, 8, 512, (size_t)-1);
    sc_tune_inv_series_newton_cutoff = d.inverse_newton.cut;

    sc_tune_series_quo_newton_cutoff = 1;
    d.series_newton = tune_series_pair(ctx, r, state,
        "series quotient classical -> Karp-Markstein (8-bit, unit constant)",
        sc_zz_poly_series_quo_classical, sc_zz_poly_series_quo_newton,
        8, 1, 8, 512, (size_t)-1);
    sc_tune_series_quo_newton_cutoff = d.series_newton.cut;

    d.quo_dc = tune_div_pair(ctx, r, state,
        "ordinary quotient classical -> divide-and-conquer (256-bit, q about 2x divisor)",
        sc_zz_poly_quo_classical, sc_zz_poly_quo_dc,
        256, 0, 0, 0, 8, 768, quo_dc0);
    sc_tune_quo_dc_cutoff = d.quo_dc.cut;
    d.divrem_dc = tune_div_pair(ctx, r, state,
        "divrem classical -> divide-and-conquer (256-bit, q about 2x divisor)",
        sc_zz_poly_divrem_classical, sc_zz_poly_divrem_dc,
        256, 0, 0, 0, 8, 768, divrem_dc0);
    sc_tune_divrem_dc_cutoff = d.divrem_dc.cut;

    sc_tune_mulders_quo_cutoff = 1;
    d.mulders_base = tune_div_pair(ctx, r, state,
        "balanced quotient classical -> recursive Mulders (256-bit coefficients)",
        sc_zz_poly_quo_classical, sc_zz_poly_quo_mulders,
        256, 1, 0, 0, 8, 1024, mulders0);
    sc_tune_mulders_quo_cutoff = d.mulders_base.cut;

    sc_tune_quo_mulders_cutoff = (size_t)-1;
    sc_tune_quo_newton_cutoff = (size_t)-1;
    lo = d.mulders_base.cut == (size_t)-1 ? 0 : d.mulders_base.cut + 1;
    if (lo != 0 && lo <= 1024)
        d.quo_mulders = tune_div_pair(ctx, r, state,
            "balanced lower quotient chain -> Mulders (256-bit coefficients)",
            sc_zz_poly_quo, sc_zz_poly_quo_mulders,
            256, 1, 0, 0, lo, 1024, (size_t)-1);
    else
        d.quo_mulders.cut = (size_t)-1;
    sc_tune_quo_mulders_cutoff = d.quo_mulders.cut;

    sc_tune_divrem_mulders_cutoff = (size_t)-1;
    sc_tune_divrem_newton_cutoff = (size_t)-1;
    if (lo != 0 && lo <= 1024)
        d.divrem_mulders = tune_div_pair(ctx, r, state,
            "balanced lower divrem chain -> Mulders (256-bit coefficients)",
            sc_zz_poly_divrem, sc_zz_poly_divrem_mulders,
            256, 1, 0, 0, lo, 1024, (size_t)-1);
    else
        d.divrem_mulders.cut = (size_t)-1;
    sc_tune_divrem_mulders_cutoff = d.divrem_mulders.cut;

    sc_tune_quo_newton_cutoff = (size_t)-1;
    d.quo_newton = tune_div_pair(ctx, r, state,
        "unit-leading lower quotient chain -> Newton (8-bit coefficients)",
        sc_zz_poly_quo, sc_zz_poly_quo_newton,
        8, 1, 1, 0, 32, 512, (size_t)-1);
    sc_tune_quo_newton_cutoff = d.quo_newton.cut;

    sc_tune_divrem_newton_cutoff = (size_t)-1;
    d.divrem_newton = tune_div_pair(ctx, r, state,
        "unit-leading lower divrem chain -> Newton (8-bit coefficients)",
        sc_zz_poly_divrem, sc_zz_poly_divrem_newton,
        8, 1, 1, 0, 32, 512, (size_t)-1);
    sc_tune_divrem_newton_cutoff = d.divrem_newton.cut;

    sc_tune_bidir_quo_cutoff = 1;
    d.bidir_base = tune_div_pair(ctx, r, state,
        "exact balanced quotient classical -> bidirectional core (256-bit coefficients)",
        sc_zz_poly_quo_classical, sc_zz_poly_quo_bidirectional,
        256, 1, 0, 1, 8, 1024, bidir0);
    sc_tune_bidir_quo_cutoff = d.bidir_base.cut;

    sc_tune_divexact_bidir_cutoff = (size_t)-1;
    d.divexact_bidir = tune_div_pair(ctx, r, state,
        "exact division lower chain -> bidirectional + middle check (256-bit coefficients)",
        sc_zz_poly_divexact, divexact_bidir_tune,
        256, 1, 0, 1, 8, 1024, divexact0);
    sc_tune_divexact_bidir_cutoff = d.divexact_bidir.cut;

    d.pseudodiv_fast = tune_pseudo_pair(ctx, r, state,
        "pseudo-division classical -> scaled-Newton fast path (32-bit coefficients)",
        sc_zz_poly_pseudodiv_impl, sc_zz_poly_pseudodiv_fast,
        32, 8, 256, (size_t)-1);
    sc_tune_pseudodiv_fast_cutoff = d.pseudodiv_fast.cut;
    d.pseudorem_fast = tune_pseudo_pair(ctx, r, state,
        "pseudo-remainder classical -> fast pseudo-division remainder (32-bit coefficients)",
        sc_zz_poly_pseudorem_classical, sc_zz_poly_pseudorem_fast,
        32, 8, 256, (size_t)-1);
    sc_tune_pseudorem_fast_cutoff = d.pseudorem_fast.cut;
    return d;
}


static gcd_tuning tune_gcd(sc_context *ctx, sc_parent *r, gmp_randstate_t state)
{
    gcd_tuning g = { 0 };
    tune_result *tr = &g.subresultant;
    const size_t lo = 4, hi = 64, fallback = sc_tune_gcd_subresultant_cutoff;
    size_t streak = 0, first_win = 0, points = 0;

    tr->cut = fallback;
    puts("\nGCD tuning: small primitive pseudo-Euclidean -> Brown subresultant.");
    printf("  bounded search: n <= %zu, at most %u sampled sizes\n",
           hi, DIV_TUNE_MAX_POINTS);
    for (size_t n = lo; n <= hi && points < DIV_TUNE_MAX_POINTS; points++) {
        sc_value *a = random_poly(ctx, r, n, 32, state);
        sc_value *b = random_poly(ctx, r, n > 1 ? n - 1 : 1, 32, state);
        double mad = 0.0, q;

        if (a == NULL || b == NULL) {
            fprintf(stderr, "out of memory while tuning gcd\n");
            exit(1);
        }
        q = ratio_bounded(ctx, sc_zz_poly_gcd_pseudo_impl,
                          sc_zz_poly_gcd_subresultant_impl, a, b, &mad);
        printf("  n=%-7zu subres/pseudo=%7.4f  MAD=%5.2f%%\n",
               n, q, 100.0 * mad);
        sc_value_free_many(2, a, b);
        if (q >= 1.05) {
            tr->loss = n;
            tr->loss_ratio = q;
            streak = 0;
        } else if (q <= 0.95) {
            if (tr->loss == 0) {
                tr->win = n;
                tr->win_ratio = q;
                tr->cut = n;
                break;
            }
            if (streak++ == 0)
                first_win = n;
            if (streak >= 2) {
                tr->win = first_win;
                tr->win_ratio = q;
                tr->cut = tr->loss + (first_win - tr->loss + 1) / 2;
                break;
            }
        } else
            streak = 0;
        if (n == hi)
            break;
        n = next_size(n) > hi ? hi : next_size(n);
    }
    if (tr->win == 0 && tr->loss != 0) {
        tr->cut = hi + 1;
        printf("  no Brown win through %zu; using pseudo gcd only below bounded cutoff %zu\n",
               hi, tr->cut);
    } else if (tr->win == 0) {
        printf("  no stable 5%% preference; keeping cutoff %zu\n", tr->cut);
    } else if (tr->loss == 0) {
        printf("  Brown already >5%% faster; cutoff %zu\n", tr->cut);
    } else {
        printf("  5%% bracket [%zu, %zu], midpoint %zu\n",
               tr->loss, tr->win, tr->cut);
    }
    sc_tune_gcd_subresultant_cutoff = tr->cut;
    return g;
}


static size_t ntt_cutoff_full(sc_context *ctx, sc_parent *r, gmp_randstate_t state,
                              const tune_result *tr)
{
    sc_value *a, *b;
    size_t np, cut = tr->cut;

    if (tr->win == 0)
        return (size_t)-1;
    a = random_poly(ctx, r, tr->cut, 48, state);
    b = random_poly(ctx, r, tr->cut, 48, state);
    np = a && b ? sc_zz_poly_ntt_nprimes(a, b) : 0;
    if (np != 0)
        cut /= np;
    sc_value_free_many(2, a, b);
    return cut;
}

static size_t ntt_cutoff_mid(sc_context *ctx, sc_parent *r, gmp_randstate_t state,
                             const tune_result *tr)
{
    sc_value *a, *b;
    size_t np, cut = tr->cut;

    if (tr->win == 0)
        return (size_t)-1;
    a = random_poly(ctx, r, 2 * tr->cut - 1, 48, state);
    b = random_poly(ctx, r, tr->cut, 48, state);
    np = a && b ? sc_zz_poly_mulmid_ntt_nprimes(a, b, tr->cut) : 0;
    if (np != 0)
        cut /= np;
    sc_value_free_many(2, a, b);
    return cut;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "PolynomialRing(ZZ)", SC_PARENT_POLY, &SC_ZZ, "x" };
    gmp_randstate_t state;
    tune_result ks, kar, low, toom, ntt, ssa, lowdc;
    tune_result low_ntt_r, low_ssa_r, high_ntt_r, high_ssa_r;
    tune_result midclass, mid63, mid_ntt_r, mid_ssa_r;
    division_tuning div;
    gcd_tuning gcd;
    size_t ntt_cut, toom_fallback, low_ntt_cut, low_ssa_cut;
    size_t high_ntt_cut, high_ssa_cut, mid_ntt_cut, mid_ssa_cut, mid63_lo;
    unsigned mfa_cut;

    sc_context_init(&ctx);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 20260923);
    sc_tune_mul_ntt_cutoff = sc_tune_mul_ssa_cutoff = (size_t)-1;
    sc_tune_mullow_ntt_cutoff = sc_tune_mullow_ssa_cutoff = (size_t)-1;
    sc_tune_mulhigh_ntt_cutoff = sc_tune_mulhigh_ssa_cutoff = (size_t)-1;
    sc_tune_mulmid_ntt_cutoff = sc_tune_mulmid_ssa_cutoff = (size_t)-1;

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
    ntt_cut = ntt_cutoff_full(&ctx, &r, state, &ntt);
    mfa_cut = tune_mfa();
    sc_tune_ssa_mfa_cutoff_log = mfa_cut;
    ssa = tune_pair(&ctx, &r, state, "lower dispatcher -> SSA (bits = length)",
                    sc_zz_poly_mul, sc_zz_poly_mul_ssa, 0, 1, 32, 2048, (size_t)-1);
    sc_tune_mul_ntt_cutoff = ntt_cut;
    sc_tune_mul_ssa_cutoff = ssa.win ? ssa.cut : (size_t)-1;

    puts("\nLow/high product tuning: lower short algorithm -> full FFT product and slice.");
    lowdc = tune_short_pair(&ctx, &r, state,
                            "mullo classical -> divide-and-conquer (256-bit coefficients)",
                            sc_zz_poly_mullow_classical, sc_zz_poly_mullow_dc,
                            256, 0, 0, 8, 256, sc_tune_mullow_dc_cutoff);
    sc_tune_mullow_dc_cutoff = lowdc.cut;
    low_ntt_r = tune_short_pair(&ctx, &r, state,
                                "mullo lower dispatcher -> full CRT-NTT (48-bit coefficients)",
                                mullow_dispatch, mullow_ntt, 48, 0, 0,
                                96, 16384, (size_t)-1);
    low_ntt_cut = ntt_cutoff_full(&ctx, &r, state, &low_ntt_r);
    sc_tune_mullow_ntt_cutoff = low_ntt_cut;
    low_ssa_r = tune_short_pair(&ctx, &r, state,
                                "mullo lower dispatcher -> full SSA (bits = length)",
                                mullow_dispatch, mullow_ssa, 0, 1, 0,
                                32, 2048, (size_t)-1);
    low_ssa_cut = low_ssa_r.win ? low_ssa_r.cut : (size_t)-1;
    sc_tune_mullow_ssa_cutoff = low_ssa_cut;
    high_ntt_r = tune_short_pair(&ctx, &r, state,
                                 "mulhi reversed mullo -> direct full CRT-NTT "
                                 "(48-bit coefficients)",
                                 mulhigh_dispatch, mulhigh_ntt, 48, 0, 0,
                                 96, 16384, (size_t)-1);
    high_ntt_cut = ntt_cutoff_full(&ctx, &r, state, &high_ntt_r);
    sc_tune_mulhigh_ntt_cutoff = high_ntt_cut;
    high_ssa_r = tune_short_pair(&ctx, &r, state,
                                 "mulhi reversed mullo -> direct full SSA (bits = length)",
                                 mulhigh_dispatch, mulhigh_ssa, 0, 1, 0,
                                 32, 2048, (size_t)-1);
    high_ssa_cut = high_ssa_r.win ? high_ssa_r.cut : (size_t)-1;
    sc_tune_mulhigh_ssa_cutoff = high_ssa_cut;

    puts("\nMiddle product tuning: classical -> Toom-4/2 -> Toom-6/3 -> FFT wraparound.");
    toom_fallback = sc_tune_mulmid_toom63_cutoff;
    sc_tune_mulmid_toom63_cutoff = (size_t)-1;
    midclass = tune_short_pair(&ctx, &r, state,
                               "mulmid classical -> Toom-4/2 (256-bit coefficients)",
                               mulmid_classical, mulmid_toom42, 256, 0, 1,
                               8, 256, sc_tune_mulmid_classical_cutoff);
    sc_tune_mulmid_classical_cutoff = midclass.cut;
    mid63_lo = midclass.cut < 24 ? 24 : midclass.cut + 1;
    mid63 = tune_short_pair(&ctx, &r, state,
                            "mulmid Toom-4/2 chain -> Toom-6/3 (256-bit coefficients)",
                            mulmid_dispatch, mulmid_toom63, 256, 0, 1,
                            mid63_lo, 768, toom_fallback);
    sc_tune_mulmid_toom63_cutoff = mid63.cut;
    mid_ntt_r = tune_short_pair(&ctx, &r, state,
                                "mulmid lower dispatcher -> CRT-NTT wraparound "
                                "(48-bit coefficients)",
                                mulmid_dispatch, sc_zz_poly_mulmid_ntt, 48, 0, 1,
                                64, 16384, (size_t)-1);
    mid_ntt_cut = ntt_cutoff_mid(&ctx, &r, state, &mid_ntt_r);
    sc_tune_mulmid_ntt_cutoff = mid_ntt_cut;
    mid_ssa_r = tune_short_pair(&ctx, &r, state,
                                "mulmid lower dispatcher -> SSA wraparound (bits = length)",
                                mulmid_dispatch, sc_zz_poly_mulmid_ssa, 0, 1, 1,
                                32, 2048, (size_t)-1);
    mid_ssa_cut = mid_ssa_r.win ? mid_ssa_r.cut : (size_t)-1;
    sc_tune_mulmid_ssa_cutoff = mid_ssa_cut;

    div = tune_division(&ctx, &r, state);
    gcd = tune_gcd(&ctx, &r, state);

    if (!write_tuning(&ks, &toom, &kar, &low, &ntt, ntt_cut, &ssa, mfa_cut,
                      &lowdc, low_ntt_cut, low_ssa_cut, high_ntt_cut, high_ssa_cut,
                      &midclass, &mid63, mid_ntt_cut, mid_ssa_cut, &div, &gcd)) {
        gmp_randclear(state);
        sc_context_clear(&ctx);
        return 1;
    }

    gmp_randclear(state);
    sc_context_clear(&ctx);
    return 0;
}
