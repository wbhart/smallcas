#include "smallcas.h"

#include <stdio.h>

static sc_value *poly_si(sc_context *ctx, sc_parent *r, size_t n, const long *c)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);
    size_t i;

    if (f == NULL)
        return NULL;
    for (i = 0; i < n; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], c[i]);
    sc_zz_poly_normalize(f);
    return f;
}

static int poly_equal(const sc_value *a, const sc_value *b)
{
    size_t i;

    if (a == NULL || b == NULL || a->data.zz_poly.length != b->data.zz_poly.length)
        return 0;
    for (i = 0; i < a->data.zz_poly.length; i++)
        if (mpz_cmp(a->data.zz_poly.coeff[i], b->data.zz_poly.coeff[i]) != 0)
            return 0;
    return 1;
}

static int zz_si(const sc_value *a, long n)
{
    return a != NULL && a->kind == SC_VALUE_ZZ && mpz_cmp_si(a->data.z, n) == 0;
}

static int test_accessors(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { 5, 0, -9, 0, 17 };
    sc_value *f = poly_si(ctx, r, 5, fc), *z = sc_value_new_zz_poly_checked(ctx, r, 0);
    sc_value *d, *dz, *lc, *lcz, *cc, *c2, *c9, *h, *bits;
    int ok;

    if (f == NULL || z == NULL)
        return 0;
    d = sc_poly_degree(ctx, f);
    dz = sc_poly_degree(ctx, z);
    lc = sc_poly_leading_coefficient(ctx, f);
    lcz = sc_poly_leading_coefficient(ctx, z);
    cc = sc_poly_constant_coefficient(ctx, f);
    c2 = sc_poly_coeff(ctx, f, 2);
    c9 = sc_poly_coeff(ctx, f, 9);
    h = sc_poly_height(ctx, f);
    bits = sc_poly_max_abs_bits(ctx, f);
    ok = zz_si(d, 4) && zz_si(dz, -1) && zz_si(lc, 17) && zz_si(lcz, 0) &&
         zz_si(cc, 5) && zz_si(c2, -9) && zz_si(c9, 0) && zz_si(h, 17) &&
         zz_si(bits, 5);
    sc_value_free_many(11, f, z, d, dz, lc, lcz, cc, c2, c9, h, bits);
    return ok;
}

static int test_rearrangements(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { 5, 0, -9, 0, 17 }, revc[] = { 17, 0, -9, 0, 5 };
    const long rev7c[] = { 0, 0, 17, 0, -9, 0, 5 }, trc[] = { 5, 0, -9 };
    const long slc[] = { 0, 0, 5, 0, -9, 0, 17 }, src[] = { -9, 0, 17 };
    sc_value *f = poly_si(ctx, r, 5, fc), *er = poly_si(ctx, r, 5, revc);
    sc_value *er7 = poly_si(ctx, r, 7, rev7c), *et = poly_si(ctx, r, 3, trc);
    sc_value *esl = poly_si(ctx, r, 7, slc), *esr = poly_si(ctx, r, 3, src);
    sc_value *rev, *rev7, *tr, *sl, *sr, *sr5;
    int ok;

    if (f == NULL || er == NULL || er7 == NULL || et == NULL || esl == NULL || esr == NULL)
        return 0;
    rev = sc_poly_reverse(ctx, f, f->data.zz_poly.length);
    rev7 = sc_poly_reverse(ctx, f, 7);
    tr = sc_poly_truncate(ctx, f, 4);
    sl = sc_poly_shift_left(ctx, f, 2);
    sr = sc_poly_shift_right(ctx, f, 2);
    sr5 = sc_poly_shift_right(ctx, f, 5);
    ok = poly_equal(rev, er) && poly_equal(rev7, er7) && poly_equal(tr, et) &&
         poly_equal(sl, esl) && poly_equal(sr, esr) && sr5 != NULL &&
         sr5->data.zz_poly.length == 0;
    sc_value_free_many(12, f, er, er7, et, esl, esr, rev, rev7, tr, sl, sr, sr5);
    return ok;
}

static int test_inflation(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { 5, 0, -9, 0, 17 }, dc[] = { 5, -9, 17 };
    const long ic[] = { 5, 0, 0, 0, -9, 0, 0, 0, 17 }, c13[] = { 13 };
    const long mc[] = { 0, 0, 0, 0, 0, 0, 3 }, kc[] = { 7 };
    sc_value *f = poly_si(ctx, r, 5, fc), *ed = poly_si(ctx, r, 3, dc);
    sc_value *ei = poly_si(ctx, r, 9, ic), *e13 = poly_si(ctx, r, 1, c13);
    sc_value *m = poly_si(ctx, r, 7, mc), *k = poly_si(ctx, r, 1, kc);
    sc_value *z = sc_value_new_zz_poly_checked(ctx, r, 0);
    sc_value *infl, *infl4, *infl0, *defl, *round, *df, *dm, *dk, *dz;
    int ok;

    if (f == NULL || ed == NULL || ei == NULL || e13 == NULL || m == NULL ||
        k == NULL || z == NULL)
        return 0;
    infl = sc_poly_inflate(ctx, ed, 2);
    infl4 = sc_poly_inflate(ctx, ed, 4);
    infl0 = sc_poly_inflate(ctx, f, 0);
    defl = sc_poly_deflate(ctx, f, 2);
    round = defl ? sc_poly_inflate(ctx, defl, 2) : NULL;
    df = sc_poly_deflation(ctx, f);
    dm = sc_poly_deflation(ctx, m);
    dk = sc_poly_deflation(ctx, k);
    dz = sc_poly_deflation(ctx, z);
    ok = poly_equal(infl, f) && poly_equal(infl4, ei) && poly_equal(infl0, e13) &&
         poly_equal(defl, ed) && poly_equal(round, f) && zz_si(df, 2) && zz_si(dm, 6) &&
         zz_si(dk, 1) && zz_si(dz, 0);
    sc_value_free_many(16, f, ed, ei, e13, m, k, z, infl, infl4, infl0, defl, round, df,
                       dm, dk, dz);
    return ok;
}

static int test_named(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { 5, 0, -9, 0, 17 }, dc[] = { 5, -9, 17 };
    sc_value *f = poly_si(ctx, r, 5, fc), *ed = poly_si(ctx, r, 3, dc);
    sc_value *two = sc_zz_from_str(ctx, "2"), *seven = sc_zz_from_str(ctx, "7");
    sc_value *d, *lc, *cc, *c2, *rev, *tr, *sl, *sr, *h, *bits, *inf, *df, *defl;
    int ok;

    if (f == NULL || ed == NULL || two == NULL || seven == NULL)
        return 0;
    d = sc_call1(ctx, "degree", f);
    lc = sc_call1(ctx, "leading_coefficient", f);
    cc = sc_call1(ctx, "constant_coefficient", f);
    c2 = sc_call2(ctx, "coeff", f, two);
    rev = sc_call1(ctx, "reverse", f);
    tr = sc_call2(ctx, "truncate", f, seven);
    sl = sc_call2(ctx, "shift_left", ed, two);
    sr = sc_call2(ctx, "shift_right", f, two);
    h = sc_call1(ctx, "height", f);
    bits = sc_call1(ctx, "max_abs_bits", f);
    inf = sc_call2(ctx, "inflate", ed, two);
    df = sc_call1(ctx, "maximal_deflation", f);
    defl = sc_call2(ctx, "deflate", f, two);
    ok = zz_si(d, 4) && zz_si(lc, 17) && zz_si(cc, 5) && zz_si(c2, -9) &&
         rev != NULL && poly_equal(tr, f) && sl != NULL && sr != NULL && zz_si(h, 17) &&
         zz_si(bits, 5) && poly_equal(inf, f) && zz_si(df, 2) && poly_equal(defl, ed);
    sc_value_free_many(17, f, ed, two, seven, d, lc, cc, c2, rev, tr, sl, sr, h, bits,
                       inf, df, defl);
    return ok;
}

static int test_invalid_deflation(sc_parent *r)
{
    const long fc[] = { 1, 1 };
    sc_context ctx;
    sc_value *f, *bad;
    int ok;

    sc_context_init(&ctx);
    f = poly_si(&ctx, r, 2, fc);
    bad = f ? sc_poly_deflate(&ctx, f, 2) : NULL;
    ok = bad == NULL && ctx.error[0] != '\0';
    sc_value_free_many(2, f, bad);
    sc_context_clear(&ctx);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "ZZ[x]", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_accessors(&ctx, &r) || !test_rearrangements(&ctx, &r) ||
        !test_inflation(&ctx, &r) || !test_named(&ctx, &r) || !test_invalid_deflation(&r)) {
        fprintf(stderr, "structural test failed: %s\n", ctx.error);
        sc_context_clear(&ctx);
        return 1;
    }
    sc_context_clear(&ctx);
    puts("structural tests passed");
    return 0;
}
