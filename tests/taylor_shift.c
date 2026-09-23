#include "smallcas.h"

#include <stdio.h>
#include <stdlib.h>

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

static sc_value *random_poly(sc_context *ctx, sc_parent *r, size_t n, long bound)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);
    size_t i;

    if (f == NULL)
        return NULL;
    for (i = 0; i < n; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], (rand() % (2 * bound + 1)) - bound);
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

static int test_fixed(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { -7, 5, -2, 3 }, ec[] = { 19, 33, 16, 3 };
    sc_value *f = poly_si(ctx, r, 4, fc), *e = poly_si(ctx, r, 4, ec);
    sc_value *c = sc_zz_from_str(ctx, "2"), *h, *d, *a;
    int ok;

    if (f == NULL || e == NULL || c == NULL)
        return 0;
    h = sc_zz_poly_taylor_shift_horner(ctx, f, c);
    d = sc_zz_poly_taylor_shift_divconquer(ctx, f, c);
    a = sc_zz_poly_taylor_shift(ctx, f, c);
    ok = poly_equal(h, e) && poly_equal(d, e) && poly_equal(a, e);
    sc_value_free_many(6, f, e, c, h, d, a);
    return ok;
}

static int test_random(sc_context *ctx, sc_parent *r)
{
    size_t n, trial;

    srand(92329);
    for (n = 0; n <= 48; n++)
        for (trial = 0; trial < 3; trial++) {
            sc_value *f = random_poly(ctx, r, n, 9), *c = sc_value_new_zz_checked(ctx);
            sc_value *g = sc_value_new_zz_poly_checked(ctx, r, 2), *h, *d, *a, *q;
            int ok;

            if (f == NULL || c == NULL || g == NULL)
                return 0;
            mpz_set_si(c->data.z, (rand() % 11) - 5);
            mpz_set(g->data.zz_poly.coeff[0], c->data.z);
            mpz_set_ui(g->data.zz_poly.coeff[1], 1);
            h = sc_zz_poly_taylor_shift_horner(ctx, f, c);
            d = sc_zz_poly_taylor_shift_divconquer(ctx, f, c);
            a = sc_zz_poly_taylor_shift(ctx, f, c);
            q = sc_zz_poly_compose(ctx, f, g);
            ok = poly_equal(h, d) && poly_equal(h, a) && poly_equal(h, q);
            sc_value_free_many(7, f, c, g, h, d, a, q);
            if (!ok)
                return 0;
        }
    return 1;
}

static int test_named(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { -7, 5, -2, 3 }, ec[] = { 19, 33, 16, 3 };
    sc_value *f = poly_si(ctx, r, 4, fc), *e = poly_si(ctx, r, 4, ec);
    sc_value *c = sc_zz_from_str(ctx, "2"), *h, *d, *a;
    int ok;

    if (f == NULL || e == NULL || c == NULL)
        return 0;
    a = sc_call2(ctx, "taylor_shift", f, c);
    h = sc_call2(ctx, "taylor_shift_horner", f, c);
    d = sc_call2(ctx, "taylor_shift_dc", f, c);
    ok = poly_equal(a, e) && poly_equal(h, e) && poly_equal(d, e);
    sc_value_free_many(6, f, e, c, h, d, a);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "ZZ[x]", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_fixed(&ctx, &r) || !test_random(&ctx, &r) || !test_named(&ctx, &r)) {
        fprintf(stderr, "Taylor-shift test failed: %s\n", ctx.error);
        sc_context_clear(&ctx);
        return 1;
    }
    sc_context_clear(&ctx);
    puts("Taylor-shift tests passed");
    return 0;
}
