#include "smallcas.h"

#include <stdio.h>
#include <stdlib.h>

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
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, 4);
    sc_value *g = sc_value_new_zz_poly_checked(ctx, r, 3);
    sc_value *h, *d, *a, *z, *lhs, *gz, *rhs;
    int ok;

    if (f == NULL || g == NULL)
        return 0;
    mpz_set_si(f->data.zz_poly.coeff[0], 1);
    mpz_set_si(f->data.zz_poly.coeff[1], 2);
    mpz_set_si(f->data.zz_poly.coeff[3], 1);
    mpz_set_si(g->data.zz_poly.coeff[0], 3);
    mpz_set_si(g->data.zz_poly.coeff[1], -1);
    mpz_set_si(g->data.zz_poly.coeff[2], 1);
    h = sc_zz_poly_compose_horner(ctx, f, g);
    d = sc_zz_poly_compose_divconquer(ctx, f, g);
    a = sc_zz_poly_compose(ctx, f, g);
    z = sc_zz_from_str(ctx, "7");
    lhs = z == NULL ? NULL : sc_zz_poly_evaluate(ctx, a, z);
    gz = z == NULL ? NULL : sc_zz_poly_evaluate(ctx, g, z);
    rhs = gz == NULL ? NULL : sc_zz_poly_evaluate(ctx, f, gz);
    ok = poly_equal(h, d) && poly_equal(h, a) && lhs != NULL && rhs != NULL &&
         mpz_cmp(lhs->data.z, rhs->data.z) == 0;
    sc_value_free_many(9, f, g, h, d, a, z, lhs, gz, rhs);
    return ok;
}

static int test_random(sc_context *ctx, sc_parent *r)
{
    size_t n, m, trial;

    srand(92326);
    for (n = 0; n <= 18; n++)
        for (m = 0; m <= 8; m++)
            for (trial = 0; trial < 4; trial++) {
                sc_value *f = random_poly(ctx, r, n, 8);
                sc_value *g = random_poly(ctx, r, m, 5);
                sc_value *z = sc_value_new_zz_checked(ctx), *h, *d, *a, *lhs, *gz, *rhs;
                int ok;

                if (f == NULL || g == NULL || z == NULL)
                    return 0;
                mpz_set_si(z->data.z, (rand() % 15) - 7);
                h = sc_zz_poly_compose_horner(ctx, f, g);
                d = sc_zz_poly_compose_divconquer(ctx, f, g);
                a = sc_zz_poly_compose(ctx, f, g);
                lhs = a == NULL ? NULL : sc_zz_poly_evaluate(ctx, a, z);
                gz = sc_zz_poly_evaluate(ctx, g, z);
                rhs = gz == NULL ? NULL : sc_zz_poly_evaluate(ctx, f, gz);
                ok = poly_equal(h, d) && poly_equal(h, a) && lhs != NULL && rhs != NULL &&
                     mpz_cmp(lhs->data.z, rhs->data.z) == 0;
                sc_value_free_many(9, f, g, z, h, d, a, lhs, gz, rhs);
                if (!ok)
                    return 0;
            }
    return 1;
}

static int test_callable(sc_context *ctx, sc_parent *r)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, 3);
    sc_value *g = sc_value_new_zz_poly_checked(ctx, r, 2);
    sc_value *y, *c, *h, *d;
    int ok;

    if (f == NULL || g == NULL)
        return 0;
    mpz_set_si(f->data.zz_poly.coeff[0], 1);
    mpz_set_si(f->data.zz_poly.coeff[2], 1);
    mpz_set_si(g->data.zz_poly.coeff[0], 1);
    mpz_set_si(g->data.zz_poly.coeff[1], 1);
    if (!sc_env_set(ctx, "f", f)) {
        sc_value_free_many(2, f, g);
        return 0;
    }
    y = sc_call1(ctx, "f", g);
    c = sc_call2(ctx, "compose", f, g);
    h = sc_call2(ctx, "compose_horner", f, g);
    d = sc_call2(ctx, "compose_dc", f, g);
    ok = poly_equal(y, c) && poly_equal(y, h) && poly_equal(y, d);
    sc_value_free_many(6, f, g, y, c, h, d);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "ZZ[x]", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_fixed(&ctx, &r) || !test_random(&ctx, &r) || !test_callable(&ctx, &r)) {
        fprintf(stderr, "composition test failed: %s\n", ctx.error);
        sc_context_clear(&ctx);
        return 1;
    }
    sc_context_clear(&ctx);
    puts("composition tests passed");
    return 0;
}
