#include "smallcas.h"

#include <stdio.h>
#include <stdlib.h>

static sc_value *random_poly(sc_context *ctx, sc_parent *r, size_t n)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);
    size_t i;

    if (f == NULL)
        return NULL;
    for (i = 0; i < n; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], (rand() % 2001) - 1000);
    sc_zz_poly_normalize(f);
    return f;
}

static int test_fixed(sc_context *ctx, sc_parent *r)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, 4);
    sc_value *x = sc_zz_from_str(ctx, "12345678901234567890");
    sc_value *h, *d, *a;
    int ok;

    if (f == NULL || x == NULL)
        return 0;
    mpz_set_si(f->data.zz_poly.coeff[0], -7);
    mpz_set_si(f->data.zz_poly.coeff[1], 5);
    mpz_set_si(f->data.zz_poly.coeff[2], -3);
    mpz_set_si(f->data.zz_poly.coeff[3], 2);
    h = sc_zz_poly_evaluate_horner(ctx, f, x);
    d = sc_zz_poly_evaluate_divconquer(ctx, f, x);
    a = sc_zz_poly_evaluate(ctx, f, x);
    ok = h != NULL && d != NULL && a != NULL && mpz_cmp(h->data.z, d->data.z) == 0 &&
         mpz_cmp(h->data.z, a->data.z) == 0;
    sc_value_free_many(5, f, x, h, d, a);
    return ok;
}

static int test_random(sc_context *ctx, sc_parent *r)
{
    size_t n, trial;

    srand(92031);
    for (n = 0; n <= 96; n++)
        for (trial = 0; trial < 10; trial++) {
            sc_value *f = random_poly(ctx, r, n);
            sc_value *x = sc_value_new_zz_checked(ctx), *h, *d, *a;
            int ok;

            if (f == NULL || x == NULL)
                return 0;
            mpz_set_si(x->data.z, (rand() % 2000001) - 1000000);
            h = sc_zz_poly_evaluate_horner(ctx, f, x);
            d = sc_zz_poly_evaluate_divconquer(ctx, f, x);
            a = sc_zz_poly_evaluate(ctx, f, x);
            ok = h != NULL && d != NULL && a != NULL &&
                 mpz_cmp(h->data.z, d->data.z) == 0 && mpz_cmp(h->data.z, a->data.z) == 0;
            sc_value_free_many(5, f, x, h, d, a);
            if (!ok)
                return 0;
        }
    return 1;
}

static int test_callable(sc_context *ctx, sc_parent *r)
{
    static const long coeff[] = { 1, 2, 3 };
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, 3);
    sc_value *x = sc_zz_from_str(ctx, "10"), *y, *e, *h, *d;
    size_t i;
    int ok;

    if (f == NULL || x == NULL)
        return 0;
    for (i = 0; i < 3; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], coeff[i]);
    if (!sc_env_set(ctx, "f", f)) {
        sc_value_free_many(2, f, x);
        return 0;
    }
    y = sc_call1(ctx, "f", x);
    e = sc_call2(ctx, "evaluate", f, x);
    h = sc_call2(ctx, "evaluate_horner", f, x);
    d = sc_call2(ctx, "evaluate_dc", f, x);
    ok = y != NULL && e != NULL && h != NULL && d != NULL &&
         mpz_cmp_ui(y->data.z, 321) == 0 && mpz_cmp(y->data.z, e->data.z) == 0 &&
         mpz_cmp(y->data.z, h->data.z) == 0 && mpz_cmp(y->data.z, d->data.z) == 0;
    sc_value_free_many(6, f, x, y, e, h, d);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "ZZ[x]", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_fixed(&ctx, &r) || !test_random(&ctx, &r) || !test_callable(&ctx, &r)) {
        fprintf(stderr, "evaluation test failed: %s\n", ctx.error);
        sc_context_clear(&ctx);
        return 1;
    }
    sc_context_clear(&ctx);
    puts("evaluation tests passed");
    return 0;
}
