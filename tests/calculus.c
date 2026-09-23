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

static int test_derivatives(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { -7, 5, -2, 0, 3 };
    const long d1c[] = { 5, -4, 0, 12 };
    const long d2c[] = { -4, 0, 36 };
    const long d4c[] = { 72 };
    sc_value *f = poly_si(ctx, r, 5, fc), *e1 = poly_si(ctx, r, 4, d1c);
    sc_value *e2 = poly_si(ctx, r, 3, d2c), *e4 = poly_si(ctx, r, 1, d4c);
    sc_value *d1, *d2, *d4, *d5, *d0, *two;
    int ok;

    if (f == NULL || e1 == NULL || e2 == NULL || e4 == NULL)
        return 0;
    d1 = sc_zz_poly_derivative(ctx, f);
    d2 = sc_zz_poly_nth_derivative(ctx, f, 2);
    d4 = sc_zz_poly_nth_derivative(ctx, f, 4);
    d5 = sc_zz_poly_nth_derivative(ctx, f, 5);
    d0 = sc_zz_poly_nth_derivative(ctx, f, 0);
    two = sc_zz_from_str(ctx, "2");
    ok = poly_equal(d1, e1) && poly_equal(d2, e2) && poly_equal(d4, e4) &&
         d5 != NULL && d5->data.zz_poly.length == 0 && poly_equal(d0, f);
    if (ok && two != NULL) {
        sc_value *named1 = sc_call1(ctx, "derivative", f);
        sc_value *named2 = sc_call2(ctx, "derivative", f, two);
        sc_value *named3 = sc_call2(ctx, "nth_derivative", f, two);

        ok = poly_equal(named1, e1) && poly_equal(named2, e2) && poly_equal(named3, e2);
        sc_value_free_many(3, named1, named2, named3);
    }
    sc_value_free_many(10, f, e1, e2, e4, d1, d2, d4, d5, d0, two);
    return ok;
}


static int test_nth_derivative_direct(sc_context *ctx, sc_parent *r)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, 25);
    size_t i, k;

    if (f == NULL)
        return 0;
    for (i = 0; i < 25; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], (long)((i * i + 3 * i + 7) % 19) - 9);
    if (mpz_sgn(f->data.zz_poly.coeff[24]) == 0)
        mpz_set_ui(f->data.zz_poly.coeff[24], 1);
    for (k = 0; k <= 25; k++) {
        sc_value *direct = sc_zz_poly_nth_derivative(ctx, f, k);
        sc_value *repeat = sc_value_copy_checked(ctx, f);
        size_t j;

        for (j = 0; repeat != NULL && j < k; j++) {
            sc_value *next = sc_zz_poly_derivative(ctx, repeat);
            sc_value_free(repeat);
            repeat = next;
        }
        if (!poly_equal(direct, repeat)) {
            sc_value_free_many(3, f, direct, repeat);
            return 0;
        }
        sc_value_free_many(2, direct, repeat);
    }
    sc_value_free(f);
    return 1;
}

static int test_discriminant(sc_context *ctx, sc_parent *r)
{
    const long qc[] = { 1, 2, 3 }, cc[] = { 1, 2, 0, 3 }, lc[] = { 5, 7 };
    const long kc[] = { 19 };
    sc_value *q = poly_si(ctx, r, 3, qc), *c = poly_si(ctx, r, 4, cc);
    sc_value *l = poly_si(ctx, r, 2, lc), *k = poly_si(ctx, r, 1, kc);
    sc_value *z = sc_value_new_zz_poly_checked(ctx, r, 0);
    sc_value *dq, *dc, *dl, *dk, *dz;
    int ok;

    if (q == NULL || c == NULL || l == NULL || k == NULL || z == NULL)
        return 0;
    dq = sc_zz_poly_discriminant(ctx, q);
    dc = sc_zz_poly_discriminant(ctx, c);
    dl = sc_zz_poly_discriminant(ctx, l);
    dk = sc_zz_poly_discriminant(ctx, k);
    dz = sc_zz_poly_discriminant(ctx, z);
    ok = dq != NULL && dc != NULL && dl != NULL && dk != NULL && dz != NULL &&
         mpz_cmp_si(dq->data.z, -8) == 0 && mpz_cmp_si(dc->data.z, -339) == 0 &&
         mpz_cmp_ui(dl->data.z, 1) == 0 && mpz_sgn(dk->data.z) == 0 &&
         mpz_sgn(dz->data.z) == 0;
    sc_value_free_many(10, q, c, l, k, z, dq, dc, dl, dk, dz);
    return ok;
}

static int test_squarefree(sc_context *ctx, sc_parent *r)
{
    const long amc[] = { -1, 1 }, bpc[] = { 2, 1 }, epc[] = { -2, 1, 1 };
    const long onepc[] = { 1, 1 }, kc[] = { -12 }, onec[] = { 1 };
    sc_value *am = poly_si(ctx, r, 2, amc), *bp = poly_si(ctx, r, 2, bpc);
    sc_value *expect = poly_si(ctx, r, 3, epc), *onep = poly_si(ctx, r, 2, onepc);
    sc_value *sq, *rep, *six, *scaled, *distinct, *twelve, *scaled_distinct;
    sc_value *sf, *ir, *is, *iz, *z, *k, *one, *sfk, *ik;
    int ok;

    if (am == NULL || bp == NULL || expect == NULL || onep == NULL)
        return 0;
    sq = sc_zz_poly_mul(ctx, am, am);
    rep = sq ? sc_zz_poly_mul(ctx, sq, bp) : NULL;
    six = sc_zz_from_str(ctx, "6");
    scaled = rep && six ? sc_zz_poly_scalar_mul(ctx, rep, six) : NULL;
    distinct = sc_zz_poly_mul(ctx, onep, bp);
    twelve = sc_zz_from_str(ctx, "12");
    scaled_distinct = distinct && twelve ? sc_zz_poly_scalar_mul(ctx, distinct, twelve) : NULL;
    sf = scaled ? sc_zz_poly_squarefree_part(ctx, scaled) : NULL;
    ir = scaled ? sc_zz_poly_is_squarefree(ctx, scaled) : NULL;
    is = scaled_distinct ? sc_zz_poly_is_squarefree(ctx, scaled_distinct) : NULL;
    z = sc_value_new_zz_poly_checked(ctx, r, 0);
    iz = z ? sc_zz_poly_is_squarefree(ctx, z) : NULL;
    k = poly_si(ctx, r, 1, kc);
    one = poly_si(ctx, r, 1, onec);
    sfk = k ? sc_zz_poly_squarefree_part(ctx, k) : NULL;
    ik = k ? sc_zz_poly_is_squarefree(ctx, k) : NULL;
    ok = poly_equal(sf, expect) && ir != NULL && mpz_sgn(ir->data.z) == 0 &&
         is != NULL && mpz_cmp_ui(is->data.z, 1) == 0 && iz != NULL &&
         mpz_sgn(iz->data.z) == 0 && poly_equal(sfk, one) && ik != NULL &&
         mpz_cmp_ui(ik->data.z, 1) == 0;
    sc_value_free_many(19, am, bp, expect, onep, sq, rep, six, scaled, distinct, twelve,
                       scaled_distinct, sf, ir, is, z, iz, k, one, sfk);
    sc_value_free(ik);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "ZZ[x]", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_derivatives(&ctx, &r) || !test_nth_derivative_direct(&ctx, &r) ||
        !test_discriminant(&ctx, &r) || !test_squarefree(&ctx, &r)) {
        fprintf(stderr, "calculus test failed: %s\n", ctx.error);
        sc_context_clear(&ctx);
        return 1;
    }
    sc_context_clear(&ctx);
    puts("calculus tests passed");
    return 0;
}
