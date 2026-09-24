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

static int same_poly(const sc_value *a, const sc_value *b)
{
    size_t i;

    if (a == NULL || b == NULL || a->data.zz_poly.length != b->data.zz_poly.length)
        return 0;
    for (i = 0; i < a->data.zz_poly.length; i++)
        if (mpz_cmp(a->data.zz_poly.coeff[i], b->data.zz_poly.coeff[i]) != 0)
            return 0;
    return 1;
}

static sc_value *random_poly(sc_context *ctx, sc_parent *r, size_t n)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);
    size_t i;

    if (f == NULL)
        return NULL;
    for (i = 0; i < n; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], (rand() % 13) - 6);
    if (n != 0 && mpz_sgn(f->data.zz_poly.coeff[n - 1]) == 0)
        mpz_set_si(f->data.zz_poly.coeff[n - 1], (rand() & 1) ? 1 : -1);
    return f;
}

static mpz_t *matrix_new(size_t n)
{
    mpz_t *a = malloc(n * n * sizeof(mpz_t));
    size_t i;

    if (a == NULL)
        return NULL;
    for (i = 0; i < n * n; i++)
        mpz_init(a[i]);
    return a;
}

static void matrix_free(mpz_t *a, size_t n)
{
    size_t i;

    for (i = 0; i < n * n; i++)
        mpz_clear(a[i]);
    free(a);
}

static void sylvester(mpz_t *a, size_t n, const sc_value *f, const sc_value *g)
{
    size_t m = f->data.zz_poly.length - 1;
    size_t d = g->data.zz_poly.length - 1;
    size_t i, j;

    for (i = 0; i < d; i++)
        for (j = 0; j <= m; j++)
            mpz_set(a[i * n + i + j], f->data.zz_poly.coeff[m - j]);
    for (i = 0; i < m; i++)
        for (j = 0; j <= d; j++)
            mpz_set(a[(d + i) * n + i + j], g->data.zz_poly.coeff[d - j]);
}

static void bareiss_det(mpz_t det, mpz_t *a, size_t n)
{
    size_t i, j, k, p;
    int sign = 1;
    mpz_t prev, t, u;

    if (n == 0) {
        mpz_set_ui(det, 1);
        return;
    }
    mpz_inits(prev, t, u, NULL);
    mpz_set_ui(prev, 1);
    for (k = 0; k + 1 < n; k++) {
        for (p = k; p < n && mpz_sgn(a[p * n + k]) == 0; p++)
            ;
        if (p == n) {
            mpz_set_ui(det, 0);
            goto done;
        }
        if (p != k) {
            for (j = k; j < n; j++)
                mpz_swap(a[k * n + j], a[p * n + j]);
            sign = -sign;
        }
        for (i = k + 1; i < n; i++)
            for (j = k + 1; j < n; j++) {
                mpz_mul(t, a[i * n + j], a[k * n + k]);
                mpz_mul(u, a[i * n + k], a[k * n + j]);
                mpz_sub(t, t, u);
                if (k != 0)
                    mpz_divexact(t, t, prev);
                mpz_set(a[i * n + j], t);
            }
        mpz_set(prev, a[k * n + k]);
    }
    mpz_set(det, a[(n - 1) * n + n - 1]);
    if (sign < 0)
        mpz_neg(det, det);
done:
    mpz_clears(prev, t, u, NULL);
}

static void resultant_reference(mpz_t det, const sc_value *f, const sc_value *g)
{
    size_t m, d, n;
    mpz_t *a;

    if (f->data.zz_poly.length == 0 || g->data.zz_poly.length == 0) {
        mpz_set_ui(det, 0);
        return;
    }
    m = f->data.zz_poly.length - 1;
    d = g->data.zz_poly.length - 1;
    n = m + d;
    a = matrix_new(n);
    if (a == NULL) {
        mpz_set_ui(det, 0);
        return;
    }
    sylvester(a, n, f, g);
    bareiss_det(det, a, n);
    matrix_free(a, n);
}

static int test_content_primitive(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { -6, 12, -18 }, pc[] = { -1, 2, -3 };
    sc_value *f = poly_si(ctx, r, 3, fc), *p = poly_si(ctx, r, 3, pc);
    sc_value *c = sc_zz_poly_content(ctx, f);
    sc_value *pp = sc_zz_poly_primitive_part(ctx, f);
    sc_value *two = sc_zz_from_str(ctx, "2");
    sc_value *half = sc_zz_poly_scalar_divexact(ctx, f, two);
    int ok = c != NULL && mpz_cmp_ui(c->data.z, 6) == 0 && same_poly(pp, p) &&
             half != NULL && mpz_cmp_si(half->data.zz_poly.coeff[0], -3) == 0;

    sc_value_free_many(6, f, p, c, pp, two, half);
    return ok;
}

static int test_gcd(sc_context *ctx, sc_parent *r)
{
    const long hc[] = { 1, -2, 0, 1 }, ac[] = { 1, 1 }, bc[] = { 2, 1 };
    sc_value *h = poly_si(ctx, r, 4, hc), *aa = poly_si(ctx, r, 2, ac);
    sc_value *bb = poly_si(ctx, r, 2, bc), *six = sc_zz_from_str(ctx, "6");
    sc_value *ten = sc_zz_from_str(ctx, "10"), *two = sc_zz_from_str(ctx, "2");
    sc_value *fa = sc_zz_poly_mul(ctx, h, aa), *fb = sc_zz_poly_mul(ctx, h, bb);
    sc_value *f = sc_zz_poly_scalar_mul(ctx, fa, six);
    sc_value *g = sc_zz_poly_scalar_mul(ctx, fb, ten);
    sc_value *expect = sc_zz_poly_scalar_mul(ctx, h, two);
    sc_value *d1 = sc_zz_poly_gcd(ctx, f, g), *d2 = sc_zz_poly_gcd_pseudo(ctx, f, g);
    int ok = same_poly(d1, expect) && same_poly(d2, expect);

    sc_value_free_many(14, h, aa, bb, six, ten, two, fa, fb, f, g, expect, d1, d2, NULL);
    return ok;
}

static int test_gcd_random(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(41891);
    for (trial = 0; trial < 50; trial++) {
        sc_value *f = random_poly(ctx, r, 2 + rand() % 5);
        sc_value *g = random_poly(ctx, r, 2 + rand() % 5);
        sc_value *d1 = sc_zz_poly_gcd(ctx, f, g);
        sc_value *d2 = sc_zz_poly_gcd_pseudo(ctx, f, g);
        sc_value *q1 = d1 ? sc_zz_poly_divexact(ctx, f, d1) : NULL;
        sc_value *q2 = d1 ? sc_zz_poly_divexact(ctx, g, d1) : NULL;
        int ok = same_poly(d1, d2) && q1 != NULL && q2 != NULL;

        sc_value_free_many(6, f, g, d1, d2, q1, q2);
        if (!ok)
            return 0;
    }
    return 1;
}


static int check_xgcd_identity(sc_context *ctx, const sc_value *f, const sc_value *g,
                               const sc_value *xg)
{
    sc_value *h, *cof, *u, *v, *uf, *vg, *sum;
    int ok;

    if (xg == NULL || xg->kind != SC_VALUE_PAIR)
        return 0;
    h = xg->data.pair.first;
    cof = xg->data.pair.second;
    if (cof == NULL || cof->kind != SC_VALUE_PAIR)
        return 0;
    u = cof->data.pair.first;
    v = cof->data.pair.second;
    uf = sc_zz_poly_mul(ctx, u, f);
    vg = sc_zz_poly_mul(ctx, v, g);
    sum = uf && vg ? sc_zz_poly_add(ctx, uf, vg) : NULL;
    ok = same_poly(sum, h);
    sc_value_free_many(3, uf, vg, sum);
    return ok;
}

static int test_fast_pseudodiv(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(27491);
    for (trial = 0; trial < 80; trial++) {
        size_t an = 3 + rand() % 18, bn = 2 + rand() % (an - 1);
        sc_value *a = random_poly(ctx, r, an), *b = random_poly(ctx, r, bn);
        sc_value *x = sc_zz_poly_pseudodiv(ctx, a, b);
        sc_value *y = sc_zz_poly_pseudodiv_fast(ctx, a, b);
        int ok = x != NULL && y != NULL &&
                 same_poly(x->data.pair.first, y->data.pair.first) &&
                 same_poly(x->data.pair.second, y->data.pair.second);

        sc_value_free_many(4, a, b, x, y);
        if (!ok)
            return 0;
    }
    {
        size_t bn = 257, an = bn + 18;
        sc_value *a = random_poly(ctx, r, an), *b = random_poly(ctx, r, bn);
        sc_value *x = sc_zz_poly_pseudodiv(ctx, a, b);
        sc_value *y = sc_zz_poly_pseudodiv_fast(ctx, a, b);
        int ok = x != NULL && y != NULL &&
                 same_poly(x->data.pair.first, y->data.pair.first) &&
                 same_poly(x->data.pair.second, y->data.pair.second);

        sc_value_free_many(4, a, b, x, y);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_hgcd_random(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(51203);
    for (trial = 0; trial < 12; trial++) {
        size_t hn = 2 + rand() % 3, an = 14 + rand() % 9, bn = 12 + rand() % 9;
        sc_value *h = random_poly(ctx, r, hn), *a = random_poly(ctx, r, an);
        sc_value *b = random_poly(ctx, r, bn), *f = sc_zz_poly_mul(ctx, h, a);
        sc_value *g = sc_zz_poly_mul(ctx, h, b);
        sc_value *hg = sc_zz_poly_gcd_hgcd(ctx, f, g);
        sc_value *ref = sc_zz_poly_gcd(ctx, f, g);
        int ok = same_poly(hg, ref);

        sc_value_free_many(7, h, a, b, f, g, hg, ref);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_hgcd_transform(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(51211);
    for (trial = 0; trial < 12; trial++) {
        size_t an = 14 + rand() % 9, bn = an - 1 - rand() % 4, mid = an / 2;
        sc_value *a = random_poly(ctx, r, an), *b = random_poly(ctx, r, bn);
        sc_value *c = NULL, *d = NULL, *p, *q, *det, *g0, *g1, *pp0, *pp1;
        sc_zz_poly_mat2 m = { 0 };
        int ok;

        if (!sc_zz_poly_hgcd_pseudo(ctx, &m, a, b) ||
            !sc_zz_poly_mat2_apply(ctx, &c, &d, &m, a, b))
            return 0;
        p = sc_zz_poly_mul(ctx, m.a00, m.a11);
        q = sc_zz_poly_mul(ctx, m.a01, m.a10);
        det = p && q ? sc_zz_poly_sub(ctx, p, q) : NULL;
        g0 = det ? sc_zz_poly_gcd(ctx, a, b) : NULL;
        g1 = g0 ? sc_zz_poly_gcd(ctx, c, d) : NULL;
        pp0 = g1 ? sc_zz_poly_primitive_part(ctx, g0) : NULL;
        pp1 = pp0 ? sc_zz_poly_primitive_part(ctx, g1) : NULL;
        ok = det != NULL && det->data.zz_poly.length == 1 &&
             mpz_sgn(det->data.zz_poly.coeff[0]) != 0 &&
             c->data.zz_poly.length > mid &&
             (d->data.zz_poly.length == 0 || d->data.zz_poly.length <= mid) &&
             same_poly(pp0, pp1);
        sc_value_free_many(11, a, b, c, d, p, q, det, g0, g1, pp0, pp1);
        sc_zz_poly_mat2_clear(&m);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_xgcd(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { 2, 1, -1 }, gc[] = { 1, -1, 1 };
    const long cc[] = { 2 }, xc[] = { 0, 1 };
    sc_value *f = poly_si(ctx, r, 3, fc), *g = poly_si(ctx, r, 3, gc);
    sc_value *xg = sc_zz_poly_xgcd(ctx, f, g), *res = sc_zz_poly_resultant(ctx, f, g);
    sc_value *h = xg ? xg->data.pair.first : NULL;
    sc_value *c = poly_si(ctx, r, 1, cc), *x = poly_si(ctx, r, 2, xc);
    sc_value *ng = sc_zz_poly_gcd(ctx, c, x), *nx = sc_zz_poly_xgcd(ctx, c, x);
    sc_value *nh = nx ? nx->data.pair.first : NULL;
    int ok = check_xgcd_identity(ctx, f, g, xg) && h != NULL && SC_VALUE_POLY == h->kind &&
             h->data.zz_poly.length == 1 && mpz_cmp_si(h->data.zz_poly.coeff[0], -3) == 0 &&
             res != NULL && mpz_cmp_ui(res->data.z, 9) == 0 &&
             mpz_divisible_p(res->data.z, h->data.zz_poly.coeff[0]) &&
             ng != NULL && ng->data.zz_poly.length == 1 &&
             mpz_cmp_ui(ng->data.zz_poly.coeff[0], 1) == 0 &&
             check_xgcd_identity(ctx, c, x, nx) && nh != NULL &&
             nh->data.zz_poly.length == 1 &&
             mpz_cmpabs_ui(nh->data.zz_poly.coeff[0], 2) == 0;

    sc_value_free_many(8, f, g, xg, res, c, x, ng, nx);
    return ok;
}

static int test_xgcd_random(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(93731);
    for (trial = 0; trial < 100; trial++) {
        sc_value *f = random_poly(ctx, r, 1 + rand() % 6);
        sc_value *g = random_poly(ctx, r, 1 + rand() % 6);
        sc_value *xg = sc_zz_poly_xgcd(ctx, f, g);
        sc_value *res = sc_zz_poly_resultant(ctx, f, g);
        sc_value *h = xg ? xg->data.pair.first : NULL;
        int ok = check_xgcd_identity(ctx, f, g, xg) && res != NULL;

        if (ok && f->data.zz_poly.length > 1 && g->data.zz_poly.length > 1 &&
            h->data.zz_poly.length == 1 && mpz_sgn(h->data.zz_poly.coeff[0]) != 0)
            ok = mpz_divisible_p(res->data.z, h->data.zz_poly.coeff[0]);
        sc_value_free_many(4, f, g, xg, res);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_resultant_fixed(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { 1, 2 }, gc[] = { 4, 3 };
    sc_value *f = poly_si(ctx, r, 2, fc), *g = poly_si(ctx, r, 2, gc);
    sc_value *res = sc_zz_poly_resultant(ctx, f, g);
    int ok = res != NULL && mpz_cmp_si(res->data.z, 5) == 0;

    sc_value_free_many(3, f, g, res);
    return ok;
}

static int test_resultant_terminal_scalar(sc_context *ctx, sc_parent *r)
{
    const long ac[] = { -2, -1, -5 }, bc[] = { -3, -1, -5 };
    sc_value *a = poly_si(ctx, r, 3, ac), *b = poly_si(ctx, r, 3, bc);
    sc_value *prs = sc_zz_poly_subres_prs_last(ctx, a, b);
    sc_value *res = sc_zz_poly_resultant(ctx, a, b);
    sc_value *last = prs ? prs->data.pair.first : NULL;
    sc_value *h = prs ? prs->data.pair.second : NULL;
    int ok = last != NULL && last->data.zz_poly.length == 1 &&
             mpz_cmp_ui(last->data.zz_poly.coeff[0], 5) == 0 && h != NULL &&
             mpz_cmp_ui(h->data.z, 25) == 0 && res != NULL &&
             mpz_cmp_ui(res->data.z, 25) == 0;

    sc_value_free_many(4, a, b, prs, res);
    return ok;
}

static int test_resultant_random(sc_context *ctx, sc_parent *r)
{
    size_t trial;
    mpz_t ref;

    srand(81673);
    mpz_init(ref);
    for (trial = 0; trial < 100; trial++) {
        size_t m = 1 + rand() % 4, n = 1 + rand() % 4;
        sc_value *f = random_poly(ctx, r, m + 1), *g = random_poly(ctx, r, n + 1);
        sc_value *res = sc_zz_poly_resultant(ctx, f, g);
        sc_value *bare = sc_zz_poly_resultant_bareiss_impl(ctx, f, g);
        sc_value *sub = sc_zz_poly_resultant_subresultant_impl(ctx, f, g);
        int ok;

        resultant_reference(ref, f, g);
        ok = res != NULL && bare != NULL && sub != NULL &&
             mpz_cmp(res->data.z, ref) == 0 && mpz_cmp(bare->data.z, ref) == 0 &&
             mpz_cmp(sub->data.z, ref) == 0;
        sc_value_free_many(5, f, g, res, bare, sub);
        if (!ok) {
            mpz_clear(ref);
            return 0;
        }
    }
    mpz_clear(ref);
    return 1;
}

static int test_resultant_common_factor(sc_context *ctx, sc_parent *r)
{
    const long hc[] = { -1, 1 }, ac[] = { 2, 1 }, bc[] = { 3, 1 };
    sc_value *h = poly_si(ctx, r, 2, hc), *a = poly_si(ctx, r, 2, ac);
    sc_value *b = poly_si(ctx, r, 2, bc), *f = sc_zz_poly_mul(ctx, h, a);
    sc_value *g = sc_zz_poly_mul(ctx, h, b), *res = sc_zz_poly_resultant(ctx, f, g);
    int ok = res != NULL && mpz_sgn(res->data.z) == 0;

    sc_value_free_many(6, h, a, b, f, g, res);
    return ok;
}

static int test_edge_cases(sc_context *ctx, sc_parent *r)
{
    const long fc[] = { -2 }, gc[] = { 1, 0, 3 }, nc[] = { 4, -2 };
    sc_value *f = poly_si(ctx, r, 1, fc), *g = poly_si(ctx, r, 3, gc);
    sc_value *n = poly_si(ctx, r, 2, nc), *z = poly_si(ctx, r, 0, NULL);
    sc_value *res = sc_zz_poly_resultant(ctx, f, g);
    sc_value *rz = sc_zz_poly_resultant(ctx, z, g);
    sc_value *d = sc_zz_poly_gcd(ctx, z, n);
    int ok = res != NULL && mpz_cmp_ui(res->data.z, 4) == 0 &&
             rz != NULL && mpz_sgn(rz->data.z) == 0 && d != NULL &&
             d->data.zz_poly.length == 2 && mpz_cmp_si(d->data.zz_poly.coeff[0], -4) == 0 &&
             mpz_cmp_si(d->data.zz_poly.coeff[1], 2) == 0;

    sc_value_free_many(7, f, g, n, z, res, rz, d);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "ZZ[x]", SC_PARENT_POLY, &SC_ZZ, "x" };
    int ok;

    sc_context_init(&ctx);
    ok = test_content_primitive(&ctx, &r) && test_gcd(&ctx, &r) &&
         test_gcd_random(&ctx, &r) && test_fast_pseudodiv(&ctx, &r) &&
         test_hgcd_random(&ctx, &r) && test_hgcd_transform(&ctx, &r) &&
         test_xgcd(&ctx, &r) &&
         test_xgcd_random(&ctx, &r) && test_resultant_fixed(&ctx, &r) &&
         test_resultant_terminal_scalar(&ctx, &r) &&
         test_resultant_random(&ctx, &r) &&
         test_resultant_common_factor(&ctx, &r) && test_edge_cases(&ctx, &r);
    sc_context_clear(&ctx);
    if (!ok) {
        fprintf(stderr, "polynomial gcd/resultant test failed\n");
        return 1;
    }
    puts("polynomial gcd/resultant tests passed");
    return 0;
}
