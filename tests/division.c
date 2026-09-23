#include "smallcas.h"

#include <stdio.h>
#include <stdlib.h>

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

static int test_nonmonic_exact(sc_context *ctx, sc_parent *r)
{
    const long bc[] = { 2, 2 }, qc[] = { 1, 2 };
    sc_value *b = poly_si(ctx, r, 2, bc), *q = poly_si(ctx, r, 2, qc);
    sc_value *a = sc_zz_poly_mul_classical(ctx, b, q);
    sc_value *qr = sc_zz_poly_divrem(ctx, a, b);
    sc_value *qo = sc_zz_poly_quo(ctx, a, b);
    sc_value *e = sc_zz_poly_divexact(ctx, a, b);
    int ok = qr != NULL && same_poly(qr->data.pair.first, q) &&
             qr->data.pair.second->data.zz_poly.length == 0 &&
             same_poly(qo, q) && same_poly(e, q);

    sc_value_free_many(6, a, b, q, qr, qo, e);
    return ok;
}

static int test_remainder(sc_context *ctx, sc_parent *r)
{
    const long ac[] = { 1, 0, 1 }, bc[] = { 1, 1 };
    const long qc[] = { -1, 1 }, rc[] = { 2 };
    sc_value *a = poly_si(ctx, r, 3, ac), *b = poly_si(ctx, r, 2, bc);
    sc_value *q = poly_si(ctx, r, 2, qc), *rem = poly_si(ctx, r, 1, rc);
    sc_value *qr = sc_zz_poly_divrem(ctx, a, b);
    sc_value *qo = sc_zz_poly_quo(ctx, a, b);
    sc_value *e = sc_zz_poly_divexact(ctx, a, b);
    int ok = qr != NULL && same_poly(qr->data.pair.first, q) &&
             same_poly(qr->data.pair.second, rem) && same_poly(qo, q) && e == NULL;

    ctx->error[0] = '\0';
    sc_value_free_many(7, a, b, q, rem, qr, qo, e);
    return ok;
}

static int test_nondivisible_lead(sc_context *ctx, sc_parent *r)
{
    const long ac[] = { 0, 1 }, bc[] = { 1, 2 };
    sc_value *a = poly_si(ctx, r, 2, ac), *b = poly_si(ctx, r, 2, bc);
    sc_value *qr = sc_zz_poly_divrem(ctx, a, b);
    int ok = qr == NULL;

    ctx->error[0] = '\0';
    sc_value_free_many(3, a, b, qr);
    return ok;
}

static int test_pseudodiv(sc_context *ctx, sc_parent *r)
{
    const long ac[] = { 1, 0, 1 }, bc[] = { 1, 2 };
    const long qc[] = { -1, 2 }, rc[] = { 5 };
    sc_value *a = poly_si(ctx, r, 3, ac), *b = poly_si(ctx, r, 2, bc);
    sc_value *q = poly_si(ctx, r, 2, qc), *rem = poly_si(ctx, r, 1, rc);
    sc_value *qr = sc_zz_poly_pseudodiv(ctx, a, b);
    sc_value *copy = sc_value_copy_checked(ctx, qr);
    int ok = qr != NULL && copy != NULL && same_poly(qr->data.pair.first, q) &&
             same_poly(qr->data.pair.second, rem) &&
             same_poly(copy->data.pair.first, q) && same_poly(copy->data.pair.second, rem);

    sc_value_free_many(6, a, b, q, rem, qr, copy);
    return ok;
}

static int test_short_pseudodiv(sc_context *ctx, sc_parent *r)
{
    const long ac[] = { 3, 4 }, bc[] = { 1, 0, 2 };
    sc_value *a = poly_si(ctx, r, 2, ac), *b = poly_si(ctx, r, 3, bc);
    sc_value *qr = sc_zz_poly_pseudodiv(ctx, a, b);
    int ok = qr != NULL && qr->data.pair.first->data.zz_poly.length == 0 &&
             same_poly(qr->data.pair.second, a);

    sc_value_free_many(3, a, b, qr);
    return ok;
}


static sc_value *random_poly(sc_context *ctx, sc_parent *r, size_t n)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);
    size_t i;

    if (f == NULL)
        return NULL;
    for (i = 0; i < n; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], (rand() % 19) - 9);
    if (n != 0 && mpz_sgn(f->data.zz_poly.coeff[n - 1]) == 0)
        mpz_set_si(f->data.zz_poly.coeff[n - 1], (rand() & 1) ? 2 : -3);
    return f;
}

static int test_random_divrem(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(9173);
    for (trial = 0; trial < 80; trial++) {
        size_t bn = 1 + rand() % 7, qn = rand() % 8, rn = rand() % bn;
        sc_value *b = random_poly(ctx, r, bn), *q = random_poly(ctx, r, qn);
        sc_value *rem = random_poly(ctx, r, rn);
        sc_value *prod = sc_zz_poly_mul_classical(ctx, b, q);
        sc_value *a = sc_zz_poly_add(ctx, prod, rem);
        sc_value *qr = sc_zz_poly_divrem(ctx, a, b);
        int ok = qr != NULL && same_poly(qr->data.pair.first, q) &&
                 same_poly(qr->data.pair.second, rem);

        sc_value_free_many(6, a, b, q, rem, prod, qr);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_random_pseudodiv(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(27181);
    for (trial = 0; trial < 80; trial++) {
        size_t an = rand() % 9, bn = 1 + rand() % 7, i;
        size_t d = an >= bn ? an - bn + 1 : 0;
        sc_value *a = random_poly(ctx, r, an), *b = random_poly(ctx, r, bn);
        sc_value *qr = sc_zz_poly_pseudodiv(ctx, a, b);
        sc_value *lhs = sc_value_copy_checked(ctx, a);
        sc_value *prod = sc_zz_poly_mul_classical(ctx, qr->data.pair.first, b);
        sc_value *rhs = sc_zz_poly_add(ctx, prod, qr->data.pair.second);
        mpz_t scale;
        int ok;

        mpz_init(scale);
        mpz_pow_ui(scale, b->data.zz_poly.coeff[bn - 1], (unsigned long)d);
        for (i = 0; i < lhs->data.zz_poly.length; i++)
            mpz_mul(lhs->data.zz_poly.coeff[i], lhs->data.zz_poly.coeff[i], scale);
        ok = same_poly(lhs, rhs) && qr->data.pair.second->data.zz_poly.length < bn;
        mpz_clear(scale);
        sc_value_free_many(6, a, b, qr, lhs, prod, rhs);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_dc_division(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(81173);
    for (trial = 0; trial < 40; trial++) {
        size_t bn = 12 + rand() % 25, qn = 32 + rand() % 50, rn = rand() % bn;
        sc_value *b = random_poly(ctx, r, bn), *q = random_poly(ctx, r, qn);
        sc_value *rem = random_poly(ctx, r, rn);
        sc_value *prod = sc_zz_poly_mul_classical(ctx, b, q);
        sc_value *a = sc_zz_poly_add(ctx, prod, rem);
        sc_value *qc = sc_zz_poly_quo_classical(ctx, a, b);
        sc_value *qd = sc_zz_poly_quo_dc(ctx, a, b);
        sc_value *qr = sc_zz_poly_divrem_dc(ctx, a, b);
        sc_value *auto_qr = sc_zz_poly_divrem(ctx, a, b);
        int ok = same_poly(qc, q) && same_poly(qd, q) && qr != NULL &&
                 same_poly(qr->data.pair.first, q) && same_poly(qr->data.pair.second, rem) &&
                 auto_qr != NULL && same_poly(auto_qr->data.pair.first, q) &&
                 same_poly(auto_qr->data.pair.second, rem);

        sc_value_free_many(10, a, b, q, rem, prod, qc, qd, qr, auto_qr, NULL);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_dc_short_quotient_long_divisor(sc_context *ctx, sc_parent *r)
{
    sc_value *b = random_poly(ctx, r, 257), *q = random_poly(ctx, r, 19);
    sc_value *rem = random_poly(ctx, r, 113);
    sc_value *prod = sc_zz_poly_mul_classical(ctx, b, q);
    sc_value *a = sc_zz_poly_add(ctx, prod, rem);
    sc_value *got = sc_zz_poly_quo_dc(ctx, a, b);
    int ok = same_poly(got, q);

    sc_value_free_many(6, a, b, q, rem, prod, got);
    return ok;
}

static int test_dc_exact(sc_context *ctx, sc_parent *r)
{
    sc_value *b = random_poly(ctx, r, 24), *q = random_poly(ctx, r, 55);
    sc_value *a = sc_zz_poly_mul_classical(ctx, b, q);
    sc_value *got = sc_zz_poly_divexact(ctx, a, b);
    int ok = same_poly(got, q);

    sc_value_free_many(4, a, b, q, got);
    return ok;
}


static int test_dc_nondivisible(sc_context *ctx, sc_parent *r)
{
    sc_value *a = random_poly(ctx, r, 70), *b = random_poly(ctx, r, 30);
    sc_value *q;
    int ok;

    mpz_set_ui(b->data.zz_poly.coeff[29], 2);
    mpz_set_ui(a->data.zz_poly.coeff[69], 1);
    q = sc_zz_poly_quo(ctx, a, b);
    ok = q == NULL;
    ctx->error[0] = '\0';
    sc_value_free_many(3, a, b, q);
    return ok;
}


static sc_value *shift_poly(sc_context *ctx, sc_parent *r, const sc_value *a, size_t k)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, a->data.zz_poly.length + k);
    size_t i;

    if (f == NULL)
        return NULL;
    for (i = 0; i < a->data.zz_poly.length; i++)
        mpz_set(f->data.zz_poly.coeff[k + i], a->data.zz_poly.coeff[i]);
    return f;
}

static int test_bidirectional_exact(sc_context *ctx, sc_parent *r)
{
    sc_value *b = random_poly(ctx, r, 70), *q = random_poly(ctx, r, 65);
    sc_value *a, *qb, *qe;
    int ok;

    if (mpz_sgn(b->data.zz_poly.coeff[0]) == 0)
        mpz_set_si(b->data.zz_poly.coeff[0], 3);
    a = sc_zz_poly_mul_classical(ctx, b, q);
    qb = sc_zz_poly_quo_bidirectional(ctx, a, b);
    qe = sc_zz_poly_divexact(ctx, a, b);
    ok = same_poly(qb, q) && same_poly(qe, q);
    sc_value_free_many(5, a, b, q, qb, qe);
    return ok;
}

static int test_bidirectional_shifted(sc_context *ctx, sc_parent *r)
{
    sc_value *c = random_poly(ctx, r, 55), *q = random_poly(ctx, r, 64);
    sc_value *b, *a, *qe;
    int ok;

    if (mpz_sgn(c->data.zz_poly.coeff[0]) == 0)
        mpz_set_si(c->data.zz_poly.coeff[0], -2);
    b = shift_poly(ctx, r, c, 3);
    a = sc_zz_poly_mul_classical(ctx, b, q);
    qe = sc_zz_poly_divexact(ctx, a, b);
    ok = same_poly(qe, q);
    sc_value_free_many(6, a, b, c, q, qe, NULL);
    return ok;
}


static int test_bidirectional_long_divisor(sc_context *ctx, sc_parent *r)
{
    sc_value *b = random_poly(ctx, r, 257), *q = random_poly(ctx, r, 65);
    sc_value *a, *got;
    int ok;

    if (mpz_sgn(b->data.zz_poly.coeff[0]) == 0)
        mpz_set_si(b->data.zz_poly.coeff[0], 7);
    a = sc_zz_poly_mul_classical(ctx, b, q);
    got = sc_zz_poly_quo_bidirectional(ctx, a, b);
    ok = same_poly(got, q);
    sc_value_free_many(4, a, b, q, got);
    return ok;
}

static int test_bidirectional_middle_reject(sc_context *ctx, sc_parent *r)
{
    size_t qn = 64, ln = (qn + 1) / 2;
    sc_value *b = random_poly(ctx, r, 80), *q = random_poly(ctx, r, qn);
    sc_value *a, *qb, *qe;
    int ok;

    if (mpz_sgn(b->data.zz_poly.coeff[0]) == 0)
        mpz_set_si(b->data.zz_poly.coeff[0], 5);
    a = sc_zz_poly_mul_classical(ctx, b, q);
    mpz_add_ui(a->data.zz_poly.coeff[ln + 5], a->data.zz_poly.coeff[ln + 5], 1);
    qb = sc_zz_poly_quo_bidirectional(ctx, a, b);
    qe = sc_zz_poly_divexact(ctx, a, b);
    ok = same_poly(qb, q) && qe == NULL;
    ctx->error[0] = '\0';
    sc_value_free_many(5, a, b, q, qb, qe);
    return ok;
}

static int test_mulders_case(sc_context *ctx, sc_parent *r, size_t bn, size_t qn)
{
    size_t rn = bn > 1 ? bn - 1 : 0;
    sc_value *b = random_poly(ctx, r, bn), *q = random_poly(ctx, r, qn);
    sc_value *rem = random_poly(ctx, r, rn);
    sc_value *prod = sc_zz_poly_mul_classical(ctx, b, q);
    sc_value *a = sc_zz_poly_add(ctx, prod, rem);
    sc_value *qm = sc_zz_poly_quo_mulders(ctx, a, b);
    sc_value *qrm = sc_zz_poly_divrem_mulders(ctx, a, b);
    sc_value *qa = sc_zz_poly_quo(ctx, a, b);
    sc_value *qra = sc_zz_poly_divrem(ctx, a, b);
    int ok = same_poly(qm, q) && qrm != NULL && same_poly(qrm->data.pair.first, q) &&
             same_poly(qrm->data.pair.second, rem) && same_poly(qa, q) && qra != NULL &&
             same_poly(qra->data.pair.first, q) && same_poly(qra->data.pair.second, rem);

    sc_value_free_many(9, a, b, q, rem, prod, qm, qrm, qa, qra);
    return ok;
}

static int test_mulders_division(sc_context *ctx, sc_parent *r)
{
    srand(93281);
    return test_mulders_case(ctx, r, 48, 48) &&
           test_mulders_case(ctx, r, 49, 49) &&
           test_mulders_case(ctx, r, 73, 73) &&
           test_mulders_case(ctx, r, 80, 55);
}

static int test_mulders_nondivisible(sc_context *ctx, sc_parent *r)
{
    sc_value *a = random_poly(ctx, r, 95), *b = random_poly(ctx, r, 48);
    sc_value *q;
    int ok;

    mpz_set_ui(b->data.zz_poly.coeff[47], 2);
    mpz_set_ui(a->data.zz_poly.coeff[94], 1);
    q = sc_zz_poly_quo_mulders(ctx, a, b);
    ok = q == NULL;
    ctx->error[0] = '\0';
    sc_value_free_many(3, a, b, q);
    return ok;
}


static int test_newton_inverse(sc_context *ctx, sc_parent *r)
{
    size_t n;

    srand(14033);
    for (n = 1; n <= 96; n++) {
        sc_value *f = random_poly(ctx, r, 1 + rand() % 35);
        sc_value *g, *p;
        int ok;

        mpz_set_si(f->data.zz_poly.coeff[0], (n & 1) ? 1 : -1);
        g = sc_zz_poly_inv_series(ctx, f, n);
        p = g == NULL ? NULL : sc_zz_poly_mullow(ctx, f, g, n);
        ok = p != NULL && p->data.zz_poly.length == 1 &&
             mpz_cmp_ui(p->data.zz_poly.coeff[0], 1) == 0;
        sc_value_free_many(3, f, g, p);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_newton_division(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(55021);
    for (trial = 0; trial < 30; trial++) {
        size_t bn = 2 + rand() % 53, qn = 64 + rand() % 55, rn = rand() % bn;
        sc_value *b = random_poly(ctx, r, bn), *q = random_poly(ctx, r, qn);
        sc_value *rem = random_poly(ctx, r, rn);
        sc_value *prod, *a, *qnw, *qrnw, *qa, *qra;
        int ok;

        mpz_set_ui(b->data.zz_poly.coeff[bn - 1], 1);
        prod = sc_zz_poly_mul_classical(ctx, b, q);
        a = sc_zz_poly_add(ctx, prod, rem);
        qnw = sc_zz_poly_quo_newton(ctx, a, b);
        qrnw = sc_zz_poly_divrem_newton(ctx, a, b);
        qa = sc_zz_poly_quo(ctx, a, b);
        qra = sc_zz_poly_divrem(ctx, a, b);
        ok = same_poly(qnw, q) && qrnw != NULL && same_poly(qrnw->data.pair.first, q) &&
             same_poly(qrnw->data.pair.second, rem) && same_poly(qa, q) && qra != NULL &&
             same_poly(qra->data.pair.first, q) && same_poly(qra->data.pair.second, rem);
        sc_value_free_many(10, b, q, rem, prod, a, qnw, qrnw, qa, qra, NULL);
        if (!ok)
            return 0;
    }
    {
        size_t bn = 257, qn = 65;
        sc_value *b = random_poly(ctx, r, bn), *q = random_poly(ctx, r, qn);
        sc_value *rem = random_poly(ctx, r, 23), *prod, *a, *got;
        int ok;

        mpz_set_ui(b->data.zz_poly.coeff[bn - 1], 1);
        prod = sc_zz_poly_mul_classical(ctx, b, q);
        a = sc_zz_poly_add(ctx, prod, rem);
        got = sc_zz_poly_quo_newton(ctx, a, b);
        ok = same_poly(got, q);
        sc_value_free_many(7, b, q, rem, prod, a, got, NULL);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_newton_preinverse(sc_context *ctx, sc_parent *r)
{
    sc_value *b = random_poly(ctx, r, 45), *pinv;
    size_t qn[2] = { 67, 93 }, i;

    mpz_set_ui(b->data.zz_poly.coeff[44], 1);
    pinv = sc_zz_poly_preinverse_newton(ctx, b, 96);
    if (pinv == NULL)
        return 0;
    for (i = 0; i < 2; i++) {
        sc_value *q = random_poly(ctx, r, qn[i]);
        sc_value *rem = random_poly(ctx, r, 17 + i);
        sc_value *prod = sc_zz_poly_mul_classical(ctx, b, q);
        sc_value *a = sc_zz_poly_add(ctx, prod, rem);
        sc_value *got = sc_zz_poly_quo_preinv(ctx, a, b, pinv);
        int ok = same_poly(got, q);

        sc_value_free_many(5, q, rem, prod, a, got);
        if (!ok) {
            sc_value_free_many(2, b, pinv);
            return 0;
        }
    }
    sc_value_free_many(2, b, pinv);
    return 1;
}

static int test_newton_exact(sc_context *ctx, sc_parent *r)
{
    sc_value *b = random_poly(ctx, r, 51), *q = random_poly(ctx, r, 88);
    sc_value *a, *got, *bad;
    int ok;

    mpz_set_ui(b->data.zz_poly.coeff[50], 1);
    a = sc_zz_poly_mul_classical(ctx, b, q);
    got = sc_zz_poly_divexact(ctx, a, b);
    mpz_add_ui(a->data.zz_poly.coeff[7], a->data.zz_poly.coeff[7], 1);
    bad = sc_zz_poly_divexact(ctx, a, b);
    ok = same_poly(got, q) && bad == NULL;
    ctx->error[0] = '\0';
    sc_value_free_many(5, b, q, a, got, bad);
    return ok;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "PolynomialRing(ZZ)", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_nonmonic_exact(&ctx, &r) || !test_remainder(&ctx, &r) ||
        !test_nondivisible_lead(&ctx, &r) || !test_pseudodiv(&ctx, &r) ||
        !test_short_pseudodiv(&ctx, &r) || !test_random_divrem(&ctx, &r) ||
        !test_random_pseudodiv(&ctx, &r) || !test_dc_division(&ctx, &r) ||
        !test_dc_exact(&ctx, &r) || !test_dc_short_quotient_long_divisor(&ctx, &r) ||
        !test_dc_nondivisible(&ctx, &r) ||
        !test_bidirectional_exact(&ctx, &r) || !test_bidirectional_shifted(&ctx, &r) ||
        !test_bidirectional_long_divisor(&ctx, &r) ||
        !test_bidirectional_middle_reject(&ctx, &r) ||
        !test_mulders_division(&ctx, &r) || !test_mulders_nondivisible(&ctx, &r) ||
        !test_newton_inverse(&ctx, &r) || !test_newton_division(&ctx, &r) ||
        !test_newton_preinverse(&ctx, &r) || !test_newton_exact(&ctx, &r)) {
        fprintf(stderr, "polynomial division test failed: %s\n", ctx.error);
        return 1;
    }
    puts("polynomial division tests passed");
    return 0;
}
