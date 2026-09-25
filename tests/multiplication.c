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

static size_t max_bits(const sc_value *a)
{
    size_t i, bits, max = 0;

    for (i = 0; i < a->data.zz_poly.length; i++) {
        bits = mpz_sgn(a->data.zz_poly.coeff[i]) == 0 ?
               0 : mpz_sizeinbase(a->data.zz_poly.coeff[i], 2);
        if (bits > max)
            max = bits;
    }
    return max;
}

static mp_bitcnt_t ks_bits(const sc_value *a, const sc_value *b)
{
    size_t n = a->data.zz_poly.length < b->data.zz_poly.length ?
               a->data.zz_poly.length : b->data.zz_poly.length;
    size_t log = 0;

    if (n != 0)
        n--;
    while (n != 0) {
        log++;
        n >>= 1;
    }
    return (mp_bitcnt_t)(max_bits(a) + max_bits(b) + log + 1);
}

static sc_value *make_poly(sc_context *ctx, sc_parent *r, size_t n, unsigned seed)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);
    size_t i;

    if (f == NULL)
        return NULL;
    srand(seed);
    for (i = 0; i < n; i++)
        mpz_set_si(f->data.zz_poly.coeff[i], (rand() % 2001) - 1000);
    if (mpz_sgn(f->data.zz_poly.coeff[n - 1]) == 0)
        mpz_set_ui(f->data.zz_poly.coeff[n - 1], 1);
    return f;
}

static int compare_karatsuba(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *c = sc_zz_poly_mul_classical(ctx, a, b);
    sc_value *k = sc_zz_poly_mul_karatsuba(ctx, a, b);
    int ok = same_poly(c, k);

    sc_value_free_many(2, c, k);
    return ok;
}

static int compare_ks(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *c = sc_zz_poly_mul_classical(ctx, a, b);
    sc_value *k = sc_zz_poly_mul_ks(ctx, a, b, ks_bits(a, b));
    int ok = same_poly(c, k);

    sc_value_free_many(2, c, k);
    return ok;
}

static sc_value *mul_toom3(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t large = a->data.zz_poly.length > b->data.zz_poly.length ?
                   a->data.zz_poly.length : b->data.zz_poly.length;
    sc_zz_poly_toom3_ws ws;

    if (!sc_zz_poly_toom3_ws_init(ctx, &ws, a, b, (large + 2) / 3))
        return NULL;
    return sc_zz_poly_mul_toom3(ctx, &ws);
}

static int compare_toom3(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *c = sc_zz_poly_mul_classical(ctx, a, b);
    sc_value *t = mul_toom3(ctx, a, b);
    int ok = same_poly(c, t);

    sc_value_free_many(2, c, t);
    return ok;
}


static sc_value *make_ssa_poly(sc_context *ctx, sc_parent *r, size_t n,
                               size_t bits, unsigned seed)
{
    sc_value *f = sc_value_new_zz_poly_checked(ctx, r, n);

    if (f == NULL)
        return NULL;
    for (size_t i = 0; i < n; i++) {
        unsigned long v = (unsigned long)(seed + 17 * i + i * i + 1);

        mpz_set_ui(f->data.zz_poly.coeff[i], v);
        if (bits != 0)
            mpz_setbit(f->data.zz_poly.coeff[i], bits - 1);
        if ((i + seed) & 1)
            mpz_neg(f->data.zz_poly.coeff[i], f->data.zz_poly.coeff[i]);
    }
    return f;
}

static int test_ssa_case(sc_context *ctx, sc_parent *r, size_t an, size_t bn,
                         size_t abits, size_t bbits, unsigned seed)
{
    sc_value *a = make_ssa_poly(ctx, r, an, abits, seed);
    sc_value *b = make_ssa_poly(ctx, r, bn, bbits, seed + 1000);
    sc_value *want = a && b ? sc_zz_poly_mul_classical(ctx, a, b) : NULL;
    sc_value *got = a && b ? sc_zz_poly_mul_ssa(ctx, a, b) : NULL;
    int ok = same_poly(want, got);

    sc_value_free_many(4, a, b, want, got);
    return ok;
}

static int test_ntt_case(sc_context *ctx, sc_parent *r, size_t an, size_t bn,
                         size_t abits, size_t bbits, unsigned seed)
{
    sc_value *a = make_ssa_poly(ctx, r, an, abits, seed);
    sc_value *b = make_ssa_poly(ctx, r, bn, bbits, seed + 1000);
    sc_value *want = a && b ? sc_zz_poly_mul_classical(ctx, a, b) : NULL;
    sc_value *got = a && b ? sc_zz_poly_mul_ntt(ctx, a, b) : NULL;
    int ok = same_poly(want, got);

    sc_value_free_many(4, a, b, want, got);
    return ok;
}

static int test_ntt(sc_context *ctx, sc_parent *r)
{
    static const size_t cases[][4] = {
        { 1, 1, 20, 19 }, { 7, 5, 90, 93 }, { 19, 23, 260, 257 },
        { 7, 5, 700, 711 }, { 33, 32, 80, 83 },
        { 64, 64, 80, 83 }, { 70, 65, 20, 23 }
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
        if (!test_ntt_case(ctx, r, cases[i][0], cases[i][1], cases[i][2],
                           cases[i][3], (unsigned)(23000 + i)))
            return 0;
    return 1;
}

static int test_ssa(sc_context *ctx, sc_parent *r)
{
    static const size_t cases[][4] = {
        { 1, 1, 20, 19 }, { 7, 5, 90, 93 }, { 19, 23, 260, 257 },
        { 40, 37, 700, 711 }, { 33, 32, 80, 83 },
        { 64, 64, 80, 83 }, { 70, 65, 20, 23 }
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
        if (!test_ssa_case(ctx, r, cases[i][0], cases[i][1], cases[i][2],
                           cases[i][3], (unsigned)(21000 + i)))
            return 0;
    {
        sc_value *z = sc_value_new_zz_poly_checked(ctx, r, 0);
        sc_value *a = make_ssa_poly(ctx, r, 9, 80, 22000);
        sc_value *p = z && a ? sc_zz_poly_mul_ssa(ctx, z, a) : NULL;
        int ok = p != NULL && p->data.zz_poly.length == 0;

        sc_value_free_many(3, z, a, p);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_power(sc_context *ctx, sc_parent *r)
{
    sc_value *f = make_poly(ctx, r, 7, 27182);
    unsigned long e;

    if (f == NULL)
        return 0;
    for (e = 0; e <= 12; e++) {
        sc_value *got = sc_zz_poly_pow_binary(ctx, f, e);
        sc_value *want = sc_value_new_zz_poly_checked(ctx, r, 1);
        unsigned long i;
        int ok;

        if (want == NULL) {
            sc_value_free_many(2, f, got);
            return 0;
        }
        mpz_set_ui(want->data.zz_poly.coeff[0], 1);
        for (i = 0; i < e; i++) {
            sc_value *t = sc_zz_poly_mul_classical(ctx, want, f);

            sc_value_free(want);
            want = t;
            if (want == NULL)
                break;
        }
        ok = same_poly(got, want);
        sc_value_free_many(2, got, want);
        if (!ok) {
            sc_value_free(f);
            return 0;
        }
    }
    sc_value_free(f);
    return 1;
}

static int test_grid(sc_context *ctx, sc_parent *r)
{
    size_t an, bn;

    for (an = 12; an <= 61; an += 7) {
        for (bn = 12; bn <= 61; bn += 5) {
            sc_value *a = make_poly(ctx, r, an, (unsigned)(1000 + an));
            sc_value *b = make_poly(ctx, r, bn, (unsigned)(2000 + bn));
            int ok = compare_karatsuba(ctx, a, b) && compare_ks(ctx, a, b) &&
                     compare_toom3(ctx, a, b);
            sc_value_free_many(2, a, b);
            if (!ok)
                return 0;
        }
    }
    return 1;
}

static int test_ks_leading_zeros(sc_context *ctx, sc_parent *r)
{
    sc_value *a = sc_value_new_zz_poly_checked(ctx, r, 19);
    sc_value *b = sc_value_new_zz_poly_checked(ctx, r, 17);
    sc_value *want, *got;
    int ok;

    if (a == NULL || b == NULL)
        return 0;
    for (size_t i = 0; i < 15; i++) {
        mpz_set_si(a->data.zz_poly.coeff[i], (long)(7 * i + 3));
        if (i & 1)
            mpz_neg(a->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i]);
    }
    for (size_t i = 0; i < 13; i++) {
        mpz_set_si(b->data.zz_poly.coeff[i], (long)(5 * i + 1));
        if ((i & 3) == 0)
            mpz_neg(b->data.zz_poly.coeff[i], b->data.zz_poly.coeff[i]);
    }
    mpz_neg(a->data.zz_poly.coeff[14], a->data.zz_poly.coeff[14]);
    want = sc_zz_poly_mul_classical(ctx, a, b);
    got = sc_zz_poly_mul_ks(ctx, a, b, ks_bits(a, b));
    ok = same_poly(want, got);
    sc_value_free_many(4, a, b, want, got);
    return ok;
}

static int test_dispatch(sc_context *ctx, sc_parent *r)
{
    sc_value *a = make_poly(ctx, r, 40, 3001);
    sc_value *b = make_poly(ctx, r, 43, 3002);
    sc_value *c = sc_zz_poly_mul_classical(ctx, a, b);
    sc_value *d = sc_zz_poly_mul(ctx, a, b);
    int ok = same_poly(c, d);

    sc_value_free_many(4, a, b, c, d);
    return ok;
}

static int test_toom3_dispatch(sc_context *ctx, sc_parent *r)
{
    sc_value *a = sc_value_new_zz_poly_checked(ctx, r, 54);
    sc_value *b = sc_value_new_zz_poly_checked(ctx, r, 57);
    sc_value *c, *d;
    size_t i;
    int ok;

    if (a == NULL || b == NULL)
        return 0;
    for (i = 0; i < 54; i++) {
        mpz_set_ui(a->data.zz_poly.coeff[i], 2 * i + 1);
        mpz_setbit(a->data.zz_poly.coeff[i], 80);
        if (i & 1)
            mpz_neg(a->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i]);
    }
    for (i = 0; i < 57; i++) {
        mpz_set_ui(b->data.zz_poly.coeff[i], 3 * i + 1);
        mpz_setbit(b->data.zz_poly.coeff[i], 80);
        if (i % 3 == 0)
            mpz_neg(b->data.zz_poly.coeff[i], b->data.zz_poly.coeff[i]);
    }
    c = sc_zz_poly_mul_classical(ctx, a, b);
    d = sc_zz_poly_mul(ctx, a, b);
    ok = same_poly(c, d);
    sc_value_free_many(4, a, b, c, d);
    return ok;
}

static int test_recursive_toom3(sc_context *ctx, sc_parent *r)
{
    sc_value *a = sc_value_new_zz_poly_checked(ctx, r, 150);
    sc_value *b = sc_value_new_zz_poly_checked(ctx, r, 153);
    sc_value *c, *d;
    size_t i;
    int ok;

    if (a == NULL || b == NULL)
        return 0;
    for (i = 0; i < 150; i++) {
        mpz_set_ui(a->data.zz_poly.coeff[i], i + 1);
        mpz_setbit(a->data.zz_poly.coeff[i], 220);
        if (i & 1)
            mpz_neg(a->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i]);
    }
    for (i = 0; i < 153; i++) {
        mpz_set_ui(b->data.zz_poly.coeff[i], 5 * i + 1);
        mpz_setbit(b->data.zz_poly.coeff[i], 220);
        if (i & 2)
            mpz_neg(b->data.zz_poly.coeff[i], b->data.zz_poly.coeff[i]);
    }
    c = sc_zz_poly_mul_classical(ctx, a, b);
    d = sc_zz_poly_mul(ctx, a, b);
    ok = same_poly(c, d);
    sc_value_free_many(4, a, b, c, d);
    return ok;
}

static int test_karatsuba_dispatch(sc_context *ctx, sc_parent *r)
{
    sc_value *a = sc_value_new_zz_poly_checked(ctx, r, 17);
    sc_value *b = sc_value_new_zz_poly_checked(ctx, r, 17);
    sc_value *c, *d;
    int ok;

    if (a == NULL || b == NULL)
        return 0;
    mpz_set_ui(a->data.zz_poly.coeff[0], 1);
    mpz_setbit(a->data.zz_poly.coeff[16], 200);
    mpz_set_si(b->data.zz_poly.coeff[0], -1);
    mpz_setbit(b->data.zz_poly.coeff[16], 200);
    c = sc_zz_poly_mul_classical(ctx, a, b);
    d = sc_zz_poly_mul(ctx, a, b);
    ok = same_poly(c, d);
    sc_value_free_many(4, a, b, c, d);
    return ok;
}

static int test_cancellation(sc_context *ctx, sc_parent *r)
{
    sc_value *a = sc_value_new_zz_poly_checked(ctx, r, 24);
    sc_value *b = sc_value_new_zz_poly_checked(ctx, r, 24);
    size_t i;
    int ok;

    if (a == NULL || b == NULL)
        return 0;
    for (i = 0; i < 12; i++) {
        mpz_set_si(a->data.zz_poly.coeff[i], (long)i + 1);
        mpz_neg(a->data.zz_poly.coeff[i + 12], a->data.zz_poly.coeff[i]);
        mpz_set_si(b->data.zz_poly.coeff[i], (long)(3 * i + 1));
        mpz_neg(b->data.zz_poly.coeff[i + 12], b->data.zz_poly.coeff[i]);
    }
    ok = compare_karatsuba(ctx, a, b) && compare_ks(ctx, a, b) &&
         compare_toom3(ctx, a, b);
    sc_value_free_many(2, a, b);
    return ok;
}

static int test_balanced_mulmid(sc_context *ctx, sc_parent *r)
{
    size_t n;

    for (n = 17; n <= 73; n++) {
        size_t an = 2 * n - 1 - n % 4, bn = n - n % 3;
        sc_value *a = make_poly(ctx, r, an, (unsigned)(8100 + n));
        sc_value *b = make_poly(ctx, r, bn, (unsigned)(9100 + n));
        sc_value *want = sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
        sc_value *got = sc_zz_poly_mulmid_balanced(ctx, a, b, n);
        int ok = same_poly(want, got);

        sc_value_free_many(4, a, b, want, got);
        if (!ok)
            return 0;
    }
    {
        const long ac[] = { 3 }, bc[] = { -5 };
        sc_value *a = sc_value_new_zz_poly_checked(ctx, r, 1);
        sc_value *b = sc_value_new_zz_poly_checked(ctx, r, 1);
        sc_value *got;
        int ok;

        mpz_set_si(a->data.zz_poly.coeff[0], ac[0]);
        mpz_set_si(b->data.zz_poly.coeff[0], bc[0]);
        got = sc_zz_poly_mulmid_balanced(ctx, a, b, 33);
        ok = got != NULL && got->data.zz_poly.length == 0;
        sc_value_free_many(3, a, b, got);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_fft_mulmid(sc_context *ctx, sc_parent *r)
{
    static const size_t cases[][3] = {
        { 7, 20, 19 }, { 8, 35, 31 }, { 9, 90, 87 }, { 31, 48, 45 },
        { 32, 130, 127 }, { 33, 260, 251 }, { 65, 75, 73 }
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        size_t n = cases[i][0];
        size_t an = 2 * n - 1 - (i % 3);
        size_t bn = n - (i % 2);
        sc_value *a = make_ssa_poly(ctx, r, an, cases[i][1],
                                    (unsigned)(18100 + i));
        sc_value *b = make_ssa_poly(ctx, r, bn, cases[i][2],
                                    (unsigned)(19100 + i));
        sc_value *want = a && b ? sc_zz_poly_mulmid_classical(ctx, a, b,
                                                               n - 1, n) : NULL;
        sc_value *ntt = a && b ? sc_zz_poly_mulmid_ntt(ctx, a, b, n) : NULL;
        sc_value *ssa = a && b ? sc_zz_poly_mulmid_ssa(ctx, a, b, n) : NULL;
        int ok = same_poly(want, ntt) && same_poly(want, ssa);

        sc_value_free_many(5, a, b, want, ntt, ssa);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_toom63_mulmid(sc_context *ctx, sc_parent *r)
{
    size_t n;

    for (n = 18; n <= 90; n += 3) {
        size_t an = 2 * n - 1 - n % 5, bn = n - n % 4;
        sc_value *a = make_poly(ctx, r, an, (unsigned)(10100 + n));
        sc_value *b = make_poly(ctx, r, bn, (unsigned)(11100 + n));
        sc_value *want = sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
        sc_zz_poly_toom63_ws ws;
        sc_value *got = NULL;
        int ok;

        if (sc_zz_poly_toom63_ws_init(ctx, &ws, a, b, n / 3))
            got = sc_zz_poly_mulmid_toom63(ctx, &ws);
        ok = same_poly(want, got);
        sc_value_free_many(4, a, b, want, got);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_toom63_tail(sc_context *ctx, sc_parent *r)
{
    size_t n;

    for (n = 49; n <= 95; n++) {
        size_t an = 2 * n - 1 - n % 7, bn = n - n % 5;
        sc_value *a, *b, *want, *got;
        int ok;

        if (n % 3 == 0)
            continue;
        a = make_poly(ctx, r, an, (unsigned)(14100 + n));
        b = make_poly(ctx, r, bn, (unsigned)(15100 + n));
        want = sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
        got = sc_zz_poly_mulmid_toom63_tail(ctx, a, b, n);
        ok = same_poly(want, got);
        sc_value_free_many(4, a, b, want, got);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_recursive_toom63_tail(sc_context *ctx, sc_parent *r)
{
    size_t n;

    for (n = 145; n <= 146; n++) {
        sc_value *a = make_poly(ctx, r, 2 * n - 1, (unsigned)(16000 + n));
        sc_value *b = make_poly(ctx, r, n, (unsigned)(17000 + n));
        sc_value *want = sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
        sc_value *got = sc_zz_poly_mulmid_balanced(ctx, a, b, n);
        int ok = same_poly(want, got);

        sc_value_free_many(4, a, b, want, got);
        if (!ok)
            return 0;
    }
    return 1;
}

static int test_recursive_toom63(sc_context *ctx, sc_parent *r)
{
    size_t n = 144;
    sc_value *a = make_poly(ctx, r, 2 * n - 1, 12144);
    sc_value *b = make_poly(ctx, r, n, 13144);
    sc_value *want = sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
    sc_value *got = sc_zz_poly_mulmid_balanced(ctx, a, b, n);
    int ok = same_poly(want, got);

    sc_value_free_many(4, a, b, want, got);
    return ok;
}


static int test_fft_short_dispatch(sc_context *ctx, sc_parent *r)
{
    static const size_t edge_cases[][2] = { { 180, 180 }, { 2400, 48 } };

    for (size_t k = 0; k < 2; k++) {
        size_t n = edge_cases[k][0], bits = edge_cases[k][1];
        sc_value *a = make_ssa_poly(ctx, r, n, bits, (unsigned)(20100 + k));
        sc_value *b = make_ssa_poly(ctx, r, n, bits, (unsigned)(20200 + k));
        sc_value *full = k == 0 ? sc_zz_poly_mul_ssa(ctx, a, b) :
                                  sc_zz_poly_mul_ntt(ctx, a, b);
        sc_value *lo = sc_zz_poly_mullow(ctx, a, b, n);
        sc_value *hi = sc_zz_poly_mulhigh(ctx, a, b, n);
        sc_value lv, hv;
        int ok;

        if (full == NULL || lo == NULL || hi == NULL) {
            sc_value_free_many(5, a, b, full, lo, hi);
            return 0;
        }
        lv = sc_zz_poly_view(full, 0, n);
        hv = sc_zz_poly_view(full, full->data.zz_poly.length - n, n);
        ok = same_poly(lo, &lv) && same_poly(hi, &hv);
        sc_value_free_many(5, a, b, full, lo, hi);
        if (!ok)
            return 0;
    }
    {
        static const size_t mid_cases[][2] = { { 100, 100 }, { 800, 48 } };

        for (size_t k = 0; k < 2; k++) {
            size_t n = mid_cases[k][0], bits = mid_cases[k][1];
            sc_value *a = make_ssa_poly(ctx, r, 2 * n - 1, bits,
                                        (unsigned)(20300 + k));
            sc_value *b = make_ssa_poly(ctx, r, n, bits, (unsigned)(20400 + k));
            sc_value *want = k == 0 ? sc_zz_poly_mulmid_ssa(ctx, a, b, n) :
                                      sc_zz_poly_mulmid_ntt(ctx, a, b, n);
            sc_value *got = sc_zz_poly_mulmid_balanced(ctx, a, b, n);
            int ok = same_poly(want, got);

            sc_value_free_many(4, a, b, want, got);
            if (!ok)
                return 0;
        }
    }
    return 1;
}

static int test_mulhigh_short_tail(sc_context *ctx, sc_parent *r)
{
    size_t edge = 7, an = 4096, bn = 3072, total = an + bn - 1;
    sc_value *a = make_poly(ctx, r, an, 17001);
    sc_value *b = make_poly(ctx, r, bn, 17002);
    sc_value *want = sc_zz_poly_mulmid_classical(ctx, a, b, total - edge, edge);
    sc_value *got = sc_zz_poly_mulhigh(ctx, a, b, edge);
    int ok = same_poly(want, got);

    sc_value_free_many(4, a, b, want, got);
    return ok;
}

static int test_mulmid(sc_context *ctx, sc_parent *r)
{
    size_t trial;

    srand(44221);
    for (trial = 0; trial < 100; trial++) {
        size_t an = 1 + rand() % 70, bn = 1 + rand() % 70;
        size_t total = an + bn - 1, start = rand() % total;
        size_t n = 1 + rand() % (total + 5);
        sc_value *a = make_poly(ctx, r, an, (unsigned)(5000 + trial));
        sc_value *b = make_poly(ctx, r, bn, (unsigned)(7000 + trial));
        sc_value *full = sc_zz_poly_mul_classical(ctx, a, b);
        sc_value *mid = sc_zz_poly_mulmid(ctx, a, b, start, n);
        size_t want = n < total - start ? n : total - start;
        size_t edge = 1 + rand() % total;
        sc_value expect = sc_zz_poly_view(full, start, want);
        sc_value low_expect = sc_zz_poly_view(full, 0, edge);
        sc_value high_expect = sc_zz_poly_view(full, total - edge, edge);
        sc_value *low = sc_zz_poly_mullow(ctx, a, b, edge);
        sc_value *high = sc_zz_poly_mulhigh(ctx, a, b, edge);
        int ok = same_poly(mid, &expect) && same_poly(low, &low_expect) &&
                 same_poly(high, &high_expect);

        sc_value_free_many(6, a, b, full, mid, low, high);
        if (!ok)
            return 0;
    }
    return 1;
}

int main(void)
{
    sc_context ctx;
    sc_parent r = { "PolynomialRing(ZZ)", SC_PARENT_POLY, &SC_ZZ, "x" };

    sc_context_init(&ctx);
    if (!test_ssa(&ctx, &r) || !test_ntt(&ctx, &r) ||
        !test_power(&ctx, &r) || !test_grid(&ctx, &r) ||
        !test_ks_leading_zeros(&ctx, &r) || !test_dispatch(&ctx, &r) ||
        !test_toom3_dispatch(&ctx, &r) || !test_recursive_toom3(&ctx, &r) ||
        !test_karatsuba_dispatch(&ctx, &r) || !test_cancellation(&ctx, &r) ||
        !test_balanced_mulmid(&ctx, &r) || !test_fft_mulmid(&ctx, &r) ||
        !test_toom63_mulmid(&ctx, &r) ||
        !test_toom63_tail(&ctx, &r) || !test_recursive_toom63(&ctx, &r) ||
        !test_recursive_toom63_tail(&ctx, &r) ||
        !test_fft_short_dispatch(&ctx, &r) ||
        !test_mulhigh_short_tail(&ctx, &r) ||
        !test_mulmid(&ctx, &r)) {
        fprintf(stderr, "polynomial multiplication comparison failed: %s\n", ctx.error);
        return 1;
    }
    puts("polynomial multiplication comparison passed");
    return 0;
}
