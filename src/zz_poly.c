#include "smallcas.h"

#include <stdlib.h>
#include <string.h>

#define SC_MPZ_ADDEQ(x, y) mpz_add((x), (x), (y))
#define SC_MPZ_SUBEQ(x, y) mpz_sub((x), (x), (y))
#define SC_MPZ_MULEQ(x, y) mpz_mul((x), (x), (y))
#define SC_MPZ_SUBMUL(x, y, z) mpz_submul((x), (y), (z))
#define SC_ZP(v, i) ((v)->data.zz_poly.coeff[i])
#define SC_ZN(v) ((v)->data.zz_poly.length)
#define SC_ZLC(v) SC_ZP((v), SC_ZN(v) - 1)
#define SC_ZP_COEFF(v, i, z) \
    ((i) < (v)->data.zz_poly.length ? (v)->data.zz_poly.coeff[i] : (z))
#define SC_T3_SIDE(w, i) ((i) / (w)->m)
#define SC_T3_POS(w, i) ((i) % (w)->m)
#define SC_T3_EVAL(w, i, p) \
    ((w)->point[SC_T3_SIDE(w, i)][p]->data.zz_poly.coeff[SC_T3_POS(w, i)])
#define SC_T3_BLOCK(w, i, b) \
    SC_ZP_COEFF(&(w)->block[SC_T3_SIDE(w, i)][b], SC_T3_POS(w, i), (w)->zero)
#define SC_T3_PROD(w, j, i) ((w)->prod[j]->data.zz_poly.coeff[i])
#define SC_T3_PROD_OR_ZERO(w, j, i) SC_ZP_COEFF((w)->prod[j], i, (w)->zero)
#define SC_T3_ADDTO_RESULT(w, j, i) \
    SC_MPZ_ADDEQ((w)->result->data.zz_poly.coeff[(i) + (j) * (w)->m], \
                 SC_T3_PROD(w, j, i))
#define SC_T63_F(w, p, i) SC_ZP((w)->fpoint[p], i)
#define SC_T63_G(w, p, i) SC_ZP((w)->gpoint[p], i)
#define SC_T63_FB(w, b, i) SC_ZP_COEFF(&(w)->fblock[b], i, (w)->zero)
#define SC_T63_GB(w, b, i) SC_ZP_COEFF(&(w)->gblock[b], i, (w)->zero)
#define SC_T63_P(w, p, i) SC_ZP_COEFF((w)->prod[p], i, (w)->zero)
#define SC_T63_OUT(w, b, i) SC_ZP((w)->result, (b) * (w)->k + (i))

sc_value *sc_zz_poly_generator(sc_context *ctx, sc_parent *parent, const char *name)
{
    sc_value *r;

    if (parent->symbol != NULL && strcmp(parent->symbol, name) != 0) {
        sc_set_error(ctx, "polynomial indeterminate is already '%s'", parent->symbol);
        return NULL;
    }
    if (parent->symbol == NULL) {
        parent->symbol = sc_string_dup(ctx, name);
        if (parent->symbol == NULL)
            return NULL;
    }
    r = sc_value_new_zz_poly_checked(ctx, parent, 2);
    if (r == NULL)
        return NULL;
    mpz_set_ui(r->data.zz_poly.coeff[1], 1);
    return r;
}

sc_value *sc_zz_poly_from_zz(sc_context *ctx, sc_parent *parent, const sc_value *a)
{
    size_t length = mpz_sgn(a->data.z) == 0 ? 0 : 1;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, parent, length);

    if (r == NULL)
        return NULL;
    if (length != 0)
        mpz_set(r->data.zz_poly.coeff[0], a->data.z);
    return r;
}

sc_value *sc_zz_poly_add_impl(sc_context *ctx,
                                  const sc_value *a, const sc_value *b)
{
    size_t i, n = a->data.zz_poly.length > b->data.zz_poly.length ?
                  a->data.zz_poly.length : b->data.zz_poly.length;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n);

    if (r == NULL)
        return NULL;
    for (i = 0; i < n; i++) {
        if (i < a->data.zz_poly.length)
            mpz_set(r->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i]);
        if (i < b->data.zz_poly.length)
            mpz_add(r->data.zz_poly.coeff[i], r->data.zz_poly.coeff[i],
                    b->data.zz_poly.coeff[i]);
    }
    while (n != 0 && mpz_sgn(r->data.zz_poly.coeff[n - 1]) == 0)
        n--;
    sc_zz_poly_truncate(r, n);
    return r;
}

sc_value *sc_zz_poly_sub_impl(sc_context *ctx,
                                  const sc_value *a, const sc_value *b)
{
    size_t i, n = a->data.zz_poly.length > b->data.zz_poly.length ?
                  a->data.zz_poly.length : b->data.zz_poly.length;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n);

    if (r == NULL)
        return NULL;
    for (i = 0; i < n; i++) {
        if (i < a->data.zz_poly.length)
            mpz_set(r->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i]);
        if (i < b->data.zz_poly.length)
            mpz_sub(r->data.zz_poly.coeff[i], r->data.zz_poly.coeff[i],
                    b->data.zz_poly.coeff[i]);
    }
    while (n != 0 && mpz_sgn(r->data.zz_poly.coeff[n - 1]) == 0)
        n--;
    sc_zz_poly_truncate(r, n);
    return r;
}

sc_value *sc_zz_poly_scalar_mul_impl(sc_context *ctx,
                                         const sc_value *a, const sc_value *b)
{
    sc_value *r;
    size_t i;

    if (mpz_sgn(b->data.z) == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, a->data.zz_poly.length);
    if (r == NULL)
        return NULL;
    for (i = 0; i < a->data.zz_poly.length; i++)
        mpz_mul(r->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i], b->data.z);
    return r;
}

sc_value *sc_zz_poly_neg_impl(sc_context *ctx, const sc_value *a)
{
    size_t i;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent,
                                               a->data.zz_poly.length);

    if (r == NULL)
        return NULL;
    for (i = 0; i < a->data.zz_poly.length; i++)
        mpz_neg(r->data.zz_poly.coeff[i], a->data.zz_poly.coeff[i]);
    return r;
}

sc_value *sc_zz_poly_mul_constant(sc_context *ctx,
                                   const sc_value *a, const sc_value *b)
{
    const sc_value *f = a->data.zz_poly.length == 1 ? b : a;
    const sc_value *c = a->data.zz_poly.length == 1 ? a : b;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent,
                                               f->data.zz_poly.length);
    size_t i;

    if (r == NULL)
        return NULL;
    for (i = 0; i < f->data.zz_poly.length; i++)
        mpz_mul(r->data.zz_poly.coeff[i], f->data.zz_poly.coeff[i],
                c->data.zz_poly.coeff[0]);
    return r;
}

sc_value *sc_zz_poly_mul_classical(sc_context *ctx,
                                    const sc_value *a, const sc_value *b)
{
    size_t i, j, n;
    sc_value *r;

    if (a->data.zz_poly.length == 0 || b->data.zz_poly.length == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    n = a->data.zz_poly.length + b->data.zz_poly.length - 1;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (r == NULL)
        return NULL;
    for (i = 0; i < a->data.zz_poly.length; i++)
        for (j = 0; j < b->data.zz_poly.length; j++)
            mpz_addmul(r->data.zz_poly.coeff[i + j], a->data.zz_poly.coeff[i],
                       b->data.zz_poly.coeff[j]);
    return r;
}

sc_value *sc_zz_poly_mullow_classical(sc_context *ctx, const sc_value *a,
                                      const sc_value *b, size_t n)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length, i, j;
    size_t total = an + bn - 1;
    sc_value *r;

    if (n > total)
        n = total;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (r == NULL)
        return NULL;
    for (i = 0; i < an && i < n; i++)
        for (j = 0; j < bn && i + j < n; j++)
            mpz_addmul(SC_ZP(r, i + j), SC_ZP(a, i), SC_ZP(b, j));
    sc_zz_poly_normalize(r);
    return r;
}

sc_value *sc_zz_poly_mullow_dc(sc_context *ctx, const sc_value *a,
                               const sc_value *b, size_t n)
{
    size_t k = (n + 1) / 2, h = n - k, i;
    sc_value a0 = sc_zz_poly_view(a, 0, k), a1 = sc_zz_poly_view(a, k, h);
    sc_value b0 = sc_zz_poly_view(b, 0, k), b1 = sc_zz_poly_view(b, k, h);
    sc_value *z0, *z1, *z2, *r;

    if (n <= SC_MULLOW_DC_CUTOFF)
        return sc_zz_poly_mullow_classical(ctx, a, b, n);
    z0 = sc_zz_poly_mul(ctx, &a0, &b0);
    z1 = sc_zz_poly_mullow(ctx, &a0, &b1, h);
    z2 = sc_zz_poly_mullow(ctx, &a1, &b0, h);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (z0 == NULL || z1 == NULL || z2 == NULL || r == NULL)
        return sc_value_free_many_null(4, z0, z1, z2, r);
    for (i = 0; i < z0->data.zz_poly.length && i < n; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(z0, i));
    for (i = 0; i < z1->data.zz_poly.length; i++)
        SC_MPZ_ADDEQ(SC_ZP(r, k + i), SC_ZP(z1, i));
    for (i = 0; i < z2->data.zz_poly.length; i++)
        SC_MPZ_ADDEQ(SC_ZP(r, k + i), SC_ZP(z2, i));
    sc_zz_poly_normalize(r);
    sc_value_free_many(3, z0, z1, z2);
    return r;
}

sc_value *sc_zz_poly_mulhigh_reverse(sc_context *ctx, const sc_value *a,
                                     const sc_value *b, size_t n)
{
    size_t an = SC_ZN(a), bn = SC_ZN(b), al = an < n ? an : n, bl = bn < n ? bn : n;
    sc_value av = sc_zz_poly_view(a, an - al, al), bv = sc_zz_poly_view(b, bn - bl, bl);
    sc_value *ar = sc_zz_poly_reverse_impl(ctx, &av, al);
    sc_value *br = sc_zz_poly_reverse_impl(ctx, &bv, bl);
    sc_value *lo = ar != NULL && br != NULL ? sc_zz_poly_mullow(ctx, ar, br, n) : NULL;
    sc_value *r = lo != NULL ? sc_zz_poly_reverse_impl(ctx, lo, n) : NULL;

    sc_value_free_many(3, ar, br, lo);
    return r;
}

sc_value *sc_zz_poly_mulmid_classical(sc_context *ctx, const sc_value *a,
                                      const sc_value *b, size_t start, size_t n)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t i, j, d, lo, hi, total = an + bn - 1;
    sc_value *r;

    if (start >= total)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (n > total - start)
        n = total - start;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (r == NULL)
        return NULL;
    for (i = 0; i < n; i++) {
        d = start + i;
        lo = d >= bn - 1 ? d - (bn - 1) : 0;
        hi = d < an ? d : an - 1;
        for (j = lo; j <= hi; j++)
            mpz_addmul(SC_ZP(r, i), SC_ZP(a, j), SC_ZP(b, d - j));
    }
    sc_zz_poly_normalize(r);
    return r;
}

sc_value *sc_zz_poly_mulmid_toom42(sc_context *ctx, const sc_value *a,
                                    const sc_value *b, size_t n)
{
    size_t k = n / 2, i;
    sc_value f0 = sc_zz_poly_view(a, 0, n - 1);
    sc_value f1 = sc_zz_poly_view(a, k, n - 1);
    sc_value f2 = sc_zz_poly_view(a, n, n - 1);
    sc_value g0 = sc_zz_poly_view(b, 0, k), g1 = sc_zz_poly_view(b, k, k);
    sc_value *s0, *d, *s2, *p0, *p1, *p2, *lo, *hi, *r;

    s0 = sc_zz_poly_add(ctx, &f0, &f1);
    d = sc_zz_poly_sub(ctx, &g0, &g1);
    s2 = sc_zz_poly_add(ctx, &f1, &f2);
    p0 = sc_zz_poly_mulmid_balanced(ctx, s0, &g1, k);
    p1 = sc_zz_poly_mulmid_balanced(ctx, &f1, d, k);
    p2 = sc_zz_poly_mulmid_balanced(ctx, s2, &g0, k);
    lo = sc_zz_poly_add(ctx, p0, p1);
    hi = sc_zz_poly_sub(ctx, p2, p1);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (lo == NULL || hi == NULL || r == NULL)
        return sc_value_free_many_null(9, s0, d, s2, p0, p1, p2, lo, hi, r);
    for (i = 0; i < lo->data.zz_poly.length; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(lo, i));
    for (i = 0; i < hi->data.zz_poly.length; i++)
        mpz_set(SC_ZP(r, k + i), SC_ZP(hi, i));
    sc_zz_poly_normalize(r);
    sc_value_free_many(8, s0, d, s2, p0, p1, p2, lo, hi);
    return r;
}

sc_value *sc_zz_poly_mulmid_toom42_odd(sc_context *ctx, const sc_value *a,
                                        const sc_value *b, size_t n)
{
    size_t m = n - 1, i, j, ai, an = a->data.zz_poly.length;
    size_t bn = b->data.zz_poly.length;
    sc_value av = sc_zz_poly_view(a, 1, 2 * m - 1);
    sc_value bv = sc_zz_poly_view(b, 0, m);
    sc_value *core = sc_zz_poly_mulmid_balanced(ctx, &av, &bv, m);
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n);

    if (core == NULL || r == NULL)
        return sc_value_free_many_null(2, core, r);
    for (i = 0; i < core->data.zz_poly.length; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(core, i));
    if (bn == n)
        for (i = 0; i < m && i < an; i++)
            mpz_addmul(SC_ZP(r, i), SC_ZP(a, i), SC_ZP(b, m));
    for (j = 0; j < bn; j++) {
        ai = 2 * n - 2 - j;
        if (ai < an)
            mpz_addmul(SC_ZP(r, m), SC_ZP(a, ai), SC_ZP(b, j));
    }
    sc_zz_poly_normalize(r);
    sc_value_free(core);
    return r;
}

sc_value *sc_zz_poly_mulmid_toom63(sc_context *ctx, sc_zz_poly_toom63_ws *w)
{
    size_t i, j;

    for (i = 0; i < 2 * w->k - 1; i++) {
        mpz_sub(SC_T63_F(w, 3, i), SC_T63_FB(w, 3, i), SC_T63_FB(w, 1, i));
        mpz_set(SC_T63_F(w, 1, i), SC_T63_F(w, 3, i));
        SC_MPZ_SUBEQ(SC_T63_F(w, 1, i), SC_T63_FB(w, 1, i));
        SC_MPZ_SUBEQ(SC_T63_F(w, 1, i), SC_T63_FB(w, 2, i));
        mpz_set(SC_T63_F(w, 2, i), SC_T63_F(w, 3, i));
        mpz_addmul_ui(SC_T63_F(w, 2, i), SC_T63_FB(w, 1, i), 3);
        mpz_submul_ui(SC_T63_F(w, 2, i), SC_T63_FB(w, 2, i), 3);
        mpz_set(SC_T63_F(w, 0, i), SC_T63_F(w, 3, i));
        mpz_addmul_ui(SC_T63_F(w, 0, i), SC_T63_FB(w, 0, i), 2);
        mpz_submul_ui(SC_T63_F(w, 0, i), SC_T63_FB(w, 2, i), 2);
        mpz_mul_si(SC_T63_F(w, 4, i), SC_T63_F(w, 3, i), -2);
        SC_MPZ_ADDEQ(SC_T63_F(w, 4, i), SC_T63_FB(w, 4, i));
        SC_MPZ_SUBEQ(SC_T63_F(w, 4, i), SC_T63_FB(w, 2, i));
    }

    for (i = 0; i < w->k; i++) {
        mpz_add(SC_T63_G(w, 1, i), SC_T63_GB(w, 0, i), SC_T63_GB(w, 2, i));
        SC_MPZ_ADDEQ(SC_T63_G(w, 1, i), SC_T63_GB(w, 1, i));
        mpz_add(SC_T63_G(w, 2, i), SC_T63_GB(w, 0, i), SC_T63_GB(w, 2, i));
        SC_MPZ_SUBEQ(SC_T63_G(w, 2, i), SC_T63_GB(w, 1, i));
        mpz_set(SC_T63_G(w, 3, i), SC_T63_GB(w, 2, i));
        mpz_addmul_ui(SC_T63_G(w, 3, i), SC_T63_GB(w, 1, i), 2);
        mpz_addmul_ui(SC_T63_G(w, 3, i), SC_T63_GB(w, 0, i), 4);
    }

    for (j = 0; j < 5; j++)
        w->prod[j] = sc_zz_poly_mulmid_balanced(ctx, w->fpoint[j], w->gpoint[j], w->k);
    if (w->prod[0] == NULL || w->prod[1] == NULL || w->prod[2] == NULL ||
        w->prod[3] == NULL || w->prod[4] == NULL)
        return sc_zz_poly_toom63_ws_abort(w);

    for (i = 0; i < w->k; i++) {
        mpz_sub(SC_T63_OUT(w, 2, i), SC_T63_P(w, 3, i), SC_T63_P(w, 2, i));
        mpz_divexact_ui(SC_T63_OUT(w, 2, i), SC_T63_OUT(w, 2, i), 3);
        mpz_sub(SC_T63_OUT(w, 1, i), SC_T63_P(w, 2, i), SC_T63_P(w, 1, i));
        mpz_fdiv_q_2exp(SC_T63_OUT(w, 1, i), SC_T63_OUT(w, 1, i), 1);
        SC_MPZ_ADDEQ(SC_T63_OUT(w, 1, i), SC_T63_OUT(w, 2, i));
        mpz_sub(SC_T63_OUT(w, 0, i), SC_T63_P(w, 0, i), SC_T63_P(w, 1, i));
        SC_MPZ_ADDEQ(SC_T63_OUT(w, 0, i), SC_T63_OUT(w, 2, i));
        mpz_fdiv_q_2exp(SC_T63_OUT(w, 0, i), SC_T63_OUT(w, 0, i), 1);
        SC_MPZ_ADDEQ(SC_T63_OUT(w, 2, i), SC_T63_OUT(w, 1, i));
        SC_MPZ_ADDEQ(SC_T63_OUT(w, 2, i), SC_T63_P(w, 4, i));
    }

    return sc_zz_poly_toom63_ws_finish(w);
}

sc_value *sc_zz_poly_mulmid_toom63_tail(sc_context *ctx, const sc_value *a,
                                                const sc_value *b, size_t n)
{
    size_t r = n % 3, m = n - r, i, j, ai, d;
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    sc_value av = sc_zz_poly_view(a, r, 2 * m - 1);
    sc_value bv = sc_zz_poly_view(b, 0, m);
    sc_value *core = sc_zz_poly_mulmid_balanced(ctx, &av, &bv, m);
    sc_value *out = sc_value_new_zz_poly_checked(ctx, a->parent, n);

    if (core == NULL || out == NULL)
        return sc_value_free_many_null(2, core, out);
    for (i = 0; i < core->data.zz_poly.length; i++)
        mpz_set(SC_ZP(out, i), SC_ZP(core, i));
    for (j = m; j < bn; j++)
        for (i = 0; i < m; i++) {
            ai = n - 1 + i - j;
            if (ai < an)
                mpz_addmul(SC_ZP(out, i), SC_ZP(a, ai), SC_ZP(b, j));
        }
    for (i = m; i < n; i++) {
        d = n - 1 + i;
        for (j = 0; j < bn; j++) {
            ai = d - j;
            if (ai < an)
                mpz_addmul(SC_ZP(out, i), SC_ZP(a, ai), SC_ZP(b, j));
        }
    }
    sc_zz_poly_normalize(out);
    sc_value_free(core);
    return out;
}

sc_value *sc_zz_poly_mul_ks(sc_context *ctx, const sc_value *a, const sc_value *b,
                                mp_bitcnt_t bits)
{
    size_t i, n = a->data.zz_poly.length + b->data.zz_poly.length - 1;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    mpz_t aa, bb, cc, t, base;

    if (r == NULL)
        return NULL;
    mpz_inits(aa, bb, cc, t, base, NULL);
    mpz_setbit(base, bits);
    mpz_set_ui(aa, 0);
    for (i = a->data.zz_poly.length; i-- > 0;) {
        mpz_mul_2exp(aa, aa, bits);
        mpz_add(aa, aa, a->data.zz_poly.coeff[i]);
    }
    mpz_set_ui(bb, 0);
    for (i = b->data.zz_poly.length; i-- > 0;) {
        mpz_mul_2exp(bb, bb, bits);
        mpz_add(bb, bb, b->data.zz_poly.coeff[i]);
    }
    mpz_mul(cc, aa, bb);
    for (i = 0; i < n; i++) {
        mpz_fdiv_r_2exp(t, cc, bits);
        if (mpz_tstbit(t, bits - 1))
            mpz_sub(t, t, base);
        mpz_set(r->data.zz_poly.coeff[i], t);
        mpz_sub(cc, cc, t);
        mpz_fdiv_q_2exp(cc, cc, bits);
    }
    mpz_clears(aa, bb, cc, t, base, NULL);
    return r;
}

sc_value *sc_zz_poly_mul_karatsuba(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t k = (an < bn ? an : bn) / 2, i, n = an + bn - 1;
    sc_value a0 = sc_zz_poly_view(a, 0, k);
    sc_value a1 = sc_zz_poly_view(a, k, an - k);
    sc_value b0 = sc_zz_poly_view(b, 0, k);
    sc_value b1 = sc_zz_poly_view(b, k, bn - k);

    sc_value *s = sc_zz_poly_add(ctx, &a0, &a1);
    sc_value *t = sc_zz_poly_add(ctx, &b0, &b1);
    sc_value *z0 = sc_zz_poly_mul(ctx, &a0, &b0);
    sc_value *z2 = sc_zz_poly_mul(ctx, &a1, &b1);
    sc_value *z1 = sc_zz_poly_mul(ctx, s, t);
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n);

    if (s == NULL || t == NULL || z0 == NULL || z1 == NULL || z2 == NULL || r == NULL)
        return sc_value_free_many_null(6, s, t, z0, z1, z2, r);

    sc_zz_poly *rp = &r->data.zz_poly;
    for (i = 0; i < z0->data.zz_poly.length; i++)
        SC_MPZ_ADDEQ(rp->coeff[i], z0->data.zz_poly.coeff[i]);
    for (i = 0; i < z1->data.zz_poly.length; i++)
        SC_MPZ_ADDEQ(rp->coeff[i + k], z1->data.zz_poly.coeff[i]);
    for (i = 0; i < z0->data.zz_poly.length; i++)
        SC_MPZ_SUBEQ(rp->coeff[i + k], z0->data.zz_poly.coeff[i]);
    for (i = 0; i < z2->data.zz_poly.length; i++)
        SC_MPZ_SUBEQ(rp->coeff[i + k], z2->data.zz_poly.coeff[i]);
    for (i = 0; i < z2->data.zz_poly.length; i++)
        SC_MPZ_ADDEQ(rp->coeff[i + 2 * k], z2->data.zz_poly.coeff[i]);

    sc_value_free_many(5, s, t, z0, z1, z2);
    return r;
}

sc_value *sc_zz_poly_mul_toom3(sc_context *ctx, sc_zz_poly_toom3_ws *w)
{
    size_t i, j;

    for (i = 0; i < 2 * w->m; i++) {
        mpz_add(SC_T3_EVAL(w, i, 2), SC_T3_BLOCK(w, i, 0), SC_T3_BLOCK(w, i, 2));
        mpz_add(SC_T3_EVAL(w, i, 1), SC_T3_EVAL(w, i, 2), SC_T3_BLOCK(w, i, 1));
        mpz_sub(SC_T3_EVAL(w, i, 2), SC_T3_EVAL(w, i, 2), SC_T3_BLOCK(w, i, 1));
        mpz_add(SC_T3_EVAL(w, i, 3), SC_T3_EVAL(w, i, 1), SC_T3_BLOCK(w, i, 1));
        mpz_addmul_ui(SC_T3_EVAL(w, i, 3), SC_T3_BLOCK(w, i, 2), 3);
    }

    for (j = 0; j < 5; j++)
        w->prod[j] = sc_zz_poly_mul(ctx, w->point[0][j], w->point[1][j]);
    if (w->prod[0] == NULL || w->prod[1] == NULL || w->prod[2] == NULL ||
        w->prod[3] == NULL || w->prod[4] == NULL)
        return sc_zz_poly_toom3_ws_abort(w);

    for (i = 0; i < w->q; i++) {
        mpz_sub(SC_T3_PROD(w, 3, i), SC_T3_PROD(w, 3, i), SC_T3_PROD(w, 1, i));
        mpz_submul_ui(SC_T3_PROD(w, 3, i), SC_T3_PROD_OR_ZERO(w, 4, i), 15);
        mpz_sub(SC_T3_PROD(w, 1, i), SC_T3_PROD(w, 1, i), SC_T3_PROD(w, 2, i));
        mpz_fdiv_q_2exp(SC_T3_PROD(w, 1, i), SC_T3_PROD(w, 1, i), 1);
        mpz_sub(SC_T3_PROD(w, 2, i), SC_T3_PROD(w, 2, i), SC_T3_PROD_OR_ZERO(w, 0, i));
        mpz_sub(SC_T3_PROD(w, 2, i), SC_T3_PROD(w, 2, i), SC_T3_PROD_OR_ZERO(w, 4, i));
        mpz_add(SC_T3_PROD(w, 2, i), SC_T3_PROD(w, 2, i), SC_T3_PROD(w, 1, i));
        mpz_sub(SC_T3_PROD(w, 3, i), SC_T3_PROD(w, 3, i), SC_T3_PROD(w, 1, i));
        mpz_submul_ui(SC_T3_PROD(w, 3, i), SC_T3_PROD(w, 2, i), 3);
        mpz_divexact_ui(SC_T3_PROD(w, 3, i), SC_T3_PROD(w, 3, i), 6);
        mpz_sub(SC_T3_PROD(w, 1, i), SC_T3_PROD(w, 1, i), SC_T3_PROD(w, 3, i));
    }

    for (j = 0; j < 5; j++)
        for (i = 0; i < w->prod[j]->data.zz_poly.length; i++)
            SC_T3_ADDTO_RESULT(w, j, i);

    return sc_zz_poly_toom3_ws_finish(w);
}

sc_value *sc_zz_poly_pow_binary(sc_context *ctx, const sc_value *a, unsigned long e)
{
    unsigned long bit = 1;
    sc_value *r, *t;

    if (e == 0) {
        r = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
        if (r != NULL)
            mpz_set_ui(SC_ZP(r, 0), 1);
        return r;
    }
    r = sc_value_copy_checked(ctx, a);
    if (r == NULL)
        return NULL;
    while (bit <= e / 2)
        bit <<= 1;
    for (bit >>= 1; bit != 0; bit >>= 1) {
        t = sc_zz_poly_mul(ctx, r, r);
        sc_value_free(r);
        if (t == NULL)
            return NULL;
        r = t;
        if (e & bit) {
            t = sc_zz_poly_mul(ctx, r, a);
            sc_value_free(r);
            if (t == NULL)
                return NULL;
            r = t;
        }
    }
    return r;
}

sc_value *sc_zz_poly_reverse_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    size_t i, j;

    if (r == NULL)
        return NULL;
    for (i = 0; i < n; i++) {
        j = n - 1 - i;
        if (j < a->data.zz_poly.length)
            mpz_set(SC_ZP(r, i), SC_ZP(a, j));
    }
    sc_zz_poly_normalize(r);
    return r;
}

sc_value *sc_zz_poly_degree_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    if (SC_ZN(a) == 0)
        mpz_set_si(r->data.z, -1);
    else
        mpz_set_ui(r->data.z, SC_ZN(a) - 1);
    return r;
}

sc_value *sc_zz_poly_leading_coefficient_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r != NULL && SC_ZN(a) != 0)
        mpz_set(r->data.z, SC_ZLC(a));
    return r;
}

sc_value *sc_zz_poly_constant_coefficient_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r != NULL && SC_ZN(a) != 0)
        mpz_set(r->data.z, SC_ZP(a, 0));
    return r;
}

sc_value *sc_zz_poly_coeff_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r != NULL && n < SC_ZN(a))
        mpz_set(r->data.z, SC_ZP(a, n));
    return r;
}

sc_value *sc_zz_poly_truncate_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    size_t i, m = n < SC_ZN(a) ? n : SC_ZN(a);
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, m);

    if (r == NULL)
        return NULL;
    for (i = 0; i < m; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(a, i));
    sc_zz_poly_normalize(r);
    return r;
}

sc_value *sc_zz_poly_shift_left_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    size_t i;
    sc_value *r;

    if (SC_ZN(a) == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (n > (size_t)-1 - SC_ZN(a)) {
        sc_set_error(ctx, "polynomial shift is too large");
        return NULL;
    }
    r = sc_value_new_zz_poly_checked(ctx, a->parent, SC_ZN(a) + n);
    if (r == NULL)
        return NULL;
    for (i = 0; i < SC_ZN(a); i++)
        mpz_set(SC_ZP(r, i + n), SC_ZP(a, i));
    return r;
}

sc_value *sc_zz_poly_shift_right_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    size_t i, m = n < SC_ZN(a) ? SC_ZN(a) - n : 0;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, m);

    if (r == NULL)
        return NULL;
    for (i = 0; i < m; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(a, i + n));
    return r;
}

sc_value *sc_zz_poly_height_impl(sc_context *ctx, const sc_value *a)
{
    size_t i;
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    for (i = 0; i < SC_ZN(a); i++)
        if (mpz_cmpabs(SC_ZP(a, i), r->data.z) > 0)
            mpz_abs(r->data.z, SC_ZP(a, i));
    return r;
}

sc_value *sc_zz_poly_max_abs_bits_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_set_ui(r->data.z, sc_zz_poly_max_abs_bits_raw(a));
    return r;
}

size_t sc_zz_poly_max_abs_bits_raw(const sc_value *a)
{
    size_t i, bits, max_bits = 0;

    for (i = 0; i < SC_ZN(a); i++) {
        if (mpz_sgn(SC_ZP(a, i)) == 0)
            continue;
        bits = mpz_sizeinbase(SC_ZP(a, i), 2);
        if (bits > max_bits)
            max_bits = bits;
    }
    return max_bits;
}

sc_value *sc_zz_poly_inflate_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    size_t i, m;
    sc_value *r;

    if (SC_ZN(a) == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (n == 0) {
        r = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
        if (r == NULL)
            return NULL;
        for (i = 0; i < SC_ZN(a); i++)
            SC_MPZ_ADDEQ(SC_ZP(r, 0), SC_ZP(a, i));
        sc_zz_poly_normalize(r);
        return r;
    }
    if (SC_ZN(a) - 1 > ((size_t)-2) / n) {
        sc_set_error(ctx, "polynomial inflation is too large");
        return NULL;
    }
    m = (SC_ZN(a) - 1) * n + 1;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, m);
    if (r == NULL)
        return NULL;
    for (i = 0; i < SC_ZN(a); i++)
        mpz_set(SC_ZP(r, i * n), SC_ZP(a, i));
    return r;
}

sc_value *sc_zz_poly_deflation_impl(sc_context *ctx, const sc_value *a)
{
    size_t i, d = 0;
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL || SC_ZN(a) == 0)
        return r;
    if (SC_ZN(a) == 1) {
        mpz_set_ui(r->data.z, 1);
        return r;
    }
    for (i = 1; i < SC_ZN(a); i++) {
        size_t x, y, t;

        if (mpz_sgn(SC_ZP(a, i)) == 0)
            continue;
        x = d;
        y = i;
        while (y != 0) {
            t = x % y;
            x = y;
            y = t;
        }
        d = x;
        if (d == 1)
            break;
    }
    mpz_set_ui(r->data.z, d);
    return r;
}

sc_value *sc_zz_poly_deflate_impl(sc_context *ctx, const sc_value *a, size_t n)
{
    size_t i, m;
    sc_value *r;

    if (n == 0) {
        sc_set_error(ctx, "polynomial deflation must be positive");
        return NULL;
    }
    for (i = 0; i < SC_ZN(a); i++)
        if (i % n != 0 && mpz_sgn(SC_ZP(a, i)) != 0) {
            sc_set_error(ctx, "polynomial is not deflatable by this factor");
            return NULL;
        }
    m = SC_ZN(a) == 0 ? 0 : (SC_ZN(a) - 1) / n + 1;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, m);
    if (r == NULL)
        return NULL;
    for (i = 0; i < m; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(a, i * n));
    return r;
}

sc_value *sc_zz_poly_quo_classical(sc_context *ctx,
                                    const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length, i, k;
    size_t qn = an >= bn ? an - bn + 1 : 0;
    sc_value *q = sc_value_new_zz_poly_checked(ctx, a->parent, qn);

    if (q == NULL)
        return NULL;
    for (i = 0; i < qn; i++)
        mpz_set(SC_ZP(q, i), SC_ZP(a, bn - 1 + i));
    for (k = qn; k-- > 0;) {
        if (!mpz_divisible_p(SC_ZP(q, k), SC_ZP(b, bn - 1))) {
            sc_set_error(ctx, "polynomial quotient is not in ZZ[x]");
            return sc_value_free_many_null(1, q);
        }
        mpz_divexact(SC_ZP(q, k), SC_ZP(q, k), SC_ZP(b, bn - 1));
        i = bn - 1 > k ? bn - 1 - k : 0;
        for (; i < bn - 1; i++)
            SC_MPZ_SUBMUL(SC_ZP(q, k + i - bn + 1), SC_ZP(q, k), SC_ZP(b, i));
    }
    return q;
}

sc_value *sc_zz_poly_series_quo_classical(sc_context *ctx, const sc_value *a,
                                           const sc_value *b, size_t n)
{
    sc_value *q = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    size_t i, j;

    if (q == NULL)
        return NULL;
    for (i = 0; i < n; i++) {
        if (i < a->data.zz_poly.length)
            mpz_set(SC_ZP(q, i), SC_ZP(a, i));
        for (j = 1; j <= i && j < b->data.zz_poly.length; j++)
            SC_MPZ_SUBMUL(SC_ZP(q, i), SC_ZP(b, j), SC_ZP(q, i - j));
        if (!mpz_divisible_p(SC_ZP(q, i), SC_ZP(b, 0))) {
            sc_set_error(ctx, "polynomial quotient is not in ZZ[x]");
            return sc_value_free_many_null(1, q);
        }
        mpz_divexact(SC_ZP(q, i), SC_ZP(q, i), SC_ZP(b, 0));
    }
    sc_zz_poly_normalize(q);
    return q;
}

sc_value *sc_zz_poly_series_quo_dc(sc_context *ctx, const sc_value *a,
                                    const sc_value *b, size_t n)
{
    size_t k = n / 2, h = n - k, i;
    sc_value av = sc_zz_poly_view(a, 0, k), bv = sc_zz_poly_view(b, 0, n);
    sc_value *q0, *p, *e, *q1, *r;

    if (n <= SC_SERIES_QUO_DC_CUTOFF)
        return sc_zz_poly_series_quo_classical(ctx, a, b, n);
    q0 = sc_zz_poly_series_quo_dc(ctx, &av, &bv, k);
    p = q0 == NULL ? NULL : sc_zz_poly_mulmid(ctx, &bv, q0, k, h);
    e = sc_value_new_zz_poly_checked(ctx, a->parent, h);
    if (q0 == NULL || p == NULL || e == NULL)
        return sc_value_free_many_null(3, q0, p, e);
    for (i = 0; i < h; i++) {
        if (k + i < a->data.zz_poly.length)
            mpz_set(SC_ZP(e, i), SC_ZP(a, k + i));
        if (i < p->data.zz_poly.length)
            SC_MPZ_SUBEQ(SC_ZP(e, i), SC_ZP(p, i));
    }
    q1 = sc_zz_poly_series_quo_dc(ctx, e, &bv, h);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (q1 == NULL || r == NULL)
        return sc_value_free_many_null(5, q0, p, e, q1, r);
    for (i = 0; i < q0->data.zz_poly.length; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(q0, i));
    for (i = 0; i < q1->data.zz_poly.length; i++)
        mpz_set(SC_ZP(r, k + i), SC_ZP(q1, i));
    sc_zz_poly_normalize(r);
    sc_value_free_many(4, q0, p, e, q1);
    return r;
}

sc_value *sc_zz_poly_inv_series_classical(sc_context *ctx, const sc_value *a,
                                               size_t n)
{
    sc_value *g = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    size_t i, j;

    if (g == NULL || n == 0)
        return g;
    mpz_set(SC_ZP(g, 0), SC_ZP(a, 0));
    for (i = 1; i < n; i++) {
        for (j = 1; j <= i && j < a->data.zz_poly.length; j++)
            mpz_addmul(SC_ZP(g, i), SC_ZP(a, j), SC_ZP(g, i - j));
        if (mpz_sgn(SC_ZP(a, 0)) > 0)
            mpz_neg(SC_ZP(g, i), SC_ZP(g, i));
    }
    sc_zz_poly_normalize(g);
    return g;
}

sc_value *sc_zz_poly_inv_series_newton(sc_context *ctx, const sc_value *a,
                                         size_t n)
{
    size_t m = (n + 1) / 2, h = n - m, i;
    sc_value *g, *e, *c, *r;

    if (n <= SC_INV_SERIES_NEWTON_CUTOFF)
        return sc_zz_poly_inv_series_classical(ctx, a, n);
    g = sc_zz_poly_inv_series_newton(ctx, a, m);
    e = g == NULL ? NULL : sc_zz_poly_mulmid(ctx, a, g, m, h);
    c = e == NULL ? NULL : sc_zz_poly_mullow(ctx, g, e, h);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (g == NULL || e == NULL || c == NULL || r == NULL)
        return sc_value_free_many_null(4, g, e, c, r);
    for (i = 0; i < g->data.zz_poly.length; i++)
        mpz_set(SC_ZP(r, i), SC_ZP(g, i));
    for (i = 0; i < c->data.zz_poly.length; i++)
        mpz_neg(SC_ZP(r, m + i), SC_ZP(c, i));
    sc_zz_poly_normalize(r);
    sc_value_free_many(3, g, e, c);
    return r;
}

sc_value *sc_zz_poly_series_quo_preinv(sc_context *ctx, const sc_value *a,
                                        const sc_value *binv, size_t n)
{
    return sc_zz_poly_mullow(ctx, a, binv, n);
}

sc_value *sc_zz_poly_series_quo_newton(sc_context *ctx, const sc_value *a,
                                        const sc_value *b, size_t n)
{
    size_t m = (n + 1) / 2, h = n - m, i;
    sc_value *g, *q0, *p, *e, *q1, *q;

    if (n <= SC_SERIES_QUO_NEWTON_CUTOFF)
        return sc_zz_poly_series_quo_classical(ctx, a, b, n);
    g = sc_zz_poly_inv_series(ctx, b, m);
    q0 = g == NULL ? NULL : sc_zz_poly_mullow(ctx, a, g, m);
    p = q0 == NULL ? NULL : sc_zz_poly_mulmid(ctx, b, q0, m, h);
    e = sc_value_new_zz_poly_checked(ctx, a->parent, h);
    if (g == NULL || q0 == NULL || p == NULL || e == NULL)
        return sc_value_free_many_null(4, g, q0, p, e);
    for (i = 0; i < h; i++) {
        if (m + i < a->data.zz_poly.length)
            mpz_set(SC_ZP(e, i), SC_ZP(a, m + i));
        if (i < p->data.zz_poly.length)
            SC_MPZ_SUBEQ(SC_ZP(e, i), SC_ZP(p, i));
    }
    q1 = sc_zz_poly_mullow(ctx, e, g, h);
    q = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (q1 == NULL || q == NULL)
        return sc_value_free_many_null(6, g, q0, p, e, q1, q);
    for (i = 0; i < q0->data.zz_poly.length; i++)
        mpz_set(SC_ZP(q, i), SC_ZP(q0, i));
    for (i = 0; i < q1->data.zz_poly.length; i++)
        mpz_set(SC_ZP(q, m + i), SC_ZP(q1, i));
    sc_zz_poly_normalize(q);
    sc_value_free_many(5, g, q0, p, e, q1);
    return q;
}

sc_value *sc_zz_poly_quo_dc(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t qn = an - bn + 1, rn = bn < qn ? bn : qn;
    sc_value av = sc_zz_poly_view(a, bn - 1, qn);
    sc_value bv = sc_zz_poly_view(b, bn - rn, rn);
    sc_value *ar = sc_zz_poly_reverse_impl(ctx, &av, qn);
    sc_value *br = sc_zz_poly_reverse_impl(ctx, &bv, rn);
    sc_value *qr = ar != NULL && br != NULL ?
                   sc_zz_poly_series_quo_dc(ctx, ar, br, qn) : NULL;
    sc_value *q = qr != NULL ? sc_zz_poly_reverse_impl(ctx, qr, qn) : NULL;

    sc_value_free_many(3, ar, br, qr);
    return q;
}

sc_value *sc_zz_poly_preinverse_newton(sc_context *ctx, const sc_value *b,
                                            size_t n)
{
    size_t bn = b->data.zz_poly.length, rn = bn < n ? bn : n;
    sc_value bv = sc_zz_poly_view(b, bn - rn, rn);
    sc_value *br = sc_zz_poly_reverse_impl(ctx, &bv, rn);
    sc_value *inv = br == NULL ? NULL : sc_zz_poly_inv_series(ctx, br, n);

    sc_value_free(br);
    return inv;
}

sc_value *sc_zz_poly_quo_preinv(sc_context *ctx, const sc_value *a,
                                 const sc_value *b, const sc_value *binv)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t qn = an >= bn ? an - bn + 1 : 0;
    sc_value av;
    sc_value *ar, *qr, *q;

    if (qn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    av = sc_zz_poly_view(a, bn - 1, qn);
    ar = sc_zz_poly_reverse_impl(ctx, &av, qn);
    qr = ar == NULL ? NULL : sc_zz_poly_series_quo_preinv(ctx, ar, binv, qn);
    q = qr == NULL ? NULL : sc_zz_poly_reverse_impl(ctx, qr, qn);
    sc_value_free_many(2, ar, qr);
    return q;
}

sc_value *sc_zz_poly_quo_newton(sc_context *ctx, const sc_value *a,
                                 const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t qn = an >= bn ? an - bn + 1 : 0, rn = bn < qn ? bn : qn;
    sc_value av, bv;
    sc_value *ar, *br, *qr, *q;

    if (qn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    av = sc_zz_poly_view(a, bn - 1, qn);
    bv = sc_zz_poly_view(b, bn - rn, rn);
    ar = sc_zz_poly_reverse_impl(ctx, &av, qn);
    br = sc_zz_poly_reverse_impl(ctx, &bv, rn);
    qr = ar != NULL && br != NULL ?
         sc_zz_poly_series_quo_newton(ctx, ar, br, qn) : NULL;
    q = qr == NULL ? NULL : sc_zz_poly_reverse_impl(ctx, qr, qn);
    sc_value_free_many(3, ar, br, qr);
    return q;
}

static sc_value *sc_zz_poly_remainder_from_quotient(sc_context *ctx,
                                                      const sc_value *a,
                                                      const sc_value *b,
                                                      const sc_value *q)
{
    size_t rn = b->data.zz_poly.length - 1, i;
    sc_value *p = sc_zz_poly_mullow(ctx, b, q, rn);
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, rn);

    if (p == NULL || r == NULL)
        return sc_value_free_many_null(2, p, r);
    for (i = 0; i < rn; i++) {
        if (i < a->data.zz_poly.length)
            mpz_set(SC_ZP(r, i), SC_ZP(a, i));
        if (i < p->data.zz_poly.length)
            SC_MPZ_SUBEQ(SC_ZP(r, i), SC_ZP(p, i));
    }
    sc_zz_poly_normalize(r);
    sc_value_free(p);
    return r;
}

sc_value *sc_zz_poly_divrem_newton(sc_context *ctx, const sc_value *a,
                                    const sc_value *b)
{
    sc_value *q = sc_zz_poly_quo_newton(ctx, a, b);
    sc_value *r = q == NULL ? NULL : sc_zz_poly_remainder_from_quotient(ctx, a, b, q);

    if (q == NULL || r == NULL)
        return sc_value_free_many_null(2, q, r);
    return sc_value_new_pair_take_checked(ctx, q, r);
}

sc_value *sc_zz_poly_quo_bidirectional(sc_context *ctx,
                                          const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t qn = an - bn + 1, ln = (qn + 1) / 2, hn = qn - ln, i;
    size_t rn = bn < hn ? bn : hn;
    sc_value av = sc_zz_poly_view(a, an - hn, hn);
    sc_value bv = sc_zz_poly_view(b, bn - rn, rn);
    sc_value *lo, *ar, *br, *hir, *hi, *q;

    if (qn <= SC_BIDIR_QUO_CUTOFF)
        return sc_zz_poly_quo_classical(ctx, a, b);
    lo = sc_zz_poly_series_quo_dc(ctx, a, b, ln);
    ar = sc_zz_poly_reverse_impl(ctx, &av, hn);
    br = sc_zz_poly_reverse_impl(ctx, &bv, rn);
    hir = lo != NULL && ar != NULL && br != NULL ?
          sc_zz_poly_series_quo_dc(ctx, ar, br, hn) : NULL;
    hi = hir != NULL ? sc_zz_poly_reverse_impl(ctx, hir, hn) : NULL;
    q = sc_value_new_zz_poly_checked(ctx, a->parent, qn);
    if (lo == NULL || ar == NULL || br == NULL || hir == NULL || hi == NULL || q == NULL)
        return sc_value_free_many_null(6, lo, ar, br, hir, hi, q);
    for (i = 0; i < lo->data.zz_poly.length; i++)
        mpz_set(SC_ZP(q, i), SC_ZP(lo, i));
    for (i = 0; i < hi->data.zz_poly.length; i++)
        mpz_set(SC_ZP(q, ln + i), SC_ZP(hi, i));
    sc_zz_poly_normalize(q);
    sc_value_free_many(5, lo, ar, br, hir, hi);
    return q;
}

sc_value *sc_zz_poly_quo_mulders_balanced(sc_context *ctx,
                                           const sc_value *a, const sc_value *b)
{
    size_t n = b->data.zz_poly.length, n1 = (n + 1) / 2, n2 = n - n1;
    size_t delta = n1 - n2, i;
    sc_value b1 = sc_zz_poly_view(b, n2, n1), b2 = sc_zz_poly_view(b, 0, n2);
    sc_value b3 = sc_zz_poly_view(b, n1, n2);
    sc_value *u, *qr1, *h, *t, *q2, *q;

    if (n <= SC_MULDERS_QUO_CUTOFF)
        return sc_zz_poly_quo_classical(ctx, a, b);
    u = sc_value_new_zz_poly_checked(ctx, a->parent, 2 * n1 - 1);
    if (u == NULL)
        return NULL;
    for (i = 0; i < n1; i++)
        mpz_set(SC_ZP(u, n1 - 1 + i), SC_ZP(a, n + n2 - 1 + i));
    qr1 = sc_zz_poly_divrem_full(ctx, u, &b1);
    if (qr1 == NULL)
        return sc_value_free_many_null(1, u);
    h = sc_zz_poly_mulmid(ctx, &b2, qr1->data.pair.first, n1 - 1, n2);
    t = sc_value_new_zz_poly_checked(ctx, a->parent, 2 * n2 - 1);
    if (h == NULL || t == NULL)
        return sc_value_free_many_null(4, u, qr1, h, t);
    for (i = 0; i < n2; i++)
        mpz_set(SC_ZP(t, n2 - 1 + i), SC_ZP(a, n - 1 + i));
    for (i = 0; i < qr1->data.pair.second->data.zz_poly.length; i++)
        SC_MPZ_ADDEQ(SC_ZP(t, n2 - delta + i), SC_ZP(qr1->data.pair.second, i));
    for (i = 0; i < h->data.zz_poly.length; i++)
        SC_MPZ_SUBEQ(SC_ZP(t, n2 - 1 + i), SC_ZP(h, i));
    sc_zz_poly_normalize(t);
    q2 = sc_zz_poly_quo_mulders(ctx, t, &b3);
    q = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (q2 == NULL || q == NULL)
        return sc_value_free_many_null(6, u, qr1, h, t, q2, q);
    for (i = 0; i < q2->data.zz_poly.length; i++)
        mpz_set(SC_ZP(q, i), SC_ZP(q2, i));
    for (i = 0; i < qr1->data.pair.first->data.zz_poly.length; i++)
        mpz_set(SC_ZP(q, n2 + i), SC_ZP(qr1->data.pair.first, i));
    sc_value_free_many(5, u, qr1, h, t, q2);
    return q;
}

sc_value *sc_zz_poly_divrem_classical(sc_context *ctx,
                                       const sc_value *a, const sc_value *b)
{
    size_t bn = b->data.zz_poly.length, k, i;
    size_t qn = a->data.zz_poly.length >= bn ? a->data.zz_poly.length - bn + 1 : 0;
    sc_value *q = sc_value_new_zz_poly_checked(ctx, a->parent, qn);
    sc_value *r = sc_value_copy_checked(ctx, a);

    if (q == NULL || r == NULL)
        return sc_value_free_many_null(2, q, r);
    while (r->data.zz_poly.length >= bn) {
        k = r->data.zz_poly.length - bn;
        if (!mpz_divisible_p(SC_ZP(r, r->data.zz_poly.length - 1), SC_ZP(b, bn - 1))) {
            sc_set_error(ctx, "polynomial quotient and remainder are not in ZZ[x]");
            return sc_value_free_many_null(2, q, r);
        }
        mpz_divexact(SC_ZP(q, k), SC_ZP(r, r->data.zz_poly.length - 1),
                     SC_ZP(b, bn - 1));
        for (i = 0; i < bn; i++)
            SC_MPZ_SUBMUL(SC_ZP(r, k + i), SC_ZP(q, k), SC_ZP(b, i));
        sc_zz_poly_normalize(r);
    }
    return sc_value_new_pair_take_checked(ctx, q, r);
}

sc_value *sc_zz_poly_divrem_dc(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *q = sc_zz_poly_quo_dc(ctx, a, b);
    sc_value *r = q == NULL ? NULL : sc_zz_poly_remainder_from_quotient(ctx, a, b, q);

    if (q == NULL || r == NULL)
        return sc_value_free_many_null(2, q, r);
    return sc_value_new_pair_take_checked(ctx, q, r);
}

sc_value *sc_zz_poly_divrem_mulders(sc_context *ctx,
                                     const sc_value *a, const sc_value *b)
{
    sc_value *q = sc_zz_poly_quo_mulders(ctx, a, b);
    sc_value *r = q == NULL ? NULL : sc_zz_poly_remainder_from_quotient(ctx, a, b, q);

    if (q == NULL || r == NULL)
        return sc_value_free_many_null(2, q, r);
    return sc_value_new_pair_take_checked(ctx, q, r);
}

sc_value *sc_zz_poly_pseudodiv_impl(sc_context *ctx,
                                     const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length, i, k;
    size_t d = an >= bn ? an - bn + 1 : 0;
    sc_value *q = sc_value_new_zz_poly_checked(ctx, a->parent, d);
    sc_value *r = sc_value_copy_checked(ctx, a);
    mpz_t scale;

    if (q == NULL || r == NULL)
        return sc_value_free_many_null(2, q, r);
    mpz_init(scale);
    while (r->data.zz_poly.length >= bn) {
        k = r->data.zz_poly.length - bn;
        for (i = 0; i < q->data.zz_poly.length; i++)
            SC_MPZ_MULEQ(SC_ZP(q, i), SC_ZP(b, bn - 1));
        mpz_set(SC_ZP(q, k), SC_ZP(r, r->data.zz_poly.length - 1));
        for (i = 0; i < r->data.zz_poly.length; i++)
            SC_MPZ_MULEQ(SC_ZP(r, i), SC_ZP(b, bn - 1));
        for (i = 0; i < bn; i++)
            SC_MPZ_SUBMUL(SC_ZP(r, k + i), SC_ZP(q, k), SC_ZP(b, i));
        sc_zz_poly_normalize(r);
        d--;
    }
    if (d != 0) {
        mpz_pow_ui(scale, SC_ZP(b, bn - 1), (unsigned long)d);
        for (i = 0; i < q->data.zz_poly.length; i++)
            SC_MPZ_MULEQ(SC_ZP(q, i), scale);
        for (i = 0; i < r->data.zz_poly.length; i++)
            SC_MPZ_MULEQ(SC_ZP(r, i), scale);
    }
    mpz_clear(scale);
    return sc_value_new_pair_take_checked(ctx, q, r);
}

sc_value *sc_zz_poly_pseudorem_classical(sc_context *ctx,
                                            const sc_value *a, const sc_value *b)
{
    size_t an = SC_ZN(a), bn = SC_ZN(b), d = an >= bn ? an - bn + 1 : 0, i, k;
    sc_value *r;
    mpz_t c, scale;

    if (bn == 0) {
        sc_set_error(ctx, "polynomial pseudo-remainder by zero");
        return NULL;
    }
    r = sc_value_copy_checked(ctx, a);
    if (r == NULL)
        return NULL;
    mpz_inits(c, scale, NULL);
    while (SC_ZN(r) >= bn) {
        k = SC_ZN(r) - bn;
        mpz_set(c, SC_ZLC(r));
        for (i = 0; i < SC_ZN(r); i++)
            SC_MPZ_MULEQ(SC_ZP(r, i), SC_ZLC(b));
        for (i = 0; i < bn; i++)
            SC_MPZ_SUBMUL(SC_ZP(r, k + i), c, SC_ZP(b, i));
        sc_zz_poly_normalize(r);
        d--;
    }
    if (d != 0) {
        mpz_pow_ui(scale, SC_ZLC(b), (unsigned long)d);
        for (i = 0; i < SC_ZN(r); i++)
            SC_MPZ_MULEQ(SC_ZP(r, i), scale);
    }
    mpz_clears(c, scale, NULL);
    return r;
}

sc_value *sc_zz_poly_pseudorem_fast(sc_context *ctx,
                                    const sc_value *a, const sc_value *b)
{
    sc_value *qr = sc_zz_poly_pseudodiv_fast(ctx, a, b);

    return sc_value_pair_take(qr, 1);
}

sc_value *sc_zz_poly_scalar_divexact_impl(sc_context *ctx, const sc_value *a,
                                           const sc_value *b)
{
    size_t i, n = a->data.zz_poly.length;
    sc_value *r;

    if (mpz_sgn(b->data.z) == 0) {
        sc_set_error(ctx, "polynomial scalar division by zero");
        return NULL;
    }
    for (i = 0; i < n; i++)
        if (!mpz_divisible_p(SC_ZP(a, i), b->data.z)) {
            sc_set_error(ctx, "polynomial scalar quotient is not exact in ZZ[x]");
            return NULL;
        }
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (r == NULL)
        return NULL;
    for (i = 0; i < n; i++)
        mpz_divexact(SC_ZP(r, i), SC_ZP(a, i), b->data.z);
    sc_zz_poly_normalize(r);
    return r;
}

sc_value *sc_zz_poly_content_impl(sc_context *ctx, const sc_value *a)
{
    size_t i, n = a->data.zz_poly.length;
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL || n == 0)
        return r;
    mpz_abs(r->data.z, SC_ZP(a, 0));
    for (i = 1; i < n && mpz_cmp_ui(r->data.z, 1) != 0; i++)
        mpz_gcd(r->data.z, r->data.z, SC_ZP(a, i));
    return r;
}

sc_value *sc_zz_poly_primitive_part_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *c = sc_zz_poly_content_impl(ctx, a), *r;

    if (c == NULL)
        return NULL;
    if (mpz_sgn(c->data.z) == 0) {
        sc_value_free(c);
        return sc_value_copy_checked(ctx, a);
    }
    r = sc_zz_poly_scalar_divexact_impl(ctx, a, c);
    sc_value_free(c);
    return r;
}

sc_value *sc_zz_poly_gcd_pseudo_impl(sc_context *ctx, const sc_value *a,
                                      const sc_value *b)
{
    sc_value *ca = sc_zz_poly_content_impl(ctx, a);
    sc_value *cb = sc_zz_poly_content_impl(ctx, b), *c, *u, *v, *w, *p, *r;

    c = ca && cb ? sc_zz_gcd(ctx, ca, cb) : NULL;
    u = c ? sc_zz_poly_primitive_part_impl(ctx, a) : NULL;
    v = u ? sc_zz_poly_primitive_part_impl(ctx, b) : NULL;
    if (v == NULL)
        return sc_value_free_many_null(5, ca, cb, c, u, v);
    if (u->data.zz_poly.length < v->data.zz_poly.length) {
        p = u;
        u = v;
        v = p;
    }
    while (v->data.zz_poly.length != 0) {
        w = sc_zz_poly_pseudorem(ctx, u, v);
        p = w ? sc_zz_poly_primitive_part_impl(ctx, w) : NULL;
        sc_value_free_many(2, u, w);
        if (p == NULL)
            return sc_value_free_many_null(4, ca, cb, c, v);
        u = v;
        v = p;
    }
    sc_value_free(v);
    if (u->data.zz_poly.length != 0 && mpz_sgn(SC_ZP(u, u->data.zz_poly.length - 1)) < 0) {
        p = sc_zz_poly_neg(ctx, u);
        sc_value_free(u);
        u = p;
    }
    r = u ? sc_zz_poly_scalar_mul_impl(ctx, u, c) : NULL;
    sc_value_free_many(4, ca, cb, c, u);
    return r;
}

sc_value *sc_zz_poly_inv_series_scaled_newton(sc_context *ctx, const sc_value *a,
                                                    size_t n)
{
    size_t m = 1;
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
    sc_value *p, *t, *u;
    mpz_t c;

    if (r == NULL || n == 0 || SC_ZN(a) == 0 || mpz_sgn(SC_ZP(a, 0)) == 0) {
        sc_value_free(r);
        sc_set_error(ctx, "scaled series inverse needs nonzero constant term");
        return NULL;
    }
    mpz_init(c);
    mpz_set_ui(SC_ZP(r, 0), 1);
    while (m < n) {
        p = sc_zz_poly_mullow(ctx, a, r, 2 * m);
        t = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
        if (p == NULL || t == NULL)
            return sc_value_free_many_null(3, r, p, t), mpz_clear(c), NULL;
        mpz_pow_ui(c, SC_ZP(a, 0), (unsigned long)m);
        mpz_mul_ui(SC_ZP(t, 0), c, 2);
        u = sc_zz_poly_sub(ctx, t, p);
        sc_value_free_many(2, t, p);
        t = u ? sc_zz_poly_mullow(ctx, r, u, 2 * m) : NULL;
        sc_value_free_many(2, r, u);
        if (t == NULL)
            return mpz_clear(c), NULL;
        r = t;
        m *= 2;
    }
    mpz_clear(c);
    return r;
}

sc_value *sc_zz_poly_pseudodiv_fast(sc_context *ctx, const sc_value *a,
                                    const sc_value *b)
{
    size_t an = SC_ZN(a), bn = SC_ZN(b), d, p = 1, brn, i;
    sc_value av, bv;
    sc_value *ar, *br, *inv, *qr, *den, *qrs, *q, *qb, *r, *pair;
    mpz_t scale;

    if (bn == 0) {
        sc_set_error(ctx, "polynomial pseudo-division by zero");
        return NULL;
    }
    if (an < bn) {
        q = sc_value_new_zz_poly_checked(ctx, a->parent, 0);
        r = sc_value_copy_checked(ctx, a);
        return q && r ? sc_value_new_pair_take_checked(ctx, q, r) :
                        sc_value_free_many_null(2, q, r);
    }
    d = an - bn + 1;
    while (p < d)
        p *= 2;
    brn = bn < p ? bn : p;
    av = sc_zz_poly_view(a, an - d, d);
    bv = sc_zz_poly_view(b, bn - brn, brn);
    ar = sc_zz_poly_reverse_impl(ctx, &av, d);
    br = sc_zz_poly_reverse_impl(ctx, &bv, brn);
    inv = br ? sc_zz_poly_inv_series_scaled_newton(ctx, br, p) : NULL;
    qr = ar && inv ? sc_zz_poly_mullow(ctx, ar, inv, d) : NULL;
    den = sc_value_new_zz_checked(ctx);
    if (qr == NULL || den == NULL)
        return sc_value_free_many_null(5, ar, br, inv, qr, den);
    mpz_pow_ui(den->data.z, SC_ZLC(b), (unsigned long)(p - d));
    qrs = sc_zz_poly_scalar_divexact_impl(ctx, qr, den);
    q = qrs ? sc_zz_poly_reverse_impl(ctx, qrs, d) : NULL;
    qb = q && bn > 1 ? sc_zz_poly_mullow(ctx, q, b, bn - 1) :
                       sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    r = qb ? sc_value_new_zz_poly_checked(ctx, a->parent, bn - 1) : NULL;
    mpz_init(scale);
    mpz_pow_ui(scale, SC_ZLC(b), (unsigned long)d);
    for (i = 0; r != NULL && i + 1 < bn; i++) {
        mpz_mul(SC_ZP(r, i), scale, SC_ZP(a, i));
        if (i < SC_ZN(qb))
            mpz_sub(SC_ZP(r, i), SC_ZP(r, i), SC_ZP(qb, i));
    }
    mpz_clear(scale);
    if (r != NULL)
        sc_zz_poly_normalize(r);
    pair = q && r ? sc_value_new_pair_take_checked(ctx, q, r) : NULL;
    if (pair != NULL) {
        q = NULL;
        r = NULL;
    }
    sc_value_free_many(9, ar, br, inv, qr, den, qrs, qb, q, r);
    return pair;
}

int sc_zz_poly_mat2_mul(sc_context *ctx, sc_zz_poly_mat2 *r,
                        const sc_zz_poly_mat2 *a, const sc_zz_poly_mat2 *b)
{
    sc_zz_poly_mat2 t = { 0 };
    sc_value *p, *q;

    p = sc_zz_poly_mul(ctx, a->a00, b->a00);
    q = sc_zz_poly_mul(ctx, a->a01, b->a10);
    t.a00 = p && q ? sc_zz_poly_add(ctx, p, q) : NULL;
    sc_value_free_many(2, p, q);

    p = sc_zz_poly_mul(ctx, a->a00, b->a01);
    q = sc_zz_poly_mul(ctx, a->a01, b->a11);
    t.a01 = p && q ? sc_zz_poly_add(ctx, p, q) : NULL;
    sc_value_free_many(2, p, q);

    p = sc_zz_poly_mul(ctx, a->a10, b->a00);
    q = sc_zz_poly_mul(ctx, a->a11, b->a10);
    t.a10 = p && q ? sc_zz_poly_add(ctx, p, q) : NULL;
    sc_value_free_many(2, p, q);

    p = sc_zz_poly_mul(ctx, a->a10, b->a01);
    q = sc_zz_poly_mul(ctx, a->a11, b->a11);
    t.a11 = p && q ? sc_zz_poly_add(ctx, p, q) : NULL;
    sc_value_free_many(2, p, q);

    if (t.a00 == NULL || t.a01 == NULL || t.a10 == NULL || t.a11 == NULL) {
        sc_zz_poly_mat2_clear(&t);
        return 0;
    }
    sc_zz_poly_mat2_move(r, &t);
    return 1;
}

int sc_zz_poly_mat2_apply(sc_context *ctx, sc_value **u, sc_value **v,
                          const sc_zz_poly_mat2 *m, const sc_value *a,
                          const sc_value *b)
{
    sc_value *p0 = sc_zz_poly_mul(ctx, m->a00, a);
    sc_value *p1 = sc_zz_poly_mul(ctx, m->a01, b);
    sc_value *q0 = sc_zz_poly_mul(ctx, m->a10, a);
    sc_value *q1 = sc_zz_poly_mul(ctx, m->a11, b);

    *u = p0 && p1 ? sc_zz_poly_add(ctx, p0, p1) : NULL;
    *v = q0 && q1 ? sc_zz_poly_add(ctx, q0, q1) : NULL;
    sc_value_free_many(4, p0, p1, q0, q1);
    if (*u != NULL && *v != NULL)
        return 1;
    sc_value_free_many(2, *u, *v);
    *u = NULL;
    *v = NULL;
    return 0;
}

int sc_zz_poly_mat2_primitive(sc_context *ctx, sc_zz_poly_mat2 *m)
{
    sc_value *e[4] = { m->a00, m->a01, m->a10, m->a11 };
    sc_value *q[4] = { NULL }, *d;
    mpz_t c;
    size_t i, j;

    mpz_init_set_ui(c, 0);
    for (i = 0; i < 4 && mpz_cmp_ui(c, 1) != 0; i++)
        for (j = 0; j < SC_ZN(e[i]) && mpz_cmp_ui(c, 1) != 0; j++)
            mpz_gcd(c, c, SC_ZP(e[i], j));
    if (mpz_cmp_ui(c, 1) <= 0) {
        mpz_clear(c);
        return 1;
    }
    d = sc_value_new_zz_checked(ctx);
    if (d == NULL) {
        mpz_clear(c);
        return 0;
    }
    mpz_set(d->data.z, c);
    mpz_clear(c);
    for (i = 0; i < 4; i++) {
        q[i] = sc_zz_poly_scalar_divexact_impl(ctx, e[i], d);
        if (q[i] == NULL) {
            sc_value_free_many(5, q[0], q[1], q[2], q[3], d);
            return 0;
        }
    }
    sc_zz_poly_mat2_clear(m);
    m->a00 = q[0];
    m->a01 = q[1];
    m->a10 = q[2];
    m->a11 = q[3];
    sc_value_free(d);
    return 1;
}

int sc_zz_poly_hgcd_pseudo(sc_context *ctx, sc_zz_poly_mat2 *m,
                            const sc_value *a, const sc_value *b)
{
    size_t ad = SC_ZN(a) - 1, mid = (ad + 1) / 2, k;
    sc_value ah, bh, dh, rh, *c = NULL, *d = NULL, *qr = NULL, *q, *rem;
    sc_zz_poly_mat2 r = { 0 }, t = { 0 }, s = { 0 }, u = { 0 };

    if (SC_ZN(b) == 0 || SC_ZN(b) - 1 < mid)
        return sc_zz_poly_mat2_identity(ctx, m, a->parent);
    ah = sc_zz_poly_view(a, mid, SC_ZN(a) - mid);
    bh = sc_zz_poly_view(b, mid, SC_ZN(b) - mid);
    if (!sc_zz_poly_hgcd_pseudo(ctx, &r, &ah, &bh) ||
        !sc_zz_poly_mat2_apply(ctx, &c, &d, &r, a, b))
        goto fail;
    if (SC_ZN(d) == 0 || SC_ZN(d) - 1 < mid) {
        sc_zz_poly_mat2_move(m, &r);
        sc_value_free_many(2, c, d);
        return 1;
    }
    qr = sc_zz_poly_pseudodiv_fast(ctx, c, d);
    if (qr == NULL || !sc_zz_poly_mat2_identity(ctx, &t, a->parent))
        goto fail;
    q = qr->data.pair.first;
    rem = qr->data.pair.second;
    sc_value_ptr_swap(&t.a00, &t.a01);
    sc_value_free(t.a10);
    t.a10 = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
    if (t.a10 == NULL)
        goto fail;
    mpz_pow_ui(SC_ZP(t.a10, 0), SC_ZLC(d),
               (unsigned long)(SC_ZN(c) - SC_ZN(d) + 1));
    sc_value_free(t.a11);
    t.a11 = sc_zz_poly_neg(ctx, q);
    if (t.a11 == NULL || !sc_zz_poly_mat2_mul(ctx, &s, &t, &r) ||
        !sc_zz_poly_mat2_primitive(ctx, &s))
        goto fail;
    if (SC_ZN(rem) == 0 || SC_ZN(rem) - 1 < mid) {
        sc_zz_poly_mat2_move(m, &s);
        sc_value_free_many(3, c, d, qr);
        sc_zz_poly_mat2_clear_many(2, &r, &t);
        return 1;
    }
    k = 2 * mid - (SC_ZN(d) - 1);
    dh = sc_zz_poly_view(d, k, SC_ZN(d) - k);
    rh = sc_zz_poly_view(rem, k, SC_ZN(rem) > k ? SC_ZN(rem) - k : 0);
    if (!sc_zz_poly_hgcd_pseudo(ctx, &u, &dh, &rh) ||
        !sc_zz_poly_mat2_mul(ctx, m, &u, &s) ||
        !sc_zz_poly_mat2_primitive(ctx, m))
        goto fail;
    sc_value_free_many(3, c, d, qr);
    sc_zz_poly_mat2_clear_many(4, &r, &t, &s, &u);
    return 1;
fail:
    sc_value_free_many(3, c, d, qr);
    sc_zz_poly_mat2_clear_many(4, &r, &t, &s, &u);
    return 0;
}

sc_value *sc_zz_poly_gcd_hgcd_impl(sc_context *ctx, const sc_value *a,
                                   const sc_value *b)
{
    sc_value *ca = sc_zz_poly_content_impl(ctx, a);
    sc_value *cb = sc_zz_poly_content_impl(ctx, b), *c, *u, *v, *p, *qr, *r;
    sc_zz_poly_mat2 m = { 0 };

    c = ca && cb ? sc_zz_gcd(ctx, ca, cb) : NULL;
    u = c ? sc_zz_poly_primitive_part_impl(ctx, a) : NULL;
    v = u ? sc_zz_poly_primitive_part_impl(ctx, b) : NULL;
    if (v == NULL)
        return sc_value_free_many_null(5, ca, cb, c, u, v);
    if (SC_ZN(u) < SC_ZN(v)) {
        p = u;
        u = v;
        v = p;
    }
    while (SC_ZN(v) != 0) {
        if (SC_ZN(u) <= SC_HGCD_BASE_CUTOFF || SC_ZN(u) == SC_ZN(v)) {
            qr = sc_zz_poly_pseudodiv_fast(ctx, u, v);
            p = sc_value_pair_take(qr, 1);
            sc_value_free(u);
            u = v;
            v = p;
        } else {
            if (!sc_zz_poly_hgcd_pseudo(ctx, &m, u, v) ||
                !sc_zz_poly_mat2_apply(ctx, &p, &r, &m, u, v))
                return sc_value_free_many_null(5, ca, cb, c, u, v);
            sc_zz_poly_mat2_clear(&m);
            sc_value_free_many(2, u, v);
            u = p;
            v = r;
            if (SC_ZN(v) != 0) {
                qr = sc_zz_poly_pseudodiv_fast(ctx, u, v);
                p = sc_value_pair_take(qr, 1);
                sc_value_free(u);
                u = v;
                v = p;
            }
        }
        if (v == NULL)
            return sc_value_free_many_null(4, ca, cb, c, u);
    }
    p = sc_zz_poly_primitive_part_impl(ctx, u);
    if (p != NULL && SC_ZN(p) != 0 && mpz_sgn(SC_ZLC(p)) < 0) {
        r = sc_zz_poly_neg(ctx, p);
        sc_value_free(p);
        p = r;
    }
    r = p ? sc_zz_poly_scalar_mul_impl(ctx, p, c) : NULL;
    sc_value_free_many(5, ca, cb, c, u, p);
    sc_value_free(v);
    return r;
}

sc_value *sc_zz_poly_subres_prs_last(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t i, d, e;
    sc_zz_poly_subres_ws s;

    if (!sc_zz_poly_subres_ws_init(ctx, &s, a, b))
        return NULL;
    d = SC_ZN(s.u) - SC_ZN(s.v);
    mpz_pow_ui(s.hp, SC_ZLC(s.v), (unsigned long)d);
    s.w = sc_zz_poly_pseudorem(ctx, s.u, s.v);
    if (SC_ZN(s.w) == 0)
        return sc_zz_poly_subres_ws_finish(ctx, &s, s.v, s.hp);
    if (((d + 1) & 1) != 0)
        for (i = 0; i < SC_ZN(s.w); i++)
            mpz_neg(SC_ZP(s.w, i), SC_ZP(s.w, i));
    e = SC_ZN(s.v) - SC_ZN(s.w);
    mpz_pow_ui(s.hc, SC_ZLC(s.w), (unsigned long)e);
    mpz_pow_ui(s.den, s.hp, (unsigned long)(e - 1));
    mpz_divexact(s.hc, s.hc, s.den);
    sc_value_free_null(&s.u);
    while (1) {
        d = SC_ZN(s.v) - SC_ZN(s.w);
        s.z = sc_zz_poly_pseudorem(ctx, s.v, s.w);
        if (SC_ZN(s.z) == 0)
            break;
        if (((d + 1) & 1) != 0)
            for (i = 0; i < SC_ZN(s.z); i++)
                mpz_neg(SC_ZP(s.z, i), SC_ZP(s.z, i));
        mpz_pow_ui(s.den, s.hp, (unsigned long)d);
        mpz_mul(s.den, s.den, SC_ZLC(s.v));
        for (i = 0; i < SC_ZN(s.z); i++)
            mpz_divexact(SC_ZP(s.z, i), SC_ZP(s.z, i), s.den);
        e = SC_ZN(s.w) - SC_ZN(s.z);
        mpz_pow_ui(s.hn, SC_ZLC(s.z), (unsigned long)e);
        mpz_pow_ui(s.den, s.hc, (unsigned long)(e - 1));
        mpz_divexact(s.hn, s.hn, s.den);
        sc_zz_poly_subres_ws_shift(&s);
        mpz_swap(s.hp, s.hc);
        mpz_swap(s.hc, s.hn);
    }
    return sc_zz_poly_subres_ws_finish(ctx, &s, s.w, s.hc);
}

sc_value *sc_zz_poly_gcd_subresultant_impl(sc_context *ctx, const sc_value *a,
                                            const sc_value *b)
{
    sc_value *ca = sc_zz_poly_content_impl(ctx, a);
    sc_value *cb = sc_zz_poly_content_impl(ctx, b), *c, *u, *v, *prs, *g, *t, *r;

    c = ca && cb ? sc_zz_gcd(ctx, ca, cb) : NULL;
    u = c ? sc_zz_poly_primitive_part_impl(ctx, a) : NULL;
    v = u ? sc_zz_poly_primitive_part_impl(ctx, b) : NULL;
    if (v == NULL)
        return sc_value_free_many_null(5, ca, cb, c, u, v);
    if (u->data.zz_poly.length == 0 || v->data.zz_poly.length == 0) {
        g = u->data.zz_poly.length == 0 ? sc_value_copy_checked(ctx, v) :
                                         sc_value_copy_checked(ctx, u);
    } else {
        prs = u->data.zz_poly.length >= v->data.zz_poly.length ?
              sc_zz_poly_subres_prs_last(ctx, u, v) :
              sc_zz_poly_subres_prs_last(ctx, v, u);
        t = sc_value_pair_take(prs, 0);
        g = t ? sc_zz_poly_primitive_part_impl(ctx, t) : NULL;
        sc_value_free(t);
    }
    if (g != NULL && g->data.zz_poly.length != 0 &&
        mpz_sgn(SC_ZP(g, g->data.zz_poly.length - 1)) < 0) {
        t = sc_zz_poly_neg(ctx, g);
        sc_value_free(g);
        g = t;
    }
    r = g ? sc_zz_poly_scalar_mul_impl(ctx, g, c) : NULL;
    sc_value_free_many(6, ca, cb, c, u, v, g);
    return r;
}

sc_value *sc_zz_poly_resultant_bareiss_impl(sc_context *ctx, const sc_value *a,
                                             const sc_value *b)
{
    size_t m, n, d, i, j, k, p;
    int sign = 1;
    mpz_t *mat = NULL;
    sc_value *res = sc_value_new_zz_checked(ctx);
    mpz_t prev, t, u;

    if (res == NULL || SC_ZN(a) == 0 || SC_ZN(b) == 0)
        return res;
    m = SC_ZN(a) - 1;
    n = SC_ZN(b) - 1;
    d = m + n;
    if (d == 0) {
        mpz_set_ui(res->data.z, 1);
        return res;
    }
    mat = malloc(d * d * sizeof(mpz_t));
    if (mat == NULL) {
        sc_set_error(ctx, "out of memory computing Sylvester determinant");
        sc_value_free(res);
        return NULL;
    }
    for (i = 0; i < d * d; i++)
        mpz_init(mat[i]);
    for (i = 0; i < n; i++)
        for (j = 0; j <= m; j++)
            mpz_set(mat[i * d + i + j], SC_ZP(a, m - j));
    for (i = 0; i < m; i++)
        for (j = 0; j <= n; j++)
            mpz_set(mat[(n + i) * d + i + j], SC_ZP(b, n - j));
    mpz_inits(prev, t, u, NULL);
    mpz_set_ui(prev, 1);
    for (k = 0; k + 1 < d; k++) {
        for (p = k; p < d && mpz_sgn(mat[p * d + k]) == 0; p++)
            ;
        if (p == d) {
            mpz_set_ui(res->data.z, 0);
            goto done;
        }
        if (p != k) {
            for (j = k; j < d; j++)
                mpz_swap(mat[k * d + j], mat[p * d + j]);
            sign = -sign;
        }
        for (i = k + 1; i < d; i++)
            for (j = k + 1; j < d; j++) {
                mpz_mul(t, mat[i * d + j], mat[k * d + k]);
                mpz_mul(u, mat[i * d + k], mat[k * d + j]);
                mpz_sub(t, t, u);
                if (k != 0)
                    mpz_divexact(t, t, prev);
                mpz_set(mat[i * d + j], t);
            }
        mpz_set(prev, mat[k * d + k]);
    }
    mpz_set(res->data.z, mat[(d - 1) * d + d - 1]);
    if (sign < 0)
        mpz_neg(res->data.z, res->data.z);
done:
    mpz_clears(prev, t, u, NULL);
    for (i = 0; i < d * d; i++)
        mpz_clear(mat[i]);
    free(mat);
    return res;
}

sc_value *sc_zz_poly_resultant_subresultant_impl(sc_context *ctx,
                                                   const sc_value *a,
                                                   const sc_value *b)
{
    size_t m, n;
    int swap;
    sc_value *ca, *cb, *u, *v, *prs, *last, *res;
    mpz_t t;

    res = sc_value_new_zz_checked(ctx);
    if (res == NULL || a->data.zz_poly.length == 0 || b->data.zz_poly.length == 0)
        return res;
    ca = sc_zz_poly_content_impl(ctx, a);
    cb = sc_zz_poly_content_impl(ctx, b);
    u = ca ? sc_zz_poly_primitive_part_impl(ctx, a) : NULL;
    v = cb ? sc_zz_poly_primitive_part_impl(ctx, b) : NULL;
    if (v == NULL)
        return sc_value_free_many_null(5, res, ca, cb, u, v);
    m = u->data.zz_poly.length - 1;
    n = v->data.zz_poly.length - 1;
    swap = m < n;
    prs = swap ? sc_zz_poly_subres_prs_last(ctx, v, u) :
                 sc_zz_poly_subres_prs_last(ctx, u, v);
    if (prs == NULL)
        return sc_value_free_many_null(5, res, ca, cb, u, v);
    sc_value_free(res);
    res = sc_value_copy_checked(ctx, prs->data.pair.second);
    last = sc_value_pair_take(prs, 0);
    if (res == NULL || last == NULL)
        return sc_value_free_many_null(6, res, ca, cb, u, v, last);
    if (last->data.zz_poly.length > 1)
        mpz_set_ui(res->data.z, 0);
    else if (swap && ((m * n) & 1) != 0)
        mpz_neg(res->data.z, res->data.z);
    mpz_init(t);
    mpz_pow_ui(t, ca->data.z, (unsigned long)n);
    mpz_mul(res->data.z, res->data.z, t);
    mpz_pow_ui(t, cb->data.z, (unsigned long)m);
    mpz_mul(res->data.z, res->data.z, t);
    mpz_clear(t);
    sc_value_free_many(5, ca, cb, u, v, last);
    return res;
}

sc_value *sc_zz_poly_xgcd_subresultant_impl(sc_context *ctx, const sc_value *a,
                                             const sc_value *b)
{
    size_t i, j, d, e;
    int first = 1, neg;
    sc_zz_poly_xgcd_ws s;
    sc_value *qr, *q, *t, *p;

    if (!sc_zz_poly_xgcd_ws_init(ctx, &s, a, b))
        return NULL;
    d = SC_ZN(s.u) - SC_ZN(s.v);
    mpz_pow_ui(s.hp, SC_ZLC(s.v), (unsigned long)d);
    while (1) {
        d = SC_ZN(s.u) - SC_ZN(s.v);
        qr = sc_zz_poly_pseudodiv(ctx, s.u, s.v);
        if (!sc_value_pair_split(qr, &q, &s.w))
            return sc_zz_poly_xgcd_ws_abort(&s);
        if (SC_ZN(s.w) == 0) {
            sc_value_free_many(2, q, s.w);
            s.w = NULL;
            return sc_zz_poly_xgcd_ws_finish(ctx, &s);
        }
        mpz_pow_ui(s.alpha, SC_ZLC(s.v), (unsigned long)(d + 1));
        if (!first) {
            mpz_pow_ui(s.den, s.hp, (unsigned long)d);
            mpz_mul(s.den, s.den, SC_ZLC(s.u));
        }
        neg = ((d + 1) & 1) != 0;
        for (i = 0; i < SC_ZN(s.w); i++) {
            if (neg)
                mpz_neg(SC_ZP(s.w, i), SC_ZP(s.w, i));
            if (!first)
                mpz_divexact(SC_ZP(s.w, i), SC_ZP(s.w, i), s.den);
        }
        for (j = 0; j < 2; j++) {
            s.cw[j] = sc_value_copy_checked(ctx, s.cu[j]);
            if (s.cw[j] == NULL) {
                sc_value_free(q);
                return sc_zz_poly_xgcd_ws_abort(&s);
            }
            for (i = 0; i < SC_ZN(s.cw[j]); i++)
                mpz_mul(SC_ZP(s.cw[j], i), SC_ZP(s.cw[j], i), s.alpha);
            t = sc_zz_poly_mul(ctx, q, s.cv[j]);
            p = t ? sc_zz_poly_sub(ctx, s.cw[j], t) : NULL;
            sc_value_free_many(2, s.cw[j], t);
            s.cw[j] = p;
            if (p == NULL) {
                sc_value_free(q);
                return sc_zz_poly_xgcd_ws_abort(&s);
            }
            for (i = 0; i < SC_ZN(p); i++) {
                if (neg)
                    mpz_neg(SC_ZP(p, i), SC_ZP(p, i));
                if (!first)
                    mpz_divexact(SC_ZP(p, i), SC_ZP(p, i), s.den);
            }
        }
        sc_value_free(q);
        e = SC_ZN(s.v) - SC_ZN(s.w);
        mpz_pow_ui(first ? s.hc : s.hn, SC_ZLC(s.w), (unsigned long)e);
        mpz_pow_ui(s.den, first ? s.hp : s.hc, (unsigned long)(e - 1));
        mpz_divexact(first ? s.hc : s.hn, first ? s.hc : s.hn, s.den);
        if (!first) {
            mpz_swap(s.hp, s.hc);
            mpz_swap(s.hc, s.hn);
        }
        first = 0;
        sc_zz_poly_xgcd_ws_shift(&s);
    }
}


sc_value *sc_zz_poly_evaluate_horner_impl(sc_context *ctx, const sc_value *a,
                                           const sc_value *b)
{
    size_t i = SC_ZN(a);
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL || i == 0)
        return r;

    mpz_set(r->data.z, SC_ZP(a, --i));

    while (i != 0) {
        mpz_mul(r->data.z, r->data.z, b->data.z);
        mpz_add(r->data.z, r->data.z, SC_ZP(a, --i));
    }

    return r;
}

sc_value *sc_zz_poly_evaluate_divconquer_impl(sc_context *ctx, const sc_value *a,
                                               const sc_value *b)
{
    size_t i, m, n = SC_ZN(a);
    sc_value *r = sc_value_new_zz_checked(ctx), *w;
    mpz_t power;

    if (r == NULL || n == 0)
        return r;

    w = sc_value_copy_checked(ctx, a);
    if (w == NULL)
        return sc_value_free_many_null(1, r);

    mpz_init_set(power, b->data.z);
    while (n > 1) {
        m = n / 2;
        mpz_addmul(SC_ZP(w, 0), power, SC_ZP(w, 1));
        for (i = 1; i < m; i++) {
            mpz_mul(SC_ZP(w, i), power, SC_ZP(w, 2 * i + 1));
            mpz_add(SC_ZP(w, i), SC_ZP(w, i), SC_ZP(w, 2 * i));
        }
        if (n & 1)
            mpz_set(SC_ZP(w, m), SC_ZP(w, n - 1));
        n = (n + 1) / 2;
        if (n > 1)
            mpz_mul(power, power, power);
    }

    mpz_set(r->data.z, SC_ZP(w, 0));
    mpz_clear(power);
    sc_value_free(w);
    return r;
}

sc_value *sc_zz_poly_compose_horner_impl(sc_context *ctx, const sc_value *a,
                                         const sc_value *b)
{
    size_t i = SC_ZN(a);
    sc_value *r, *t;

    if (i == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);

    r = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
    if (r == NULL)
        return NULL;
    mpz_set(SC_ZP(r, 0), SC_ZP(a, --i));

    while (i != 0) {
        t = sc_zz_poly_mul(ctx, r, b);
        sc_value_free(r);
        if (t == NULL)
            return NULL;
        r = t;
        i--;
        if (SC_ZN(r) == 0 && mpz_sgn(SC_ZP(a, i)) != 0) {
            sc_value_free(r);
            r = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
            if (r == NULL)
                return NULL;
        }
        if (SC_ZN(r) != 0)
            mpz_add(SC_ZP(r, 0), SC_ZP(r, 0), SC_ZP(a, i));
        sc_zz_poly_normalize(r);
    }

    return r;
}

sc_value *sc_zz_poly_compose_divconquer_impl(sc_context *ctx, const sc_value *a,
                                             const sc_value *b)
{
    size_t i, m, n = SC_ZN(a), count = n;
    sc_value **w;
    sc_value *power, *t, *u, *r;

    if (n <= 1)
        return sc_value_copy_checked(ctx, a);

    w = sc_value_array_new_checked(ctx, n);
    if (w == NULL)
        return NULL;
    for (i = 0; i < n; i++) {
        size_t len = mpz_sgn(SC_ZP(a, i)) != 0;

        w[i] = sc_value_new_zz_poly_checked(ctx, a->parent, len);
        if (w[i] == NULL) {
            sc_value_array_free(count, w);
            return NULL;
        }
        if (len != 0)
            mpz_set(SC_ZP(w[i], 0), SC_ZP(a, i));
    }

    power = sc_value_copy_checked(ctx, b);
    if (power == NULL) {
        sc_value_array_free(count, w);
        return NULL;
    }
    while (n > 1) {
        m = n / 2;
        for (i = 0; i < m; i++) {
            t = sc_zz_poly_mul(ctx, power, w[2 * i + 1]);
            u = t == NULL ? NULL : sc_zz_poly_add(ctx, w[2 * i], t);
            sc_value_free(t);
            sc_value_free_null(&w[2 * i]);
            sc_value_free_null(&w[2 * i + 1]);
            if (u == NULL) {
                sc_value_free(power);
                sc_value_array_free(count, w);
                return NULL;
            }
            w[i] = u;
        }
        if (n & 1) {
            w[m] = w[n - 1];
            w[n - 1] = NULL;
        }
        n = (n + 1) / 2;
        if (n > 1) {
            t = sc_zz_poly_mul(ctx, power, power);
            sc_value_free(power);
            power = t;
            if (power == NULL) {
                sc_value_array_free(count, w);
                return NULL;
            }
        }
    }

    r = w[0];
    w[0] = NULL;
    sc_value_free(power);
    sc_value_array_free(count, w);
    return r;
}

sc_value *sc_zz_poly_taylor_shift_horner_impl(sc_context *ctx, const sc_value *a,
                                              const sc_value *b)
{
    size_t n = SC_ZN(a), i, j, d = 0;
    sc_value *r;

    if (n == 0 || mpz_sgn(b->data.z) == 0)
        return sc_value_copy_checked(ctx, a);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (r == NULL)
        return NULL;
    mpz_set(SC_ZP(r, 0), SC_ZP(a, n - 1));
    for (i = n - 1; i-- != 0; d++) {
        mpz_set(SC_ZP(r, d + 1), SC_ZP(r, d));
        for (j = d; j != 0; j--) {
            mpz_mul(SC_ZP(r, j), SC_ZP(r, j), b->data.z);
            mpz_add(SC_ZP(r, j), SC_ZP(r, j), SC_ZP(r, j - 1));
        }
        mpz_mul(SC_ZP(r, 0), SC_ZP(r, 0), b->data.z);
        mpz_add(SC_ZP(r, 0), SC_ZP(r, 0), SC_ZP(a, i));
    }
    return r;
}

sc_value *sc_zz_poly_taylor_shift_divconquer_impl(sc_context *ctx,
                                                  const sc_value *a,
                                                  const sc_value *b)
{
    sc_value *g, *r;

    if (SC_ZN(a) == 0 || mpz_sgn(b->data.z) == 0)
        return sc_value_copy_checked(ctx, a);
    g = sc_value_new_zz_poly_checked(ctx, a->parent, 2);
    if (g == NULL)
        return NULL;
    mpz_set(SC_ZP(g, 0), b->data.z);
    mpz_set_ui(SC_ZP(g, 1), 1);
    r = sc_zz_poly_compose_divconquer_impl(ctx, a, g);
    sc_value_free(g);
    return r;
}


sc_value *sc_zz_poly_derivative_impl(sc_context *ctx, const sc_value *a)
{
    size_t i, n = SC_ZN(a);
    sc_value *r = sc_value_new_zz_poly_checked(ctx, a->parent, n > 1 ? n - 1 : 0);

    if (r == NULL)
        return NULL;
    for (i = 1; i < n; i++)
        mpz_mul_ui(SC_ZP(r, i - 1), SC_ZP(a, i), (unsigned long)i);
    return r;
}

sc_value *sc_zz_poly_nth_derivative_impl(sc_context *ctx, const sc_value *a, size_t k)
{
    size_t i, n = SC_ZN(a);
    sc_value *r;
    mpz_t fall;

    if (k == 0)
        return sc_value_copy_checked(ctx, a);
    if (k >= n)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n - k);
    if (r == NULL)
        return NULL;
    mpz_init(fall);
    mpz_fac_ui(fall, (unsigned long)k);
    for (i = k; i < n; i++) {
        mpz_mul(SC_ZP(r, i - k), SC_ZP(a, i), fall);
        if (i + 1 < n) {
            mpz_mul_ui(fall, fall, (unsigned long)(i + 1));
            mpz_divexact_ui(fall, fall, (unsigned long)(i + 1 - k));
        }
    }
    mpz_clear(fall);
    return r;
}

sc_value *sc_zz_poly_discriminant_impl(sc_context *ctx, const sc_value *a)
{
    size_t n;
    sc_value *d, *r;

    if (SC_ZN(a) <= 1)
        return sc_value_new_zz_checked(ctx);
    n = SC_ZN(a) - 1;
    d = sc_zz_poly_derivative_impl(ctx, a);
    r = d ? sc_zz_poly_resultant(ctx, a, d) : NULL;
    sc_value_free(d);
    if (r == NULL)
        return NULL;
    mpz_divexact(r->data.z, r->data.z, SC_ZLC(a));
    if ((n & 3) == 2 || (n & 3) == 3)
        mpz_neg(r->data.z, r->data.z);
    return r;
}

sc_value *sc_zz_poly_is_squarefree_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx), *d, *g;

    if (r == NULL || SC_ZN(a) <= 2) {
        if (r != NULL)
            mpz_set_ui(r->data.z, SC_ZN(a) != 0);
        return r;
    }
    d = sc_zz_poly_derivative_impl(ctx, a);
    g = d ? sc_zz_poly_gcd(ctx, a, d) : NULL;
    sc_value_free(d);
    if (g == NULL) {
        sc_value_free(r);
        return NULL;
    }
    mpz_set_ui(r->data.z, SC_ZN(g) <= 1);
    sc_value_free(g);
    return r;
}

sc_value *sc_zz_poly_squarefree_part_impl(sc_context *ctx, const sc_value *a)
{
    sc_value *d, *g, *q, *r;

    if (SC_ZN(a) == 0)
        return sc_value_copy_checked(ctx, a);
    d = sc_zz_poly_derivative_impl(ctx, a);
    g = d ? sc_zz_poly_gcd(ctx, a, d) : NULL;
    q = g ? sc_zz_poly_divexact(ctx, a, g) : NULL;
    sc_value_free_many(2, d, g);
    if (q == NULL || SC_ZN(q) == 0 || mpz_sgn(SC_ZLC(q)) > 0)
        return q;
    r = sc_zz_poly_neg(ctx, q);
    sc_value_free(q);
    return r;
}
