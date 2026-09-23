#include "smallcas.h"

typedef struct sc_zz_poly_qpoly {
    sc_value *poly;
    mpq_t scale;
    int initialized;
} sc_zz_poly_qpoly;

typedef struct sc_zz_poly_qmat2 {
    sc_zz_poly_qpoly e[4];
} sc_zz_poly_qmat2;

#define QP(v, i) ((v)->poly->data.zz_poly.coeff[i])
#define QN(v) ((v)->poly->data.zz_poly.length)
#define QLC(v) QP((v), QN(v) - 1)

static void qpoly_init(sc_zz_poly_qpoly *a)
{
    a->poly = NULL;
    mpq_init(a->scale);
    mpq_set_ui(a->scale, 1, 1);
    a->initialized = 1;
}

static void qpoly_clear(sc_zz_poly_qpoly *a)
{
    if (!a->initialized)
        return;
    sc_value_free(a->poly);
    mpq_clear(a->scale);
    a->poly = NULL;
    a->initialized = 0;
}

static void qmat_clear(sc_zz_poly_qmat2 *a)
{
    size_t i;

    for (i = 0; i < 4; i++)
        qpoly_clear(&a->e[i]);
}

static int qpoly_normalize(sc_context *ctx, sc_zz_poly_qpoly *a)
{
    size_t i;
    mpz_t c;

    if (a->poly == NULL)
        return 0;
    if (QN(a) == 0) {
        mpq_set_ui(a->scale, 1, 1);
        return 1;
    }
    mpz_init_set_ui(c, 0);
    for (i = 0; i < QN(a) && mpz_cmp_ui(c, 1) != 0; i++)
        mpz_gcd(c, c, QP(a, i));
    if (mpz_cmp_ui(c, 1) > 0) {
        for (i = 0; i < QN(a); i++)
            mpz_divexact(QP(a, i), QP(a, i), c);
        mpz_mul(mpq_numref(a->scale), mpq_numref(a->scale), c);
        mpq_canonicalize(a->scale);
    }
    mpz_clear(c);
    (void)ctx;
    return 1;
}

static int qpoly_set(sc_context *ctx, sc_zz_poly_qpoly *r, const sc_value *a)
{
    qpoly_init(r);
    r->poly = sc_value_copy_checked(ctx, a);
    return r->poly != NULL && qpoly_normalize(ctx, r);
}

static int qpoly_copy(sc_context *ctx, sc_zz_poly_qpoly *r,
                      const sc_zz_poly_qpoly *a)
{
    qpoly_init(r);
    r->poly = sc_value_copy_checked(ctx, a->poly);
    if (r->poly == NULL)
        return 0;
    mpq_set(r->scale, a->scale);
    return 1;
}

static int qpoly_set_si(sc_context *ctx, sc_zz_poly_qpoly *r,
                        sc_parent *parent, long x)
{
    qpoly_init(r);
    r->poly = sc_value_new_zz_poly_checked(ctx, parent, x == 0 ? 0 : 1);
    if (r->poly == NULL)
        return 0;
    if (x != 0)
        mpz_set_ui(QP(r, 0), 1);
    mpq_set_si(r->scale, x == 0 ? 1 : x, 1);
    return 1;
}

static int qpoly_add(sc_context *ctx, sc_zz_poly_qpoly *r,
                     const sc_zz_poly_qpoly *a, const sc_zz_poly_qpoly *b)
{
    size_t i, n = QN(a) > QN(b) ? QN(a) : QN(b);
    mpz_t den, ma, mb, t;

    qpoly_init(r);
    r->poly = sc_value_new_zz_poly_checked(ctx, a->poly->parent, n);
    if (r->poly == NULL)
        return 0;
    mpz_inits(den, ma, mb, t, NULL);
    mpz_lcm(den, mpq_denref(a->scale), mpq_denref(b->scale));
    mpz_divexact(ma, den, mpq_denref(a->scale));
    mpz_mul(ma, ma, mpq_numref(a->scale));
    mpz_divexact(mb, den, mpq_denref(b->scale));
    mpz_mul(mb, mb, mpq_numref(b->scale));
    for (i = 0; i < n; i++) {
        if (i < QN(a))
            mpz_mul(QP(r, i), ma, QP(a, i));
        if (i < QN(b)) {
            mpz_mul(t, mb, QP(b, i));
            mpz_add(QP(r, i), QP(r, i), t);
        }
    }
    mpq_set_z(r->scale, den);
    mpq_inv(r->scale, r->scale);
    mpz_clears(den, ma, mb, t, NULL);
    sc_zz_poly_normalize(r->poly);
    return qpoly_normalize(ctx, r);
}

static int qpoly_neg(sc_context *ctx, sc_zz_poly_qpoly *r,
                     const sc_zz_poly_qpoly *a)
{
    if (!qpoly_copy(ctx, r, a))
        return 0;
    mpq_neg(r->scale, r->scale);
    return 1;
}

static int qpoly_mul(sc_context *ctx, sc_zz_poly_qpoly *r,
                     const sc_zz_poly_qpoly *a, const sc_zz_poly_qpoly *b)
{
    qpoly_init(r);
    r->poly = sc_zz_poly_mul(ctx, a->poly, b->poly);
    if (r->poly == NULL)
        return 0;
    mpq_mul(r->scale, a->scale, b->scale);
    return qpoly_normalize(ctx, r);
}

static int qpoly_high(sc_context *ctx, sc_zz_poly_qpoly *r,
                      const sc_zz_poly_qpoly *a, size_t start)
{
    size_t i, n = QN(a) > start ? QN(a) - start : 0;

    qpoly_init(r);
    r->poly = sc_value_new_zz_poly_checked(ctx, a->poly->parent, n);
    if (r->poly == NULL)
        return 0;
    for (i = 0; i < n; i++)
        mpz_set(QP(r, i), QP(a, start + i));
    mpq_set(r->scale, a->scale);
    sc_zz_poly_normalize(r->poly);
    return qpoly_normalize(ctx, r);
}

static int qpoly_divrem(sc_context *ctx, sc_zz_poly_qpoly *q,
                        sc_zz_poly_qpoly *r, const sc_zz_poly_qpoly *a,
                        const sc_zz_poly_qpoly *b)
{
    size_t delta;
    sc_value *pair = NULL, *qp = NULL, *rp = NULL;
    mpz_t lc;
    mpq_t t;

    if (QN(b) == 0) {
        sc_set_error(ctx, "LR division by zero polynomial");
        return 0;
    }
    if (QN(a) < QN(b)) {
        if (!qpoly_set_si(ctx, q, a->poly->parent, 0) || !qpoly_copy(ctx, r, a))
            goto fail;
        return 1;
    }
    pair = sc_zz_poly_pseudodiv_fast(ctx, a->poly, b->poly);
    if (pair == NULL || !sc_value_pair_split(pair, &qp, &rp))
        goto fail;
    pair = NULL;
    delta = QN(a) - QN(b) + 1;
    mpz_init(lc);
    mpq_init(t);
    mpz_pow_ui(lc, QLC(b), (unsigned long)delta);
    qpoly_init(q);
    q->poly = qp;
    qp = NULL;
    mpq_div(q->scale, a->scale, b->scale);
    mpq_set_z(t, lc);
    mpq_div(q->scale, q->scale, t);
    qpoly_init(r);
    r->poly = rp;
    rp = NULL;
    mpq_div(r->scale, a->scale, t);
    mpq_clear(t);
    mpz_clear(lc);
    if (!qpoly_normalize(ctx, q) || !qpoly_normalize(ctx, r))
        goto fail;
    return 1;
fail:
    sc_value_free_many(3, pair, qp, rp);
    qpoly_clear(q);
    qpoly_clear(r);
    return 0;
}

static int qmat_identity(sc_context *ctx, sc_zz_poly_qmat2 *m, sc_parent *parent)
{
    return qpoly_set_si(ctx, &m->e[0], parent, 1) &&
           qpoly_set_si(ctx, &m->e[1], parent, 0) &&
           qpoly_set_si(ctx, &m->e[2], parent, 0) &&
           qpoly_set_si(ctx, &m->e[3], parent, 1);
}

static int qmat_step(sc_context *ctx, sc_zz_poly_qmat2 *m,
                     const sc_zz_poly_qpoly *q)
{
    return qpoly_set_si(ctx, &m->e[0], q->poly->parent, 0) &&
           qpoly_set_si(ctx, &m->e[1], q->poly->parent, 1) &&
           qpoly_set_si(ctx, &m->e[2], q->poly->parent, 1) &&
           qpoly_neg(ctx, &m->e[3], q);
}

static int qmat_mul(sc_context *ctx, sc_zz_poly_qmat2 *r,
                    const sc_zz_poly_qmat2 *a, const sc_zz_poly_qmat2 *b)
{
    size_t i, j;
    sc_zz_poly_qmat2 t = { 0 };
    sc_zz_poly_qpoly p = { 0 }, q = { 0 };

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            if (!qpoly_mul(ctx, &p, &a->e[2 * i], &b->e[j]) ||
                !qpoly_mul(ctx, &q, &a->e[2 * i + 1], &b->e[2 + j]) ||
                !qpoly_add(ctx, &t.e[2 * i + j], &p, &q))
                goto fail;
            qpoly_clear(&p);
            qpoly_clear(&q);
        }
    }
    *r = t;
    return 1;
fail:
    qpoly_clear(&p);
    qpoly_clear(&q);
    qmat_clear(&t);
    return 0;
}

static int qmat_apply(sc_context *ctx, sc_zz_poly_qpoly *u,
                      sc_zz_poly_qpoly *v, const sc_zz_poly_qmat2 *m,
                      const sc_zz_poly_qpoly *a, const sc_zz_poly_qpoly *b)
{
    sc_zz_poly_qpoly p = { 0 }, q = { 0 }, s = { 0 }, t = { 0 };

    if (!qpoly_mul(ctx, &p, &m->e[0], a) || !qpoly_mul(ctx, &q, &m->e[1], b) ||
        !qpoly_add(ctx, &s, &p, &q))
        goto fail;
    qpoly_clear(&p);
    qpoly_clear(&q);
    if (!qpoly_mul(ctx, &p, &m->e[2], a) || !qpoly_mul(ctx, &q, &m->e[3], b) ||
        !qpoly_add(ctx, &t, &p, &q))
        goto fail;
    qpoly_clear(&p);
    qpoly_clear(&q);
    *u = s;
    *v = t;
    return 1;
fail:
    qpoly_clear(&p);
    qpoly_clear(&q);
    qpoly_clear(&s);
    qpoly_clear(&t);
    return 0;
}

static int qhgcd(sc_context *ctx, sc_zz_poly_qmat2 *m,
                 const sc_zz_poly_qpoly *a, const sc_zz_poly_qpoly *b)
{
    size_t mid, shift;
    sc_zz_poly_qpoly ah = { 0 }, bh = { 0 }, c = { 0 }, d = { 0 };
    sc_zz_poly_qpoly q = { 0 }, rem = { 0 }, dh = { 0 }, rh = { 0 };
    sc_zz_poly_qmat2 r = { 0 }, step = { 0 }, s = { 0 }, t = { 0 };

    if (QN(a) == 0 || QN(b) == 0)
        return qmat_identity(ctx, m, a->poly->parent);
    mid = QN(a) / 2;
    if (QN(b) <= mid)
        return qmat_identity(ctx, m, a->poly->parent);
    if (QN(a) == QN(b)) {
        if (!qpoly_divrem(ctx, &q, &rem, a, b) || !qmat_step(ctx, &step, &q))
            goto fail;
        if (QN(&rem) <= mid) {
            *m = step;
            step = (sc_zz_poly_qmat2){ 0 };
            goto done;
        }
        if (!qhgcd(ctx, &t, b, &rem) || !qmat_mul(ctx, m, &t, &step))
            goto fail;
        goto done;
    }
    if (!qpoly_high(ctx, &ah, a, mid) || !qpoly_high(ctx, &bh, b, mid) ||
        !qhgcd(ctx, &r, &ah, &bh) || !qmat_apply(ctx, &c, &d, &r, a, b))
        goto fail;
    if (QN(&d) <= mid) {
        *m = r;
        r = (sc_zz_poly_qmat2){ 0 };
        goto done;
    }
    if (!qpoly_divrem(ctx, &q, &rem, &c, &d) || !qmat_step(ctx, &step, &q) ||
        !qmat_mul(ctx, &s, &step, &r))
        goto fail;
    if (QN(&rem) <= mid) {
        *m = s;
        s = (sc_zz_poly_qmat2){ 0 };
        goto done;
    }
    shift = 2 * mid - (QN(&d) - 1);
    if (!qpoly_high(ctx, &dh, &d, shift) || !qpoly_high(ctx, &rh, &rem, shift) ||
        !qhgcd(ctx, &t, &dh, &rh) || !qmat_mul(ctx, m, &t, &s))
        goto fail;
done:
    qpoly_clear(&ah);
    qpoly_clear(&bh);
    qpoly_clear(&c);
    qpoly_clear(&d);
    qpoly_clear(&q);
    qpoly_clear(&rem);
    qpoly_clear(&dh);
    qpoly_clear(&rh);
    qmat_clear(&r);
    qmat_clear(&step);
    qmat_clear(&s);
    qmat_clear(&t);
    return 1;
fail:
    qmat_clear(m);
    goto done_fail;
done_fail:
    qpoly_clear(&ah);
    qpoly_clear(&bh);
    qpoly_clear(&c);
    qpoly_clear(&d);
    qpoly_clear(&q);
    qpoly_clear(&rem);
    qpoly_clear(&dh);
    qpoly_clear(&rh);
    qmat_clear(&r);
    qmat_clear(&step);
    qmat_clear(&s);
    qmat_clear(&t);
    return 0;
}

sc_value *sc_zz_poly_gcd_lr_impl(sc_context *ctx, const sc_value *a,
                                 const sc_value *b)
{
    sc_value *ca = NULL, *cb = NULL, *cont = NULL, *pa = NULL, *pb = NULL;
    sc_value *g = NULL, *out = NULL, *tmp;
    sc_zz_poly_qpoly u = { 0 }, v = { 0 }, x = { 0 }, y = { 0 };
    sc_zz_poly_qpoly q = { 0 }, rem = { 0 };
    sc_zz_poly_qmat2 m = { 0 };

    ca = sc_zz_poly_content_impl(ctx, a);
    cb = sc_zz_poly_content_impl(ctx, b);
    cont = ca && cb ? sc_zz_gcd(ctx, ca, cb) : NULL;
    pa = cont ? sc_zz_poly_primitive_part_impl(ctx, a) : NULL;
    pb = pa ? sc_zz_poly_primitive_part_impl(ctx, b) : NULL;
    if (pb == NULL)
        goto fail;
    if (pa->data.zz_poly.length < pb->data.zz_poly.length) {
        tmp = pa;
        pa = pb;
        pb = tmp;
    }
    if (!qpoly_set(ctx, &u, pa) || !qpoly_set(ctx, &v, pb))
        goto fail;
    while (QN(&v) != 0) {
        if (QN(&u) <= 12 || QN(&u) == QN(&v)) {
            if (!qpoly_divrem(ctx, &q, &rem, &u, &v))
                goto fail;
            qpoly_clear(&u);
            u = v;
            v = rem;
            v.initialized = 1;
            rem = (sc_zz_poly_qpoly){ 0 };
            qpoly_clear(&q);
            continue;
        }
        if (!qhgcd(ctx, &m, &u, &v) || !qmat_apply(ctx, &x, &y, &m, &u, &v))
            goto fail;
        qmat_clear(&m);
        if (QN(&y) >= QN(&v)) {
            qpoly_clear(&x);
            qpoly_clear(&y);
            if (!qpoly_divrem(ctx, &q, &rem, &u, &v))
                goto fail;
            qpoly_clear(&u);
            u = v;
            v = rem;
            v.initialized = 1;
            rem = (sc_zz_poly_qpoly){ 0 };
            qpoly_clear(&q);
        } else {
            qpoly_clear(&u);
            qpoly_clear(&v);
            u = x;
            v = y;
            x = (sc_zz_poly_qpoly){ 0 };
            y = (sc_zz_poly_qpoly){ 0 };
        }
    }
    g = sc_value_copy_checked(ctx, u.poly);
    if (g == NULL)
        goto fail;
    if (g->data.zz_poly.length != 0 &&
        mpz_sgn(g->data.zz_poly.coeff[g->data.zz_poly.length - 1]) < 0) {
        tmp = sc_zz_poly_neg(ctx, g);
        sc_value_free(g);
        g = tmp;
    }
    out = g ? sc_zz_poly_scalar_mul_impl(ctx, g, cont) : NULL;
fail:
    sc_value_free_many(6, ca, cb, cont, pa, pb, g);
    qpoly_clear(&u);
    qpoly_clear(&v);
    qpoly_clear(&x);
    qpoly_clear(&y);
    qpoly_clear(&q);
    qpoly_clear(&rem);
    qmat_clear(&m);
    return out;
}
