#include "smallcas.h"
#include "smallcas_fft.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

static size_t sc_ssa_log2ceil(size_t n)
{
    size_t k = 0;

    if (n != 0)
        n--;
    while (n != 0)
        k++, n >>= 1;
    return k;
}

static int sc_ssa_params_len(size_t *k, unsigned *depth, unsigned *logn,
                             const sc_value *a, const sc_value *b, size_t len)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, need, n;
    size_t ba = sc_zz_poly_max_abs_bits_raw(a);
    size_t bb = sc_zz_poly_max_abs_bits_raw(b);

    if ((GMP_NUMB_BITS & (GMP_NUMB_BITS - 1)) != 0 || len == 0)
        return 0;
    *logn = (unsigned)sc_ssa_log2ceil(len);
    if (*logn >= sizeof(size_t) * 8)
        return 0;
    n = (size_t)1 << *logn;
    if (ba > SIZE_MAX - bb || ba + bb > SIZE_MAX - sc_ssa_log2ceil(small) - 1)
        return 0;
    need = ba + bb + sc_ssa_log2ceil(small) + 1;
    if (need < (n >> 1))
        need = n >> 1;
    *k = GMP_NUMB_BITS;
    while (*k < need) {
        if (*k > SIZE_MAX / 2)
            return 0;
        *k <<= 1;
    }
    n = 2 * *k;
    *depth = 0;
    while (n > 1)
        (*depth)++, n >>= 1;
    return 1;
}

static void sc_ssa_import(mp_ptr r, mpz_srcptr z, const sc_fft_mod *m)
{
    size_t count = 0;

    mpn_zero(r, m->n);
    mpz_export(r, &count, -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS, z);
    if (mpz_sgn(z) < 0)
        sc_fft_neg(r, r, m);
}

static int sc_ssa_negative(mp_srcptr a, const sc_fft_mod *m)
{
    mp_size_t l = m->n - 1;
    mp_limb_t half = (mp_limb_t)1 << (GMP_NUMB_BITS - 1);

    if (a[l] != 0 || a[l - 1] > half)
        return 1;
    if (a[l - 1] < half)
        return 0;
    for (mp_size_t i = 0; i + 1 < l; i++)
        if (a[i] != 0)
            return 1;
    return 0;
}

static void sc_ssa_export(mpz_ptr z, mp_srcptr a, const sc_fft_mod *m, mp_ptr t)
{
    int neg = sc_ssa_negative(a, m);

    if (neg) {
        sc_fft_neg(t, a, m);
        a = t;
    }
    mpz_import(z, (size_t)m->n, -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS, a);
    if (neg)
        mpz_neg(z, z);
}

static void sc_ssa_forward(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                           mp_ptr scratch)
{
    if (p->logn >= SC_SSA_MFA_CUTOFF_LOG)
        sc_fft_forward_mfa(a, p, m, scratch);
    else
        sc_fft_forward(a, p, m, scratch);
}

static void sc_ssa_inverse(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                           mp_ptr scratch)
{
    if (p->logn >= SC_SSA_MFA_CUTOFF_LOG)
        sc_fft_inverse_mfa(a, p, m, scratch);
    else
        sc_fft_inverse(a, p, m, scratch);
}

static sc_value *sc_zz_poly_mul_ssa_core(sc_context *ctx, const sc_value *a,
                                           const sc_value *b, size_t len,
                                           size_t start, size_t outn, int cyclic,
                                           int square)
{
    size_t an = a->data.zz_poly.length;
    size_t bn = square ? an : b->data.zz_poly.length;
    size_t k, limbs, words, pointn;
    unsigned depth, logn;
    sc_fft_mod m;
    sc_fft_plan p;
    mp_ptr av = NULL, bv, work = NULL;
    sc_value *r = NULL;

    if (!sc_ssa_params_len(&k, &depth, &logn, a, square ? a : b, len) ||
        !sc_fft_mod_init(&m, 1, (mp_bitcnt_t)k, 2, depth)) {
        sc_set_error(ctx, "SSA parameter setup failed");
        return NULL;
    }
    if (!sc_fft_plan_init(&p, logn, &m))
        goto fail_mod;
    limbs = (size_t)m.n;
    if (p.len > SIZE_MAX / limbs / (square ? 1 : 2))
        goto fail_plan;
    words = (square ? 1 : 2) * p.len * limbs;
    av = calloc(words, sizeof(mp_limb_t));
    work = calloc(4 * limbs + 1, sizeof(mp_limb_t));
    r = sc_value_new_zz_poly_checked(ctx, a->parent, outn);
    if (av == NULL || work == NULL || r == NULL)
        goto fail;
    bv = square ? av : av + p.len * limbs;
    for (size_t i = 0; i < an; i++)
        sc_ssa_import(sc_fft_entry(av, i, &m), a->data.zz_poly.coeff[i], &m);
    if (!square)
        for (size_t i = 0; i < bn; i++)
            sc_ssa_import(sc_fft_entry(bv, i, &m), b->data.zz_poly.coeff[i], &m);
    if (!cyclic && outn < p.len) {
        sc_fft_forward_tft(av, an, outn, &p, &m, work);
        if (!square)
            sc_fft_forward_tft(bv, bn, outn, &p, &m, work);
    } else {
        sc_ssa_forward(av, &p, &m, work);
        if (!square)
            sc_ssa_forward(bv, &p, &m, work);
    }
    pointn = cyclic ? p.len : outn;
    for (size_t i = 0; i < pointn; i++)
        if (square)
            sc_fft_sqr(sc_fft_entry(av, i, &m), sc_fft_entry(av, i, &m), &m, work);
        else
            sc_fft_mul(sc_fft_entry(av, i, &m), sc_fft_entry(av, i, &m),
                       sc_fft_entry(bv, i, &m), &m, work);
    if (!cyclic && outn < p.len)
        sc_fft_inverse_tft(av, outn, &p, &m, work);
    else
        sc_ssa_inverse(av, &p, &m, work);
    for (size_t i = 0; i < outn; i++)
        sc_ssa_export(r->data.zz_poly.coeff[i],
                      sc_fft_entry_const(av, start + i, &m), &m, work);
    free(work);
    free(av);
    sc_fft_plan_clear(&p);
    sc_fft_mod_clear(&m);
    sc_zz_poly_normalize(r);
    return r;
fail:
    sc_value_free(r);
    free(work);
    free(av);
fail_plan:
    sc_fft_plan_clear(&p);
fail_mod:
    sc_fft_mod_clear(&m);
    if (ctx->error[0] == '\0')
        sc_set_error(ctx, "out of memory in SSA multiplication");
    return NULL;
}

sc_value *sc_zz_poly_mulmid_ssa(sc_context *ctx, const sc_value *a,
                                const sc_value *b, size_t n)
{
    size_t an, bn;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (n == 0 || an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (n > SIZE_MAX / 2 + 1 || an > 2 * n - 1 || bn > n) {
        sc_set_error(ctx, "SSA middle-product shape unsupported");
        return NULL;
    }
    return sc_zz_poly_mul_ssa_core(ctx, a, b, 2 * n - 1, n - 1, n, 1, 0);
}

sc_value *sc_zz_poly_mul_ssa(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, outn;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (an > SIZE_MAX - bn + 1) {
        sc_set_error(ctx, "SSA parameter setup failed");
        return NULL;
    }
    outn = an + bn - 1;
    return sc_zz_poly_mul_ssa_core(ctx, a, b, outn, 0, outn, 0, 0);
}

sc_value *sc_zz_poly_sqr_ssa(sc_context *ctx, const sc_value *a)
{
    size_t n, outn;

    if (a == NULL)
        return NULL;
    n = a->data.zz_poly.length;
    if (n == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (n > SIZE_MAX / 2 + 1) {
        sc_set_error(ctx, "SSA square parameter setup failed");
        return NULL;
    }
    outn = 2 * n - 1;
    return sc_zz_poly_mul_ssa_core(ctx, a, a, outn, 0, outn, 0, 1);
}

static int sc_ssa_taylor_params(size_t *k, size_t *need_out, unsigned *depth,
                                unsigned *logn, const sc_value *a,
                                const sc_value *c)
{
    size_t n = a->data.zz_poly.length, need, span, N, bs, bq;
    mpz_t s, q, t;

    if (n == 0 || n - 1 > ULONG_MAX ||
        (GMP_NUMB_BITS & (GMP_NUMB_BITS - 1)) != 0 ||
        n > SIZE_MAX / 2 + 1)
        return 0;
    span = 2 * n - 1;
    *logn = (unsigned)sc_ssa_log2ceil(span);
    if (*logn >= sizeof(size_t) * 8)
        return 0;
    N = (size_t)1 << *logn;
    mpz_inits(s, q, t, NULL);
    mpz_set_ui(s, 0);
    for (size_t i = 0; i < n; i++) {
        mpz_abs(t, a->data.zz_poly.coeff[i]);
        mpz_add(s, s, t);
    }
    mpz_abs(t, c->data.z);
    mpz_add_ui(t, t, 1);
    mpz_pow_ui(q, t, (unsigned long)(n - 1));
    bs = mpz_sgn(s) == 0 ? 0 : mpz_sizeinbase(s, 2);
    bq = mpz_sizeinbase(q, 2);
    mpz_clears(s, q, t, NULL);
    if (bs > SIZE_MAX - bq - 1)
        return 0;
    need = bs + bq + 1;
    if (need < (N >> 1))
        need = N >> 1;
    if (need_out != NULL)
        *need_out = need;
    *k = GMP_NUMB_BITS;
    while (*k < need) {
        if (*k > SIZE_MAX / 2)
            return 0;
        *k <<= 1;
    }
    N = 2 * *k;
    *depth = 0;
    while (N > 1)
        (*depth)++, N >>= 1;
    return 1;
}

static int sc_ssa_taylor_fill(mp_ptr av, mp_ptr bv, mp_ptr rt, mp_ptr work,
                              size_t n, const sc_value *a, const sc_value *c,
                              const sc_fft_mod *m, mpz_srcptr mod, mpz_ptr invtop)
{
    mpz_t fact, invfact, cpow, t;

    mpz_inits(fact, invfact, cpow, t, NULL);
    mpz_set_ui(fact, 1);
    for (size_t i = 0; i < n; i++) {
        mpz_mod(t, a->data.zz_poly.coeff[i], mod);
        mpz_mul(t, t, fact);
        mpz_mod(t, t, mod);
        sc_ssa_import(sc_fft_entry(av, n - 1 - i, m), t, m);
        if (i + 1 < n)
            mpz_mul_ui(fact, fact, (unsigned long)(i + 1)), mpz_mod(fact, fact, mod);
    }
    if (mpz_invert(invtop, fact, mod) == 0) {
        mpz_clears(fact, invfact, cpow, t, NULL);
        return 0;
    }
    mpz_set_ui(cpow, 1);
    for (size_t i = 0; i < n; i++) {
        sc_ssa_import(sc_fft_entry(bv, i, m), cpow, m);
        if (i + 1 < n)
            mpz_mul(cpow, cpow, c->data.z), mpz_mod(cpow, cpow, mod);
    }
    mpz_set(invfact, invtop);
    for (size_t i = n; i-- != 0;) {
        sc_ssa_import(rt, invfact, m);
        sc_fft_mul(sc_fft_entry(bv, i, m), sc_fft_entry(bv, i, m), rt, m, work);
        if (i != 0)
            mpz_mul_ui(invfact, invfact, (unsigned long)i), mpz_mod(invfact, invfact, mod);
    }
    mpz_clears(fact, invfact, cpow, t, NULL);
    return 1;
}

static void sc_ssa_taylor_extract(sc_value *r, mp_ptr av, mp_ptr rt, mp_ptr work,
                                  size_t n, const sc_fft_mod *m, mpz_srcptr mod,
                                  mpz_srcptr invtop)
{
    mpz_t invfact;

    mpz_init_set(invfact, invtop);
    for (size_t j = n; j-- != 0;) {
        sc_ssa_import(rt, invfact, m);
        sc_fft_mul(sc_fft_entry(av, n - 1 - j, m),
                   sc_fft_entry(av, n - 1 - j, m), rt, m, work);
        sc_ssa_export(r->data.zz_poly.coeff[j],
                      sc_fft_entry_const(av, n - 1 - j, m), m, work);
        if (j != 0)
            mpz_mul_ui(invfact, invfact, (unsigned long)j), mpz_mod(invfact, invfact, mod);
    }
    mpz_clear(invfact);
}

int sc_zz_poly_taylor_shift_convolution_preferred(const sc_value *a,
                                                     const sc_value *c)
{
    size_t n = a->data.zz_poly.length, k, need, p, lo, hi;
    unsigned depth, logn, logp;

    if (n < SC_TAYLOR_CONV_CUTOFF || n < 16)
        return 0;
    logp = sc_ssa_log2ceil(n);
    if (logp >= sizeof(size_t) * 8)
        return 0;
    p = (size_t)1 << logp;
    lo = 5 * (p >> 3);
    hi = 13 * (p >> 4);
    if (n < lo || n > hi)
        return 0;
    if (!sc_ssa_taylor_params(&k, &need, &depth, &logn, a, c))
        return 0;
    (void)depth;
    (void)logn;
    return need >= 5 * (k >> 3);
}

sc_value *sc_zz_poly_taylor_shift_convolution_impl(sc_context *ctx,
                                                   const sc_value *a,
                                                   const sc_value *c)
{
    size_t n = a->data.zz_poly.length, k, limbs, words;
    unsigned depth, logn;
    sc_fft_mod m;
    sc_fft_plan p;
    mp_ptr av = NULL, bv, work = NULL, rt;
    sc_value *r = NULL;
    mpz_t mod, invtop;

    if (n == 0 || mpz_sgn(c->data.z) == 0)
        return sc_value_copy_checked(ctx, a);
    if (!sc_ssa_taylor_params(&k, NULL, &depth, &logn, a, c) ||
        !sc_fft_mod_init(&m, 1, (mp_bitcnt_t)k, 2, depth)) {
        sc_set_error(ctx, "SSA Taylor-shift parameter setup failed");
        return NULL;
    }
    if (!sc_fft_plan_init(&p, logn, &m))
        goto fail_mod;
    limbs = (size_t)m.n;
    if (p.len > SIZE_MAX / limbs / 2 || limbs > (SIZE_MAX - 1) / 5)
        goto fail_plan;
    words = 2 * p.len * limbs;
    av = calloc(words, sizeof(mp_limb_t));
    work = calloc(5 * limbs + 1, sizeof(mp_limb_t));
    r = sc_value_new_zz_poly_checked(ctx, a->parent, n);
    if (av == NULL || work == NULL || r == NULL)
        goto fail;
    bv = av + p.len * limbs;
    rt = work + 4 * limbs + 1;
    mpz_inits(mod, invtop, NULL);
    mpz_set_ui(mod, 1);
    mpz_mul_2exp(mod, mod, (mp_bitcnt_t)k);
    mpz_add_ui(mod, mod, 1);
    if (!sc_ssa_taylor_fill(av, bv, rt, work, n, a, c, &m, mod, invtop)) {
        sc_set_error(ctx, "Taylor-shift factorial is not invertible");
        goto fail_mpz;
    }
    sc_ssa_forward(av, &p, &m, work);
    sc_ssa_forward(bv, &p, &m, work);
    for (size_t i = 0; i < p.len; i++)
        sc_fft_mul(sc_fft_entry(av, i, &m), sc_fft_entry(av, i, &m),
                   sc_fft_entry(bv, i, &m), &m, work);
    sc_ssa_inverse(av, &p, &m, work);
    sc_ssa_taylor_extract(r, av, rt, work, n, &m, mod, invtop);
    mpz_clears(mod, invtop, NULL);
    free(work);
    free(av);
    sc_fft_plan_clear(&p);
    sc_fft_mod_clear(&m);
    sc_zz_poly_normalize(r);
    return r;
fail_mpz:
    mpz_clears(mod, invtop, NULL);
fail:
    sc_value_free(r);
    free(work);
    free(av);
fail_plan:
    sc_fft_plan_clear(&p);
fail_mod:
    sc_fft_mod_clear(&m);
    if (ctx->error[0] == '\0')
        sc_set_error(ctx, "out of memory in SSA Taylor shift");
    return NULL;
}
