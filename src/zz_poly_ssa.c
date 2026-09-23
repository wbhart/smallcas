#include "smallcas.h"
#include "smallcas_fft.h"

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

static int sc_ssa_params(size_t *k, unsigned *depth, unsigned *logn,
                         const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, total, need, n;
    size_t ba = sc_zz_poly_max_abs_bits_raw(a);
    size_t bb = sc_zz_poly_max_abs_bits_raw(b);

    if ((GMP_NUMB_BITS & (GMP_NUMB_BITS - 1)) != 0 || an > SIZE_MAX - bn + 1)
        return 0;
    total = an + bn - 1;
    *logn = (unsigned)sc_ssa_log2ceil(total);
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

sc_value *sc_zz_poly_mul_ssa(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, outn, k, limbs, words;
    unsigned depth, logn;
    sc_fft_mod m;
    sc_fft_plan p;
    mp_ptr av = NULL, bv, work = NULL;
    sc_value *r = NULL;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (!sc_ssa_params(&k, &depth, &logn, a, b) ||
        !sc_fft_mod_init(&m, 1, (mp_bitcnt_t)k, 2, depth)) {
        sc_set_error(ctx, "SSA parameter setup failed");
        return NULL;
    }
    if (!sc_fft_plan_init(&p, logn, &m))
        goto fail_mod;
    limbs = (size_t)m.n;
    if (p.len > SIZE_MAX / limbs / 2)
        goto fail_plan;
    words = 2 * p.len * limbs;
    av = calloc(words, sizeof(mp_limb_t));
    work = calloc(4 * limbs + 1, sizeof(mp_limb_t));
    outn = an + bn - 1;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, outn);
    if (av == NULL || work == NULL || r == NULL)
        goto fail;
    bv = av + p.len * limbs;
    for (size_t i = 0; i < an; i++)
        sc_ssa_import(sc_fft_entry(av, i, &m), a->data.zz_poly.coeff[i], &m);
    for (size_t i = 0; i < bn; i++)
        sc_ssa_import(sc_fft_entry(bv, i, &m), b->data.zz_poly.coeff[i], &m);
    sc_ssa_forward(av, &p, &m, work);
    sc_ssa_forward(bv, &p, &m, work);
    for (size_t i = 0; i < p.len; i++)
        sc_fft_mul(sc_fft_entry(av, i, &m), sc_fft_entry(av, i, &m),
                   sc_fft_entry(bv, i, &m), &m, work);
    sc_ssa_inverse(av, &p, &m, work);
    for (size_t i = 0; i < outn; i++)
        sc_ssa_export(r->data.zz_poly.coeff[i], sc_fft_entry_const(av, i, &m),
                      &m, work);
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
