#include "smallcas_fft.h"

#include <stdlib.h>
#include <string.h>

static void sc_fft_set_u64(mpz_t z, uint64_t x)
{
    mpz_import(z, 1, 1, sizeof(x), 0, 0, &x);
}

static void sc_fft_export(mp_ptr r, mp_size_t n, mpz_srcptr z)
{
    size_t count = 0;
    mpn_zero(r, n);
    mpz_export(r, &count, -1, sizeof(mp_limb_t), 0, GMP_NAIL_BITS, z);
}

int sc_fft_mod_init(sc_fft_mod *m, uint64_t c, mp_bitcnt_t bits,
                    uint64_t root, unsigned depth)
{
    mpz_t p, w, t;
    int ok = 0;

    memset(m, 0, sizeof(*m));
    mpz_inits(p, w, t, NULL);
    if (c == 0 || bits == 0 || depth == 0)
        goto done;
    sc_fft_set_u64(p, c);
    mpz_mul_2exp(p, p, bits);
    mpz_add_ui(p, p, 1);
    sc_fft_set_u64(w, root);
    mpz_mod(w, w, p);
    mpz_set(t, w);
    for (unsigned i = 1; i < depth; i++)
        mpz_mul(t, t, t), mpz_mod(t, t, p);
    mpz_add_ui(t, t, 1);
    if (mpz_cmp(t, p) != 0)
        goto done;
    m->n = mpz_size(p);
    m->mod = calloc((size_t)3 * m->n, sizeof(mp_limb_t));
    if (m->mod == NULL)
        goto done;
    m->root = m->mod + m->n;
    m->root_inv = m->root + m->n;
    sc_fft_export(m->mod, m->n, p);
    sc_fft_export(m->root, m->n, w);
    if (mpz_invert(t, w, p) == 0)
        goto done;
    sc_fft_export(m->root_inv, m->n, t);
    m->c = c;
    m->bits = bits;
    m->depth = depth;
    m->fermat = c == 1 && root == 2 && bits % GMP_NUMB_BITS == 0 &&
                bits <= (mp_bitcnt_t)((size_t)-1 / 2);
    ok = 1;
done:
    mpz_clears(p, w, t, NULL);
    if (!ok)
        sc_fft_mod_clear(m);
    return ok;
}

void sc_fft_mod_clear(sc_fft_mod *m)
{
    free(m->mod);
    memset(m, 0, sizeof(*m));
}

static void sc_fft_reduce_fermat(mp_ptr r, mp_srcptr lo, mp_srcptr hi,
                                 const sc_fft_mod *m)
{
    mp_size_t l = m->n - 1;
    mp_limb_t cy = mpn_sub_n(r, lo, hi, l);

    r[l] = 0;
    if (cy)
        r[l] = mpn_add_1(r, r, l, 1);
}

void sc_fft_mul(mp_ptr r, mp_srcptr a, mp_srcptr b,
                const sc_fft_mod *m, mp_ptr scratch)
{
    mp_size_t n = m->n, pn = 2 * n;
    mp_ptr prod = scratch, q = scratch + 2 * n;

    if (m->fermat) {
        mp_size_t l = n - 1;
        if (a[l]) {
            if (b[l])
                sc_fft_set_ui(r, 1, m);
            else
                sc_fft_neg(r, b, m);
        } else if (b[l]) {
            sc_fft_neg(r, a, m);
        } else {
            mpn_mul_n(prod, a, b, l);
            sc_fft_reduce_fermat(r, prod, prod + l, m);
        }
        return;
    }
    mpn_mul_n(prod, a, b, n);
    while (pn > 0 && prod[pn - 1] == 0)
        pn--;
    if (pn < n || (pn == n && mpn_cmp(prod, m->mod, n) < 0)) {
        mpn_zero(r, n);
        if (pn != 0)
            mpn_copyi(r, prod, pn);
    } else {
        mpn_tdiv_qr(q, r, 0, prod, pn, m->mod, n);
    }
}

void sc_fft_sqr(mp_ptr r, mp_srcptr a, const sc_fft_mod *m, mp_ptr scratch)
{
    mp_size_t n = m->n, pn = 2 * n;
    mp_ptr prod = scratch, q = scratch + 2 * n;

    if (m->fermat) {
        mp_size_t l = n - 1;
        if (a[l])
            sc_fft_set_ui(r, 1, m);
        else {
            mpn_sqr(prod, a, l);
            sc_fft_reduce_fermat(r, prod, prod + l, m);
        }
        return;
    }
    mpn_sqr(prod, a, n);
    while (pn > 0 && prod[pn - 1] == 0)
        pn--;
    if (pn < n || (pn == n && mpn_cmp(prod, m->mod, n) < 0)) {
        mpn_zero(r, n);
        if (pn != 0)
            mpn_copyi(r, prod, pn);
    } else {
        mpn_tdiv_qr(q, r, 0, prod, pn, m->mod, n);
    }
}

void sc_fft_mul_2exp(mp_ptr r, mp_srcptr a, size_t e,
                     const sc_fft_mod *m, mp_ptr scratch)
{
    size_t bits = (size_t)m->bits, ord = 2 * bits;
    mp_size_t l = m->n - 1;
    int neg;

    e %= ord;
    neg = e >= bits;
    if (neg)
        e -= bits;
    if (e == 0) {
        if (r != a)
            mpn_copyi(r, a, m->n);
    } else {
        size_t q = e / GMP_NUMB_BITS, sh = e % GMP_NUMB_BITS;
        mp_ptr t = scratch;
        mp_limb_t cy = 0;
        mpn_zero(t, 2 * l + 1);
        if (sh != 0)
            cy = mpn_lshift(t + q, a, m->n, sh);
        else
            mpn_copyi(t + q, a, m->n);
        t[q + m->n] = cy;
        sc_fft_reduce_fermat(r, t, t + l, m);
    }
    if (neg)
        sc_fft_neg(r, r, m);
}

void sc_fft_div_2exp(mp_ptr a, size_t k, const sc_fft_mod *m)
{
    while (k-- != 0) {
        mp_limb_t cy = 0;
        if (a[0] & 1)
            cy = mpn_add_n(a, a, m->mod, m->n);
        mpn_rshift(a, a, m->n, 1);
        if (cy)
            a[m->n - 1] |= (mp_limb_t)1 << (GMP_NUMB_BITS - 1);
    }
}
