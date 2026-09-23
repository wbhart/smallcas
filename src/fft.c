#include "smallcas_fft.h"
#include "tuning.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

int sc_fft_plan_init(sc_fft_plan *p, unsigned logn, const sc_fft_mod *m)
{
    mp_ptr s, w, wi, ms;
    size_t half;

    memset(p, 0, sizeof(*p));
    if (logn >= sizeof(size_t) * CHAR_BIT ||
        m->depth >= sizeof(size_t) * CHAR_BIT || logn > m->depth)
        return 0;
    p->len = (size_t)1 << logn;
    p->logn = logn;
    p->root_stride = (size_t)1 << (m->depth - logn);
    if (logn == 0 || m->fermat)
        return 1;
    half = p->len >> 1;
    p->fwd = calloc(2 * half * (size_t)m->n, sizeof(mp_limb_t));
    s = calloc((size_t)5 * m->n + 1, sizeof(mp_limb_t));
    if (p->fwd == NULL || s == NULL)
        goto fail;
    p->inv = p->fwd + half * (size_t)m->n;
    w = s;
    wi = w + m->n;
    ms = wi + m->n;
    mpn_copyi(w, m->root, m->n);
    mpn_copyi(wi, m->root_inv, m->n);
    for (unsigned i = logn; i < m->depth; i++) {
        sc_fft_mul(w, w, w, m, ms);
        sc_fft_mul(wi, wi, wi, m, ms);
    }
    sc_fft_set_ui(p->fwd, 1, m);
    sc_fft_set_ui(p->inv, 1, m);
    for (size_t i = 1; i < half; i++) {
        sc_fft_mul(sc_fft_entry(p->fwd, i, m),
                   sc_fft_entry(p->fwd, i - 1, m), w, m, ms);
        sc_fft_mul(sc_fft_entry(p->inv, i, m),
                   sc_fft_entry(p->inv, i - 1, m), wi, m, ms);
    }
    free(s);
    return 1;
fail:
    free(s);
    sc_fft_plan_clear(p);
    return 0;
}

void sc_fft_plan_clear(sc_fft_plan *p)
{
    free(p->fwd);
    memset(p, 0, sizeof(*p));
}

static void sc_fft_root_pow(mp_ptr a, size_t e, int inv,
                            const sc_fft_plan *p, const sc_fft_mod *m,
                            mp_ptr scratch)
{
    size_t half;
    int neg = 0;

    if (p->len == 1)
        return;
    e &= p->len - 1;
    if (m->fermat) {
        size_t k = e * p->root_stride;
        if (inv && k != 0)
            k = ((size_t)1 << m->depth) - k;
        sc_fft_mul_2exp(a, a, k, m, scratch);
        return;
    }
    half = p->len >> 1;
    if (e >= half)
        neg = 1, e -= half;
    if (e != 0) {
        mp_srcptr w = sc_fft_entry_const(inv ? p->inv : p->fwd, e, m);
        sc_fft_mul(a, a, w, m, scratch);
    }
    if (neg)
        sc_fft_neg(a, a, m);
}

static void sc_fft_twiddle(mp_ptr a, size_t e, int inv,
                           const sc_fft_plan *p, const sc_fft_mod *m,
                           mp_ptr scratch)
{
    if (m->fermat) {
        size_t k = e * p->root_stride;

        if (inv && k != 0)
            k = ((size_t)1 << m->depth) - k;
        sc_fft_mul_2exp(a, a, k, m, scratch);
    } else {
        mp_srcptr w = sc_fft_entry_const(inv ? p->inv : p->fwd, e, m);

        sc_fft_mul(a, a, w, m, scratch);
    }
}

void sc_fft_forward(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                    mp_ptr scratch)
{
    mp_ptr t = scratch, ms = scratch + m->n;
    size_t stride = 1;

    for (size_t len = p->len; len > 1; len >>= 1, stride <<= 1) {
        size_t h = len >> 1;
        for (size_t b = 0; b < p->len; b += len)
            for (size_t j = 0; j < h; j++) {
                mp_ptr x = sc_fft_entry(a, b + j, m);
                mp_ptr y = sc_fft_entry(a, b + j + h, m);
                sc_fft_addsub(x, y, t, m);
                sc_fft_twiddle(y, j * stride, 0, p, m, ms);
            }
    }
}

void sc_fft_inverse(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                    mp_ptr scratch)
{
    mp_ptr t = scratch, ms = scratch + m->n;
    size_t stride = p->len >> 1;

    for (size_t len = 2; len <= p->len && len != 0; len <<= 1, stride >>= 1) {
        size_t h = len >> 1;
        for (size_t b = 0; b < p->len; b += len)
            for (size_t j = 0; j < h; j++) {
                mp_ptr x = sc_fft_entry(a, b + j, m);
                mp_ptr y = sc_fft_entry(a, b + j + h, m);
                sc_fft_twiddle(y, j * stride, 1, p, m, ms);
                sc_fft_addsub(x, y, t, m);
            }
        if (len == p->len)
            break;
    }
    for (size_t i = 0; i < p->len; i++)
        sc_fft_div_2exp(sc_fft_entry(a, i, m), p->logn, m);
}

static void sc_fft_forward_base(mp_ptr a, unsigned logn, size_t stride,
                                size_t root_stride, size_t zeta,
                                const sc_fft_plan *p, const sc_fft_mod *m,
                                mp_ptr scratch)
{
    mp_ptr t = scratch, ms = scratch + m->n;
    size_t len = (size_t)1 << logn, h = len >> 1;

    if (logn == 0)
        return;
    for (size_t j = 0; j < h; j++) {
        mp_ptr x = sc_fft_entry(a, j * stride, m);
        mp_ptr y = sc_fft_entry(a, (j + h) * stride, m);

        sc_fft_addsub(x, y, t, m);
        sc_fft_root_pow(y, zeta + j * root_stride, 0, p, m, ms);
    }
    if (logn > 1) {
        zeta <<= 1;
        sc_fft_forward_base(a, logn - 1, stride, root_stride << 1,
                            zeta, p, m, scratch);
        sc_fft_forward_base(sc_fft_entry(a, h * stride, m), logn - 1,
                            stride, root_stride << 1, zeta, p, m, scratch);
    }
}

static void sc_fft_inverse_base(mp_ptr a, unsigned logn, size_t stride,
                                size_t root_stride, size_t zeta,
                                const sc_fft_plan *p, const sc_fft_mod *m,
                                mp_ptr scratch)
{
    mp_ptr t = scratch, ms = scratch + m->n;
    size_t len = (size_t)1 << logn, h = len >> 1;

    if (logn == 0)
        return;
    if (logn > 1) {
        size_t z2 = zeta << 1;

        sc_fft_inverse_base(a, logn - 1, stride, root_stride << 1,
                            z2, p, m, scratch);
        sc_fft_inverse_base(sc_fft_entry(a, h * stride, m), logn - 1,
                            stride, root_stride << 1, z2, p, m, scratch);
    }
    for (size_t j = 0; j < h; j++) {
        mp_ptr x = sc_fft_entry(a, j * stride, m);
        mp_ptr y = sc_fft_entry(a, (j + h) * stride, m);

        sc_fft_root_pow(y, zeta + j * root_stride, 1, p, m, ms);
        sc_fft_addsub(x, y, t, m);
    }
}

static void sc_fft_forward_mfa_rec(mp_ptr a, unsigned logn, size_t stride,
                                   size_t root_stride, size_t zeta,
                                   const sc_fft_plan *p, const sc_fft_mod *m,
                                   mp_ptr scratch)
{
    unsigned log1, log2;
    size_t n1, n2;

    if (logn <= SC_FFT_MFA_BASE_LOG) {
        sc_fft_forward_base(a, logn, stride, root_stride, zeta, p, m, scratch);
        return;
    }
    log1 = logn >> 1;
    log2 = logn - log1;
    n1 = (size_t)1 << log1;
    n2 = (size_t)1 << log2;
    for (size_t u = 0; u < n2; u++)
        sc_fft_forward_mfa_rec(sc_fft_entry(a, u * stride, m), log1,
                               stride * n2, root_stride * n2,
                               zeta + u * root_stride, p, m, scratch);
    zeta *= n1;
    for (size_t u = 0; u < n1; u++)
        sc_fft_forward_mfa_rec(sc_fft_entry(a, u * n2 * stride, m), log2,
                               stride, root_stride * n1, zeta, p, m, scratch);
}

static void sc_fft_inverse_mfa_rec(mp_ptr a, unsigned logn, size_t stride,
                                   size_t root_stride, size_t zeta,
                                   const sc_fft_plan *p, const sc_fft_mod *m,
                                   mp_ptr scratch)
{
    unsigned log1, log2;
    size_t n1, n2, row_zeta;

    if (logn <= SC_FFT_MFA_BASE_LOG) {
        sc_fft_inverse_base(a, logn, stride, root_stride, zeta, p, m, scratch);
        return;
    }
    log1 = logn >> 1;
    log2 = logn - log1;
    n1 = (size_t)1 << log1;
    n2 = (size_t)1 << log2;
    row_zeta = zeta * n1;
    for (size_t u = 0; u < n1; u++)
        sc_fft_inverse_mfa_rec(sc_fft_entry(a, u * n2 * stride, m), log2,
                               stride, root_stride * n1, row_zeta,
                               p, m, scratch);
    for (size_t u = 0; u < n2; u++)
        sc_fft_inverse_mfa_rec(sc_fft_entry(a, u * stride, m), log1,
                               stride * n2, root_stride * n2,
                               zeta + u * root_stride, p, m, scratch);
}

void sc_fft_forward_mfa(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                        mp_ptr scratch)
{
    sc_fft_forward_mfa_rec(a, p->logn, 1, 1, 0, p, m, scratch);
}

void sc_fft_inverse_mfa(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                        mp_ptr scratch)
{
    sc_fft_inverse_mfa_rec(a, p->logn, 1, 1, 0, p, m, scratch);
    for (size_t i = 0; i < p->len; i++)
        sc_fft_div_2exp(sc_fft_entry(a, i, m), p->logn, m);
}

static void sc_fft_forward_tft_rec(mp_ptr a, unsigned logn, size_t stride,
                                   size_t root_stride, size_t zeta,
                                   size_t z, size_t n,
                                   const sc_fft_plan *p, const sc_fft_mod *m,
                                   mp_ptr scratch)
{
    size_t len = (size_t)1 << logn;
    unsigned log1, log2;
    size_t n1, n2, n1p, z1, z2, z2p, l1, l2;

    if (z == len && n == len) {
        sc_fft_forward_mfa_rec(a, logn, stride, root_stride, zeta,
                               p, m, scratch);
        return;
    }
    if (logn == 0)
        return;
    if (logn == 1) {
        mp_ptr x = a, y = sc_fft_entry(a, stride, m);

        if (n == 2 && z == 2) {
            sc_fft_addsub(x, y, scratch, m);
            sc_fft_root_pow(y, zeta, 0, p, m, scratch + m->n);
        } else if (n == 2) {
            mpn_copyi(y, x, m->n);
            sc_fft_root_pow(y, zeta, 0, p, m, scratch);
        } else if (z == 2) {
            sc_fft_add(x, x, y, m);
        }
        return;
    }
    log1 = logn >> 1;
    log2 = logn - log1;
    l1 = (size_t)1 << log1;
    l2 = (size_t)1 << log2;
    n1 = n / l2;
    n2 = n % l2;
    n1p = n1 + (n2 != 0);
    z1 = z / l2;
    z2 = z % l2;
    z2p = z1 != 0 ? l2 : z2;

    for (size_t u = 0; u < z2; u++)
        sc_fft_forward_tft_rec(sc_fft_entry(a, u * stride, m), log1,
                               stride * l2, root_stride * l2,
                               zeta + u * root_stride, z1 + 1, n1p,
                               p, m, scratch);
    for (size_t u = z2; u < z2p; u++)
        sc_fft_forward_tft_rec(sc_fft_entry(a, u * stride, m), log1,
                               stride * l2, root_stride * l2,
                               zeta + u * root_stride, z1, n1p,
                               p, m, scratch);
    zeta *= l1;
    for (size_t u = 0; u < n1; u++)
        sc_fft_forward_tft_rec(sc_fft_entry(a, u * l2 * stride, m), log2,
                               stride, root_stride * l1, zeta, z2p, l2,
                               p, m, scratch);
    if (n2 != 0)
        sc_fft_forward_tft_rec(sc_fft_entry(a, n1 * l2 * stride, m), log2,
                               stride, root_stride * l1, zeta, z2p, n2,
                               p, m, scratch);
}

static void sc_fft_itft_rec(mp_ptr a, unsigned logn, size_t stride,
                            size_t root_stride, size_t zeta,
                            size_t z, size_t n, int f,
                            const sc_fft_plan *p, const sc_fft_mod *m,
                            mp_ptr scratch)
{
    size_t len = (size_t)1 << logn;
    unsigned log1, log2;
    size_t n1, n2, z1, z2, z2p, l1, l2, lo, hi;

    if (n == len) {
        sc_fft_inverse_mfa_rec(a, logn, stride, root_stride, zeta,
                               p, m, scratch);
        return;
    }
    if (logn == 0)
        return;
    if (logn == 1) {
        mp_ptr x = a, y = sc_fft_entry(a, stride, m), t = scratch;

        if (n == 1) {
            mpn_copyi(t, x, m->n);
            if (z == 2)
                sc_fft_sub(x, x, y, m);
            if (f) {
                if (z == 2)
                    sc_fft_sub(y, t, y, m);
                else
                    mpn_copyi(y, t, m->n);
                sc_fft_root_pow(y, zeta, 0, p, m, t + m->n);
            }
            sc_fft_add(x, x, t, m);
        } else {
            if (z == 2)
                sc_fft_add(x, x, y, m);
            sc_fft_div_2exp(x, 1, m);
        }
        return;
    }
    log1 = logn >> 1;
    log2 = logn - log1;
    l1 = (size_t)1 << log1;
    l2 = (size_t)1 << log2;
    n1 = n / l2;
    n2 = n % l2;
    z1 = z / l2;
    z2 = z % l2;
    z2p = z1 != 0 ? l2 : z2;
    lo = n2 < z2 ? n2 : z2;
    hi = n2 > z2 ? n2 : z2;

    for (size_t u = 0; u < n1; u++)
        sc_fft_itft_rec(sc_fft_entry(a, u * l2 * stride, m), log2,
                        stride, root_stride * l1, zeta * l1,
                        l2, l2, 0, p, m, scratch);
    for (size_t u = n2; u < hi; u++)
        sc_fft_itft_rec(sc_fft_entry(a, u * stride, m), log1,
                        stride * l2, root_stride * l2,
                        zeta + u * root_stride, z1 + 1, n1,
                        n2 + (size_t)f != 0, p, m, scratch);
    for (size_t u = hi; u < z2p; u++)
        sc_fft_itft_rec(sc_fft_entry(a, u * stride, m), log1,
                        stride * l2, root_stride * l2,
                        zeta + u * root_stride, z1, n1,
                        n2 + (size_t)f != 0, p, m, scratch);
    if (n2 + (size_t)f != 0)
        sc_fft_itft_rec(sc_fft_entry(a, n1 * l2 * stride, m), log2,
                        stride, root_stride * l1, zeta * l1,
                        z2p, n2, f, p, m, scratch);
    for (size_t u = 0; u < lo; u++)
        sc_fft_itft_rec(sc_fft_entry(a, u * stride, m), log1,
                        stride * l2, root_stride * l2,
                        zeta + u * root_stride, z1 + 1, n1 + 1, 0,
                        p, m, scratch);
    for (size_t u = lo; u < n2; u++)
        sc_fft_itft_rec(sc_fft_entry(a, u * stride, m), log1,
                        stride * l2, root_stride * l2,
                        zeta + u * root_stride, z1, n1 + 1, 0,
                        p, m, scratch);
}

void sc_fft_forward_tft(mp_ptr a, size_t z, size_t n, const sc_fft_plan *p,
                        const sc_fft_mod *m, mp_ptr scratch)
{
    if (z == 0 || n == 0 || z > p->len || n > p->len)
        return;
    sc_fft_forward_tft_rec(a, p->logn, 1, 1, 0, z, n, p, m, scratch);
}

void sc_fft_inverse_tft(mp_ptr a, size_t n, const sc_fft_plan *p,
                        const sc_fft_mod *m, mp_ptr scratch)
{
    if (n == 0 || n > p->len)
        return;
    sc_fft_itft_rec(a, p->logn, 1, 1, 0, n, n, 0, p, m, scratch);
    for (size_t i = 0; i < n; i++)
        sc_fft_div_2exp(sc_fft_entry(a, i, m), p->logn, m);
}
