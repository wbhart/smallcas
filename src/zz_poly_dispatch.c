#include "smallcas.h"

#include <stdint.h>

static size_t sc_zz_poly_valuation(const sc_value *a)
{
    size_t i = 0;

    while (i < a->data.zz_poly.length && mpz_sgn(a->data.zz_poly.coeff[i]) == 0)
        i++;
    return i;
}

static int sc_zz_poly_zero_prefix(const sc_value *a, size_t n)
{
    size_t i;

    for (i = 0; i < n && i < a->data.zz_poly.length; i++)
        if (mpz_sgn(a->data.zz_poly.coeff[i]) != 0)
            return 0;
    return 1;
}

static int sc_zz_poly_unit_lead(const sc_value *a)
{
    size_t n = a->data.zz_poly.length;

    return n != 0 && mpz_cmpabs_ui(a->data.zz_poly.coeff[n - 1], 1) == 0;
}

static int sc_zz_poly_window_equal(const sc_value *a, size_t start,
                                   const sc_value *b, size_t n)
{
    size_t i;

    for (i = 0; i < n; i++) {
        int az = start + i >= a->data.zz_poly.length;
        int bz = i >= b->data.zz_poly.length;

        if (az != bz && (az ? mpz_sgn(b->data.zz_poly.coeff[i]) :
                         mpz_sgn(a->data.zz_poly.coeff[start + i])) != 0)
            return 0;
        if (!az && !bz && mpz_cmp(a->data.zz_poly.coeff[start + i],
                                  b->data.zz_poly.coeff[i]) != 0)
            return 0;
    }
    return 1;
}

static size_t sc_ceil_log2_size(size_t n)
{
    size_t bits = 0;

    if (n != 0)
        n--;
    while (n != 0) {
        bits++;
        n >>= 1;
    }
    return bits;
}

static mp_bitcnt_t sc_zz_poly_ks_bits(const sc_value *a, const sc_value *b)
{
    size_t small = a->data.zz_poly.length < b->data.zz_poly.length ?
                   a->data.zz_poly.length : b->data.zz_poly.length;
    size_t bits_a = sc_zz_poly_max_abs_bits_raw(a);
    size_t bits_b = sc_zz_poly_max_abs_bits_raw(b);

    return (mp_bitcnt_t)(bits_a + bits_b + sc_ceil_log2_size(small) + 1);
}

static int sc_zz_poly_use_sqr_ks(const sc_value *a)
{
    size_t n = a->data.zz_poly.length;
    size_t bits = sc_zz_poly_max_abs_bits_raw(a);

    return n >= SC_SQR_KS_CUTOFF && bits < n;
}

static int sc_zz_poly_use_sqr_toom3(const sc_value *a)
{
    return a->data.zz_poly.length >= SC_SQR_TOOM3_CUTOFF;
}

static int sc_zz_poly_use_sqr_karatsuba(const sc_value *a)
{
    return a->data.zz_poly.length >= SC_SQR_KARATSUBA_CUTOFF;
}

static int sc_zz_poly_use_ks(const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn;
    size_t bits_a = sc_zz_poly_max_abs_bits_raw(a);
    size_t bits_b = sc_zz_poly_max_abs_bits_raw(b);
    size_t bits = bits_a > bits_b ? bits_a : bits_b;

    return small >= SC_MUL_KS_CUTOFF && bits < small;
}

static int sc_zz_poly_use_toom3(const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, large = an > bn ? an : bn;

    return small >= SC_MUL_TOOM3_CUTOFF && large <= small + small / 2;
}

static int sc_zz_poly_use_karatsuba(const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, large = an > bn ? an : bn;
    size_t bits_a, bits_b, bits;

    if (large > 2 * small)
        return 0;
    bits_a = sc_zz_poly_max_abs_bits_raw(a);
    bits_b = sc_zz_poly_max_abs_bits_raw(b);
    bits = bits_a > bits_b ? bits_a : bits_b;
    if (bits < SC_MUL_KARATSUBA_LOW_BITS)
        return small >= SC_MUL_KARATSUBA_LOW_BITS_CUTOFF;
    return small >= SC_MUL_KARATSUBA_CUTOFF;
}

static int sc_zz_poly_balanced(const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, large = an > bn ? an : bn;

    return large <= small + small / 2;
}

static int sc_zz_poly_use_ntt_cutoff(const sc_value *a, const sc_value *b,
                                     size_t cutoff)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, np, ba, bb, bits;

    if (small < cutoff || !sc_zz_poly_balanced(a, b))
        return 0;
    ba = sc_zz_poly_max_abs_bits_raw(a);
    bb = sc_zz_poly_max_abs_bits_raw(b);
    bits = ba > bb ? ba : bb;
    if (bits >= small)
        return 0;
    np = sc_zz_poly_ntt_nprimes(a, b);
    return np != 0 && small / np >= cutoff;
}

static int sc_zz_poly_use_ssa_cutoff(const sc_value *a, const sc_value *b,
                                     size_t cutoff)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, ba, bb, bits;

    if (small < cutoff || !sc_zz_poly_balanced(a, b))
        return 0;
    ba = sc_zz_poly_max_abs_bits_raw(a);
    bb = sc_zz_poly_max_abs_bits_raw(b);
    bits = ba > bb ? ba : bb;
    return bits >= small;
}

static int sc_zz_poly_use_ntt(const sc_value *a, const sc_value *b)
{
    return sc_zz_poly_use_ntt_cutoff(a, b, SC_MUL_NTT_CUTOFF);
}

static int sc_zz_poly_use_ssa(const sc_value *a, const sc_value *b)
{
    return sc_zz_poly_use_ssa_cutoff(a, b, SC_MUL_SSA_CUTOFF);
}

static int sc_zz_poly_use_short_ntt(const sc_value *a, const sc_value *b,
                                    size_t n, size_t cutoff)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, ba, bb, bits, np;

    if (n < cutoff || !sc_zz_poly_balanced(a, b))
        return 0;
    ba = sc_zz_poly_max_abs_bits_raw(a);
    bb = sc_zz_poly_max_abs_bits_raw(b);
    bits = ba > bb ? ba : bb;
    if (bits >= small)
        return 0;
    np = sc_zz_poly_ntt_nprimes(a, b);
    return np != 0 && n / np >= cutoff;
}

static int sc_zz_poly_use_short_ssa(const sc_value *a, const sc_value *b,
                                    size_t n, size_t cutoff)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, ba, bb, bits;

    if (n < cutoff || !sc_zz_poly_balanced(a, b))
        return 0;
    ba = sc_zz_poly_max_abs_bits_raw(a);
    bb = sc_zz_poly_max_abs_bits_raw(b);
    bits = ba > bb ? ba : bb;
    return bits >= small;
}

static int sc_zz_poly_use_mulmid_ntt(const sc_value *a, const sc_value *b, size_t n)
{
    size_t ba, bb, bits, np;

    if (n < SC_MULMID_NTT_CUTOFF || n > SIZE_MAX / 2 + 1 ||
        a->data.zz_poly.length > 2 * n - 1 || b->data.zz_poly.length > n)
        return 0;
    ba = sc_zz_poly_max_abs_bits_raw(a);
    bb = sc_zz_poly_max_abs_bits_raw(b);
    bits = ba > bb ? ba : bb;
    if (bits >= n)
        return 0;
    np = sc_zz_poly_mulmid_ntt_nprimes(a, b, n);
    return np != 0 && n / np >= SC_MULMID_NTT_CUTOFF;
}

static int sc_zz_poly_use_mulmid_ssa(const sc_value *a, const sc_value *b, size_t n)
{
    size_t ba, bb, bits;

    if (n < SC_MULMID_SSA_CUTOFF || n > SIZE_MAX / 2 + 1 ||
        a->data.zz_poly.length > 2 * n - 1 || b->data.zz_poly.length > n)
        return 0;
    ba = sc_zz_poly_max_abs_bits_raw(a);
    bb = sc_zz_poly_max_abs_bits_raw(b);
    bits = ba > bb ? ba : bb;
    return bits >= n;
}

static sc_value *sc_zz_poly_fft_window(sc_context *ctx, const sc_value *a,
                                       const sc_value *b, size_t start, size_t n,
                                       int ntt)
{
    sc_value *p = ntt ? sc_zz_poly_mul_ntt(ctx, a, b) :
                        sc_zz_poly_mul_ssa(ctx, a, b);
    sc_value view, *r;

    if (p == NULL)
        return NULL;
    view = sc_zz_poly_view(p, start, n);
    r = sc_value_copy_checked(ctx, &view);
    sc_value_free(p);
    return r;
}

sc_value *sc_zz_poly_add(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_add_impl(ctx, a, b);
}

sc_value *sc_zz_poly_sub(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_sub_impl(ctx, a, b);
}

sc_value *sc_zz_poly_sqr(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    if (a->data.zz_poly.length <= 1)
        return sc_zz_poly_sqr_classical(ctx, a);
    if (sc_zz_poly_use_ntt_cutoff(a, a, SC_SQR_NTT_CUTOFF))
        return sc_zz_poly_sqr_ntt(ctx, a);
    if (sc_zz_poly_use_ssa_cutoff(a, a, SC_SQR_SSA_CUTOFF))
        return sc_zz_poly_sqr_ssa(ctx, a);
    if (sc_zz_poly_use_sqr_ks(a))
        return sc_zz_poly_sqr_ks(ctx, a, sc_zz_poly_ks_bits(a, a));
    if (sc_zz_poly_use_sqr_toom3(a))
        return sc_zz_poly_sqr_toom3(ctx, a);
    if (sc_zz_poly_use_sqr_karatsuba(a))
        return sc_zz_poly_sqr_karatsuba(ctx, a);
    return sc_zz_poly_sqr_classical(ctx, a);
}

sc_value *sc_zz_poly_mul(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (a == b)
        return sc_zz_poly_sqr(ctx, a);
    if (a->data.zz_poly.length == 0 || b->data.zz_poly.length == 0)
        return sc_zz_poly_mul_classical(ctx, a, b);
    if (a->data.zz_poly.length == 1 || b->data.zz_poly.length == 1)
        return sc_zz_poly_mul_constant(ctx, a, b);
    if (sc_zz_poly_use_ntt(a, b))
        return sc_zz_poly_mul_ntt(ctx, a, b);
    if (sc_zz_poly_use_ssa(a, b))
        return sc_zz_poly_mul_ssa(ctx, a, b);
    if (sc_zz_poly_use_ks(a, b))
        return sc_zz_poly_mul_ks(ctx, a, b, sc_zz_poly_ks_bits(a, b));
    if (sc_zz_poly_use_toom3(a, b)) {
        size_t large = a->data.zz_poly.length > b->data.zz_poly.length ?
                       a->data.zz_poly.length : b->data.zz_poly.length;
        sc_zz_poly_toom3_ws ws;

        if (!sc_zz_poly_toom3_ws_init(ctx, &ws, a, b, (large + 2) / 3))
            return NULL;
        return sc_zz_poly_mul_toom3(ctx, &ws);
    }
    if (sc_zz_poly_use_karatsuba(a, b))
        return sc_zz_poly_mul_karatsuba(ctx, a, b);
    return sc_zz_poly_mul_classical(ctx, a, b);
}

sc_value *sc_zz_poly_mullow(sc_context *ctx, const sc_value *a, const sc_value *b,
                            size_t n)
{
    size_t an, bn, total;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (n == 0 || an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    total = an + bn - 1;
    if (n >= total)
        return sc_zz_poly_mul(ctx, a, b);
    if (sc_zz_poly_use_short_ntt(a, b, n, SC_MULLOW_NTT_CUTOFF))
        return sc_zz_poly_fft_window(ctx, a, b, 0, n, 1);
    if (sc_zz_poly_use_short_ssa(a, b, n, SC_MULLOW_SSA_CUTOFF))
        return sc_zz_poly_fft_window(ctx, a, b, 0, n, 0);
    if (n <= SC_MULLOW_DC_CUTOFF)
        return sc_zz_poly_mullow_classical(ctx, a, b, n);
    return sc_zz_poly_mullow_dc(ctx, a, b, n);
}

sc_value *sc_zz_poly_mulhigh(sc_context *ctx, const sc_value *a, const sc_value *b,
                             size_t n)
{
    size_t an, bn, total;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (n == 0 || an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    total = an + bn - 1;
    if (n >= total)
        return sc_zz_poly_mul(ctx, a, b);
    if (sc_zz_poly_use_short_ntt(a, b, n, SC_MULHIGH_NTT_CUTOFF))
        return sc_zz_poly_fft_window(ctx, a, b, total - n, n, 1);
    if (sc_zz_poly_use_short_ssa(a, b, n, SC_MULHIGH_SSA_CUTOFF))
        return sc_zz_poly_fft_window(ctx, a, b, total - n, n, 0);
    return sc_zz_poly_mulhigh_reverse(ctx, a, b, n);
}

sc_value *sc_zz_poly_mulmid_balanced(sc_context *ctx, const sc_value *a,
                                     const sc_value *b, size_t n)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (n == 0 || a->data.zz_poly.length == 0 || b->data.zz_poly.length == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (sc_zz_poly_use_mulmid_ntt(a, b, n))
        return sc_zz_poly_mulmid_ntt(ctx, a, b, n);
    if (sc_zz_poly_use_mulmid_ssa(a, b, n))
        return sc_zz_poly_mulmid_ssa(ctx, a, b, n);
    if (n <= SC_MULMID_CLASSICAL_CUTOFF)
        return sc_zz_poly_mulmid_classical(ctx, a, b, n - 1, n);
    if (n >= SC_MULMID_TOOM63_CUTOFF) {
        sc_zz_poly_toom63_ws ws;

        if (n % 3 != 0)
            return sc_zz_poly_mulmid_toom63_tail(ctx, a, b, n);
        if (!sc_zz_poly_toom63_ws_init(ctx, &ws, a, b, n / 3))
            return NULL;
        return sc_zz_poly_mulmid_toom63(ctx, &ws);
    }
    if (n & 1)
        return sc_zz_poly_mulmid_toom42_odd(ctx, a, b, n);
    return sc_zz_poly_mulmid_toom42(ctx, a, b, n);
}

sc_value *sc_zz_poly_mulmid(sc_context *ctx, const sc_value *a, const sc_value *b,
                            size_t start, size_t n)
{
    size_t an, bn, total, count, small, shift;
    sc_value *p, *r;
    sc_value view;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (n == 0 || an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    total = an + bn - 1;
    if (start >= total)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    count = n < total - start ? n : total - start;
    if (start == 0)
        return sc_zz_poly_mullow(ctx, a, b, count);
    if (start + count == total)
        return sc_zz_poly_mulhigh(ctx, a, b, count);
    if (start >= count - 1 && bn <= count) {
        shift = start - count + 1;
        view = sc_zz_poly_view(a, shift, 2 * count - 1);
        return sc_zz_poly_mulmid_balanced(ctx, &view, b, count);
    }
    if (start >= count - 1 && an <= count) {
        shift = start - count + 1;
        view = sc_zz_poly_view(b, shift, 2 * count - 1);
        return sc_zz_poly_mulmid_balanced(ctx, &view, a, count);
    }
    small = an < bn ? an : bn;
    if (count <= SC_MULMID_GENERAL_CLASSICAL_CUTOFF ||
        small <= SC_MULMID_GENERAL_CLASSICAL_CUTOFF)
        return sc_zz_poly_mulmid_classical(ctx, a, b, start, count);
    p = sc_zz_poly_mul(ctx, a, b);
    if (p == NULL)
        return NULL;
    view = sc_zz_poly_view(p, start, count);
    r = sc_value_copy_checked(ctx, &view);
    sc_value_free(p);
    return r;
}

sc_value *sc_zz_poly_inv_series(sc_context *ctx, const sc_value *a, size_t n)
{
    if (a == NULL)
        return NULL;
    if (n == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (a->data.zz_poly.length == 0 ||
        mpz_cmpabs_ui(a->data.zz_poly.coeff[0], 1) != 0) {
        sc_set_error(ctx, "series inverse requires constant coefficient +/-1");
        return NULL;
    }
    if (n <= SC_INV_SERIES_NEWTON_CUTOFF)
        return sc_zz_poly_inv_series_classical(ctx, a, n);
    return sc_zz_poly_inv_series_newton(ctx, a, n);
}

sc_value *sc_zz_poly_scalar_mul(sc_context *ctx,
                                const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_scalar_mul_impl(ctx, a, b);
}

sc_value *sc_zz_poly_neg(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_neg_impl(ctx, a);
}

sc_value *sc_zz_poly_pow(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (mpz_sgn(b->data.z) < 0 || !mpz_fits_ulong_p(b->data.z)) {
        sc_set_error(ctx, "exponent must be a nonnegative machine integer");
        return NULL;
    }
    return sc_zz_poly_pow_binary(ctx, a, mpz_get_ui(b->data.z));
}

sc_value *sc_zz_poly_divrem_full(sc_context *ctx,
                                  const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t qn = an >= bn ? an - bn + 1 : 0;

    if (qn >= SC_DIVREM_DC_CUTOFF && bn > 1)
        return sc_zz_poly_divrem_dc(ctx, a, b);
    return sc_zz_poly_divrem_classical(ctx, a, b);
}

sc_value *sc_zz_poly_quo_mulders(sc_context *ctx,
                                  const sc_value *a, const sc_value *b)
{
    size_t an, bn, qn, start;
    sc_value av, bv;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (bn == 0) {
        sc_set_error(ctx, "polynomial division by zero");
        return NULL;
    }
    qn = an >= bn ? an - bn + 1 : 0;
    if (qn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (qn > bn)
        return sc_zz_poly_quo_dc(ctx, a, b);
    start = an - (2 * qn - 1);
    av = sc_zz_poly_view(a, start, 2 * qn - 1);
    bv = sc_zz_poly_view(b, bn - qn, qn);
    return sc_zz_poly_quo_mulders_balanced(ctx, &av, &bv);
}

sc_value *sc_zz_poly_quo(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, qn;

    if (a == NULL || b == NULL)
        return NULL;
    bn = b->data.zz_poly.length;
    if (bn == 0) {
        sc_set_error(ctx, "polynomial division by zero");
        return NULL;
    }
    an = a->data.zz_poly.length;
    qn = an >= bn ? an - bn + 1 : 0;
    if (qn >= SC_QUO_NEWTON_CUTOFF && bn > 1 && sc_zz_poly_unit_lead(b))
        return sc_zz_poly_quo_newton(ctx, a, b);
    if (qn >= SC_QUO_MULDERS_CUTOFF && qn <= bn && bn > 1)
        return sc_zz_poly_quo_mulders(ctx, a, b);
    if (qn >= SC_QUO_DC_CUTOFF && bn > 1)
        return sc_zz_poly_quo_dc(ctx, a, b);
    return sc_zz_poly_quo_classical(ctx, a, b);
}

sc_value *sc_zz_poly_divrem(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, qn;

    if (a == NULL || b == NULL)
        return NULL;
    bn = b->data.zz_poly.length;
    if (bn == 0) {
        sc_set_error(ctx, "polynomial division by zero");
        return NULL;
    }
    an = a->data.zz_poly.length;
    qn = an >= bn ? an - bn + 1 : 0;
    if (qn >= SC_DIVREM_NEWTON_CUTOFF && bn > 1 && sc_zz_poly_unit_lead(b))
        return sc_zz_poly_divrem_newton(ctx, a, b);
    if (qn >= SC_DIVREM_MULDERS_CUTOFF && qn <= bn && bn > 1)
        return sc_zz_poly_divrem_mulders(ctx, a, b);
    if (qn >= SC_DIVREM_DC_CUTOFF && bn > 1)
        return sc_zz_poly_divrem_dc(ctx, a, b);
    return sc_zz_poly_divrem_classical(ctx, a, b);
}

static sc_value *sc_zz_poly_divexact_via_divrem(sc_context *ctx,
                                                 const sc_value *a,
                                                 const sc_value *b)
{
    sc_value *qr = sc_zz_poly_divrem(ctx, a, b), *q;

    if (qr == NULL)
        return NULL;
    if (qr->data.pair.second->data.zz_poly.length != 0) {
        sc_set_error(ctx, "polynomial quotient is not exact in ZZ[x]");
        sc_value_free(qr);
        return NULL;
    }
    q = qr->data.pair.first;
    qr->data.pair.first = NULL;
    sc_value_free(qr);
    return q;
}

sc_value *sc_zz_poly_divexact(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, qn, v, ln, midn;
    sc_value av, bv, *q, *mid;

    if (a == NULL || b == NULL)
        return NULL;
    bn = b->data.zz_poly.length;
    if (bn == 0) {
        sc_set_error(ctx, "polynomial division by zero");
        return NULL;
    }
    an = a->data.zz_poly.length;
    qn = an >= bn ? an - bn + 1 : 0;
    if (qn < SC_DIVEXACT_BIDIR_CUTOFF || sc_zz_poly_unit_lead(b))
        return sc_zz_poly_divexact_via_divrem(ctx, a, b);
    v = sc_zz_poly_valuation(b);
    if (!sc_zz_poly_zero_prefix(a, v)) {
        sc_set_error(ctx, "polynomial quotient is not exact in ZZ[x]");
        return NULL;
    }
    av = sc_zz_poly_view(a, v, an - v);
    bv = sc_zz_poly_view(b, v, bn - v);
    q = sc_zz_poly_quo_bidirectional(ctx, &av, &bv);
    if (q == NULL)
        return NULL;
    ln = (qn + 1) / 2;
    midn = bv.data.zz_poly.length - 1;
    mid = sc_zz_poly_mulmid(ctx, &bv, q, ln, midn);
    if (mid == NULL || !sc_zz_poly_window_equal(&av, ln, mid, midn)) {
        sc_set_error(ctx, "polynomial quotient is not exact in ZZ[x]");
        return sc_value_free_many_null(2, q, mid);
    }
    sc_value_free(mid);
    return q;
}

sc_value *sc_zz_poly_pseudodiv(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, qn;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (bn == 0) {
        sc_set_error(ctx, "polynomial pseudo-division by zero");
        return NULL;
    }
    qn = an >= bn ? an - bn + 1 : 0;
    if (qn >= SC_PSEUDODIV_FAST_CUTOFF && bn > 1)
        return sc_zz_poly_pseudodiv_fast(ctx, a, b);
    return sc_zz_poly_pseudodiv_impl(ctx, a, b);
}

sc_value *sc_zz_poly_pseudorem(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, qn;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (bn == 0) {
        sc_set_error(ctx, "polynomial pseudo-remainder by zero");
        return NULL;
    }
    qn = an >= bn ? an - bn + 1 : 0;
    if (qn >= SC_PSEUDOREM_FAST_CUTOFF && bn > 1)
        return sc_zz_poly_pseudorem_fast(ctx, a, b);
    return sc_zz_poly_pseudorem_classical(ctx, a, b);
}

sc_value *sc_zz_poly_scalar_divexact(sc_context *ctx, const sc_value *a,
                                     const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_scalar_divexact_impl(ctx, a, b);
}

sc_value *sc_zz_poly_content(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_content_impl(ctx, a);
}

sc_value *sc_zz_poly_primitive_part(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_primitive_part_impl(ctx, a);
}

sc_value *sc_zz_poly_gcd_pseudo(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_gcd_pseudo_impl(ctx, a, b);
}

sc_value *sc_zz_poly_gcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t n;

    if (a == NULL || b == NULL)
        return NULL;
    n = a->data.zz_poly.length > b->data.zz_poly.length ?
        a->data.zz_poly.length : b->data.zz_poly.length;
    if (n < SC_GCD_SUBRESULTANT_CUTOFF)
        return sc_zz_poly_gcd_pseudo_impl(ctx, a, b);
    return sc_zz_poly_gcd_subresultant_impl(ctx, a, b);
}

sc_value *sc_zz_poly_gcd_hgcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_gcd_hgcd_impl(ctx, a, b);
}

sc_value *sc_zz_poly_gcd_lr(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_gcd_lr_impl(ctx, a, b);
}

sc_value *sc_zz_poly_resultant(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t d;

    if (a == NULL || b == NULL)
        return NULL;
    if (a->data.zz_poly.length == 0 || b->data.zz_poly.length == 0)
        return sc_zz_poly_resultant_subresultant_impl(ctx, a, b);
    d = a->data.zz_poly.length + b->data.zz_poly.length - 2;
    if (d < SC_RESULTANT_SUBRESULTANT_CUTOFF)
        return sc_zz_poly_resultant_bareiss_impl(ctx, a, b);
    return sc_zz_poly_resultant_subresultant_impl(ctx, a, b);
}


sc_value *sc_zz_poly_xgcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r, *h, *cof, *u, *v, *z, *o;

    if (a == NULL || b == NULL)
        return NULL;
    if (a->data.zz_poly.length == 0 || b->data.zz_poly.length == 0) {
        h = sc_value_copy_checked(ctx, a->data.zz_poly.length == 0 ? b : a);
        z = sc_value_new_zz_poly_checked(ctx, a->parent, 0);
        o = sc_value_new_zz_poly_checked(ctx, a->parent, 1);
        if (o != NULL)
            mpz_set_ui(o->data.zz_poly.coeff[0], 1);
        if (h == NULL || z == NULL || o == NULL)
            return sc_value_free_many_null(3, h, z, o);
        u = a->data.zz_poly.length == 0 ? z : o;
        v = a->data.zz_poly.length == 0 ? o : z;
        cof = sc_value_new_pair_take_checked(ctx, u, v);
        return cof == NULL ? sc_value_free_many_null(1, h) :
                             sc_value_new_pair_take_checked(ctx, h, cof);
    }
    if (a->data.zz_poly.length >= b->data.zz_poly.length)
        return sc_zz_poly_xgcd_subresultant_impl(ctx, a, b);
    r = sc_zz_poly_xgcd_subresultant_impl(ctx, b, a);
    if (r == NULL)
        return NULL;
    cof = r->data.pair.second;
    u = cof->data.pair.first;
    cof->data.pair.first = cof->data.pair.second;
    cof->data.pair.second = u;
    return r;
}

sc_value *sc_zz_poly_evaluate_horner(sc_context *ctx, const sc_value *a,
                                     const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_evaluate_horner_impl(ctx, a, b);
}

sc_value *sc_zz_poly_evaluate_divconquer(sc_context *ctx, const sc_value *a,
                                         const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_evaluate_divconquer_impl(ctx, a, b);
}

sc_value *sc_zz_poly_evaluate(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (a->data.zz_poly.length < SC_EVALUATE_DC_CUTOFF ||
        mpz_cmpabs_ui(b->data.z, 1) <= 0)
        return sc_zz_poly_evaluate_horner_impl(ctx, a, b);
    return sc_zz_poly_evaluate_divconquer_impl(ctx, a, b);
}

sc_value *sc_zz_poly_compose_horner(sc_context *ctx, const sc_value *a,
                                    const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_compose_horner_impl(ctx, a, b);
}

sc_value *sc_zz_poly_compose_divconquer(sc_context *ctx, const sc_value *a,
                                        const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_compose_divconquer_impl(ctx, a, b);
}

static int sc_zz_poly_compose_dc_phase(size_t n)
{
    size_t p = 1;

    while (p < n && p <= SIZE_MAX / 2)
        p <<= 1;
    return n >= p - p / 8;
}

sc_value *sc_zz_poly_compose(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (a->data.zz_poly.length < SC_COMPOSE_DC_CUTOFF ||
        b->data.zz_poly.length <= 1 ||
        !sc_zz_poly_compose_dc_phase(a->data.zz_poly.length))
        return sc_zz_poly_compose_horner_impl(ctx, a, b);
    return sc_zz_poly_compose_divconquer_impl(ctx, a, b);
}

sc_value *sc_zz_poly_taylor_shift_horner(sc_context *ctx, const sc_value *a,
                                         const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_taylor_shift_horner_impl(ctx, a, b);
}

sc_value *sc_zz_poly_taylor_shift_divconquer(sc_context *ctx, const sc_value *a,
                                             const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_taylor_shift_divconquer_impl(ctx, a, b);
}

sc_value *sc_zz_poly_taylor_shift_convolution(sc_context *ctx,
                                              const sc_value *a,
                                              const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    return sc_zz_poly_taylor_shift_convolution_impl(ctx, a, b);
}

sc_value *sc_zz_poly_taylor_shift(sc_context *ctx, const sc_value *a,
                                  const sc_value *b)
{
    if (a == NULL || b == NULL)
        return NULL;
    if (mpz_cmpabs_ui(b->data.z, 1) <= 0)
        return sc_zz_poly_taylor_shift_horner_impl(ctx, a, b);
    if (sc_zz_poly_taylor_shift_convolution_preferred(a, b))
        return sc_zz_poly_taylor_shift_convolution_impl(ctx, a, b);
    if (a->data.zz_poly.length < SC_TAYLOR_DC_CUTOFF)
        return sc_zz_poly_taylor_shift_horner_impl(ctx, a, b);
    return sc_zz_poly_taylor_shift_divconquer_impl(ctx, a, b);
}


sc_value *sc_zz_poly_derivative(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_derivative_impl(ctx, a);
}

sc_value *sc_zz_poly_nth_derivative(sc_context *ctx, const sc_value *a, size_t n)
{
    if (a == NULL)
        return NULL;
    if (n == 1)
        return sc_zz_poly_derivative_impl(ctx, a);
    return sc_zz_poly_nth_derivative_impl(ctx, a, n);
}

sc_value *sc_zz_poly_discriminant(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_discriminant_impl(ctx, a);
}

sc_value *sc_zz_poly_is_squarefree(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_is_squarefree_impl(ctx, a);
}

sc_value *sc_zz_poly_squarefree_part(sc_context *ctx, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    return sc_zz_poly_squarefree_part_impl(ctx, a);
}

sc_value *sc_zz_poly_degree(sc_context *ctx, const sc_value *a)
{
    return a == NULL ? NULL : sc_zz_poly_degree_impl(ctx, a);
}

sc_value *sc_zz_poly_leading_coefficient(sc_context *ctx, const sc_value *a)
{
    return a == NULL ? NULL : sc_zz_poly_leading_coefficient_impl(ctx, a);
}

sc_value *sc_zz_poly_constant_coefficient(sc_context *ctx, const sc_value *a)
{
    return a == NULL ? NULL : sc_zz_poly_constant_coefficient_impl(ctx, a);
}

sc_value *sc_zz_poly_coeff(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_coeff_impl(ctx, a, n);
}

sc_value *sc_zz_poly_reverse(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_reverse_impl(ctx, a, n);
}

sc_value *sc_zz_poly_truncate_copy(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_truncate_impl(ctx, a, n);
}

sc_value *sc_zz_poly_shift_left(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_shift_left_impl(ctx, a, n);
}

sc_value *sc_zz_poly_shift_right(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_shift_right_impl(ctx, a, n);
}

sc_value *sc_zz_poly_height(sc_context *ctx, const sc_value *a)
{
    return a == NULL ? NULL : sc_zz_poly_height_impl(ctx, a);
}

sc_value *sc_zz_poly_max_abs_bits(sc_context *ctx, const sc_value *a)
{
    return a == NULL ? NULL : sc_zz_poly_max_abs_bits_impl(ctx, a);
}

sc_value *sc_zz_poly_inflate(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_inflate_impl(ctx, a, n);
}

sc_value *sc_zz_poly_deflation(sc_context *ctx, const sc_value *a)
{
    return a == NULL ? NULL : sc_zz_poly_deflation_impl(ctx, a);
}

sc_value *sc_zz_poly_deflate(sc_context *ctx, const sc_value *a, size_t n)
{
    return a == NULL ? NULL : sc_zz_poly_deflate_impl(ctx, a, n);
}
