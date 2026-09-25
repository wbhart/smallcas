#ifndef SMALLCAS_FFT_H
#define SMALLCAS_FFT_H

#include <gmp.h>
#include <stddef.h>
#include <stdint.h>

typedef struct sc_fft_mod {
    uint64_t c;
    mp_bitcnt_t bits;
    mp_size_t n;
    mp_ptr mod, root, root_inv;
    unsigned depth;
    int fermat;
} sc_fft_mod;

typedef struct sc_fft_plan {
    size_t len;
    unsigned logn;
    mp_ptr fwd, inv;
    size_t root_stride;
} sc_fft_plan;

int sc_fft_mod_init(sc_fft_mod *m, uint64_t c, mp_bitcnt_t bits,
                    uint64_t root, unsigned depth);
void sc_fft_mod_clear(sc_fft_mod *m);
int sc_fft_plan_init(sc_fft_plan *p, unsigned logn, const sc_fft_mod *m);
void sc_fft_plan_clear(sc_fft_plan *p);
void sc_fft_forward(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                    mp_ptr scratch);
void sc_fft_inverse(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                    mp_ptr scratch);
void sc_fft_forward_mfa(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                        mp_ptr scratch);
void sc_fft_inverse_mfa(mp_ptr a, const sc_fft_plan *p, const sc_fft_mod *m,
                        mp_ptr scratch);
void sc_fft_forward_tft(mp_ptr a, size_t z, size_t n, const sc_fft_plan *p,
                        const sc_fft_mod *m, mp_ptr scratch);
void sc_fft_inverse_tft(mp_ptr a, size_t n, const sc_fft_plan *p,
                        const sc_fft_mod *m, mp_ptr scratch);
void sc_fft_mul(mp_ptr r, mp_srcptr a, mp_srcptr b,
                const sc_fft_mod *m, mp_ptr scratch);
void sc_fft_sqr(mp_ptr r, mp_srcptr a, const sc_fft_mod *m, mp_ptr scratch);
void sc_fft_div_2exp(mp_ptr a, size_t k, const sc_fft_mod *m);
void sc_fft_mul_2exp(mp_ptr r, mp_srcptr a, size_t e,
                     const sc_fft_mod *m, mp_ptr scratch);

static inline mp_ptr sc_fft_entry(mp_ptr a, size_t i, const sc_fft_mod *m)
{
    return a + i * (size_t)m->n;
}

static inline mp_srcptr sc_fft_entry_const(mp_srcptr a, size_t i,
                                            const sc_fft_mod *m)
{
    return a + i * (size_t)m->n;
}

static inline void sc_fft_set_ui(mp_ptr r, mp_limb_t a, const sc_fft_mod *m)
{
    mpn_zero(r, m->n);
    r[0] = m->n == 1 ? a % m->mod[0] : a;
}

static inline int sc_fft_is_zero(mp_srcptr a, const sc_fft_mod *m)
{
    for (mp_size_t i = 0; i < m->n; i++)
        if (a[i] != 0)
            return 0;
    return 1;
}

static inline void sc_fft_neg(mp_ptr r, mp_srcptr a, const sc_fft_mod *m)
{
    if (sc_fft_is_zero(a, m))
        mpn_zero(r, m->n);
    else
        mpn_sub_n(r, m->mod, a, m->n);
}

static inline int sc_fft_equal(mp_srcptr a, mp_srcptr b, const sc_fft_mod *m)
{
    return mpn_cmp(a, b, m->n) == 0;
}

static inline int sc_fft_equal_ui(mp_srcptr a, mp_limb_t b,
                                  const sc_fft_mod *m)
{
    if (a[0] != b)
        return 0;
    for (mp_size_t i = 1; i < m->n; i++)
        if (a[i] != 0)
            return 0;
    return 1;
}

static inline void sc_fft_add(mp_ptr r, mp_srcptr a, mp_srcptr b,
                              const sc_fft_mod *m)
{
    mp_limb_t cy = mpn_add_n(r, a, b, m->n);
    if (cy || mpn_cmp(r, m->mod, m->n) >= 0)
        mpn_sub_n(r, r, m->mod, m->n);
}

static inline void sc_fft_sub(mp_ptr r, mp_srcptr a, mp_srcptr b,
                              const sc_fft_mod *m)
{
    mp_limb_t cy = mpn_sub_n(r, a, b, m->n);
    if (cy)
        mpn_add_n(r, r, m->mod, m->n);
}

static inline void sc_fft_addsub(mp_ptr a, mp_ptr b, mp_ptr t,
                                 const sc_fft_mod *m)
{
    mpn_copyi(t, a, m->n);
    sc_fft_add(a, a, b, m);
    sc_fft_sub(b, t, b, m);
}

#endif
