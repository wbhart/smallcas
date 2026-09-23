#ifndef SMALLCAS_FFT_H
#define SMALLCAS_FFT_H

#include <gmp.h>
#include <stddef.h>

typedef struct sc_fft_mod {
    mpz_t mod;
    mpz_t root;
    mpz_t root_inv;
    unsigned depth;
} sc_fft_mod;

int sc_fft_mod_init(sc_fft_mod *m, mpz_srcptr mod, mpz_srcptr root,
                    unsigned depth);
void sc_fft_mod_clear(sc_fft_mod *m);
void sc_fft_root_power(mpz_t r, size_t e, int inverse, const sc_fft_mod *m);
void sc_fft_div_2exp(mpz_t a, size_t k, const sc_fft_mod *m);

static inline void sc_fft_set_z(mpz_t r, mpz_srcptr a, const sc_fft_mod *m)
{
    mpz_mod(r, a, m->mod);
}

static inline void sc_fft_add(mpz_t r, mpz_srcptr a, mpz_srcptr b,
                              const sc_fft_mod *m)
{
    mpz_add(r, a, b);
    if (mpz_cmp(r, m->mod) >= 0)
        mpz_sub(r, r, m->mod);
}

static inline void sc_fft_sub(mpz_t r, mpz_srcptr a, mpz_srcptr b,
                              const sc_fft_mod *m)
{
    mpz_sub(r, a, b);
    if (mpz_sgn(r) < 0)
        mpz_add(r, r, m->mod);
}

static inline void sc_fft_mul(mpz_t r, mpz_srcptr a, mpz_srcptr b,
                              const sc_fft_mod *m)
{
    mpz_mul(r, a, b);
    mpz_mod(r, r, m->mod);
}

static inline void sc_fft_addsub(mpz_t a, mpz_t b, mpz_t t,
                                 const sc_fft_mod *m)
{
    mpz_set(t, a);
    sc_fft_add(a, a, b, m);
    sc_fft_sub(b, t, b, m);
}

#define sc_fft_mul_root(r, a, w, m) sc_fft_mul((r), (a), (w), (m))

#endif
