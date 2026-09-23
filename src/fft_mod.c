#include "smallcas_fft.h"

static void sc_fft_set_size_t(mpz_t z, size_t a)
{
    mpz_import(z, 1, 1, sizeof(a), 0, 0, &a);
}

int sc_fft_mod_init(sc_fft_mod *m, mpz_srcptr mod, mpz_srcptr root,
                    unsigned depth)
{
    mpz_t e, t;
    int ok = 0;

    mpz_inits(m->mod, m->root, m->root_inv, e, t, NULL);
    if (depth == 0 || mpz_cmp_ui(mod, 2) <= 0 || mpz_even_p(mod))
        goto done;
    mpz_set(m->mod, mod);
    mpz_mod(m->root, root, mod);
    mpz_set_ui(e, 0);
    mpz_setbit(e, depth - 1);
    mpz_powm(t, m->root, e, m->mod);
    mpz_sub_ui(e, m->mod, 1);
    if (mpz_cmp(t, e) != 0 || mpz_invert(m->root_inv, m->root, m->mod) == 0)
        goto done;
    m->depth = depth;
    ok = 1;
done:
    mpz_clears(e, t, NULL);
    if (!ok)
        sc_fft_mod_clear(m);
    return ok;
}

void sc_fft_mod_clear(sc_fft_mod *m)
{
    mpz_clears(m->mod, m->root, m->root_inv, NULL);
    m->depth = 0;
}

void sc_fft_root_power(mpz_t r, size_t e, int inverse, const sc_fft_mod *m)
{
    mpz_t z;

    mpz_init(z);
    sc_fft_set_size_t(z, e);
    mpz_powm(r, inverse ? m->root_inv : m->root, z, m->mod);
    mpz_clear(z);
}

void sc_fft_div_2exp(mpz_t a, size_t k, const sc_fft_mod *m)
{
    while (k-- != 0) {
        if (mpz_odd_p(a))
            mpz_add(a, a, m->mod);
        mpz_fdiv_q_2exp(a, a, 1);
    }
}
