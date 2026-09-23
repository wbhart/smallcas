#include "smallcas_fft.h"

#include <stdio.h>

static int test_fermat_mod(void)
{
    sc_fft_mod m;
    mpz_t p, root, a, b, c, d, t;
    int ok = 1;

    mpz_inits(p, root, a, b, c, d, t, NULL);
    mpz_set_ui(p, 1);
    mpz_setbit(p, 32);
    mpz_set_ui(root, 2);
    if (!sc_fft_mod_init(&m, p, root, 6))
        ok = 0;
    if (ok) {
        mpz_sub_ui(a, p, 2);
        mpz_set_ui(b, 5);
        sc_fft_add(c, a, b, &m);
        ok &= mpz_cmp_ui(c, 3) == 0;
        sc_fft_sub(c, a, b, &m);
        mpz_sub_ui(t, p, 7);
        ok &= mpz_cmp(c, t) == 0;
        sc_fft_mul(c, a, b, &m);
        mpz_mul(t, a, b);
        mpz_mod(t, t, p);
        ok &= mpz_cmp(c, t) == 0;
        mpz_set(c, a);
        mpz_set(d, b);
        sc_fft_addsub(c, d, t, &m);
        ok &= mpz_cmp_ui(c, 3) == 0;
        mpz_sub_ui(a, p, 7);
        ok &= mpz_cmp(d, a) == 0;
        sc_fft_root_power(c, 17, 0, &m);
        sc_fft_root_power(d, 17, 1, &m);
        sc_fft_mul(c, c, d, &m);
        ok &= mpz_cmp_ui(c, 1) == 0;
        mpz_set_ui(a, 123456789);
        mpz_mul_2exp(b, a, 11);
        mpz_mod(b, b, p);
        sc_fft_div_2exp(b, 11, &m);
        ok &= mpz_cmp(a, b) == 0;
        sc_fft_mod_clear(&m);
    }
    mpz_clears(p, root, a, b, c, d, t, NULL);
    return ok;
}

static int test_crt_prime(void)
{
    sc_fft_mod m;
    mpz_t p, g, root, t;
    int ok;

    mpz_inits(p, g, root, t, NULL);
    mpz_set_ui(p, 2013265921UL);
    mpz_set_ui(g, 31);
    mpz_powm_ui(root, g, 15, p);
    ok = sc_fft_mod_init(&m, p, root, 27);
    if (ok) {
        sc_fft_root_power(t, (size_t)1 << 26, 0, &m);
        mpz_add_ui(t, t, 1);
        ok = mpz_cmp(t, p) == 0;
        sc_fft_mod_clear(&m);
    }
    mpz_clears(p, g, root, t, NULL);
    return ok;
}

static int test_reject_bad_root(void)
{
    sc_fft_mod m;
    mpz_t p, root;
    int ok;

    mpz_inits(p, root, NULL);
    mpz_set_ui(p, 97);
    mpz_set_ui(root, 1);
    ok = !sc_fft_mod_init(&m, p, root, 5);
    mpz_clears(p, root, NULL);
    return ok;
}

int main(void)
{
    if (!test_fermat_mod() || !test_crt_prime() || !test_reject_bad_root()) {
        fprintf(stderr, "fft modulus tests failed\n");
        return 1;
    }
    printf("fft modulus tests passed\n");
    return 0;
}
