#include "smallcas_fft.h"

#include <stdio.h>
#include <stdlib.h>

static int test_moduli(void)
{
    sc_fft_mod a, b;
    int ok;

    ok = sc_fft_mod_init(&a, 1, 64, 2, 7);
    ok &= sc_fft_mod_init(&b, 15, 27, 440564289, 27);
    if (ok) {
        ok &= a.n == 2 && a.fermat && b.n == 1 && !b.fermat;
        ok &= a.mod[0] == 1 && a.mod[1] == 1;
        ok &= b.mod[0] == 2013265921UL;
    }
    if (ok) {
        sc_fft_mod_clear(&a);
        sc_fft_mod_clear(&b);
    }
    return ok;
}

static int test_arithmetic(void)
{
    sc_fft_mod m;
    mp_ptr a, b, c, d, s;
    int ok = sc_fft_mod_init(&m, 1, 64, 2, 7);

    if (!ok)
        return 0;
    s = calloc((size_t)9 * m.n + 1, sizeof(mp_limb_t));
    a = s;
    b = a + m.n;
    c = b + m.n;
    d = c + m.n;
    sc_fft_set_ui(a, 12345, &m);
    sc_fft_set_ui(b, 6789, &m);
    sc_fft_mul(c, a, b, &m, d + m.n);
    ok &= sc_fft_equal_ui(c, 12345UL * 6789UL, &m);
    sc_fft_sqr(d, a, &m, d + m.n);
    ok &= sc_fft_equal_ui(d, 12345UL * 12345UL, &m);
    sc_fft_div_2exp(c, 5, &m);
    for (int i = 0; i < 5; i++)
        sc_fft_add(c, c, c, &m);
    ok &= sc_fft_equal_ui(c, 12345UL * 6789UL, &m);
    mpn_copyi(a, m.mod, m.n);
    mpn_sub_1(a, a, m.n, 123);
    for (size_t e = 0; e < 128; e++) {
        mpn_copyi(b, a, m.n);
        for (size_t j = 0; j < e; j++)
            sc_fft_add(b, b, b, &m);
        sc_fft_mul_2exp(c, a, e, &m, d);
        ok &= sc_fft_equal(b, c, &m);
    }
    free(s);
    sc_fft_mod_clear(&m);
    return ok;
}


static int test_prime_square(void)
{
    sc_fft_mod m;
    mp_ptr a, c, work, s;
    int ok = sc_fft_mod_init(&m, 15, 27, 440564289, 27);

    if (!ok)
        return 0;
    s = calloc((size_t)7 * m.n + 1, sizeof(mp_limb_t));
    a = s;
    c = a + m.n;
    work = c + m.n;
    sc_fft_set_ui(a, 12345, &m);
    sc_fft_sqr(c, a, &m, work);
    ok &= sc_fft_equal_ui(c, 12345UL * 12345UL, &m);
    free(s);
    sc_fft_mod_clear(&m);
    return ok;
}

int main(void)
{
    if (!test_moduli() || !test_arithmetic() || !test_prime_square()) {
        fprintf(stderr, "fft modulus tests failed\n");
        return 1;
    }
    printf("fft modulus tests passed\n");
    return 0;
}
