#include "smallcas_fft.h"

#include <stdio.h>
#include <stdlib.h>

static int test_roundtrips(uint64_t c, mp_bitcnt_t bits, uint64_t root,
                           unsigned depth, unsigned maxlog)
{
    sc_fft_mod m;
    int ok = sc_fft_mod_init(&m, c, bits, root, depth);

    if (!ok)
        return 0;
    for (unsigned logn = 0; logn <= maxlog && ok; logn++) {
        sc_fft_plan p;
        mp_ptr a, b, d, work;

        ok = sc_fft_plan_init(&p, logn, &m);
        if (!ok)
            break;
        a = calloc(3 * p.len * (size_t)m.n, sizeof(mp_limb_t));
        b = a + p.len * (size_t)m.n;
        d = b + p.len * (size_t)m.n;
        work = calloc((size_t)4 * m.n + 1, sizeof(mp_limb_t));
        if (a == NULL || work == NULL) {
            free(work);
            free(a);
            sc_fft_plan_clear(&p);
            ok = 0;
            break;
        }
        for (size_t i = 0; i < p.len; i++) {
            mp_limb_t x = (mp_limb_t)(i * 6364136223846793005ULL +
                                      1442695040888963407ULL);
            sc_fft_set_ui(sc_fft_entry(a, i, &m), x, &m);
            mpn_copyi(sc_fft_entry(b, i, &m), sc_fft_entry(a, i, &m), m.n);
            mpn_copyi(sc_fft_entry(d, i, &m), sc_fft_entry(a, i, &m), m.n);
        }
        sc_fft_forward(b, &p, &m, work);
        sc_fft_forward_mfa(d, &p, &m, work);
        for (size_t i = 0; i < p.len; i++)
            ok &= sc_fft_equal(sc_fft_entry(b, i, &m),
                               sc_fft_entry(d, i, &m), &m);
        sc_fft_inverse(b, &p, &m, work);
        sc_fft_inverse_mfa(d, &p, &m, work);
        for (size_t i = 0; i < p.len; i++) {
            ok &= sc_fft_equal(sc_fft_entry(b, i, &m), sc_fft_entry(a, i, &m), &m);
            ok &= sc_fft_equal(sc_fft_entry(d, i, &m), sc_fft_entry(a, i, &m), &m);
        }
        free(work);
        free(a);
        sc_fft_plan_clear(&p);
    }
    sc_fft_mod_clear(&m);
    return ok;
}

static int test_convolution(void)
{
    static const mp_limb_t av[] = {1, 2, 3};
    static const mp_limb_t bv[] = {4, 5, 6};
    static const mp_limb_t cv[] = {4, 13, 28, 27, 18, 0, 0, 0};
    sc_fft_mod m;
    sc_fft_plan p;
    mp_ptr a, b, s, work;
    int ok = sc_fft_mod_init(&m, 15, 27, 440564289, 27);

    if (!ok || !sc_fft_plan_init(&p, 3, &m))
        return 0;
    a = calloc(2 * p.len * (size_t)m.n, sizeof(mp_limb_t));
    b = a + p.len * (size_t)m.n;
    s = calloc((size_t)3 * m.n + 1, sizeof(mp_limb_t));
    work = calloc((size_t)4 * m.n + 1, sizeof(mp_limb_t));
    for (size_t i = 0; i < 3; i++) {
        sc_fft_set_ui(sc_fft_entry(a, i, &m), av[i], &m);
        sc_fft_set_ui(sc_fft_entry(b, i, &m), bv[i], &m);
    }
    sc_fft_forward(a, &p, &m, work);
    sc_fft_forward(b, &p, &m, work);
    for (size_t i = 0; i < p.len; i++)
        sc_fft_mul(sc_fft_entry(a, i, &m), sc_fft_entry(a, i, &m),
                   sc_fft_entry(b, i, &m), &m, s);
    sc_fft_inverse(a, &p, &m, work);
    for (size_t i = 0; i < p.len; i++)
        ok &= sc_fft_equal_ui(sc_fft_entry(a, i, &m), cv[i], &m);
    free(work);
    free(s);
    free(a);
    sc_fft_plan_clear(&p);
    sc_fft_mod_clear(&m);
    return ok;
}

int main(void)
{
    if (!test_roundtrips(1, 64, 2, 7, 7) ||
        !test_roundtrips(1, 512, 2, 10, 10) ||
        !test_roundtrips(15, 27, 440564289, 27, 12) || !test_convolution()) {
        fprintf(stderr, "fft tests failed\n");
        return 1;
    }
    printf("fft tests passed\n");
    return 0;
}
