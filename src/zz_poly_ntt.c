#include "smallcas.h"
#include "smallcas_fft.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define SC_NTT_DEPTH 24
#define SC_NTT_PRIME_BITS 60

typedef struct {
    uint64_t c, root;
} sc_ntt_prime;

static const sc_ntt_prime sc_ntt_primes[] = {
    { UINT64_C(68719476771), UINT64_C(792893038599829801) },
    { UINT64_C(68719476775), UINT64_C(91800153599965448) },
    { UINT64_C(68719476793), UINT64_C(274471583722122878) },
    { UINT64_C(68719476813), UINT64_C(269433890925539316) },
    { UINT64_C(68719476831), UINT64_C(426217368599154363) },
    { UINT64_C(68719476835), UINT64_C(1147363667100351485) },
    { UINT64_C(68719476943), UINT64_C(659064345051291461) },
    { UINT64_C(68719476961), UINT64_C(557327111208626070) },
    { UINT64_C(68719476963), UINT64_C(268223570106781950) },
    { UINT64_C(68719477057), UINT64_C(921751764709828673) },
    { UINT64_C(68719477083), UINT64_C(1100673306306755326) },
    { UINT64_C(68719477167), UINT64_C(525652500323346389) },
    { UINT64_C(68719477185), UINT64_C(536204626827999747) },
    { UINT64_C(68719477225), UINT64_C(741734686530937815) },
    { UINT64_C(68719477395), UINT64_C(205184732444534124) },
    { UINT64_C(68719477417), UINT64_C(798929101234008049) },
    { UINT64_C(68719477425), UINT64_C(492848167428983656) },
    { UINT64_C(68719477447), UINT64_C(800278120038747653) },
    { UINT64_C(68719477485), UINT64_C(343868008004814888) },
    { UINT64_C(68719477563), UINT64_C(240469417569708092) },
    { UINT64_C(68719477713), UINT64_C(1131754979111235552) },
    { UINT64_C(68719477753), UINT64_C(174468147219459643) },
    { UINT64_C(68719477765), UINT64_C(165783880346690137) },
    { UINT64_C(68719477831), UINT64_C(791899134593205249) },
    { UINT64_C(68719477857), UINT64_C(249978664075842080) },
    { UINT64_C(68719477863), UINT64_C(504539099892706417) },
    { UINT64_C(68719477887), UINT64_C(991637250584924236) },
    { UINT64_C(68719477975), UINT64_C(396354790989170163) },
    { UINT64_C(68719477981), UINT64_C(1041509218677101781) },
    { UINT64_C(68719478077), UINT64_C(1045369832772303646) },
    { UINT64_C(68719478095), UINT64_C(631341474070207279) },
    { UINT64_C(68719478121), UINT64_C(770637029477773296) }
};

static size_t sc_ntt_log2ceil(size_t n)
{
    size_t k = 0;

    if (n != 0)
        n--;
    while (n != 0)
        k++, n >>= 1;
    return k;
}

static int sc_ntt_params(size_t *np, unsigned *logn,
                         const sc_value *a, const sc_value *b)
{
    size_t an = a->data.zz_poly.length, bn = b->data.zz_poly.length;
    size_t small = an < bn ? an : bn, need;
    size_t ba = sc_zz_poly_max_abs_bits_raw(a);
    size_t bb = sc_zz_poly_max_abs_bits_raw(b);

    if (an > SIZE_MAX - bn + 1 || ba > SIZE_MAX - bb)
        return 0;
    *logn = (unsigned)sc_ntt_log2ceil(an + bn - 1);
    if (*logn > SC_NTT_DEPTH || ba + bb > SIZE_MAX - sc_ntt_log2ceil(small) - 1)
        return 0;
    need = ba + bb + sc_ntt_log2ceil(small) + 1;
    *np = (need + SC_NTT_PRIME_BITS - 1) / SC_NTT_PRIME_BITS;
    return *np <= sizeof(sc_ntt_primes) / sizeof(sc_ntt_primes[0]);
}

size_t sc_zz_poly_ntt_nprimes(const sc_value *a, const sc_value *b)
{
    size_t np;
    unsigned logn;

    return sc_ntt_params(&np, &logn, a, b) ? np : 0;
}

static mp_limb_t sc_ntt_mulmod(mp_limb_t a, mp_limb_t b, mp_limb_t p)
{
    mp_limb_t t[2];

    t[1] = mpn_mul_1(t, &a, 1, b);
    return mpn_mod_1(t, 2, p);
}

static mp_limb_t sc_ntt_inverse(mpz_t t, mpz_t u, mpz_srcptr M, mp_limb_t p)
{
    mp_limb_t a = (mp_limb_t)mpz_fdiv_ui(M, (unsigned long)p);

    mpz_set_ui(t, (unsigned long)a);
    mpz_set_ui(u, (unsigned long)p);
    return mpz_invert(t, t, u) ? (mp_limb_t)mpz_get_ui(t) : 0;
}

static void sc_ntt_crt(sc_value *r, mp_srcptr v, size_t n, mpz_t M,
                       mp_limb_t p, mp_limb_t inv)
{
    for (size_t i = 0; i < n; i++) {
        mp_limb_t x = (mp_limb_t)mpz_fdiv_ui(r->data.zz_poly.coeff[i], (unsigned long)p);
        mp_limb_t d = v[i] >= x ? v[i] - x : p - (x - v[i]);
        mp_limb_t q = sc_ntt_mulmod(d, inv, p);

        mpz_addmul_ui(r->data.zz_poly.coeff[i], M, (unsigned long)q);
    }
}

static int sc_ntt_prime_pass(sc_value *r, const sc_value *a, const sc_value *b,
                             unsigned logn, size_t pi, mpz_t M, mpz_t t, mpz_t u,
                             mp_ptr v, mp_ptr work)
{
    const sc_ntt_prime *q = sc_ntt_primes + pi;
    sc_fft_mod m;
    sc_fft_plan plan;
    size_t n = (size_t)1 << logn;
    mp_limb_t p, inv;
    mp_ptr w = v + n;

    if (!sc_fft_mod_init(&m, q->c, SC_NTT_DEPTH, q->root, SC_NTT_DEPTH))
        return 0;
    if (m.n != 1 || !sc_fft_plan_init(&plan, logn, &m)) {
        sc_fft_mod_clear(&m);
        return 0;
    }
    p = m.mod[0];
    inv = sc_ntt_inverse(t, u, M, p);
    if (inv == 0) {
        sc_fft_plan_clear(&plan);
        sc_fft_mod_clear(&m);
        return 0;
    }
    memset(v, 0, 2 * n * sizeof(mp_limb_t));
    for (size_t i = 0; i < a->data.zz_poly.length; i++)
        v[i] = (mp_limb_t)mpz_fdiv_ui(a->data.zz_poly.coeff[i], (unsigned long)p);
    for (size_t i = 0; i < b->data.zz_poly.length; i++)
        w[i] = (mp_limb_t)mpz_fdiv_ui(b->data.zz_poly.coeff[i], (unsigned long)p);
    sc_fft_forward(v, &plan, &m, work);
    sc_fft_forward(w, &plan, &m, work);
    for (size_t i = 0; i < n; i++)
        sc_fft_mul(v + i, v + i, w + i, &m, work);
    sc_fft_inverse(v, &plan, &m, work);
    sc_ntt_crt(r, v, r->data.zz_poly.length, M, p, inv);
    mpz_mul_ui(M, M, (unsigned long)p);
    sc_fft_plan_clear(&plan);
    sc_fft_mod_clear(&m);
    return 1;
}

sc_value *sc_zz_poly_mul_ntt(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    size_t an, bn, outn, np, n;
    unsigned logn;
    sc_value *r = NULL;
    mp_ptr v = NULL, work = NULL;
    mpz_t M, half, t, u;

    if (a == NULL || b == NULL)
        return NULL;
    an = a->data.zz_poly.length;
    bn = b->data.zz_poly.length;
    if (an == 0 || bn == 0)
        return sc_value_new_zz_poly_checked(ctx, a->parent, 0);
    if (GMP_NUMB_BITS != 64 || ULONG_MAX < UINT64_MAX ||
        !sc_ntt_params(&np, &logn, a, b)) {
        sc_set_error(ctx, "NTT parameters unsupported");
        return NULL;
    }
    outn = an + bn - 1;
    n = (size_t)1 << logn;
    r = sc_value_new_zz_poly_checked(ctx, a->parent, outn);
    v = calloc(2 * n, sizeof(mp_limb_t));
    work = calloc(5, sizeof(mp_limb_t));
    if (r == NULL || v == NULL || work == NULL)
        goto fail;
    mpz_inits(M, half, t, u, NULL);
    mpz_set_ui(M, 1);
    for (size_t i = 0; i < np; i++)
        if (!sc_ntt_prime_pass(r, a, b, logn, i, M, t, u, v, work))
            goto fail_mpz;
    mpz_fdiv_q_2exp(half, M, 1);
    for (size_t i = 0; i < outn; i++)
        if (mpz_cmp(r->data.zz_poly.coeff[i], half) > 0)
            mpz_sub(r->data.zz_poly.coeff[i], r->data.zz_poly.coeff[i], M);
    mpz_clears(M, half, t, u, NULL);
    free(work);
    free(v);
    sc_zz_poly_normalize(r);
    return r;
fail_mpz:
    mpz_clears(M, half, t, u, NULL);
fail:
    free(work);
    free(v);
    sc_value_free(r);
    if (ctx->error[0] == '\0')
        sc_set_error(ctx, "out of memory in NTT multiplication");
    return NULL;
}
