#include "smallcas.h"

sc_value *sc_zz_from_str(sc_context *ctx, const char *text)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    if (mpz_set_str(r->data.z, text, 10) != 0) {
        sc_value_free(r);
        sc_set_error(ctx, "invalid integer literal");
        return NULL;
    }
    return r;
}

sc_value *sc_zz_add(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_add(r->data.z, a->data.z, b->data.z);
    return r;
}

sc_value *sc_zz_sub(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_sub(r->data.z, a->data.z, b->data.z);
    return r;
}

sc_value *sc_zz_mul(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_mul(r->data.z, a->data.z, b->data.z);
    return r;
}

sc_value *sc_zz_divexact(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r;

    if (mpz_sgn(b->data.z) == 0) {
        sc_set_error(ctx, "division by zero");
        return NULL;
    }
    if (!mpz_divisible_p(a->data.z, b->data.z)) {
        sc_set_error(ctx, "quotient is not in ZZ");
        return NULL;
    }
    r = sc_value_new_zz_checked(ctx);
    if (r == NULL)
        return NULL;
    mpz_divexact(r->data.z, a->data.z, b->data.z);
    return r;
}

sc_value *sc_zz_mod(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r;

    if (mpz_sgn(b->data.z) == 0) {
        sc_set_error(ctx, "modulo by zero");
        return NULL;
    }
    r = sc_value_new_zz_checked(ctx);
    if (r == NULL)
        return NULL;
    mpz_mod(r->data.z, a->data.z, b->data.z);
    return r;
}

sc_value *sc_zz_neg(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_neg(r->data.z, a->data.z);
    return r;
}

sc_value *sc_zz_abs(sc_context *ctx, const sc_value *a)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_abs(r->data.z, a->data.z);
    return r;
}

sc_value *sc_zz_pow(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r;

    if (mpz_sgn(b->data.z) < 0 || !mpz_fits_ulong_p(b->data.z)) {
        sc_set_error(ctx, "exponent must be a nonnegative machine integer");
        return NULL;
    }
    r = sc_value_new_zz_checked(ctx);
    if (r == NULL)
        return NULL;
    mpz_pow_ui(r->data.z, a->data.z, mpz_get_ui(b->data.z));
    return r;
}

sc_value *sc_zz_gcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_gcd(r->data.z, a->data.z, b->data.z);
    return r;
}

sc_value *sc_zz_lcm(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *r = sc_value_new_zz_checked(ctx);

    if (r == NULL)
        return NULL;
    mpz_lcm(r->data.z, a->data.z, b->data.z);
    return r;
}
