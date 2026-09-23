#include "smallcas.h"

#include <string.h>

static int sc_is_zz(const sc_value *a)
{
    return a != NULL && a->kind == SC_VALUE_ZZ;
}

static int sc_is_poly(const sc_value *a)
{
    return a != NULL && a->kind == SC_VALUE_POLY;
}

static int sc_get_size(const sc_value *a, size_t *n)
{
    if (!sc_is_zz(a) || mpz_sgn(a->data.z) < 0 || !mpz_fits_ulong_p(a->data.z))
        return 0;
    *n = (size_t)mpz_get_ui(a->data.z);
    return 1;
}

static sc_value *sc_pair_take(sc_value *pair, int first)
{
    sc_value *value;

    if (pair == NULL)
        return NULL;
    value = first ? pair->data.pair.first : pair->data.pair.second;
    if (first)
        pair->data.pair.first = NULL;
    else
        pair->data.pair.second = NULL;
    sc_value_free(pair);
    return value;
}

sc_value *sc_coerce(sc_context *ctx, sc_parent *parent, const sc_value *a)
{
    if (a == NULL)
        return NULL;
    if (a->kind != SC_VALUE_PARENT && a->parent == parent)
        return sc_value_copy_checked(ctx, a);
    if (parent == &SC_ZZ && sc_is_zz(a))
        return sc_value_copy_checked(ctx, a);
    if (parent->kind == SC_PARENT_POLY && a->kind != SC_VALUE_PARENT &&
        a->parent == parent->base)
        return sc_poly_from_base(ctx, parent, a);
    sc_set_error(ctx, "cannot coerce value into %s", parent->name);
    return NULL;
}

sc_value *sc_add(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *c, *r;

    if (sc_is_zz(a) && sc_is_zz(b))
        return sc_zz_add(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent)
        return sc_poly_add(ctx, a, b);
    if (sc_is_poly(a) && b != NULL && b->kind != SC_VALUE_PARENT &&
        b->parent == a->parent->base) {
        c = sc_coerce(ctx, a->parent, b);
        r = c ? sc_poly_add(ctx, a, c) : NULL;
        sc_value_free(c);
        return r;
    }
    if (sc_is_poly(b) && a != NULL && a->kind != SC_VALUE_PARENT &&
        a->parent == b->parent->base) {
        c = sc_coerce(ctx, b->parent, a);
        r = c ? sc_poly_add(ctx, c, b) : NULL;
        sc_value_free(c);
        return r;
    }
    sc_set_error(ctx, "no addition method for these parents");
    return NULL;
}

sc_value *sc_sub(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    sc_value *c, *r;

    if (sc_is_zz(a) && sc_is_zz(b))
        return sc_zz_sub(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent)
        return sc_poly_sub(ctx, a, b);
    if (sc_is_poly(a) && b != NULL && b->kind != SC_VALUE_PARENT &&
        b->parent == a->parent->base) {
        c = sc_coerce(ctx, a->parent, b);
        r = c ? sc_poly_sub(ctx, a, c) : NULL;
        sc_value_free(c);
        return r;
    }
    if (sc_is_poly(b) && a != NULL && a->kind != SC_VALUE_PARENT &&
        a->parent == b->parent->base) {
        c = sc_coerce(ctx, b->parent, a);
        r = c ? sc_poly_sub(ctx, c, b) : NULL;
        sc_value_free(c);
        return r;
    }
    sc_set_error(ctx, "no subtraction method for these parents");
    return NULL;
}

sc_value *sc_mul(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (sc_is_zz(a) && sc_is_zz(b))
        return sc_zz_mul(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent)
        return sc_poly_mul(ctx, a, b);
    if (sc_is_poly(a) && b != NULL && b->kind != SC_VALUE_PARENT &&
        b->parent == a->parent->base)
        return sc_poly_scalar_mul(ctx, a, b);
    if (sc_is_poly(b) && a != NULL && a->kind != SC_VALUE_PARENT &&
        a->parent == b->parent->base)
        return sc_poly_scalar_mul(ctx, b, a);
    sc_set_error(ctx, "no multiplication method for these parents");
    return NULL;
}

sc_value *sc_div(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (sc_is_zz(a) && sc_is_zz(b))
        return sc_zz_divexact(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent)
        return sc_poly_divexact(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && b->parent == a->parent->base)
        return sc_poly_scalar_divexact(ctx, a, b);
    sc_set_error(ctx, "no division method for these parents");
    return NULL;
}

sc_value *sc_mod(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (sc_is_zz(a) && sc_is_zz(b))
        return sc_zz_mod(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent)
        return sc_pair_take(sc_poly_divrem(ctx, a, b), 0);
    sc_set_error(ctx, "no modulo method for these parents");
    return NULL;
}

sc_value *sc_neg(sc_context *ctx, const sc_value *a)
{
    if (sc_is_zz(a))
        return sc_zz_neg(ctx, a);
    if (sc_is_poly(a))
        return sc_poly_neg(ctx, a);
    sc_set_error(ctx, "no negation method for this parent");
    return NULL;
}

sc_value *sc_pow(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (sc_is_zz(a) && sc_is_zz(b))
        return sc_zz_pow(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b))
        return sc_poly_pow(ctx, a, b);
    sc_set_error(ctx, "no power method for these parents");
    return NULL;
}

sc_value *sc_call1(sc_context *ctx, const char *name, const sc_value *a)
{
    sc_value *callee, *r;

    if (strcmp(name, "PolynomialRing") == 0 && a != NULL &&
        a->kind == SC_VALUE_PARENT)
        return sc_polynomial_ring(ctx, a->data.parent_value);
    if (strcmp(name, "abs") == 0 && sc_is_zz(a))
        return sc_zz_abs(ctx, a);
    if (strcmp(name, "content") == 0 && sc_is_poly(a))
        return sc_poly_content(ctx, a);
    if (strcmp(name, "primitive_part") == 0 && sc_is_poly(a))
        return sc_poly_primitive_part(ctx, a);
    if (strcmp(name, "derivative") == 0 && sc_is_poly(a))
        return sc_poly_derivative(ctx, a);
    if (strcmp(name, "discriminant") == 0 && sc_is_poly(a))
        return sc_poly_discriminant(ctx, a);
    if (strcmp(name, "is_squarefree") == 0 && sc_is_poly(a))
        return sc_poly_is_squarefree(ctx, a);
    if (strcmp(name, "squarefree_part") == 0 && sc_is_poly(a))
        return sc_poly_squarefree_part(ctx, a);
    if (strcmp(name, "degree") == 0 && sc_is_poly(a))
        return sc_poly_degree(ctx, a);
    if (strcmp(name, "leading_coefficient") == 0 && sc_is_poly(a))
        return sc_poly_leading_coefficient(ctx, a);
    if (strcmp(name, "constant_coefficient") == 0 && sc_is_poly(a))
        return sc_poly_constant_coefficient(ctx, a);
    if (strcmp(name, "reverse") == 0 && sc_is_poly(a))
        return sc_poly_reverse(ctx, a, a->data.zz_poly.length);
    if (strcmp(name, "height") == 0 && sc_is_poly(a))
        return sc_poly_height(ctx, a);
    if (strcmp(name, "max_abs_bits") == 0 && sc_is_poly(a))
        return sc_poly_max_abs_bits(ctx, a);
    if ((strcmp(name, "deflation") == 0 || strcmp(name, "maximal_deflation") == 0) &&
        sc_is_poly(a))
        return sc_poly_deflation(ctx, a);
    if (strcmp(name, "ZZ") == 0)
        return sc_coerce(ctx, &SC_ZZ, a);
    callee = sc_env_get(ctx, name);
    if (callee == NULL)
        return NULL;
    if (callee->kind == SC_VALUE_POLY) {
        r = a != NULL && a->kind == SC_VALUE_POLY ?
            sc_poly_compose(ctx, callee, a) : sc_poly_evaluate(ctx, callee, a);
        sc_value_free(callee);
        return r;
    }
    if (callee->kind != SC_VALUE_PARENT) {
        sc_value_free(callee);
        sc_set_error(ctx, "'%s' is not callable", name);
        return NULL;
    }
    r = sc_coerce(ctx, callee->data.parent_value, a);
    sc_value_free(callee);
    return r;
}

sc_value *sc_call_string(sc_context *ctx, const char *name, const char *text)
{
    sc_value *callee, *r;

    callee = sc_env_get(ctx, name);
    if (callee == NULL)
        return NULL;
    if (callee->kind != SC_VALUE_PARENT) {
        sc_value_free(callee);
        sc_set_error(ctx, "'%s' is not a parent", name);
        return NULL;
    }
    r = sc_poly_generator(ctx, callee->data.parent_value, text);
    sc_value_free(callee);
    return r;
}

sc_value *sc_call2(sc_context *ctx, const char *name,
                   const sc_value *a, const sc_value *b)
{
    if (sc_is_zz(a) && sc_is_zz(b) && strcmp(name, "gcd") == 0)
        return sc_zz_gcd(ctx, a, b);
    if (sc_is_zz(a) && sc_is_zz(b) && strcmp(name, "lcm") == 0)
        return sc_zz_lcm(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "gcd") == 0)
        return sc_poly_gcd(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "gcd_pseudo") == 0)
        return sc_poly_gcd_pseudo(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "gcd_hgcd") == 0)
        return sc_poly_gcd_hgcd(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "gcd_lr") == 0)
        return sc_poly_gcd_lr(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "resultant") == 0)
        return sc_poly_resultant(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "xgcd") == 0)
        return sc_poly_xgcd(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "quo") == 0)
        return sc_poly_quo(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "divrem") == 0)
        return sc_poly_divrem(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "pseudodiv") == 0)
        return sc_poly_pseudodiv(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) &&
        (strcmp(name, "derivative") == 0 || strcmp(name, "nth_derivative") == 0)) {
        size_t n;

        if (!sc_get_size(b, &n)) {
            sc_set_error(ctx, "derivative order must be a nonnegative machine integer");
            return NULL;
        }
        return sc_poly_nth_derivative(ctx, a, n);
    }
    if (sc_is_poly(a) && sc_is_zz(b) &&
        (strcmp(name, "coeff") == 0 || strcmp(name, "reverse") == 0 ||
         strcmp(name, "truncate") == 0 || strcmp(name, "shift_left") == 0 ||
         strcmp(name, "shift_right") == 0 || strcmp(name, "inflate") == 0 ||
         strcmp(name, "deflate") == 0)) {
        size_t n;

        if (!sc_get_size(b, &n)) {
            sc_set_error(ctx, "polynomial index/factor must be a nonnegative machine integer");
            return NULL;
        }
        if (strcmp(name, "coeff") == 0)
            return sc_poly_coeff(ctx, a, n);
        if (strcmp(name, "reverse") == 0)
            return sc_poly_reverse(ctx, a, n);
        if (strcmp(name, "truncate") == 0)
            return sc_poly_truncate(ctx, a, n);
        if (strcmp(name, "shift_left") == 0)
            return sc_poly_shift_left(ctx, a, n);
        if (strcmp(name, "shift_right") == 0)
            return sc_poly_shift_right(ctx, a, n);
        if (strcmp(name, "inflate") == 0)
            return sc_poly_inflate(ctx, a, n);
        return sc_poly_deflate(ctx, a, n);
    }
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "evaluate") == 0)
        return sc_poly_evaluate(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "evaluate_horner") == 0)
        return sc_poly_evaluate_horner(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "evaluate_dc") == 0)
        return sc_poly_evaluate_divconquer(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "taylor_shift") == 0)
        return sc_poly_taylor_shift(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "taylor_shift_horner") == 0)
        return sc_poly_taylor_shift_horner(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "taylor_shift_dc") == 0)
        return sc_poly_taylor_shift_divconquer(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "compose") == 0)
        return sc_poly_compose(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "compose_horner") == 0)
        return sc_poly_compose_horner(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "compose_dc") == 0)
        return sc_poly_compose_divconquer(ctx, a, b);
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "inv_series") == 0) {
        size_t n;

        if (!sc_get_size(b, &n)) {
            sc_set_error(ctx, "series precision must be a nonnegative machine integer");
            return NULL;
        }
        return sc_zz_poly_inv_series(ctx, a, n);
    }
    if (sc_is_poly(a) && sc_is_zz(b) && strcmp(name, "preinverse") == 0) {
        size_t n;

        if (!sc_get_size(b, &n)) {
            sc_set_error(ctx, "preinverse precision must be a nonnegative machine integer");
            return NULL;
        }
        return sc_zz_poly_preinverse_newton(ctx, a, n);
    }
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "quo_newton") == 0)
        return sc_zz_poly_quo_newton(ctx, a, b);
    if (sc_is_poly(a) && sc_is_poly(b) && a->parent == b->parent &&
        strcmp(name, "divrem_newton") == 0)
        return sc_zz_poly_divrem_newton(ctx, a, b);
    if (strcmp(name, "divexact") == 0)
        return sc_div(ctx, a, b);
    if (strcmp(name, "pow") == 0)
        return sc_pow(ctx, a, b);
    sc_set_error(ctx, "unknown binary function '%s'", name);
    return NULL;
}


sc_value *sc_call3(sc_context *ctx, const char *name,
                   const sc_value *a, const sc_value *b, const sc_value *c)
{
    if (strcmp(name, "quo_preinv") == 0 && sc_is_poly(a) && sc_is_poly(b) &&
        sc_is_poly(c) && a->parent == b->parent && a->parent == c->parent)
        return sc_zz_poly_quo_preinv(ctx, a, b, c);
    sc_set_error(ctx, "unknown ternary function '%s'", name);
    return NULL;
}
