#include "smallcas.h"


static int sc_poly_pair(sc_context *ctx, const sc_value *a, const sc_value *b,
                        const char *operation)
{
    if (a != NULL && b != NULL && a->kind == SC_VALUE_POLY &&
        b->kind == SC_VALUE_POLY && a->parent == b->parent)
        return 1;
    sc_set_error(ctx, "%s requires polynomials with a common parent", operation);
    return 0;
}

static int sc_poly_scalar_pair(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (a != NULL && b != NULL && a->kind == SC_VALUE_POLY &&
        b->kind != SC_VALUE_PARENT && b->parent == a->parent->base)
        return 1;
    sc_set_error(ctx, "scalar is not in the polynomial coefficient ring");
    return 0;
}

static int sc_poly_value(sc_context *ctx, const sc_value *a, const char *operation)
{
    if (a != NULL && a->kind == SC_VALUE_POLY)
        return 1;
    sc_set_error(ctx, "%s operand is not a polynomial", operation);
    return 0;
}

sc_value *sc_polynomial_ring(sc_context *ctx, sc_parent *base)
{
    sc_parent *parent;

    if (base != &SC_ZZ) {
        sc_set_error(ctx, "only PolynomialRing(ZZ) is supported");
        return NULL;
    }
    if (ctx->nparents == 32) {
        sc_set_error(ctx, "too many dynamic parents");
        return NULL;
    }
    parent = sc_parent_new_checked(ctx);
    if (parent == NULL)
        return NULL;
    parent->name = "PolynomialRing(ZZ)";
    parent->kind = SC_PARENT_POLY;
    parent->base = base;
    ctx->parents[ctx->nparents++] = parent;
    return sc_value_new_parent_checked(ctx, parent);
}

sc_value *sc_poly_generator(sc_context *ctx, sc_parent *parent, const char *name)
{
    if (parent == NULL || parent->kind != SC_PARENT_POLY) {
        sc_set_error(ctx, "parent is not a polynomial ring");
        return NULL;
    }
    if (parent->base == &SC_ZZ)
        return sc_zz_poly_generator(ctx, parent, name);
    sc_set_error(ctx, "no polynomial generator method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_from_base(sc_context *ctx, sc_parent *parent, const sc_value *a)
{
    if (parent == NULL || parent->kind != SC_PARENT_POLY || a == NULL ||
        a->kind == SC_VALUE_PARENT || a->parent != parent->base) {
        sc_set_error(ctx, "value is not in the polynomial coefficient ring");
        return NULL;
    }
    if (parent->base == &SC_ZZ && a->kind == SC_VALUE_ZZ)
        return sc_zz_poly_from_zz(ctx, parent, a);
    sc_set_error(ctx, "no polynomial coercion method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_add(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "addition"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_add(ctx, a, b);
    sc_set_error(ctx, "no polynomial addition method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_sub(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "subtraction"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_sub(ctx, a, b);
    sc_set_error(ctx, "no polynomial subtraction method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_mul(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "multiplication"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_mul(ctx, a, b);
    sc_set_error(ctx, "no polynomial multiplication method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_scalar_mul(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_scalar_mul(ctx, a, b);
    sc_set_error(ctx, "no scalar multiplication method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_neg(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "negation"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_neg(ctx, a);
    sc_set_error(ctx, "no polynomial negation method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_pow(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_value(ctx, a, "power") || b == NULL || b->kind != SC_VALUE_ZZ) {
        sc_set_error(ctx, "polynomial power requires an integer exponent");
        return NULL;
    }
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_pow(ctx, a, b);
    sc_set_error(ctx, "no polynomial power method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_divexact(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "division"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_divexact(ctx, a, b);
    sc_set_error(ctx, "no exact polynomial division method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_quo(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "quo"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_quo(ctx, a, b);
    sc_set_error(ctx, "no polynomial quotient method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_divrem(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "divrem"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_divrem(ctx, a, b);
    sc_set_error(ctx, "no polynomial divrem method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_pseudodiv(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "pseudodiv"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_pseudodiv(ctx, a, b);
    sc_set_error(ctx, "no polynomial pseudodiv method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_scalar_divexact(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_scalar_divexact(ctx, a, b);
    sc_set_error(ctx, "no scalar division method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_content(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "content"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_content(ctx, a);
    sc_set_error(ctx, "no content method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_primitive_part(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "primitive_part"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_primitive_part(ctx, a);
    sc_set_error(ctx, "no primitive-part method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_gcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "gcd"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_gcd(ctx, a, b);
    sc_set_error(ctx, "no polynomial gcd method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_gcd_pseudo(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "gcd_pseudo"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_gcd_pseudo(ctx, a, b);
    sc_set_error(ctx, "no polynomial pseudo-gcd method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_gcd_hgcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "gcd_hgcd"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_gcd_hgcd(ctx, a, b);
    sc_set_error(ctx, "no polynomial HGCD method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_gcd_lr(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "gcd_lr"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_gcd_lr(ctx, a, b);
    sc_set_error(ctx, "no LR gcd method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_resultant(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "resultant"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_resultant(ctx, a, b);
    sc_set_error(ctx, "no polynomial resultant method for this coefficient ring");
    return NULL;
}


sc_value *sc_poly_xgcd(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "xgcd"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_xgcd(ctx, a, b);
    sc_set_error(ctx, "no polynomial xgcd method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_evaluate(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_evaluate(ctx, a, b);
    sc_set_error(ctx, "no evaluation method for this polynomial coefficient ring");
    return NULL;
}

sc_value *sc_poly_evaluate_horner(sc_context *ctx, const sc_value *a,
                                  const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_evaluate_horner(ctx, a, b);
    sc_set_error(ctx, "no Horner evaluation method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_evaluate_divconquer(sc_context *ctx, const sc_value *a,
                                      const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_evaluate_divconquer(ctx, a, b);
    sc_set_error(ctx, "no D&C evaluation method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_compose(sc_context *ctx, const sc_value *a, const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "composition"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_compose(ctx, a, b);
    sc_set_error(ctx, "no composition method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_compose_horner(sc_context *ctx, const sc_value *a,
                                 const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "composition"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_compose_horner(ctx, a, b);
    sc_set_error(ctx, "no Horner composition method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_compose_divconquer(sc_context *ctx, const sc_value *a,
                                     const sc_value *b)
{
    if (!sc_poly_pair(ctx, a, b, "composition"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_compose_divconquer(ctx, a, b);
    sc_set_error(ctx, "no D&C composition method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_taylor_shift(sc_context *ctx, const sc_value *a,
                               const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_taylor_shift(ctx, a, b);
    sc_set_error(ctx, "no Taylor-shift method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_taylor_shift_horner(sc_context *ctx, const sc_value *a,
                                      const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_taylor_shift_horner(ctx, a, b);
    sc_set_error(ctx, "no Horner Taylor-shift method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_taylor_shift_divconquer(sc_context *ctx, const sc_value *a,
                                          const sc_value *b)
{
    if (!sc_poly_scalar_pair(ctx, a, b))
        return NULL;
    if (a->parent->base == &SC_ZZ && b->kind == SC_VALUE_ZZ)
        return sc_zz_poly_taylor_shift_divconquer(ctx, a, b);
    sc_set_error(ctx, "no D&C Taylor-shift method for this coefficient ring");
    return NULL;
}


sc_value *sc_poly_derivative(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "derivative"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_derivative(ctx, a);
    sc_set_error(ctx, "no derivative method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_nth_derivative(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "derivative"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_nth_derivative(ctx, a, n);
    sc_set_error(ctx, "no nth derivative method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_discriminant(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "discriminant"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_discriminant(ctx, a);
    sc_set_error(ctx, "no discriminant method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_is_squarefree(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "is_squarefree"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_is_squarefree(ctx, a);
    sc_set_error(ctx, "no squarefree test for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_squarefree_part(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "squarefree_part"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_squarefree_part(ctx, a);
    sc_set_error(ctx, "no squarefree-part method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_degree(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "degree"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_degree(ctx, a);
    sc_set_error(ctx, "no degree method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_leading_coefficient(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "leading_coefficient"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_leading_coefficient(ctx, a);
    sc_set_error(ctx, "no leading-coefficient method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_constant_coefficient(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "constant_coefficient"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_constant_coefficient(ctx, a);
    sc_set_error(ctx, "no constant-coefficient method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_coeff(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "coeff"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_coeff(ctx, a, n);
    sc_set_error(ctx, "no coefficient-access method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_reverse(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "reverse"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_reverse(ctx, a, n);
    sc_set_error(ctx, "no reverse method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_truncate(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "truncate"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_truncate_copy(ctx, a, n);
    sc_set_error(ctx, "no truncation method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_shift_left(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "shift_left"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_shift_left(ctx, a, n);
    sc_set_error(ctx, "no left-shift method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_shift_right(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "shift_right"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_shift_right(ctx, a, n);
    sc_set_error(ctx, "no right-shift method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_height(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "height"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_height(ctx, a);
    sc_set_error(ctx, "no height method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_max_abs_bits(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "max_abs_bits"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_max_abs_bits(ctx, a);
    sc_set_error(ctx, "no coefficient bit-size method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_inflate(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "inflate"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_inflate(ctx, a, n);
    sc_set_error(ctx, "no inflation method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_deflation(sc_context *ctx, const sc_value *a)
{
    if (!sc_poly_value(ctx, a, "deflation"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_deflation(ctx, a);
    sc_set_error(ctx, "no deflation method for this coefficient ring");
    return NULL;
}

sc_value *sc_poly_deflate(sc_context *ctx, const sc_value *a, size_t n)
{
    if (!sc_poly_value(ctx, a, "deflate"))
        return NULL;
    if (a->parent->base == &SC_ZZ)
        return sc_zz_poly_deflate(ctx, a, n);
    sc_set_error(ctx, "no deflate method for this coefficient ring");
    return NULL;
}
