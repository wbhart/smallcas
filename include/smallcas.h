#ifndef SMALLCAS_H
#define SMALLCAS_H

#include <gmp.h>
#include <stddef.h>

#include "tuning.h"

typedef enum {
    SC_PARENT_ZZ = 1,
    SC_PARENT_POLY
} sc_parent_kind;

typedef enum {
    SC_VALUE_PARENT = 1,
    SC_VALUE_ZZ,
    SC_VALUE_POLY,
    SC_VALUE_PAIR
} sc_value_kind;

typedef enum {
    SC_PATTERN_NAME = 1,
    SC_PATTERN_PAIR
} sc_pattern_kind;

typedef struct sc_parent {
    const char *name;
    sc_parent_kind kind;
    struct sc_parent *base;
    char *symbol;
} sc_parent;

typedef struct sc_zz_poly {
    size_t length;
    mpz_t *coeff;
} sc_zz_poly;

typedef struct sc_value {
    sc_value_kind kind;
    sc_parent *parent;
    union {
        mpz_t z;
        sc_zz_poly zz_poly;
        sc_parent *parent_value;
        struct {
            struct sc_value *first;
            struct sc_value *second;
        } pair;
    } data;
} sc_value;

typedef struct sc_zz_poly_toom3_ws {
    size_t m;
    size_t q;
    size_t length;
    sc_value block[2][3];
    sc_value *point[2][5];
    sc_value *prod[5];
    sc_value *result;
    mpz_t zero;
} sc_zz_poly_toom3_ws;

typedef struct sc_zz_poly_toom63_ws {
    size_t k;
    sc_value fblock[5];
    sc_value gblock[3];
    sc_value *fpoint[5];
    sc_value *gpoint[5];
    sc_value *prod[5];
    sc_value *result;
    mpz_t zero;
} sc_zz_poly_toom63_ws;

typedef struct sc_zz_poly_subres_ws {
    sc_value *u;
    sc_value *v;
    sc_value *w;
    sc_value *z;
    mpz_t hp;
    mpz_t hc;
    mpz_t hn;
    mpz_t den;
} sc_zz_poly_subres_ws;


typedef struct sc_zz_poly_mat2 {
    sc_value *a00;
    sc_value *a01;
    sc_value *a10;
    sc_value *a11;
} sc_zz_poly_mat2;

typedef struct sc_zz_poly_xgcd_ws {
    sc_value *u;
    sc_value *v;
    sc_value *w;
    sc_value *cu[2];
    sc_value *cv[2];
    sc_value *cw[2];
    mpz_t hp;
    mpz_t hc;
    mpz_t hn;
    mpz_t den;
    mpz_t alpha;
} sc_zz_poly_xgcd_ws;

typedef struct sc_pattern {
    sc_pattern_kind kind;
    union {
        char *name;
        struct {
            struct sc_pattern *first;
            struct sc_pattern *second;
        } pair;
    } data;
} sc_pattern;

typedef struct sc_binding {
    char *name;
    sc_value *value;
} sc_binding;

typedef struct sc_context {
    sc_binding vars[128];
    size_t nvars;
    sc_parent *parents[32];
    size_t nparents;
    sc_value *result;
    int show_result;
    int timing;
    char error[256];
} sc_context;

extern sc_parent SC_ZZ;
extern sc_context *sc_parse_context;

void sc_context_init(sc_context *ctx);
void sc_context_clear(sc_context *ctx);
void sc_set_error(sc_context *ctx, const char *fmt, ...);

char *sc_string_dup(sc_context *ctx, const char *text);
char *sc_string_ndup(sc_context *ctx, const char *text, size_t n);
sc_parent *sc_parent_new_checked(sc_context *ctx);

sc_value *sc_value_new_zz(void);
sc_value *sc_value_new_zz_checked(sc_context *ctx);
sc_value *sc_value_new_zz_poly(sc_parent *parent, size_t length);
sc_value *sc_value_new_zz_poly_checked(sc_context *ctx, sc_parent *parent,
                                       size_t length);
sc_value *sc_value_new_parent(sc_parent *parent);
sc_value *sc_value_new_parent_checked(sc_context *ctx, sc_parent *parent);
sc_value *sc_value_new_pair_take_checked(sc_context *ctx, sc_value *first,
                                         sc_value *second);
sc_value *sc_value_copy(const sc_value *value);
sc_value *sc_value_copy_checked(sc_context *ctx, const sc_value *value);
void sc_value_free(sc_value *value);
void sc_value_free_null(sc_value **value);
void sc_value_ptr_swap(sc_value **a, sc_value **b);
void sc_value_free_many(size_t count, ...);
sc_value **sc_value_array_new_checked(sc_context *ctx, size_t count);
void sc_value_array_free(size_t count, sc_value **values);
sc_value *sc_value_free_many_null(size_t count, ...);
sc_value *sc_value_pair_take(sc_value *pair, int second);
int sc_value_pair_split(sc_value *pair, sc_value **first, sc_value **second);
sc_pattern *sc_pattern_new_name_take_checked(sc_context *ctx, char *name);
sc_pattern *sc_pattern_new_pair_take_checked(sc_context *ctx, sc_pattern *first,
                                             sc_pattern *second);
void sc_pattern_free(sc_pattern *pattern);
void sc_zz_poly_truncate(sc_value *value, size_t length);
void sc_zz_poly_normalize(sc_value *value);
sc_value sc_zz_poly_view(const sc_value *value, size_t start, size_t length);
int sc_zz_poly_toom3_ws_init(sc_context *ctx, sc_zz_poly_toom3_ws *ws,
                             const sc_value *a, const sc_value *b, size_t m);
sc_value *sc_zz_poly_toom3_ws_abort(sc_zz_poly_toom3_ws *ws);
sc_value *sc_zz_poly_toom3_ws_finish(sc_zz_poly_toom3_ws *ws);
int sc_zz_poly_toom63_ws_init(sc_context *ctx, sc_zz_poly_toom63_ws *ws,
                              const sc_value *a, const sc_value *b, size_t k);
sc_value *sc_zz_poly_toom63_ws_abort(sc_zz_poly_toom63_ws *ws);
sc_value *sc_zz_poly_toom63_ws_finish(sc_zz_poly_toom63_ws *ws);
int sc_zz_poly_subres_ws_init(sc_context *ctx, sc_zz_poly_subres_ws *ws,
                              const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_subres_ws_finish(sc_context *ctx, sc_zz_poly_subres_ws *ws,
                                      sc_value *keep, mpz_srcptr h);
void sc_zz_poly_subres_ws_shift(sc_zz_poly_subres_ws *ws);
int sc_zz_poly_xgcd_ws_init(sc_context *ctx, sc_zz_poly_xgcd_ws *ws,
                            const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_xgcd_ws_finish(sc_context *ctx, sc_zz_poly_xgcd_ws *ws);
sc_value *sc_zz_poly_xgcd_ws_abort(sc_zz_poly_xgcd_ws *ws);
void sc_zz_poly_xgcd_ws_shift(sc_zz_poly_xgcd_ws *ws);
int sc_zz_poly_mat2_identity(sc_context *ctx, sc_zz_poly_mat2 *m, sc_parent *parent);
void sc_zz_poly_mat2_clear(sc_zz_poly_mat2 *m);
void sc_zz_poly_mat2_clear_many(size_t count, ...);
void sc_zz_poly_mat2_move(sc_zz_poly_mat2 *dst, sc_zz_poly_mat2 *src);
void sc_value_print(const sc_value *value);

int sc_env_set(sc_context *ctx, const char *name, const sc_value *value);
sc_value *sc_env_get(sc_context *ctx, const char *name);
int sc_env_set_pattern(sc_context *ctx, const sc_pattern *pattern,
                       const sc_value *value);

sc_value *sc_zz_from_str(sc_context *ctx, const char *text);
sc_value *sc_zz_add(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_sub(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_mul(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_divexact(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_mod(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_neg(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_abs(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_pow(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_gcd(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_lcm(sc_context *ctx, const sc_value *a, const sc_value *b);

sc_value *sc_polynomial_ring(sc_context *ctx, sc_parent *base);
sc_value *sc_poly_generator(sc_context *ctx, sc_parent *parent, const char *name);
sc_value *sc_poly_from_base(sc_context *ctx, sc_parent *parent, const sc_value *a);
sc_value *sc_poly_add(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_sub(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_mul(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_scalar_mul(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_neg(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_pow(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_divexact(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_quo(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_divrem(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_pseudodiv(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_scalar_divexact(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_content(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_primitive_part(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_gcd(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_gcd_pseudo(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_gcd_hgcd(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_gcd_lr(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_resultant(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_xgcd(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_evaluate(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_evaluate_horner(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_evaluate_divconquer(sc_context *ctx, const sc_value *a,
                                      const sc_value *b);
sc_value *sc_poly_compose(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_poly_compose_horner(sc_context *ctx, const sc_value *a,
                                 const sc_value *b);
sc_value *sc_poly_compose_divconquer(sc_context *ctx, const sc_value *a,
                                     const sc_value *b);
sc_value *sc_poly_taylor_shift(sc_context *ctx, const sc_value *a,
                               const sc_value *b);
sc_value *sc_poly_taylor_shift_horner(sc_context *ctx, const sc_value *a,
                                      const sc_value *b);
sc_value *sc_poly_taylor_shift_divconquer(sc_context *ctx, const sc_value *a,
                                          const sc_value *b);
sc_value *sc_poly_derivative(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_nth_derivative(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_discriminant(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_is_squarefree(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_squarefree_part(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_degree(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_leading_coefficient(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_constant_coefficient(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_coeff(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_reverse(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_truncate(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_shift_left(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_shift_right(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_height(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_max_abs_bits(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_inflate(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_poly_deflation(sc_context *ctx, const sc_value *a);
sc_value *sc_poly_deflate(sc_context *ctx, const sc_value *a, size_t n);

sc_value *sc_zz_poly_generator(sc_context *ctx, sc_parent *parent, const char *name);
sc_value *sc_zz_poly_from_zz(sc_context *ctx, sc_parent *parent, const sc_value *a);
sc_value *sc_zz_poly_add_impl(sc_context *ctx,
                                  const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_sub_impl(sc_context *ctx,
                                  const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_scalar_mul_impl(sc_context *ctx,
                                         const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_neg_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_mul_constant(sc_context *ctx,
                                  const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mul_classical(sc_context *ctx,
                                   const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mullow_classical(sc_context *ctx, const sc_value *a,
                                      const sc_value *b, size_t n);
sc_value *sc_zz_poly_mullow_dc(sc_context *ctx, const sc_value *a,
                               const sc_value *b, size_t n);
sc_value *sc_zz_poly_mulhigh_reverse(sc_context *ctx, const sc_value *a,
                                     const sc_value *b, size_t n);
sc_value *sc_zz_poly_mulmid_classical(sc_context *ctx, const sc_value *a,
                                      const sc_value *b, size_t start, size_t n);
sc_value *sc_zz_poly_mulmid_toom42(sc_context *ctx, const sc_value *a,
                                   const sc_value *b, size_t n);
sc_value *sc_zz_poly_mulmid_toom42_odd(sc_context *ctx, const sc_value *a,
                                       const sc_value *b, size_t n);
sc_value *sc_zz_poly_mulmid_toom63(sc_context *ctx, sc_zz_poly_toom63_ws *ws);
sc_value *sc_zz_poly_mulmid_toom63_tail(sc_context *ctx, const sc_value *a,
                                          const sc_value *b, size_t n);
sc_value *sc_zz_poly_mul_ssa(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mul_ntt(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mulmid_ssa(sc_context *ctx, const sc_value *a,
                                const sc_value *b, size_t n);
sc_value *sc_zz_poly_mulmid_ntt(sc_context *ctx, const sc_value *a,
                                const sc_value *b, size_t n);
size_t sc_zz_poly_ntt_nprimes(const sc_value *a, const sc_value *b);
size_t sc_zz_poly_mulmid_ntt_nprimes(const sc_value *a, const sc_value *b, size_t n);
sc_value *sc_zz_poly_mul_ks(sc_context *ctx, const sc_value *a, const sc_value *b,
                                mp_bitcnt_t bits);
sc_value *sc_zz_poly_mul_karatsuba(sc_context *ctx,
                                   const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mul_toom3(sc_context *ctx, sc_zz_poly_toom3_ws *ws);
sc_value *sc_zz_poly_pow_binary(sc_context *ctx, const sc_value *a, unsigned long e);
sc_value *sc_zz_poly_reverse_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_degree_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_leading_coefficient_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_constant_coefficient_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_coeff_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_truncate_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_shift_left_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_shift_right_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_height_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_max_abs_bits_impl(sc_context *ctx, const sc_value *a);
size_t sc_zz_poly_max_abs_bits_raw(const sc_value *a);
sc_value *sc_zz_poly_inflate_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_deflation_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_deflate_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_quo_classical(sc_context *ctx,
                                   const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_series_quo_classical(sc_context *ctx, const sc_value *a,
                                          const sc_value *b, size_t n);
sc_value *sc_zz_poly_series_quo_dc(sc_context *ctx, const sc_value *a,
                                   const sc_value *b, size_t n);
sc_value *sc_zz_poly_inv_series_classical(sc_context *ctx, const sc_value *a,
                                              size_t n);
sc_value *sc_zz_poly_inv_series(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_inv_series_newton(sc_context *ctx, const sc_value *a,
                                        size_t n);
sc_value *sc_zz_poly_series_quo_preinv(sc_context *ctx, const sc_value *a,
                                       const sc_value *binv, size_t n);
sc_value *sc_zz_poly_series_quo_newton(sc_context *ctx, const sc_value *a,
                                       const sc_value *b, size_t n);
sc_value *sc_zz_poly_quo_dc(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_preinverse_newton(sc_context *ctx, const sc_value *b,
                                           size_t n);
sc_value *sc_zz_poly_quo_preinv(sc_context *ctx, const sc_value *a,
                                const sc_value *b, const sc_value *binv);
sc_value *sc_zz_poly_quo_newton(sc_context *ctx, const sc_value *a,
                                const sc_value *b);
sc_value *sc_zz_poly_divrem_newton(sc_context *ctx, const sc_value *a,
                                   const sc_value *b);
sc_value *sc_zz_poly_quo_bidirectional(sc_context *ctx,
                                           const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_quo_mulders_balanced(sc_context *ctx,
                                          const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_quo_mulders(sc_context *ctx,
                                 const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_divrem_classical(sc_context *ctx,
                                      const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_divrem_dc(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_divrem_full(sc_context *ctx,
                                 const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_divrem_mulders(sc_context *ctx,
                                    const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_pseudodiv_impl(sc_context *ctx,
                                    const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_pseudorem_classical(sc_context *ctx,
                                            const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_pseudorem_fast(sc_context *ctx,
                                    const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_add(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_sub(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mul(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_mullow(sc_context *ctx, const sc_value *a, const sc_value *b,
                            size_t n);
sc_value *sc_zz_poly_mulhigh(sc_context *ctx, const sc_value *a, const sc_value *b,
                             size_t n);
sc_value *sc_zz_poly_mulmid(sc_context *ctx, const sc_value *a, const sc_value *b,
                            size_t start, size_t n);
sc_value *sc_zz_poly_mulmid_balanced(sc_context *ctx, const sc_value *a,
                                     const sc_value *b, size_t n);
sc_value *sc_zz_poly_scalar_mul(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_neg(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_pow(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_divexact(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_quo(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_divrem(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_pseudodiv(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_pseudorem(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_scalar_divexact_impl(sc_context *ctx, const sc_value *a,
                                          const sc_value *b);
sc_value *sc_zz_poly_content_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_primitive_part_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_gcd_pseudo_impl(sc_context *ctx, const sc_value *a,
                                     const sc_value *b);
sc_value *sc_zz_poly_inv_series_scaled_newton(sc_context *ctx, const sc_value *a,
                                                    size_t n);
sc_value *sc_zz_poly_pseudodiv_fast(sc_context *ctx, const sc_value *a,
                                    const sc_value *b);
int sc_zz_poly_mat2_mul(sc_context *ctx, sc_zz_poly_mat2 *r,
                        const sc_zz_poly_mat2 *a, const sc_zz_poly_mat2 *b);
int sc_zz_poly_mat2_apply(sc_context *ctx, sc_value **u, sc_value **v,
                          const sc_zz_poly_mat2 *m, const sc_value *a,
                          const sc_value *b);
int sc_zz_poly_mat2_primitive(sc_context *ctx, sc_zz_poly_mat2 *m);
int sc_zz_poly_hgcd_pseudo(sc_context *ctx, sc_zz_poly_mat2 *m,
                            const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_gcd_hgcd_impl(sc_context *ctx, const sc_value *a,
                                   const sc_value *b);
sc_value *sc_zz_poly_gcd_lr_impl(sc_context *ctx, const sc_value *a,
                                 const sc_value *b);
sc_value *sc_zz_poly_subres_prs_last(sc_context *ctx, const sc_value *a,
                                           const sc_value *b);
sc_value *sc_zz_poly_gcd_subresultant_impl(sc_context *ctx, const sc_value *a,
                                           const sc_value *b);
sc_value *sc_zz_poly_resultant_bareiss_impl(sc_context *ctx,
                                             const sc_value *a,
                                             const sc_value *b);
sc_value *sc_zz_poly_resultant_subresultant_impl(sc_context *ctx,
                                                  const sc_value *a,
                                                  const sc_value *b);
sc_value *sc_zz_poly_xgcd_subresultant_impl(sc_context *ctx, const sc_value *a,
                                             const sc_value *b);
sc_value *sc_zz_poly_evaluate_horner_impl(sc_context *ctx, const sc_value *a,
                                           const sc_value *b);
sc_value *sc_zz_poly_evaluate_divconquer_impl(sc_context *ctx, const sc_value *a,
                                               const sc_value *b);
sc_value *sc_zz_poly_compose_horner_impl(sc_context *ctx, const sc_value *a,
                                         const sc_value *b);
sc_value *sc_zz_poly_compose_divconquer_impl(sc_context *ctx, const sc_value *a,
                                             const sc_value *b);
sc_value *sc_zz_poly_taylor_shift_horner_impl(sc_context *ctx, const sc_value *a,
                                              const sc_value *b);
sc_value *sc_zz_poly_taylor_shift_divconquer_impl(sc_context *ctx,
                                                  const sc_value *a,
                                                  const sc_value *b);
sc_value *sc_zz_poly_derivative_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_nth_derivative_impl(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_discriminant_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_is_squarefree_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_squarefree_part_impl(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_scalar_divexact(sc_context *ctx, const sc_value *a,
                                     const sc_value *b);
sc_value *sc_zz_poly_content(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_primitive_part(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_gcd(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_gcd_pseudo(sc_context *ctx, const sc_value *a,
                                const sc_value *b);
sc_value *sc_zz_poly_gcd_hgcd(sc_context *ctx, const sc_value *a,
                              const sc_value *b);
sc_value *sc_zz_poly_gcd_lr(sc_context *ctx, const sc_value *a,
                            const sc_value *b);
sc_value *sc_zz_poly_resultant(sc_context *ctx, const sc_value *a,
                               const sc_value *b);
sc_value *sc_zz_poly_xgcd(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_evaluate(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_evaluate_horner(sc_context *ctx, const sc_value *a,
                                     const sc_value *b);
sc_value *sc_zz_poly_evaluate_divconquer(sc_context *ctx, const sc_value *a,
                                         const sc_value *b);
sc_value *sc_zz_poly_compose(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_zz_poly_derivative(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_nth_derivative(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_discriminant(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_is_squarefree(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_squarefree_part(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_degree(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_leading_coefficient(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_constant_coefficient(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_coeff(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_reverse(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_truncate_copy(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_shift_left(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_shift_right(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_height(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_max_abs_bits(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_inflate(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_deflation(sc_context *ctx, const sc_value *a);
sc_value *sc_zz_poly_deflate(sc_context *ctx, const sc_value *a, size_t n);
sc_value *sc_zz_poly_compose_horner(sc_context *ctx, const sc_value *a,
                                    const sc_value *b);
sc_value *sc_zz_poly_compose_divconquer(sc_context *ctx, const sc_value *a,
                                        const sc_value *b);
sc_value *sc_zz_poly_taylor_shift(sc_context *ctx, const sc_value *a,
                                  const sc_value *b);
sc_value *sc_zz_poly_taylor_shift_horner(sc_context *ctx, const sc_value *a,
                                         const sc_value *b);
sc_value *sc_zz_poly_taylor_shift_divconquer(sc_context *ctx, const sc_value *a,
                                             const sc_value *b);

sc_value *sc_coerce(sc_context *ctx, sc_parent *parent, const sc_value *a);
sc_value *sc_add(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_sub(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_mul(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_div(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_mod(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_neg(sc_context *ctx, const sc_value *a);
sc_value *sc_pow(sc_context *ctx, const sc_value *a, const sc_value *b);
sc_value *sc_call1(sc_context *ctx, const char *name, const sc_value *a);
sc_value *sc_call_string(sc_context *ctx, const char *name, const char *text);
sc_value *sc_call2(sc_context *ctx, const char *name,
                   const sc_value *a, const sc_value *b);
sc_value *sc_call3(sc_context *ctx, const char *name,
                   const sc_value *a, const sc_value *b, const sc_value *c);

void sc_lexer_set_input(const char *input);
void sc_lexer_clear_input(void);
int sc_parse_line(sc_context *ctx, const char *line);
int sc_repl_line(sc_context *ctx, const char *line);

#endif
