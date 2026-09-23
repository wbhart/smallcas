#include "smallcas.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

char *sc_string_dup(sc_context *ctx, const char *text)
{
    size_t n = strlen(text) + 1;
    char *copy = malloc(n);

    if (copy == NULL) {
        sc_set_error(ctx, "out of memory");
        return NULL;
    }
    memcpy(copy, text, n);
    return copy;
}

char *sc_string_ndup(sc_context *ctx, const char *text, size_t n)
{
    char *copy = malloc(n + 1);

    if (copy == NULL) {
        sc_set_error(ctx, "out of memory");
        return NULL;
    }
    memcpy(copy, text, n);
    copy[n] = '\0';
    return copy;
}

sc_parent *sc_parent_new_checked(sc_context *ctx)
{
    sc_parent *parent = calloc(1, sizeof(*parent));

    if (parent == NULL)
        sc_set_error(ctx, "out of memory");
    return parent;
}

sc_value *sc_value_new_zz(void)
{
    sc_value *value = malloc(sizeof(*value));

    if (value == NULL)
        return NULL;
    value->kind = SC_VALUE_ZZ;
    value->parent = &SC_ZZ;
    mpz_init(value->data.z);
    return value;
}

sc_value *sc_value_new_zz_checked(sc_context *ctx)
{
    sc_value *value = sc_value_new_zz();

    if (value == NULL)
        sc_set_error(ctx, "out of memory");
    return value;
}

sc_value *sc_value_new_zz_poly(sc_parent *parent, size_t length)
{
    sc_value *value;
    size_t i;

    if (length > (size_t)-1 / sizeof(mpz_t))
        return NULL;
    value = malloc(sizeof(*value));
    if (value == NULL)
        return NULL;
    value->kind = SC_VALUE_POLY;
    value->parent = parent;
    value->data.zz_poly.length = length;
    value->data.zz_poly.coeff = length ? malloc(length * sizeof(mpz_t)) : NULL;
    if (length != 0 && value->data.zz_poly.coeff == NULL) {
        free(value);
        return NULL;
    }
    for (i = 0; i < length; i++)
        mpz_init(value->data.zz_poly.coeff[i]);
    return value;
}

sc_value *sc_value_new_zz_poly_checked(sc_context *ctx, sc_parent *parent,
                                       size_t length)
{
    sc_value *value = sc_value_new_zz_poly(parent, length);

    if (value == NULL)
        sc_set_error(ctx, "out of memory");
    return value;
}

sc_value *sc_value_new_parent(sc_parent *parent)
{
    sc_value *value = malloc(sizeof(*value));

    if (value == NULL)
        return NULL;
    value->kind = SC_VALUE_PARENT;
    value->parent = NULL;
    value->data.parent_value = parent;
    return value;
}

sc_value *sc_value_new_parent_checked(sc_context *ctx, sc_parent *parent)
{
    sc_value *value = sc_value_new_parent(parent);

    if (value == NULL)
        sc_set_error(ctx, "out of memory");
    return value;
}

sc_value *sc_value_new_pair_take_checked(sc_context *ctx, sc_value *first,
                                         sc_value *second)
{
    sc_value *value = malloc(sizeof(*value));

    if (value == NULL) {
        sc_value_free_many(2, first, second);
        sc_set_error(ctx, "out of memory");
        return NULL;
    }
    value->kind = SC_VALUE_PAIR;
    value->parent = NULL;
    value->data.pair.first = first;
    value->data.pair.second = second;
    return value;
}

sc_value *sc_value_copy(const sc_value *value)
{
    sc_value *copy;
    size_t i;

    if (value == NULL)
        return NULL;
    if (value->kind == SC_VALUE_PARENT)
        return sc_value_new_parent(value->data.parent_value);
    if (value->kind == SC_VALUE_PAIR) {
        sc_value *first = sc_value_copy(value->data.pair.first);
        sc_value *second = sc_value_copy(value->data.pair.second);

        if (first == NULL || second == NULL)
            return sc_value_free_many_null(2, first, second);
        copy = malloc(sizeof(*copy));
        if (copy == NULL)
            return sc_value_free_many_null(2, first, second);
        copy->kind = SC_VALUE_PAIR;
        copy->parent = NULL;
        copy->data.pair.first = first;
        copy->data.pair.second = second;
        return copy;
    }
    if (value->kind == SC_VALUE_ZZ) {
        copy = sc_value_new_zz();
        if (copy != NULL)
            mpz_set(copy->data.z, value->data.z);
        return copy;
    }
    if (value->parent->base != &SC_ZZ)
        return NULL;
    copy = sc_value_new_zz_poly(value->parent, value->data.zz_poly.length);
    if (copy == NULL)
        return NULL;
    for (i = 0; i < value->data.zz_poly.length; i++)
        mpz_set(copy->data.zz_poly.coeff[i], value->data.zz_poly.coeff[i]);
    return copy;
}

sc_value *sc_value_copy_checked(sc_context *ctx, const sc_value *value)
{
    sc_value *copy = sc_value_copy(value);

    if (copy == NULL)
        sc_set_error(ctx, "out of memory");
    return copy;
}

void sc_value_free(sc_value *value)
{
    size_t i;

    if (value == NULL)
        return;
    if (value->kind == SC_VALUE_ZZ)
        mpz_clear(value->data.z);
    if (value->kind == SC_VALUE_PAIR)
        sc_value_free_many(2, value->data.pair.first, value->data.pair.second);
    if (value->kind == SC_VALUE_POLY && value->parent->base == &SC_ZZ) {
        for (i = 0; i < value->data.zz_poly.length; i++)
            mpz_clear(value->data.zz_poly.coeff[i]);
        free(value->data.zz_poly.coeff);
    }
    free(value);
}

void sc_value_free_null(sc_value **value)
{
    sc_value_free(*value);
    *value = NULL;
}

void sc_value_free_many(size_t count, ...)
{
    va_list ap;
    size_t i;

    va_start(ap, count);
    for (i = 0; i < count; i++)
        sc_value_free(va_arg(ap, sc_value *));
    va_end(ap);
}

sc_value **sc_value_array_new_checked(sc_context *ctx, size_t count)
{
    sc_value **values = calloc(count ? count : 1, sizeof(*values));

    if (values == NULL)
        sc_set_error(ctx, "out of memory");
    return values;
}

void sc_value_array_free(size_t count, sc_value **values)
{
    size_t i;

    if (values == NULL)
        return;
    for (i = 0; i < count; i++)
        sc_value_free(values[i]);
    free(values);
}

sc_value *sc_value_free_many_null(size_t count, ...)
{
    va_list ap;
    size_t i;

    va_start(ap, count);
    for (i = 0; i < count; i++)
        sc_value_free(va_arg(ap, sc_value *));
    va_end(ap);
    return NULL;
}

sc_value *sc_value_pair_take(sc_value *pair, int second)
{
    sc_value *value;

    if (pair == NULL || pair->kind != SC_VALUE_PAIR)
        return NULL;
    value = second ? pair->data.pair.second : pair->data.pair.first;
    if (second)
        pair->data.pair.second = NULL;
    else
        pair->data.pair.first = NULL;
    sc_value_free(pair);
    return value;
}

int sc_value_pair_split(sc_value *pair, sc_value **first, sc_value **second)
{
    if (pair == NULL || pair->kind != SC_VALUE_PAIR)
        return 0;
    *first = pair->data.pair.first;
    *second = pair->data.pair.second;
    pair->data.pair.first = NULL;
    pair->data.pair.second = NULL;
    sc_value_free(pair);
    return 1;
}

sc_pattern *sc_pattern_new_name_take_checked(sc_context *ctx, char *name)
{
    sc_pattern *pattern = malloc(sizeof(*pattern));

    if (pattern == NULL) {
        free(name);
        sc_set_error(ctx, "out of memory");
        return NULL;
    }
    pattern->kind = SC_PATTERN_NAME;
    pattern->data.name = name;
    return pattern;
}

sc_pattern *sc_pattern_new_pair_take_checked(sc_context *ctx, sc_pattern *first,
                                             sc_pattern *second)
{
    sc_pattern *pattern = malloc(sizeof(*pattern));

    if (pattern == NULL) {
        sc_pattern_free(first);
        sc_pattern_free(second);
        sc_set_error(ctx, "out of memory");
        return NULL;
    }
    pattern->kind = SC_PATTERN_PAIR;
    pattern->data.pair.first = first;
    pattern->data.pair.second = second;
    return pattern;
}

void sc_pattern_free(sc_pattern *pattern)
{
    if (pattern == NULL)
        return;
    if (pattern->kind == SC_PATTERN_NAME)
        free(pattern->data.name);
    else {
        sc_pattern_free(pattern->data.pair.first);
        sc_pattern_free(pattern->data.pair.second);
    }
    free(pattern);
}

void sc_zz_poly_truncate(sc_value *value, size_t length)
{
    size_t i;

    for (i = length; i < value->data.zz_poly.length; i++)
        mpz_clear(value->data.zz_poly.coeff[i]);
    value->data.zz_poly.length = length;
}

void sc_zz_poly_normalize(sc_value *value)
{
    size_t length = value->data.zz_poly.length;

    while (length != 0 && mpz_sgn(value->data.zz_poly.coeff[length - 1]) == 0)
        length--;
    sc_zz_poly_truncate(value, length);
}

sc_value sc_zz_poly_view(const sc_value *value, size_t start, size_t length)
{
    sc_value view = *value;
    size_t n = value->data.zz_poly.length;

    if (start >= n) {
        start = n;
        length = 0;
    } else if (length > n - start) {
        length = n - start;
    }
    while (length != 0 &&
           mpz_sgn(value->data.zz_poly.coeff[start + length - 1]) == 0)
        length--;
    view.data.zz_poly.length = length;
    view.data.zz_poly.coeff = value->data.zz_poly.coeff;
    if (view.data.zz_poly.coeff != NULL)
        view.data.zz_poly.coeff += start;
    return view;
}

static void sc_zz_poly_toom3_ws_clear(sc_zz_poly_toom3_ws *ws, int keep_result)
{
    size_t i, j;

    for (i = 0; i < 2; i++)
        for (j = 1; j <= 3; j++)
            sc_value_free(ws->point[i][j]);
    for (i = 0; i < 5; i++)
        sc_value_free(ws->prod[i]);
    if (!keep_result)
        sc_value_free(ws->result);
    mpz_clear(ws->zero);
}

int sc_zz_poly_toom3_ws_init(sc_context *ctx, sc_zz_poly_toom3_ws *ws,
                             const sc_value *a, const sc_value *b, size_t m)
{
    const sc_value *src[2] = { a, b };
    size_t i, j;

    memset(ws, 0, sizeof(*ws));
    ws->m = m;
    ws->q = 2 * m - 1;
    ws->length = a->data.zz_poly.length + b->data.zz_poly.length - 1;
    mpz_init(ws->zero);
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++)
            ws->block[i][j] = sc_zz_poly_view(src[i], j * m, m);
        ws->point[i][0] = &ws->block[i][0];
        ws->point[i][4] = &ws->block[i][2];
        for (j = 1; j <= 3; j++) {
            ws->point[i][j] = sc_value_new_zz_poly_checked(ctx, a->parent, m);
            if (ws->point[i][j] == NULL) {
                sc_zz_poly_toom3_ws_abort(ws);
                return 0;
            }
        }
    }
    ws->result = sc_value_new_zz_poly_checked(ctx, a->parent, 6 * m - 1);
    if (ws->result == NULL) {
        sc_zz_poly_toom3_ws_abort(ws);
        return 0;
    }
    return 1;
}

sc_value *sc_zz_poly_toom3_ws_abort(sc_zz_poly_toom3_ws *ws)
{
    sc_zz_poly_toom3_ws_clear(ws, 0);
    return NULL;
}

sc_value *sc_zz_poly_toom3_ws_finish(sc_zz_poly_toom3_ws *ws)
{
    sc_value *result = ws->result;

    sc_zz_poly_truncate(result, ws->length);
    sc_zz_poly_toom3_ws_clear(ws, 1);
    return result;
}

static void sc_zz_poly_toom63_ws_clear(sc_zz_poly_toom63_ws *ws, int keep_result)
{
    size_t i;

    for (i = 0; i < 5; i++) {
        sc_value_free(ws->fpoint[i]);
        sc_value_free(ws->prod[i]);
    }
    for (i = 1; i <= 3; i++)
        sc_value_free(ws->gpoint[i]);
    if (!keep_result)
        sc_value_free(ws->result);
    mpz_clear(ws->zero);
}

int sc_zz_poly_toom63_ws_init(sc_context *ctx, sc_zz_poly_toom63_ws *ws,
                              const sc_value *a, const sc_value *b, size_t k)
{
    size_t i;

    memset(ws, 0, sizeof(*ws));
    ws->k = k;
    mpz_init(ws->zero);
    for (i = 0; i < 5; i++) {
        ws->fblock[i] = sc_zz_poly_view(a, i * k, 2 * k - 1);
        ws->fpoint[i] = sc_value_new_zz_poly_checked(ctx, a->parent, 2 * k - 1);
        if (ws->fpoint[i] == NULL) {
            sc_zz_poly_toom63_ws_abort(ws);
            return 0;
        }
    }
    for (i = 0; i < 3; i++)
        ws->gblock[i] = sc_zz_poly_view(b, i * k, k);
    ws->gpoint[0] = &ws->gblock[2];
    ws->gpoint[4] = &ws->gblock[0];
    for (i = 1; i <= 3; i++) {
        ws->gpoint[i] = sc_value_new_zz_poly_checked(ctx, a->parent, k);
        if (ws->gpoint[i] == NULL) {
            sc_zz_poly_toom63_ws_abort(ws);
            return 0;
        }
    }
    ws->result = sc_value_new_zz_poly_checked(ctx, a->parent, 3 * k);
    if (ws->result == NULL) {
        sc_zz_poly_toom63_ws_abort(ws);
        return 0;
    }
    return 1;
}

sc_value *sc_zz_poly_toom63_ws_abort(sc_zz_poly_toom63_ws *ws)
{
    sc_zz_poly_toom63_ws_clear(ws, 0);
    return NULL;
}

sc_value *sc_zz_poly_toom63_ws_finish(sc_zz_poly_toom63_ws *ws)
{
    sc_value *result = ws->result;

    sc_zz_poly_normalize(result);
    sc_zz_poly_toom63_ws_clear(ws, 1);
    return result;
}

int sc_zz_poly_subres_ws_init(sc_context *ctx, sc_zz_poly_subres_ws *ws,
                              const sc_value *a, const sc_value *b)
{
    memset(ws, 0, sizeof(*ws));
    mpz_inits(ws->hp, ws->hc, ws->hn, ws->den, NULL);
    ws->u = sc_value_copy_checked(ctx, a);
    ws->v = sc_value_copy_checked(ctx, b);
    if (ws->u == NULL || ws->v == NULL) {
        sc_value_free_many(2, ws->u, ws->v);
        mpz_clears(ws->hp, ws->hc, ws->hn, ws->den, NULL);
        return 0;
    }
    return 1;
}

sc_value *sc_zz_poly_subres_ws_finish(sc_context *ctx, sc_zz_poly_subres_ws *ws,
                                      sc_value *keep, mpz_srcptr h)
{
    sc_value *hz = sc_value_new_zz_checked(ctx);

    if (ws->u != keep)
        sc_value_free(ws->u);
    if (ws->v != keep)
        sc_value_free(ws->v);
    if (ws->w != keep)
        sc_value_free(ws->w);
    if (ws->z != keep)
        sc_value_free(ws->z);
    if (hz != NULL)
        mpz_set(hz->data.z, h);
    mpz_clears(ws->hp, ws->hc, ws->hn, ws->den, NULL);
    if (hz == NULL) {
        sc_value_free(keep);
        return NULL;
    }
    return sc_value_new_pair_take_checked(ctx, keep, hz);
}

void sc_zz_poly_subres_ws_shift(sc_zz_poly_subres_ws *ws)
{
    sc_value_free(ws->v);
    ws->v = ws->w;
    ws->w = ws->z;
    ws->z = NULL;
}

static sc_value *sc_zz_poly_constant_checked(sc_context *ctx, sc_parent *parent,
                                             unsigned long value)
{
    sc_value *r = sc_value_new_zz_poly_checked(ctx, parent, value != 0 ? 1 : 0);

    if (r != NULL && value != 0)
        mpz_set_ui(r->data.zz_poly.coeff[0], value);
    return r;
}

int sc_zz_poly_xgcd_ws_init(sc_context *ctx, sc_zz_poly_xgcd_ws *ws,
                            const sc_value *a, const sc_value *b)
{
    memset(ws, 0, sizeof(*ws));
    mpz_inits(ws->hp, ws->hc, ws->hn, ws->den, ws->alpha, NULL);
    ws->u = sc_value_copy_checked(ctx, a);
    ws->v = sc_value_copy_checked(ctx, b);
    ws->cu[0] = sc_zz_poly_constant_checked(ctx, a->parent, 1);
    ws->cu[1] = sc_zz_poly_constant_checked(ctx, a->parent, 0);
    ws->cv[0] = sc_zz_poly_constant_checked(ctx, a->parent, 0);
    ws->cv[1] = sc_zz_poly_constant_checked(ctx, a->parent, 1);
    if (ws->u != NULL && ws->v != NULL && ws->cu[0] != NULL && ws->cu[1] != NULL &&
        ws->cv[0] != NULL && ws->cv[1] != NULL)
        return 1;
    sc_value_free_many(6, ws->u, ws->v, ws->cu[0], ws->cu[1], ws->cv[0], ws->cv[1]);
    mpz_clears(ws->hp, ws->hc, ws->hn, ws->den, ws->alpha, NULL);
    return 0;
}

sc_value *sc_zz_poly_xgcd_ws_finish(sc_context *ctx, sc_zz_poly_xgcd_ws *ws)
{
    sc_value *cof = sc_value_new_pair_take_checked(ctx, ws->cv[0], ws->cv[1]);
    sc_value *r = NULL;

    ws->cv[0] = NULL;
    ws->cv[1] = NULL;
    if (cof != NULL) {
        r = sc_value_new_pair_take_checked(ctx, ws->v, cof);
        ws->v = NULL;
    }
    sc_value_free_many(6, ws->u, ws->v, ws->w, ws->cu[0], ws->cu[1],
                       ws->cw[0]);
    sc_value_free(ws->cw[1]);
    mpz_clears(ws->hp, ws->hc, ws->hn, ws->den, ws->alpha, NULL);
    return r;
}

sc_value *sc_zz_poly_xgcd_ws_abort(sc_zz_poly_xgcd_ws *ws)
{
    sc_value_free_many(9, ws->u, ws->v, ws->w, ws->cu[0], ws->cu[1],
                       ws->cv[0], ws->cv[1], ws->cw[0], ws->cw[1]);
    mpz_clears(ws->hp, ws->hc, ws->hn, ws->den, ws->alpha, NULL);
    return NULL;
}

void sc_zz_poly_xgcd_ws_shift(sc_zz_poly_xgcd_ws *ws)
{
    size_t j;

    sc_value_free(ws->u);
    ws->u = ws->v;
    ws->v = ws->w;
    ws->w = NULL;
    for (j = 0; j < 2; j++) {
        sc_value_free(ws->cu[j]);
        ws->cu[j] = ws->cv[j];
        ws->cv[j] = ws->cw[j];
        ws->cw[j] = NULL;
    }
}

int sc_zz_poly_mat2_identity(sc_context *ctx, sc_zz_poly_mat2 *m, sc_parent *parent)
{
    memset(m, 0, sizeof(*m));
    m->a00 = sc_zz_poly_constant_checked(ctx, parent, 1);
    m->a01 = sc_zz_poly_constant_checked(ctx, parent, 0);
    m->a10 = sc_zz_poly_constant_checked(ctx, parent, 0);
    m->a11 = sc_zz_poly_constant_checked(ctx, parent, 1);
    if (m->a00 != NULL && m->a01 != NULL && m->a10 != NULL && m->a11 != NULL)
        return 1;
    sc_zz_poly_mat2_clear(m);
    return 0;
}

void sc_zz_poly_mat2_clear(sc_zz_poly_mat2 *m)
{
    if (m == NULL)
        return;
    sc_value_free_many(4, m->a00, m->a01, m->a10, m->a11);
    memset(m, 0, sizeof(*m));
}

void sc_zz_poly_mat2_clear_many(size_t count, ...)
{
    va_list ap;
    size_t i;

    va_start(ap, count);
    for (i = 0; i < count; i++)
        sc_zz_poly_mat2_clear(va_arg(ap, sc_zz_poly_mat2 *));
    va_end(ap);
}

void sc_zz_poly_mat2_move(sc_zz_poly_mat2 *dst, sc_zz_poly_mat2 *src)
{
    sc_zz_poly_mat2_clear(dst);
    *dst = *src;
    memset(src, 0, sizeof(*src));
}

void sc_value_ptr_swap(sc_value **a, sc_value **b)
{
    sc_value *t = *a;

    *a = *b;
    *b = t;
}
