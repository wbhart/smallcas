#include "smallcas.h"

#include <string.h>

int sc_env_set(sc_context *ctx, const char *name, const sc_value *value)
{
    size_t i;
    sc_value *copy;
    char *copy_name;

    if (value == NULL)
        return 0;
    if (strcmp(name, "ZZ") == 0) {
        sc_set_error(ctx, "cannot reassign built-in parent ZZ");
        return 0;
    }
    copy = sc_value_copy_checked(ctx, value);
    if (copy == NULL)
        return 0;
    for (i = 0; i < ctx->nvars; i++) {
        if (strcmp(ctx->vars[i].name, name) == 0) {
            sc_value_free(ctx->vars[i].value);
            ctx->vars[i].value = copy;
            return 1;
        }
    }
    if (ctx->nvars == 128) {
        sc_value_free(copy);
        sc_set_error(ctx, "too many variables");
        return 0;
    }
    copy_name = sc_string_dup(ctx, name);
    if (copy_name == NULL) {
        sc_value_free(copy);
        return 0;
    }
    ctx->vars[ctx->nvars].name = copy_name;
    ctx->vars[ctx->nvars].value = copy;
    ctx->nvars++;
    return 1;
}

sc_value *sc_env_get(sc_context *ctx, const char *name)
{
    size_t i;

    if (strcmp(name, "ZZ") == 0)
        return sc_value_new_parent_checked(ctx, &SC_ZZ);
    for (i = 0; i < ctx->nvars; i++) {
        if (strcmp(ctx->vars[i].name, name) == 0)
            return sc_value_copy_checked(ctx, ctx->vars[i].value);
    }
    sc_set_error(ctx, "unknown identifier '%s'", name);
    return NULL;
}

static int sc_env_pattern_validate(sc_context *ctx, const sc_pattern *pattern,
                                   const sc_value *value)
{
    if (pattern->kind == SC_PATTERN_NAME) {
        if (strcmp(pattern->data.name, "ZZ") != 0)
            return 1;
        sc_set_error(ctx, "cannot reassign built-in parent ZZ");
        return 0;
    }
    if (value == NULL || value->kind != SC_VALUE_PAIR) {
        sc_set_error(ctx, "cannot unpack non-pair value");
        return 0;
    }
    return sc_env_pattern_validate(ctx, pattern->data.pair.first,
                                   value->data.pair.first)
        && sc_env_pattern_validate(ctx, pattern->data.pair.second,
                                   value->data.pair.second);
}

static int sc_env_pattern_bind(sc_context *ctx, const sc_pattern *pattern,
                               const sc_value *value)
{
    if (pattern->kind == SC_PATTERN_NAME)
        return sc_env_set(ctx, pattern->data.name, value);
    return sc_env_pattern_bind(ctx, pattern->data.pair.first,
                               value->data.pair.first)
        && sc_env_pattern_bind(ctx, pattern->data.pair.second,
                               value->data.pair.second);
}

int sc_env_set_pattern(sc_context *ctx, const sc_pattern *pattern,
                       const sc_value *value)
{
    if (pattern == NULL || value == NULL)
        return 0;
    if (!sc_env_pattern_validate(ctx, pattern, value))
        return 0;
    return sc_env_pattern_bind(ctx, pattern, value);
}
