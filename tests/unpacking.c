#include "smallcas.h"

#include <stdio.h>
#include <stdlib.h>

static int zz_equal_ui(const sc_value *value, unsigned long n)
{
    return value != NULL && value->kind == SC_VALUE_ZZ
        && mpz_cmp_ui(value->data.z, n) == 0;
}

static sc_pattern *name_pattern(sc_context *ctx, const char *name)
{
    return sc_pattern_new_name_take_checked(ctx, sc_string_dup(ctx, name));
}

int main(void)
{
    sc_context ctx;
    sc_value *one;
    sc_value *two;
    sc_value *three;
    sc_value *inner;
    sc_value *value;
    sc_value *h;
    sc_value *u;
    sc_value *v;
    sc_pattern *pattern;
    sc_pattern *uv;

    sc_context_init(&ctx);
    one = sc_zz_from_str(&ctx, "1");
    two = sc_zz_from_str(&ctx, "2");
    three = sc_zz_from_str(&ctx, "3");
    inner = sc_value_new_pair_take_checked(&ctx, two, three);
    value = sc_value_new_pair_take_checked(&ctx, one, inner);
    uv = sc_pattern_new_pair_take_checked(&ctx, name_pattern(&ctx, "u"),
                                          name_pattern(&ctx, "v"));
    pattern = sc_pattern_new_pair_take_checked(&ctx, name_pattern(&ctx, "h"), uv);
    if (!sc_env_set_pattern(&ctx, pattern, value))
        return 1;
    h = sc_env_get(&ctx, "h");
    u = sc_env_get(&ctx, "u");
    v = sc_env_get(&ctx, "v");
    if (!zz_equal_ui(h, 1) || !zz_equal_ui(u, 2) || !zz_equal_ui(v, 3))
        return 1;
    sc_value_free_many(4, value, h, u, v);
    sc_pattern_free(pattern);
    sc_context_clear(&ctx);
    puts("tuple unpacking tests passed");
    return 0;
}
