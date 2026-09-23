#define _POSIX_C_SOURCE 200809L

#include "smallcas.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void sc_repl_reset(sc_context *ctx)
{
    sc_value_free(ctx->result);
    ctx->result = NULL;
    ctx->show_result = 0;
    ctx->error[0] = '\0';
}

static const char *sc_skip_space(const char *s)
{
    while (isspace((unsigned char)*s))
        s++;
    return s;
}

static int sc_repl_word(const char *s, const char *word)
{
    size_t n = strlen(word);

    return strncmp(s, word, n) == 0 &&
           (s[n] == '\0' || isspace((unsigned char)s[n]));
}

static const char *sc_repl_sep(const char *s, int ch)
{
    int depth = 0;
    int quote = 0;

    for (; *s != '\0'; s++) {
        if (quote != 0) {
            if (*s == quote)
                quote = 0;
        } else if (*s == '\'' || *s == '"') {
            quote = *s;
        } else if (*s == '(') {
            depth++;
        } else if (*s == ')' && depth > 0) {
            depth--;
        } else if (*s == ch && depth == 0) {
            return s;
        }
    }
    return NULL;
}

static int sc_repl_bound(sc_context *ctx, const char *s, size_t n, mpz_t z)
{
    char *text = sc_string_ndup(ctx, s, n);
    int ok;

    if (text == NULL)
        return 0;
    ok = sc_parse_line(ctx, text);
    if (ok && ctx->show_result && ctx->result != NULL &&
        ctx->result->kind == SC_VALUE_ZZ)
        mpz_set(z, ctx->result->data.z);
    else if (ok)
        sc_set_error(ctx, "for-loop bound must be an integer expression");
    free(text);
    return ok && ctx->error[0] == '\0';
}

static int sc_repl_for(sc_context *ctx, const char *line)
{
    const char *p = sc_skip_space(line + 3);
    const char *name = p, *colon, *semi;
    char *var;
    mpz_t lo, hi, i;
    sc_value *index;
    int ok = 0;

    while (isalnum((unsigned char)*p) || *p == '_')
        p++;
    if (p == name || (!isalpha((unsigned char)*name) && *name != '_'))
        goto syntax;
    var = sc_string_ndup(ctx, name, (size_t)(p - name));
    if (var == NULL)
        return 0;
    p = sc_skip_space(p);
    if (*p++ != '=') {
        free(var);
        goto syntax;
    }
    colon = sc_repl_sep(p, ':');
    semi = colon ? sc_repl_sep(colon + 1, ';') : NULL;
    if (colon == NULL || semi == NULL || *sc_skip_space(semi + 1) == '\0') {
        free(var);
        goto syntax;
    }
    if (strcmp(var, "ZZ") == 0) {
        sc_set_error(ctx, "cannot reassign built-in parent ZZ");
        free(var);
        return 0;
    }
    mpz_inits(lo, hi, i, NULL);
    if (!sc_repl_bound(ctx, p, (size_t)(colon - p), lo) ||
        !sc_repl_bound(ctx, colon + 1, (size_t)(semi - colon - 1), hi))
        goto done;
    index = sc_value_new_zz_checked(ctx);
    if (index == NULL)
        goto done;
    for (mpz_set(i, lo); mpz_cmp(i, hi) <= 0; mpz_add_ui(i, i, 1)) {
        mpz_set(index->data.z, i);
        if (!sc_env_set(ctx, var, index) || !sc_parse_line(ctx, semi + 1))
            goto loop_done;
    }
    ok = 1;
loop_done:
    sc_value_free(index);
done:
    free(var);
    mpz_clears(lo, hi, i, NULL);
    sc_value_free(ctx->result);
    ctx->result = NULL;
    ctx->show_result = 0;
    return ok;
syntax:
    sc_set_error(ctx, "use: for name = lo:hi; command");
    return 0;
}

static int sc_repl_time(sc_context *ctx, const char *line)
{
    const char *p = sc_skip_space(line + 5);

    if (sc_repl_word(p, "on") && *sc_skip_space(p + 2) == '\0')
        ctx->timing = 1;
    else if (sc_repl_word(p, "off") && *sc_skip_space(p + 3) == '\0')
        ctx->timing = 0;
    else {
        sc_set_error(ctx, "use: @time on or @time off");
        return 0;
    }
    return 1;
}

int sc_repl_line(sc_context *ctx, const char *line)
{
    struct timespec t0, t1;
    const char *p = sc_skip_space(line);
    int timed, ok;

    if (sc_repl_word(p, "@time")) {
        sc_repl_reset(ctx);
        return sc_repl_time(ctx, p);
    }
    timed = ctx->timing;
    if (timed)
        clock_gettime(CLOCK_MONOTONIC, &t0);
    if (sc_repl_word(p, "for")) {
        sc_repl_reset(ctx);
        ok = sc_repl_for(ctx, p);
    } else {
        ok = sc_parse_line(ctx, line);
    }
    if (timed) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        printf("time: %.9f s\n", (double)(t1.tv_sec - t0.tv_sec) +
               1e-9 * (double)(t1.tv_nsec - t0.tv_nsec));
    }
    return ok;
}
