#include <stdio.h>

#include "smallcas.h"

#include <stdarg.h>
#include <stdlib.h>

sc_parent SC_ZZ = { "ZZ", SC_PARENT_ZZ, NULL, NULL };

void sc_context_init(sc_context *ctx)
{
    ctx->nvars = 0;
    ctx->nparents = 0;
    ctx->result = NULL;
    ctx->show_result = 0;
    ctx->timing = 0;
    ctx->error[0] = '\0';
}

void sc_context_clear(sc_context *ctx)
{
    size_t i;

    sc_value_free(ctx->result);
    for (i = 0; i < ctx->nvars; i++) {
        free(ctx->vars[i].name);
        sc_value_free(ctx->vars[i].value);
    }
    for (i = 0; i < ctx->nparents; i++) {
        free(ctx->parents[i]->symbol);
        free(ctx->parents[i]);
    }
    ctx->nvars = 0;
    ctx->nparents = 0;
}

void sc_set_error(sc_context *ctx, const char *fmt, ...)
{
    va_list ap;

    if (ctx->error[0] != '\0')
        return;
    va_start(ap, fmt);
    vsnprintf(ctx->error, sizeof(ctx->error), fmt, ap);
    va_end(ap);
}

static void sc_zz_poly_fprint(FILE *out, const sc_value *value)
{
    const sc_zz_poly *p = &value->data.zz_poly;
    const char *x = value->parent->symbol ? value->parent->symbol : "x";
    mpz_t mag;
    size_t i;
    int first = 1;

    if (p->length == 0) {
        fputc('0', out);
        return;
    }
    mpz_init(mag);
    for (i = p->length; i-- > 0;) {
        int sign = mpz_sgn(p->coeff[i]);

        if (sign == 0)
            continue;
        if (!first)
            fputs(sign < 0 ? " - " : " + ", out);
        else if (sign < 0)
            fputc('-', out);
        mpz_abs(mag, p->coeff[i]);
        if (i == 0 || mpz_cmp_ui(mag, 1) != 0)
            gmp_fprintf(out, "%Zd", mag);
        if (i > 0) {
            if (mpz_cmp_ui(mag, 1) != 0)
                fputc('*', out);
            fputs(x, out);
            if (i > 1)
                fprintf(out, "^%zu", i);
        }
        first = 0;
    }
    mpz_clear(mag);
}

static void sc_value_fprint(FILE *out, const sc_value *value)
{
    if (value == NULL) {
        fputs("<null>", out);
        return;
    }
    if (value->kind == SC_VALUE_PARENT) {
        fputs(value->data.parent_value->name, out);
        return;
    }
    if (value->kind == SC_VALUE_ZZ) {
        gmp_fprintf(out, "%Zd", value->data.z);
        return;
    }
    if (value->kind == SC_VALUE_POLY && value->parent->base == &SC_ZZ) {
        sc_zz_poly_fprint(out, value);
        return;
    }
    if (value->kind == SC_VALUE_PAIR) {
        fputc('(', out);
        sc_value_fprint(out, value->data.pair.first);
        fputs(", ", out);
        sc_value_fprint(out, value->data.pair.second);
        fputc(')', out);
        return;
    }
    fputs("<unprintable value>", out);
}

void sc_value_print(const sc_value *value)
{
    sc_value_fprint(stdout, value);
    fputc('\n', stdout);
}
