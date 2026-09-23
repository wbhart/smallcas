%{
#include "smallcas.h"

#include <stdlib.h>

int yylex(void);
void yyerror(const char *message);

sc_context *sc_parse_context;
%}

%code requires {
#include "smallcas.h"
}

%union {
    char *text;
    sc_value *value;
    sc_pattern *pattern;
}

%token <text> INTEGER IDENT STRING
%type <value> expr
%type <pattern> pattern

%destructor { free($$); } <text>
%destructor { sc_value_free($$); } <value>
%destructor { sc_pattern_free($$); } <pattern>

%left '+' '-'
%left '*' '/' '%'
%precedence UMINUS
%right '^'

%%

line:
      expr
      {
          sc_parse_context->result = $1;
          sc_parse_context->show_result = 1;
      }
    | pattern '=' expr
      {
          sc_env_set_pattern(sc_parse_context, $1, $3);
          sc_pattern_free($1);
          sc_value_free($3);
          sc_parse_context->show_result = 0;
      }
    | pattern ',' pattern '=' expr
      {
          sc_pattern *pair = sc_pattern_new_pair_take_checked(sc_parse_context,
                                                               $1, $3);

          sc_env_set_pattern(sc_parse_context, pair, $5);
          sc_pattern_free(pair);
          sc_value_free($5);
          sc_parse_context->show_result = 0;
      }
    ;

pattern:
      IDENT
      {
          $$ = sc_pattern_new_name_take_checked(sc_parse_context, $1);
      }
    | '(' pattern ',' pattern ')'
      {
          $$ = sc_pattern_new_pair_take_checked(sc_parse_context, $2, $4);
      }
    ;

expr:
      INTEGER
      {
          $$ = sc_zz_from_str(sc_parse_context, $1);
          free($1);
      }
    | IDENT
      {
          $$ = sc_env_get(sc_parse_context, $1);
          free($1);
      }
    | '(' expr ')'
      {
          $$ = $2;
      }
    | '-' expr %prec UMINUS
      {
          $$ = sc_neg(sc_parse_context, $2);
          sc_value_free($2);
      }
    | expr '+' expr
      {
          $$ = sc_add(sc_parse_context, $1, $3);
          sc_value_free($1);
          sc_value_free($3);
      }
    | expr '-' expr
      {
          $$ = sc_sub(sc_parse_context, $1, $3);
          sc_value_free($1);
          sc_value_free($3);
      }
    | expr '*' expr
      {
          $$ = sc_mul(sc_parse_context, $1, $3);
          sc_value_free($1);
          sc_value_free($3);
      }
    | expr '/' expr
      {
          $$ = sc_div(sc_parse_context, $1, $3);
          sc_value_free($1);
          sc_value_free($3);
      }
    | expr '%' expr
      {
          $$ = sc_mod(sc_parse_context, $1, $3);
          sc_value_free($1);
          sc_value_free($3);
      }
    | expr '^' expr
      {
          $$ = sc_pow(sc_parse_context, $1, $3);
          sc_value_free($1);
          sc_value_free($3);
      }
    | IDENT '(' expr ')'
      {
          $$ = sc_call1(sc_parse_context, $1, $3);
          free($1);
          sc_value_free($3);
      }
    | IDENT '(' STRING ')'
      {
          $$ = sc_call_string(sc_parse_context, $1, $3);
          free($1);
          free($3);
      }
    | IDENT '(' expr ',' expr ')'
      {
          $$ = sc_call2(sc_parse_context, $1, $3, $5);
          free($1);
          sc_value_free($3);
          sc_value_free($5);
      }
    | IDENT '(' expr ',' expr ',' expr ')'
      {
          $$ = sc_call3(sc_parse_context, $1, $3, $5, $7);
          free($1);
          sc_value_free($3);
          sc_value_free($5);
          sc_value_free($7);
      }
    ;

%%

void yyerror(const char *message)
{
    sc_set_error(sc_parse_context, "%s", message);
}

int sc_parse_line(sc_context *ctx, const char *line_text)
{
    int status;

    sc_value_free(ctx->result);
    ctx->result = NULL;
    ctx->show_result = 0;
    ctx->error[0] = '\0';
    sc_parse_context = ctx;
    sc_lexer_set_input(line_text);
    status = yyparse();
    sc_lexer_clear_input();
    if (ctx->error[0] != '\0')
        return 0;
    return status == 0;
}
