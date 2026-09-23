%{
#include "smallcas.h"
#include "parser.h"
%}

%option noyywrap noinput nounput nodefault never-interactive

%%

[ \t\r\n]+              ;
[0-9]+                  {
                            yylval.text = sc_string_dup(sc_parse_context, yytext);
                            if (yylval.text == NULL)
                                return 0;
                            return INTEGER;
                        }
[A-Za-z_][A-Za-z0-9_]*  {
                            yylval.text = sc_string_dup(sc_parse_context, yytext);
                            if (yylval.text == NULL)
                                return 0;
                            return IDENT;
                        }
\'[^\'\n]*\'             {
                            yylval.text = sc_string_ndup(sc_parse_context,
                                                         yytext + 1, yyleng - 2);
                            if (yylval.text == NULL)
                                return 0;
                            return STRING;
                        }
\"[^\"\n]*\"             {
                            yylval.text = sc_string_ndup(sc_parse_context,
                                                         yytext + 1, yyleng - 2);
                            if (yylval.text == NULL)
                                return 0;
                            return STRING;
                        }
[()+*/%^,=\-]           { return (unsigned char) yytext[0]; }
.                       { return (unsigned char) yytext[0]; }

%%

void sc_lexer_set_input(const char *input)
{
    yy_scan_string(input);
}

void sc_lexer_clear_input(void)
{
    yylex_destroy();
}
