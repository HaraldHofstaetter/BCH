%option noyywrap
%option nounput
%option noinput

natural  [0-9]+
rational {natural}|({natural}\/{natural})

%{
#include"bch.h"

typedef struct {
    int num;
    int den;
} rational_t;

#include"parser.tab.h"

static rational_t parse_rational(char *s);

%}

%%

[a-zA-Z]   { yylval.gen = *yytext; return GEN; }
"log"      { return LOG;}
"exp"      { return EXP;}
[\+\-\*\(\)\[\],]  { return yytext[0]; }
{rational} { yylval.rat = parse_rational(yytext); return RAT; } 
[ \t]      { /* ignore whitespace */ }
<<EOF>>    { return END; } 

%%

static rational_t parse_rational(char *s)
{
    rational_t r;
    char *t = s;
    while ((*t!='/')&&(*t!='\0')) { t++; }
    if (*t=='/') {
        *t = '\0';
        r.num = atoi(s);
        r.den = atoi(t+1);
        *t = '/';
    }
    else {
        r.num = atoi(s);
        r.den = 1;
    }
    return r;
}




    


