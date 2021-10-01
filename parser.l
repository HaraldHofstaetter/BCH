%option noyywrap
%option nounput
%option noinput

natural  [0-9]+
rational {natural}|({natural}\/{natural})

%{
#include"bch.h"

#include"parser.tab.h"

static rat_t parse_rational(char *s);

%}

%%

[a-zA-Z]   { yylval.gen = *yytext; return GEN; }
"Id"       { return ID;}
"Zero"     { return ZERO;}
"log"      { return LOG;}
"exp"      { return EXP;}
[\+\-\*\(\)\[\],]  { return yytext[0]; }
{rational} { yylval.rat = parse_rational(yytext); return RAT; } 
[ \t]      { /* ignore whitespace */ }
<<EOF>>    { return END; } 

%%

static rat_t parse_rational(char *s)
{
    rat_t r;
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




    


