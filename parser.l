%option noyywrap
%option nounput
%option noinput

natural  [0-9]+
rational {natural}|({natural}\/{natural})
float    (([0-9]+[.][0-9]*)|([0-9]*[.][0-9]+))([eE][-+]?[0-9]+)?

 /* (E(+|-)?[0-9]+)? */

%{
#include"bch.h"

#include"parser.tab.h"

static rat_t parse_RATIONAL(char *s);

%}

%%

[a-zA-Z]   { yylval.gen = *yytext; return GEN; }
"Id"       { return ID;}
"Zero"     { return ZERO;}
"log"      { return LOG;}
"exp"      { return EXP;}
[\+\-\*\(\)\[\],]  { return yytext[0]; }
{rational} { yylval.rat = parse_RATIONAL(yytext); return RAT; } 
{float}    { yylval.flt = parse_FLOAT(yytext); return FLT; }  
[ \t]      { /* ignore whitespace */ }
<<EOF>>    { return END; } 

%%

static rat_t parse_RATIONAL(char *s)
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

#include <stdio.h>

/* generic n-th power of a FLOAT number
 * based on npwr() from src/dd_real.cpp in
 * https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz 
 * by David H. Bailey
 */
static FLOAT npwr_f(FLOAT a, int n) {  
  if (n == 0) {
    if (is_zero_f(a)) {
      fprintf(stderr, "PANIC: npwr_f: invalid argument");
      abort(); 
    }
    return i2f(1);
  }

  FLOAT r = a;
  FLOAT s = i2f(1);
  int N = abs(n);

  if (N > 1) {
    /* Use binary exponentiation */
    while (N > 0) {
      if (N % 2 == 1) {
        s = mul_f(s, r);
      }
      N /= 2;
      if (N > 0)
        r = mul_f(r, r);
    }
  } else {
    s = r;
  }

  /* Compute the reciprocal if n is negative. */
  if (n < 0)
    return div_f(i2f(1), s);
  
  return s;
}

/* generic parser for FLOAT numbers 
 * based on dd_real::read from src/dd_real.cpp in
 * https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz 
 * by David H. Bailey
 */
FLOAT generic_parse_FLOAT(char *s) {
  char *p = s;
  char ch;
  int sign = 0;
  int point = -1;
  int nd = 0;
  int e = 0;
  int done = 0;
  FLOAT r = zero_f();
  int nread;
  
  /* Skip any leading spaces */
  while (*p == ' ')
    p++;

  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') {
      int d = ch - '0';
      r = mul_f(r, i2f(10));
      r = add_f(r, i2f(d));
      nd++;
    } else {

      switch (ch) {

      case '.':
        if (point >= 0) {
          fprintf(stderr, "PANIC: generic_parse_FLOAT: more than one '.'");
          abort();
        }
        point = nd;
        break;

      case '-':
      case '+':
        if (sign != 0 || nd > 0) {
          fprintf(stderr, "PANIC: generic_parse_FLOAT: more than one '+' or '-'");
          abort();
        }
        sign = (ch == '-') ? -1 : 1;
        break;

      case 'E':
      case 'e':
        nread = sscanf(p+1, "%d", &e);
        done = 1;
        if (nread != 1) {
          fprintf(stderr, "PANIC: generic_parse_FLOAT: bad exponent");
          abort();
        }
        break;

      default:
        fprintf(stderr, "PANIC: generic_parse_FLOAT: unknwown character");
        abort();
      }
    }
    
    p++;
  }

  if (point >= 0) {
    e -= (nd - point);
  }

  if (e != 0) {
    r = mul_f(r, npwr_f(i2f(10), e));
  }

  if (sign == -1) {
    r = neg_f(r);
  }
  
  return r;
}


    


