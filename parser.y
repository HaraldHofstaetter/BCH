%{
#include<stdio.h>
#include <stdlib.h> /* qsort */

#include"bch.h" 

int yylex(void); 
void yyerror(char *);

static char *gens;
static int gens_tab[256];
static int num_gens;
static expr_t* result;

%}

%union {     
    char gen;
    rat_t rat;
    expr_t * expr;
};

%token <rat> RAT 
%token <gen> GEN 
%token LOG EXP ID ZERO
%token '+' '-' '*' 
%token '(' ')' '[' ']' ','
%token END 

%left '+' '-'
%left '*' 

%type <expr> expr 

%%

line: expr END { result = $1; YYACCEPT; }

expr: GEN  { $$ = generator(gens_tab[(size_t) $1]); }
    | ZERO { $$ = zero_element(); }
    | ID   { $$ = identity(); }
    | expr '+' expr  { $$ = sum($1, $3); }
    | expr '-' expr  { $$ = difference($1, $3); } 
    | expr '*' expr  { $$ = product($1, $3); } 
 /* | '-' RAT '*' expr   { $$ = term_from_rat(rat_neg($2), $4); } */
    | RAT '*' expr   { $$ = term_from_rat($1, $3); } 
    | '-' expr       { $$ = negation($2); } 
    | '+' expr       { $$ = $2; }
    | LOG '(' expr ')'   { $$ = logarithm($3); } 
    | EXP '(' expr ')'   { $$ = exponential($3); } 
    | '[' expr ',' expr ']'          { $$ = commutator($2, $4); }
    | '(' expr ')'       { $$ = $2; }


%%


void yyerror(char *s) 
{
    fprintf(stderr, "ERROR: while parsing expression: %s\n", s); 
    result = 0; 
} 


typedef struct yy_buffer_state * YY_BUFFER_STATE;
/* extern int yyparse(); */
extern YY_BUFFER_STATE yy_scan_string(char * str);
extern void yy_delete_buffer(YY_BUFFER_STATE buffer);


static int compare_chars(const void *a, const void *b)
{
  const char *da = (const char *) a;
  const char *db = (const char *) b;
  return (*da > *db) - (*da < *db);
}


expr_t* parse(char *inp, char *generators, int *num_generators) {
    for (int i=0; i<256; i++) {
        gens_tab[i] = -1;
    }
    num_gens = 0;
    result = 0;
    gens = generators;

    /* Pass 1: determine all generator symbols in input string*/
    YY_BUFFER_STATE buffer0 = yy_scan_string(inp);
    int tok;
    while ((tok = yylex())) {
        if (tok == END) {
            break;
        }
        if (tok == GEN) {
            if (gens_tab[(size_t) yylval.gen]==-1) {
                 gens_tab[(size_t) yylval.gen] = num_gens;
                 gens[num_gens] = yylval.gen;
                 num_gens++;
            }
        }
    }
    yy_delete_buffer(buffer0);

    /* sort generator symbols */
    qsort(gens, num_gens, sizeof(char), compare_chars);
    for (int i=0; i<num_gens; i++) {
        gens_tab[(size_t) gens[i]] = i;
    }
    
    /* Pass 2: parse input string */
    YY_BUFFER_STATE buffer = yy_scan_string(inp);
    yyparse();     
    yy_delete_buffer(buffer);

    *num_generators = num_gens;
    return result;
}

