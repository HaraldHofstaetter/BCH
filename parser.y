%{
#include<stdio.h>

#include"bch.h" 

int yylex(void); 
void yyerror(char *);

typedef struct {
    int num;
    int den;
} rational_t;

static char *gens;
static int gens_tab[256];
static int num_gens;
static expr_t* result;

%}

%union {     
    char gen;
    rational_t rat;
    expr_t * expr;
};

%token <rat> RAT 
%token <gen> GEN 
%token ID LOG EXP
%token '+' '-' '*' 
%token '(' ')' '[' ']' ','
%token END 

%left '+' '-'
%left '*'

%type <expr> expr

%%

line: expr END { result = $1; YYACCEPT; }

expr: GEN { if (gens_tab[(size_t) $1]==-1) {
                gens_tab[(size_t) $1] = num_gens;
                gens[num_gens] = $1;
                num_gens++;
            }
            $$ = generator(gens_tab[(size_t) $1]); 
          }
    | ID                    { $$ = identity(); } 
    | expr '+' expr         { $$ = sum($1, $3); }
    | expr '-' expr         { $$ = difference($1, $3); } 
    | expr '*' expr         { $$ = product($1, $3); }
    | '-' expr              { $$ = negation($2); }
    | '+' expr              { $$ = $2; }
    | RAT '*' expr          { $$ = term($1.num, $1.den, $3); }
    | EXP '(' expr ')'      { $$ = exponential($3); }
    | LOG '(' expr ')'      { $$ = logarithm($3); }
    | '[' expr ',' expr ']' { $$ = commutator($2, $4); }
    | '(' expr ')'          { $$ = $2; }

%%


void yyerror(char *s) 
{
    fprintf(stderr, "ERROR: while parsing expression: %s\n", s); 
    if (result!=0) {
        free_expr(result);
    }
    result = 0; 
} 

typedef struct yy_buffer_state * YY_BUFFER_STATE;
/* extern int yyparse(); */
extern YY_BUFFER_STATE yy_scan_string(char * str);
extern void yy_delete_buffer(YY_BUFFER_STATE buffer);

expr_t* parse(char *inp, char *generators, int *num_generators) {
    for (int i=0; i<256; i++) {
        gens_tab[i] = -1;
    }
    num_gens = 0;
    result = 0;
    gens = generators;

    YY_BUFFER_STATE buffer = yy_scan_string(inp);
    yyparse();     
    yy_delete_buffer(buffer);

    if (result==0) { 
        return 0;
    }

    gens[num_gens] = '\0';
    *num_generators = num_gens;
    return result;
}

/*
int main(void)  // main() for testing only
{  
    while(1) {
    
        char *inp = 0;
        size_t len = 0;
        getline(&inp, &len, stdin);
        
        char generators[56];
        int num_generators;
        expr_t* expr = parse(inp, generators, &num_generators);
        free(inp);
    
        if (expr!=0) { 
            printf("generators = %s # = %d\n", generators, num_generators); 
            printf("expression = "); print_expr(expr, generators); printf("\n");
            free_expr(expr);
         } 
    }

    return 0; 
} 
*/


