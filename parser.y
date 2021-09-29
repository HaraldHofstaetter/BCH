%{
#include<stdio.h>

#include"bch.h" 

int yylex(void); 
void yyerror(char *);

typedef struct {
    int num;
    int den;
} rational_t;

static int gcd(int a, int b);

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
%token LOG EXP 
%token '+' '-' '*' 
%token '(' ')' '[' ']' ','
%token END 

%left '+' '-'
%left '*' 

%type <expr> lie_expr prod_of_exp_of_lie_exprs

%%

line: lie_expr END { result = $1; YYACCEPT; }

lie_expr: GEN { if (gens_tab[(size_t) $1]==-1) {
                    gens_tab[(size_t) $1] = num_gens;
                    gens[num_gens] = $1;
                    num_gens++;
                 }
                 $$ = generator(gens_tab[(size_t) $1]); 
              }
    | lie_expr '+' lie_expr  { $$ = sum($1, $3); }
    | lie_expr '-' lie_expr  { $$ = difference($1, $3); } 
 /* | '-' RAT '*' lie_expr   { if ($2.den==0) { 
                                   yyerror("zero denominator");
                                   YYABORT;
                               }
                               int d = gcd($2.num, $2.den);
                               $$ = term(-$2.num/d, $2.den/d, $4); 
                             } */
    | RAT '*' lie_expr       { if ($1.den==0) { 
                                   yyerror("zero denominator");
                                   YYABORT;
                               }
                               int d = gcd($1.num, $1.den);
                               $$ = term($1.num/d, $1.den/d, $3); 
                             } 
    | '-' lie_expr           { $$ = negation($2); } 
    | '+' lie_expr           { $$ = $2; }
    | LOG '(' prod_of_exp_of_lie_exprs ')'   { $$ = logarithm($3); } 
    | '[' lie_expr ',' lie_expr ']'          { $$ = commutator($2, $4); }
    | '(' lie_expr ')'       { $$ = $2; }


prod_of_exp_of_lie_exprs: prod_of_exp_of_lie_exprs  '*' prod_of_exp_of_lie_exprs {
                                                 $$ = product($1, $3);
                                               } 
                        | EXP '(' lie_expr ')' { $$ = exponential($3); }


%%


void yyerror(char *s) 
{
    fprintf(stderr, "ERROR: while parsing Lie expression: %s\n", s); 
    if (result!=0) {
        free_expr(result);
    }
    result = 0; 
} 


static int gcd(int a, int b) {
    /* computes greatest common divisor of a and b
     * METHOD: Euclid's classical algorithm
     */
    while (b!=0) {
       INTEGER t = b; 
       b = a%b; 
       a = t; 
    }
    return a>=0 ? a : -a;
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


