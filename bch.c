#include"bch.h"
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
static void call_onMainFinished(void) 
{
    EM_ASM(onMainFinished(););
}    
#endif

static
long long int get_arg(int argc, char*argv[], char* varname, 
                      long long int default_value, 
                      long long int min, 
                      long long int max)
{
    for (int k=1; k<argc; k++) {
        char* sep = strchr(argv[k], '=');
        if (sep) {
            *sep ='\0';
            if (strcmp(argv[k], varname)==0) {
                char *endptr;
                long long int value = strtoll(sep+1, &endptr, 10); 
                if ((*endptr != '\0')||(endptr == sep+1)) {
                    fprintf(stderr, "ERROR: expected %s=integer\n", argv[k]);
                    exit(EXIT_FAILURE);
                }
                if (value<min) {
                    fprintf(stderr, "ERROR: expected %s>=%lli, got %lli\n", argv[k], min, value);
                    exit(EXIT_FAILURE);
                }
                if (value>max) {
                    fprintf(stderr, "ERROR: expected %s<=%lli, got %lli\n", argv[k], max, value);
                    exit(EXIT_FAILURE);
                }
                *sep = '=';
                return value;
            }
            *sep = '=';
        }   
    }
    return default_value;
}

static
char* get_string_arg(int argc, char*argv[], char* varname, char* default_value) {
    for (int k=1; k<argc; k++) {
        char* sep = strchr(argv[k], '=');
        if (sep) {
            *sep ='\0';
            if (strcmp(argv[k], varname)==0) {
                *sep = '=';
                return sep+1;
            }
            *sep = '=';
        }   
    }
    return default_value;
}

static char default_generators[57] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

int main(int argc, char*argv[]) {

#ifdef __EMSCRIPTEN__
    atexit(call_onMainFinished);
#endif

#ifdef USE_INT128_T
    /* Note: N>30 almost certainly causes overflow */
    size_t N = get_arg(argc, argv, "N", 5, 1, 30);
#else
    size_t N = get_arg(argc, argv, "N", 5, 1, 16);
#endif 
    set_verbosity_level(get_arg(argc, argv ,"verbosity_level", 0, 0, 9));

    if (get_arg(argc, argv, "goldberg_coefficients", 0, 0, 1)) {
        goldberg_t *G = goldberg(N);
        print_goldberg(G);
        free_goldberg(G);
        return EXIT_SUCCESS;
    }

    int basis = get_arg(argc, argv, "basis", 0, 0, 5);
    int table_output = get_arg(argc, argv, "table_output", -1, 0, 1);
    char *generators = get_string_arg(argc, argv, "generators", default_generators);

    expr_t *A = generator(0);
    expr_t *B = generator(1);
    expr_t *C = generator(2);
    expr_t *ex = NULL;
    lie_series_t *LS = NULL;

    char *inp = get_string_arg(argc, argv, "expression", NULL);
    if ((inp==0) || (strlen(inp)==1 && inp[0]>='0' && inp[0]<='9')) {
      /* predefined expressions */
        switch(get_arg(argc, argv, "expression", 0, 0, 8)) {
        case 0:  /* log(exp(A)*exp(B)), with optimizations spezific for this expression */ 
            LS = BCH(N, basis);
            break;
        case 1: 
            LS = symBCH(N, basis); /* log(exp(A/2)*exp(B)*exp(A/2), with optimizations spezific for this expression */  
            break;
        case 2: /* log(exp(A)*exp(B)*exp(A)) */
            ex = logarithm(product(product(exponential(A), exponential(B)), 
                                   exponential(A)));
            LS = lie_series(2, ex, N, basis);
            break;
        case 3: /* log(exp(A)*exp(B)*exp(C)), 3 generators */
            ex = logarithm(product(product(exponential(A), exponential(B)), exponential(C)));
            LS = lie_series(3, ex, N, basis);
            break;
        case 4: /* log(exp(A)*exp(B)*exp(-A)*exp(-B)) */
            ex = logarithm(product(product(exponential(A), exponential(B)),
                           product(exponential(negation(A)),exponential(negation(B)))));
            LS = lie_series(2, ex, N, basis);
            break;
        case 5: /* log(exp(B/6)*exp(A/2)*exp(2/3*B+1/72*[B,[A,B]])*exp(A/2)*exp(B/6)) */
            ex = logarithm(product(product(product(product(
                    exponential(term(1, 6, B)), exponential(term(1, 2, A))),
                    exponential(sum(term(2, 3, B), term(1, 72, commutator(B, commutator(A, B)))))), 
                    exponential(term(1, 2, A))), exponential(term(1, 6, B))));
            LS = lie_series(2, ex, N, basis); 
            break;
        case 6: /* log(exp(A)*exp(B)) computed in Lie algebra over 3 generators */
            ex = logarithm(product(exponential(A), exponential(B)));
            LS = lie_series(3, ex, N, basis); /* SIC! K=3 */
            break;
        case 7: /* same as case 0 but without specific optimizations */
            ex = logarithm(product(exponential(A), exponential(B)));
            LS = lie_series(2, ex, N, basis); 
            break;
        case 8: /* same as case 1 but without specific optimizations */
            ex = logarithm(product(product(exponential(term(1, 2, A)), exponential(B)), 
                                   exponential(term(1, 2, A))));
            LS = lie_series(2, ex, N, basis); 
            break;
        }
    }
    else {
        /* parse expression */
        int num_generators;
        ex = parse(inp, default_generators, &num_generators);
        if (ex==0) {
            /* parse error, function parse() should have already printed 
             * an error message
             */
            exit(EXIT_FAILURE);
        }
        print_expr(ex, default_generators);
        printf("\n");
        if (!is_lie_element(ex)) {
            fprintf(stderr, "ERROR: expression is not a Lie element\n");
            exit(EXIT_FAILURE);
        }

        LS = lie_series(num_generators, ex, N, basis);
    }

    if (strlen(generators)<LS->K) {
        fprintf(stderr, "WARNING: expected generators of length %i, got \"%s\"\n", LS->K, generators);
        generators = default_generators;
    }

    if (get_verbosity_level()>0) {
        print_statistics(LS);
        print_statistics_n(LS, N);
    }

    if (table_output==-1) {
        table_output = LS->dim<=200 ? 0 : 1;
    }

    /* output result: */
    switch(table_output) {
        case 0:
            print_lie_series(LS, generators);
            printf("\n");
            break;
        case 1: {
            int what = 
                PRINT_INDEX           *get_arg(argc, argv, "print_index",         1, 0, 1) |
                PRINT_DEGREE          *get_arg(argc, argv, "print_degree",        1, 0, 1) |
                PRINT_MULTI_DEGREE    *get_arg(argc, argv, "print_multi_degree",  0, 0, 1) |
                PRINT_FACTORS         *get_arg(argc, argv, "print_factors",       1, 0, 1) |
                PRINT_FOLIAGE         *get_arg(argc, argv, "print_foliage",       0, 0, 1) |
                PRINT_BASIS_ELEMENT   *get_arg(argc, argv, "print_basis_element", 0, 0, 1) |
                PRINT_COEFFICIENT     *get_arg(argc, argv, "print_coefficient",   1, 0, 1); 
            print_table(LS, what, generators);
            break;
        }
    }

    free_lie_series(LS);
    free_all_expressions();

    return EXIT_SUCCESS;
}
