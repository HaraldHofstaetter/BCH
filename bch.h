#ifndef BCH_H
#define BCH_H

#include<stdint.h>
#include<stddef.h>

#define USE_INT128_T 1

#ifdef USE_INT128_T
typedef __int128_t INTEGER; 
#else
typedef int64_t INTEGER;
#endif

enum expr_type { UNDEFINED, IDENTITY, GENERATOR, SUM, DIFFERENCE, PRODUCT, 
                 NEGATION, TERM, EXPONENTIAL, LOGARITHM };

typedef struct expr_t {
    enum expr_type type;
    struct expr_t *arg1;
    struct expr_t *arg2;
    int num;
    int den;
} expr_t;

expr_t* identity(void);
expr_t* generator(uint8_t n);
expr_t* sum(expr_t* arg1, expr_t* arg2);
expr_t* difference(expr_t* arg1, expr_t* arg2);
expr_t* product(expr_t* arg1, expr_t* arg2);
expr_t* negation(expr_t* arg);
expr_t* term(int num, int den, expr_t* arg);
expr_t* exponential(expr_t* arg);
expr_t* logarithm(expr_t* arg);
expr_t* commutator(expr_t* arg1, expr_t* arg2);

void print_expr(expr_t* ex);
void free_expr(expr_t* ex);

int phi(INTEGER y[], int m, uint8_t w[], expr_t* ex, INTEGER v[]);
INTEGER common_denominator(int n, expr_t* ex);
void print_INTEGER(INTEGER x);
void print_RATIONAL(INTEGER p, INTEGER q);


typedef struct goldberg_t {
    size_t N;
    size_t n_partitions;
    uint8_t **P;
    INTEGER denom;
    INTEGER *c;
} goldberg_t;

goldberg_t goldberg(size_t N);
INTEGER goldberg_coefficient(int n, uint8_t w[], goldberg_t *G);
void print_goldberg(goldberg_t *G);
void free_goldberg(goldberg_t G);


typedef struct lie_series_t {
    uint8_t K;     /* number of generators */
    uint8_t N;     /* maximum length of Lyndon words (=maximum order of Lie series expansion) */
    uint32_t dim;  /* number of Lyndon words of length<=N */
    uint8_t **W;   /* W[i] ... ith Lyndon word, ordered primarily by length and 
                               secondarily by lexicographical order */
    uint32_t *p1;  /* standard factorization of W[i] is W[p1[i]]*W[p2[i]] */
    uint32_t *p2;
    uint8_t *nn;   /* nn[i] = length of W[i] */
    uint32_t *ii;  /* W[ii[n-1]] = first Lyndon word of length n; 
                      W[ii[n]-1] = last Lyndon word of length n; 
                      ii[N] = dim */
    uint8_t **R;   /* R[i] ... ith rightnormed basis element corresponding to ith Lyndon word */
    INTEGER denom;
    INTEGER *c;
} lie_series_t;

lie_series_t lie_series(size_t K, expr_t* expr, size_t N, int rightnormed);
lie_series_t BCH(size_t N, int rightnormed);
lie_series_t symBCH(size_t N, int rightnormed);

void set_verbosity_level(unsigned int verbosity_level);
unsigned int get_verbosity_level(void);


void print_lie_series(lie_series_t *LS);
void print_lie_series_statistics(lie_series_t *LS);

enum {
    PRINT_INDEX =            1 << 0, 
    PRINT_DEGREE =           1 << 1, 
    PRINT_MULTI_DEGREE =     1 << 2, 
    PRINT_FACTORS =          1 << 3, 
    PRINT_LYNDON_WORD =      1 << 4, 
    PRINT_RIGHTNORMED_WORD = 1 << 5, 
    PRINT_BASIS_ELEMENT =    1 << 6, 
    PRINT_COEFFICIENT =      1 << 7
};


void print_lists(lie_series_t *LS, unsigned int what);
int get_degree(lie_series_t *LS, size_t i);
int get_degree_of_generator(lie_series_t *LS, size_t i, uint8_t g);
void print_lyndon_word(lie_series_t *LS,  size_t i);
void print_rightnormed_word(lie_series_t *LS,  size_t i);
void print_basis_element(lie_series_t *LS,  size_t i);

void free_lie_series(lie_series_t LS);

void lyndon2rightnormed(int lw, uint8_t w[], uint8_t r[]);


/**********************************************/
double tic(void); 
double toc(double t0);
size_t get_right_factors(size_t i, size_t J[], size_t kmax, uint32_t *p1, uint32_t *p2);
void init_lyndon_words(lie_series_t *LS);
uint32_t* multi_degree_indices(size_t K, size_t dim,  uint8_t **W, uint8_t *nn);
void convert_to_lie_series(lie_series_t *LS, int N);
void compute_BCH_terms_of_even_order_N(lie_series_t *LS);
void init_rightnormed(lie_series_t *LS);
void compute_rightnormed_BCH_terms_of_even_orders(lie_series_t *LS);
void convert_to_rightnormed_lie_series(lie_series_t *LS, int N, int odd_orders_only);

#endif /*BCH_H */
