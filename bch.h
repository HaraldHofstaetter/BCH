#ifndef BCH_H
#define BCH_H

#include<stdint.h>
#include<stddef.h>

#ifdef USE_INT128_T
typedef __int128_t INTEGER; 
#else
typedef int64_t INTEGER;
#endif

typedef uint64_t generator_t;

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
expr_t* generator(generator_t n);
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

void init_factorial(int n);
void free_factorial(void);
int phi(INTEGER y[], int m, generator_t w[], expr_t* ex, INTEGER v[]);
INTEGER common_denominator(int n);
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
INTEGER goldberg_coefficient(int n, generator_t w[], goldberg_t *G);
void print_goldberg(goldberg_t *G);
void free_goldberg(goldberg_t G);


typedef struct lie_series_t {
    size_t K;
    size_t N;
    size_t n_lyndon;
    uint32_t *p1;
    uint32_t *p2;
    INTEGER denom;
    INTEGER *c;
} lie_series_t;

lie_series_t lie_series(size_t K, expr_t* expr, size_t N, int64_t fac, size_t M);
lie_series_t BCH(size_t N, size_t M);
lie_series_t symBCH(size_t N, size_t M);

void set_verbosity_level(unsigned int verbosity_level);
unsigned int get_verbosity_level(void);
void set_max_lookup_size(size_t M);
double tic(void); 
double toc(double t0);


void print_lie_series(lie_series_t *LS);

enum {
    PRINT_INDEX =         1 << 0, 
    PRINT_DEGREE =        1 << 1, 
    PRINT_MULTI_DEGREE =  1 << 2, 
    PRINT_FACTORS =       1 << 3, 
    PRINT_WORD =          1 << 4, 
    PRINT_BASIS_ELEMENT = 1 << 5, 
    PRINT_COEFFICIENT =   1 << 6
};


void print_lists(lie_series_t *LS, unsigned int what);
int get_degree(lie_series_t *LS, size_t i);
int get_degree_of_generator(lie_series_t *LS, size_t i, uint8_t g);
void print_word(lie_series_t *LS,  size_t i);
void print_basis_element(lie_series_t *LS,  size_t i);

void free_lie_series(lie_series_t LS);

#endif /*BCH_H */
