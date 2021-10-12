#ifndef BCH_H
#define BCH_H

#include<stdint.h>
#include<stddef.h>

#define USE_INT128_T 1
//#define USE_QUADMATH 1
#define USE_DOUBLEDOUBLE 1

#ifdef USE_INT128_T
typedef __int128_t INTEGER; 
#else
typedef int64_t INTEGER;
#endif

#if defined(USE_QUADMATH)
#include <quadmath.h>
typedef __float128 FLOAT;
#elif defined(USE_DOUBLEDOUBLE)
typedef struct {
    double hi;
    double lo;
} FLOAT;
#else
typedef double FLOAT;
#endif

int str_INTEGER(char *out, INTEGER x);
int str_RATIONAL(char *out, INTEGER p, INTEGER q);
void print_INTEGER(INTEGER x);
void print_RATIONAL(INTEGER p, INTEGER q);
INTEGER gcd_INTEGER(INTEGER a, INTEGER b);


typedef struct lie_series_t {
    uint8_t K;     /* number of generators */
    uint8_t N;     /* maximum degree of Lie series */ 
    uint32_t dim;  /* dimension of Lie algebra */ 
    uint8_t **W;   /* W[i], 0<=i<dim ... ith Lyndon word, ordered primarily by length and 
                                         secondarily by lexicographical order */
    uint8_t *nn;   /* nn[i] = length of W[i] */
    uint32_t *p1;  /* standard factorization of W[i] is W[p1[i]]*W[p2[i]] */
    uint32_t *p2;
    uint32_t *ii;  /* W[ii[n-1]] = first Lyndon word of length n; 
                      W[ii[n]-1] = last Lyndon word of length n; 
                      ii[N] = dim */
    uint8_t **R;   /* R[i] ... ith rightnormed basis element corresponding to ith Lyndon word */
    INTEGER denom;
    INTEGER *c;
    FLOAT *c_f;
} lie_series_t;


enum {
    LYNDON_BASIS = 0,
    RIGHTNORMED_BASIS = 1,
    HALL_BASIS = 2,
};     


lie_series_t* BCH(int N, int basis);
lie_series_t* symBCH(int N, int basis);

void free_lie_series(lie_series_t *LS);

void set_verbosity_level(int verbosity_level);
int get_verbosity_level(void);
void set_float_output_digits(int digits);
void set_float_output_threshold(double threshold);

int dimension(lie_series_t *LS);
int maximum_degree(lie_series_t *LS);
int number_of_generators(lie_series_t *LS);
INTEGER denominator(lie_series_t *LS);
INTEGER numerator_of_coefficient(lie_series_t *LS,  int i);
int degree(lie_series_t *LS, int i);
int degree_of_generator(lie_series_t *LS, int i, uint8_t g);
int left_factor(lie_series_t *LS, int i);
int right_factor(lie_series_t *LS, int i);
int str_foliage(char *out, lie_series_t *LS,  int i, char *generators);
int str_basis_element(char *out, lie_series_t *LS,  int i, char *generators);
int str_coefficient(char *out, lie_series_t *LS,  int i);
void print_foliage(lie_series_t *LS,  int i, char *generators);
void print_basis_element(lie_series_t *LS,  int i, char *generators);
void print_coefficient(lie_series_t *LS,  int i);

void print_lie_series(lie_series_t *LS, char *generators);
void print_statistics(lie_series_t *LS);
void print_statistics_n(lie_series_t *LS, int n);

enum {
    PRINT_INDEX =            1 << 0, 
    PRINT_DEGREE =           1 << 1, 
    PRINT_MULTI_DEGREE =     1 << 2, 
    PRINT_FACTORS =          1 << 3, 
    PRINT_FOLIAGE =          1 << 4, 
    PRINT_BASIS_ELEMENT =    1 << 5, 
    PRINT_COEFFICIENT =      1 << 6
};

void print_table(lie_series_t *LS, int what, char *generators);


enum expr_type { ZERO_ELEMENT, IDENTITY, GENERATOR, SUM, DIFFERENCE, PRODUCT, 
                 NEGATION, TERM, TERM_F, EXPONENTIAL, LOGARITHM };

typedef struct rat_t {
    int num;
    int den;
} rat_t;

typedef struct expr_t {
    enum expr_type type;
    struct expr_t *arg1;  /* pointer to subexpression */
    struct expr_t *arg2;  /* pointer to subexpression */
    int gen;              /* for expr_type GENERATOR */
    rat_t factor;         /* for expr_type TERM */
    FLOAT factor_f;       /* for expr_type TERM */
    int mindeg;           /* minimal degree */
} expr_t;

expr_t* zero_element(void);
expr_t* identity(void);
expr_t* identity(void);
expr_t* generator(uint8_t n);
expr_t* sum(expr_t* arg1, expr_t* arg2);
expr_t* difference(expr_t* arg1, expr_t* arg2);
expr_t* product(expr_t* arg1, expr_t* arg2);
expr_t* negation(expr_t* arg);
expr_t* term(int num, int den, expr_t* arg);
expr_t* term_f(FLOAT x, expr_t* arg);
expr_t* exponential(expr_t* arg);
expr_t* logarithm(expr_t* arg);
expr_t* commutator(expr_t* arg1, expr_t* arg2);

int str_expr(char *out, expr_t* ex, char *gens);
void print_expr(expr_t* ex, char *gens);
void free_all_expressions(void);

expr_t* parse(char *inp, char *generators, int *num_generators);

int phi(INTEGER y[], int m, uint8_t w[], expr_t* ex, INTEGER v[]);
INTEGER common_denominator(int n, expr_t* ex);
int is_lie_element(expr_t* ex);

int phi_f(FLOAT y[], int m, uint8_t w[], expr_t* ex, FLOAT v[]);
void init_phi_f(int n);
static inline int is_contaminated_with_floats(expr_t* ex) {
    return ex->factor.den==0;
}


lie_series_t* lie_series(int K, expr_t* expr, int N, int basis);


typedef struct goldberg_t {
    size_t N;
    size_t n_partitions;
    uint8_t **P;
    INTEGER denom;
    INTEGER *c;
} goldberg_t;

goldberg_t* goldberg(size_t N);
INTEGER goldberg_coefficient(int n, uint8_t w[], goldberg_t *G);
void print_goldberg(goldberg_t *G);
void free_goldberg(goldberg_t *G);


/**********************************************/
/* floating point operations: */

FLOAT generic_parse_FLOAT(char *in);

#if defined(USE_QUADMATH)

/* https://gcc.gnu.org/onlinedocs/libquadmath/ */

static inline FLOAT i2f(int x) {return ((FLOAT) x);}
static inline FLOAT i64_to_f(int64_t x) {return ((FLOAT) x);}
static inline FLOAT r2f(rat_t x) {return ((FLOAT) x.num)/((FLOAT) x.den);}
static inline FLOAT add_f(FLOAT a, FLOAT b) {return a+b;}
static inline FLOAT sub_f(FLOAT a, FLOAT b) {return a-b;}
static inline FLOAT mul_f(FLOAT a, FLOAT b) {return a*b;}
static inline FLOAT div_f(FLOAT a, FLOAT b) {return a/b;}
static inline FLOAT neg_f(FLOAT a) {return -a;}
static inline FLOAT abs_f(FLOAT x) {return fabsq(x);}

static inline FLOAT zero_f(void) { return 0.0q; }
static inline FLOAT one_f(void) { return 1.0q; }
static inline FLOAT eps_f(void) { return 1e-30q; }

static inline int is_zero_f(FLOAT x) { return fabsq(x)<eps_f(); }
static inline int is_one_f(FLOAT x)  { return fabsq(x-1.0q)<eps_f(); }
static inline int lt_f(FLOAT x, FLOAT y) { return x < y; }

#include <stdio.h>
extern char *FLOAT_OUTPUT_FORMAT;
static inline int str_FLOAT(char *out, FLOAT x) {
    /* TODO use quadmath_snprintf */
    return out==NULL ? snprintf(NULL, 0, FLOAT_OUTPUT_FORMAT, (double) x) : 
                        sprintf(out,     FLOAT_OUTPUT_FORMAT, (double) x) ; 
}

static inline void print_FLOAT(FLOAT x) { 
    char buf[128];
    str_FLOAT(buf, x);
    printf("%s", buf); 
}
static inline FLOAT parse_FLOAT(char *in) { return strtoflt128(in, NULL); }

#elif defined(USE_DOUBLEDOUBLE)
/* extended-precision arithmetic using pairs of floating-point numbers based on:
 * https://www.davidhbailey.com/dhbsoftware/ 
 * https://github.com/JuliaMath/DoubleDouble.jl/blob/master/src/DoubleDouble.jl 
 */
static inline FLOAT doubledouble(double u, double v) { 
    double w = u + v;
    FLOAT x;
    x.hi = w;
    x.lo = (u-w) + v;
    return x;
}

static inline double halfprec(double x) {
    double p = x*1.34217729e8; /* half64 */
    return (x-p)+p;
}

#include<math.h> /* fabs */

static inline FLOAT add_f(FLOAT x, FLOAT y) {
    double r = x.hi + y.hi;
    double s = fabs(x.hi) > fabs(y.hi) ? (((x.hi - r) + y.hi) + y.lo) + x.lo : 
                                         (((y.hi - r) + x.hi) + x.lo) + y.lo ;
    return doubledouble(r, s);
}

static inline FLOAT sub_f(FLOAT x, FLOAT y) {
    double r = x.hi - y.hi;
    double s = fabs(x.hi) > fabs(y.hi) ? (((x.hi - r) - y.hi) - y.lo) + x.lo  : 
                                         (((-y.hi - r) + x.hi) + x.lo) - y.lo ;
    return doubledouble(r, s);
}

static inline FLOAT mul_f(FLOAT x, FLOAT y) {
    double hx = halfprec(x.hi);
    double lx = x.hi - hx;
    double hy = halfprec(y.hi);
    double ly = y.hi - hy;
    double z = x.hi*y.hi;
    FLOAT  c = doubledouble(z, ((hx*hy-z) + hx*ly + lx*hy) + lx*ly);
    double cc = (x.hi*y.lo + x.lo*y.hi) + c.lo;
    return doubledouble(c.hi, cc);
}

static inline FLOAT div_f(FLOAT x, FLOAT y) {
    double c = x.hi/y.hi;
    double hx = halfprec(c);
    double lx = c - hx;
    double hy = halfprec(y.hi);
    double ly = y.hi - hy;
    double z = c*y.hi;
    FLOAT  u = doubledouble(z, ((hx*hy-z) + hx*ly + lx*hy) + lx*ly);
    double cc = ((((x.hi - u.hi) - u.lo) + x.lo) - c*y.lo)/y.hi;
    return doubledouble(c, cc);
}

static inline FLOAT neg_f(FLOAT x) {return doubledouble(-x.hi, -x.lo); }
static inline FLOAT abs_f(FLOAT x) { return x.hi>0 ? x : neg_f(x); }

static inline FLOAT i2f(int x) {return doubledouble((double) x, 0.0);}
static inline FLOAT i64_to_f(int64_t x) {return doubledouble((double) x, 0.0);}
static inline FLOAT d2f(double x) {return doubledouble(x, 0.0);}
static inline FLOAT r2f(rat_t x) {return div_f(i2f(x.num), i2f(x.den));}

static inline FLOAT zero_f(void) { return doubledouble(0.0, 0.0); }
static inline FLOAT one_f(void) { return doubledouble(1.0, 0.0); }
static inline FLOAT eps_f(void) { return doubledouble(1e-30, 0.0); }

static inline int lt_f(FLOAT x, FLOAT y) { return x.hi + x.lo < y.hi + y.lo; }
static inline int is_zero_f(FLOAT x) { return lt_f(abs_f(x), eps_f()); }
static inline int is_one_f(FLOAT x)  { return lt_f(abs_f(sub_f(x, i2f(1))), eps_f()); }


#include <stdio.h>
extern char *FLOAT_OUTPUT_FORMAT;
static inline int str_FLOAT(char *out, FLOAT x) { 
    return out==NULL ? snprintf(NULL, 0, FLOAT_OUTPUT_FORMAT, x.hi) : 
                        sprintf(out,     FLOAT_OUTPUT_FORMAT, x.hi) ; 
}

/* TODO: better implementation */
static inline void print_FLOAT(FLOAT x) { 
    char buf[128];
    str_FLOAT(buf, x);
    printf("%s", buf); 
}

static inline FLOAT parse_FLOAT(char *in) { return generic_parse_FLOAT(in); }

#else

static inline FLOAT i2f(int x) {return ((FLOAT) x);}
static inline FLOAT i64_to_f(int64_t x) {return ((FLOAT) x);}
static inline FLOAT r2f(rat_t x) {return ((FLOAT) x.num)/((FLOAT) x.den);}
static inline FLOAT add_f(FLOAT a, FLOAT b) {return a+b;}
static inline FLOAT sub_f(FLOAT a, FLOAT b) {return a-b;}
static inline FLOAT mul_f(FLOAT a, FLOAT b) {return a*b;}
static inline FLOAT div_f(FLOAT a, FLOAT b) {return a/b;}
static inline FLOAT neg_f(FLOAT a) {return -a;}
static inline FLOAT abs_f(FLOAT x) {return fabs(x);}

static inline FLOAT zero_f(void) { return 0.0; }
static inline FLOAT one_f(void) { return 1.0; }

#include <math.h>
static inline int is_zero_f(FLOAT x) { return fabs(x)<1e-14; }
static inline int is_one_f(FLOAT x)  { return fabs(x-1.0)<1e-14; }
static inline int lt_f(FLOAT x, FLOAT y) { return x < y; }

#include <stdio.h>
extern char *FLOAT_OUTPUT_FORMAT;
/* TODO: better implementation */
static inline int str_FLOAT(char *out, FLOAT x) { 
    return out==NULL ? snprintf(NULL, 0, FLOAT_OUTPUT_FORMAT, x) : 
                       sprintf(out,      FLOAT_OUTPUT_FORMAT, x) ; 
}

static inline void print_FLOAT(FLOAT x) { printf("%g", x); }
static inline FLOAT parse_FLOAT(char *in) { FLOAT d; sscanf(in, "%lf", &d ); return d; }

#endif

/* parse.y, expr.c: */
rat_t rat(int num, int den);
rat_t add_r(rat_t a, rat_t b);
rat_t sub_r(rat_t a, rat_t b);
rat_t mul_r(rat_t a, rat_t b);
rat_t div_r(rat_t a, rat_t b);
rat_t neg_r(rat_t a);
expr_t* term_r(rat_t factor, expr_t* arg);

/* lie_series.c: */
double tic(void); 
double toc(double t0);
size_t get_right_factors(size_t i, size_t J[], size_t kmax, uint32_t *p1, uint32_t *p2);

/* lyndon.c: */
void init_lyndon_words(lie_series_t *LS);
uint32_t* multi_degree_indices(size_t K, size_t dim,  uint8_t **W, uint8_t *nn);

/* convert_lyndon.c: */ 
void convert_to_lie_series(lie_series_t *LS, int N);
void compute_BCH_terms_of_even_degree_N(lie_series_t *LS);

/* convert_rightnormed.c: */ 
void init_rightnormed(lie_series_t *LS);
void compute_rightnormed_BCH_terms_of_even_degrees(lie_series_t *LS);
void convert_to_rightnormed_lie_series(lie_series_t *LS, int N, int odd_degrees_only);

/* rightnormed.c : */
void lyndon2rightnormed(int lw, uint8_t w[], uint8_t r[]);

/* convert_hall.c: */
void convert_lyndon_to_hall_lie_series(lie_series_t *LS, lie_series_t *HS, int basis);


#endif /*BCH_H */
