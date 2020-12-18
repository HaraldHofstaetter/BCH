#include"bch.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static uint8_t K;             /* number of generators */
static uint8_t N;             /* maximum length of Lyndon words (=maximum order of Lie series expansion) */
static generator_t **W=NULL; /* W[i] ... ith Lyndon word, ordered primarily by length and 
                                secondarily by lexicographical order */
static generator_t **R=NULL; /* R[i] ... ith rightnormed basis element corresponding to ith Lyndon word */
static uint32_t *p1=NULL;    /* standard factorization of W[i] is W[p1[i]]*W[p2[i]] */
static uint32_t *p2=NULL;
static uint8_t  *nn=NULL;    /* nn[i] = length of W[i] */
static uint32_t *ii=NULL;    /* W[ii[n-1]] = first Lyndon word of length n; 
                                W[ii[n]-1] = last Lyndon word of length n */
static uint32_t *WI=NULL;    /* WI[i] = word index of W[i] */
static uint32_t *DI=NULL;    /* DI[i] = multi degree index of W[i] */

static size_t N_LYNDON;      /* number of Lyndon words of length <=N, N_LYNDON = ii[N] */

static unsigned int VERBOSITY_LEVEL = 0;


double tic(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timespec tt;
    clock_gettime(CLOCK_MONOTONIC, &tt);	
    return tt.tv_sec + ((double) tt.tv_nsec)*1e-9;
#endif
}

double toc(double t0) {
#ifdef _OPENMP
    double t1 = omp_get_wtime();
#else
    struct timespec tt;
    clock_gettime(CLOCK_MONOTONIC, &tt);	
    double t1 = tt.tv_sec + ((double) tt.tv_nsec)*1e-9;
#endif
    return t1-t0;
}


/*********************************************/
#include"lyndon.c"
/*********************************************/


static void compute_goldberg_coefficients(int N, INTEGER c[], INTEGER denom) {
    if (VERBOSITY_LEVEL>=1) {
        printf("#expression=log(exp(A)*exp(A))\n");
        printf("#denominator="); print_INTEGER(denom); printf("\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    double t0 = tic();
    
    goldberg_t G = goldberg(N);        
    int f = denom/G.denom;
    if (f!=1) {
        #pragma omp for schedule(dynamic,256) 
        for (int i=ii[N]-1; i>=0; i--) {
            c[i] = f*goldberg_coefficient(nn[i], W[i], &G); 
        }       
    }
    else {
        #pragma omp for schedule(dynamic,256) 
        for (int i=ii[N]-1; i>=0; i--) {
            c[i] = goldberg_coefficient(nn[i], W[i], &G); 
        }       
    }

    free_goldberg(G);

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#compute coeffs of words: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


static void compute_word_coefficients(int N, expr_t* ex, INTEGER c[], INTEGER denom) {
    if (VERBOSITY_LEVEL>=1) {
        printf("#expression="); print_expr(ex); printf("\n"); 
        printf("#denominator="); print_INTEGER(denom); printf("\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    double t0 = tic();

    size_t i1 = ii[N-1];
    size_t i2 = ii[N]-1;

    INTEGER e[N+1];

    /* c[0] needs special handling */
    INTEGER t1[2];
    e[0] = 0;
    e[1] = denom;
    int  m = phi(t1, 2, W[0], ex, e);
    c[0] = m>0 ? t1[0] : 0;

    /* now the other coeffs */
    for (int j=0; j<N; j++){
        e[j] = 0;
    }
    e[N] = denom;

    #pragma omp parallel 
    {

    size_t JW[N];
    INTEGER t[N+1];

    #pragma omp for schedule(dynamic,256) 
    for (int i=i1; i<=i2; i++) {
            generator_t *w = W[i];
            int m = phi(t, N+1, w, ex, e);
            size_t kW = get_right_factors(i, JW, N);
            for (int k=0; k<=kW; k++) {
                c[JW[k]] = k<m ? t[k] : 0;
            }
    }
    }

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#compute coeffs of words: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}    

/*********************************************/
#include"convert_lyndon.c"
#include"convert_rightnormed.c"
/*********************************************/



static void init_all(size_t number_of_generators, size_t order, int rightnormed) {
    K = number_of_generators;
    N = order;
    init_lyndon_words(rightnormed);
}


static void free_all(void) {
    free_lyndon_words();
}


static lie_series_t gen_result(INTEGER *c, INTEGER denom) {
    lie_series_t LS;
    LS.K = K;
    LS.N = N;
    LS.dim = N_LYNDON;
    LS.p1 = p1;
    LS.p2 = p2;
    LS.nn = nn;
    LS.R = R;
    LS.denom = denom;
    LS.c = c;
    return LS;
}


lie_series_t lie_series(size_t K, expr_t* expr, size_t N, int rightnormed) {
    double t0 = tic();
    init_all(K, N, rightnormed);
    INTEGER *c = malloc(N_LYNDON*sizeof(INTEGER));
    INTEGER denom = common_denominator(N, expr);
    compute_word_coefficients(N, expr, c, denom);
    if (rightnormed) {
        convert_to_rightnormed_lie_series(N, c, 0);
    }
    else {
        convert_to_lie_series(N, c);
    }
    lie_series_t LS = gen_result(c, denom);
    free_all();
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


lie_series_t BCH(size_t N, int rightnormed) {
    double t0 = tic();
    init_all(2, N,rightnormed);
    INTEGER *c = malloc(N_LYNDON*sizeof(INTEGER));
    INTEGER denom = common_denominator(N, 0);
    if (rightnormed) {
        compute_goldberg_coefficients(N, c, denom);
        convert_to_rightnormed_lie_series(N, c, 1);
        compute_rightnormed_BCH_terms_of_even_orders(c);
    }
    else {
        if (N%2) {
            compute_goldberg_coefficients(N, c, denom);
            convert_to_lie_series(N, c);
        }
        else {
            compute_goldberg_coefficients(N-1, c, denom);
            convert_to_lie_series(N-1, c);
            compute_BCH_terms_of_even_order_N(c);
        }
    }
    lie_series_t LS = gen_result(c, denom);
    free_all();
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


lie_series_t symBCH(size_t N, int rightnormed) {
    double t0 = tic();
    expr_t *halfA = generator(0);
    expr_t *B = generator(1);
    expr_t *expr = logarithm(product(product(exponential(halfA), exponential(B)), 
                                     exponential(halfA)));
    init_all(2, N, rightnormed);
    INTEGER *c = calloc(N_LYNDON, sizeof(INTEGER)); /* calloc initializes to zero */
    INTEGER denom = common_denominator(N, 0);
    if (VERBOSITY_LEVEL>=1) {
        printf("#NOTE: in the following expression, A stands for A/2\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    int N1 = N%2 ? N : N-1;
    compute_word_coefficients(N1, expr, c, denom);
    if (rightnormed) {
        convert_to_rightnormed_lie_series(N1, c, 1);
    }
    else {
        convert_to_lie_series(N1, c);
    }
    lie_series_t LS = gen_result(c, denom);
    for (int i=0; i<N_LYNDON; i++) {
        int nA = get_degree_of_generator(&LS, i, 0);
        LS.c[i] <<= N-1-nA; /* c[i] = c[i]*2^(N-1-nA) */
    }
    LS.denom <<= N-1; /* denom = denom*2^(N-1) */
    if (VERBOSITY_LEVEL>=1) {
        printf("#denominator changed to "); print_INTEGER(LS.denom); printf("\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    free_all();
    free_expr(halfA);
    free_expr(B);
    free_expr(expr);
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


void free_lie_series(lie_series_t LS) {
    free(LS.p1);
    free(LS.p2);
    free(LS.nn);
    free(LS.R);
    free(LS.c);
}


void set_verbosity_level(unsigned int level) {
    VERBOSITY_LEVEL = level;
}

unsigned int get_verbosity_level(void) {
    return VERBOSITY_LEVEL;
}


void print_lyndon_word(lie_series_t *LS,  size_t i) {
    if (i<LS->K) {
        printf("%c", (char) ('A'+i));
    }
    else {
        print_lyndon_word(LS, LS->p1[i]);
        print_lyndon_word(LS, LS->p2[i]);
    }
}   

void print_rightnormed_word(lie_series_t *LS,  size_t i) {
    if (LS->R) {
        for (int j=0; j < LS->nn[i]; j++) {
            printf("%c", (char) ('A'+LS->R[i][j]));
        }
    }
}

void print_basis_element(lie_series_t *LS,  size_t i) {
    if (i<LS->K) {
        printf("%c", (char) ('A'+i));
    }
    else {
        if (LS->R) { /* rightnormed basis element */
            for (int j=0; j < LS->nn[i]-1; j++) {
                printf("[%c,", (char) ('A'+LS->R[i][j]));
            }
            printf("%c", (char) ('A'+LS->R[i][LS->nn[i]-1]));
            for (int j=0; j < LS->nn[i]-1; j++) {
                printf("]");
            }
        }
        else { /* Lyndon basis element */
            printf("[");
            print_basis_element(LS, LS->p1[i]);
            printf(",");
            print_basis_element(LS, LS->p2[i]);
            printf("]");
        }
    }
}

void print_lie_series(lie_series_t *LS) {
    for (int i=0; i<LS->dim; i++) {
        if (LS->c[i]!=0) {
            if (LS->c[i]>0) {
                printf("+");
            }
            print_RATIONAL(LS->c[i], LS->denom);
            printf("*");
            print_basis_element(LS, i);
        }
    }
}


void print_lie_series_statistics(lie_series_t *LS) {
    int n = 1;
    int dim_n = 0;
    int dim = 0;
    int nonzero_n = 0;
    int nonzero = 0;
    printf("# degree         dim    #nonzero   dim(cum.)   #nz(cum.)\n");
    for (int i=0; i<LS->dim; i++) {
        int nn = get_degree(LS,i);
        if (nn > n) { 
            dim += dim_n;
            nonzero += nonzero_n;
            printf("#  %5i  %10i  %10i  %10i  %10i\n", n, dim_n, nonzero_n, dim, nonzero);
            dim_n = 0;
            nonzero_n = 0;
            n++;
        }
        dim_n++;
        if (LS->c[i]!=0) {
            nonzero_n++;
        }
    }
    dim += dim_n;
    nonzero += nonzero_n;
    printf("#  %5i  %10i  %10i  %10i  %10i\n", n, dim_n, nonzero_n, dim, nonzero);
    printf("#\n");
}

int get_degree(lie_series_t *LS, size_t i) {
    return LS->nn[i];
/*    
    if (i<LS->K) {
        return 1;
    }
    else {
        return get_degree(LS, LS->p1[i])+get_degree(LS, LS->p2[i]);
    }
*/
}

int get_degree_of_generator(lie_series_t *LS, size_t i, uint8_t g) {
    if (i<LS->K) {
        return i==g ? 1 : 0;
    }
    else {
        return get_degree_of_generator(LS, LS->p1[i], g)
              +get_degree_of_generator(LS, LS->p2[i], g);
    }
}

void print_lists(lie_series_t *LS, unsigned int what) {
    if (VERBOSITY_LEVEL>=1) {
        printf("# ");
        if (what & PRINT_INDEX) printf("i");
        if (what & PRINT_DEGREE) printf("\t|i|");
        if (what & PRINT_MULTI_DEGREE) printf("\tmulti degree"); 
        if (what & PRINT_FACTORS) printf("\ti'\ti\"");
        if (what & PRINT_LYNDON_WORD) printf("\tLyndon word");
        if ((LS->R) && (what & PRINT_RIGHTNORMED_WORD)) printf("\trightnormed word");
        if (what & PRINT_BASIS_ELEMENT) printf("\tbasis element");
        if (what & PRINT_COEFFICIENT) printf("\tcoefficient"); 
        printf("\n");
    }
    for (int i=0; i<LS->dim; i++) {
        if (what & PRINT_INDEX) printf("%i", i);
        if (what & PRINT_DEGREE) printf("\t%i", get_degree(LS, i));
        if (what & PRINT_MULTI_DEGREE) {
            printf("\t(%i", get_degree_of_generator(LS, i, 0));
            for (int g=1; g<LS->K; g++) {
                printf(",%i", get_degree_of_generator(LS, i, g));
            }
            printf(")");
        }
        if (what & PRINT_FACTORS) printf("\t%i\t%i", LS->p1[i], LS->p2[i]);
        if (what & PRINT_LYNDON_WORD) {
            printf("\t");
            print_lyndon_word(LS, i);
        }
        if ((LS->R) && (what & PRINT_RIGHTNORMED_WORD)) {
            printf("\t");
            print_rightnormed_word(LS, i);
        }
        if (what & PRINT_BASIS_ELEMENT) {
            printf("\t");
            print_basis_element(LS, i);

        }
        if (what & PRINT_COEFFICIENT) {
            printf("\t");
            print_RATIONAL(LS->c[i], LS->denom);
        }
        printf("\n");
    }
}



