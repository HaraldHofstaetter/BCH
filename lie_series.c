#include"bch.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

unsigned int VERBOSITY_LEVEL = 0;

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


size_t get_right_factors(size_t i, size_t J[], size_t kmax, uint32_t *p1, uint32_t *p2) {
    size_t k = 0;
    J[0] = i;
    size_t l = i;
    while ((k<kmax) && (p1[l]==0)) {
        k++;
        l = p2[l];
        J[k] = l;
    }
    return k;
}

static int degree_of_generator(lie_series_t *LS, size_t i, uint8_t g) {
    if (LS->nn[i]==1) {
        return i==g ? 1 : 0;
    }
    else {
        return degree_of_generator(LS, LS->p1[i], g)
              +degree_of_generator(LS, LS->p2[i], g);
    }
}


static void compute_goldberg_coefficients(lie_series_t *LS, int N) {
    /* computes coefficients of Lyndon words LS->W[] up to length N<=LS->N
     * in log(exp(A)exp(B)) and stores them in LS->c[]
     */
    if (VERBOSITY_LEVEL>=1) {
        printf("#expression=log(exp(A)*exp(A))\n");
        printf("#denominator="); print_INTEGER(LS->denom); printf("\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    double t0 = tic();
    
    goldberg_t G = goldberg(N);        
    int f = LS->denom/G.denom;
    if (f!=1) {
        #pragma omp for schedule(dynamic,256) 
        for (int i=LS->ii[N]-1; i>=0; i--) {
            LS->c[i] = f*goldberg_coefficient(LS->nn[i], LS->W[i], &G); 
        }       
    }
    else {
        #pragma omp for schedule(dynamic,256) 
        for (int i=LS->ii[N]-1; i>=0; i--) {
            LS->c[i] = goldberg_coefficient(LS->nn[i], LS->W[i], &G); 
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


static void compute_word_coefficients(lie_series_t *LS, int N, expr_t* ex) {
    /* computes coefficients of Lyndon words LS->W[] up to length N<=LS->N
     * in expression ex and stores them in LS->c[]
     */
    if (VERBOSITY_LEVEL>=1) {
        printf("#expression="); print_expr(ex); printf("\n"); 
        printf("#denominator="); print_INTEGER(LS->denom); printf("\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    double t0 = tic();

    size_t i1 = LS->ii[N-1];
    size_t i2 = LS->ii[N]-1;

    INTEGER e[N+1];

    /* c[0] needs special handling */
    INTEGER t1[2];
    e[0] = 0;
    e[1] = LS->denom;
    int  m = phi(t1, 2, LS->W[0], ex, e);
    LS->c[0] = m>0 ? t1[0] : 0;

    /* now the other coeffs */
    for (int j=0; j<N; j++){
        e[j] = 0;
    }
    e[N] = LS->denom;

    #pragma omp parallel 
    {

    size_t JW[N];
    INTEGER t[N+1];

    #pragma omp for schedule(dynamic,256) 
    for (int i=i1; i<=i2; i++) {
            uint8_t *w = LS->W[i];
            int m = phi(t, N+1, w, ex, e);
            size_t kW = get_right_factors(i, JW, N, LS->p1, LS->p2);
            for (int k=0; k<=kW; k++) {
                LS->c[JW[k]] = k<m ? t[k] : 0;
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



lie_series_t lie_series(size_t K, expr_t* expr, size_t N, int basis) {
    double t0 = tic();
    lie_series_t LS;
    LS.K = K;
    LS.N = N;
    init_lyndon_words(&LS);
    LS.c = malloc(LS.dim*sizeof(INTEGER));
    LS.denom = common_denominator(N, expr);
    compute_word_coefficients(&LS, N, expr);
    if (basis==RIGHTNORMED_BASIS) {
        init_rightnormed(&LS);
        convert_to_rightnormed_lie_series(&LS, N, 0);
    }
    else {
        convert_to_lie_series(&LS, N);
        LS.R = 0;
        if (basis>=HALL_BASIS) {
            lie_series_t HS;
            convert_lyndon_to_hall_lie_series(&LS, &HS, basis);
            free_lie_series(LS);
            LS = HS;
        }
    }
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


lie_series_t BCH(size_t N, int basis) {
    double t0 = tic();
    lie_series_t LS;
    LS.K = 2;
    LS.N = N;
    init_lyndon_words(&LS);
    LS.c = malloc(LS.dim*sizeof(INTEGER));
    LS.denom = common_denominator(N, 0);
    if (basis==RIGHTNORMED_BASIS) {
        init_rightnormed(&LS);
        compute_goldberg_coefficients(&LS, N);
        convert_to_rightnormed_lie_series(&LS, N, 1);
        compute_rightnormed_BCH_terms_of_even_orders(&LS);
    }
    else {
        if (N%2) {
            compute_goldberg_coefficients(&LS, N);
            convert_to_lie_series(&LS, N);
        }
        else {
            compute_goldberg_coefficients(&LS, N-1);
            convert_to_lie_series(&LS, N-1);
            compute_BCH_terms_of_even_order_N(&LS);
        }
        if (basis>=HALL_BASIS) {
            lie_series_t HS;
            convert_lyndon_to_hall_lie_series(&LS, &HS, basis);
            free_lie_series(LS);
            LS = HS;
        }
    }
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


lie_series_t symBCH(size_t N, int basis) {
    double t0 = tic();
    lie_series_t LS;
    LS.K = 2;
    LS.N = N;
    expr_t *halfA = generator(0);
    expr_t *B = generator(1);
    expr_t *expr = logarithm(product(product(exponential(halfA), exponential(B)), 
                                     exponential(halfA)));
    init_lyndon_words(&LS);
    LS.c = calloc(LS.dim, sizeof(INTEGER)); /* calloc initializes to zero */
    LS.denom = common_denominator(N, 0);
    if (VERBOSITY_LEVEL>=1) {
        printf("#NOTE: in the following expression, A stands for A/2\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    int N1 = N%2 ? N : N-1;
    compute_word_coefficients(&LS, N1, expr);
    if (basis==RIGHTNORMED_BASIS) {
        init_rightnormed(&LS);
        convert_to_rightnormed_lie_series(&LS, N1, 1);
    }
    else {
        convert_to_lie_series(&LS, N1);
        if (basis>=HALL_BASIS) {
            lie_series_t HS;
            convert_lyndon_to_hall_lie_series(&LS, &HS, basis);
            free_lie_series(LS);
            LS = HS;
        }
    }
    for (int i=0; i<LS.dim; i++) {
        int nA = degree_of_generator(&LS, i, 0);
        LS.c[i] <<= N-1-nA; /* c[i] = c[i]*2^(N-1-nA) */
    }
    LS.denom <<= N-1; /* denom = denom*2^(N-1) */
    if (VERBOSITY_LEVEL>=1) {
        printf("#denominator changed to "); print_INTEGER(LS.denom); printf("\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
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
    if (LS.W) {
        free(LS.W[0]);
        free(LS.W);
    }
    free(LS.nn);
    free(LS.p1);
    free(LS.p2);
    free(LS.ii);
    if (LS.R) {
        free(LS.R[0]);
        free(LS.R);
    }
    free(LS.c);
}


void set_verbosity_level(unsigned int level) {
    VERBOSITY_LEVEL = level;
}

unsigned int get_verbosity_level(void) {
    return VERBOSITY_LEVEL;
}


void print_lyndon_word(lie_series_t *LS,  size_t i, char *g) {
    if (LS->nn[i]==1) {
        printf("%c", g[i]);
    }
    else {
        print_lyndon_word(LS, LS->p1[i], g);
        print_lyndon_word(LS, LS->p2[i], g);
    }
}   

void print_rightnormed_word(lie_series_t *LS,  size_t i, char *g) {
    if (LS->R) {
        for (int j=0; j < LS->nn[i]; j++) {
            printf("%c", g[LS->R[i][j]]);
        }
    }
}

void print_basis_element(lie_series_t *LS,  size_t i, char *g) {
    if (LS->nn[i]==1) {
        printf("%c", g[i]);
    }
    else {
        if (LS->R) { /* rightnormed basis element */
            for (int j=0; j < LS->nn[i]-1; j++) {
                printf("[%c,", g[LS->R[i][j]]);
            }
            printf("%c", g[LS->R[i][LS->nn[i]-1]]);
            for (int j=0; j < LS->nn[i]-1; j++) {
                printf("]");
            }
        }
        else { /* Lyndon or Hall basis element */
            printf("[");
            print_basis_element(LS, LS->p1[i], g);
            printf(",");
            print_basis_element(LS, LS->p2[i], g);
            printf("]");
        }
    }
}

void print_lie_series(lie_series_t *LS, char *g) {
    for (int i=0; i<LS->dim; i++) {
        if (LS->c[i]!=0) {
            if (LS->c[i]>0) {
                printf("+");
            }
            print_RATIONAL(LS->c[i], LS->denom);
            printf("*");
            print_basis_element(LS, i, g);
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
        if (LS->nn[i] > n) { 
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


void print_lists(lie_series_t *LS, unsigned int what, char* g) {
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
        if (what & PRINT_DEGREE) printf("\t%i", LS->nn[i]);
        if (what & PRINT_MULTI_DEGREE) {
            printf("\t(%i", degree_of_generator(LS, i, 0));
            for (int g=1; g<LS->K; g++) {
                printf(",%i", degree_of_generator(LS, i, g));
            }
            printf(")");
        }
        if (what & PRINT_FACTORS) printf("\t%i\t%i", LS->p1[i], LS->p2[i]);
        if (what & PRINT_LYNDON_WORD) {
            printf("\t");
            print_lyndon_word(LS, i, g);
        }
        if ((LS->R) && (what & PRINT_RIGHTNORMED_WORD)) {
            printf("\t");
            print_rightnormed_word(LS, i, g);
        }
        if (what & PRINT_BASIS_ELEMENT) {
            printf("\t");
            print_basis_element(LS, i, g);

        }
        if (what & PRINT_COEFFICIENT) {
            printf("\t");
            print_RATIONAL(LS->c[i], LS->denom);
        }
        printf("\n");
    }
}



