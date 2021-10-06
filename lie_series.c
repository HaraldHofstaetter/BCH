#include"bch.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int VERBOSITY_LEVEL = 0;

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


static void compute_goldberg_coefficients(lie_series_t *LS, int N) {
    /* computes coefficients of Lyndon words LS->W[] up to length N<=LS->N
     * in log(exp(A)exp(B)) and stores them in LS->c[]
     */
    if (get_verbosity_level()>=1) {
        printf("#expression=log(exp(A)*exp(B))\n");
        printf("#denominator="); print_INTEGER(LS->denom); printf("\n");
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
    double t0 = tic();
    
    goldberg_t *G = goldberg(N);        
    int f = LS->denom/G->denom;
    if (f!=1) {
        #pragma omp for schedule(dynamic,256) 
        for (int i=LS->ii[N]-1; i>=0; i--) {
            LS->c[i] = f*goldberg_coefficient(LS->nn[i], LS->W[i], G); 
        }       
    }
    else {
        #pragma omp for schedule(dynamic,256) 
        for (int i=LS->ii[N]-1; i>=0; i--) {
            LS->c[i] = goldberg_coefficient(LS->nn[i], LS->W[i], G); 
        }       
    }

    free_goldberg(G);

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#compute coefficients of Lyndon words: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}


static void compute_word_coefficients(lie_series_t *LS, int N, expr_t* ex) {
    /* computes coefficients of Lyndon words LS->W[] up to length N<=LS->N
     * in expression ex and stores them in LS->c[]
     */
    if (get_verbosity_level()>=1) {
        printf("#expression="); print_expr(ex, 0); printf("\n"); 
        printf("#denominator="); print_INTEGER(LS->denom); printf("\n");
        if (get_verbosity_level()>=2) {
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

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#compute coeffs of words: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}    



lie_series_t* lie_series(int K, expr_t* expr, int N, int basis) {
    double t0 = tic();
    lie_series_t *LS = malloc(sizeof(lie_series_t));
    LS->K = K;
    LS->N = N;
    init_lyndon_words(LS);
    LS->c = malloc(LS->dim*sizeof(INTEGER));
    LS->denom = common_denominator(N, expr);
    compute_word_coefficients(LS, N, expr);
    if (basis==RIGHTNORMED_BASIS) {
        init_rightnormed(LS);
        convert_to_rightnormed_lie_series(LS, N, 0);
    }
    else {
        convert_to_lie_series(LS, N);
        LS->R = 0;
        if (basis>=HALL_BASIS) {
            lie_series_t *HS = malloc(sizeof(lie_series_t));
            convert_lyndon_to_hall_lie_series(LS, HS, basis);
            free_lie_series(LS);
            LS = HS;
        }
    }
    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


lie_series_t* BCH(int N, int basis) {
    double t0 = tic();
    lie_series_t *LS = malloc(sizeof(lie_series_t));
    LS->K = 2;
    LS->N = N;
    init_lyndon_words(LS);
    LS->c = malloc(LS->dim*sizeof(INTEGER));
    LS->denom = common_denominator(N, 0);
    if (basis==RIGHTNORMED_BASIS) {
        init_rightnormed(LS);
        compute_goldberg_coefficients(LS, N);
        convert_to_rightnormed_lie_series(LS, N, 1);
        compute_rightnormed_BCH_terms_of_even_degrees(LS);
    }
    else {
        if (N%2) {
            compute_goldberg_coefficients(LS, N);
            convert_to_lie_series(LS, N);
        }
        else {
            compute_goldberg_coefficients(LS, N-1);
            convert_to_lie_series(LS, N-1);
            compute_BCH_terms_of_even_degree_N(LS);
        }
        if (basis>=HALL_BASIS) {
            lie_series_t *HS = malloc(sizeof(lie_series_t));
            convert_lyndon_to_hall_lie_series(LS, HS, basis);
            free_lie_series(LS);
            LS = HS;
        }
    }
    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


lie_series_t* symBCH(int N, int basis) {
    double t0 = tic();
    lie_series_t *LS = malloc(sizeof(lie_series_t));
    LS->K = 2;
    LS->N = N;
    expr_t *halfA = generator(0);
    expr_t *B = generator(1);
    expr_t *expr = logarithm(product(product(exponential(halfA), exponential(B)), 
                                     exponential(halfA)));
    init_lyndon_words(LS);
    LS->c = calloc(LS->dim, sizeof(INTEGER)); /* calloc initializes to zero */
    LS->denom = common_denominator(N, 0);
    if (get_verbosity_level()>=1) {
        printf("#NOTE: in the following expression, A stands for A/2\n");
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
    int N1 = N%2 ? N : N-1;
    compute_word_coefficients(LS, N1, expr);
    if (basis==RIGHTNORMED_BASIS) {
        init_rightnormed(LS);
        convert_to_rightnormed_lie_series(LS, N1, 1);
    }
    else {
        convert_to_lie_series(LS, N1);
        if (basis>=HALL_BASIS) {
            lie_series_t *HS = malloc(sizeof(lie_series_t));
            convert_lyndon_to_hall_lie_series(LS, HS, basis);
            free_lie_series(LS);
            LS = HS;
        }
    }
    for (int i=0; i<LS->dim; i++) {
        int nA = degree_of_generator(LS, i, 0);
        LS->c[i] *= (1<<(N-1-nA)); /* c[i] = c[i]*2^(N-1-nA) */
    }
    LS->denom <<= N-1; /* denom = denom*2^(N-1) */
    if (get_verbosity_level()>=1) {
        printf("#denominator changed to "); print_INTEGER(LS->denom); printf("\n");
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
    free_expr(halfA);
    free_expr(B);
    free_expr(expr);
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#total time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
    return LS;
}


void free_lie_series(lie_series_t *LS) {
    if (LS->W) {
        free(LS->W[0]);
        free(LS->W);
    }
    free(LS->nn);
    free(LS->p1);
    free(LS->p2);
    free(LS->ii);
    if (LS->R) {
        free(LS->R[0]);
        free(LS->R);
    }
    free(LS->c);
    free(LS);
}


void set_verbosity_level(int level) {
    VERBOSITY_LEVEL = level;
}

int get_verbosity_level(void) {
    return VERBOSITY_LEVEL;
}


int dimension(lie_series_t *LS) {
    return LS->dim;
}

int maximum_degree(lie_series_t *LS) {
    return LS->N;
}

int number_of_generators(lie_series_t *LS) {
    return LS->K;
}

INTEGER denominator(lie_series_t *LS){
    return LS->denom;
}

INTEGER numerator_of_coefficient(lie_series_t *LS,  int i) {
    return LS->c[i];
}

int degree(lie_series_t *LS, int i) {
    return LS->nn[i];
}


int degree_of_generator(lie_series_t *LS, int i, uint8_t g) {
    if (LS->nn[i]==1) {
        return i==g ? 1 : 0;
    }
    else {
        return degree_of_generator(LS, LS->p1[i], g)
              +degree_of_generator(LS, LS->p2[i], g);
    }
}


int left_factor(lie_series_t *LS, int i) {
    return LS->p1[i];
}


int right_factor(lie_series_t *LS, int i) {
    return LS->p2[i];
}


int str_foliage(char *out, lie_series_t *LS,  int i, char *g) {
    if (!out) {
        return LS->nn[i];
    }
    int pos = 0;
    if (LS->nn[i]==1) {
        out[pos++] = g[LS->p1[i]];
    }
    else {
        pos += str_foliage(out+pos, LS, LS->p1[i], g);
        pos += str_foliage(out+pos, LS, LS->p2[i], g);
    }
    out[pos] = '\0';
    return pos;
}


void print_foliage(lie_series_t *LS,  int i, char *g) {
    int n = LS->nn[i];
    char out[n+1];
    str_foliage(out, LS, i , g);
    printf("%s", out);
}


int str_basis_element(char *out, lie_series_t *LS,  int i, char *g) {
    int n = LS->nn[i];
    if (!out) {
        return n+3*(n-1);
    }
    int pos = 0;
    if (n==1) {
        out[pos++] = g[LS->p1[i]];
    }
    else { 
        out[pos++] = '[';
        pos += str_basis_element(out+pos, LS, LS->p1[i], g);
        out[pos++] = ',';
        pos += str_basis_element(out+pos, LS, LS->p2[i], g);
        out[pos++] = ']';
    }
    out[pos] = '\0';
    return pos;
}


void print_basis_element(lie_series_t *LS,  int i, char *g) {
    int n = LS->nn[i];
    char out[n+3*(n-1)+1];
    str_basis_element(out, LS, i , g);
    printf("%s", out);
}


int str_coefficient(char *out, lie_series_t *LS,  int i) {
    return str_RATIONAL(out, LS->c[i], LS->denom);
}


void print_coefficient(lie_series_t *LS,  int i) {
    print_RATIONAL(LS->c[i], LS->denom);
}


void print_lie_series(lie_series_t *LS, char *g) {
    for (int i=0; i<dimension(LS); i++) {
        INTEGER num = numerator_of_coefficient(LS, i);
        if (num!=0) {
            if (num>0) {
                printf("+");
            }
            print_coefficient(LS, i);
            printf("*");
            print_basis_element(LS, i, g);
        }
    }
}


void print_table(lie_series_t *LS, int what, char* g) {
    if (get_verbosity_level()>=1) {
        printf("# ");
        if (what & PRINT_INDEX) printf("i");
        if (what & PRINT_DEGREE) printf("\t|i|");
        if (what & PRINT_MULTI_DEGREE) printf("\tmulti degree"); 
        if (what & PRINT_FACTORS) printf("\ti'\ti\"");
        if (what & PRINT_FOLIAGE) printf("\tfoliage");
        if (what & PRINT_BASIS_ELEMENT) printf("\tbasis element");
        if (what & PRINT_COEFFICIENT) printf("\tcoefficient"); 
        printf("\n");
    }
    for (int i=0; i<dimension(LS); i++) {
        if (what & PRINT_INDEX) printf("%i", i);
        if (what & PRINT_DEGREE) printf("\t%i", degree(LS, i));
        if (what & PRINT_MULTI_DEGREE) {
            printf("\t(%i", degree_of_generator(LS, i, 0));
            for (int g=1; g<number_of_generators(LS); g++) {
                printf(",%i", degree_of_generator(LS, i, g));
            }
            printf(")");
        }
        if (what & PRINT_FACTORS) printf("\t%i\t%i", left_factor(LS, i), right_factor(LS, i));
        if (what & PRINT_FOLIAGE) {
            printf("\t");
            print_foliage(LS, i, g);
        }
        if (what & PRINT_BASIS_ELEMENT) {
            printf("\t");
            print_basis_element(LS, i, g);

        }
        if (what & PRINT_COEFFICIENT) {
            printf("\t");
            print_coefficient(LS, i);
        }
        printf("\n");
    }
}


void print_statistics(lie_series_t *LS) {
    int N = maximum_degree(LS);
    int dim[N];
    int nonzero[N];
    for (int i=0; i<N; i++) {
        dim[i] = 0;
        nonzero[i] = 0;
    }
    for (int i=0; i<dimension(LS); i++) {
        int n = degree(LS, i);
        dim[n-1]++;
        if (numerator_of_coefficient(LS, i)!=0) {
            nonzero[n-1]++;
        }
    }
    printf("# degree         dim    #nonzero   dim(cum.)   #nz(cum.)\n");
    int dim_cum = 0;
    int nonzero_cum = 0;
    for (int n=1; n<=N; n++) {
        dim_cum += dim[n-1];
        nonzero_cum += nonzero[n-1];
        printf("#  %5i  %10i  %10i  %10i  %10i\n", n, dim[n-1], nonzero[n-1], dim_cum, nonzero_cum);
    }
    printf("#\n");
}


void print_statistics_n(lie_series_t *LS, int n) {
    assert(n<=maximum_degree(LS));
    int K = number_of_generators(LS);
    int m = 1;
    for (int k=0; k<K-1; k++) {
        m *= n;
    }
    int* dim = calloc(m, sizeof(int));
    int* nonzero = calloc(m, sizeof(int));
    for (int i=0; i<dimension(LS); i++) {
        if (degree(LS, i)==n) {
            int j = 0;
            int d = 1;
            for (int k=0; k<K-1; k++) {
                j += d*degree_of_generator(LS, i, k);
                d *= n;
            }
            dim[j]++;
            if (numerator_of_coefficient(LS, i)!=0) {
                nonzero[j]++;
            }
        }
    }
    printf("# multi-degree\tdim\t#nonzero\n");
    for (int j=0; j<m; j++) {
        if (dim[j]>0) {
        printf("# (");
        int jj = j;
        int s = 0;
        for (int k=0; k<K-1; k++) {
            int d = jj % n;
            printf("%2i,", d);
            s += d;

            jj /= n;
        }
        printf("%2i)\t%i\t%i\n", n-s, dim[j], nonzero[j]);
        }
    }
    printf("#\n");

    free(dim);
    free(nonzero);
}

