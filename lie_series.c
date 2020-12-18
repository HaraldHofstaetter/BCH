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

static int ipow(int base, unsigned int exp) {
    /* computes base^exp 
     * METHOD: see https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
     */
    if (base==2) {
        return 2<<(exp-1);
    }
    else {
        int result = 1;
        for (;;)
        {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            if (!exp)
                break;
            base *= base;
        }
        return result;
    }
}


static void moebius_mu(size_t N, int mu[N]) {
    /* INPUT: N
     * OUTPUT: mu[n] = Moebius mu function of n+1, n=0,...,N-1
     * METHOD: see https://mathoverflow.net/questions/99473/calculating-m%C3%B6bius-function
     */
    for (int i=0; i<N; i++) {
        mu[i] = 0;
    }
    mu[0] = 1;
    for (int n=1; n<=N/2; n++) {
        int mu_n = mu[n-1];
        for(int i=2*n-1; i<N; i+=n) {
            mu[i] -= mu_n;
        }
    }
}

static void number_of_lyndon_words(generator_t K, size_t N, size_t nLW[N]) {
    /* INPUT: K ... number of letters
     *        N ... maximum lenght of lyndon words
     * OUTPUT: nLW[n] ... number of lyndon words with K letters of length n+1, n=0,...,N-1
     * METHOD: Witt's formula
     */
    int mu[N];
    moebius_mu(N, mu);

    for (int n=1; n<=N; n++) {
        int d = 1;
        int h = 0;
        while (d*d < n) {
            div_t d1r = div(n, d);
            if (d1r.rem==0) {
               int d1 = d1r.quot; 
               h += mu[d-1]*ipow(K, d1)+mu[d1-1]*ipow(K, d);
            }
            d++;
        }
        if (d*d == n) {
            h += mu[d-1]*ipow(K, d);
        }
        nLW[n-1] = h/n;
    }
}

static size_t word_index(size_t K, generator_t w[], size_t l, size_t r) {
    /* computes the index of the subword w[l:r] of w starting at position l and
     * ending at position r. The index is given as w[l:r] interpreted as a K-adic
     * number plus the number (K^n-1)/(K-1)-1 of words of length < n, where 
     * n = r-l+1 = length of w[l:r] 
     */
    size_t x = 0;
    size_t y = 1;
    if (K==2) {
        for (int j=r; j>= (signed) l; j--) { /* CAUTION! comparison between signed and unsigned */
            x += w[j]*y;
            y <<= 1;
        }
        return x + y - 2;
    }
    else {
        for (int j=r; j>= (signed) l; j--) { /* CAUTION! comparison between signed and unsigned */
            x += w[j]*y;
            y *= K;
        }
        return x + (y-1)/(K-1) - 1;
    }
}

static size_t find_lyndon_word_index(size_t l, size_t r, size_t wi) {
    /* finds index wi in the sorted list of indices WI. Start search at position l 
     * and stop * it at position r. This function is only applied in situations where 
     * the search will not fail.
     * METHOD: binary search
     */
    while (l<=r) {
        size_t m = l + (r-l)/2;
        if (WI[m]==wi) {
            return m;
        }
        if (WI[m]<wi) {
            l = m+1;
        }
        else {
            r = m-1;
        }
    }
    fprintf(stderr, "ERROR: Lyndon word index not found: %li\n", wi);
    exit(EXIT_FAILURE);
}

static unsigned int binomial(unsigned int n, unsigned int k) {
    /* computes binomial coefficient n over k
     * METHOD: from Julia base library, see
     * https://github.com/JuliaLang/julia/blob/master/base/intfuncs.jl     
     */ 
    if (k < 0 || k > n ) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    if (k == 1) {
        return n;
    }
    if (k > (n>>1)) {
        k = (n - k);
    }
    uint64_t x = n - k +1;
    uint64_t nn = x;
    nn++;
    uint64_t rr = 2;
    while (rr <= k) {
        x = (x*nn) / rr;  
        rr++;
        nn++;
    }
    return x;
}

static size_t tuple(size_t K, size_t h[]) {
    if (K==2) {
        int s = h[0]+h[1];
        return ((s*(s+1))>>1)+h[1];
    }
    else {
        size_t index = 0;
        size_t n = 0;
        for (int k=0; k<K; k++) {
            n += h[K-k-1];
            index += binomial(k+n, n-1);
        }
        return index;
    }
}

static size_t multi_degree_index(size_t K, generator_t w[], size_t l, size_t r) {
    size_t h[K];
    for (int j=0; j<K; j++) {
        h[j] = 0;
    }
    for (int j=l; j<=r; j++) {
        h[w[j]]++;
    }
    return tuple(K, h); 
}


static int longest_right_lyndon_factor(generator_t w[], size_t l, size_t r) {
/* returns starting position of the longest right Lyndon factor of the subword w[l:r]
 * METHOD: based on the algorithm MaxLyn from
 *   F. Franek, A. S. M. S. Islam, M. S. Rahman, W. F. Smyth: Algorithms to Compute the Lyndon Array. 
 *   Stringology 2016: 172-184
 */
    for (int j=l+1; j<r; j++) {        
        int i = j+1;
        while (i <= r) {
            int k = 0;
            while ((i+k <= r) && (w[j+k]==w[i+k])) {
                k += 1;
            } 
            if ((i+k > r) || (w[j+k] >= w[i+k])) {
                break;
            }
            else {
                i += k + 1;
            }
        }
        if (i==r+1) {
            return j;
        }
    }
    return r;
}

/* The following two functions are for the generation of Lyndon words.
 * METHOD: Algorithm 2.1 from
 *   K. Cattell, F. Ruskey, J. Sawada, M. Serra, C.R. Miers, Fast algorithms 
 *   to generate necklaces, unlabeled necklaces and irreducible polynomials over GF(2), 
 *   J. Algorithms 37 (2) (2000) 267–282
 */

static void genLW(size_t K, size_t n, size_t t, size_t p, generator_t a[], size_t wp[]) {
    if (t>n) {
        if (p==n) {
            int H = 0;
            size_t j2 = 0;
            while ((longest_right_lyndon_factor(a, H+1, n)==H+2) && (a[H+1]==0)) { 
                H++;
            }
            for (int h=H; h>=0; h--)  {
                size_t n0 = n-h;
                size_t j = wp[n0-1];
                for (int i=0; i<n0; i++) {
                    W[j][i] = a[i+h+1];
                }
                WI[j] = word_index(K, a, h+1, n);
                DI[j] = multi_degree_index(K, a, h+1, n);
                if (n0>1) {
                    if (h<H) {
                        p1[j] = 0;
                        p2[j] = j2;
                    }
                    else {
                        size_t m = longest_right_lyndon_factor(a, h+1, n);
                        size_t wi1 = word_index(K, a, h+1, m-1);
                        size_t wi2 = word_index(K, a, m, n);
                        int n1 = m-h-1;
                        int n2 = n0-n1;
                        p1[j] = find_lyndon_word_index(ii[n1-1], wp[n1-1], wi1);
                        p2[j] = find_lyndon_word_index(ii[n2-1], wp[n2-1], wi2);
                    }
                }
                j2 = j;
                wp[n0-1]++;
            } 
        }
    }
    else {
        a[t] = a[t-p];
        genLW(K, n, t+1, p, a, wp); 
        for (int j=a[t-p]+1; j<K; j++) {
             a[t] = j;
             genLW(K, n, t+1, t, a, wp);
        }
    }
}


static void init_lyndon_words(int rightnormed) {
    double t0 = tic();
    size_t nLW[N];
    number_of_lyndon_words(K, N, nLW);
    size_t mem_len = 0;
    N_LYNDON = 0;
    for (int n=1; n<=N; n++) {
        N_LYNDON += nLW[n-1];
        mem_len += n*nLW[n-1];
    }
    W = malloc(N_LYNDON*sizeof(generator_t *));
    R = NULL;
    if (rightnormed) {
        R = malloc(N_LYNDON*sizeof(generator_t *));
    }
    p1 = malloc(N_LYNDON*sizeof(uint32_t)); 
    p2 = malloc(N_LYNDON*sizeof(uint32_t)); 
    nn = malloc(N_LYNDON*sizeof(uint8_t)); 
    WI = malloc(N_LYNDON*sizeof(uint32_t));
    DI = malloc(N_LYNDON*sizeof(uint32_t));
    ii = malloc((N+1)*sizeof(uint32_t)); 
    W[0] = malloc(mem_len*sizeof(generator_t)); 
    if (rightnormed) {
        R[0] = malloc(mem_len*sizeof(generator_t)); 
    }
    ii[0] = 0;
    int m=0;
    for (int n=1; n<=N; n++) {
        ii[n] = ii[n-1] + nLW[n-1];
        for (int k=0; k<nLW[n-1]; k++) {            
            if (m<N_LYNDON-1) { /* avoiding illegal W[N_LYNDON] */
                W[m+1] = W[m]+n;
                if (rightnormed) {
                    R[m+1] = R[m]+n;
                }
            }
            nn[m] = n;
            m++;
        }
    }
    assert(m==N_LYNDON);
    for (int i=0; i<K; i++) {
        p1[i] = i;
        p2[i] = 0;
    }

    generator_t a[N+1];
    size_t wp[N];
    for (int i=0; i<N; i++) {
        wp[i] = ii[i]; 
    }
    wp[0] = 1;
    W[0][0] = 0;
    WI[0] = 0;
    DI[0] = 1;

    for (int i=0; i<=N; i++) {
        a[i] = 0; 
    }
    
    genLW(K, N, 1, 1, a, wp);

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#number of Lyndon words of length<=%i over set of %i letters: %li\n", N, K, N_LYNDON);
        printf("#init Lyndon words: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }

    if (rightnormed) {
        double t0 = tic();

        #pragma omp for schedule(dynamic,256) 
        for (int i=0; i<N_LYNDON; i++) {
            lyndon2rightnormed(nn[i], W[i], R[i]);
        }
        
        if (VERBOSITY_LEVEL>=1) {
            double t1 = toc(t0);
            printf("#compute rightnormed basis elements: time=%g sec\n", t1);
            if (VERBOSITY_LEVEL>=2) {
                fflush(stdout);
            }
        }
    }
}

static void free_lyndon_words(void) {
    free(W[0]);
    free(W);
    free(ii);
    free(WI);
    free(DI);
    /* Note: p1, p2, and nn are taken over by a lie_series_t struct
       and are eventually freed by free_lie_series */
}


static inline size_t get_right_factors(size_t i, size_t J[], size_t kmax) {
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



