/* indent -npsl -npcs -br -i4 bch.c */

#include"bch.h"
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static size_t K;             /* number of generators */
static size_t N;             /* maximum length of Lyndon words (=maximum order of Lie series expansion) */
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

static size_t M = 0;         /* maximum lookup length */ 

typedef int32_t TINT_t ;
static TINT_t **T = NULL;    /* precomputed lookup table: word with index i has coefficient 
                                T[i][T_P[j]]  in basis element with number j.  */
static uint32_t *T_P = NULL;

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

static void gen_D(size_t K, size_t N, generator_t w[], size_t D[]) {
    size_t h[K];
    for (int r=N-1; r>=0; r--) {
        for (int j=0; j<K; j++) {
            h[j] = 0;
        }
        for (int l=r; l>=0; l--) {
            h[w[l]]++;
            D[l + r*N] = tuple(K, h);
        }
    }
}


static void gen_TWI(size_t K, size_t N, size_t M, generator_t w[], TINT_t **TWI) {
    for (int r=N-1; r>=0; r--) {
        int x = 0;
        int y = 1;
        if (K==2) {
            for (int l=r; l>=0 && l>r-(signed) M; l--) {
                x += w[l]*y;
                y <<= 1;
                TWI[l + r*N] = T[x + y - 2]; 
            }
        }
        else {
            int os = 0;
            for (int l=r; l>=0 && l>r-(signed) M; l--) {
                x += w[l]*y;
                y *= K;
                TWI[l + r*N] = T[x + os]; 
                os += y;
            }
        }
    }
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
        printf("#number of Lyndon words of length<=%li over set of %li letters: %li\n", N, K, N_LYNDON);
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
            printf("#compute righnormed basis elements: time=%g sec\n", t1);
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





static void gen_ith_word_of_length_n(size_t i, size_t n, generator_t w[]) {
    /* METHOD: compute base K expansion of i */
    for (int j=0; j<n; j++) {
        w[j] = 0;
    }
    size_t k=n-1;
    if (K==2) {
        while (i>0) {
            w[k] = i & 1;
            i >>= 1;
            k--;
        }
    }
    else {
        while (i>0) {
            w[k] = i%K;
            i/=K;
            k--;
        }
    }
}

static void init_lookup_table() {
/* Define T, T_P such that T[i]=T0[d]+T_P2[i]*T_D1[d] where d is multi degree index of 
 * the word with index i * and T_P=T_P1.
 * Here:
 * T_D1[i] = number of Lyndon words (Lyndon basis elements) which have multi degree index i. 
 * T D2[i] = number of all words of length <= M which have multi degree index i. 
 * T_P1[i] such that Lyndon word W[i] is the T_P1[i]-th Lyndon word having multi degree 
 *         index DI[i]. 
 * T_P2[i] such that word with index i is the T_P2[i]-th word in the list of all words 
 *          having the same multi degree index as the given word with index i 
 * Then: word with index i and multi-degree index d has coefficient
 *      T[i][T_P[j]] = T0[d][[T_P1[j] + T_P2[i]*T_D1[d]]
 * in basis element with number j. 
 * Basis element j and word i are assumed to have the same multi-degree  and length <= M. 
 * Note: transposing rows and columns such that the coefficient is given by 
 * T0[d][[T_P1[j]*T_D2[d] + T_P2[i]] results in a  significant loss of performance. 
 */
    if (M==0) {
        return;
    }    
    
    double t0 = tic();
    size_t H = DI[ii[M]-1]+1; 
    uint32_t *T_D1 = calloc(H, sizeof(uint32_t));
    uint32_t *T_P1 = calloc(ii[M], sizeof(uint32_t));
    for (int i=0; i<ii[M]; i++) {
        T_P1[i] = T_D1[DI[i]];
        T_D1[DI[i]]++;
    }
    uint32_t *T_D2 = calloc(H, sizeof(uint32_t));
    uint32_t *T_P2 = calloc((ipow(K, M+1)-1)/(K-1)-1, sizeof(uint32_t));
    uint32_t *FWD = calloc(H, sizeof(uint32_t)); 
    uint32_t *WDI = calloc((ipow(K, M+1)-1)/(K-1)-1, sizeof(uint32_t));
    generator_t w[M];
    for (int n=1; n<=M; n++) {
        int os = (ipow(K, n)-1)/(K-1)-1;  
        for (int i=0; i<ipow(K, n); i++) {
            gen_ith_word_of_length_n(i, n, w);
            size_t wi = word_index(K, w, 0, n-1);
            size_t di = multi_degree_index(K, w, 0, n-1);
            if (di<H) { /* this holds for all di except the last one */
               T_P2[wi] = T_D2[di];
               T_D2[di]++;
               WDI[wi] = di; 
               if (FWD[di]==0) {
                    FWD[di] = wi - os;
               }
            }
        }
    }
    TINT_t **T0 = calloc(H, sizeof(TINT_t*));
    for (int h=0; h<H; h++) {
        size_t d = T_D1[h]*T_D2[h];
        if (d>0) {
            T0[h] = calloc(d, sizeof(TINT_t)); // !!! TODO: eventually free this memory
        }
    }

    T = calloc((ipow(K, M+1)-1)/(K-1)-1, sizeof(int*));
    for (int wi=0; wi<(ipow(K, M+1)-1)/(K-1)-1; wi++) {
        int di = WDI[wi]; 
        T[wi] = T0[di] + T_P2[wi]*T_D1[di];
    }
    
    /* case n=1: */
    for (int j=0; j<K; j++) {
        uint32_t di =  DI[j];
        T0[di][T_P1[j] + T_P2[j]*T_D1[di]] = 1;
    }

    /* n>=2: */
    for (int n=2; n<=M; n++) {
        size_t ii1 = ii[n-1];
        size_t ii2 = ii[n]-1;
        int os = (ipow(K, n)-1)/(K-1)-1;
        #pragma omp parallel for schedule(dynamic,256)
        for (int j=ii1; j<=ii2; j++) {
            uint32_t j1 = p1[j];
            uint32_t j2 = p2[j];
            uint8_t n1 = nn[j1];
            uint8_t n2 = nn[j2];
            int Kn1 = ipow(K, n1);
            int Kn2 = ipow(K, n2);
            int os1 = (Kn1-1)/(K-1)-1;
            int os2 = (Kn2-1)/(K-1)-1;
            uint32_t di = DI[j];
            uint32_t di1 = DI[j1];
            uint32_t di2 = DI[j2];
            int y1 = T_D2[di1];
            int y2 = T_D2[di2];
            int x = T_D1[di];
            int x1 = T_D1[di1];
            int x2 = T_D1[di2];
            TINT_t *L = T0[di]+T_P1[j];
            TINT_t *L1 = T0[di1]+T_P1[j1];
            TINT_t *L2 = T0[di2]+T_P1[j2];
            for (int i1=FWD[di1]; i1<=WI[j1]; i1++) {
                if (WDI[i1+os1]==di1) {
                    int i = T_P2[i1*Kn2+FWD[di2]+os];
                    int c1 = L1[T_P2[i1+os1]*x1];
                    if (c1!=0) {
                        int k = i*x;
                        int k2 = 0;
                        for (int i2=0; i2<y2; i2++) {
                            int c2 = L2[k2];
                            L[k] = c1*c2;  
                            k += x;
                            k2 += x2;
                            i++;
                        }
                    }
                }
            }
            for (int i2=FWD[di2]; i2<=WI[j2]; i2++) {
                if (WDI[i2+os2]==di2) {
                    int i = T_P2[i2*Kn1+FWD[di1]+os];
                    int c2 = L2[T_P2[i2+os2]*x2];
                    if (c2!=0) {
                        int k = i*x;
                        int k1 = 0;
                        for (int i1=0; i1<y1; i1++) {
                            int c1 = L1[k1];
                            L[k] -= c1*c2;  
                            k += x;
                            k1 += x1;
                            i++;
                        }
                    }
                }
            }
        }
    }
    T_P = T_P1;
    free(T_P2);
    free(T_D1);
    free(T_D2);
    free(FWD);
    free(WDI);
    free(T0);

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#lookup table for word lengths<=%li\n", M);
        printf("#init lookup table: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}

static void free_lookup_table(void) {
    free(T);
    free(T_P);
    // TODO: free memory as indicated above !!!
}


static int coeff_word_in_basis_element(/* generator_t w[], */ size_t l, size_t r, size_t j, size_t N1, size_t D[], TINT_t **TWI) {  
    /* computes the coefficient of the word with index wi=W2I[l+r*N] in the basis element
     * with number j.
     * W2I is a table of indices such that W2I[l'+r'*N] is the index of the subword w[l':r'] 
     * of a word w which is given only implicitely.
     * D is a table of multi degree indices such that D[l'+r'*N] is the multi degree index
     * of w[l':r']. 
     */
    int n=r-l+1;

    if (n==1) {
        return DI[j]==D[l + r*N1];
    }

    if (n<=M) {  /* use lookup table */
        return TWI[l + r*N1][T_P[j]]; 
    }


    size_t j1 = p1[j];
    size_t j2 = p2[j];

    size_t m1 = nn[j1];
    size_t m2 = r-l+1-m1;

    int mi = DI[j1];
    int c2 = 0;
    if (D[l+m2 + r*N1] == mi) {
        c2 = coeff_word_in_basis_element(/* w, */ l+m2, r, j1, N1, D, TWI); 
        if (c2!=0) {
            c2 *= coeff_word_in_basis_element(/* w, */ l, l+m2-1, j2, N1, D, TWI); 
        }
    }

    int c1 = 0;
    if (D[l + (l+m1-1)*N1] == mi) {
        c1 = coeff_word_in_basis_element(/* w, */ l+m1, r, j2, N1, D, TWI); 
        if (c1!=0) {
            c1 *= coeff_word_in_basis_element(/* w, */ l, l+m1-1, j1, N1, D, TWI); 
        }
    }

    return c1 - c2;
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

static void compute_word_coefficients(int N, expr_t* ex, INTEGER c[], INTEGER denom, int classical_bch) {
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

    if (classical_bch) {
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
    }
    else {
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
    }
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#compute coeffs of words: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}    


static void convert_to_lie_series(int N, INTEGER c[]) {
    if (VERBOSITY_LEVEL>=2) {
#ifdef _OPENMP
        printf("# degree     #basis        time thread\n");
#else
        printf("# degree     #basis        time\n");
#endif
    }
    double t0 = tic();

    size_t i1 = ii[N-1];
    size_t i2 = ii[N]-1;

    size_t h1 = DI[i1];
    size_t h2 = DI[i2];

    double h_time[h2-h1+1];
    int h_n[h2-h1+1];
#ifdef _OPENMP
    int h_thread[h2-h1+1];
#endif
    int hh[h2-h1+1];
    if (K==2 && (N&1)==0) {
        /* generate loop order corresponding to decreasing running times:
         * N/2, N/2-1, N/2+1, N/2-1, N/2+1, ..., 1, N-1
         */ 
        int n1 = N>>1;
        hh[0]=n1-1;
        for (int k=1; k<n1; k++) {
            hh[N-2*k-1] = k-1;
            hh[2*k+1-1] = n1+k-1;
        }
    }
    else {
        for (int k=0; k<=h2-h1; k++) {
            hh[k] = k;
        }
    }

    #pragma omp parallel 
    {
    int *jj = calloc(N_LYNDON, sizeof(int));  // N_LYNDON far too large upper bound
    size_t D[N*N];
    TINT_t *TWI[N*N];
    
    size_t JW[N];
    size_t JB[N];

    /* Note: We choose schedule(dynamic, 1) because each
     * iteration of the loop is associated with a specific 
     * multi degree index, and the work to be done varies widely
     * for different multi degree indices. 
     */
    #pragma omp for schedule(dynamic,1) 
    for (int k=0; k<=h2-h1; k++) {
        int h = h1+hh[k];
        h_time[k] = tic();
        h_n[k] = 0;
#ifdef _OPENMP
        h_thread[k] = omp_get_thread_num();
#endif
        size_t kW1 = N+1; 

        int jj_max = 0; 
        for (int i=i1; i<=i2; i++) {
            if (DI[i]==h) {
                jj[jj_max] = i;
                jj_max++;
            }
        }

        for (int x=0; x<jj_max; x++) {
            int i = jj[x];
            {
                h_n[k]++;

                size_t kW = get_right_factors(i, JW, N);
                generator_t *w = W[i];

                kW1 = kW<kW1 ? kW : kW1;
                size_t N1 = N-kW1;

                gen_D(K, N1, w+kW1, D);
                gen_TWI(K, N1, M, w+kW1, TWI);

                for (int y=0; y<=x-1; y++) {
                    int j = jj[y];
                    size_t kB = get_right_factors(j, JB, N);
                    if (D[kB-kW1 + (N1-1)*N1] == DI[JB[kB]]) { /* check if multi degrees match */
                        int d = coeff_word_in_basis_element(/* w+kW1, */ kB-kW1, N1-1, JB[kB], N1, D, TWI);
                        if (d!=0) {
                            for (int k=0; k<=kB && k<=kW; k++) {
                                c[JW[k]] -= d*c[JB[k]];
                            }
                        }
                    }
                }
            }
        }
        h_time[k] = toc(h_time[k]); 
        if (VERBOSITY_LEVEL>=2) {
#ifdef _OPENMP
            printf("#%7i %10i %11.2f   %4i\n", hh[k]+1, h_n[k], h_time[k], h_thread[k]);
#else
            printf("#%7i %10i %11.2f\n", hh[k]+1, h_n[k], h_time[k]);
#endif
            fflush(stdout);
        }
    }
    free(jj);
    }

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#convert to lie series: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}

/* tables beta_num_h, beta_num_l, beta_den_h, beta_den_l defined such that the 
   rational numbers
      beta[k] = (H*beta_num_h[k]+beta_num_l[k])/(H*beta_den_h[k]+beta_den_l[k])
   (where H = 1000000000000000000) are the coefficients of the power series 
   of the function f(x)=tanh(x/2), i.e., 
   beta[] = {           1/2,                                   
                       -1/24,                                  
                        1/240,                                 
                      -17/40320,                               
                       31/725760,                              
                     -691/159667200,                           
                     5461/12454041600,                         
                  -929569/20922789888000,                      
                  3202291/711374856192000,                     
               -221930581/486580401635328000,                  
               4722116521/102181884343418880000,               
             -56963745931/12165654935945871360000,             
           14717667114151/31022420086661971968000000,          
        -2093660879252671/43555477801673408643072000000,       
        86125672563201181/17683523987479403909087232000000,    
   -129848163681107301953/263130836933693530167218012160000000}

   This data was computed with the following Julia code:
   
   n=16
   A = zeros(Rational{BigInt}, 2*n+1)
   B = similar(A)
   for k = 0:2*n 
       A[k+1] = 1//(k+1)
       for j = k:-1:1
           A[j] = j*(A[j] - A[j+1])
       end
       B[k+1] =  A[1]
   end
   beta = [2*(2^(2*k)-1)*B[2*k+1]/factorial(BigInt(2*k)) for k=1:n]
*/   

static long long int beta_num_h[16] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -129};
    
static long long int beta_num_l[16] = {
    1, -1, 1, -17, 31, -691, 5461, -929569, 3202291, -221930581, 4722116521, -56963745931, 14717667114151,
    -2093660879252671, 86125672563201181, -848163681107301953};

static long long int beta_den_h[16] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 102, 12165, 31022420, 43555477801, 17683523987479, 263130836933693530};

static long long int beta_den_l[16] = {
    2, 24, 240, 40320, 725760, 159667200, 12454041600, 20922789888000, 711374856192000, 
    486580401635328000, 181884343418880000, 654935945871360000, 86661971968000000, 673408643072000000, 
    403909087232000000, 167218012160000000};


static void  compute_BCH_terms_of_order_N(INTEGER c[], INTEGER denom) {
    double t0 = tic();
    assert(!(N&1));
    int N2 = N/2;
    INTEGER beta_num[N2];
    INTEGER beta_den[N2];
    INTEGER H = 1000000000000000000;
    for (int k=0; k<N2; k++) {
        beta_num[k] = beta_num_h[k]*H + beta_num_l[k];
        beta_den[k] = beta_den_h[k]*H + beta_den_l[k];
    }

    #pragma omp parallel for schedule(dynamic,256)
    for (int i=ii[N-1]; i<=ii[N]-1; i++) {
        c[i] = 0;
        int k = 0;
        int l = 0;
        int q = i;
        while (p1[q]==0) {
            k += 1;
            q = p2[q];
            if (k&1) {
                INTEGER d = c[q]/beta_den[l];
                if (d*beta_den[l]!=c[q]) {
                    fprintf(stderr, "ERROR: divisibility check failed in compute_BCH_terms_of_order_N");
                    exit(EXIT_FAILURE);
                }
                c[i] += beta_num[l]*d; 
                l += 1;
            }
        }
    } 

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#compute terms of order %li: time=%g sec\n", N, t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


static void init_all(size_t number_of_generators, size_t order, 
                     size_t max_lookup_length, int rightnormed) {
    K = number_of_generators;
    N = order;
    init_factorial(N);
    init_lyndon_words(rightnormed);
    //if (rightnormed) {
    //    M = 0;
    //}
    //else {
        M = max_lookup_length;
        init_lookup_table();
    //}
}

static void free_all(void) {
    free_factorial();
    free_lookup_table();
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

lie_series_t lie_series(size_t K, expr_t* expr, size_t N, int64_t fac, size_t M, int rightnormed) {
    double t0 = tic();
    init_all(K, N, M, rightnormed);
    INTEGER *c = malloc(N_LYNDON*sizeof(INTEGER));
    INTEGER denom = common_denominator(N)*fac;
    compute_word_coefficients(N, expr, c, denom, 0);
    convert_to_lie_series(N, c);
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

lie_series_t BCH(size_t N, size_t M, int rightnormed) {
    double t0 = tic();
    expr_t *A = generator(0);
    expr_t *B = generator(1);
    expr_t *expr = logarithm(product(exponential(A), exponential(B)));
    init_all(2, N, M, rightnormed);
    INTEGER *c = malloc(N_LYNDON*sizeof(INTEGER));
    INTEGER denom = common_denominator(N);
    if (N%2) {
        compute_word_coefficients(N, expr, c, denom, 1);
        convert_to_lie_series(N, c);
    }
    else {
        compute_word_coefficients(N-1, expr, c, denom, 1);
        convert_to_lie_series(N-1, c);
        compute_BCH_terms_of_order_N(c, denom);
    }
    lie_series_t LS = gen_result(c, denom);
    free_all();
    free_expr(A);
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

lie_series_t symBCH(size_t N, size_t M, int rightnormed) {
    double t0 = tic();
    expr_t *halfA = generator(0);
    expr_t *B = generator(1);
    expr_t *expr = logarithm(product(product(exponential(halfA), exponential(B)), 
                                     exponential(halfA)));
    init_all(2, N, M, rightnormed);
    INTEGER *c = malloc(N_LYNDON*sizeof(INTEGER));
    INTEGER denom = common_denominator(N);
    if (VERBOSITY_LEVEL>=1) {
        printf("#NOTE: in the following expression, A stands for A/2\n");
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
    compute_word_coefficients(N, expr, c, denom, 0);
    convert_to_lie_series(N, c);
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


void print_word(lie_series_t *LS,  size_t i) {
    if (i<LS->K) {
        printf("%c", (char) ('A'+i));
    }
    else {
        print_word(LS, LS->p1[i]);
        print_word(LS, LS->p2[i]);
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
        if (what & PRINT_WORD) printf("\tword");
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
        if (what & PRINT_WORD) {
            printf("\t");
            print_word(LS, i);
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



