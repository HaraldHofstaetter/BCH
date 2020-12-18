#include"bch.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

extern unsigned int VERBOSITY_LEVEL;


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

static void number_of_lyndon_words(uint8_t K, size_t N, size_t nLW[N]) {
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

static size_t word_index(size_t K, uint8_t w[], size_t l, size_t r) {
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

static size_t find_lyndon_word_index(uint32_t *WI, size_t l, size_t r, size_t wi) {
    /* finds index wi in the sorted list of indices WI. Start search at position l 
     * and stop it at position r. This function is only applied in situations where 
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



static int longest_right_lyndon_factor(uint8_t w[], size_t l, size_t r) {
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

static void genLW(size_t K, size_t n, size_t t, size_t p, uint8_t a[], uint8_t **W, 
        size_t wp[], uint32_t *WI, uint32_t *p1, uint32_t *p2, uint32_t *ii) {
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
                        p1[j] = find_lyndon_word_index(WI, ii[n1-1], wp[n1-1], wi1);
                        p2[j] = find_lyndon_word_index(WI, ii[n2-1], wp[n2-1], wi2);
                    }
                }
                j2 = j;
                wp[n0-1]++;
            } 
        }
    }
    else {
        a[t] = a[t-p];
        genLW(K, n, t+1, p, a, W, wp, WI, p1, p2, ii); 
        for (int j=a[t-p]+1; j<K; j++) {
             a[t] = j;
             genLW(K, n, t+1, t, a, W, wp, WI, p1, p2, ii);
        }
    }
}


void init_lyndon_words(lie_series_t *LS) {
    double t0 = tic();
    size_t nLW[LS->N];
    number_of_lyndon_words(LS->K, LS->N, nLW);
    size_t mem_len = 0;
    size_t N_LYNDON = 0;
    for (int n=1; n<=LS->N; n++) {
        N_LYNDON += nLW[n-1];
        mem_len += n*nLW[n-1];
    }
    LS->dim = N_LYNDON;
    LS->W = malloc(N_LYNDON*sizeof(uint8_t *));
    LS->p1 = malloc(N_LYNDON*sizeof(uint32_t)); 
    LS->p2 = malloc(N_LYNDON*sizeof(uint32_t)); 
    LS->nn = malloc(N_LYNDON*sizeof(uint8_t)); 
    LS->ii = malloc((LS->N+1)*sizeof(uint32_t)); 
    LS->W[0] = malloc(mem_len*sizeof(uint8_t)); 
    LS->ii[0] = 0;
    int m=0;
    for (int n=1; n<=LS->N; n++) {
        LS->ii[n] = LS->ii[n-1] + nLW[n-1];
        for (int k=0; k<nLW[n-1]; k++) {            
            if (m<N_LYNDON-1) { /* avoiding illegal W[N_LYNDON] */
                LS->W[m+1] = LS->W[m]+n;
            }
            LS->nn[m] = n;
            m++;
        }
    }
    assert(m==N_LYNDON);
    for (int i=0; i<LS->K; i++) {
        LS->p1[i] = i;
        LS->p2[i] = 0;
    }

    uint8_t a[LS->N+1];
    size_t wp[LS->N];
    for (int i=0; i<LS->N; i++) {
        wp[i] = LS->ii[i]; 
    }
    uint32_t *WI = malloc(N_LYNDON*sizeof(uint32_t));
    wp[0] = 1;
    LS->W[0][0] = 0;
    WI[0] = 0;

    for (int i=0; i<=LS->N; i++) {
        a[i] = 0; 
    }
    
    genLW(LS->K, LS->N, 1, 1, a, LS->W, wp, WI, LS->p1, LS->p2, LS->ii);

    free(WI);

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#number of Lyndon words of length<=%i over set of %i letters: %i\n", 
                LS->N, LS->K, LS->dim);
        printf("#init Lyndon words: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
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

static size_t multi_degree_index(size_t K, uint8_t w[], size_t l, size_t r) {
    size_t h[K];
    for (int j=0; j<K; j++) {
        h[j] = 0;
    }
    for (int j=l; j<=r; j++) {
        h[w[j]]++;
    }
    return tuple(K, h); 
}



uint32_t* multi_degree_indices(size_t K, size_t dim,  uint8_t **W, uint8_t *nn) {
    uint32_t *DI = malloc(dim*sizeof(uint32_t));
    for (int i=0; i<dim; i++) {
        DI[i] = multi_degree_index(K, W[i], 0, nn[i-1]);
    }
    return DI;
}

