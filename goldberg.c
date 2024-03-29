#include"bch.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

/* #define USE_PHI_FOR_GOLDBERG 1 */


/* n_partitions[n] = number of partitions of given number n, hardcoded for simplicity.
 * See: https://oeis.org/A000041
 */
static int n_partitions[33] = {1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 
           101, 135, 176, 231, 297, 385, 490, 627, 792, 1002, 1255, 1575, 1958, 
           2436, 3010, 3718, 4565, 5604, 6842, 8349};


static void partitions(int n, uint8_t **P) {
    /* INPUT: n
     * OUTPUT: P[j][i]=a[i] for the jth partition 
     *     a[0]+...+a[m-1]=n, a[0]>=a[1]>=...>=a[m-1]>=1
     *     of n, where  0<=j<n_partitions[n], P[j][m]=0.
     * METHOD: Algorithm P in Section 7.2.14 of
     *    D.E. Knuth: The Art of Computer Programming, Volume 4A
     */
    /* P1: Initialize */
    int np = 0;
    int a[n+1]; 
    a[0] = 0; 
    for (int m=2; m<=n; m++) {
        a[m] = 1;
    }
    int m = 1;
    while (1) { 
        /* P2: Store the final part */
        a[m] = n;
        int q = m - (n==1 ? 1 : 0);
        while (1) {
            /* # P3: Visit */
            for (int j=0; j<m; j++) {
                P[np][j] = a[j+1];
            }
            P[np][m] = 0;
            np++;        
            if (a[q]!=2) {
                break;
            }
            /* P4: Change 2 to 1+1 */
            a[q] = 1;
            q--;
            m++;
        }
        /* P5: Decrease a[q] */
        if (q==0) {
            break;
        }
        int x = a[q]-1;
        a[q] = x;
        n = m-q+1;
        m = q+1;
        /* P6: Copy x if necessary */
        while (n>x) {
            a[m] = x;
            m++;
            n -= x;
        }
    }  
}

/* n!*goldberg_denominator[n] = common denominator for all goldberg coefficients of degree <=n.
 * See http://oeis.org/A338025
 * Not static because also needed in phi.c 
 */
int goldberg_denominator[33] = {1, 1, 1, 2, 1, 6, 2, 6, 3, 10, 2, 6, 2, 210, 30, 12, 3, 30, 10, 
                               210, 42, 330, 30, 60, 30, 546, 42, 28, 2, 60, 4, 924, 231};

extern INTEGER FACTORIAL[]; /* defined in phi.c */


#ifndef USE_PHI_FOR_GOLDBERG    

static void compute_goldberg_coeffs(INTEGER y[], uint8_t q[], int Afirst, INTEGER d, INTEGER *C) {
    /* INPUT: q[0],...,q[m-1]>=1 ... exponents in word of the form w=A^q[0]B^q[1]... or w=B^q[0]A^q[1]... 
     *           depending on Afirst; q[m] = 0 is marker for first after last entry;
     *        d ... common denominator
     *        C ... auxilliary array of size at least N*N where N=q[0]+...+q[m-1]
     * OUTPUT: y[k]/d is the coefficient of tight subword w[end:end-k] of w in log(exp(A)exp(B))
     * METHOD: Algorithm 2 from the appendix of: 
     *   Harald Hofstätter, Smallest common denominators for the homogeneous components
     *   of the Baker-Campbell-Hausdorff series, https://arxiv.org/abs/2012.03818
     */
    int N = 0;
    int m = 0;
    for (; q[m]!=0; m++) {
        N += q[m];
    }
    for (int n=0; n<N*N; n++) {
        C[n] = 0;
    }
    int Acurrent = Afirst; 
    if (!(m%2)) {
        Acurrent = !Afirst;
    }
    int n = 0;
    for (int i=m-1; i>=0; i--) {
        for (int r=1; r<=q[i]; r++) {
            n += 1;
            INTEGER h = 0;
            if (i==m-1) {
                h = d/FACTORIAL[n];
            }
            else if (Acurrent && (i==m-2)) {
                h = d/(FACTORIAL[r]*FACTORIAL[q[i+1]]);
            }
            C[0 + (n-1)*N] = h;
            for (int k=2; k<n; k++) { /* case k>=2 */
                INTEGER h = 0;
                for (int j=1; j<=r; j++) {
                    if ((n>j) && C[k-2 + (n-j-1)*N]!=0) {
                        h += C[k-2 + (n-j-1)*N]/FACTORIAL[j];
                    }
                }
                if (Acurrent && (i<=m-2)) {
                    for (int j=1; j<=q[i+1]; j++) {
                        if ((n>r+j) && C[k-2 + (n-r-j-1)*N]!=0){
                            h += C[k-2 + (n-r-j-1)*N]/(FACTORIAL[r]*FACTORIAL[j]);
                        }
                    }
                }
                C[k-1 + (n-1)*N] = h;
            }
            C[n-1 + (n-1)*N] = d;
        }
        Acurrent = !Acurrent;
    }
    y[N] = 0;
    for (n=1;n<=N; n++) {
        INTEGER h = 0;
        for (int k=1; k<=N; k++) {
            h += C[k-1 + (n-1)*N]/(k%2 ? +k : -k);
        }
        y[N-n] = h;
    }
}
#endif


goldberg_t* goldberg(size_t n) {
    double t0 = tic();

    goldberg_t *G = malloc(sizeof(goldberg_t));
    int jj[n];
    for (int j=0; j<n; j++) {
        jj[j] = 0;
    }
    int ii[n+1];
    ii[0] = 0;
    int mem_len = 0;
    for (int j=1; j<=n; j++) {
        ii[j] = ii[j-1] + n_partitions[j];
        mem_len += n_partitions[j]*(j+1);
    }

    G->N = n;
    int np = ii[n];
    G->n_partitions = np;

    G->c = malloc(np*sizeof(INTEGER));     
    G->P = malloc(np*sizeof(uint8_t *));     
    G->P[0] = malloc(2*mem_len*sizeof(uint8_t));
    int i=0;
    for (int j=1; j<=n; j++) {
        for (int k=0; k<n_partitions[j]; k++) {
            if (i>0) {
                G->P[i] = G->P[i-1]+j+1;
            }
            i++;
        }
    }
    
    uint8_t **Pn = G->P+ii[n-1];
    partitions(n, Pn);

    G->denom = FACTORIAL[n]*goldberg_denominator[n];

#ifdef USE_PHI_FOR_GOLDBERG    
    expr_t *A = generator(0);
    expr_t *B = generator(1);
    expr_t *ex = logarithm(product(exponential(A), exponential(B)));

    INTEGER e[n+1];
    for (int j=0; j<n; j++){
        e[j] = 0;
    }
    e[n] = G->denom;

    uint8_t w[n];
#else
    INTEGER *C = malloc(n*n*sizeof(INTEGER));     
#endif
    INTEGER t[n+1];
    for (int k=0; k<n_partitions[n]; k++) {
#ifdef USE_PHI_FOR_GOLDBERG    
        int m=0;
        int l=0;
        for ( ; Pn[k][l]!=0; l++) {
            for (int i=0; i<Pn[k][l]; i++) {
                w[m] = l & 1;
                m++;
            }
        }

        phi(t, n+1, w, ex, e);
#else
        int l=0;
        for ( ; Pn[k][l]!=0; l++) { }
        compute_goldberg_coeffs(t, Pn[k], 1, G->denom, C);
#endif

        int j = n-1;
        int q = Pn[k][0];
        while (1) {
            G->c[ii[j]+jj[j]] = t[n-j-1];
            if (j<n-1) {
                G->P[ii[j]+jj[j]][0] = q;
                for (int i=1; i<l; i++) {
                    G->P[ii[j]+jj[j]][i] = Pn[k][i];
                }
                G->P[ii[j]+jj[j]][l] = 0;
            }
            jj[j]++;
            if (!((l==1 && q>1) || (l>1 && q>Pn[k][1]))) {
                break;
            }
            j--;
            q--;
        }
    }
#ifndef USE_PHI_FOR_GOLDBERG    
    free(C);
#endif

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#compute Goldberg coefficients: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }

    return G;
}


static int cumsum_partitions[33] = {1, 2, 4, 7, 12, 19, 30, 45, 67, 97, 139, 195, 272, 
    373, 508, 684, 915, 1212, 1597, 2087, 2714, 3506, 4508, 5763, 7338, 9296, 11732, 
    14742, 18460, 23025, 28629, 35471, 43820};

INTEGER goldberg_coefficient(int n, uint8_t w[], goldberg_t *G) {
    /* INPUT: word w of length n, 
     *        G ... struct as computed by function goldberg();
     * OUTPUT: returns coefficient of word w in log(exp(A)exp(B));
     * METHOD: binary search for entry in G corresponding to w.
     */
    int x[n];
    for (int i=0; i<n; i++) {
        x[i] = 0;
    }
    int k = 0;
    for (int i=1; i<n; i++) {
        if (w[i-1]==w[i]) {
            k++;
        }
        else {
            x[k]++;
            k = 0;
        }
    }
    x[k]++;
    
    uint8_t p[n+1];
    k = 0;
    for (int i=n-1; i>=0; i--) {
        for (int j=0; j<x[i]; j++) {
            p[k] = i+1;
            k++;
        }
    }
    for (; k<n+1; k++) {
        p[k] = 0;
    }

    int a = cumsum_partitions[n-1]-1;
    int b = a + n_partitions[n] - 1;
    while (a<=b) {
        int c = a + (b-a)/2;
        int i = 0;
        while (p[i]==G->P[c][i]) {
            if (p[i]==0) {
                return ((w[0]==1) && !(n&1)) ? -G->c[c] : G->c[c];
            }
            i+=1;
        }
        if (p[i]>G->P[c][i]) {
            b = c-1;
        }
        else {
            a = c+1;
        }
    }
    fprintf(stderr, "PANIC: partition not found"); 
    abort();
}


void print_goldberg(goldberg_t *G) {
    for (int i=0; i<G->n_partitions; i++) {
        for (int j=0; G->P[i][j]!=0; j++) {
            printf("%3i", G->P[i][j]);
        }
        printf("\t");
        print_RATIONAL(G->c[i], G->denom);
        printf("\n");
    }
}


void free_goldberg(goldberg_t *G) {
    free(G->P[0]);
    free(G->P);
    free(G->c);
    free(G);
}

