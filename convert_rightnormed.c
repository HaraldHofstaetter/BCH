#include"bch.h"
#include <stdio.h>
#include <stdlib.h>

extern unsigned int VERBOSITY_LEVEL;


static int coeff_word_in_rightnormed(uint8_t w[], uint8_t c[], int l1, int r1, int l2) {
    if (l1==r1) {
        return w[l1]==c[l2] ? 1 : 0;
    }
    else {
        return (w[l1]==c[l2] ? coeff_word_in_rightnormed(w, c, l1+1, r1, l2+1) : 0) -
               (w[r1]==c[l2] ? coeff_word_in_rightnormed(w, c, l1, r1-1, l2+1) : 0);
    }
}


static void integer_lu(int n, int64_t *A) {    
    /* LU factorization */
    for (int k=0; k<n; k++) {
        int64_t s = A[k+n*k];
        if ((s!=1) && (s!=-1)) {
            fprintf(stderr, "ERROR: integer LU factorization does not exist"); 
            exit(EXIT_FAILURE);
        }
        #pragma omp parallel for schedule(static, 32) 
        for (int i=k+1; i<n; i++) {
            if (A[i+n*k]!=0) {
                A[i+n*k] *= s;
                for (int j=k+1; j<n; j++) { 
                    A[i+n*j] -= A[i+n*k]*A[k+n*j];
                }
            }
        }
    }
}


static void integer_lu_solve(int n, int64_t *A, INTEGER *x) {    
    /* forward substitution */
    for (int i=1; i<n; i++) {
        INTEGER s=0;
        for (int j=0; j<i; j++) {
            s += A[i+n*j]*x[j];
        }
        x[i] -= s;
    }

    /* back substitution */
    x[n-1] *= A[n-1 +n*(n-1)];
    for (int i=n-2; i>=0; i--) {
        INTEGER s=0;
        for (int j=i+1; j<n; j++) {
            s += A[i+n*j] * x[j];
        }
        x[i] = (x[i] - s)*A[i+n*i];
    }
}


void convert_to_rightnormed_lie_series(lie_series_t *LS, int N, int odd_orders_only) {
    double t0 = tic();
    for (int n=2; n<=N; n++) { /* over all word sizes */
    if ((!odd_orders_only)||(n&1)) {
        size_t i1 = LS->ii[n-1];
        size_t i2 = LS->ii[n]-1;
        uint32_t *DI  = multi_degree_indices( LS->K, LS->dim, LS->W, LS->nn);
        size_t h1 = DI[i1];
        size_t h2 = DI[i2];
        for (int h=h1; h<=h2; h++) { /* over all multi-degrees */
            /* get dimension */
            int m=0;
            for (int j=i1; j<=i2; j++) {
                if (DI[j]==h) {
                    m++;
                }
            }
            if (m==0) {
                continue;
            }

            /* set up matrix and right-hand side */
            INTEGER *x = calloc(m, sizeof(INTEGER));
            int64_t *A = calloc(m*m, sizeof(int64_t));
            int jj=0;
            for (int j=i1; j<=i2; j++) {
                if (DI[j]==h) {
                    int ii=0;
                    for (int i=i1; i<=i2; i++) {
                        if (DI[i]==h) {
                            A[ii+jj*m] = coeff_word_in_rightnormed(LS->W[i], LS->R[j], 0, n-1, 0);
                            ii++;
                        }
                    }
                    x[jj] = LS->c[j];
                    jj++;
                }
            }

            integer_lu(m, A);
            integer_lu_solve(m, A, x);

            /* copy result */
            jj=0;
            for (int j=i1; j<=i2; j++) {
                if (DI[j]==h) {
                    LS->c[j] = x[jj];
                    jj++;
                }
            }

            free(x);
            free(A);
        }
        free(DI);
    }
    else {
        for (int j=LS->ii[n-1]; j<=LS->ii[n]-1; j++) {
            LS->c[j] = 0;
        }
    }
    }
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#convert to rightnormed lie series: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


extern INTEGER beta_num[];  /* defined in convert_lyndon.c */
extern INTEGER beta_den[];  /* defined in convert_lyndon.c */


void compute_rightnormed_BCH_terms_of_even_orders(lie_series_t *LS) {
    double t0 = tic();

    for(int n=2; n<=LS->N; n+=2) {
        #pragma omp parallel for schedule(dynamic,256)
        for (int i=LS->ii[n-1]; i<=LS->ii[n]-1; i++) {
            LS->c[i] = 0;
            int k=0;
            int l=0;
            while (LS->R[i][k]==LS->K-1) {
                k += 1;
                if (k&1) {
                    int q = LS->ii[n-1-k];
                    for (; q<=LS->ii[n-k]-1; q++) {
                        int m=0;
                        for (; (m<n-k) && (LS->R[q][m]==LS->R[i][k+m]) ; m++) {}
                        if (m==n-k) {
                            break;
                        }
                    }
                    if (q>LS->ii[n-k]) {
                        fprintf(stderr, "ERROR: basis element not found in compute_rightnormed_BCH_terms_of_even_orders");
                        exit(EXIT_FAILURE);
                    }
                    INTEGER d = LS->c[q]/beta_den[l];
                    if (d*beta_den[l]!=LS->c[q]) {
                        fprintf(stderr, "ERROR: divisibility check failed in compute_rightnormed_BCH_terms_of_even_orders");
                        exit(EXIT_FAILURE);
                    }
                    LS->c[i] -= beta_num[l]*d; 
                    l += 1;
                }
            }
        }
    }

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#compute terms of even orders: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


void init_rightnormed(lie_series_t *LS) {
    double t0 = tic();
    size_t mem_len = 0;
    for (int n=1; n<=LS->N; n++) {
        mem_len += n*(LS->ii[n]-LS->ii[n-1]);
    }

    LS->R = malloc(LS->dim*sizeof(uint8_t *));
    LS->R[0] = malloc(mem_len*sizeof(uint8_t)); 
    for (int i=1; i<LS->dim; i++) {
        LS->R[i] = LS->R[i-1] + LS->nn[i-1];
    }

    #pragma omp for schedule(dynamic,256) 
    for (int i=0; i<LS->dim; i++) {
        lyndon2rightnormed(LS->nn[i], LS->W[i], LS->R[i]);
    }

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#init rightnormed basis elements: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


