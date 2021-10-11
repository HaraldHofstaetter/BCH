#include"bch.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


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
            fprintf(stderr, "PANIC: integer LU factorization does not exist\n"); 
            abort();
        }
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


static void integer_lu_solve_f(int n, int64_t *A, FLOAT *x) {    
    /* forward substitution */
    for (int i=1; i<n; i++) {
        FLOAT s=0;
        for (int j=0; j<i; j++) {
            s = add_f(s, mul_f(i64_to_f(A[i+n*j]), x[j]));
        }
        x[i] = sub_f(x[i], s);
    }

    /* back substitution */
    x[n-1] = mul_f(x[n-1], i64_to_f(A[n-1+n*(n-1)]));
    for (int i=n-2; i>=0; i--) {
        FLOAT s=0;
        for (int j=i+1; j<n; j++) {
            s = add_f(s, mul_f(i64_to_f(A[i+n*j]), x[j]));
        }
        x[i] = mul_f(sub_f(x[i], s), i64_to_f(A[i+n*i]));
    }
}



void convert_to_rightnormed_lie_series(lie_series_t *LS, int N, int odd_degrees_only) {
    double t0 = tic();
    uint32_t *DI  = multi_degree_indices( LS->K, LS->dim, LS->W, LS->nn);
    for (int n=2; n<=N; n++) { /* over all word sizes */
    if ((!odd_degrees_only)||(n&1)) {
        size_t i1 = LS->ii[n-1];
        size_t i2 = LS->ii[n]-1;
        size_t h1 = DI[i1];
        size_t h2 = DI[i2];
        #pragma omp parallel for schedule(dynamic,1)
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

            /* set up matrix */
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
                    jj++;
                }
            }
            integer_lu(m, A);

            if (LS->c) {
                /* set up right-hand side */
                INTEGER *x = calloc(m, sizeof(INTEGER));
                jj=0;
                for (int j=i1; j<=i2; j++) {
                    if (DI[j]==h) {
                        x[jj] = LS->c[j];
                        jj++;
                    }
                }
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
            }
            else {
                /* set up right-hand side */
                FLOAT *x = calloc(m, sizeof(FLOAT));
                jj=0;
                for (int j=i1; j<=i2; j++) {
                    if (DI[j]==h) {
                        x[jj] = LS->c_f[j];
                        jj++;
                    }
                }
                integer_lu_solve_f(m, A, x);

                /* copy result */
                jj=0;
                for (int j=i1; j<=i2; j++) {
                    if (DI[j]==h) {
                        LS->c_f[j] = x[jj];
                        jj++;
                    }
                }
                free(x);
            }
            free(A);
        }
    }
    else {
        if (LS->c) {
            for (int j=LS->ii[n-1]; j<=LS->ii[n]-1; j++) {
                LS->c[j] = 0;
            }
        }
        else {
            for (int j=LS->ii[n-1]; j<=LS->ii[n]-1; j++) {
                LS->c_f[j] = 0;
            }
        }
    }
    }
    free(DI);

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#convert to rightnormed Lie series: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}


extern INTEGER beta_num[];  /* defined in convert_lyndon.c */
extern INTEGER beta_den[];  /* defined in convert_lyndon.c */


void compute_rightnormed_BCH_terms_of_even_degrees(lie_series_t *LS) {
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
                        fprintf(stderr, "PANIC: basis element not found in compute_rightnormed_BCH_terms_of_even_degrees");
                        abort();
                    }
                    INTEGER d = LS->c[q]/beta_den[l];
                    if (d*beta_den[l]!=LS->c[q]) {
                        fprintf(stderr, "PANIC: divisibility check failed in compute_rightnormed_BCH_terms_of_even_degrees");
                        abort();
                    }
                    LS->c[i] -= beta_num[l]*d; 
                    l += 1;
                }
            }
        }
    }

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#compute terms of even degrees: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}


static void adjust_p1_p2(lie_series_t *LS) {
    
    for (int k=0; k<LS->K; k++) {
        LS->p1[k] = k;
        LS->p2[k] = 0;
    }

    for (int n=2; n<=LS->N; n++) {
        for (int i=LS->ii[n-1]; i<LS->ii[n]; i++) {
            int l = LS->R[i][0];
            LS->p1[i] = l;
            int c=LS->ii[n-2];
            for (; c<LS->ii[n-1]; c++) {
                /* Test if R[i][2:end]==R[c] */
                int j=0;
                for ( ; j<LS->nn[i]-1; j++) {
                    if (LS->R[i][j+1]!=LS->R[c][j]) {
                        break;
                    }
                }
                if (j==LS->nn[i]-1) {
                    LS->p2[i] = c;
                    break;
                }
            }
            assert(c<LS->ii[n-1]);
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

    adjust_p1_p2(LS);

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#initialize rightnormed basis elements: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}


