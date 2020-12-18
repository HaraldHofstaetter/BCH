static int coeff_word_in_rightnormed(generator_t w[], generator_t c[], int l1, int r1, int l2) {
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


static void convert_to_rightnormed_lie_series(int N, INTEGER c[], int odd_orders_only) {
    double t0 = tic();
    for (int n=2; n<=N; n++) { /* over all word sizes */
    if ((!odd_orders_only)||(n&1)) {
        size_t i1 = ii[n-1];
        size_t i2 = ii[n]-1;
        uint32_t *DI  = multi_degree_indices( K, N_LYNDON, W, nn);
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
                            A[ii+jj*m] = coeff_word_in_rightnormed(W[i], R[j], 0, n-1, 0);
                            ii++;
                        }
                    }
                    x[jj] = c[j];
                    jj++;
                }
            }

            integer_lu(m, A);
            integer_lu_solve(m, A, x);

            /* copy result */
            jj=0;
            for (int j=i1; j<=i2; j++) {
                if (DI[j]==h) {
                    c[j] = x[jj];
                    jj++;
                }
            }

            free(x);
            free(A);
        }
        free(DI);
    }
    else {
        for (int j=ii[n-1]; j<=ii[n]-1; j++) {
            c[j] = 0;
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


static void compute_rightnormed_BCH_terms_of_even_orders(INTEGER c[]) {
    double t0 = tic();

    for(int n=2; n<=N; n+=2) {
        #pragma omp parallel for schedule(dynamic,256)
        for (int i=ii[n-1]; i<=ii[n]-1; i++) {
            c[i] = 0;
            int k=0;
            int l=0;
            while (R[i][k]==K-1) {
                k += 1;
                if (k&1) {
                    int q = ii[n-1-k];
                    for (; q<=ii[n-k]-1; q++) {
                        int m=0;
                        for (; (m<n-k) && (R[q][m]==R[i][k+m]) ; m++) {}
                        if (m==n-k) {
                            break;
                        }
                    }
                    if (q>ii[n-k]) {
                        fprintf(stderr, "ERROR: basis element not found in compute_rightnormed_BCH_terms_of_even_orders");
                        exit(EXIT_FAILURE);
                    }
                    INTEGER d = c[q]/beta_den[l];
                    if (d*beta_den[l]!=c[q]) {
                        fprintf(stderr, "ERROR: divisibility check failed in compute_rightnormed_BCH_terms_of_even_orders");
                        exit(EXIT_FAILURE);
                    }
                    c[i] -= beta_num[l]*d; 
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


static void init_rightnormed(void) {
    double t0 = tic();
    size_t mem_len = 0;
    for (int n=1; n<=N; n++) {
        mem_len += n*(ii[n]-ii[n-1]);
    }

    R = malloc(N_LYNDON*sizeof(generator_t *));
    R[0] = malloc(mem_len*sizeof(generator_t)); 
    for (int i=1; i<N_LYNDON; i++) {
        R[i] = R[i-1] + nn[i-1];
    }

    #pragma omp for schedule(dynamic,256) 
    for (int i=0; i<N_LYNDON; i++) {
        lyndon2rightnormed(nn[i], W[i], R[i]);
    }

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#init rightnormed basis elements: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


