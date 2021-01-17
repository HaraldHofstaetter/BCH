#include"bch.h"
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

extern unsigned int VERBOSITY_LEVEL;


void init_hall(lie_series_t *LS, int basis) {
    /* METHOD: Algorithm 1 in Section II.B of
     * F. Casas, A. Murua, An efficient algoritzhm for computing the 
     * Baker-Campbell-Hausdorff series and some of its applications,
     * J. Math. Phys. 50, 033513 (2009).
     */
    if (basis==LYNDON_AS_HALL_BASIS) {
        /* leave p1, p2 as they are */
        return;
    }
    if (basis==REVERSE_LYNDON_AS_HALL_BASIS) {
        /* swap p1, p2 */
        uint32_t *h = LS->p1;
        LS->p1 = LS->p2;
        LS->p2 = h;
        return;
    }
    uint32_t *p1, *p2;
    if (basis==REVERSE_HALL_BASIS) {
        p1 = LS->p2;
        p2 = LS->p1;
    }
    else {
        p1 = LS->p1;
        p2 = LS->p2;
    }
    uint8_t *nn= LS->nn;
    for (int i=0; i<LS->K; i++) {
        p1[i] = i;
        p2[i] = 0;
        assert(nn[i]==1);
    }
    int i = LS->K;
    for (int n=2; n<=LS->N; n++) {
        for (int j=0; j<i; j++) {
            for (int k=j+1; k<i; k++) {
                if ((nn[j]+nn[k]==n) && j>=p2[k]) {
                    p1[i] = k;
                    p2[i] = j;
                    assert(nn[i]==n);
                    i++;
                }
            }
        }
    }       
    assert(i==LS->dim);
}


static uint32_t* hall_multi_degree_indices(size_t K, size_t dim, uint32_t *p1, uint32_t *p2) {
    uint32_t *DI = malloc(dim*sizeof(uint32_t));
    uint8_t *h = malloc(K*dim*sizeof(uint8_t));
    for (int i=0; i<K; i++) {
        for (int j=0; j<K; j++) {
            h[i*K+j] = 0;
        }
        h[i*K+i] = 1;
        DI[i] = tuple_index(K, h+K*i);
    }
    for (int i=K; i<dim; i++) {
        for (int j=0; j<K; j++) {
            h[i*K+j] = h[p1[i]*K+j] + h[p2[i]*K+j];
        }
        DI[i] = tuple_index(K, h+K*i);
    }
    free(h);
    return DI;
}


static int compare_w1w2_w2w1(int n1, uint8_t w1[], int n2, uint8_t w2[]) {
    int i = 0;
    while (i<n1+n2) {
        uint8_t c12 = i<n1 ? w1[i] : w2[i-n1];
        uint8_t c21 = i<n2 ? w2[i] : w1[i-n2];
        if (c12<c21) {
            return -1;
        }
        else if (c12>c21) {
            return +1;
        }
        i++;
    }
    return 0;
}

static void leading_word(int K, int i, uint8_t w[], uint8_t *nn, uint32_t *p1, uint32_t *p2) {
    if (i<K) {
        w[0] = i;
        return;
    }
    int n1 =nn[p1[i]];
    int n2 =nn[p2[i]];
    uint8_t w1[n1];
    uint8_t w2[n2];
    leading_word(K, p1[i], w1, nn, p1, p2);
    leading_word(K, p2[i], w2, nn, p1, p2);
    int c = compare_w1w2_w2w1(n1, w1, n2, w2);
    int k=0;
    if (c<0) {
        for (int j=0; j<n1; j++) {
            w[k] = w1[j];
            k++;
        }
        for (int j=0; j<n2; j++) {
            w[k] = w2[j];
            k++;
        }
    }
    else {
        for (int j=0; j<n2; j++) {
            w[k] = w2[j];
            k++;
        }
        for (int j=0; j<n1; j++) {
            w[k] = w1[j];
            k++;
        }
    }
    // assert(k==nn[i]);
}


static size_t find_smallest_lyndon_word_index(uint32_t *WI, size_t l, size_t r, size_t wi) {
    /* finds smallest index not less than wi in the sorted list of indices WI. 
     * Start search at position l and stop it at position r. 
     * METHOD: binary search, 
     * See also: lyndon.c/find_lyndon_word_index
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
    return l;
}


static void sortperm(int n, uint32_t a[], uint32_t p[]) {
    /* METHOD: bubble sort */
    for (int i=0; i<n; i++) {
        p[i] = i;
    }
    int swapped = 1;
    while (swapped) {
        swapped = 0;
        for (int i=0; i<n-1; i++) {
            if (a[p[i]]>a[p[i+1]]) {
                uint32_t h = p[i];
                p[i] = p[i+1];
                p[i+1] = h;
                swapped = 1;
            }
        }
        n--;
    }
}


static inline INT_FF_LU_T IABS(INT_FF_LU_T x) {
    return x>=0 ? x : -x;
}


/* W. Zhou, D.J. Jeffrey, Fraction-free matrix factors: new forms for 
 * LU and QR factors, Front. Comput. Sci. China 2 (1) (2008)1–13.
 *
 * David Dureisseix. Generalized fraction-free LU factorization for 
 * singular systems with kernel extraction. Linear Algebra and its 
 * Applications, Elsevier, 2012, 436 (1), pp.27-40.
 */


static void fraction_free_lu(int n, INT_FF_LU_T *A, uint32_t *p) {
    for (int i=0; i<n; i++) {
        p[i] = i;
    }
    INT_FF_LU_T oldpivot = 1;
    for (int k=0; k<n; k++) {
        INT_FF_LU_T pivot = INT_FF_LU_MAX;
        int kpivot = n;
        for (int i=k; i<n; i++) { /* search for smallest nonzero pivot */
            if ((A[i+n*k]!=0) && (IABS(A[i+n*k])<IABS(pivot))) {
                kpivot = i;
                pivot = A[i+n*k];
                if (IABS(pivot)==1) { /* pivot already as small as possible */
                    break;
                }
            }
        }
        if (kpivot>=n) {
            fprintf(stderr, "ERROR: fraction-free LU factorization does not exist\n"); 
            exit(EXIT_FAILURE);
        }
        if (kpivot!=k) { /* swap k-th and kpivot-th row */
            for (int i=0; i<n; i++) {
                INT_FF_LU_T h = A[k+n*i];
                A[k+n*i] = A[kpivot+n*i];
                A[kpivot+n*i] = h;
            }
            uint32_t h = p[k];
            p[k] = p[kpivot];
            p[kpivot] = h;
        }
        for (int i=k+1; i<n; i++) {
            INT_FF_LU_T Aik = A[i+n*k];
            if (!((Aik==0) && (pivot==oldpivot))) {
                if (oldpivot==1) { /* avoid expensive div operation */
                    for (int j=k+1; j<n; j++) {
                        A[i+n*j] = pivot*A[i+n*j] - A[k+n*j]*Aik;
                    }
                }
                else if (oldpivot==-1) { /* avoid expensive div operation */
                    for (int j=k+1; j<n; j++) {
                        A[i+n*j] = -pivot*A[i+n*j] + A[k+n*j]*Aik;
                    }
                }
                else {
                    for (int j=k+1; j<n; j++) {
                        A[i+n*j] = (pivot*A[i+n*j] - A[k+n*j]*Aik)/oldpivot; /* division is exact without remainder */
                    }
                }
            }
        }
        oldpivot = pivot;
    }
}


static INT_FF_LU_T fraction_free_lu_solve(int n, INT_FF_LU_T *A, INTEGER *x) {
    INT_FF_LU_T oldpivot = 1;
    for (int k=0; k<n-1; k++) {
        INT_FF_LU_T pivot = A[k+n*k];
        for (int i=k+1; i<n; i++) {
            x[i] = (pivot*x[i] - A[i+n*k]*x[k])/oldpivot;
        }
        oldpivot = pivot;
    }

    INT_FF_LU_T d = A[n-1+n*(n-1)];
    for (int i=n-1; i>=0; i--) {
        INTEGER h = 0;
        for (int j=i+1; j<n; j++) {
            h += A[i+j*n]*x[j];
        }
        x[i] = (d*x[i]-h)/A[i+n*i];
    }
    return d;
}

void convert_to_hall_lie_series(lie_series_t *LS, int N, int odd_orders_only) {
    double t0 = tic();
    uint32_t *WI = malloc(LS->dim*sizeof(uint32_t));
    for (int i=0; i<LS->dim; i++) {
        WI[i] = word_index(LS->K, LS->W[i], 0, LS->nn[i]-1);
    }
    uint32_t *DI_lyndon  = multi_degree_indices(LS->K, LS->dim, LS->W, LS->nn);
    uint32_t *DI_hall = hall_multi_degree_indices(LS->K, LS->dim, LS->p1, LS->p2);
    INTEGER *c_hall =malloc(LS->dim*sizeof(INTEGER));
    for (int i=0; i<=LS->K; i++) {
        c_hall[i] = LS->c[i];
    }

    for (int n=2; n<=N; n++) { /* over all word sizes */
    if ((!odd_orders_only)||(n&1)) {
        size_t i1 = LS->ii[n-1];
        size_t i2 = LS->ii[n]-1;
        size_t h1 = DI_lyndon[i1];
        size_t h2 = DI_lyndon[i2];
        #pragma omp parallel for schedule(dynamic,1)
        for (int h=h1; h<=h2; h++) { /* over all multi-degrees */
            /* get dimension */
            int m=0;
            for (int j=i1; j<=i2; j++) {
                if (DI_lyndon[j]==h) {
                    m++;
                }
            }
            if (m==0) {
                continue;
            }

            uint32_t *WI1 = malloc(m*sizeof(uint32_t));
            uint32_t *I = malloc(m*sizeof(uint32_t));
            uint32_t *J = malloc(m*sizeof(uint32_t));
            int k=0;
            int l=0;
            for (int j=i1; j<=i2; j++) {
                if (DI_lyndon[j]==h) {
                    I[k] = j;
                    WI1[k] = WI[j];
                    k++;
                }
                if (DI_hall[j]==h) {
                    J[l] = j;
                    l++;
                }
            }
            assert(l==m);
            uint32_t* LWI = malloc(m*sizeof(uint32_t));
            uint8_t w[n];
            for (int j=0; j<m; j++) {
                leading_word(LS->K, J[j], w, LS->nn, LS->p1, LS->p2);
                int i = word_index(LS->K, w, 0, n-1);
                LWI[j] = find_smallest_lyndon_word_index(WI1, 0, m-1, i);
            }
            uint32_t* p0 = malloc(m*sizeof(uint32_t));
            sortperm(m, LWI, p0);

            P_t *P = P_init(LS->K, N, 2*m);
            uint32_t *r = malloc(m*sizeof(uint32_t));
            int stop = 0;
            for (int j=0; j<m; j++) {
                r[j] = P_append(P, J[j], 0, LS->p1, LS->p2, LS->nn);
                stop = r[j] > stop ? r[j] : stop;
            }
            int32_t *X = malloc((P->len)*sizeof(int32_t));

            /* set up matrix */
            INTEGER *x = calloc(m, sizeof(INTEGER));
            INT_FF_LU_T *A = calloc(m*m, sizeof(INT_FF_LU_T));
            for (int i=0; i<m; i++) {
                //stop = r[i] > stop ? r[i] : stop;
                P_run(X, P, LS->W[I[m-i-1]], stop);
                for (int j=0; j<m; j++) {
                    A[i+j*m] = X[r[p0[m-j-1]]];
                }
            }

            uint32_t* p1 = malloc(m*sizeof(uint32_t));
            fraction_free_lu(m, A, p1);

            /* set up righthand side */
            for (int i=0; i<m; i++) {
                x[i] = LS->c[I[m-p1[i]-1]];
            }

            INT_FF_LU_T det = fraction_free_lu_solve(m, A, x);
            assert(IABS(det)==1);

            /* copy result */
            if (det==1) {
                for (int j=0; j<m; j++) {
                     c_hall[J[p0[m-j-1]]] = x[j];
                }
            }
            else { /* det==-1 */
                for (int j=0; j<m; j++) {
                     c_hall[J[p0[m-j-1]]] = -x[j];
                }
            }

            free(WI1);
            free(I);
            free(J);
            free(LWI);
            free(p0);
            free(p1);
            free(r);
            free(X);
            free(x);
            free(A);
        }
    }
    else {
        for (int j=LS->ii[n-1]; j<=LS->ii[n]-1; j++) {
            LS->c[j] = 0;
        }
    }
    }
    free(LS->c);
    LS->c = c_hall;
    free(WI);
    free(DI_lyndon);
    free(DI_hall);
    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#convert to Lie series: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


