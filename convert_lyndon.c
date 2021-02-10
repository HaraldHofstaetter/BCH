#include"bch.h"
#include <stdio.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __GNUC__ 
#define SIMD_VECTORIZED 1
// #define USE_SIMD_INTRINSICS 1
#endif


#include"khash.h"
KHASH_MAP_INIT_INT64(P_Dict, uint32_t)      // instantiate structs and methods
    
typedef khash_t(P_Dict) PH;


typedef struct P_line_t {
    uint32_t a11;
    uint32_t a12;
    uint32_t a21;
    uint32_t a22;
} P_line_t;


typedef struct P_t {
    P_line_t *L;
    uint32_t len;
    uint32_t maxlen;
    uint8_t n;
    uint8_t K;
    void *H;
} P_t;


static P_t *P_init(uint8_t K, uint8_t n, uint32_t len) {
    P_t *P = malloc(sizeof(P_t));
    P->H = kh_init(P_Dict);  // allocate hash table
    P->n = n;
    P->K = K;
    if (len<K*n) {
        len = K*n;
    }
    P->L = malloc(sizeof(P_line_t)*len);
    P->maxlen = len;
    khint_t k;
    int l=0;
    for (int i=0; i<K; i++) {
        for (int j=0; j<n; j++) {
            P->L[l].a11 = 0;
            P->L[l].a12 = 0;
            P->L[l].a21 = 0;
            P->L[l].a22 = 0;
            uint64_t key = (((uint64_t) i)<< 8) | j;
            int absent;
            k = kh_put(P_Dict, (PH*) P->H, key, &absent);  // insert a key to the hash table
            kh_val((PH*) P->H, k) = l;
            l++;
        }    
    }
    P->len = l;
    return P;
}


static void P_free(P_t *P) {
    free(P->L);
    kh_destroy(P_Dict, (PH*) P->H);  // deallocate hash table
    free(P);
}


static uint32_t P_append(P_t *P, uint32_t i, uint8_t l, uint32_t *p1, uint32_t *p2, uint8_t* nn) {
    uint64_t key = (((uint64_t) i)<< 8) | l;
    khint_t k = kh_get(P_Dict, (PH*) P->H, key);  // query the hash table
    if (k == kh_end((PH*) P->H)) {                // test if the key is missing
        if (nn[i]>1) {
            uint32_t        a11 = P_append(P, p1[i], l,           p1, p2, nn);
            uint32_t        a12 = P_append(P, p1[i], l+nn[p2[i]], p1, p2, nn);
            uint32_t        a21 = P_append(P, p2[i], l,           p1, p2, nn);
            uint32_t        a22 = P_append(P, p2[i], l+nn[p1[i]], p1, p2, nn);
            if (P->len>=P->maxlen) {
                P->maxlen *= 2;
                P->L = realloc(P->L, sizeof(P_line_t)*P->maxlen);
            }
            P->L[P->len].a11 = a11;
            P->L[P->len].a12 = a12;
            P->L[P->len].a21 = a21;
            P->L[P->len].a22 = a22;
        }
        int absent;
        k = kh_put(P_Dict, (PH*) P->H, key, &absent);  // insert a key to the hash table
        kh_val((PH*) P->H, k) = P->len;
        P->len++;
        return (P->len-1);
    }
    else {
        return kh_val((PH*) P->H, k);
    }
}


#ifndef SIMD_VECTORIZED
static void P_run(int32_t *X, P_t *P, uint8_t w[], uint32_t stop) { 
    if (stop>=P->len) {
        stop = P->len-1;
    }
    P_line_t *L = P->L;
    for (int k=0; k<P->K; k++) {
        for (int i=0; i<P->n; i++) {
            X[k*P->n+i] = w[i]==k ? 1 : 0;
        }
    }
    for (int p=P->K*P->n; p<=stop; p++) {
        X[p] = X[L[p].a11]*X[L[p].a22] - X[L[p].a12]*X[L[p].a21];
    }
}
#endif



#ifdef SIMD_VECTORIZED
#ifdef USE_SIMD_INTRINSICS    
#include <smmintrin.h>
#endif

typedef int32_t v4int32_t __attribute__ ((vector_size(16), aligned(16))); 

static void  P_run_4(v4int32_t* X0, P_t *P, uint8_t w0[], uint8_t w1[], uint8_t w2[], uint8_t w3[], uint32_t stop) {
    v4int32_t *X = __builtin_assume_aligned (X0, 16);
    if (stop>=P->len) {
        stop = P->len-1;
    }
    P_line_t *L = P->L;
    for (int k=0; k<P->K; k++) {
        for (int i=0; i<P->n; i++) {
            X[k*P->n+i][0] = w0[i]==k ? 1 : 0;
            X[k*P->n+i][1] = w1[i]==k ? 1 : 0;
            X[k*P->n+i][2] = w2[i]==k ? 1 : 0;
            X[k*P->n+i][3] = w3[i]==k ? 1 : 0;
        }
    }
#ifndef USE_SIMD_INTRINSICS    
    for (int p=P->K*P->n; p<=stop; p++) {
       X[p] = X[L[p].a11]*X[L[p].a22] - X[L[p].a12]*X[L[p].a21];
    }
#else
    for (int p=P->K*P->n; p<=stop; p++) {
        __m128i x11 = _mm_load_si128( (__m128i*) X + L[p].a11 );
        __m128i x12 = _mm_load_si128( (__m128i*) X + L[p].a12 );
        __m128i x21 = _mm_load_si128( (__m128i*) X + L[p].a21 );
        __m128i x22 = _mm_load_si128( (__m128i*) X + L[p].a22 );
        __m128i x11x22 = _mm_mullo_epi32(x11, x22);
        __m128i x12x21 = _mm_mullo_epi32(x12, x21);
        __m128i y = _mm_sub_epi32(x11x22, x12x21);
        _mm_store_si128( (__m128i*) X +p, y);  
    }
    
#endif 
}
#endif


void convert_to_lie_series(lie_series_t *LS, int N) {
    if (get_verbosity_level()>=2) {
#ifdef _OPENMP
        printf("# degree     #basis        time thread\n");
#else
        printf("# degree     #basis        time\n");
#endif
    }
    double t0 = tic();

    size_t i1 = LS->ii[N-1];
    size_t i2 = LS->ii[N]-1;

    uint32_t *DI  = multi_degree_indices( LS->K, LS->dim, LS->W, LS->nn);

    size_t h1 = DI[i1];
    size_t h2 = DI[i2];

    double h_time[h2-h1+1];
    double h_time1[h2-h1+1];
    double h_time2[h2-h1+1];
    int h_n[h2-h1+1];
#ifdef _OPENMP
    int h_thread[h2-h1+1];
#endif
    #pragma omp parallel 
    {
    int *jj = calloc(LS->dim, sizeof(int));  // LS->dim far too large upper bound
    size_t JW[N];
    size_t JB[N];

    /* Note: We choose schedule(dynamic, 1) because each
     * iteration of the loop is associated with a specific 
     * multi degree index, and the work to be done varies widely
     * for different multi degree indices. 
     */
    #pragma omp for schedule(dynamic,1) 
    for (int h=h1; h<=h2; h++) { /* over all multi-degrees */
        int k = h-h1;
        h_time[k] = tic();
        h_n[k] = 0;
#ifdef _OPENMP
        h_thread[k] = omp_get_thread_num();
#endif

        int jj_max = 0; 
        for (int i=i1; i<=i2; i++) {
            if (DI[i]==h) {
                jj[jj_max] = i;
                jj_max++;
            }
        }

        P_t *P = P_init(LS->K, N, 2*jj_max);
        uint32_t *r = malloc(jj_max*sizeof(uint32_t));

        for (int y=0; y<jj_max; y++) {
            int j = jj[y];
            size_t kB = get_right_factors(j, JB, N, LS->p1, LS->p2);
            r[y] = P_append(P, JB[kB], kB, LS->p1, LS->p2, LS->nn);
        }
#ifdef SIMD_VECTORIZED
        h_time1[k] = tic();
        v4int32_t *X = aligned_alloc(16, (P->len)*sizeof(v4int32_t));
        h_time1[k] = toc(h_time1[k]); 

        h_time2[k] = tic();
        int stop = 0;
        for (int x=0; x<jj_max; x+=4) {
            int i[4];
            for (int s=0; s<4; s++) {
                i[s] = x+s<jj_max ? jj[x+s] : jj[jj_max-1];
                stop = (x+s < jj_max) && (r[x+s] > stop) ? r[x+s] : stop;
            }
            
            h_n[k]++; // TODO: adapt this to vectorized version

            P_run_4(X, P, LS->W[i[0]], LS->W[i[1]], LS->W[i[2]], LS->W[i[3]], stop);

            for (int s=0; s<4 && x+s<jj_max; s++) {
                size_t kW =  get_right_factors(i[s], JW, N, LS->p1, LS->p2);
                int lA = 0; for (;LS->W[i[s]][lA]==0; lA++);
            
                for (int y=0; y<=x+s-1; y++) {
                    int j = jj[y];
                    size_t kB = get_right_factors(j, JB, N, LS->p1, LS->p2);
                    if (lA>=kB) {
                        int d = X[r[y]][s];
                        if (d!=0) {
                            for (int k=0; k<=kB && k<=kW; k++) {
                                LS->c[JW[k]] -= d*LS->c[JB[k]];
                            }
                        }
                    }
                }
            }
            
        }
        h_time2[k] = toc(h_time2[k]); 

#else        
        h_time1[k] = tic();
        int32_t *X = malloc((P->len)*sizeof(int32_t));
        h_time1[k] = toc(h_time1[k]); 

        h_time2[k] = tic();
        int stop = 0;
        for (int x=0; x<jj_max; x++) {
            int i = jj[x];
            stop = r[x] > stop ? r[x] : stop;
            h_n[k]++;

            size_t kW = get_right_factors(i, JW, N, LS->p1, LS->p2);
            uint8_t *w = LS->W[i];

            P_run(X, P, w, stop);

            int lA = 0; 
            for (;w[lA]==0; lA++) ;

            for (int y=0; y<=x-1; y++) {
                int j = jj[y];
                size_t kB = get_right_factors(j, JB, N, LS->p1, LS->p2);
                if (lA>=kB) {
                    int d = X[r[y]];
                    if (d!=0) {
                        for (int k=0; k<=kB && k<=kW; k++) {
                            LS->c[JW[k]] -= d*LS->c[JB[k]];
                        }
                    }
                }
            }
        }
        h_time2[k] = toc(h_time2[k]); 
#endif
        free(r);
        int len = P->len;
        P_free(P);
        free(X);
        LS->R = 0;
        h_time[k] = toc(h_time[k]); 
        if (get_verbosity_level()>=2) {
#ifdef _OPENMP
            printf("#%7i %10i %11.2f %11.2f %11.2f %10i %4i\n", k+1, h_n[k], h_time1[k], h_time2[k], h_time[k], len, h_thread[k]);
#else
            printf("#%7i %10i %11.2f\n", k+1, h_n[k], h_time[k]);
#endif
            fflush(stdout);
        }
    }
    free(jj);
    }
    free(DI);

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#convert to Lie series: time=%g sec\n", t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}


/* tables beta_num[] and beta_den[]: numerators and denominators of
   the coefficients of the power series of the function f(x)=tanh(x/2), 
   i.e., beta[] = {     1/2,                                   
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

static const INTEGER H =  1000000000000000000;

/* beta_num, beta_den not static because also needed in convert_rightnormed.c */

INTEGER beta_num[16] = {1, -1, 1, -17, 31, -691, 5461, -929569, 3202291, -221930581, 4722116521, 
    -56963745931, 14717667114151, -2093660879252671, 86125672563201181, -129*H-848163681107301953};

INTEGER beta_den[16] = {2, 24, 240, 40320, 725760, 159667200, 12454041600, 20922789888000, 
    711374856192000, 486580401635328000, 102*H+181884343418880000, 12165*H+654935945871360000, 
    31022420*H+86661971968000000, 43555477801*H+673408643072000000, 17683523987479*H+403909087232000000, 
    263130836933693530*H+167218012160000000};


void compute_BCH_terms_of_even_degree_N(lie_series_t *LS) {
    double t0 = tic();
    assert(!(LS->N&1));

    #pragma omp parallel for schedule(dynamic,256)
    for (int i=LS->ii[LS->N-1]; i<=LS->ii[LS->N]-1; i++) {
        LS->c[i] = 0;
        int k = 0;
        int l = 0;
        int q = i;
        while (LS->p1[q]==0) {
            k += 1;
            q = LS->p2[q];
            if (k&1) {
                INTEGER d = LS->c[q]/beta_den[l];
                if (d*beta_den[l]!=LS->c[q]) {
                    fprintf(stderr, "ERROR: divisibility check failed in compute_BCH_terms_of_degree_N");
                    exit(EXIT_FAILURE);
                }
                LS->c[i] += beta_num[l]*d; 
                l += 1;
            }
        }
    } 

    if (get_verbosity_level()>=1) {
        double t1 = toc(t0);
        printf("#compute terms of degree %i: time=%g sec\n", LS->N, t1);
        if (get_verbosity_level()>=2) {
            fflush(stdout);
        }
    }
}

