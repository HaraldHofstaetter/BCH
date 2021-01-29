#include"bch.h"
#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<stdbool.h>
#include<assert.h>

extern unsigned int VERBOSITY_LEVEL;


#include"khash.h"
KHASH_MAP_INIT_STR(str_int, int)
KHASH_MAP_INIT_INT(LinComb, int64_t)
KHASH_MAP_INIT_STR(str_LinComb, khash_t(LinComb) *)


typedef struct magma_element_t {
    uint8_t deg;
    uint8_t g;
    struct magma_element_t *l;
    struct magma_element_t *r;
} magma_element_t;


static magma_element_t* gen(uint8_t g) {
    magma_element_t *m = malloc(sizeof(magma_element_t));
    m->deg = 1;
    m->g = g;
    m->l = NULL;
    m->r = NULL;
    return m;
}

static magma_element_t* bracket(magma_element_t *l, magma_element_t *r) {
    magma_element_t *m = malloc(sizeof(magma_element_t));
    m->deg = l->deg + r->deg;
    m->g = 0;
    m->l = l;
    m->r = r;
    return m;
}


static uint8_t m2str(const magma_element_t *m, uint8_t p, char *s) {
    if (m->deg==1) {
        s[p] = '0' + m->g;
        s[p+1] = '\0';
        return p+1;
    }
    else {
        //s[p] = '(';
        //p++;
        p = m2str(m->l, p, s);
        p = m2str(m->r, p, s);
        s[p] = ')';
        s[p+1] = '\0';
        return p+1;
    }
}


static void data2m(int dim, uint8_t *nn, uint32_t *p1, uint32_t *p2,  magma_element_t **H) {
    for (int i=0; i<dim; i++) {
        if (nn[i]==1) {
            H[i] = gen(i);
        }
        else {
            H[i] = bracket(H[p1[i]], H[p2[i]]);
        }
    }
}

static khash_t(str_int) *hall_inverse_table(int dim, magma_element_t **H) {
    khash_t(str_int) *HT = kh_init(str_int);
    for(int i=0; i<dim; i++) {
        magma_element_t *h = H[i];
        char s[3*h->deg-2+1];
        m2str(h, 0, s);
        int absent;
        khint_t k = kh_put(str_int, HT, s, &absent);
        kh_value(HT, k) = i;
        if (absent) kh_key(HT, k) = strdup(s);
    }
    return HT;
}



static int hall_data(int K, int N, int size, uint8_t **_nn, uint32_t **_p1, uint32_t **_p2,  
        magma_element_t **H, khash_t(str_int) **HT) {
    /* METHOD: Algorithm 1 in Section II.B of
     * F. Casas, A. Murua, An efficient algoritzhm for computing the 
     * Baker-Campbell-Hausdorff series and some of its applications,
     * J. Math. Phys. 50, 033513 (2009).
     */
    uint8_t *nn = malloc(size*sizeof(uint8_t));
    uint32_t *p1 = malloc(size*sizeof(uint32_t));
    uint32_t *p2 = malloc(size*sizeof(uint32_t));
    for (int i=0; i<K; i++) {
        p1[i] = i;
        p2[i] = 0;
        nn[i] = 1;
    }
    int i = K;
    for (int n=2; n<=N; n++) {
        for (int j=0; (j<i) && (nn[j]<n); j++) {
            for (int k=j+1; (k<i) && (nn[k]<=n-nn[j]); k++) {
                if ((nn[j]+nn[k]==n) && j>=p2[k]) {
                    if (i>=size) {
                        size *= 2;
                        nn = realloc(nn, size*sizeof(uint8_t));
                        p1 = realloc(p1, size*sizeof(uint32_t));
                        p2 = realloc(p2, size*sizeof(uint32_t));
                    }
                    p1[i] = k;
                    p2[i] = j;
                    nn[i] = n;
                    i++;
                }
            }
        }
    }
    data2m(i, nn, p1, p2, H);
    *HT = hall_inverse_table(i, H);
    *_nn = nn;
    *_p1 = p1;
    *_p2 = p2;
    return i;
}


static void data_from_table(int dim, magma_element_t **H, khash_t(str_int) *HT,
        uint8_t **_nn, uint32_t **_p1, uint32_t **_p2) {
    uint8_t *nn = malloc(dim*sizeof(uint8_t));
    uint32_t *p1 = malloc(dim*sizeof(uint32_t));
    uint32_t *p2 = malloc(dim*sizeof(uint32_t));
    for (int i=0; i<dim; i++) {
        magma_element_t *h = H[i];
        nn[i] = h->deg;
        if (nn[i]==1) {
            p1[i] = h->g;
            p2[i] = 0;
        }
        else {
            char s[3*h->deg-2+1];
            m2str(h->l, 0, s);
            khint_t k = kh_get(str_int, HT, s);
            assert (k != kh_end(HT));
            p1[i] = kh_value(HT, k); 
            m2str(h->r, 0, s);
            k = kh_get(str_int, HT, s);
            assert (k != kh_end(HT));
            p2[i] = kh_value(HT, k); 
        }
    }
    *_nn = nn;
    *_p1 = p1;
    *_p2 = p2;
}


static void foliage(const magma_element_t *m, uint8_t p, char *f) {
    if (m->deg==1) {
        f[p] = '0' + m->g;
        f[p+1] = '\0';
    }
    else {
        foliage(m->l, p, f);
        foliage(m->r, p+m->l->deg, f);
    }
}

/* hcmp_0, hcmp_1, hcmp_2 are examples of Hall orders. 
 * n1, n2 ... degrees of first and second hall elements 
 * f1, f2 ... foliages of first and second hall elements 
 * Note that it shold not be expected that f1[n1]=='\0'
 * or f2[n2]=='\0', so that comparisons of these strings
 * should be done with strncmp rather than strcmp.
 */

int hcmp_0(int n1, const char *f1, int n2, const char *f2) {
    int n = n1<n2 ? n1 : n2;
    int r = strncmp(f2, f1, n); 
    if (r==0) {
        if (n2<n1) {
            return -1;
        }
        else if (n2>n1) {
            return +1;
        }
        else {
            return 0;
        }
    }
    else {
        return r;
    }
}


int hcmp_1(int n1, const char *f1, int n2, const char *f2) {
    if (n1<n2) {
        return -1;
    }
    else if (n1>n2) {
        return +1;
    }
    else { /* n1==n2 */
        return strncmp(f1, f2, n1);
    }
}


int hcmp_2(int n1, const char *f1, int n2, const char *f2) {
    if (n1<n2) {
        return -1;
    }
    else if (n1>n2) {
        return +1;
    }
    else {
        int c1 = 0;
        int c2 = 0;
        for (int i=0; i<n1; i++) {
            if (f1[i]=='0') {
                c1++;
            }
            if (f2[i]=='0') {
                c2++;
            }
        }
        if (c1<c2) {
            return -1;
        }
        else if (c1>c2) {
            return +1;
        }
        else {
            return strcmp(f1, f2);
        }
    }
    assert(0); /* never reach this place */
}


static bool ishall(magma_element_t *m, char *f, int (*hcmp)(int n1, const char *f1, int n2, const char *f2)) {
    if (m->deg==1) {
        return true;
    }
    else if (!ishall(m->l, f, hcmp) 
          || !ishall(m->r, f+m->l->deg, hcmp) 
          || (hcmp(m->l->deg, f, m->r->deg, f+m->l->deg)<=0)) {
        return  false;
    }
    else if (m->l->deg==1) {
        return true;
    }
    else {
        return hcmp(m->l->r->deg, f+m->l->l->deg, m->r->deg, f+m->l->deg)<=0;
    }
}


void qsort_r(void *base, size_t nmemb, size_t size,
                  int (*compar)(const void *, const void *, void *),
                  void *arg);


static int cmp_for_qsort_r(const void *_m1, const void *_m2, void *_hcmp) {
    int (*hcmp)(int n1, const char *f1, int n2, const char *f2) = _hcmp;
    const magma_element_t *m1 = *(magma_element_t * const *) _m1;
    const magma_element_t *m2 = *(magma_element_t * const *) _m2;
    int n1 = m1->deg; 
    int n2 = m2->deg;
    char f1[n1+1]; 
    char f2[n2+1]; 
    foliage(m1, 0, f1);
    foliage(m2, 0, f2);
    return  hcmp(n1, f1, n2, f2);
}


static int hall_data_from_hall_order(int K, int N, int size, 
        int (*hcmp)(int n1, const char *f1, int n2, const char *f2),
        uint8_t **_nn, uint32_t **_p1, uint32_t **_p2,
        magma_element_t **H, khash_t(str_int) **HT) {
    for (int i=0; i<K; i++) {
        H[i] = gen(i);
    }
    int n = 2;
    int k = K;
    char f[N+1];
    while (n<=N) {
        int k0 = k;
        for(int i=0; i<k0; i++) {
            foliage(H[i], 0, f);
            for(int j=0; j<k0; j++) {
                if (H[i]->deg + H[j]->deg==n) {
                    magma_element_t *h = bracket(H[i], H[j]);
                    foliage(H[j], H[i]->deg, f);
                    if (ishall(h, f, hcmp)) {
                        H[k] = h;
                        k++;
                    }
                    else {
                        free(h);
                    }
                }
            }
        }
        n++;
    }
    qsort_r(H, k, sizeof(magma_element_t*), cmp_for_qsort_r, (void *) hcmp);
    *HT = hall_inverse_table(k, H);
    data_from_table(k, H, *HT, _nn, _p1, _p2);
    return k;
}


static khash_t(LinComb)* rewrite_magma_element(magma_element_t *m, 
                      magma_element_t **H, khash_t(str_int) *HT,
                      khash_t(str_LinComb) *LT) {
    khash_t(LinComb) *R;
    char *s = malloc((3*m->deg-2+1)*sizeof(char));
    m2str(m, 0, s);
    if (LT!=NULL) {
        /* Take result from lookup table if it already contains m */
        khint_t k = kh_get(str_LinComb, LT, s);
        if (k != kh_end(LT)) {
           R = kh_value(LT, k);
           free(s);
           return R;
        }
    }
    khint_t k = kh_get(str_int, HT, s);
    /* If m is hall element, then R = 1*m */
    if (k != kh_end(HT)) {
        int64_t i = kh_value(HT, k);
        R = kh_init(LinComb);
        int absent;
        khint_t k = kh_put(LinComb, R, i, &absent);
        kh_value(R, k) = 1;
        goto exit;
    }
    assert(m->deg>1);
    khint_t kl;
    {
    char sl[3*m->l->deg-2+1];
    m2str(m->l, 0, sl);
    kl = kh_get(str_int, HT, sl);
    }
    if (kl != kh_end(HT)) {
        khint_t kr;
        {
        char sr[3*m->r->deg-2+1];
        m2str(m->r, 0, sr);
        kr = kh_get(str_int, HT, sr);
        }
        if (kr != kh_end(HT)) {
            /* m->l and m->r are both Hall elements */
            int il = kh_value(HT, kl);
            int ir = kh_value(HT, kr);
            if (ir>il) {
                /* m = (l,r) = -(r, l) = -mrev */
                magma_element_t *mrev = bracket(m->r, m->l);
                khash_t(LinComb) *S = rewrite_magma_element(mrev, H, HT, LT);
                free(mrev);
                /* compute R = -S */
                R = kh_init(LinComb); 
                for (khint_t ks=kh_begin(S); ks!=kh_end(S); ks++) {
                    if (kh_exist(S, ks)) {
                        int absent;
                        khint_t kr = kh_put(LinComb, R, kh_key(S, ks), &absent); 
                        kh_value(R, kr) = -kh_value(S, ks);
                    }
                }
                if (LT==NULL) {
                    /* If a lookup table is available, auxiliary results in the form of 
                     * linear combinations must not be destroyed! On the other hand, if
                     * no lookup table is available, they should be destroyed for a 
                     * proper cleaning up. 
                     */
                    kh_destroy(LinComb, S);
                }
                goto exit;
            }
            else if (il>ir) {
                magma_element_t *a = m->l->l;
                magma_element_t *b = m->l->r;
                magma_element_t *n1 = bracket(bracket(a, m->r), b);
                magma_element_t *n2 = bracket(a, bracket(b, m->r));
                khash_t(LinComb) *S1 = rewrite_magma_element(n1, H, HT, LT);
                khash_t(LinComb) *S2 = rewrite_magma_element(n2, H, HT, LT);
                free(n1->l);
                free(n1);
                free(n2->r);
                free(n2);
                /* compute R = S1 + S2 */
                R = kh_init(LinComb); 
                for (khint_t ks=kh_begin(S1); ks!=kh_end(S1); ks++) { /* copy S1 to R */
                    if (kh_exist(S1, ks)) {
                        int absent;
                        khint_t kr = kh_put(LinComb, R, kh_key(S1, ks), &absent); 
                        kh_value(R, kr) = kh_value(S1, ks);
                    }
                }
                for (khint_t ks=kh_begin(S2); ks!=kh_end(S2); ks++) { /* add S2 to R */
                    if (kh_exist(S2, ks)) {
                        int absent;
                        khint_t kr = kh_put(LinComb, R, kh_key(S2, ks), &absent);
                        if (absent) {
                            kh_value(R, kr) = kh_value(S2, ks);
                        }
                        else {
                            kh_value(R, kr) += kh_value(S2, ks);
                        }
                        if (kh_value(R, kr)==0) {
                            kh_del(LinComb, R, kr);
                        }
                    }
                }
                if (LT==NULL) {
                    kh_destroy(LinComb, S1);
                    kh_destroy(LinComb, S2);
                }
                goto exit;
            }
            else { 
                R = kh_init(LinComb); /* zero linear combination */
                goto exit;
            }
        }
    }
    /* at least one of m->l, m->r is not a Hall element */
    khash_t(LinComb) *U = rewrite_magma_element(m->l, H, HT, LT);
    khash_t(LinComb) *V = rewrite_magma_element(m->r, H, HT, LT);
    R = kh_init(LinComb); 
    for (khint_t ku=kh_begin(U); ku!=kh_end(U); ku++) {
        if (kh_exist(U, ku)) {
            int64_t cu = kh_value(U, ku);
            if (cu!=0) {
                for (khint_t kv=kh_begin(V); kv!=kh_end(V); kv++) {
                    if (kh_exist(V, kv)) {
                        int64_t cv = kh_value(V, kv);
                        if (cv!=0) {
                            int64_t cuv = cu*cv;
                            magma_element_t *muv = bracket(H[kh_key(U, ku)], H[kh_key(V, kv)]);
                            khash_t(LinComb) *S = rewrite_magma_element(muv, H, HT, LT);
                            free(muv);
                            for (khint_t ks=kh_begin(S); ks!=kh_end(S); ks++) {
                                if (kh_exist(S, ks)) {
                                    int absent;
                                    khint_t kr = kh_put(LinComb, R, kh_key(S, ks), &absent);
                                    int64_t cuvs = cuv*kh_value(S, ks);
                                    if (absent) {
                                        kh_value(R, kr) = cuvs;
                                    }
                                    else {
                                        kh_value(R, kr) += cuvs;
                                    }
                                    if (kh_value(R, kr)==0) {
                                        kh_del(LinComb, R, kr);
                                    }
                                }
                            }
                            if (LT==NULL) {
                                kh_destroy(LinComb, S);
                            }
                        }
                    }
                }
            }
        }
    }
    if (LT==NULL) {
         kh_destroy(LinComb, U);
         kh_destroy(LinComb, V);
    }
exit: ;
    if (LT!=NULL) { /* put result into lookup table */
        int absent;
        k = kh_put(str_LinComb, LT, s, &absent);
        kh_value(LT, k) = R;
    } else {
        free(s);
    }
    return R;
}



void convert_lyndon_to_hall_lie_series(lie_series_t *LS, lie_series_t *HS,
        int (*hcmp)(int n1, const char *f1, int n2, const char *f2)) {
    double t0 = tic();
    int N = LS->N;
    HS->N = LS->N;
    HS->K = LS->K;
    HS->denom = LS->denom;
    HS->W = NULL;
    HS->R = NULL;
    HS->ii = NULL;
    HS->c = calloc(LS->dim, sizeof(INTEGER));

    magma_element_t **H = malloc(LS->dim*sizeof(magma_element_t*));
    khash_t(str_int) *HT;

    if  (hcmp==NULL) {
        HS->dim = hall_data(LS->K, LS->N, LS->dim, &HS->nn, &HS->p1, &HS->p2, H, &HT);
    }
    else {
        HS->dim = hall_data_from_hall_order(LS->K, LS->N, LS->dim, hcmp, &HS->nn, &HS->p1, &HS->p2, H, &HT);
    }
    assert(HS->dim==LS->dim);

    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#init Hall basis: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }

    {
        khint_t k = kh_get(str_int, HT, "0");
        assert(k != kh_end(HT));
        int i = kh_value(HT, k);
        HS->c[i] = LS->c[0];
    }

    magma_element_t **L = malloc(LS->dim*sizeof(magma_element_t*));
    data2m(LS->dim, LS->nn, LS->p1, LS->p2, L);

    size_t i1 = LS->ii[N-1];
    size_t i2 = LS->ii[N]-1;
    uint32_t *DI  = multi_degree_indices(LS->K, LS->dim, LS->W, LS->nn);
    size_t h1 = DI[i1];
    size_t h2 = DI[i2];
    #pragma omp parallel for schedule(dynamic,1)
    for (int h=h1; h<=h2; h++) { /* over all multi-degrees */
         khash_t(str_LinComb) *LUT = kh_init(str_LinComb);
         //kh_resize(str_LinComb, LUT, 300000);
         for (int i=i1; i<=i2; i++) {
            if (DI[i]==h) {
                size_t I[N];
                size_t nr = get_right_factors(i, I, N, LS->p1, LS->p2);
                for (int r=0; r<=nr; r++) {
                    int ii = I[r];
                    if (LS->c[ii]!=0) {
                        khash_t(LinComb) *R = rewrite_magma_element(L[ii], H, HT, LUT);
                        for (khint_t k=kh_begin(R); k!=kh_end(R); k++) {
                            if (kh_exist(R, k)) {
                                int64_t v = kh_value(R, k);
                                if (v!=0) {
                                    int j = kh_key(R, k);
                                    HS->c[j] += v*LS->c[ii];
                                }
                            }
                        }
                    }
                }
            }
        }
        for (khint_t k=kh_begin(LUT); k!=kh_end(LUT); k++) {
            if (kh_exist(LUT, k)) {
                free((char*) kh_key(LUT, k));
                kh_destroy(LinComb, kh_value(LUT, k));
            }
        }
        kh_destroy(str_LinComb, LUT);
    }

    free(DI);
    for (int i=0; i<LS->dim; i++) {
        free(H[i]);
        free(L[i]);
    }
    free(H);
    free(L);
    for (khint_t k=kh_begin(HT); k!=kh_end(HT); k++) {
        if (kh_exist(HT, k)) {
            free((char*) kh_key(HT, k));
        }
    }
    kh_destroy(str_int, HT);


    if (VERBOSITY_LEVEL>=1) {
        double t1 = toc(t0);
        printf("#convert from Lyndon to Hall Lie series: time=%g sec\n", t1);
        if (VERBOSITY_LEVEL>=2) {
            fflush(stdout);
        }
    }
}


