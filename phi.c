#include"bch.h"
#include<stdlib.h>
#include<stdio.h>
#include<inttypes.h>

static const INTEGER H =  1000000000000000000;

/* FACTORIAL not static because also needed in goldberg.c */
INTEGER FACTORIAL[33] =  {             1,
                                       1,
                                       2,
                                       6,
                                      24,
                                     120,
                                     720,
                                    5040,
                                   40320,
                                  362880,
                                 3628800,
                                39916800,
                               479001600,
                              6227020800,
                             87178291200,
                           1307674368000,
                          20922789888000,
                         355687428096000,
                        6402373705728000,
                      121645100408832000,
                  2*H+432902008176640000,
                 51*H+ 90942171709440000,
               1124*H+   727777607680000,
              25852*H+ 16738884976640000,
             620448*H+401733239439360000,
           15511210*H+ 43330985984000000,
          403291461*H+126605635584000000,
        10888869450*H+418352160768000000,
       304888344611*H+713860501504000000,
      8841761993739*H+701954543616000000,
    265252859812191*H+ 58636308480000000,
   8222838654177922*H+817725562880000000,
 263130836933693530*H+167218012160000000 };


static inline void check_for_divisibility_by_int(INTEGER p, int q, INTEGER d, char *s) {
    if (q*d!=p) {
        int q1 = (q>0?q:-q)/gcd_INTEGER(p,q);
        fprintf(stderr, "PANIC: dividend not divisble by %i %s\n", q1, s);
        abort();
    }
}

static inline void check_for_divisibility_by_int64(INTEGER p, int64_t q, INTEGER d, char *s) {
    if (q*d!=p) {
        int64_t q1 = (q>0?q:-q)/gcd_INTEGER(p,q);
        fprintf(stderr, "PANIC: dividend not divisble by %" PRId64 " %s\n", q1, s);
        abort();
    }
}

static inline void check_for_divisibility_by_INTEGER(INTEGER p, INTEGER q, INTEGER d, char *s) {
    if (q*d!=p) {
        int64_t q1 = (q>0?q:-q)/gcd_INTEGER(p,q);
        fprintf(stderr, "PANIC: dividend not divisble by %" PRId64 " %s\n", q1, s);
        abort();
    }
}

int phi(INTEGER y[], int m, uint8_t w[], expr_t* ex, INTEGER v[]) {
    if (m==0) {
        return 0;
    }
    switch (ex->type) {
        case ZERO_ELEMENT: 
            return 0;
        case IDENTITY: 
            for (int j=0; j<m; j++) {
                y[j] = v[j];
            }
            return m;
        case GENERATOR: {
            int m1=0;
            for (int j=0; j<m-1; j++) {
                if (w[j]==ex->gen) {
                    y[j] = v[j+1];
                    if (y[j]!=0) {
                        m1 = j+1;
                    }
                }
                else {
                    y[j] = 0;
                }
            }
            return m1;
            }
        case SUM: { 
            INTEGER y2[m];
            for (int j=0; j<m; j++) {
                y2[j] = v[j];
            }
            int m1 = phi(y, m, w, ex->arg1, v);
            int m2 = phi(y2, m, w, ex->arg2, y2);
            if (m1<m2) {
                for (int j=0; j<m1; j++) {
                    y[j] += y2[j];
                }
                for (int j=m1; j<m2; j++) {
                    y[j] = y2[j];
                }
                return m2;
            }
            else {
                for (int j=0; j<m2; j++) {
                    y[j] += y2[j];
                }
                return m1;
            }
            } 
        case DIFFERENCE: {
            INTEGER y2[m];
            for (int j=0; j<m; j++) {
                y2[j] = v[j];
            }
            int m1 = phi(y, m, w, ex->arg1, v);
            int m2 = phi(y2, m, w, ex->arg2, y2);
            if (m1<m2) {
                for (int j=0; j<m1; j++) {
                    y[j] -= y2[j];
                }
                for (int j=m1; j<m2; j++) {
                    y[j] = -y2[j];
                }
                return m2;
            }
            else {
                for (int j=0; j<m2; j++) {
                    y[j] -= y2[j];
                }
                return m1;
            }
            } 
        case PRODUCT: {
            int md;
            if (ex->arg1->const_term.num!=0) {
                md = 0;
            }
            else {
                md = ex->arg1->mindeg;
            }
            if (md>=m) {
                return 0;
            }
            int m1 = phi(y+md, m-md, w+md, ex->arg2, v+md);
            if (m1==0) {
                return 0;
            }
            for (int j=0; j<md; j++) {
                y[j] = 0;
            }
            return phi(y, m1+md, w, ex->arg1, y);
            }
        case NEGATION: { 
            int m1 = phi(y, m, w, ex->arg1, v);
            for (int j=0; j<m1; j++) {
                y[j] = -y[j];
            }
            return m1;
            }
        case TERM: { 
            int p = ex->factor.num;
            for (int j=0; j<m; j++) {
                y[j] = p*v[j];
            }
            int m1 = phi(y, m, w, ex->arg1, y);
            int q = ex->factor.den;
            if (q!=1) {
                for (int j=0; j<m1; j++) {
                    INTEGER h = y[j];
                    INTEGER d = h/q;
                    check_for_divisibility_by_int(h, q, d, "in phi()/TERM");
                    y[j] = d;
                }
            }
            return m1;
            }
        case EXPONENTIAL: {
            INTEGER z[m];
            for (int j=0; j<m; j++) {
                z[j] = v[j];
                y[j] = v[j];
            }
            int m1 = m;
            for (int k=1; k<m; k++) {
                m1 = phi(z, m1, w, ex->arg1, z);
                if (m1==0) {
                    return m;
                }
                if (k<=20) {
                    int64_t f = FACTORIAL[k]; /* fits into int64 => faster execution expected */
                    for (int j=0; j<m1; j++) {
                        INTEGER d = z[j]/f;
                        check_for_divisibility_by_int64(z[j], f, d, "in phi()/EXPONENTIAL");
                        y[j] += d;
                    }
                }
                else {
                    INTEGER f = FACTORIAL[k];
                    for (int j=0; j<m1; j++) {
                        INTEGER d = z[j]/f;
                        check_for_divisibility_by_INTEGER(z[j], f, d, "in phi()/EXPONENTIAL");
                        y[j] += d;
                    }
                }
            }
            return m;
            }
        case LOGARITHM: {
            INTEGER z[m];
            for (int j=0; j<m; j++) {
                z[j] = v[j];
                y[j] = 0;                     
            } 
            INTEGER h[m];
            int m1 = m; 
            int mret = m;
            for (int k=1; k<m; k++) {
                for (int j=0; j<m1; j++) {
                    h[j] = z[j];
                }
                int m2 = phi(z, m1, w, ex->arg1, z);
                int m3 = 0;
                for (int j=0; j<m2; j++) {
                    z[j] -= h[j];
                    if (z[j]!=0) {
                        m3 = j+1;
                    }
                }
                for (int j=m2; j<m1; j++) {
                    z[j] = -h[j];
                    if (z[j]!=0) {
                        m3 = j+1;
                    }
                }
                if (k==1) {
                    mret = m3;
                }
                if (m3==0) {
                    return mret;
                }
                m1 = m3;
                int f = k%2 ? +k : -k; /* f = (-1)^(k+1)*k */ 
                for (int j=0; j<m1; j++) {
                    INTEGER d = z[j]/f;
                    check_for_divisibility_by_int(z[j], f, d, "in phi()/LOGARITHM");
                    y[j] += d;
                }
            }
            return mret;
            }
        default:
            fprintf(stderr, "PANIC: unknown expr type %i\n", ex->type);
            abort();
    }
}


static inline INTEGER lcm_INTEGER(INTEGER a, INTEGER b) {
    /* computes least common multiple of a and b */
    return (a/gcd_INTEGER(a,b))*b;
}


static inline INTEGER lcm1_INTEGER(INTEGER a, INTEGER b) {
    if (a==0) {
        return b;
    }
    else if (b==0) {
        return a;
    }
    else {
        return lcm_INTEGER(a, b);
    }
}


static void delta(INTEGER d[], int N, expr_t* ex) {
    switch (ex->type) {
        case ZERO_ELEMENT:
            for (int n=0; n<=N; n++) {
                d[n] = 0;
            }
            return;
        case GENERATOR:
            d[0] = 0;
            if (N>=1) {
                d[1] = 1;
            }
            for (int n=2; n<=N; n++) {
                d[n] = 0;
            }
            return;
        case IDENTITY:
            d[0] = 1;
            for (int n=1; n<=N; n++) {
                d[n] = 0;
            }
            return;
        case NEGATION:
            delta(d, N, ex->arg1);
            return;
        case SUM:
        case DIFFERENCE: {
            INTEGER h[N+1];
            delta(d, N, ex->arg1);
            delta(h, N, ex->arg2);
            for (int n=0; n<=N; n++) {
                d[n] = lcm1_INTEGER(d[n], h[n]);
            }
            return;
            }
        case TERM:
            if (ex->factor.num==0) {
                for (int n=0; n<=N; n++) {
                    d[n] = 0;
                }
                return;
            }
            delta(d, N, ex->arg1);
            for (int n=0; n<=N; n++) {
                INTEGER x = ex->factor.den*d[n];
                d[n] = x/gcd_INTEGER(x, ex->factor.num);
            }
            return;
        case PRODUCT: {
            INTEGER h1[N+1];
            INTEGER h2[N+1];
            delta(h1, N, ex->arg1);
            delta(h2, N, ex->arg2);
            for (int n=0; n<=N; n++) {
                INTEGER m = 0;
                for (int k=0; k<=n; k++) {
                    m = lcm1_INTEGER(m, h1[k]*h2[n-k]);
                }
                d[n] = m;
            }
            return;
            }
        case EXPONENTIAL: {
            INTEGER a[N+1];
            INTEGER h[N+1];
            INTEGER h1[N+1];
            delta(a, N, ex->arg1);
            if (a[0]!=0) {
                fprintf(stderr, "PANIC: delta(): exponential expects argument with no constant term\n");
                abort();
            }
            for (int n=0; n<=N; n++) {
                h[n] = a[n];
                d[n] = a[n];
            }
            d[0] = 1;
            for (int j=2; j<=N; j++) {
                for (int n=j; n<=N; n++) {
                    INTEGER m = 0;
                    for (int k=j-1; k<=n; k++) {
                        m = lcm1_INTEGER(m, h[k]*a[n-k]);
                    }
                    h1[n] = m;
                }
                for (int n=j; n<=N; n++) {
                    h[n] = h1[n];
                    d[n] = lcm1_INTEGER(d[n], FACTORIAL[j]*h[n]);
                }
            }
            return; 
            }
        case LOGARITHM: {
            INTEGER a[N+1];
            INTEGER h[N+1];
            INTEGER h1[N+1];
            delta(a, N, ex->arg1);
            if (a[0]!=1) {
                fprintf(stderr, "PANIC: delta(): logarithm expects argument with constant term 1\n");
                abort();
            }
            a[0] = 0;
            for (int n=0; n<=N; n++) {
                h[n] = a[n];
                d[n] = a[n];
            }
            for (int j=2; j<=N; j++) {
                for (int n=j; n<=N; n++) {
                    INTEGER m = 0;
                    for (int k=j-1; k<=n; k++) {
                        m = lcm1_INTEGER(m, h[k]*a[n-k]);
                    }
                    h1[n] = m;
                }
                for (int n=j; n<=N; n++) {
                    h[n] = h1[n];
                    d[n] = lcm1_INTEGER(d[n], j*h[n]);
                }
            }
            return;
            }                          
        default:
            fprintf(stderr, "PANIC: unknown expr type %i\n", ex->type);
            abort();
    }
}


extern int goldberg_denominator[]; /* defined in goldberg.c */


INTEGER common_denominator(int n, expr_t* ex) {
    if (n==0) {
        return 1;
    }
    if (ex==NULL) {
        return FACTORIAL[n]*goldberg_denominator[n];
    }
    INTEGER d[n+1];
    delta(d, n, ex);
    INTEGER h = d[0];
    for (int j=1; j<=n; j++) {
        h = lcm1_INTEGER(h, d[j]);
    }
    return h;
}



