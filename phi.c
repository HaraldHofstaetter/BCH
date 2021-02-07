#include"bch.h"
#include<stdlib.h>
#include<stdio.h>

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



#ifdef USE_INT128_T
int str_INTEGER(char *out, __int128_t x) {
    int pos = 0;
    int s = 1;
    if (x<0) {
        s = -1;
        x = -x;
    }
    uint64_t F = 100000000000000000ULL;
    int64_t x1 = x % F;
    x /=F;
    if (x>0) {
        int64_t x2 = x % F;
        x /= F;
        if (x>0) {
            int64_t x3 = x;
            if (out) {
                pos += sprintf(out+pos, "%li%017li%017li", s*x3,x2,x1);
            }
            else {
                pos += snprintf(NULL, 0, "%li%017li%017li", s*x3,x2,x1);
            }
        }
        else {
            if (out) {
                pos += sprintf(out+pos, "%li%017li", s*x2,x1);
            }
            else {
                pos += snprintf(NULL, 0, "%li%017li", s*x2,x1);
            }
        }
    }
    else {
        if (out) {
            pos += sprintf(out+pos, "%li", s*x1);
        }
        else {
            pos += snprintf(NULL, 0,  "%li", s*x1);
        }
    }
    return pos;
}

void print_INTEGER(__int128_t x) {
    char out[40];  /* log10(2^128) = 38.53... */
    str_INTEGER(out, x);
    printf("%s", out);
}
#else
int str_INTEGER(char *out, int64_t x) {
    if (out) {
        return sprintf(out, "%li",x);
    }
    else {
        return snprintf(NULL, 0, "%li",x);
    }
}

void print_INTEGER(int64_t x) {
    printf("%li",x);
}
#endif


static INTEGER gcd(INTEGER a, INTEGER b) {
    /* computes greatest common divisor of a and b
     * METHOD: Euclid's classical algorithm
     */
    while (b!=0) {
       INTEGER t = b; 
       b = a%b; 
       a = t; 
    }
    return a>=0 ? a : -a;
}

int str_RATIONAL(char *out, INTEGER p, INTEGER q) {
    INTEGER d = gcd(p, q);
    int pos = 0;
    if (out) {
        pos += str_INTEGER(out+pos, p/d);
        pos += sprintf(out+pos, "/");
        pos += str_INTEGER(out+pos, q/d);
    }
    else {
        pos += 1 + str_INTEGER(NULL, p/d) + str_INTEGER(NULL, q/d);
    }
    return pos;
}

void print_RATIONAL(INTEGER p, INTEGER q) {
    INTEGER d = gcd(p, q);
    print_INTEGER(p/d);
    printf("/");
    print_INTEGER(q/d);
}

static expr_t* undefined_expr(void) {
    expr_t *ex = malloc(sizeof(expr_t));
    ex->type = UNDEFINED;
    ex->arg1 = NULL;
    ex->arg2 = NULL;
    ex->num = 0;
    ex->den = 0;
    return ex;
}

expr_t* identity(void) {
    expr_t *ex = undefined_expr();
    ex->type = IDENTITY;
    return ex;
}

expr_t* generator(uint8_t n) {
    expr_t *ex = undefined_expr();
    ex->type = GENERATOR;
    ex->num = n;
    return ex;
}

expr_t* sum(expr_t* arg1, expr_t* arg2) {
    expr_t *ex = undefined_expr();
    ex->type = SUM;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    return ex;
}

expr_t* difference(expr_t* arg1, expr_t* arg2) {
    expr_t *ex = undefined_expr();
    ex->type = DIFFERENCE;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    return ex;
}

expr_t* product(expr_t* arg1, expr_t* arg2) {
    expr_t *ex = undefined_expr();
    ex->type = PRODUCT;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    return ex;
}

expr_t* negation(expr_t* arg) {
    expr_t *ex = undefined_expr();
    ex->type = NEGATION;
    ex->arg1 = arg;
    return ex;
}

expr_t* term(int num, int den, expr_t* arg) {
    expr_t *ex = undefined_expr();
    ex->type = TERM;
    ex->arg1 = arg;
    ex->num = num;
    ex->den = den;
    return ex;
}

expr_t* exponential(expr_t* arg) {
    expr_t *ex = undefined_expr();
    ex->type = EXPONENTIAL;
    ex->arg1 = arg;
    return ex;
}

expr_t* logarithm(expr_t* arg) {
    expr_t *ex = undefined_expr();
    ex->type = LOGARITHM;
    ex->arg1 = arg;
    return ex;
}

expr_t* commutator(expr_t* arg1, expr_t* arg2) {
    return difference(product(arg1, arg2), 
                      product(arg2, arg1));
}


void free_expr(expr_t* ex) {
    if (ex) {
        free(ex->arg1);
        free(ex->arg2);
        free(ex);
    }
}   


int str_expr(char *out, expr_t* ex) {
    int pos = 0;
    switch(ex->type) {
        case IDENTITY:
            if (out) {
                pos += sprintf(out+pos, "Id");
            }
            else {
                pos += 2;
            }
            break;
        case GENERATOR: 
            if (out) {
                pos += sprintf(out+pos, "%c", 'A'+ex->num);
            }
            else {
                pos += 1;
            }
            break;
        case SUM:
            if (out) {
                pos += sprintf(out+pos, "(");
                pos += str_expr(out+pos, ex->arg1);
                pos += sprintf(out+pos, "+");
                pos += str_expr(out+pos, ex->arg2);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 3 + str_expr(NULL, ex->arg1)+ str_expr(NULL, ex->arg2);
            }
            break;
        case DIFFERENCE:
            if (out) {
                pos += sprintf(out+pos, "(");
                pos += str_expr(out+pos, ex->arg1);
                pos += sprintf(out+pos, "-");
                pos += str_expr(out+pos, ex->arg2);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 3 + str_expr(NULL, ex->arg1)+ str_expr(NULL, ex->arg2);
            }
            break;
        case PRODUCT: 
            if (out) {
                pos += str_expr(out+pos, ex->arg1);
                pos += sprintf(out+pos, "*");
                pos += str_expr(out+pos, ex->arg2);
            }
            else {
                pos += 1 + str_expr(NULL, ex->arg1)+ str_expr(NULL, ex->arg2);
            }
            break;
        case NEGATION: 
            if (out) {
                pos += sprintf(out+pos, "(-1)*");
                pos += str_expr(out+pos, ex->arg1);
            }
            else {
                pos += 5 + str_expr(NULL, ex->arg1);
            }
            break;
        case TERM: 
            if (out) {
                pos += sprintf(out+pos, "(%i/%i)*", ex->num, ex->den);
                pos += str_expr(out+pos, ex->arg1);
            }
            else {
                pos += snprintf(NULL, 0, "(%i/%i)*", ex->num, ex->den)
                       + str_expr(NULL, ex->arg1);
            }
            break;
        case EXPONENTIAL:
            if (out) {
                pos += sprintf(out+pos, "exp(");
                pos += str_expr(out+pos, ex->arg1);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 5 + str_expr(NULL, ex->arg1);
            }
            break;
        case LOGARITHM: 
            if (out) {
                pos += sprintf(out+pos, "log(");
                pos += str_expr(out+pos, ex->arg1);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 5 + str_expr(NULL, ex->arg1);
            }
            break;
        default:
            fprintf(stderr, "ERROR: unknown expr type %i\n", ex->type);
            exit(EXIT_FAILURE);
    }
    return pos;
}

void print_expr(expr_t* ex) {
    int n = str_expr(NULL, ex) + 1;
    char *s = malloc(n*sizeof(char));
    str_expr(s, ex);
    printf("%s", s);
    free(s);
}


static inline void check_for_divisibility_by_int(INTEGER p, int q, INTEGER d) {
    if (q*d!=p) {
        int q1 = (q>0?q:-q)/gcd(p,q);
        fprintf(stderr, "ERROR: dividend not divisble by %i\n", q1);
        exit(EXIT_FAILURE);
    }
}

static inline void check_for_divisibility_by_long_int(INTEGER p, long int q, INTEGER d) {
    if (q*d!=p) {
        long int q1 = (q>0?q:-q)/gcd(p,q);
        fprintf(stderr, "ERROR: dividend not divisble by %li\n", q1);
        exit(EXIT_FAILURE);
    }
}

static inline void check_for_divisibility_by_INTEGER(INTEGER p, INTEGER q, INTEGER d) {
    if (q*d!=p) {
        long int q1 = (q>0?q:-q)/gcd(p,q);
        fprintf(stderr, "ERROR: dividend not divisble by %li\n", q1);
        exit(EXIT_FAILURE);
    }
}

int phi(INTEGER y[], int m, uint8_t w[], expr_t* ex, INTEGER v[]) {
    if (m==0) {
        return 0;
    }
    switch (ex->type) {
        case IDENTITY: 
            for (int j=0; j<m; j++) {
                y[j] = v[j];
            }
            return m;
        case GENERATOR: {
            int m1=0;
            for (int j=0; j<m-1; j++) {
                if (w[j]==ex->num) {
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
            int m1 = phi(y, m, w, ex->arg2, v);
            if (m1==0) {
                return 0;
            }
            return phi(y, m1, w, ex->arg1, y);
            }
        case NEGATION: { 
            int m1 = phi(y, m, w, ex->arg1, v);
            for (int j=0; j<m1; j++) {
                y[j] = -y[j];
            }
            return m1;
            }
        case TERM: { 
            int m1 = phi(y, m, w, ex->arg1, v);
            int p = ex->num;
            int q = ex->den;
            for (int j=0; j<m1; j++) {
                INTEGER h = y[j]*p;
                INTEGER d = h/q;
                check_for_divisibility_by_int(h, q, d);
                y[j] = d;
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
                    long int f = FACTORIAL[k]; /* fits into long int => faster execution expected */
                    for (int j=0; j<m1; j++) {
                        INTEGER d = z[j]/f;
                        check_for_divisibility_by_long_int(z[j], f, d);
                        y[j] += d;
                    }
                }
                else {
                    INTEGER f = FACTORIAL[k];
                    for (int j=0; j<m1; j++) {
                        INTEGER d = z[j]/f;
                        check_for_divisibility_by_INTEGER(z[j], f, d);
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
                y[j] = v[j];                    
            } 
            INTEGER h[m];
            int m1 = m; 
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
                if (m3==0) {
                    return m;
                }
                m1 = m3;
                int f = k%2 ? +k : -k; /* f = (-1)^(k+1)*k */ 
                for (int j=0; j<m1; j++) {
                    INTEGER d = z[j]/f;
                    check_for_divisibility_by_int(z[j], f, d);
                    y[j] += d;
                }
            }
            return m;
            }
        default:
            fprintf(stderr, "ERROR: unknown expr type %i\n", ex->type);
            exit(EXIT_FAILURE);
    }
}


static inline INTEGER lcm(INTEGER a, INTEGER b) {
    /* computes least common multiple of a and b */
    return (a/gcd(a,b))*b;
}


static inline INTEGER lcm1(INTEGER a, INTEGER b) {
    if (a==0) {
        return b;
    }
    else if (b==0) {
        return a;
    }
    else {
        return lcm(a, b);
    }
}


static void delta(INTEGER d[], int N, expr_t* ex) {
    switch (ex->type) {
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
        case NEGATION:
            delta(d, N, ex->arg1);
            return;
        case SUM:
        case DIFFERENCE: {
            INTEGER h[N+1];
            delta(d, N, ex->arg1);
            delta(h, N, ex->arg2);
            for (int n=0; n<=N; n++) {
                d[n] = lcm1(d[n], h[n]);
            }
            return;
            }
        case TERM:
            delta(d, N, ex->arg1);
            for (int n=0; n<=N; n++) {
                INTEGER x = ex->den*d[n];
                d[n] = x/gcd(x, ex->num);
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
                    m = lcm1(m, h1[k]*h2[n-k]);
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
                fprintf(stderr, "ERROR: Logarithm expects argument with no constant\n");
                exit(EXIT_FAILURE);
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
                        m = lcm1(m, h[k]*a[n-k]);
                    }
                    h1[n] = m;
                }
                for (int n=j; n<=N; n++) {
                    h[n] = h1[n];
                    d[n] = lcm1(d[n], FACTORIAL[j]*h[n]);
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
                fprintf(stderr, "ERROR: Logarithm expects argument with constant term 1\n");
                exit(EXIT_FAILURE);
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
                        m = lcm1(m, h[k]*a[n-k]);
                    }
                    h1[n] = m;
                }
                for (int n=j; n<=N; n++) {
                    h[n] = h1[n];
                    d[n] = lcm1(d[n], j*h[n]);
                }
            }
            return;
            }                          
        default:
            fprintf(stderr, "ERROR: unknown expr type %i\n", ex->type);
            exit(EXIT_FAILURE);
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
        h = lcm1(h, d[j]);
    }
    return h;
}


