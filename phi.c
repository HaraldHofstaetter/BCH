#include"bch.h"
#include<stdlib.h>
#include<stdio.h>
#include <inttypes.h>

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


static const int HUGE_NUMBER = 1<<30;


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
                pos += sprintf(out+pos, "%" PRId64 "%017" PRId64 "%017" PRId64, s*x3,x2,x1);
            }
            else {
                pos += snprintf(NULL, 0, "%" PRId64 "%017" PRId64 "%017" PRId64, s*x3,x2,x1);
            }
        }
        else {
            if (out) {
                pos += sprintf(out+pos, "%" PRId64 "%017" PRId64, s*x2,x1);
            }
            else {
                pos += snprintf(NULL, 0, "%" PRId64 "%017" PRId64, s*x2,x1);
            }
        }
    }
    else {
        if (out) {
            pos += sprintf(out+pos, "%" PRId64, s*x1);
        }
        else {
            pos += snprintf(NULL, 0,  "%" PRId64, s*x1);
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


static inline int minimum(int a, int b) {
    return a<=b ? a : b;
}

static inline int maximum(int a, int b) {
    return a>=b ? a : b;
}

rat_t rat(int num, int den) {
    if (den==0) {
        fprintf(stderr, "ERROR: rat(): zero denominator\n");
        exit(EXIT_FAILURE);
    }
    rat_t r;
    if (den<0) {
        num = -num;
        den = -den;
    }
    int d = gcd(num, den); 
    r.num = num/d;
    r.den = den/d;
    return r;
}

rat_t rat_add(rat_t a, rat_t b) {
    return rat(a.num*b.den+b.num*a.den, a.den*b.den);
}
rat_t rat_sub(rat_t a, rat_t b) {
    return rat(a.num*b.den-b.num*a.den, a.den*b.den);
}
rat_t rat_mul(rat_t a, rat_t b) {
    return rat(a.num*b.num, a.den*b.den);
}
rat_t rat_div(rat_t a, rat_t b) {
    return rat(a.den*b.den, a.num*b.num);
}
rat_t rat_neg(rat_t a) {
    return rat(-a.num, a.den);
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

expr_t* zero_element(void) {
    expr_t *ex = malloc(sizeof(expr_t));
    ex->type = ZERO_ELEMENT;
    ex->arg1 = NULL;
    ex->arg2 = NULL;
    ex->factor = rat(0,1);
    ex->gen = -1;
    ex->const_term = rat(0,1);
    ex->mindeg = HUGE_NUMBER;
    return ex;
}

expr_t* identity(void) {
    expr_t *ex = zero_element();
    ex->type = IDENTITY;
    ex->const_term = rat(1,1);
    ex->mindeg = HUGE_NUMBER;
    return ex;
}

expr_t* generator(uint8_t n) {
    expr_t *ex = zero_element();
    ex->type = GENERATOR;
    ex->gen = n;
    ex->const_term = rat(0,1);
    ex->mindeg = 1;
    return ex;
}

expr_t* sum(expr_t* arg1, expr_t* arg2) {
    expr_t *ex = zero_element();
    ex->type = SUM;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    ex->const_term = rat_add(arg1->const_term, arg2->const_term);
    ex->mindeg = minimum(arg1->mindeg, arg2->mindeg);
    return ex;
}

expr_t* difference(expr_t* arg1, expr_t* arg2) {
    expr_t *ex = zero_element();
    ex->type = DIFFERENCE;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    ex->const_term = rat_sub(arg1->const_term, arg2->const_term);
    ex->mindeg = minimum(arg1->mindeg, arg2->mindeg);
    return ex;
}

expr_t* product(expr_t* arg1, expr_t* arg2) {
    expr_t *ex = zero_element();
    ex->type = PRODUCT;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    ex->const_term = rat_mul(arg1->const_term, arg2->const_term);
    if ((arg1->const_term.num!=0)&&(arg2->const_term.num!=0)) {
        ex->mindeg = minimum(arg1->mindeg, arg2->mindeg);
    }
    else if ((arg1->const_term.num==0)&&(arg2->const_term.num!=0)) {
        ex->mindeg = arg2->mindeg;
    }
    else if ((arg2->const_term.num==0)&&(arg1->const_term.num!=0)) {
        ex->mindeg = arg1->mindeg;
    }
    else {
        ex->mindeg = minimum(HUGE_NUMBER, arg1->mindeg + arg2->mindeg);
    }
    return ex;
}

expr_t* negation(expr_t* arg) {
    expr_t *ex = zero_element();
    ex->type = NEGATION;
    ex->arg1 = arg;
    ex->const_term = rat_neg(arg->const_term);
    ex->mindeg = arg->mindeg;
    return ex;
}

expr_t* term_from_rat(rat_t factor, expr_t* arg) {
    if (factor.den==0) { 
        fprintf(stderr, "ERROR: zero denominator\n");
        exit(EXIT_FAILURE);
    }
    expr_t *ex = zero_element();
    ex->type = TERM;
    ex->arg1 = arg;
    ex->factor = factor;
    ex->const_term = rat_mul(factor, arg->const_term);
    if (factor.num==0) {
         ex->mindeg = HUGE_NUMBER;
    }
    else {
         ex->mindeg = arg->mindeg;
    }
    return ex;
}

expr_t* term(int num, int den, expr_t* arg) {
    rat_t r;
    r.num = num;
    r.den = den;
    return term_from_rat(r, arg);
}

expr_t* exponential(expr_t* arg) {
    if (arg->const_term.num!=0) {
        fprintf(stderr, "ERROR: exponential expects argument with no constant term\n");
        exit(EXIT_FAILURE);
    }
    expr_t *ex = zero_element();
    ex->type = EXPONENTIAL;
    ex->arg1 = arg;
    ex->const_term = rat(1, 1);
    ex->mindeg = arg->mindeg; 
    return ex;
}

expr_t* logarithm(expr_t* arg) {
    if (!((arg->const_term.num==1) && (arg->const_term.den==1))) {
        fprintf(stderr, "ERROR: logarithm expects argument with constant term == 1\n");
        exit(EXIT_FAILURE);
    }
    expr_t *ex = zero_element();
    ex->type = LOGARITHM;
    ex->arg1 = arg;
    ex->const_term = rat(0, 1);
    ex->mindeg = arg->mindeg; 
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

int str_expr(char *out, expr_t* ex, char* gens) {
    int pos = 0;
    if (gens == 0) {
        gens = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    }
    switch(ex->type) {
        case ZERO_ELEMENT:
            if (out) {
                pos += sprintf(out+pos, "0");
            }
            else {
                pos += 1;
            }
            break;
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
                pos += sprintf(out+pos, "%c", gens[ex->gen]);
            }
            else {
                pos += 1;
            }
            break;
        case SUM:
            if (out) {
                pos += sprintf(out+pos, "(");
                pos += str_expr(out+pos, ex->arg1, gens);
                pos += sprintf(out+pos, "+");
                pos += str_expr(out+pos, ex->arg2, gens);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 3 + str_expr(NULL, ex->arg1, gens)+str_expr(NULL, ex->arg2, gens);
            }
            break;
        case DIFFERENCE:
            if ((ex->arg1->type==PRODUCT)&&(ex->arg1->type==PRODUCT)
                    &&(ex->arg1->arg1==ex->arg2->arg2)
                    &&(ex->arg1->arg2==ex->arg2->arg1)) { /* Commutator */
                if (out) {
                    pos += sprintf(out+pos, "[");
                    pos += str_expr(out+pos, ex->arg1->arg1, gens);
                    pos += sprintf(out+pos, ",");
                    pos += str_expr(out+pos, ex->arg1->arg2, gens);
                    pos += sprintf(out+pos, "]");
                }
                else {
                    pos += 3 + str_expr(NULL, ex->arg1->arg1, gens)+str_expr(NULL, ex->arg1->arg2, gens);
                }
            }
            else {
                if (out) {
                    pos += sprintf(out+pos, "(");
                    pos += str_expr(out+pos, ex->arg1, gens);
                    pos += sprintf(out+pos, "-");
                    pos += str_expr(out+pos, ex->arg2, gens);
                    pos += sprintf(out+pos, ")");
                }
                else {
                    pos += 3 + str_expr(NULL, ex->arg1, gens)+str_expr(NULL, ex->arg2, gens);
                }
            }
            break;
        case PRODUCT: 
            if (out) {
                pos += str_expr(out+pos, ex->arg1, gens);
                pos += sprintf(out+pos, "*");
                pos += str_expr(out+pos, ex->arg2, gens);
            }
            else {
                pos += 1 + str_expr(NULL, ex->arg1, gens)+str_expr(NULL, ex->arg2, gens);
            }
            break;
        case NEGATION: 
            if (out) {
                pos += sprintf(out+pos, "(-1)*");
                pos += str_expr(out+pos, ex->arg1, gens);
            }
            else {
                pos += 5 + str_expr(NULL, ex->arg1, gens);
            }
            break;
        case TERM: 
            if (out) {
                pos += sprintf(out+pos, "(%i/%i)*", ex->factor.num, ex->factor.den);
                pos += str_expr(out+pos, ex->arg1, gens);
            }
            else {
                pos += snprintf(NULL, 0, "(%i/%i)*", ex->factor.num, ex->factor.den)
                       + str_expr(NULL, ex->arg1, gens);
            }
            break;
        case EXPONENTIAL:
            if (out) {
                pos += sprintf(out+pos, "exp(");
                pos += str_expr(out+pos, ex->arg1, gens);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 5 + str_expr(NULL, ex->arg1, gens);
            }
            break;
        case LOGARITHM: 
            if (out) {
                pos += sprintf(out+pos, "log(");
                pos += str_expr(out+pos, ex->arg1, gens);
                pos += sprintf(out+pos, ")");
            }
            else {
                pos += 5 + str_expr(NULL, ex->arg1, gens);
            }
            break;
        default:
            fprintf(stderr, "ERROR: unknown expr type %i\n", ex->type);
            exit(EXIT_FAILURE);
    }
    return pos;
}

void print_expr(expr_t* ex, char* gens) {
    int n = str_expr(NULL, ex, gens) + 1;
    char *s = malloc(n*sizeof(char));
    str_expr(s, ex, gens);
    printf("%s", s);
    free(s);
}



static inline void check_for_divisibility_by_int(INTEGER p, int q, INTEGER d, char *s) {
    if (q*d!=p) {
        int q1 = (q>0?q:-q)/gcd(p,q);
        fprintf(stderr, "ERROR: dividend not divisble by %i %s\n", q1, s);
        exit(EXIT_FAILURE);
    }
}

static inline void check_for_divisibility_by_long_int(INTEGER p, long int q, INTEGER d, char *s) {
    if (q*d!=p) {
        long int q1 = (q>0?q:-q)/gcd(p,q);
        fprintf(stderr, "ERROR: dividend not divisble by %li %s\n", q1, s);
        exit(EXIT_FAILURE);
    }
}

static inline void check_for_divisibility_by_INTEGER(INTEGER p, INTEGER q, INTEGER d, char *s) {
    if (q*d!=p) {
        long int q1 = (q>0?q:-q)/gcd(p,q);
        fprintf(stderr, "ERROR: dividend not divisble by %li %s\n", q1, s);
        exit(EXIT_FAILURE);
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
                    long int f = FACTORIAL[k]; /* fits into long int => faster execution expected */
                    for (int j=0; j<m1; j++) {
                        INTEGER d = z[j]/f;
                        check_for_divisibility_by_long_int(z[j], f, d, "in phi()/EXPONENTIAL");
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
                y[j] = v[j];                    
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
                d[n] = lcm1(d[n], h[n]);
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
                d[n] = x/gcd(x, ex->factor.num);
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
                fprintf(stderr, "ERROR: Exponential expects argument with no constant term\n");
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


static int is_product_of_exponentials_of_lie_elements(expr_t *ex) {
    return ((ex->type==EXPONENTIAL)&&is_lie_element(ex->arg1))
        || ((ex->type==PRODUCT)&&is_product_of_exponentials_of_lie_elements(ex->arg1)
                               &&is_product_of_exponentials_of_lie_elements(ex->arg1));
}


int is_lie_element(expr_t* ex) {
    switch (ex->type) {
        case ZERO_ELEMENT:
        case GENERATOR:
                return 1;
        case SUM:
        case DIFFERENCE: /* Commutator */
                if ((ex->arg1->type==PRODUCT)&&(ex->arg1->type==PRODUCT)
                    &&(ex->arg1->arg1==ex->arg2->arg2)
                    &&(ex->arg1->arg2==ex->arg2->arg1)
                    &&is_lie_element(ex->arg1->arg1)
                    &&is_lie_element(ex->arg1->arg2))
                    return 1;
                else {
                    return is_lie_element(ex->arg1) && is_lie_element(ex->arg2);
                }
        case LOGARITHM:
                return is_product_of_exponentials_of_lie_elements(ex->arg1);
        case NEGATION:
        case TERM:
                return is_lie_element(ex->arg1);
        case IDENTITY:
        case PRODUCT:
        case EXPONENTIAL:
        default:
            return 0;
    }
}

