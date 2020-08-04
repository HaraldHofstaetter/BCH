#include"bch.h"
#include<stdlib.h>
#include<stdio.h>


static INTEGER *FACTORIAL=NULL;

/* table den_fac obtained with the following Julia code:
n = 33
F = [factorial(Int128(k)) for k=0:n-1]
M = zeros(Int128,n,n)
M[:,1] = F
for m = 2:n
    M[m+1:end,m] = [lcm([F[k]*M[n-k+1,m-1] for k=2:n-m+1]) for n=m+1:n]  
end
using LinearAlgebra # for diagm
M *= diagm(1:n)
D = [lcm(M[k,1:k-1]) for k=1:n]
den_fac = [div(D[i],F[i]) for i=1:n]
 */

static int den_fac[33] = {1, 1, 1, 2, 1, 6, 2, 6, 3, 10, 2, 6, 2, 210, 30, 12, 3, 30, 10, 
                          210, 42, 330, 30, 60, 30, 546, 42, 28, 2, 60, 4, 924, 231};


void init_factorial(int n) {
    FACTORIAL = malloc((n+1)*sizeof(INTEGER)); 
    FACTORIAL[0] = 1;
    for (int k=1; k<=n; k++) {
        FACTORIAL[k] = k*FACTORIAL[k-1];
    }
}

void free_factorial(void) {
    free(FACTORIAL);
}

INTEGER common_denominator(int n) {
    return FACTORIAL[n]*den_fac[n];
}

#ifdef USE_INT128_T
void print_INTEGER(__int128_t x) {
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
            printf("%li%017li%017li",s*x3,x2,x1);
        }
        else {
            printf("%li%017li",s*x2,x1);
        }
    }
    else {
        printf("%li",s*x1);
    }
}
#else
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

expr_t* generator(generator_t n) {
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

void print_expr(expr_t* ex) {
    switch(ex->type) {
        case IDENTITY:
            printf("Id");
            break;
        case GENERATOR: 
            printf("%c", 'A'+ex->num);
            break;
        case SUM:
            printf("(");
            print_expr(ex->arg1);
            printf("+");
            print_expr(ex->arg2);
            printf(")");
            break;
        case DIFFERENCE:
            printf("(");
            print_expr(ex->arg1);
            printf("-");
            print_expr(ex->arg2);
            printf(")");
            break;
        case PRODUCT: 
            print_expr(ex->arg1);
            printf("*");
            print_expr(ex->arg2);
            break;
        case NEGATION: 
            printf("(-1)*");
            print_expr(ex->arg1);
            break;
        case TERM: 
            printf("(%i/%i)*", ex->num, ex->den);
            print_expr(ex->arg1);
            break;
        case EXPONENTIAL:
            printf("exp(");
            print_expr(ex->arg1);
            printf(")");
            break;
        case LOGARITHM: 
            printf("log(");
            print_expr(ex->arg1);
            printf(")");
            break;
        default:
            fprintf(stderr, "ERROR: unknown expr type %i\n", ex->type);
            exit(EXIT_FAILURE);
    }
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

int phi(INTEGER y[], int m, generator_t w[], expr_t* ex, INTEGER v[]) {
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
            int m1 = phi(y, m, w, ex->arg1, v);
            int m2 = phi(y2, m, w, ex->arg2, v);
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
            int m1 = phi(y, m, w, ex->arg1, v);
            int m2 = phi(y2, m, w, ex->arg2, v);
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
