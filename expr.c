#include"bch.h"
#include<stdlib.h>
#include<stdio.h>
#include<inttypes.h>


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


static int gcd(int a, int b) {
    /* computes greatest common divisor of a and b
     * METHOD: Euclid's classical algorithm
     */
    while (b!=0) {
       int t = b; 
       b = a%b; 
       a = t; 
    }
    return a>=0 ? a : -a;
}


INTEGER gcd_INTEGER(INTEGER a, INTEGER b) {
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
        fprintf(stderr, "PANIC: rat(): zero denominator\n");
        abort();
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
    INTEGER d = gcd_INTEGER(p, q);
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
    INTEGER d = gcd_INTEGER(p, q);
    print_INTEGER(p/d);
    printf("/");
    print_INTEGER(q/d);
}


typedef struct dealloc_list_t {
    void *p;
    struct dealloc_list_t *next;
} dealloc_list_t;


static dealloc_list_t *expr_dealloc_list = NULL;
static int free_all_expressions_already_registered = 0;


void free_all_expressions(void) {
    dealloc_list_t *h = expr_dealloc_list;
    while (h) {
        dealloc_list_t *h0 = h;
        h = h->next;
        free(h0->p);
        free(h0);
    }
    expr_dealloc_list = NULL;
}


expr_t* zero_element(void) {
    expr_t *ex = malloc(sizeof(expr_t));

    if (free_all_expressions_already_registered==0) {
        atexit(free_all_expressions);
        free_all_expressions_already_registered = 1;
    }
    dealloc_list_t *h = malloc(sizeof(dealloc_list_t));
    h->p = ex;
    h->next = expr_dealloc_list;
    expr_dealloc_list = h;

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
    if ((arg1==NULL)||(arg2==NULL)) return NULL;
    expr_t *ex = zero_element();
    ex->type = SUM;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    ex->const_term = rat_add(arg1->const_term, arg2->const_term);
    ex->mindeg = minimum(arg1->mindeg, arg2->mindeg);
    return ex;
}

expr_t* difference(expr_t* arg1, expr_t* arg2) {
    if ((arg1==NULL)||(arg2==NULL)) return NULL;
    expr_t *ex = zero_element();
    ex->type = DIFFERENCE;
    ex->arg1 = arg1;
    ex->arg2 = arg2;
    ex->const_term = rat_sub(arg1->const_term, arg2->const_term);
    ex->mindeg = minimum(arg1->mindeg, arg2->mindeg);
    return ex;
}

expr_t* product(expr_t* arg1, expr_t* arg2) {
    if ((arg1==NULL)||(arg2==NULL)) return NULL;
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
    if (arg==NULL) return NULL;
    expr_t *ex = zero_element();
    ex->type = NEGATION;
    ex->arg1 = arg;
    ex->const_term = rat_neg(arg->const_term);
    ex->mindeg = arg->mindeg;
    return ex;
}

expr_t* term_from_rat(rat_t factor, expr_t* arg) {
    if (arg==NULL) return NULL;
    if (factor.den==0) { 
        fprintf(stderr, "ERROR: zero denominator\n");
        return NULL;
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
    if (arg==NULL) return NULL;
    if (arg->const_term.num!=0) {
        fprintf(stderr, "ERROR: exponential expects argument with no constant term\n");
        return NULL;
    }
    expr_t *ex = zero_element();
    ex->type = EXPONENTIAL;
    ex->arg1 = arg;
    ex->const_term = rat(1, 1);
    ex->mindeg = arg->mindeg; 
    return ex;
}

expr_t* logarithm(expr_t* arg) {
    if (arg==NULL) return NULL;
    if (!((arg->const_term.num==1) && (arg->const_term.den==1))) {
        fprintf(stderr, "ERROR: logarithm expects argument with constant term == 1\n");
        return NULL;
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
            fprintf(stderr, "PANIC: unknown expr type %i\n", ex->type);
            abort();
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
                return is_lie_element(ex->arg1) && is_lie_element(ex->arg2);
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

