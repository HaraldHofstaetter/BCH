#include"bch.h"
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

//#include<stdint.h>
//typedef uint8_t generator_t;


static void invperm(int n, generator_t a[], generator_t b[]) {
    for (int j=0; j<n; j++) {
        b[a[j]-1] = j+1;
    }
}

static void reverse(int n, generator_t a[], generator_t b[]) {
    for (int j=0; j<n; j++) {
        b[j] = a[n-j-1];
    }
}    

static void copy(int n, generator_t a[], generator_t b[]) {
    for (int j=0; j<n; j++) {
        b[j] = a[j];
    }
}    

static int lexcmp(generator_t a[], int la, generator_t b[], int lb) {
    int l = la<lb ? la : lb;
    int j = 0;
    while ((j<l) && (a[j]==b[j])) {
        j++;
    }
    if (j==l) {
        if (la==lb) {
            return 0;
        }
        else {
            return la < lb ? -1 : +1;
        }
    }
    else {
        return a[j]<b[j] ?  -1 : +1;
    }
}

/*
static void print_word(generator_t w[], int lw) {
    printf("[");
    for (int j=0; j<lw; j++) {
        printf("%i", w[j]);
        if (j<lw-1) {
            printf(", ");
        }
    }
    printf("]");
}
*/

static generator_t digits[] =  {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};


typedef struct list_t {
    generator_t* w;
    int lw;
    int value;
    struct list_t *next;
} list_t;


list_t* init_list(int q) {
    list_t *L = NULL;
    for (int j=q; j>=1; j--) {
        list_t *x = malloc(sizeof(list_t));
        x->w = digits+j;
        x->lw = 1;
        x->value = j;
        x->next = L;
        L = x;
    }
    return L;
}

/*
static void print_list(list_t *x) {
    while (x!=NULL) {
        print_word(x->w, x->lw);
        printf(": %i; ", x->value);
        x = x->next;
    }
}
*/

static int list2arrays(list_t *x, generator_t **w, int lw[], generator_t v[]) {
    int i=0;
    while (x!=NULL) {
        w[i] = x->w;
        lw[i] = x->lw;
        v[i] = x->value;
        i++;
        list_t *x0 = x;
        x = x->next;
        free(x0);
    }
    return i;
}

static int get_list(list_t **L, generator_t w[], int lw, int value) {
    list_t *h = *L;
    list_t *h0 = NULL;
    int c; 
    while ((h!=NULL) && ((c=lexcmp(h->w, h->lw, w, lw))<0)) {
        h0 = h;
        h = h->next;
    }
    if ((h!=NULL) && (c==0)) {
        return h->value;
    }
    else {
        list_t *x = malloc(sizeof(list_t));
        x->w = w;
        x->lw = lw;
        x->value = value;
        x->next = h;
        if (h0==NULL) {
            *L = x;
        }
        else {
            h0->next = x;
        }
        return value;
    }
}


static void analyze_lyndon_word(generator_t w[], int lw, generator_t w2[], int *lw2, generator_t **tt, int *ltt, int *n ) {
    generator_t a = w[0];
    generator_t q = w[0];
    for (int i=1; i<lw; i++) {
        if (w[i]<a) {
            a = w[i];
        }
        if (w[i]>q) {
            q = w[i];
        }
    }
    list_t *L = init_list(q);

    generator_t w1[lw];
    int lw1 = 0;
    int s = 1;
    int m1 = 1;
    int m2 = 0;

    //assert(w[s-1]==a);
    s++;
    while (s<=lw && w[s-1]!=a) {
        s++;
    }
    generator_t* v = w+1;
    int lv = s-2;
    int lav = s-1;
    while (s<=lw) {
        if (m2!=0) {
            m1 = s;
        }
        while ((s+lav<=lw) && (w[s+lav-1]==a) && 
               (w[s-1]==a) && lexcmp(w+s, lv, v, lv)==0) {
            s += lav;
        }
        //assert(w[s-1]==a);
        s++;

        generator_t* uu = w+s-1;
        int luu=0;

        while ((s<lw) && (w[s-1]!=a)) { 
            luu++;
            s++;
        }
        if ((s==lw) && (w[s-1]!=a)) {
            luu++;
            s++;
        }
        int j=1;
        while (!((lexcmp(v, lv, uu, j)<0)&&(lexcmp(v, lv, uu, j-1)>=0))) {
            j++;
        }

        generator_t* u1 = uu+j;
        int lu1 = luu-j;
        m2 = s-lu1-1;
        int x = get_list(&L, w+m1-1, m2-m1+1, q+1);
        if (x==q+1) {
            q++;
        }

        w1[lw1] = x; /* w1 = concat(w1,x,u1) */
        lw1++;
        for(int j=0; j<lu1; j++) {
            w1[lw1] = u1[j];
            lw1++;
        }
    }
    generator_t vv[lw];
    int lvv = list2arrays(L, tt, ltt, vv);
    generator_t pp[lw];
    invperm(lvv, vv, pp);
    for (int j=0; j<lw1; j++) {
        w2[j] = pp[w1[j]-1];
    }
    *lw2 = lw1;
    *n = lvv;
}

void lyndon2rightnormed(int lw, generator_t w[], generator_t r[]) {
    generator_t aa = w[0]; // aa = minimum(w)
    int k = 1; // number of occurences of aa in w
    for (int i=1; i<lw; i++) {
        if (w[i]<aa) {
            aa = w[i];
            k = 1;
        } 
        else if (w[i]==aa) {
            k++;
        }
    }
    if (k==1) {
        reverse(lw, w, r);
        return;
    }

    generator_t w1[lw];
    int lw1;
    generator_t* tt[2*lw]; //TODO: Check that this is large enough...
    int ltt[2*lw];
    int m;
    analyze_lyndon_word(w, lw, w1, &lw1, tt, ltt, &m);

    generator_t r1[lw1];
    lyndon2rightnormed(lw1, w1, r1);

    generator_t* y = tt[r1[lw1-1]-1];
    int ly = ltt[r1[lw1-1]-1];
    generator_t a = y[0];
    
    int k0 = 1; // index of first a in y[1:end]
    for( ; y[k0]!=a; k0++) {} 
    int k1 = ly-1; // index of last a in y
    for( ; y[k1]!=a; k1--) {} 

    generator_t *v = y+1;
    int lv = k0-1;

    generator_t *avn = y+k0;
    int lavn = k1-k0;

    generator_t *u1 = y+k1+1;
    int lu1 = ly-k1-1;

    int lr = 0;
    for (int i=0; i<lw1-1; i++) {
        copy(ltt[r1[i]-1], tt[r1[i]-1], r+lr);
        lr += ltt[r1[i]-1];
    }
    copy(lavn, avn, r+lr);
    lr += lavn;
    r[lr] = a;
    lr++;
    copy(lu1, u1, r+lr);
    lr += lu1;
    reverse(lv, v, r+lr);
    lr += lv;
    r[lr] = a;
    lr++;
    //assert(lw==lr);
}


/*
int main(void) {
    generator_t w[] = {0, 0, 1, 2, 0, 2, 0, 1, 1};
    int n = 9;
    generator_t w1[n];
    int lw1;
    generator_t* tt[n];
    int ltt[n];
    int m;
    
    analyze_lyndon_word(w, n, w1, &lw1, tt, ltt, &m);
    printf("w1="); print_word(w1, lw1); printf("\n");
    printf("tt=["); 
    for (int j=0; j<m; j++) {
        print_word(tt[j], ltt[j]);
        if (j<m-1) {
            printf(", ");
        }
    }
    printf("]\n");

    generator_t r[n];

    lyndon2rightnormed(n, w, r);
    printf("r="); print_word(r, n); printf("\n");
}
*/

