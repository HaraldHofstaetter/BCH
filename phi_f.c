#include"bch.h"
#include <stdlib.h>

static FLOAT RECIP_FACTORIAL[33];
static FLOAT RECIP[33];


void init_phi_f(int n) {
    RECIP_FACTORIAL[0] = one_f();
    for (int k=1; k<=n; k++) {
        RECIP_FACTORIAL[k] = div_f(RECIP_FACTORIAL[k-1], i2f(k));
    }
    RECIP[0] = zero_f(); /* should be undefined (1/0) */
    for (int k=1; k<=n; k++) {
        RECIP[k] = div_f(one_f(), i2f(k));
    }

}


int phi_f(FLOAT y[], int m, uint8_t w[], expr_t* ex, FLOAT v[]) {
    if (m==0) {
        return 0;;
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
                    if (!is_zero_f(y[j])) {
                        m1 = j+1;
                    }
                }
                else {
                    y[j] = zero_f();
                }
            }
            return m1;
            }
        case SUM: { 
            FLOAT y2[m];
            for (int j=0; j<m; j++) {
                y2[j] = v[j];
            }
            int m1 = phi_f(y, m, w, ex->arg1, v);
            int m2 = phi_f(y2, m, w, ex->arg2, y2);
            if (m1<m2) {
                for (int j=0; j<m1; j++) {
                    y[j] = add_f(y[j], y2[j]);
                }
                for (int j=m1; j<m2; j++) {
                    y[j] = y2[j];
                }
                return m2;
            }
            else {
                for (int j=0; j<m2; j++) {
                    y[j] = add_f(y[j], y2[j]);
                }
                return m1;
            }
            } 
        case DIFFERENCE: {
            FLOAT y2[m];
            for (int j=0; j<m; j++) {
                y2[j] = v[j];
            }
            int m1 = phi_f(y, m, w, ex->arg1, v);
            int m2 = phi_f(y2, m, w, ex->arg2, y2);
            if (m1<m2) {
                for (int j=0; j<m1; j++) {
                    y[j] = sub_f(y[j], y2[j]);
                }
                for (int j=m1; j<m2; j++) {
                    y[j] = neg_f(y2[j]);
                }
                return m2;
            }
            else {
                for (int j=0; j<m2; j++) {
                    y[j] = sub_f(y[j], y2[j]);
                }
                return m1;
            }
            } 
        case PRODUCT: {
            int  md = ex->arg1->mindeg;
            if (md>=m) {
                return 0;
            }
            int m1 = phi_f(y+md, m-md, w+md, ex->arg2, v+md);
            if (m1==0) {
                return 0;
            }
            for (int j=0; j<md; j++) {
                y[j] = zero_f();
            }
            return phi_f(y, m1+md, w, ex->arg1, y);
            }
        case NEGATION: { 
            int m1 = phi_f(y, m, w, ex->arg1, v);
            for (int j=0; j<m1; j++) {
                y[j] = neg_f(y[j]);
            }
            return m1;
            }
        case TERM_F: 
        case TERM: { 
            int m1 = phi_f(y, m, w, ex->arg1, v);
            for (int j=0; j<m1; j++) {
                y[j] = mul_f(ex->factor_f, y[j]);
            }
            return m1;
            }
        case EXPONENTIAL: {
            FLOAT z[m];
            for (int j=0; j<m; j++) {
                z[j] = v[j];
                y[j] = v[j];
            }
            int m1 = m;
            for (int k=1; k<m; k++) {
                m1 = phi_f(z, m1, w, ex->arg1, z);
                if (m1==0) {
                    return m;
                }
                FLOAT rf = RECIP_FACTORIAL[k];
                for (int j=0; j<m1; j++) {
                    y[j] = add_f(y[j], mul_f(z[j], rf));
                }
            }
            return m;
            }
        case LOGARITHM: {
            FLOAT z[m];
            for (int j=0; j<m; j++) {
                z[j] = v[j];
                y[j] = zero_f();                     
            } 
            FLOAT h[m];
            int m1 = m; 
            int mret = m;
            for (int k=1; k<m; k++) {
                for (int j=0; j<m1; j++) {
                    h[j] = z[j];
                }
                int m2 = phi_f(z, m1, w, ex->arg1, z);
                int m3 = 0;
                for (int j=0; j<m2; j++) {
                    z[j] = sub_f(z[j], h[j]);
                    if (!is_zero_f(z[j])) {
                        m3 = j+1;
                    }
                }
                for (int j=m2; j<m1; j++) {
                    z[j] = neg_f(h[j]);
                    if (!is_zero_f(z[j])) {
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
                FLOAT sr = k%2 ? RECIP[k] : neg_f(RECIP[k]);
                for (int j=0; j<m1; j++) {
                    y[j] = add_f(y[j], mul_f(z[j], sr));
                }
            }
            return mret;
            }
        default:
            fprintf(stderr, "PANIC: unknown expr type %i\n", ex->type);
            abort();
    }
}

