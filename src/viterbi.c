#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <tgmath.h>
#include <string.h>
#define __DEBUG__
#ifdef __DEBUG__
#define DBGMSG(FORMAT, ARGS) printf(FORMAT, ARGS)
#else
#define DBGMSG(x)
#endif

#define LOGZERO -INFINITY
//#define LOGZERO 99999
#define my_log(x) log(x) // detect_nan_log(x, __FUNCTION__, __LINE__)//log(x)
#define FLOAT_APROX 0.000000000000001

typedef long double my_float;
static const my_float EPSILON = 1.0;

static double detect_nan_log(double x, const char *func, int line) {
    if (x < 0) {
        printf("[%s:%d] nan: %lf\n", func, line, x);
        exit(1);
    }
    return log(x);
}

static unsigned long long ncr(unsigned long long n, unsigned long long r) {
    if(r > n - r) r = n - r; // because C(n, r) == C(n, n - r)
    unsigned long long ans = 1;
    unsigned long long i;

    for(i = 0; i < r; i++) {
        ans *= (n - i);
        ans /= (i + 1);
    }

    return ans;
}

static int fact(int z) {
    int f = 1;
    if (z == 0) {
        return f;
    }
    else {
        for (int i = 1; i <= z; i++) {
            f = f * i;
        }
    }
    return f;
}
/*
   static int ncr(int n, int r) {
   return fact(n) / (fact(r) * fact(n - r));
   }*/

static int cal_m(my_float *o, int n) {
    my_float max = o[0];
    for (int i = 1; i < n; i++) {
        if (o[i] > max) {
            max = o[i];
        }
    }
    //to speed up, much smaller than maxIValue but at least 20 larger than
    //max(outputs) if outputs = true + laplace
    //return max + 1; //max + 10
    //return max + 2; //why this cause 2,2, 3???
    return max + 50;
}

static my_float cal_init_prob(int m, int n, int k) {
    return my_log((my_float)ncr(m - k + n - 1, n - 1)) - my_log((my_float)ncr(m + n, n));
}

static my_float (*pi_k)(int m, int n, int k) = cal_init_prob;

static my_float cal_trans_prop_cache_lvy(int m, int n, int i, int j, int t, my_float **cache) {
    if (i > j) {
        return LOGZERO;
    }
    else if (cache[i][j] > 0.0000001) { //??
        return cache[i][j];
    }
    else if (i == j) {
        //cache[i][j] = my_log((n - t) / (my_float)(m - i + n - t));
        cache[i][j] = my_log((n - t)) - my_log(m - i + n - t);
    }
    else if (t == n - 1) {
        //cache[i][j] = my_log(1 / (my_float)(m - i + 1));
        cache[i][j] = my_log(1) - my_log(m - i + 1);
    }
    else {
        //cache[i][j] = (cal_trans_prop_cache_lvy(m, n, i, j - 1, t, cache) + my_log((m - j + 1) /(my_float)(m - j + n - t)));
        cache[i][j] = cal_trans_prop_cache_lvy(m, n, i, j - 1, t, cache) + my_log(m - j + 1) - my_log(m - j + n - t);
    }
    return cache[i][j];
}

static my_float cal_trans_prob_cache(int m, int n, int i, int j, int t) {
    static my_float **trans_prob_cache;
    static int cached_t;
    if (!trans_prob_cache) {
        trans_prob_cache = malloc(sizeof(my_float*) * (m + 1));
        for (int ii = 0; ii <= m; ii++) {
            trans_prob_cache[ii] = malloc(sizeof(my_float) * (m + 1));
            memset(trans_prob_cache[ii], 0, sizeof(my_float) * (m + 1));
        }
        cached_t = t;
    }

    if (cached_t != t) {
        for (int ii = 0; ii <= m; ii++) {
            free(trans_prob_cache[ii]);
        }
        free(trans_prob_cache);
        trans_prob_cache = malloc(sizeof(my_float*) * (m + 1));
        for (int ii = 0; ii <= m; ii++) {
            trans_prob_cache[ii] = malloc(sizeof(my_float) * (m + 1));
            memset(trans_prob_cache[ii], 0, sizeof(my_float) * (m + 1));
        }
        cached_t = t;
    }
    return cal_trans_prop_cache_lvy(m, n, i, j, t, trans_prob_cache);
}

static my_float cal_trans_prob(int m, int n, int i, int j, int t) {
    //return ncr(m -j + n - t - 1, n - t -  1) / (my_float)ncr(m - i + n - t, n - t);
    if (i > j) {
        return LOGZERO;
    }
    else if (i == j) {
        return my_log((n - t) / (my_float)(m - i + n - t));
    }
    else if (t == n - 1) {
        return my_log(1 / (my_float)(m - i + 1));
    }
    else {
        return cal_trans_prob(m, n, i, j - 1, t) + my_log((m - j + 1) / (my_float)(m - j + n - t));
    }
}

static my_float (*a_i_j_t)(int m, int n, int i, int j, int t) = cal_trans_prob_cache;//cal_trans_prob;//cal_trans_prob_cache;

//use my_log version?
static my_float cal_emission_prob(my_float *o, int t, int j) {
    //return my_log(EPSILON / 2 * exp(-EPSILON * abs(o[t] - j)));
    return my_log(EPSILON) - my_log(2) - EPSILON * abs(o[t] - j);
}

static my_float (*bj)(my_float *o, int t, int j) = cal_emission_prob;

static my_float maximum(my_float a[], my_float n) {
    int c, idx;
    my_float max;

    max = a[0];
    idx = 0;

    for (idx = 1; idx < n; idx++) {
        if (a[idx] > max) {
            max = a[idx];
        }
    }

    return max;
}

static my_float mimimum(my_float a[], my_float n) {
    int c, idx;
    my_float min;

    min = a[0];
    idx = 0;

    for (idx = 1; idx < n; idx++) {
        if (a[idx] < min) {
            min = a[idx];
        }
    }

    return min;
}

static void adjust_to_prevent_underflow(my_float a[], int n) {
    my_float c_t = maximum(a, n);
    //my_float c_t = mimimum(a, n);

    for (int i = 0; i < n; i++) {
        a[i] -= c_t;
    }
}

void viterbi(my_float o[], int q_star[], int n) {
    /* initilization */
    int m = cal_m(o, n);
    my_float **delta; //delta_t(j) for 0 <= j <= m, 0 <= t <= (T-1)=(n-1)
    int **psi;

    delta = malloc(sizeof(my_float*) * n);// malloc memory from heap
    psi = malloc(sizeof(int*) * n);
    for (int i = 0; i < n; i++) {
        delta[i] = malloc(sizeof(my_float) * (m + 1));
        memset(delta[i], 0, sizeof(my_float) * (m + 1));
        psi[i] = malloc(sizeof(int) * (m + 1));
        memset(psi[i], 0, sizeof(int) * (m + 1));
    }

    printf("m = %d\n", m); // m = 4
    // compute t = 0, delta_0(i) = pi_i * bi(o_0)
    // my_log delta_0(i) = my_log pi_i + my_log bi(o_0)
    for (int i = 0; i <= m; i++) {
        delta[0][i] = pi_k(m, n, i) + bj(o, 0, i);
        psi[0][i] = 0;
    }

    // adjust_to_prevent_underflow(delta[0], m + 1);
    // delta[0][i] stores my_log delta[0][i]'

    /* recursion */
    for (int t = 1; t < n; t++) {
        printf("t = %d\n", t);
        for (int j = 0; j <= m; j++) {
            delta[t][j] = delta[t - 1][j] + a_i_j_t(m, n, 0, j, t);
            psi[t][j] = 0;
            for (int i = 0; i <= m; i++) {
                if (delta[t - 1][i] + a_i_j_t(m, n, i, j, t) > delta[t][j]) {
                    delta[t][j] = delta[t - 1][i] + a_i_j_t(m, n, i, j, t);
                    psi[t][j] = i;
                }
            }
            delta[t][j] += bj(o, t, j);

        }
    // adjust_to_prevent_underflow(delta[t], m + 1);
    }

    for (int t = 0; t < n; t++) {
        for (int i = 0; i <= m; i++) {
            printf("delta[%d][%d]=%LF,psi[%d][%d]=%d\n", t, i, delta[t][i], t, i, psi[t][i]);
        }
    }
    /* termination */
    {
        int num = n - 1;
        my_float p_star = delta[num][0];

        for (int i = 1; i <= m; i++) {
            if (delta[num][i] > p_star) {
                p_star = delta[num][i];
                q_star[num] = i;
            }
        }
    }

    /* backstrace*/
    for (int t = n - 2; t >= 0; t--) {
        q_star[t] = psi[t + 1][q_star[t + 1]];
    }
    return;
}

void viterbi_with_log_init_prob(my_float o[], int q_star[], int n, my_float log_init_prop[], int m) {
    /* initilization */
    my_float **delta; //delta_t(j) for 0 <= j <= m, 0 <= t <= (T-1)=(n-1)
    int **psi;

    delta = malloc(sizeof(my_float*) * n);// malloc memory from heap
    psi = malloc(sizeof(int*) * n);
    for (int i = 0; i < n; i++) {
        delta[i] = malloc(sizeof(my_float) * (m + 1));
        memset(delta[i], 0, sizeof(my_float) * (m + 1));
        psi[i] = malloc(sizeof(int) * (m + 1));
        memset(psi[i], 0, sizeof(int) * (m + 1));
    }

    for (int i = 0; i <= m; i++) {
        delta[0][i] = log_init_prop[i] + bj(o, 0, i);
    }

    printf("m = %d\n", m); // m = 4

    // adjust_to_prevent_underflow(delta[0], m + 1);
    // delta[0][i] stores my_log delta[0][i]'

    /* recursion */
    for (int t = 1; t < n; t++) {
        for (int j = 0; j <= m; j++) {
            delta[t][j] = delta[t - 1][j] + a_i_j_t(m, n, 0, j, t);
            psi[t][j] = 0;
            for (int i = 0; i <= m; i++) {
                if (delta[t - 1][i] + a_i_j_t(m, n, i, j, t) - delta[t][j] > FLOAT_APROX) {
                    delta[t][j] = delta[t - 1][i] + a_i_j_t(m, n, i, j, t);
                    psi[t][j] = i;
                }
            }
            delta[t][j] += bj(o, t, j);
        }
//        adjust_to_prevent_underflow(delta[t], m + 1);
    }

    for (int t = 0; t < n; t++) {
        for (int i = 0; i <= m; i++) {
            printf("delta[%d][%d]=%LF,psi[%d][%d]=%d\n", t, i, delta[t][i], t, i, psi[t][i]);
        }
    }
    /* termination */
    {
        int num = n - 1;
        my_float p_star = delta[num][0];

        for (int i = 1; i <= m; i++) {
            if (delta[num][i] > p_star) {
                p_star = delta[num][i];
                q_star[num] = i;
            }
        }

    /* backstrace*/
        for (int t = num - 1; t >= 0; t--) {
            q_star[t] = psi[t + 1][q_star[t + 1]];
        }
    }
    return;
}


