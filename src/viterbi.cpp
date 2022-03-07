#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "viterbi.hpp"
#include <limits>

using namespace std;
using namespace pcyh;

#define LOGZERO -INFINITY
#define my_log(x) log(x) // detect_nan_log(x, __FUNCTION__, __LINE__)//log(x)
#define FLOAT_APPROX 0.000000000000001

typedef long double my_float;
static const my_float EPSILON = 1.0;

extern "C" {
    Viterbi* viterbi_new() { return new Viterbi(); }

    ViterbiResult* viterbi_calculate(Viterbi* viterbi, my_float o[],
            int n, const my_float log_init_prob[], int m, int path_type) {
        return new ViterbiResult(viterbi->calculate(o, n, log_init_prob, m,
                (Viterbi::path_type_e)path_type));
    }

    void viterbi_delete(Viterbi* viterbi) { free(viterbi); }

    void viterbi_set_theta(Viterbi* viterbi, int theta) {
        viterbi->set_theta(theta);
    }

    bool viterbi_result_next(ViterbiResult *result, int q_star[]) {
        return result->next(q_star);
    }
}

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

static int cal_m(const my_float *o, int n) {
    my_float max = o[0];
    for (int i = 1; i < n; i++) {
        if (o[i] > max) {
            max = o[i];
        }
    }
    //to speed up, much smaller than maxIValue but at least 20 larger than
    //max(outputs) if outputs = true + laplace
    return max + 50;
}

static my_float cal_init_prob(int m, int n, int k) {
    return my_log((my_float)ncr(m - k + n - 1, n - 1))
        - my_log((my_float)ncr(m + n, n));
}

static my_float (*pi_k)(int m, int n, int k) = cal_init_prob;

my_float Viterbi::cal_trans_prob_cache_lvy(int m, int n, int i, int j, int t,
        my_float *cache) {
    if (i > j) {
        return LOGZERO;
    }
    else if (cache[j] > 0.0000001) { //??
        return cache[j];
    }
    else if (i == j) {
        cache[j] = my_log((n - t)) - my_log(m - i + n - t);
    }
    else if (t == n - 1) {
        cache[j] = my_log(1) - my_log(m - i + 1);
    }
    else {
        cache[j] = cal_trans_prob_cache_lvy(m, n, i, j - 1, t, cache)
            + my_log(m - j + 1) - my_log(m - j + n - t);
    }
    return cache[j];
}

my_float Viterbi::cal_trans_prob_cache(int m, int n, int i, int j, int t) {
    if (!trans_prob_cache) {
        trans_prob_cache = new my_float[m + 1];
        memset(trans_prob_cache, 0, sizeof(my_float) * (m + 1));
        cached_t = t;
        cached_i = i;
    }

    if (cached_t != t || cached_i != i) {
        memset(trans_prob_cache, 0, sizeof(my_float) * (m + 1));
        cached_t = t;
        cached_i = i;
    }
    return cal_trans_prob_cache_lvy(m, n, i, j, t, trans_prob_cache);
}

static my_float cal_trans_prob(int m, int n, int i, int j, int t) {
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
        return cal_trans_prob(m, n, i, j - 1, t) +
            my_log((m - j + 1) / (my_float)(m - j + n - t));
    }
}

#define a_i_j_t(...) cal_trans_prob_cache(__VA_ARGS__)

//use my_log version?
static my_float cal_emission_prob(const my_float *o, int t, int j) {
    return my_log(EPSILON) - my_log(2) - EPSILON * fabs(o[t] - j);
}

static my_float (*bj)(const my_float *o, int t, int j) = cal_emission_prob;

static my_float maximum(const my_float a[], my_float n) {
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

static my_float mimimum(const my_float a[], my_float n) {
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

    for (int i = 0; i < n; i++) {
        a[i] -= c_t;
    }
}

static int get_start_j(int data, int theta, int min) {
    if (theta < 0) return min;
    return (data - theta) > min ? (data - theta) : min;
}

static int get_stop_j(int data, int theta, int max) {
    if (theta < 0) return max;
    return (data + theta) < max ? (data + theta) : max;
}

#if 0
static void dump_delta_psi(int t, int m, vector<my_float> delta,
        vector< vector<int> > psi) {
    for (int i = 0; i < m + 1; i++) {
        cout << "delta[" << t << "][" << i << "]=" << delta[i];

        cout << ", psi: ";
        for (vector<int>::const_iterator iter = psi[i].begin(); iter != psi[i].end(); ++iter)
            cout << *iter << ' ';
        cout << endl;
    }
}
#endif

ViterbiResult::ViterbiResult() {
}

ViterbiResult::ViterbiResult(int n, int m) {
    psi.resize(n, vector< vector<int> >(m + 1, vector<int>(0)));
    q_star.resize(n);
    stack.push((ViterbiResult::step_t){.t = n, .i = 0, .nth = 0});
}

ViterbiResult::ViterbiResult(const ViterbiResult &result) {
    ViterbiResult::start = result.start; /* max at T */
    ViterbiResult::psi = result.psi;
    int num = ViterbiResult::psi.size();
    q_star.resize(num);
    stack.push((ViterbiResult::step_t){.t = num, .i = 0, .nth = 0});
}

ViterbiResult::ViterbiResult(vector< vector< vector<int> > > psi,
        const vector<int> &start) {
    ViterbiResult::start = start;
    ViterbiResult::psi = psi;
    int num = ViterbiResult::psi.size();
    q_star.resize(num);
    stack.push((ViterbiResult::step_t){.t = num, .i = 0, .nth = 0});
}

bool ViterbiResult::next(int ret[]) {
    int num = psi.size();
    while (!stack.empty()) {
        ViterbiResult::step_t step = stack.top();
        stack.pop();
        if (step.t <= 0) {
            memcpy(ret, &q_star[0], sizeof(int) * num);
            return true;
        }
        else if (step.t >= num) {
            if (step.nth < start.size()) {
                q_star[step.t - 1] = start[step.nth];
                step.nth++;
                stack.push(step);
                step.i = start[step.nth - 1];
                step.nth = 0;
                step.t--;
                stack.push(step);
            }
        }
        else {
            if (step.nth < psi[step.t][step.i].size()) {
                q_star[step.t - 1] = psi[step.t][step.i][step.nth];
                step.nth++;
                stack.push(step);
                step.i = psi[step.t][step.i][step.nth - 1];
                step.nth = 0;
                step.t--;
                stack.push(step);
            }
        }
    }
    return false;
}

void ViterbiResult::reset() {
}

Viterbi::Viterbi() {
    theta = -1;
    trans_prob_cache = NULL;
    path_type = PATH_TYPE_ALL;
    cached_t = 0;
}

Viterbi::~Viterbi() {
    delete_cache();
}

void Viterbi::delete_cache() {
    if (trans_prob_cache) {
        delete[] trans_prob_cache;
        trans_prob_cache = NULL;
    }
}

ViterbiResult Viterbi::calculate( my_float o[], int n,
        const my_float log_init_prob[], int m, path_type_e path_type) {
#define SWICH_DELTA_IDX(x) ((x) ^= 1)
    this->n = n;
    this->m = m;
    this->path_type = path_type;
    int delta_len = theta >= 0 ? theta * 2 + 1 : m + 1;
    /* initilization */
    vector< vector<my_float> > delta(2, vector<my_float>(delta_len, -INFINITY));
    ViterbiResult viterbi_result(n, m);
    vector< vector< vector<int> > > &psi = viterbi_result.psi;
    this->o = o;
    int cnt_delta_buf_idx = 0;

    for (int i = get_start_j(o[0], theta, 0); i <= get_stop_j(o[0], theta, m); i++) {
        delta[cnt_delta_buf_idx][i - get_start_j(o[0], theta, 0)] =
            log_init_prob[i] + bj(o, 0, i);
    }

#if 0
    dump_delta_psi(0, m, delta[0], psi[0]);
#endif

	//adjust_to_prevent_underflow(delta[cnt_delta_buf_idx], m + 1);
    SWICH_DELTA_IDX(cnt_delta_buf_idx);

    /* recursion */
    for (int t = 1; t < n; t++) {
#ifdef __DEBUG__
        cout << "Calculating t=" << t << endl;
#endif
        int start_j = get_start_j(o[t], theta, 0);
        int stop_j = get_stop_j(o[t], theta, m);
        for (int j = start_j; j <= stop_j; j++) {
            int start_i = get_start_j(o[t - 1], theta, 0);
            int stop_i = get_stop_j(o[t - 1], theta, m);
            delta[cnt_delta_buf_idx][j - start_j] = -INFINITY;
            for (int i = start_i; i <= stop_i; i++) {
                if (delta[cnt_delta_buf_idx ^ 1][i - start_i]
                        + a_i_j_t(m, n, i, j, t) -
                        delta[cnt_delta_buf_idx][j - start_j] > FLOAT_APPROX) {
                    /* if A larger than B */
                    delta[cnt_delta_buf_idx][j - start_j] =
                        delta[cnt_delta_buf_idx ^ 1][i - start_i] +
                        a_i_j_t(m, n, i, j, t);
                    psi[t][j].clear();
                    psi[t][j].push_back(i);
                }
                else if (path_type != PATH_TYPE_FIRST
                        && fabs(delta[cnt_delta_buf_idx ^ 1][i - start_i]
                            + a_i_j_t(m, n, i, j, t) -
                            delta[cnt_delta_buf_idx][j - start_j])
                        < FLOAT_APPROX) {
                     /* if A equal B */
                    if (path_type == PATH_TYPE_LAST) {
                        psi[t][j].clear();
                    }
                    psi[t][j].push_back(i);
                }
            }
            delta[cnt_delta_buf_idx][j - start_j] += bj(o, t, j);
        }
#if 0
        dump_delta_psi(t, m, delta[cnt_delta_buf_idx], psi[t]);
#endif
        for (int p = 0; p < delta_len; p++) {
            delta[cnt_delta_buf_idx ^ 1][p] = -std::numeric_limits<my_float>::infinity();
        }
        //adjust_to_prevent_underflow(delta[cnt_delta_buf_idx], m + 1);
        SWICH_DELTA_IDX(cnt_delta_buf_idx);
    }

    delete_cache();

    /* termination */
    {
        int num = n;
        vector<int> &max = viterbi_result.start;
        my_float p_star = -std::numeric_limits<my_float>::infinity();

        int start_i = get_start_j(o[num - 1], theta, 0);
        int stop_i = get_stop_j(o[num - 1], theta, m);
        for (int i = start_i; i <= stop_i; i++) {
            if (p_star == -INFINITY || delta[cnt_delta_buf_idx ^ 1][i - start_i]
                    - p_star > FLOAT_APPROX) {
                p_star = delta[cnt_delta_buf_idx ^ 1][i - start_i];
                max.clear();
                max.push_back(i);
            }
            else if (path_type != PATH_TYPE_FIRST
                    && fabs(delta[cnt_delta_buf_idx ^ 1][i - start_i] - p_star)
                    < FLOAT_APPROX) {
                if (path_type == PATH_TYPE_LAST) {
                    max.clear();
                }
                max.push_back(i);
            }
        }

        return viterbi_result;
    }

}

void Viterbi::set_theta(int theta) {
    Viterbi::theta = theta;
}
