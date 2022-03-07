#ifndef __VITERBI_HPP__
#define __VITERBI_HPP__

#include <vector>
#include <stack>

namespace pcyh {

typedef long double my_float;

class ViterbiResult {
 public:
    ViterbiResult();

    ViterbiResult(int n, int m);

    ViterbiResult(const ViterbiResult &result);

    ViterbiResult(std::vector< std::vector< std::vector<int> > > psi,
            const std::vector<int> &start);

    bool next(int ret[]);

    void reset();

    std::vector< std::vector< std::vector<int> > > psi;
    std::vector<int> start;

 private:
    struct step_t {
        int t;
        int i;
        int nth;
    };
    std::stack<step_t> stack;
    std::vector<int> q_star;
};

class Viterbi {
 public:
    typedef enum {
        PATH_TYPE_ALL = 0,
        PATH_TYPE_FIRST,
        PATH_TYPE_LAST,
    } path_type_e;

    Viterbi();

    ~Viterbi();

    ViterbiResult calculate(my_float o[], int n,
            const my_float log_init_prob[], int m,
            path_type_e path_type);

    void set_theta(int theta);

 private:
    int m;

    int n;

    int theta;

    path_type_e path_type;

    my_float *trans_prob_cache;

    my_float* o;

    int cached_t;

    int cached_i;

    my_float cal_trans_prob_cache(int m, int n, int i, int j, int t);

    my_float cal_trans_prob_cache_lvy(int m, int n, int i, int j, int t,
            my_float *cache);

    void delete_cache();
};

}; /* namespace pcyh */

#endif  /* _VITERBI_HPP_ */
