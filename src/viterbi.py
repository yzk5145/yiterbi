import logging
import math
from ctypes import *
import operator as op
import os 

def ncr(n, r):
    from functools import reduce
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom

def cal_m(data, bias=5):
    return int(max(data) + bias)

def cal_init_prob(m, n, k):
    return math.log(ncr(m - k + n -1, n -1)) - math.log(ncr(m + n, n))

class ViterbiResult:
    def __init__(self, obj, n):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.c_viterbi = CDLL(dir_path + '/../build/libviterbi.so')
        self.obj = obj
        self.n = n

    def __iter__(self):
        q_star = (c_int * self.n)()
        while self.c_viterbi.viterbi_result_next(self.obj, q_star):
            yield [int (i) for i in q_star]

(PATH_TYPE_ALL, PATH_TYPE_FIRST, PATH_TYPE_LAST) = (0, 1, 2)
class Viterbi:
    def __init__(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.c_viterbi = CDLL(dir_path + '/../build/libviterbi.so')
        self.c_obj = self.c_viterbi.viterbi_new()

    def __del__(self):
        self.c_viterbi.viterbi_delete(self.c_obj)

    def calculate(self, data, bias=1, theta=-1, path_type=PATH_TYPE_ALL):
        self.m = cal_m(data, bias)
        self.n = len(data)

        self.c_viterbi.viterbi_set_theta(self.c_obj, (c_int)(theta))
        logging.debug("m = {}".format(self.m))
        def is_in_theta_range(v):
            if theta < 0: return True
            #data[0] is float
            #else: return v in range(max(data[0] - theta, 0), min(data[0] + theta 3, self.m + 1))
            else: return v in range(max(int(data[0]) - theta, 0), min(int(data[0]) + theta + 1, self.m + 1))
        log_init_prob_list = [
                cal_init_prob(self.m, self.n, i)
                if is_in_theta_range(i)
                else 0
                for i in range(self.m + 1)
                ]
        o =  (c_longdouble * len(data))(*data)
        n = (c_int)(len(data))

        c_log_init_prob_list = (c_longdouble * len(log_init_prob_list))(*log_init_prob_list)
        del log_init_prob_list

        c_result_obj = self.c_viterbi.viterbi_calculate(self.c_obj, o, n, 
                c_log_init_prob_list, self.m, (c_int)(path_type))

        return ViterbiResult(c_result_obj, self.n)

def main():
    viterbi = Viterbi();
    result = viterbi.calculate([2,1,3], bias=2, theta=2,
            path_type=PATH_TYPE_ALL)
    for q_star in result:
        print(q_star)
        

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
