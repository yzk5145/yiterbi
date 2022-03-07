from ctypes import *

def veterbi(data):
    c_veterbi = CDLL('../build/libviterbi.so')
    o =  (c_longdouble * len(data))(*data)
    q_star =  (c_int * len(data))()
    n = (c_int)(len(data))
    c_veterbi.viterbi(o, q_star, n)
    return q_star

data = [2.0, 1.0, 3.0]

for n in veterbi(data):
    print(n)

