#include <iostream>

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "viterbi.hpp"
using namespace std;
using namespace pcyh;
static void handler(int sig) {
    void *array[10];
    size_t size;

    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}


int main(void) {
    Viterbi viterbi;
    ViterbiResult result;
    signal(SIGSEGV, handler);

    my_float o[] = {2, 1, 3};
    my_float init_prob[] = {-0.8472978603872034, -1.2527629684953676, -1.7635885922613586, -2.456735772821304, -3.5553480614894135};
    int q_star[3] = {0};

    result = viterbi.calculate(o, 3, init_prob, 4, Viterbi::PATH_TYPE_ALL);
    while (result.next(q_star)) {
        for (int i = 0; i < 3; i++) {
            cout << q_star[i] << "";
        }
        cout << endl;
    }
    return 0;
}
