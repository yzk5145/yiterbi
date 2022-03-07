#include <stdio.h>
#include <string.h>

#include "viterbi.h"

int main() {
    long double o[] = {2, 1, 3};
    // int o[] = {5, 4, 3, 2, 1};
    // n = T = len(o)
    int n = sizeof(o) / sizeof(long double);
    printf("n = %d\n", n); //n = 3 
    int q_star[n];
    memset(q_star, 0, sizeof(q_star));
    //starting address of memory to be filled, value to be filled, # of bytes to be filled 
    viterbi(o, q_star, n);
    for (int i = 0; i < n; i++) {
        printf(">%d: %d\n", i, q_star[i]);
    }
    return 0;
}
