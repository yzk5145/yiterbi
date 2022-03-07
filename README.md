# YITERBI
HMM viterbi decoding implementation written in C++ and wrapped by python.
Refer to Section 6.1.2.1 HMM.viterbi and Algorithm 3: HMM.viterbi in 
```
@book{kuo2019singling,
  title={Singling-Out vs. Blending-In: Outlier Detection and Differential Privacy in Data},
  author={Kuo, Yu-Hsuan},
  year={2019},
  publisher={The Pennsylvania State University}
}
```

## Build Shared Library

generate build/libviterbi.so:

```bash
$ make
```


## Run the Test
compile build/test:
```bash
$ make test
```
compile and run build/test:
```bash
$ make test_run
```
compile and run build/test with memory analyzer (Valgrind):
```bash
$ make test_valgrind
```
