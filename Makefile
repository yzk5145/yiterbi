OUT_DIR := build
.PHONY: all

all: $(OUT_DIR)/ libs

libs: src/viterbi.c
	g++ -O3 -shared -o $(OUT_DIR)/libviterbi.so -fPIC src/viterbi.cpp -lm

.PHONY: test
test: src/viterbi.cpp test/test.cpp
	mkdir -p build
	g++ test/test.cpp src/viterbi.cpp -o $(OUT_DIR)/test -Isrc/ -Lbuild -g -rdynamic

test_run: test
	./build/test

test_valgrind: test
	valgrind --leak-check=yes --leak-check=full --show-leak-kinds=all --show-reachable=yes ./build/test 

$(OUT_DIR)/:
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(OUT_DIR)
	rm -rf logs
