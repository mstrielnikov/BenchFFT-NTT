CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c11 -I include -mavx2 -mfma -funroll-loops
CFLAGS_SCALAR = -Wall -Wextra -O3 -std=c11 -I include -fno-tree-vectorize
CFLAGS_AVX = -Wall -Wextra -O3 -std=c11 -I include -mavx2 -mfma -funroll-loops -DBUILD_AVX
LDFLAGS = -lm

SRC = src
OBJ = obj
TEST = test
LIB = lib

.PHONY: all clean build build-scalar build-avx build-mersenne build-mersenne-avx test bench bench-scalar bench-avx bench-mersenne bench-mersenne-avx

all: build

build: $(OBJ) $(LIB)/libfft.a

build-scalar: $(OBJ) $(LIB)/libfft_scalar.a

build-avx: $(OBJ) $(LIB)/libfft_avx.a

build-mersenne: $(OBJ) $(LIB)/libfft_mersenne.a

build-mersenne-avx: $(OBJ) $(LIB)/libfft_mersenne_avx.a

$(OBJ):
	mkdir -p $(OBJ)

$(LIB):
	mkdir -p $(LIB)

$(OBJ)/%.o: $(SRC)/%.c | $(OBJ)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ)/%_scalar.o: $(SRC)/%.c | $(OBJ)
	$(CC) $(CFLAGS_SCALAR) -c -o $@ $<

$(OBJ)/%_avx.o: $(SRC)/%_avx.c | $(OBJ)
	$(CC) $(CFLAGS_AVX) -c -o $@ $<

$(OBJ)/bigint.o $(OBJ)/fft_split.o $(OBJ)/ntt_mont.o: $(SRC)/bigint.c $(SRC)/fft_split.c $(SRC)/ntt_mont.c

$(LIB)/libfft.a: $(OBJ)/bigint.o $(OBJ)/fft_split.o $(OBJ)/ntt_mont.o
	mkdir -p $(LIB)
	ar rcs $@ $^

$(LIB)/libfft_scalar.a: $(OBJ)/bigint.o $(OBJ)/fft_split_scalar.o $(OBJ)/ntt_mont_scalar.o
	mkdir -p $(LIB)
	ar rcs $@ $^

$(LIB)/libfft_avx.a: $(OBJ)/bigint.o $(OBJ)/fft_split_avx.o $(OBJ)/ntt_mont_avx.o
	mkdir -p $(LIB)
	ar rcs $@ $^

$(LIB)/libfft_mersenne.a: $(OBJ)/bigint.o $(OBJ)/fft_split.o $(OBJ)/ntt_mersenne.o
	mkdir -p $(LIB)
	ar rcs $@ $^

$(LIB)/libfft_mersenne_avx.a: $(OBJ)/bigint.o $(OBJ)/fft_split_avx.o $(OBJ)/ntt_mersenne_avx.o
	mkdir -p $(LIB)
	ar rcs $@ $^

test: build
	$(CC) $(CFLAGS) $(TEST)/test.c -L $(LIB) -lfft -o $(TEST)/test $(LDFLAGS)
	./$(TEST)/test

test-scalar: build-scalar
	$(CC) $(CFLAGS_SCALAR) $(TEST)/test.c -L $(LIB) -lfft_scalar -o $(TEST)/test_scalar $(LDFLAGS)
	./$(TEST)/test_scalar

bench: build
	$(CC) $(CFLAGS) $(TEST)/bench.c -L $(LIB) -lfft -o $(TEST)/bench $(LDFLAGS)
	./$(TEST)/bench

bench-scalar: build-scalar
	$(CC) $(CFLAGS_SCALAR) $(TEST)/bench.c -L $(LIB) -lfft_scalar -o $(TEST)/bench_scalar $(LDFLAGS)
	./$(TEST)/bench_scalar

bench-avx: build-avx
	$(CC) $(CFLAGS_AVX) $(TEST)/bench.c -L $(LIB) -lfft_avx -o $(TEST)/bench_avx $(LDFLAGS)
	./$(TEST)/bench_avx

bench-mersenne: build-mersenne
	$(CC) $(CFLAGS) $(TEST)/bench.c -L $(LIB) -lfft_mersenne -o $(TEST)/bench_mersenne $(LDFLAGS)
	./$(TEST)/bench_mersenne

bench-mersenne-avx: build-mersenne-avx
	$(CC) $(CFLAGS_AVX) $(TEST)/bench.c -L $(LIB) -lfft_mersenne_avx -o $(TEST)/bench_mersenne_avx $(LDFLAGS)
	./$(TEST)/bench_mersenne_avx

clean:
	rm -rf $(OBJ) $(TEST)/test $(TEST)/test_scalar $(TEST)/bench $(TEST)/bench_scalar $(TEST)/bench_avx $(TEST)/bench_mersenne $(LIB)
