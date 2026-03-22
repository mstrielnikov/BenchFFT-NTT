CC      = gcc
CFLAGS        = -Wall -Wextra -O3 -std=c11 -I include -mavx2 -mfma -funroll-loops
CFLAGS_SCALAR = -Wall -Wextra -O3 -std=c11 -I include -fno-tree-vectorize
LDFLAGS = -lm

SRC_DIR   = src
OBJ_DIR   = obj
LIB_DIR   = lib
TEST_DIR  = test
BENCH_DIR = bench

# ── Source files ─────────────────────────────────────────────────────────────
# Only real implementation sources (copy files and other/ are intentionally excluded)
SRC_BIGINT      = $(SRC_DIR)/bigint.c
SRC_FFT_SCALAR  = $(SRC_DIR)/fft_split.c
SRC_FFT_AVX     = $(SRC_DIR)/fft_split_avx.c

# ── Test entry points ────────────────────────────────────────────────────────
# test.c includes fft_split_test.c (scalar + avx via #if HAS_AVX)
TEST_MAIN   = $(TEST_DIR)/test.c
TEST_BIN    = $(TEST_DIR)/test

# ── Bench entry points ───────────────────────────────────────────────────────
# bench.c includes fft_split_bench.c and ../src/fft_split_avx.c internally
BENCH_MAIN  = $(BENCH_DIR)/bench.c
BENCH_BIN   = $(BENCH_DIR)/bench

.PHONY: all build build-scalar build-avx test bench clean

all: build test bench

# ── Library targets ──────────────────────────────────────────────────────────

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# Scalar objects
$(OBJ_DIR)/bigint_scalar.o: $(SRC_BIGINT) | $(OBJ_DIR)
	$(CC) $(CFLAGS_SCALAR) -c -o $@ $<

$(OBJ_DIR)/fft_split_scalar.o: $(SRC_FFT_SCALAR) | $(OBJ_DIR)
	$(CC) $(CFLAGS_SCALAR) -c -o $@ $<

# AVX objects
$(OBJ_DIR)/bigint_avx.o: $(SRC_BIGINT) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR)/fft_split_avx.o: $(SRC_FFT_AVX) | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Static libraries (output into lib/, not obj/)
$(LIB_DIR)/libfft_scalar.a: $(OBJ_DIR)/bigint_scalar.o $(OBJ_DIR)/fft_split_scalar.o | $(LIB_DIR)
	ar rcs $@ $^

$(LIB_DIR)/libfft_avx.a: $(OBJ_DIR)/bigint_avx.o $(OBJ_DIR)/fft_split_avx.o | $(LIB_DIR)
	ar rcs $@ $^

build-scalar: $(LIB_DIR)/libfft_scalar.a
build-avx:    $(LIB_DIR)/libfft_avx.a
build: build-scalar build-avx

# ── Test ─────────────────────────────────────────────────────────────────────
# test.c textually includes fft_split_test.c (which includes both src .c files
# directly), so we compile everything in one translation unit with -mavx2 so
# that the #if HAS_AVX branch is active.

$(TEST_BIN): $(TEST_MAIN) $(SRC_BIGINT)
	$(CC) $(CFLAGS) -o $@ $(TEST_MAIN) $(SRC_BIGINT) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

# ── Bench ────────────────────────────────────────────────────────────────────
# bench.c textually includes fft_split_bench.c and ../src/fft_split_avx.c,
# so we only need bigint.c as an explicit extra source.

$(BENCH_BIN): $(BENCH_MAIN) $(SRC_BIGINT)
	$(CC) $(CFLAGS) -o $@ $(BENCH_MAIN) $(SRC_BIGINT) $(LDFLAGS)

bench: $(BENCH_BIN)
	./$(BENCH_BIN)

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OBJ_DIR) $(LIB_DIR) $(TEST_BIN) $(BENCH_BIN)
