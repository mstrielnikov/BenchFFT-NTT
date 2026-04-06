CC      = gcc
CFLAGS        = -Wall -Wextra -O3 -std=c11 -I include -mavx2 -mfma -funroll-loops -fPIC
CFLAGS_SCALAR = -Wall -Wextra -O3 -std=c11 -I include -fno-tree-vectorize -fPIC
LDFLAGS = -lm

SRC_DIR   = src
OBJ_DIR   = obj
LIB_DIR   = lib
TEST_DIR  = test
BENCH_DIR = bench

# ── Source files ─────────────────────────────────────────────────────────────
# Automatically collect all .c files in src/
ALL_SRCS       := $(wildcard $(SRC_DIR)/*.c)

# We identify AVX and ASM sources by their suffixes
SRCS_AVX       := $(filter %_avx.c, $(ALL_SRCS))
SRCS_ASM       := $(filter %_asm.c, $(ALL_SRCS))

# Everything else is considered scalar (including bigint.c, fft_mersenne.c, etc.)
SRCS_SCALAR    := $(filter-out $(SRCS_AVX) $(SRCS_ASM), $(ALL_SRCS)) src/nussbaumer.c

# Object files
OBJS_SCALAR    := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%_scalar.o, $(SRCS_SCALAR))
OBJS_AVX       := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%_avx.o, $(SRCS_AVX))
OBJS_ASM       := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%_asm.o, $(SRCS_ASM))

# ── Test entry points ────────────────────────────────────────────────────────
TEST_MAIN   = $(TEST_DIR)/test.c
TEST_BIN    = $(TEST_DIR)/test
TEST_DEPS   = $(wildcard $(TEST_DIR)/*.h) $(wildcard $(SRC_DIR)/*.c) $(wildcard include/*.h)

# ── Bench entry points ───────────────────────────────────────────────────────
BENCH_MAIN  = $(BENCH_DIR)/bench.c
BENCH_BIN   = $(BENCH_DIR)/bench
BENCH_DEPS  = $(wildcard $(BENCH_DIR)/*.h) $(wildcard $(SRC_DIR)/*.c) $(wildcard include/*.h)

.PHONY: all build build-scalar build-avx build-shared test bench clean verify

all: build test bench

# ── Library targets ──────────────────────────────────────────────────────────

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# Scalar pattern rule
$(OBJ_DIR)/%_scalar.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS_SCALAR) -c -o $@ $<

# AVX pattern rule
$(OBJ_DIR)/%_avx.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# ASM pattern rule (compiled with default CFLAGS for now, as AVX won't hurt scalar ASM)
$(OBJ_DIR)/%_asm.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Static libraries
$(LIB_DIR)/libfft_scalar.a: $(OBJS_SCALAR) | $(LIB_DIR)
	ar rcs $@ $^

$(LIB_DIR)/libfft_avx.a: $(OBJS_SCALAR) $(OBJS_AVX) $(OBJS_ASM) | $(LIB_DIR)
	ar rcs $@ $^

$(LIB_DIR)/libfft_both.so: $(OBJS_SCALAR) $(OBJS_AVX) $(OBJS_ASM) | $(LIB_DIR)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

build-scalar: $(LIB_DIR)/libfft_scalar.a
build-avx:    $(LIB_DIR)/libfft_avx.a
build-shared: $(LIB_DIR)/libfft_both.so
build: build-scalar build-avx build-shared

# ── Test ─────────────────────────────────────────────────────────────────────

$(TEST_BIN): $(TEST_MAIN) $(TEST_DEPS) $(SRC_DIR)/bigint.c
	$(CC) $(CFLAGS) -o $@ $(TEST_MAIN) $(SRC_DIR)/bigint.c $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

# ── Bench ────────────────────────────────────────────────────────────────────

$(BENCH_BIN): $(BENCH_MAIN) $(BENCH_DEPS) $(SRC_DIR)/bigint.c
	$(CC) $(CFLAGS) -o $@ $(BENCH_MAIN) $(SRC_DIR)/bigint.c $(LDFLAGS)

bench: $(BENCH_BIN)
	./$(BENCH_BIN)

# ── Lean Verification ────────────────────────────────────────────────────────

verify: build-shared
	python3 verify_runtime.py

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OBJ_DIR) $(LIB_DIR) $(TEST_BIN) $(BENCH_BIN)
