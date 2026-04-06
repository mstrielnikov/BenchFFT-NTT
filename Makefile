# BenchFFT-NTT - FFT/NTT Multiplication Library
#
# Structure:
#   c/           - C implementation
#   c/src/       - Implementation: bigint.c, fft_split.c, fft_split_avx.c
#   c/test/      - Tests: test_scalar.c, test_avx.c
#   c/bench/     - Benchmarks: bench_scalar.c, bench_avx.c
#   c/obj/       - Object files: *_scalar.o, *_avx.o
#   rust/        - Rust implementation
#   lean/        - Lean4 formal verification
#
# Targets:
#   build        - Build both scalar and AVX libraries
#   test         - Run both scalar and AVX tests
#   bench        - Run both scalar and AVX benchmarks
#   bench-rust   - Build and run Rust benchmarks
#   bench-all    - Run all benchmarks (C + Rust)
#   clean        - Clean build artifacts
#   formal-build - Build Lean4 formal verification
#   formal-test  - Run Lean4 formal verification

C_DIR = c

.PHONY: all build test bench bench-rust bench-all clean formal-build formal-test

all: build test bench

build:
	$(MAKE) -C $(C_DIR) build

test:
	$(MAKE) -C $(C_DIR) test

bench:
	$(MAKE) -C $(C_DIR) bench

bench-rust:
	cd rust && cargo build --release --quiet
	@echo "=== RUST ==="
	cd rust && ./target/release/bench

bench-all:
	$(MAKE) -C $(C_DIR) bench
	$(MAKE) bench-rust

clean:
	$(MAKE) -C $(C_DIR) clean
	cd rust && cargo clean
	cd lean && rm -rf build .lake

formal-build:
	cd lean && lake build

formal-test:
	cd lean && lake exe Formal
