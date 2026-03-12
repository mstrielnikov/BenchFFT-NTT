# BenchFFT-NTT - FFT/NTT Multiplication Library
# 
# Targets:
#   build     - Build C library and tests
#   test      - Run C tests
#   bench     - Run C benchmarks
#   bench-rust - Run Rust benchmarks
#   bench-all - Run all benchmarks (C + Rust)
#   clean     - Clean build artifacts
#   formal    - Build and run Lean4 formal verification

.PHONY: all build test bench bench-rust bench-all clean formal formal-build formal-test

all: build test bench

build:
	$(MAKE) -C c build

test:
	$(MAKE) -C c test

bench:
	$(MAKE) -C c bench

bench-rust:
	cd rust && cargo build --release --quiet
	@echo "=== RUST ==="
	cd rust && ./target/release/bench

bench-all:
	$(MAKE) -C c bench-both
	$(MAKE) bench-rust

clean:
	$(MAKE) -C c clean
	cd rust && cargo clean
	cd lean && rm -rf build .lake

formal-build:
	cd lean && lake build

formal-test:
	cd lean && lake exe Formal
