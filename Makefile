# BenchFFT-NTT - FFT/NTT Multiplication Library
# 
# Targets:
#   build     - Build C library and tests
#   test      - Run C tests
#   bench     - Run C benchmarks
#   clean     - Clean build artifacts
#   formal    - Build and run Lean4 formal verification

.PHONY: all build test bench clean formal formal-build formal-test

all: build test bench

build:
	$(MAKE) -C c build

test:
	$(MAKE) -C c test

bench:
	$(MAKE) -C c bench

clean:
	$(MAKE) -C c clean
	cd lean && rm -rf build .lake

formal-build:
	cd lean && lake build

formal-test:
	cd lean && lake exe Formal
