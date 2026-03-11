CC = gcc
CFLAGS = -Wall -Wextra -O3 -std=c11 -I include
LDFLAGS = -lm

SRC = src
OBJ = obj
TEST = test
LIB = lib

.PHONY: all clean build test bench

all: build

build: $(OBJ) $(LIB)/libfft.a

$(OBJ):
	mkdir -p $(OBJ)

$(LIB):
	mkdir -p $(LIB)

$(OBJ)/%.o: $(SRC)/%.c | $(OBJ)
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB)/libfft.a: $(OBJ)/bigint.o $(OBJ)/fft_split.o $(OBJ)/ntt_mont.o
	mkdir -p $(LIB)
	ar rcs $@ $^

test: build
	$(CC) $(CFLAGS) $(TEST)/test.c -L $(LIB) -lfft -o $(TEST)/test $(LDFLAGS)
	./$(TEST)/test

bench: build
	$(CC) $(CFLAGS) $(TEST)/bench.c -L $(LIB) -lfft -o $(TEST)/bench $(LDFLAGS)
	./$(TEST)/bench

clean:
	rm -rf $(OBJ) $(TEST)/test $(TEST)/bench $(LIB)
