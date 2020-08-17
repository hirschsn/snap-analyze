
.SUFFIXES:
.SUFFIXES: .cpp

CXXFLAGS?=-O3 -march=native -mtune=native -Wall -Wextra

BIN=sana-bond df
OBJ=$(BIN:=.o)

all: $(BIN)

clean:
	rm -f $(BIN)

.cpp:
	$(CXX) -std=c++17 $(CXXFLAGS) $< -o $@

.PHONY: all clean
