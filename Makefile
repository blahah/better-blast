SRC=./src
CXXFLAGS+=-I./seqan/core/include
CXXFLAGS+=-I./seqan/extras/include

default: all
all: main

main: bblast.o
	$(CXX) $(LDFLAGS) -o bblast bblast.o

bblast.o: $(SRC)/bb.cpp
	$(CXX) $(CXXFLAGS) -c -o bblast.o $(SRC)/bb.cpp

clean:
	rm -f bblast.o bblast

.PHONY: default all clean