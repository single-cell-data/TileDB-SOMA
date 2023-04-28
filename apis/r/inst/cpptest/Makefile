
## By using a Makefile we should get CXX etc set automagically

INCL = $(shell pkg-config --cflags tiledbsoma tiledb)
LIBS = $(shell pkg-config --libs tiledbsoma tiledb)

all:	min

min:	minimal.cpp
	$(CXX) -o $@ ${INCL} ${LIBS}
