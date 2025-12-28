TARGET=../../bin/eq3nr
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
LIBRARIES=../../lib/libeqlib.a ../../lib/libeqlibu.a ../../lib/libeqlibg.a
INCLUDES=-I../../include

include ../Makefile.exe.in
