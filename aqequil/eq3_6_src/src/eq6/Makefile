TARGET=../../bin/eq6
MODSOURCES=src/mod6pt.f src/mod6xf.f
SOURCES=$(wildcard src/*.f)
OBJECTS=$(subst src/,obj/,$(SOURCES:.f=.o))
MODOBJECTS=$(subst src/,obj/,$(MODSOURCES:.f=.o))
LIBRARIES=../../lib/libeqlib.a ../../lib/libeqlibg.a ../../lib/libeqlibu.a 
INCLUDES=-I../../include

include ../Makefile.exe.in
