include ../tacs-dev/Makefile.in
include ../tacs-dev/TACS_Common.mk

CXX=mpicxx

FLAGS=-O3
DEBUG_FLAGS=-g

default: TMROctant.o TMROctree.o tmr_demo.o
	${CXX} TMROctant.o TMROctree.o tmr_demo.o -o tmr_demo

debug: FLAGS=${DEBUG_FLAGS}
debug: default

%.o: %.c
	${CXX} ${FLAGS} -c $< -o $*.o

clean:
	rm *.o
