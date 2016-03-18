include ../tacs-dev/Makefile.in
include ../tacs-dev/TACS_Common.mk

OBJS = TMROctant.o TMROctree.o tmr_demo.o

default: ${OBJS}
	${CXX} -o tmr_demo tmr_demo.o TMROctant.o TMROctree.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	rm *.o
