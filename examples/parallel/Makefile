TMR_DIR=${HOME}/hg/tmr

include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = parallel.o \
	block_parallel.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} parallel.o ${TMR_LD_FLAGS} -o parallel
	${CXX} block_parallel.o ${TMR_LD_FLAGS} -o block_parallel

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf block_parallel parallel *.o
