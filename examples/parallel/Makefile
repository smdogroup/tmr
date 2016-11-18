TMR_DIR=${HOME}/hg/tmr
PAROPT_DIR=${HOME}/hg/paropt

include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk
# include ${PAROPT_DIR}/Makefile.in
# include ${PAROPT_DIR}/ParOpt_Common.mk

OBJS = octant_test.o
#       parallel.o \
# 	block_parallel.o \
# 	TMRTopoProblem.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${PAROPT_CC_FLAGS} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} octant_test.o ${TMR_LD_FLAGS} -o octant_test
# 	${CXX} parallel.o ${TMR_LD_FLAGS} -o parallel
# 	${CXX} block_parallel.o TMRTopoProblem.o ${PAROPT_LD_FLAGS} ${TMR_LD_FLAGS} -o block_parallel

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: PAROPT_CC_FLAGS=${PAROPT_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf block_parallel parallel *.o
