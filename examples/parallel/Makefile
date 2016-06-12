TMR_DIR=${HOME}/hg/tmr
TACS_DIR=${HOME}/hg/tacs
PAROPT_DIR=${HOME}/hg/paropt

include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk
include ${TACS_DIR}/Makefile.in
include ${TACS_DIR}/TACS_Common.mk
include ${PAROPT_DIR}/Makefile.in
include ${PAROPT_DIR}/ParOpt_Common.mk

OBJS = parallel.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} ${TACS_CC_FLAGS} ${PAROPT_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} ${OBJS} ${TACS_LD_FLAGS} ${PAROPT_LD_FLAGS} ${TMR_LD_FLAGS} -o parallel

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: PAROPT_CC_FLAGS=${PAROPT_DEBUG_CC_FLAGS}
debug: default

clean:
	rm parallel *.o
