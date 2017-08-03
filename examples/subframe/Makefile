TMR_DIR=${HOME}/hg/tmr
include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

EXEC_OBJS = subframe_opt_multi.o

OBJS = TMRMultiLoadTopoProblem.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS} ${EXEC_OBJS}
	${CXX} -o subframe_multi subframe_opt_multi.o ${OBJS} ${TMR_LD_FLAGS}

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf subframe_opt subframe_multi *.o
