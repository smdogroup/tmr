TMR_DIR=${HOME}/hg/tmr

include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = tfi.o \
	quartersphere.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} ${TACS_CC_FLAGS} ${PAROPT_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} tfi.o ${TACS_LD_FLAGS} ${PAROPT_LD_FLAGS} ${TMR_LD_FLAGS} -o tfi
	${CXX} quartersphere.o ${TACS_LD_FLAGS} ${PAROPT_LD_FLAGS} ${TMR_LD_FLAGS} -o quartersphere

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm tfi quartersphere *.o
