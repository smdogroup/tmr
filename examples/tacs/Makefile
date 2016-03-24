TMR_DIR=${HOME}/hg/tmr
TACS_DIR=${HOME}/hg/tacs-dev

include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk
include ${TACS_DIR}/Makefile.in
include ${TACS_DIR}/TACS_Common.mk

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} ${TACS_CC_FLAGS} -c $< -o $*.o

default: tmr_demo.o
	${CXX} tmr_demo.o ${TACS_LD_FLAGS} ${TMR_LD_FLAGS}  -o tmr_demo

clean:
	rm tmr_demo *.o
