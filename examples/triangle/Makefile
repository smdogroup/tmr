TMR_DIR=${HOME}/hg/tmr
include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = armtest.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
#	${CXX} tritest.o ${TMR_LD_FLAGS} -o tritest
	${CXX} armtest.o -pg ${TMR_LD_FLAGS} -o armtest

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf tritest armtest *.o
