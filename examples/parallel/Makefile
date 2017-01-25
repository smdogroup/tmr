TMR_DIR=${HOME}/hg/tmr
include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = octant_test.o \
	quadrant_test.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} octant_test.o ${TMR_LD_FLAGS} -o octant_test
	${CXX} quadrant_test.o ${TMR_LD_FLAGS} -o quadrant_test

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf octant_test quadrant_test *.o
