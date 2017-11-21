include ../../Makefile.in
include ../../TMR_Common.mk

OBJS = tritest.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} tritest.o ${TMR_LD_FLAGS} -o tritest

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf tritest *.o
