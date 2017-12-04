include ../../Makefile.in
include ../../TMR_Common.mk

OBJS = cylinder.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} cylinder.o ${TMR_LD_FLAGS} -o cylinder

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf *.o
