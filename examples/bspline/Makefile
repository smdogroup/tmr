include ../../Makefile.in
include ../../TMR_Common.mk

OBJS = bspline.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} bspline.o ${TMR_LD_FLAGS} -o bspline

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf bspline *.o
