include ../../Makefile.in
include ../../TMR_Common.mk

OBJS = bracket.o

%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} bracket.o ${TMR_LD_FLAGS} -o bracket

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf bracket *.o
