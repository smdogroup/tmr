include ../../Makefile.in
include ../../TMR_Common.mk

OBJS = bracket.o \
	bracket_opt.o \
	bracket_opt_noshells.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} bracket.o ${TMR_LD_FLAGS} -o bracket
	${CXX} bracket_opt.o ${TMR_LD_FLAGS} -o bracket_opt
	${CXX} bracket_opt_noshells.o ${TMR_LD_FLAGS} -o bracket_noshell

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf bracket bracket_opt bracket_noshell *.o
