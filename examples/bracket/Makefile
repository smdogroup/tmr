TMR_DIR=${HOME}/hg/tmr
include ${TMR_DIR}/Makefile.in
include ${TMR_DIR}/TMR_Common.mk

OBJS = bracket.o \
	bracket_opt.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o

default: ${OBJS}
	${CXX} bracket.o ${TMR_LD_FLAGS} -o bracket
	${CXX} bracket_opt.o ${TMR_LD_FLAGS} -o bracket_opt

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf bracket bracket_opt *.o
