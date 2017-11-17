include ../../Makefile.in
include ../../TMR_Common.mk

NETGEN_INCLUDE=-I${HOME}/packages/netgen/include/
NG_LD_PATH=${HOME}/packages/netgen/lib
NETGEN_LD_FLAGS=-L${NG_LD_PATH} -Wl,-rpath,${NG_LD_PATH} -lnglib 

OBJS = step_import.o step_to_netgen.o

# Create a new rule for the code that requires both TACS and TMR
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} ${NETGEN_INCLUDE} -c $< -o $*.o

default: ${OBJS}
	${CXX} step_import.o ${TMR_LD_FLAGS} -o step_import
	${CXX} step_to_netgen.o ${TMR_LD_FLAGS} ${NETGEN_LD_FLAGS} -o step_to_netgen

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	rm -rf step_import *.o
