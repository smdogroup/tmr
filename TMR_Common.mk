include ${TACS_DIR}/Makefile.in
include ${TACS_DIR}/TACS_Common.mk

TMR_LIB = ${TMR_DIR}/lib/libtmr.a

TMR_INCLUDE = -I${TMR_DIR}/src -I${TMR_DIR}/src/interfaces ${TACS_INCLUDE} 

# Set the compiler flags for TMR
TMR_CC_FLAGS = ${TMR_FLAGS} ${TMR_INCLUDE} ${TRIANGLE_INCLUDE} ${BLOSSOM_INCLUDE} ${METIS_INCLUDE} ${TACS_OPT_CC_FLAGS}
TMR_DEBUG_CC_FLAGS = ${TMR_DEBUG_FLAGS} ${TMR_INCLUDE} ${TRIANGLE_INCLUDE} ${BLOSSOM_INCLUDE} ${METIS_INCLUDE} ${TACS_DEBUG_CC_FLAGS}

# Set the compiler flags
TMR_LD_FLAGS = ${TMR_LD_CMD} ${TRIANGLE_LIB} ${BLOSSOM_LIB} ${METIS_LIB} ${TACS_LD_FLAGS}

# This is the one rule that is used to compile all the source
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
