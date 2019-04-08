include ${TACS_DIR}/Makefile.in
include ${TACS_DIR}/TACS_Common.mk
include ${PAROPT_DIR}/Makefile.in
include ${PAROPT_DIR}/ParOpt_Common.mk
include ${EGADS_DIR}/Makefile.in
include ${EGADS_DIR}/EGADS_Common.mk

TMR_LIB = ${TMR_DIR}/lib/libtmr.a

TMR_INCLUDE = -I${TMR_DIR}/src -I${TMR_DIR}/src/interfaces -I${TMR_DIR}/src/topology \
	${EGADS_INCLUDE} ${TACS_INCLUDE} ${PAROPT_INCLUDE} ${OPENCASCADE_INCLUDE} ${NETGEN_INCLUDE}

# Set the compiler flags for TMR
TMR_CC_FLAGS = ${TMR_FLAGS} ${TMR_INCLUDE} ${BLOSSOM_INCLUDE} ${TACS_OPT_CC_FLAGS} ${EGADS_CC_FLAGS}
TMR_DEBUG_CC_FLAGS = ${TMR_DEBUG_FLAGS} ${TMR_INCLUDE} ${BLOSSOM_INCLUDE} ${TACS_DEBUG_CC_FLAGS} ${EGADS_DEBUG_CC_FLAGS}

# Set the compiler flags
TMR_EXTERN_LIBS = ${BLOSSOM_LIB} ${TACS_LD_FLAGS} ${PAROPT_LD_FLAGS} ${EGADS_LD_FLAGS} ${OPENCASCADE_LIB_PATH} ${OPENCASCADE_LIBS} ${NETGEN_LD_FLAGS}
TMR_LD_FLAGS = ${TMR_LD_CMD} ${TMR_EXTERN_LIBS}

# This is the one rule that is used to compile all the source
%.o: %.cpp
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.cpp successfully ---"
	@echo
