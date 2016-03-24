TMR_LIB = ${TMR_DIR}/lib/libtmr.a

TMR_INCLUDE = -I${TMR_DIR}/src

# Set the compiler flags for TMR
TMR_CC_FLAGS = ${TMR_OPT_FLAGS} ${TMR_INCLUDE}
TMR_CC_DEBUG_FLAGS = ${TMR_CC_FLAGS} ${TMR_INCLUDE}

# Set the compiler flags
TMR_LD_FLAGS = -L${TMR_DIR}/lib/ -Wl,-rpath,${TMR_DIR}/lib -ltmr

# This is the one rule that is used to compile all the source
%.o: %.c
	${CXX} ${TMR_CC_FLAGS} -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo


