# ============================================
#
# Make file for TMR_DIR/src/interfaces
#
# ============================================

include ../../Makefile.in
include ../../TMR_Common.mk

CXX_OBJS = TMRTriangleInterface.o

DIR=${TMR_DIR}/src/interfaces

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TMR_LIB} ${CXX_OBJS}

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
