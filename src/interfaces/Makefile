# ============================================
#
# Make file for TMR_DIR/src/interfaces
#
# ============================================

include ../../Makefile.in
include ../../TMR_Common.mk

# CXX_OBJS = TMRTriangleInterface.o \
# 	TMRPerfectMatchInterface.o \
# 	TMRTriangularize.o

CXX_OBJS = TMRTriangularize.o

DIR=${TMR_DIR}/src/interfaces

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS} predicates.o
	${AR} ${AR_FLAGS} ${TMR_LIB} ${CXX_OBJS} predicates.o

predicates.o: predicates.c
	${CC} -c predicates.c -o predicates.o
	
debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
