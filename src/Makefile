# ============================================
#
# Make file for TMR_DIR/src/
#
# ============================================

include ../Makefile.in
include ../TMR_Common.mk

CXX_OBJS = TMRBase.o \
	TMROctant.o \
	TMROctForest.o \
	TMRQuadrant.o \
	TMRQuadForest.o \
	TMRGeometry.o \
	TMRTriangularize.o \
	TMRMesh.o \
	TMRMeshSmoothing.o \
	TMROctStiffness.o \
	TMRBspline.o \
	TMRTopology.o \
	TMRNativeTopology.o \
	TMR_STLTools.o \
	TMR_RefinementTools.o \
	TMR_TACSCreator.o

DIR=${TMR_DIR}/src

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TMR_LIB} ${CXX_OBJS}

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
