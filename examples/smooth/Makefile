include ../../Makefile.in
include ../../TMR_Common.mk

#Flags for python compilation
PY_INCLUDE = ${PYTHON_INCLUDE} -I${NUMPY_DIR} -I${MPI4PY_DIR}
CYTHON_INCLUDE = -I${NUMPY_DIR} -I${MPI4PY_DIR}

# locate shared objects
CYTHON_SO = locate.${SO_EXT}

default: ${CYTHON_SO}

debug: TMR_CC_FLAGS=${TMR_DEBUG_CC_FLAGS}
debug: default

%.cpp:%.pyx
	cython --cplus $*.pyx

%.o:%.cpp
	${CXX} ${TMR_CC_FLAGS} ${PY_INCLUDE} -c $< -o $@

locate.${SO_EXT}: locate.cpp LocatePoint.o
	${CXX} ${SO_LINK_FLAGS} $< LocatePoint.o -L/usr/lib -llapack ${TMR_LD_FLAGS} ${PY_INCLUDE} -o $@

clean:
	${RM} *.${SO_EXT}
	${RM} *.o
