TACS_DIR=${HOME}/git/tacs
include ${TACS_DIR}/Makefile.in
include ${TACS_DIR}/TACS_Common.mk
include ../../Makefile.in

#Flags for python compilation
PY_INCLUDE = ${PYTHON_INCLUDE} -I${NUMPY_DIR} -I${MPI4PY_DIR}
CYTHON_INCLUDE = -I${NUMPY_DIR} -I${MPI4PY_DIR}

# ParOpt shared objects
CYTHON_SO = ksFSDT.so

default: ${CYTHON_SO}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

%.cpp:%.pyx
	cython -I${TACS_DIR}/tacs/ --cplus $*.pyx --include-dir ${TACS_DIR}/tacs/

%.o:%.cpp
	${CXX} -I${TACS_DIR}/tacs/ ${TACS_CC_FLAGS} ${PY_INCLUDE} -c $< -o $@

ksFSDT.so: ksFSDT.o
	${CXX} ${SO_LINK_FLAGS} ksFSDT.o ${TACS_LD_FLAGS} -o $@

clean:
	${RM} *.so *.o 
