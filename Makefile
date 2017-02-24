include Makefile.in
include TMR_Common.mk

TMR_SUBDIRS = src \
	src/interfaces
TMR_OBJS := $(addsuffix /*.o, ${TMR_SUBDIRS})

default:
	echo "Building Real TMR"
	@for subdir in $(TMR_SUBDIRS) ; do \
	    echo "making $@ in $$subdir"; \
	    echo; (cd $$subdir && $(MAKE) TMR_DIR=${TMR_DIR}) || exit 1; \
	done
	${CXX} ${SO_LINK_FLAGS} ${TMR_OBJS} ${TMR_EXTERN_LIBS} -o ${TMR_DIR}/lib/libtmr.${SO_EXT} 
	python setup.py build_ext --inplace

debug:
	echo "Building Real TMR"
	@for subdir in $(TMR_SUBDIRS) ; do \
	    echo "making $@ in $$subdir"; \
	    echo; (cd $$subdir && $(MAKE) debug TMR_DIR=${TMR_DIR}) || exit 1; \
	done
	${CXX} ${SO_LINK_FLAGS} ${TMR_OBJS} ${TMR_EXTERN_LIBS} -o ${TMR_DIR}/lib/libtmr.${SO_EXT} 
	python setup.py build_ext --inplace

clean:
	${RM} lib/*.a lib/*.so
	@for subdir in $(TMR_SUBDIRS) ; do \
	    echo "making $@ in $$subdir"; \
	    echo; (cd $$subdir && $(MAKE) clean TMR_DIR=${TMR_DIR}) || exit 1; \
	done
	${RM} tmr/*.so
	${RM} tmr/*.cpp

