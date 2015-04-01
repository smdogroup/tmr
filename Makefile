
%.o: %.c
	${CXX} -g -c $< -o $*.o

default: TMROctant.o TMROctree.o tmr_demo.o
	${CXX} TMROctant.o TMROctree.o tmr_demo.o -o tmr_demo

clean:
	rm *.o
