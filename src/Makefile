PROG=QMC
OBJECTS=Support.o Classes/QMC_state.o Classes/System_Model.o

MPICPP=mpic++
CPP_FLAGS=-std=c++11 -O0 -g3 -Wall -c 
LD_FLAGS=
INCLUDE=
LIBS=-lgsl -lgslcblas -lm

all: ${OBJECTS} ${PROG}
	
${PROG}: ${OBJECTS} QMC_main.cpp
	${MPICPP} ${LD_FLAGS} ${INCLUDE} QMC_main.cpp ${OBJECTS} ${LIBS} -o $@

%.o: %.cpp
	${MPICPP} ${INCLUDE} ${CPP_FLAGS} -o $@ $<

clean:
	rm -f ${OBJECTS} ${PROG}

