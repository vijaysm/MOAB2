MOAB_DIR = ${HOME}/MOAB
CHACO_DIR = ${HOME}/Chaco-2.2

include ${MOAB_DIR}/MB.lin
include ${MOAB_DIR}/MB.local

INCLUDES = -I${CHACO_DIR}/code/main -I${MOAB_DIR}

mbchaco:	main.o ComputePartition.o makefile
		${CXX} ${DEBUG_FLAG} main.o ComputePartition.o -L${MOAB_DIR} -lMOAB ${NETCDF_LIBS} -L${CHACO_DIR}/code -lchaco -lm -o mbchaco

lint:
		lint ${IFLAG} ${FILES.c} -lm

alint:
		alint ${IFLAG} ${FILES.c} -lm

clean:
		rm -f */*.o ${DEST_DIR}/core

.cpp.o:
		${CXX} -c ${DEBUG_FLAG} ${INCLUDES} -o $*.o $*.cpp
