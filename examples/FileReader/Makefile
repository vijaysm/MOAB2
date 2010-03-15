MOAB_LIB_DIR = /homes/iulian/lib/MOAB/lib
CXX=/homes/sharp/3rdparty/mpich2/mpich2-1.1.1p1/gcc/bin/mpicxx
include ${MOAB_LIB_DIR}/moab.make

fread : FileRead.o
	${CXX} -g $< ${MOAB_LIBS_LINK} -o $@

.cpp.o : 
	${CXX} -g ${MOAB_INCLUDES} -c $<

