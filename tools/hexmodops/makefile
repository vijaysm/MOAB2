#
MOAB_DIR = ${HOME}/MOAB
NETCDF_DIR = /cubit/netcdf/lib
DEBUG_FLAG = -g
#

testops: testops.o
	${CXX} ${DEBUG_FLAG} -o testops testops.o -L${MOAB_DIR} -lMOAB -L${NETCDF_DIR} -lnetcdf_c++ -lnetcdf

makeops: makeops.o
	${CXX} ${DEBUG_FLAG} -o makeops makeops.o -L${MOAB_DIR} -lMOAB -L${NETCDF_DIR} -lnetcdf_c++ -lnetcdf

.SUFFIXES: .cpp .cxx .o

.cpp.o:
	${CXX} -I${MOAB_DIR} -c ${DEBUG_FLAG} $<


