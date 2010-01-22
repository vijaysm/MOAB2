# define the source: SRC=trunk or SRC=${MOAB_DIR}/branches/perf, etc
SRC=${MOAB_DIR}/branches/perf
cd ${SRC}
make distclean
rm -f config.log make.log
autoreconf -fi
#
# define the build: BUILD=gcc, etc
BUILD=gcc
mkdir -p ${SRC}/build/${BUILD}
cd ${SRC}/build/${BUILD}
#
export MPI_PREFIX=${SHARP_DIR}/mpich2/mpich2-1.2.1/gcc 
export HDF5_PREFIX=${SHARP_DIR}/hdf5/hdf5-1.8.3/parallel/gcc
export ZLIB_PREFIX=${SHARP_DIR}/zlib/zlib-1.2.3/gcc
export SZIP_PREFIX=${SHARP_DIR}/szip/szip-2.1/gcc
CC=${MPI_PREFIX}/bin/mpicc CXX=${MPI_PREFIX}/bin/mpicxx ../../configure --with-mpi=${MPI_PREFIX} --with-hdf5=${HDF5_PREFIX} --with-zlib=${ZLIB_PREFIX} --with-szip=${SZIP_PREFIX} --disable-debug 2>&1 | tee configure.log
#
make 2>&1 | tee make.log
make check 2>&1 | tee make.check.log

