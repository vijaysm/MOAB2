LDADD = $(top_builddir)/src/libMOAB.la
AM_CPPFLAGS += -DIS_BUILDING_MB \
               -DSRCDIR=$(srcdir) \
               -I$(builddir) \
               -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_builddir)/src/moab \
               -I$(top_srcdir)/src/moab \
               -I$(top_srcdir)/src/parallel/moab

# Run parallel tests in parallel
if USE_MPIEXEC
  OLD_ENVIRONMENT := ${TESTS_ENVIRONMENT}
  TESTS_ENVIRONMENT = ${MPIEXEC} ${MPIEXEC_NP} ${NP} ${OLD_ENVIRONMENT} 
endif

TESTS = pcomm_unit \
        parallel_unit_tests \
        uber_parallel_test \
        scdtest \
        pcomm_serial \
        $(HDF5_TESTS)

if PARALLEL_HDF5
  AM_CPPFLAGS += -I$(top_srcdir)/mhdf/include
  HDF5_TESTS = parallel_hdf5_test mhdf_parallel parallel_write_test
else
  HDF5_TESTS = 
endif

check_PROGRAMS = $(TESTS) mbparallelcomm_test

pcomm_unit_SOURCES = pcomm_unit.cpp
parallel_hdf5_test_SOURCES = parallel_hdf5_test.cc
mhdf_parallel_SOURCES = mhdf_parallel.c
parallel_unit_tests_SOURCES = parallel_unit_tests.cpp
parallel_write_test_SOURCES = parallel_write_test.cc
uber_parallel_test_SOURCES = uber_parallel_test.cpp
pcomm_serial_SOURCES = pcomm_serial.cpp
mbparallelcomm_test_SOURCES = mbparallelcomm_test.cpp
scdtest_SOURCES = scdtest.cpp

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = mhdf_ll.h5m

EXTRA_DIST = ptest.cub \
             64bricks_512hex.h5m 
             
