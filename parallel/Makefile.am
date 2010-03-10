# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libMOABpar.la
         
# Some variables
DEFS = $(DEFINES) -DIS_BUILDING_MB -DSRCDIR=$(srcdir)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(srcdir)/../mhdf/inc

# The directory in which to install headers
libMOABpar_la_includedir = $(includedir)

# Run parallel tests in parallel
if USE_MPIEXEC
  OLD_ENVIRONMENT := ${TESTS_ENVIRONMENT}
  TESTS_ENVIRONMENT = ${MPIEXEC} ${MPIEXEC_NP} ${NP} ${OLD_ENVIRONMENT}
endif

# Conditional sources
MOAB_PARALLEL_SRCS =
MOAB_PARALLEL_HDRS =
MOAB_PARALLEL_TEST = 
MOAB_PARALLEL_CHECK = 
if PARALLEL
# The list of source files, and any header files that do not need to be installed
  MOAB_PARALLEL_SRCS += \
     MBParallelComm.cpp \
     MBProcConfig.cpp \
     MBParallelData.cpp \
     ReadParallel.cpp \
     crystal.c errmem.c \
     transfer.c gs.c gs.h tuple_list.c \
     tuple_list.h \
     sort.c sort.h minmax.h transfer.h sort_imp.c

# The list of header files which are to be installed
  MOAB_PARALLEL_HDRS += \
     MBmpi.h MBmpi_config.h \
     MBParallelComm.hpp \
     MBProcConfig.hpp \
     MBParallelData.hpp \
     ReadParallel.hpp \
     crystal.h errmem.h types.h

  MOAB_PARALLEL_TEST += pcomm_unit parallel_unit_tests uber_parallel_test scdtest pcomm_serial

if PARALLEL_HDF5
#  libMOABpar_la_LIBADD = $(top_builddir)/mhdf/libmhdf.la
  INCLUDES += -I$(top_srcdir)/mhdf/include
  MOAB_PARALLEL_SRCS += WriteHDF5Parallel.cpp 
  MOAB_PARALLEL_HDRS += WriteHDF5Parallel.hpp
  
  MOAB_PARALLEL_TEST += parallel_hdf5_test mhdf_parallel parallel_write_test
endif

  MOAB_PARALLEL_CHECK += mbparallelcomm_test
  mbparallelcomm_test_SOURCES = mbparallelcomm_test.cpp
  mbparallelcomm_test_LDADD = $(top_builddir)/libMOAB.la
  mbparallelcomm_test_DEPENDENCIES = $(mbparallelcomm_test_LDADD)

endif

# The list of source files, and any header files that do not need to be installed
libMOABpar_la_SOURCES = \
   $(MOAB_PARALLEL_SRCS)

# The list of header files which are to be installed
libMOABpar_la_include_HEADERS = \
  $(MOAB_PARALLEL_HDRS)

# Tests and such

check_PROGRAMS = $(MOAB_PARALLEL_CHECK) $(MOAB_PARALLEL_TEST)
TESTS = $(MOAB_PARALLEL_TEST)
pcomm_unit_SOURCES = pcomm_unit.cpp
pcomm_unit_LDADD = $(top_builddir)/libMOAB.la
parallel_hdf5_test_SOURCES = parallel_hdf5_test.cc
parallel_hdf5_test_LDADD = ../libMOAB.la
mhdf_parallel_SOURCES = mhdf_parallel.c
mhdf_parallel_LDADD = ../mhdf/libmhdf.la
parallel_unit_tests_SOURCES = parallel_unit_tests.cpp
parallel_unit_tests_LDADD = ../libMOAB.la
parallel_write_test_SOURCES = parallel_write_test.cc
parallel_write_test_LDADD = ../libMOAB.la
uber_parallel_test_SOURCES = uber_parallel_test.cpp
uber_parallel_test_LDADD = ../libMOAB.la
pcomm_serial_SOURCES = pcomm_serial.cpp
pcomm_serial_LDADD = ../libMOAB.la

scdtest_SOURCES = scdtest.cpp
scdtest_LDADD = ../libMOAB.la


# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = mhdf_ll.h5m

EXTRA_DIST = ptest.cub

