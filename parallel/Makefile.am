# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libMOABpar.la
         
# Some variables
DEFS = $(DEFINES) -DIS_BUILDING_MB
INCLUDES += -I$(top_srcdir) -I$(top_builddir) 

# The directory in which to install headers
libMOABpar_la_includedir = $(includedir)


# Conditional sources
MOAB_PARALLEL_SRCS =
MOAB_PARALLEL_HDRS =
MOAB_PARALLEL_TEST = 
MOAB_PARALLEL_CHECK = 
if PARALLEL
  MOAB_PARALLEL_SRCS += \
     MBParallelComm.cpp \
     MBProcConfig.cpp \
     MBParallelData.cpp \
     ReadParallel.cpp \
     crystal.c crystal.h errmem.h errmem.c \
     transfer.c gs.c gs.h tuple_list.c \
     tuple_list.h types.h \
     sort.c sort.h

  MOAB_PARALLEL_HDRS += \
     MBParallelComm.hpp \
     MBProcConfig.hpp \
     MBParallelData.hpp \
     ReadParallel.hpp 

if PARALLEL_HDF5
#  libMOABpar_la_LIBADD = $(top_builddir)/mhdf/libmhdf.la
  INCLUDES += -I$(top_srcdir)/mhdf/include
  MOAB_PARALLEL_SRCS += WriteHDF5Parallel.cpp 
  MOAB_PARALLEL_HDRS += WriteHDF5Parallel.hpp
  
  MOAB_PARALLEL_TEST += parallel_hdf5_test
endif

  MOAB_PARALLEL_CHECK += mbparallelcomm_test
  mbparallelcomm_test_SOURCES = mbparallelcomm_test.cpp
  mbparallelcomm_test_LDADD = $(top_builddir)/libMOAB.la
  mbparallelcomm_test_DEPENDENCIES = $(mbparallelcomm_test_LDADD)

  MOAB_PARALLEL_TEST += pcomm_unit
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
parallel_hdf5_test_CPPFLAGS = -DSRCDIR=$(srcdir)

