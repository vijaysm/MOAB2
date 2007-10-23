# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libMOABpar.la
         
# Some variables
DEFS = $(DEFINES) -DIS_BUILDING_MB
INCLUDES += -I$(top_builddir) 

# The directory in which to install headers
libMOABpar_la_includedir = $(includedir)


# Conditional sources
MOAB_PARALLEL_SRCS =
MOAB_PARALLEL_HDRS =
if USE_MPI
  MOAB_PARALLEL_SRCS += \
     MBParallelComm.cpp \
     MBParallelConventions.h \
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
endif

endif

# The list of source files, and any header files that do not need to be installed
libMOABpar_la_SOURCES = \
   $(MOAB_PARALLEL_SRCS)

# The list of header files which are to be installed
libMOABpar_la_include_HEADERS = \
  $(MOAB_PARALLEL_HDRS)

# Tests and such

#moab_test_SOURCES = MBTest.cpp
#moab_test_LDADD = $(top_builddir)/libMOAB.la
#moab_test_DEPENDENCIES = test/mb_big_test.g test/cell1.gen test/cell2.gen $(top_builddir)/libMOAB.la

