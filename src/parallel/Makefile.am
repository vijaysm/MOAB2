# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Subdirectories to build
SUBDIRS = 

# Things to build
noinst_LTLIBRARIES = libMOABpar.la
         
# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB -DSRCDIR=$(srcdir) \
               -I$(srcdir)/.. -I.. -I$(srcdir)/../io \
               -I$(srcdir)/../io/mhdf/include  \
               -I$(builddir)/../moab -I../moab  \
               -I$(srcdir)/moab -Imoab

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
if PARALLEL
# The list of source files, and any header files that do not need to be installed
  MOAB_PARALLEL_SRCS += \
     MBParallelComm.cpp \
     MBProcConfig.cpp \
     MBParallelData.cpp \
     ReadParallel.cpp \
     ReadParallel.hpp \
     crystal.h errmem.h types.h \
     crystal.c errmem.c \
     transfer.c gs.c gs.h tuple_list.c \
     tuple_list.h \
     sort.c sort.h minmax.h transfer.h sort_imp.c

# The list of header files which are to be installed
  MOAB_PARALLEL_HDRS += \
     moab/MBmpi.h \
     moab/MBmpi_config.h \
     moab/MBParallelComm.hpp \
     moab/MBProcConfig.hpp \
     moab/MBParallelData.hpp \
     moab/MBParallelConventions.h 

if PARALLEL_HDF5
#  libMOABpar_la_LIBADD = $(top_builddir)/mhdf/libmhdf.la
  MOAB_PARALLEL_SRCS += WriteHDF5Parallel.cpp 
  MOAB_PARALLEL_HDRS += WriteHDF5Parallel.hpp
endif

endif

# The list of source files, and any header files that do not need to be installed
libMOABpar_la_SOURCES = \
   $(MOAB_PARALLEL_SRCS)

# The list of header files which are to be installed
nobase_libMOABpar_la_include_HEADERS = \
  $(MOAB_PARALLEL_HDRS)
