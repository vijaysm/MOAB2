DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_srcdir)/mhdf/include -I$(top_builddir)
check_PROGRAMS = h5test h5legacy
if PARALLEL_HDF5
  check_PROGRAMS += parallel
  INCLUDES += -I$(top_srcdir)/parallel
endif
TESTS = $(check_PROGRAMS)

LDADD = $(top_builddir)/libMOAB.la
h5test_SOURCES = h5file_test.cpp
h5legacy_SOURCES = h5legacy.cpp
parallel_SOURCES = parallel.cpp
