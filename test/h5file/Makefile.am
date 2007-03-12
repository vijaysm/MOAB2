DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_srcdir)/mhdf/include -I$(top_builddir)
check_PROGRAMS = h5test
if PARALLEL_HDF5
  check_PROGRAMS += parallel
endif
h5test_SOURCES = h5file_test.cpp
h5test_LDADD = $(top_builddir)/libMOAB.la
h5test_DEPENDENCIES = $(top_builddir)/libMOAB.la
parallel_SOURCES = parallel.cpp
parallel_LDADD = $(top_builddir)/libMOAB.la
parallel_DEPENDENCIES = $(top_builddir)/libMOAB.la
TESTS = $(check_PROGRAMS)
