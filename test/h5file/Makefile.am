DEFS = $(DEFINES)
INCLUDES += -I$(top_builddir)
check_PROGRAMS = h5file_test
if PARALLEL
  check_PROGRAMS += parallel
endif
h5file_test_SOURCES = h5file_test.cpp
h5file_test_LDADD = $(top_builddir)/libMOAB.la
parallel_SOURCES = parallel.cpp
parallel_LDADD = $(top_builddir)/libMOAB.la
TESTS = $(check_PROGRAMS)
