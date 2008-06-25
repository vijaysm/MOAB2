DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_srcdir)/mhdf/include -I$(top_builddir)
check_PROGRAMS = h5test h5legacy h5varlen h5sets_test
parallel_programs = parallel varlen_ll mhdf_parallel

if PARALLEL_HDF5
  check_PROGRAMS += $(parallel_programs)
  INCLUDES += -I$(top_srcdir)/parallel
endif

TESTS = $(check_PROGRAMS)
LDADD = $(top_builddir)/libMOAB.la

h5test_SOURCES = h5file_test.cpp
h5legacy_SOURCES = h5legacy.cpp
h5varlen_SOURCES = h5varlen.cpp
parallel_SOURCES = parallel.cpp
varlen_ll_SOURCES = varlen_ll.cpp
mhdf_parallel_SOURCES = mhdf_parallel.c
mhdf_parallel_LDADD = ../../mhdf/libmhdf.la
h5sets_test_SOURCES = h5sets_test.cpp
