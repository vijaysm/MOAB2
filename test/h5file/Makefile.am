DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_srcdir)/mhdf/include -I$(top_builddir)
check_PROGRAMS = h5test h5legacy h5varlen h5sets_test

TESTS = $(check_PROGRAMS)
LDADD = $(top_builddir)/libMOAB.la

h5test_SOURCES = h5file_test.cpp
h5legacy_SOURCES = h5legacy.cpp
h5varlen_SOURCES = h5varlen.cpp
h5sets_test_SOURCES = h5sets_test.cpp
