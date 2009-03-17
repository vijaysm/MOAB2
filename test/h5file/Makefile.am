
@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_srcdir)/mhdf/include -I$(top_builddir)
check_PROGRAMS = h5test h5legacy h5varlen h5sets_test h5regression

TESTS = $(check_PROGRAMS)
LDADD = $(top_builddir)/libMOAB.la

if HAVE_CGM
  LDADD += $(CGM_LTFLAGS) 
endif

h5test_SOURCES = h5file_test.cpp
h5legacy_SOURCES = h5legacy.cpp
h5varlen_SOURCES = h5varlen.cpp
h5sets_test_SOURCES = h5sets_test.cpp
h5regression_SOURCES = h5regression.cpp
