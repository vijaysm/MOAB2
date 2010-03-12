EXTRA_DIST = v3_dodec.h5m

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -I.. -I$(srcdir)/.. \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
	       -I$(top_builddir)/src/parallel 
            
check_PROGRAMS = h5test h5legacy h5varlen h5sets_test h5regression h5partial

TESTS = $(check_PROGRAMS)
LDADD = $(top_builddir)/src/libMOAB.la

h5test_SOURCES = h5file_test.cpp
h5legacy_SOURCES = h5legacy.cpp
h5varlen_SOURCES = h5varlen.cpp
h5sets_test_SOURCES = h5sets_test.cpp
h5regression_SOURCES = h5regression.cpp
h5partial_SOURCES = h5partial.cpp
