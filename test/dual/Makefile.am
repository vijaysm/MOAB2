
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
check_PROGRAMS = dual_test
dual_test_SOURCES = dual_test.cpp
LDADD = $(top_builddir)/libMOAB.la
dual_test_DEPENDENCIES = $(top_builddir)/libMOAB.la
TESTS = dual_test 
