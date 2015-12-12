
MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src $(BGDEF)-DMESHDIR=$(MESHDIR) 
check_PROGRAMS = dual_test
dual_test_SOURCES = dual_test.cpp
LDADD = $(top_builddir)/src/libMOAB.la
dual_test_DEPENDENCIES = $(LDADD)
TESTS = dual_test 
MOSTLYCLEANFILES = dual.vtk
