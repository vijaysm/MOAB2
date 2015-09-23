@TEST_MAKE_INCLUDE@

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR)

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = jacobi_1
jacobi_1_SOURCES = jacobi_1.cpp
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt
