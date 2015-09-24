@TEST_MAKE_INCLUDE@

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR)

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = par_hex_smooth_laplace
par_hex_smooth_laplace_SOURCES = par_hex_smooth_laplace.cpp
par_hex_smooth_laplace_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt

MOSTLYCLEANFILES = *.*.vtk
