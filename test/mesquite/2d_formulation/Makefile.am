@TEST_MAKE_INCLUDE@
MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR)

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = main 
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt initial.vtk 
MOSTLYCLEANFILES = 1-1.results.vtk \
                   1-2.results.vtk \
                   1-4.results.vtk \
                   1-8.results.vtk \
                   2-NE.results.vtk \
                   2-NW.results.vtk \
                   2-SE.results.vtk \
                   2-SW.results.vtk \
                   3-Left.results.vtk \
                   3-Mid.results.vtk \
                   3-Right.results.vtk

