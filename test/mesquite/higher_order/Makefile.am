include $(top_srcdir)/src/mesquite/make.inc

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS  = -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(MSQ_INCLUDES) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/test

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = main voshell
TESTS = $(check_PROGRAMS)

main_SOURCES = main.cpp
voshell_SOURCES = voshell.cpp

EXTRA_DIST = CMakeLists.txt \
             expected_linear_output.vtk \
             expected_quadratic_output.vtk \
             hour-quad8.vtk \
	     linear_input.vtk \
	     quadratic_input.vtk \
	     cutCube.vtk \
	     hex3Dpart.vtk \
	     homogeneousPart.vtk \
	     quadPart.vtk \
	     sphereCube.vtk \
	     sphereCylinder_1194_inv.vtk \
	     triPart.vtk

MOSTLYCLEANFILES = smooth_ho.vtk
