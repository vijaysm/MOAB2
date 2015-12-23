include $(top_srcdir)/src/mesquite/make.inc
MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS  = -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(MSQ_INCLUDES) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/test

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = exo_to_vtk vtk_to_exo
exo_to_vtk_SOURCES = exo_to_vtk.cpp
vtk_to_exo_SOURCES = vtk_to_exo.cpp
exo_to_vtk_DEPENDENCIES = $(LDADD)
vtk_to_exo_DEPENDENCIES = $(LDADD)
EXTRA_DIST = CMakeLists.txt
