include $(top_srcdir)/src/mesquite/make.inc

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(MSQ_INCLUDES) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/test 

LDADD = $(top_builddir)/src/mesquite/libmesquite.la
check_PROGRAMS = main setup
main_SOURCES = main.cpp
main_DEPENDENCIES = $(LDADD)
setup_SOURCES = setup.cpp
setup_DEPENDENCIES = $(LDADD)
TESTS = main
EXTRA_DIST = CMakeLists.txt input.vtk

MOSTLYCLEANFILES = 

#input.vtk: mesh.vtk setup$(EXEEXT)
#	./setup(EXEEXT) 1 $(srcdir)/mesh.vtk input.vtk
