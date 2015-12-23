include $(top_srcdir)/src/mesquite/make.inc

LDADD = $(top_builddir)/itaps/mesquite/libmsqitaps.la \
        $(top_builddir)/src/mesquite/libmesquite.la \
        $(top_builddir)/itaps/imesh/libiMesh.la \
        ${IGEOM_LDFLAGS} \
        ${CGM_LDFLAGS} ${CGM_LIBS}

check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/src/mesquite/libmesquite.la \
                    $(top_builddir)/itaps/mesquite/libmsqitaps.la
TESTS = main
EXTRA_DIST = CMakeLists.txt

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(MSQ_INCLUDES) \
							 -I$(top_srcdir)/src \
							 -I$(top_srcdir)/test \
               -I$(top_srcdir)/src/parallel \
							 -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_srcdir)/itaps/mesquite \
               -I$(top_builddir)/itaps/imesh

if ENABLE_imesh

else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	feasible-newton-result.vtk \
			original.vtk \
			smart-laplacian-result.vtk
