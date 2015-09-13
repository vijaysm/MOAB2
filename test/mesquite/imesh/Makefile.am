@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/itaps/mesquite/libmsqimesh.la \
        $(top_builddir)/src/mesquite/libmesquite.la \
        $(IMESH_LIBS)
check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/src/mesquite/libmesquite.la \
                    $(top_builddir)/itaps/mesquite/libmsqimesh.la
TESTS = main
EXTRA_DIST = CMakeLists.txt

AM_CPPFLAGS += -I$(top_srcdir)/src/parallel \
							 -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps/imesh

if ENABLE_imesh

else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	feasible-newton-result.vtk \
			original.vtk \
			smart-laplacian-result.vtk
