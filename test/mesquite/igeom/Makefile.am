@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/itaps/mesquite/libmsqitaps.la \
        $(top_builddir)/src/mesquite/libmesquite.la \
        ${IGEOM_LDFLAGS} ${IGEOM_LIBS} \
        ${CGM_LDFLAGS} ${CGM_LIBS}

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/itaps/mesquite/libmsqitaps.la \
                    $(top_builddir)/src/mesquite/libmesquite.la 

TESTS = main
EXTRA_DIST = CMakeLists.txt 

AM_CPPFLAGS += -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(IGEOM_DEFINES) $(IGEOM_INCLUDES) \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/itaps/mesquite \
               -I$(top_srcdir)/src/parallel

if HAVE_CGM
else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	
