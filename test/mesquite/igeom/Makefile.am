@TEST_MAKE_INCLUDE@
LDADD = $(top_builddir)/itaps/mesquite/libmsqigeom.la \
        $(top_builddir)/src/mesquite/libmesquite.la \
        ${IGEOM_LDFLAGS} \
        ${CGM_LDFLAGS} ${CGM_LIBS}

check_PROGRAMS = main
main_SOURCES = main.cpp
main_DEPENDENCIES = $(top_builddir)/itaps/mesquite/libmsqigeom.la \
										$(top_builddir)/src/mesquite/libmesquite.la 

TESTS = main
EXTRA_DIST = CMakeLists.txt 

AM_CPPFLAGS += $(IGEOM_DEFINES) $(IGEOM_INCLUDES) \
							 -I$(top_srcdir)/itaps/imesh \
							 -I$(top_builddir)/itaps/imesh \
							 -I$(top_srcdir)/src/parallel

if HAVE_CGM
else
main$(EXEEXT):
	echo "ITAPS support disabled.  Reconfigure to build this test."
endif

MOSTLYCLEANFILES =	
