BUILT_SOURCES = iGeom_protos.h

MESHDIR = $(top_srcdir)/MeshFiles/unittest/iGeom

AM_CPPFLAGS = -DSRCDIR=$(srcdir) \
	      -DMESHDIR=$(MESHDIR) \
              -I$(top_srcdir)/src \
              -I$(top_builddir)/src \
              -I.. -I$(srcdir)/.. \
              -I$(top_srcdir)/itaps/imesh \
              -I$(top_builddir)/itaps/imesh 


lib_LTLIBRARIES = libiGeomMOAB.la

libiGeomMOAB_la_includedir = $(includedir)

#libiGeomMOAB_la_LIBADD = $(top_builddir)/libMOAB.la

libiGeomMOAB_la_SOURCES = \
        iGeom_MOAB.hpp \
	iGeom_MOAB.cpp \
	SmoothCurveEval.hpp SmoothCurveEval.cpp \
	SmoothFaceEval.hpp SmoothFaceEval.cpp 

libiGeomMOAB_la_include_HEADERS = \
	iGeom.h \
	iGeom_protos.h \
	iGeom_f.h

# stuff for test program
check_PROGRAMS = testgeom  testSmoothGeom 
testgeom_SOURCES = testgeom.cc
testgeom_DEPENDENCIES = libiGeomMOAB.la \
                        $(top_builddir)/itaps/imesh/libiMesh.la \
                        $(top_builddir)/src/libMOAB.la 
testgeom_LDADD = libiGeomMOAB.la \
                 $(top_builddir)/itaps/imesh/libiMesh.la \
                 $(top_builddir)/src/libMOAB.la \
                 ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}

testSmoothGeom_SOURCES = testSmooth.cc 
testSmoothGeom_DEPENDENCIES = libiGeomMOAB.la \
                        $(top_builddir)/itaps/imesh/libiMesh.la \
                        $(top_builddir)/src/libMOAB.la 
testSmoothGeom_LDADD = libiGeomMOAB.la \
                 $(top_builddir)/itaps/imesh/libiMesh.la \
                 $(top_builddir)/src/libMOAB.la \
                 ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}
                 
TESTS = $(check_PROGRAMS)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iGeom-Defs.inc

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	$(AM_V_at)echo "IGEOM_DIR=${cfgdir}/.." >> $(DESTDIR)$(cfgdir)/iGeom-SIDL-Defs.inc
	$(AM_V_at)echo "IGEOM_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	$(AM_V_at)echo "IGEOM_INCLUDEDIR2=" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	$(AM_V_at)echo "IGEOM_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	$(AM_V_at)echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	$(AM_V_at)echo "IGEOM_LIBS=-L${libdir} -liGeomMOAB" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc

iGeom_protos.h: iGeom.h $(top_srcdir)/itaps/mkprotos.sh Makefile.am
	$(AM_V_GEN)$(top_srcdir)/itaps/mkprotos.sh iGeom MOAB $< $@ MOAB_FCDefs.h
