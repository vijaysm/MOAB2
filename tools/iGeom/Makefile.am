AUTOMAKE_OPTIONS = foreign
BUILT_SOURCES = iGeom_protos.h

IMESH_DIR = $(top_srcdir)/tools/iMesh
INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I. -I$(IMESH_DIR)

DEFS += $(DEFINES) -DSRCDIR=$(srcdir)

lib_LTLIBRARIES = libiGeomMOAB.la

libiGeomMOAB_la_includedir = $(includedir)

#libiGeomMOAB_la_LIBADD = $(top_builddir)/libMOAB.la

libiGeomMOAB_la_SOURCES = \
	iGeom_MOAB.cpp 

libiGeomMOAB_la_include_HEADERS = \
	iBase.h \
	iGeom.h \
	iGeom_protos.h \
	iBase_f.h \
	iGeom_f.h

# stuff for c test program
check_PROGRAMS = testgeom
testgeom_SOURCES = testgeom.cc
testgeom_DEPENDENCIES = libiGeomMOAB.la $(IMESH_DIR)/libiMesh.la $(top_builddir)/libMOAB.la
testgeom_LDADD = libiGeomMOAB.la $(IMESH_DIR)/libiMesh.la $(top_builddir)/libMOAB.la ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}

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
	echo "IGEOM_DIR=${cfgdir}/.." >> $(DESTDIR)$(cfgdir)/iGeom-SIDL-Defs.inc
	echo "IGEOM_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "IGEOM_INCLUDEDIR2=" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "IGEOM_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc
	echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iGeom-Defs.inc

iGeom_protos.h: iGeom.h mkprotos.sh
	$(srcdir)/mkprotos.sh iGeom MBCN $< $@ MBCN_FCDefs.h

iGeom_FCDefs.h: iBase_FCDefs.h
	cd ../.. && ./config.status tools/iGeom/iGeom_FCDefs.h