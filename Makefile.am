@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@
@FBIGEOM_CONFIG_OPTIONS@

AUTOMAKE_OPTIONS = foreign subdir-objects
BUILT_SOURCES = itaps/irel/iRel_protos.h

AM_CPPFLAGS = -DSRCDIR=$(srcdir)

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = test

AM_CPPFLAGS += -I$(top_srcdir)/src -I$(top_srcdir)/itaps/irel -I$(top_builddir)/itaps/irel -I$(top_builddir) -I$(top_srcdir)/src/lasso
DEFS += $(DEFINES)
LIBS =

libiRel_la_SOURCES = \
	src/lasso/ArrayManager.hpp \
	src/AssocPair.cpp \
	src/lasso/AssocPair.hpp \
	src/lasso/AssocPairSide.hpp \
	src/iRel_Lasso.cpp \
	src/lasso/iRel_Lasso.hpp \
	src/Lasso.cpp \
	src/lasso/Lasso.hpp

libiRel_la_include_HEADERS = \
	itaps/irel/iRel.h \
	itaps/irel/iRel_f.h \
	itaps/irel/iRel_protos.h \
	itaps/irel/iRel_FCDefs.h

# The directory in which to install headers
libiRel_la_includedir = $(includedir)
libiRel_la_LIBADD =

lib_LTLIBRARIES = libiRel.la

if ENABLE_IGEOM
    AM_CPPFLAGS += -DENABLE_IGEOM -I${IGEOM_DIR}/include
    libiRel_la_SOURCES += \
        src/GeomAssocPairSide.cpp \
        src/lasso/GeomAssocPairSide.hpp
    libiRel_la_LIBADD += \
        ${IGEOM_LTFLAGS} ${IGEOM_LDFLAGS} ${IGEOM_DIR}/lib/libiGeom.la
    LIBS += ${IGEOM_LIBS}
endif

if ENABLE_FBIGEOM
    AM_CPPFLAGS += -DENABLE_FBIGEOM -I${FBIGEOM_DIR}/include
    libiRel_la_SOURCES += \
        src/FBGeomAssocPairSide.cpp \
        src/lasso/FBGeomAssocPairSide.hpp
    libiRel_la_LIBADD += \
        ${FBIGEOM_LTFLAGS} ${FBIGEOM_LDFLAGS} ${FBIGEOM_DIR}/lib/libFBiGeomMOAB.la
    LIBS += ${FBIGEOM_LIBS}
endif

if ENABLE_IMESH
    AM_CPPFLAGS += -DENABLE_IMESH -I${IMESH_DIR}/include
    libiRel_la_SOURCES += \
        src/MeshAssocPairSide.cpp \
        src/lasso/MeshAssocPairSide.hpp
    libiRel_la_LIBADD += \
        ${IMESH_LTFLAGS} ${IMESH_LDFLAGS} ${IMESH_DIR}/lib/libiMesh.la
    LIBS += ${IMESH_LIBS}
endif

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = itaps/irel/iRel-Defs.inc

# By default, iRel-Defs.inc will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "IREL_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iRel-Defs.inc
	echo "IREL_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iRel-Defs.inc

CLEANFILES=itaps/irel/iRel_protos.h
FPFX=iRel
SEDEXPR = 's/^[[:space:]]*void[[:space:]][[:space:]]*$(FPFX)_\([a-z][_a-zA-Z0-9]*\)[[:space:]]*(.*$$/$(FPFX)_\1/p'
itaps/irel/iRel_protos.h: itaps/irel/iRel.h
	echo '#include "iRel_FCDefs.h"' > $@
	echo '#ifdef IREL_FC_FUNC_' >> $@
	echo >> $@
	for func in `$(SED) -n $(SEDEXPR) $<`; do \
	  lower=`echo $$func | tr '[:upper:]' '[:lower:]'`; \
          upper=`echo $$func | tr '[:lower:]' '[:upper:]'`; \
	  echo "#define $$func IREL_FC_FUNC_( $$lower, $$upper )" >> $@; \
        done
	echo >> $@
	echo "#endif" >> $@

itaps/irel/iRel_FCDefs.h: itaps/iBase_FCDefs.h
	./config.status itaps/irel/iRel_FCDefs.h

DISTCLEANFILES = itaps/irel/iRel_FCDefs.h
