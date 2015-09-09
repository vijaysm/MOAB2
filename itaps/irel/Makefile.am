@IGEOM_CONFIG_OPTIONS@

AUTOMAKE_OPTIONS = foreign
BUILT_SOURCES = $(top_builddir)/itaps/irel/iRel_protos.h

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest/irel
AM_CPPFLAGS = -DSRCDIR=$(srcdir)

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS += -I$(top_builddir) $(IGEOM_DEFINES) $(IGEOM_INCLUDES) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/itaps -I$(top_builddir)/itaps/fbigeom \
               -I$(top_srcdir)/itaps/imesh -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/itaps/irel -I$(top_builddir)/itaps/irel
DEFS += 
LIBS =

libiRel_la_SOURCES = \
	ArrayManager.hpp \
	AssocPair.cpp \
	AssocPair.hpp \
	AssocPairSide.hpp \
	iRel_Lasso.cpp \
	iRel_Lasso.hpp \
	Lasso.cpp \
	Lasso.hpp

libiRel_la_include_HEADERS = \
	$(top_srcdir)/itaps/irel/iRel.h \
	$(top_srcdir)/itaps/irel/iRel_f.h \
	$(top_builddir)/itaps/irel/iRel_protos.h

# The directory in which to install headers
libiRel_la_includedir = $(includedir)
libiRel_la_LIBADD =

lib_LTLIBRARIES = libiRel.la

if HAVE_CGM
    AM_CPPFLAGS += -DENABLE_IGEOM -I${CGM_DIR}/include
    libiRel_la_SOURCES += \
        GeomAssocPairSide.cpp \
        GeomAssocPairSide.hpp
    libiRel_la_LIBADD += \
        ${IGEOM_LTFLAGS} ${IGEOM_LDFLAGS} ${CGM_DIR}/lib/libiGeom.la
    LIBS += ${CGM_LIBS}
endif

if ENABLE_fbigeom
    AM_CPPFLAGS += -DENABLE_FBIGEOM -I$(top_srcdir)/itaps/fbigeom
    libiRel_la_SOURCES += \
        FBGeomAssocPairSide.cpp \
        FBGeomAssocPairSide.hpp
    libiRel_la_LIBADD += $(top_builddir)/itaps/fbigeom/libFBiGeomMOAB.la
endif

if ENABLE_imesh
    AM_CPPFLAGS += -DENABLE_IMESH -I$(top_srcdir)/itaps/imesh
    libiRel_la_SOURCES += \
        MeshAssocPairSide.cpp \
        MeshAssocPairSide.hpp
    libiRel_la_LIBADD += $(top_builddir)/itaps/imesh/libiMesh.la
endif


TESTS = test_entset test_both test_status 
if ENABLE_fbigeom
  TESTS+= test_entsetFB
endif
check_PROGRAMS = ${TESTS}

# stuff for c test programs
test_entset_SOURCES = test_entset.cpp
test_entset_DEPENDENCIES = libiRel.la
test_entset_LDADD = libiRel.la
test_entset_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/test

test_both_SOURCES = test_both.cpp
test_both_DEPENDENCIES = libiRel.la
test_both_LDADD = libiRel.la
test_both_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/test

test_status_SOURCES = test_status.cpp
test_status_DEPENDENCIES = libiRel.la
test_status_LDADD = libiRel.la
test_status_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/test

if ENABLE_fbigeom
test_entsetFB_SOURCES = test_entsetFB.cpp
test_entsetFB_DEPENDENCIES = libiRel.la
test_entsetFB_LDADD = libiRel.la
test_entsetFB_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/test
endif

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iRel-Defs.inc

# By default, iRel-Defs.inc will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "IREL_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iRel-Defs.inc
	echo "IREL_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iRel-Defs.inc

CLEANFILES=$(top_builddir)/itaps/irel/iRel_protos.h

#itaps/irel/iRel_FCDefs.h: itaps/iBase_FCDefs.h
#	./config.status itaps/irel/iRel_FCDefs.h

# DISTCLEANFILES = $(top_srcdir)/itaps/irel/iRel_FCDefs.h

$(top_builddir)/itaps/irel/iRel_protos.h: iRel.h $(top_srcdir)/itaps/mkprotos.sh Makefile.am
	$(AM_V_GEN)$(top_srcdir)/itaps/mkprotos.sh iRel MOAB $< $@ MOAB_FCDefs.h

