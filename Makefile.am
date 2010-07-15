@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

AUTOMAKE_OPTIONS = foreign
BUILT_SOURCES = iRel_protos.h

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I. -I${IMESH_DIR}/include -I${IGEOM_DIR}/include

DEFS += $(DEFINES)

LIBS = 

# stuff for c test program
rtestc_SOURCES = rtestc.c

LDADD = libiRel.la 

TESTS = rtestc
check_PROGRAMS = ${TESTS}

rtestc_DEPENDENCIES = libiRel.la ${IGEOM_DIR}/lib/libiGeom.la ${IMESH_DIR}/lib/libiMesh.la 
rtestc_LDADD = libiRel.la #${IGEOM_DIR}/lib/libiGeom.la ${IMESH_DIR}/lib/libiMesh.la 
#rtestc_LDFLAGS = ${IGEOM_LTFLAGS} ${IGEOM_LDFLAGS} ${IMESH_LTFLAGS} ${IMESH_LDFLAGS} 

libiRel_la_SOURCES = \
	ArrayManager.hpp \
	AssocPair.cpp \
	AssocPair.hpp \
	AssocPairC.cpp \
	AssocPairC.hpp \
	iRel_Lasso.cpp \
	iRel_Lasso.hpp \
	Lasso.cpp \
	Lasso.hpp

libiRel_la_include_HEADERS = \
	iRel.h \
	iRel_f.h \
	iRel_protos.h \
	iRel_FCDefs.h

# The directory in which to install headers
libiRel_la_includedir = $(includedir)
libiRel_la_LIBADD = ${IGEOM_LTFLAGS} ${IGEOM_LDFLAGS} ${IMESH_LTFLAGS} \
                    ${IMESH_LDFLAGS} ${IGEOM_DIR}/lib/libiGeom.la ${IMESH_DIR}/lib/libiMesh.la 

lib_LTLIBRARIES = libiRel.la

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

if USE_BABEL
  SIDL_SUBDIR = SIDL
else
  SIDL_SUBDIR =
endif
SUBDIRS = . $(SIDL_SUBDIR)

CLEANFILES=iRel_protos.h
FPFX=iRel
SEDEXPR = 's/^[[:space:]]*void[[:space:]][[:space:]]*$(FPFX)_\([a-z][_a-zA-Z0-9]*\)[[:space:]]*(.*$$/$(FPFX)_\1/p'
iRel_protos.h: iRel.h
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

iRel_FCDefs.h: iBase_FCDefs.h
	./config.status iRel_FCDefs.h

DISTCLEANFILES = iRel_FCDefs.h
