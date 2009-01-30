AUTOMAKE_OPTIONS = foreign
BUILT_SOURCES = iRel_protos.h

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I. ${IMESH_INCLUDES} ${IGEOM_INCLUDES}

DEFS += $(DEFINES)

LIBS = $(IMESH_CXX_LDFLAGS) $(IMESH_CXX_LIBS)

# stuff for c test program
bin_PROGRAMS = rtestc
rtestc_SOURCES = rtestc.c

LDADD = libiRel.la ${IGEOM_LIBS} ${IMESH_LIBS}

rtestc_DEPENDENCIES = libiRel.la 

libiRel_la_SOURCES = \
	iRel_Lasso.cpp \
	AssocPair.cpp \
	AssocPairC.cpp 
                   
libiRel_la_include_HEADERS = iRel.h \
                             iRel_f.h \
                             iRel_protos.h 

# The directory in which to install headers
libiRel_la_includedir = $(includedir)

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
	echo '#include "iBase_FCDefs.h"' > $@
	echo '#ifdef FC_FUNC_' >> $@
	echo >> $@
	for func in `$(SED) -n $(SEDEXPR) $<`; do \
	  lower=`echo $$func | tr '[:upper:]' '[:lower:]'`; \
          upper=`echo $$func | tr '[:lower:]' '[:upper:]'`; \
	  echo "#define $$func FC_FUNC_( $$lower, $$upper )" >> $@; \
        done
	echo >> $@
	echo "#endif" >> $@
