AUTOMAKE_OPTIONS = foreign

if USE_BABEL
  babel_DIR = SIDL
else
  babel_DIR = 
endif 

SUBDIRS = . $(babel_DIR)

INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I.

DEFS += $(DEFINES)

# stuff for c test program
bin_PROGRAMS = testc_cbind
testc_cbind_SOURCES = testc_cbind.c

LDADD = libiMesh.la $(top_builddir)/libMOAB.la
testc_cbind_DEPENDENCIES = libiMesh.la $(top_builddir)/libMOAB.la

libiMesh_la_SOURCES = \
	iMesh_MOAB.cpp
                   
libMOAB_la_include_HEADERS = iMesh.h \
                             iMesh_f.h \
                             iMesh_protos.h \
                             iBase.h \
                             iBase_f.h \
			     iBase_FCDefs.h 

lib_LTLIBRARIES = libiMesh.la

libMOAB_la_includedir = $(includedir)

libiMesh_la_LIBADD = $(top_builddir)/libMOAB.la

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iMesh-Defs.inc

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "iMesh_DIR=${DESTDIR}${cfgdir}/.." >> $(DESTDIR)$(cfgdir)/iMesh-SIDL-Defs.inc
	echo "iMesh_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iMesh-Defs.inc
	echo "iMesh_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iMesh-Defs.inc
	echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iMesh-Defs.inc

