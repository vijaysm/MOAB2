AUTOMAKE_OPTIONS = foreign

@IGEOM_CONFIG_OPTIONS@
@IMESH_CONFIG_OPTIONS@

INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I. ${iMesh_INCLUDES} ${iGeom_INCLUDES}

DEFS += $(DEFINES)

LIBS = $(iMesh_CXX_LDFLAGS) $(iMesh_CXX_LIBS)

# stuff for c test program
bin_PROGRAMS = rtestc
rtestc_SOURCES = rtestc.c

LDADD = libiRel.la ${iMesh_DIR}/lib/libiMesh.la ${iGeom_DIR}/lib/libiGeom.la ${iGeom_LIBS} ${iMesh_LIBS}

rtestc_DEPENDENCIES = libiRel.la ${iMesh_DIR}/lib/libiMesh.la ${iGeom_DIR}/lib/libiGeom.la

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
	echo "iRel_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iRel-Defs.inc
	echo "iRel_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/iRel-Defs.inc
