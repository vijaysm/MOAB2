AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/itaps \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_builddir)/src/parallel

LDADD = $(top_builddir)/src/libMOAB.la
if ENABLE_imesh
 LDADD+= $(top_builddir)/itaps/imesh/libiMesh.la
endif

if HAVE_CGM
  AM_CPPFLAGS += $(CGM_CPPFLAGS)
  LDADD+=  ${CGM_DIR}/lib/libcgm.la
endif

bin_PROGRAMS += mbmesquite
AM_CPPFLAGS += $(MESQUITE_INC_FLAGS)
AM_LDFLAGS += $(MESQUITE_LIB_FLAGS)
mbmesquite_SOURCES = mbmesquite.cpp
mbmesquite_LDADD = $(top_builddir)/src/libMOAB.la $(MESQUITE_LIBS)

libmbmesquite_la_SOURCES = \
   MsqIMesh.cpp MsqIMesh.hpp MsqIBase.cpp MsqIBase.hpp

if PARALLEL
    libmbmesquite_la_SOURCES += MsqIMeshP.cpp MsqIMeshP.hpp
endif

if HAVE_CGM
    libmbmesquite_la_SOURCES += MsqIGeom.hpp MsqIGeom.hpp
endif

lib_LTLIBRARIES = libmbmesquite.la
LDADD = libmbmesquite.la  $(top_builddir)/src/libMOAB.la

mbmesquite_SOURCES = mbmesquite.cpp
