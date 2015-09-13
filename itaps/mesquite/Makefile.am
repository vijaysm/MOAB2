if HAVE_CGM
  IGEOM_LIB = libmsqigeom.la
else
  IGEOM_LIB = 
endif

if ENABLE_irel
  IREL_LIB = libmsqirel.la
else
  IREL_LIB = 
endif

if ENABLE_imesh
  IMESH_LIB = libmsqimesh.la
if HAVE_MPI
  IMESHP_LIB = libmsqimeshp.la
else
  IMESHP_LIB =
endif
else
  IMESH_LIB = 
  IMESHP_LIB = 
endif

lib_LTLIBRARIES = $(IMESH_LIB) $(IGEOM_LIB) $(IREL_LIB) $(IMESHP_LIB)

libmsqimesh_la_includedir = $(includedir)
libmsqimesh_la_SOURCES = MsqIMesh.cpp 
libmsqimesh_la_include_HEADERS = MsqIMesh.hpp 
libmsqimesh_la_LIBADD = $(top_builddir)/src/mesquite/libmesquite.la $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la

libmsqigeom_la_includedir = $(includedir)
libmsqigeom_la_SOURCES = MsqIGeom.cpp 
libmsqigeom_la_include_HEADERS = MsqIGeom.hpp 
libmsqigeom_la_LIBADD = $(top_builddir)/src/mesquite/libmesquite.la ${IGEOM_LTFLAGS} ${IGEOM_LDFLAGS} ${CGM_DIR}/lib/libiGeom.la

libmsqirel_la_includedir = $(includedir)
libmsqirel_la_SOURCES = MsqIRel.cpp 
libmsqirel_la_include_HEADERS = MsqIRel.hpp 
libmsqirel_la_LIBADD = libmsqigeom.la libmsqimesh.la $(top_builddir)/itaps/irel/libiRel.la $(top_builddir)/src/libMOAB.la

if HAVE_MPI
libmsqimeshp_la_includedir = $(includedir)
libmsqimeshp_la_SOURCES = MsqIMeshP.cpp 
libmsqimeshp_la_include_HEADERS = MsqIMeshP.hpp 
libmsqimeshp_la_LIBADD = libmsqimesh.la $(top_builddir)/src/mesquite/libmesquite.la $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la
endif

AM_CPPFLAGS = $(CPPFLAGS) \
              -I$(top_builddir)/src/mesquite/include \
              -I$(top_srcdir)/src/mesquite/include \
              -I$(top_srcdir)/src/mesquite/Mesh \
              -I$(top_srcdir)/src/mesquite/Misc \
              -I$(top_srcdir)/src/mesquite/Control \
              -I$(top_srcdir)/src/mesquite/MappingFunction \
              $(IGEOM_DEFINES) $(IGEOM_INCLUDES) \
              -I$(top_srcdir)/src \
              -I$(top_srcdir)/src/parallel \
              -I$(top_srcdir)/itaps \
              -I$(top_srcdir)/itaps/imesh \
              -I$(top_builddir)/itaps/imesh \
              -I$(top_srcdir)/itaps/fbigeom \
              -I$(top_srcdir)/itaps/irel \
              -I$(top_builddir)/itaps/irel

EXTRA_DIST = 

