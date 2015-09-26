MSQ_ITAPS_LIB = libmsqitaps.la
MSQ_ITAPS_HDRS =
MSQ_ITAPS_SRCS =
MSQ_LIBDEPS = $(top_builddir)/src/mesquite/libmesquite.la

if ENABLE_irel
  MSQ_ITAPS_SRCS += MsqIRel.cpp
  MSQ_ITAPS_HDRS += MsqIRel.hpp
  MSQ_LIBDEPS += $(top_builddir)/itaps/irel/libiRel.la
endif

if ENABLE_imesh
  MSQ_ITAPS_SRCS += MsqIMesh.cpp
  MSQ_ITAPS_HDRS += MsqIMesh.hpp
  MSQ_LIBDEPS += $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la
if HAVE_MPI
  MSQ_ITAPS_SRCS += MsqIMeshP.cpp
  MSQ_ITAPS_HDRS += MsqIMeshP.hpp
endif
endif

if HAVE_CGM
  MSQ_ITAPS_SRCS += MsqIGeom.cpp
  MSQ_ITAPS_HDRS += MsqIGeom.hpp
  MSQ_LIBDEPS += ${IGEOM_LTFLAGS} ${IGEOM_LDFLAGS} ${CGM_DIR}/lib/libiGeom.la
endif

lib_LTLIBRARIES = $(MSQ_ITAPS_LIB)

libmsqitaps_la_includedir = $(includedir)
libmsqitaps_la_SOURCES = $(MSQ_ITAPS_SRCS)
libmsqitaps_la_include_HEADERS = $(MSQ_ITAPS_HDRS) 
libmsqitaps_la_LIBADD = $(MSQ_LIBDEPS)

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

