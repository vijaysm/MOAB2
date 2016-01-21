AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS =  -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/itaps \
               -I$(top_builddir)/itaps/imesh

LDADD = $(top_builddir)/src/libMOAB.la
if ENABLE_imesh
 LDADD+= $(top_builddir)/itaps/imesh/libiMesh.la
endif

if HAVE_CGM
  AM_CPPFLAGS += $(CGM_CPPFLAGS)
  LDADD+=  ${CGM_DIR}/lib/libcgm.la
endif 

# For old (pre 1.10.x) versions of Automake
docdir = @docdir@

bin_PROGRAMS =
dist_man1_MANS = 
doc_DATA = 

EXTRA_DIST = README.hexmodops \
             README.mbconvert \
             README.mbsize \
             README.mbskin \
             README.mbsurfplot \
             README.mbtagprop \
             README.spheredecomp

if ENABLE_mbconvert
  bin_PROGRAMS += mbconvert
  dist_man1_MANS += mbconvert.man
  doc_DATA += README.mbconvert
endif

if ENABLE_hexmodops
  bin_PROGRAMS += hexmodops
  doc_DATA += README.hexmodops
endif

if ENABLE_mbsize
  bin_PROGRAMS += mbsize
  doc_DATA += README.mbsize
endif

if ENABLE_mbskin
  bin_PROGRAMS += mbskin
  doc_DATA += README.mbskin
endif

if ENABLE_mbdepth
  bin_PROGRAMS += mbdepth
endif

if ENABLE_mbtagprop
  bin_PROGRAMS += mbtagprop
  doc_DATA += README.mbtagprop
endif

if ENABLE_mbmem
  bin_PROGRAMS += mbmem
endif

if ENABLE_spheredecomp
  bin_PROGRAMS += spheredecomp
  doc_DATA += README.spheredecomp
endif

if ENABLE_mbsurfplot
  bin_PROGRAMS += mbsurfplot
  doc_DATA += README.mbsurfplot
endif

if ENABLE_gsets
  bin_PROGRAMS += mbgsets
endif 

if ENABLE_refiner
  refiner_DIR = refiner
else
  refiner_DIR=
endif

if ENABLE_mbcoupler
  mbcoupler_DIR = mbcoupler
else
  mbcoupler_DIR=
endif

if ENABLE_mbpart
  bin_PROGRAMS += mbpart
  AM_CPPFLAGS += $(ZOLTAN_INC_FLAGS)
  AM_LDFLAGS += $(ZOLTAN_LIB_FLAGS)
  mbpart_SOURCES = mbpart.cpp
  mbpart_LDADD = $(top_builddir)/src/libMOAB.la $(ZOLTAN_LIBS)
if HAVE_CGM
  mbpart_LDADD += $(CGM_LIBS)
endif 
else
if ENABLE_metis
  bin_PROGRAMS += mbpart
  AM_CPPFLAGS += $(METIS_IN_FLAGS)
  AM_LDFLAGS += $(METIS_LIB_FLAGS)
  mbpart_SOURCES = mbpart.cpp
  mbpart_LDADD = $(top_builddir)/src/libMOAB.la $(METIS_LIBS)
endif
endif

if ENABLE_mcnpmit
  mcnpmit_DIR = mcnpmit
else
  mcnpmit_DIR =
endif

if ENABLE_dagmc
  dagmc_DIR = dagmc
else
  dagmc_DIR = 
endif 

if PARALLEL
  parread = parread
else
  parread =
endif

if ENABLE_mbcslam
  mbcslam_DIR = mbcslam
else
  mbcslam_DIR = 
endif

if ENABLE_mbquality
  AM_CPPFLAGS += -I$(top_srcdir)/src/verdict
  bin_PROGRAMS += mbquality
endif

if ENABLE_mbmeshtogeom
if HAVE_CGM
  bin_PROGRAMS += mbmeshtogeom
endif
endif

# The directory in which to install headers

TEST_DEPS =
if ENABLE_imesh
    AM_CPPFLAGS += -DENABLE_IMESH -I$(top_srcdir)/itaps/imesh
if HAVE_MPI
    AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif
    TEST_DEPS += $(top_builddir)/itaps/imesh/libiMesh.la
endif

if ENABLE_fbigeom
    AM_CPPFLAGS += -DENABLE_FBIGEOM -I$(top_srcdir)/itaps/fbigeom
    TEST_DEPS += $(top_builddir)/itaps/fbigeom/libFBiGeomMOAB.la
endif
TEST_DEPS += $(top_builddir)/src/libMOAB.la

if HAVE_CGM
    AM_CPPFLAGS += -DENABLE_IGEOM -I${CGM_DIR}/include
    LIBS += ${CGM_LIBS}
    TEST_DEPS += ${CGM_DIR}/lib/libiGeom.la
endif
TEST_DEPS += ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}

if ENABLE_mbumr
  AM_CPPFLAGS += -I$(top_srcdir)/src/RefineMesh
  bin_PROGRAMS += mbumr
endif

if ENABLE_mbmerge
  bin_PROGRAMS += mbmerge
endif

noinst_PROGRAMS = $(parread)

SUBDIRS = . \
          $(refiner_DIR) \
          $(mbcoupler_DIR) \
          $(mbcslam_DIR) \
          $(mbzoltan_DIR) \
          $(mcnpmit_DIR) \
          $(vtkMOABReader_DIR) \
          $(dagmc_DIR) 

AM_CPPFLAGS += -I$(top_builddir) $(IGEOM_DEFINES) $(IGEOM_INCLUDES) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/itaps -I$(top_builddir)/itaps/fbigeom \
               -I$(top_srcdir)/itaps/imesh -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/itaps/irel -I$(top_builddir)/itaps/irel

mbconvert_SOURCES = convert.cpp
mbmerge_SOURCES = merge.cpp
mbgsets_SOURCES = gsets.cc
hexmodops_SOURCES = makeops.cpp
mbsize_SOURCES = measure.cpp measure.hpp size.cpp
mbskin_SOURCES = skin.cpp
mbdepth_SOURCES = depth.cpp
spheredecomp_SOURCES = SphereDecomp.hpp SphereDecomp.cpp sphere_decomp.cpp
mbsurfplot_SOURCES = surfplot.cpp
mbtagprop_SOURCES = parse.cpp parse.hpp propagate_tags.cpp
mbmem_SOURCES = mbmem.cpp
parread_SOURCES = parread.cpp
mbmeshtogeom_SOURCES = meshtogeom.cpp
mbmeshtogeom_LDADD = $(TEST_DEPS)
mbquality_SOURCES =  quality.cpp
mbumr_SOURCES = umr.cpp

