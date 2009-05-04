AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir)
LDADD = $(top_builddir)/libMOAB.la

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
             README.spheredecomp \
             mbzoltan/README


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

if ENABLE_mbtagprop
  bin_PROGRAMS += mbtagprop
  doc_DATA += README.mbtagprop
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

if ENABLE_mbperf
  mbperf_DIR = mbperf
else
  mbperf_DIR =
endif

if ENABLE_mbchaco
  mbchaco_DIR = mbchaco
else
  mbchaco_DIR=
endif

if ENABLE_mbcoupler
  mbcoupler_DIR = mbcoupler
else
  mbcoupler_DIR=
endif

if ENABLE_mbzoltan
  mbzoltan_DIR = mbzoltan
else
  mbzoltan_DIR=
endif

if ENABLE_qvdual
  qvdual_DIR = qvdual
else
  qvdual_DIR=
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

if ENABLE_imesh
  imesh_DIR = iMesh
else
  imesh_DIR = 
endif 

SUBDIRS = $(imesh_DIR) \
          $(mbchaco_DIR) \
          $(mbcoupler_DIR) \
          $(mbzoltan_DIR) \
          $(mbperf_DIR) \
	  $(mcnpmit_DIR) \
          $(qvdual_DIR) \
          $(dagmc_DIR)

mbconvert_SOURCES = convert.cpp
mbconvert_DEPENDENCIES = $(top_builddir)/libMOAB.la

mbgsets_SOURCES = gsets.cc
mbgsets_DEPENDENCIES = $(top_builddir)/libMOAB.la

hexmodops_SOURCES = makeops.cpp
hexmodops_DEPENDENCIES = $(top_builddir)/libMOAB.la

mbsize_SOURCES = measure.cpp measure.hpp size.cpp
mbsize_DEPENDENCIES = $(top_builddir)/libMOAB.la

mbskin_SOURCES = skin.cpp
mbskin_DEPENDENCIES = $(top_builddir)/libMOAB.la

spheredecomp_SOURCES = SphereDecomp.hpp SphereDecomp.cpp main.cpp
spheredecomp_DEPENDENCIES = $(top_builddir)/libMOAB.la

mbsurfplot_SOURCES = surfplot.cpp
mbsurfplot_DEPENDENCIES = $(top_builddir)/libMOAB.la

mbtagprop_SOURCES = parse.cpp parse.hpp propagate_tags.cpp
mbtagprop_DEPENDENCIES = $(top_builddir)/libMOAB.la
