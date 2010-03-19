AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_builddir)/src/parallel

LDADD = $(top_builddir)/src/libMOAB.la

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

if ENABLE_mbdepth
  bin_PROGRAMS += mbdepth
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

SUBDIRS = $(mbcoupler_DIR) \
          $(mbzoltan_DIR) \
          $(mbperf_DIR) \
	  $(mcnpmit_DIR) \
          $(qvdual_DIR) \
          $(dagmc_DIR)

mbconvert_SOURCES = convert.cpp
mbconvert_DEPENDENCIES = $(LDADD)

mbgsets_SOURCES = gsets.cc
mbgsets_DEPENDENCIES = $(LDADD)

hexmodops_SOURCES = makeops.cpp
hexmodops_DEPENDENCIES = $(LDADD)

mbsize_SOURCES = measure.cpp measure.hpp size.cpp
mbsize_DEPENDENCIES = $(LDADD)

mbskin_SOURCES = skin.cpp
mbskin_DEPENDENCIES = $(LDADD)

mbdepth_SOURCES = depth.cpp
mbdepth_DEPENDENCIES = $(LDADD)

spheredecomp_SOURCES = SphereDecomp.hpp SphereDecomp.cpp sphere_decomp.cpp
spheredecomp_DEPENDENCIES = $(LDADD)

mbsurfplot_SOURCES = surfplot.cpp
mbsurfplot_DEPENDENCIES = $(LDADD)

mbtagprop_SOURCES = parse.cpp parse.hpp propagate_tags.cpp
mbtagprop_DEPENDENCIES = $(LDADD)
