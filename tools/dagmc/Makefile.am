AM_CPPFLAGS = -I$(top_srcdir)/src \
              -I$(top_builddir)/src \
	            -I$(top_srcdir)/tools

lib_LTLIBRARIES = libdagmc.la
libdagmc_la_SOURCES = DagMC.cpp
libdagmc_la_includedir = $(includedir)
libdagmc_la_include_HEADERS = DagMC.hpp

AM_CXXFLAGS =
if HAVE_CGM
  AM_CPPFLAGS += $(CGM_INCLUDES) $(CGM_CPPFLAGS) 
  AM_CXXFLAGS += $(CGM_CXXFLAGS)

if HAVE_NETCDF
  netcdf_progs = update_coords 
else
  netcdf_progs = 
endif

endif

if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

bin_PROGRAMS = dagmc_preproc
noinst_PROGRAMS = quads_to_tris $(netcdf_progs) 

CGM_CPPFLAGS = @CGM_CPPFLAGS@
CGM_LDFLAGS = @CGM_LDFLAGS@
CGM_LTFLAGS = @CGM_LTFLAGS@
CGM_LIBS = @CGM_LIBS@

LDADD = libdagmc.la $(top_builddir)/src/libMOAB.la   
if HAVE_CGM
 LDADD += $(CGM_DIR)/lib/libcgm.la 
endif

TESTS = 
if HAVE_HDF5
   TESTS+=test_geom 
endif
check_PROGRAMS = $(TESTS) pt_vol_test ray_fire_test

quads_to_tris_SOURCES = quads_to_tris_driver.cpp quads_to_tris.cpp
quads_to_tris_LDADD = $(top_builddir)/src/libMOAB.la

test_geom_SOURCES = test_geom.cc

update_coords_SOURCES = cub2h5m.cc quads_to_tris.cpp quads_to_tris.hpp

pt_vol_test_SOURCES = pt_vol_test.cc

ray_fire_test_SOURCES = ray_fire_test.cc

dagmc_preproc_SOURCES = dagmc_preproc.cpp obb_analysis.cpp dagmc_preproc.hpp

