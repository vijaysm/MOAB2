
AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(top_builddir)/ \
               -I$(top_builddir)/include \
               -I$(top_srcdir)/src/LocalDiscretization/ \
               -I$(top_srcdir)/tools/mbcoupler \
               -I$(top_srcdir)/src/moab/point_locator/lotte \
               -I$(top_srcdir)/src/LocalDiscretization/moab/ \
               -DIS_BUILDING_MB

# might need parallel here because of timer
if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

LDADD = $(top_builddir)/src/libMOAB.la $(top_builddir)/tools/mbcoupler/libmbcoupler.la

check_PROGRAMS = point_location elem_eval_time tree_searching_perf sploc_searching_perf
noinst_PROGRAMS =

point_location_SOURCES = point_location.cpp
elem_eval_time_SOURCES = elem_eval_time.cpp
tree_searching_perf_SOURCES = tree_searching_perf.cpp
sploc_searching_perf_SOURCES = sploc_searching_perf.cpp
