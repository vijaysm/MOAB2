
AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io/mhdf/include \
               -I$(top_builddir)/ \
               -I$(top_builddir)/include \
               -I$(top_srcdir)/src/LocalDiscretization/ \
               -I$(top_srcdir)/tools/mbcoupler \
               -I$(top_srcdir)/src/moab/point_locator/lotte \
               -DIS_BUILDING_MB
            
LDADD = $(top_builddir)/src/libMOAB.la

check_PROGRAMS = 
noinst_PROGRAMS = point_location elem_eval_time

point_location_SOURCES = point_location.cpp
elem_eval_time_SOURCES = elem_eval_time.cpp

