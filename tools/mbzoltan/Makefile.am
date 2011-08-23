
AM_CPPFLAGS += -DIS_BUILDING_MB\
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_builddir)/src \
	       -I$(srcdir)/.. \
               $(ZOLTAN_INC_FLAGS)

if HAVE_CGM
AM_CPPFLAGS += $(CGM_CPPFLAGS)
endif

AM_LDFLAGS += $(ZOLTAN_LIB_FLAGS)
lib_LTLIBRARIES = libmbzoltan.la
libmbzoltan_la_LIBADD = $(top_builddir)/src/libMOAB.la $(ZOLTAN_LIBS)
libmbzoltan_la_includedir = $(includedir)
libmbzoltan_la_SOURCES = MBZoltan.cpp
libmbzoltan_la_include_HEADERS = MBZoltan.hpp

bin_PROGRAMS = mbpart
LDADD = libmbzoltan.la
mbpart_SOURCES = mbpart.cpp
mbpart_LDADD = libmbzoltan.la 

