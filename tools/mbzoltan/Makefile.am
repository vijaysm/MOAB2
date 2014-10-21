
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
bin_PROGRAMS = mbpart
mbpart_SOURCES = mbpart.cpp
mbpart_LDADD = $(top_builddir)/src/libMOAB.la $(ZOLTAN_LIBS)

