AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB

@CGM_CONFIG_OPTIONS@

LDADD = $(top_builddir)/libMOAB.la $(ZOLTAN_LIB_FLAGS) $(ZOLTAN_LIBS)
if HAVE_CGM
  LDADD += $(CGM_LTFLAGS) 
endif


INCLUDES += -I$(top_srcdir) \
            -I$(top_srcdir)/parallel \
            -I$(top_builddir) \
            $(ZOLTAN_INC_FLAGS)
bin_PROGRAMS = mbzoltan
mbzoltan_SOURCES = MBZoltan.cpp MBZoltan.hpp main.cpp
mbzoltan_DEPENDENCIES = $(top_builddir)/libMOAB.la


