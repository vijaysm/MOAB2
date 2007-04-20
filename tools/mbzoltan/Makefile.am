AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB
INCLUDES += -I$(top_srcdir) -I$(top_builddir) $(ZOLTAN_INCLUDES)
bin_PROGRAMS = mbzoltan
mbzoltan_SOURCES = MBZoltan.cpp MBZoltan.hpp main.cpp
LDADD = $(top_builddir)/libMOAB.la $(ZOLTAN_LIBS)
mbzoltan_DEPENDENCIES = $(top_builddir)/libMOAB.la


