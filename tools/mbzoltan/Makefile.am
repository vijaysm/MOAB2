AUTOMAKE_OPTIONS = foreign
DEFS = $(DEFINES) -DIS_BUILDING_MB
INCLUDES += -I$(top_srcdir) \
            -I$(top_srcdir)/parallel \
            -I$(top_builddir) \
            $(ZOLTAN_INC_FLAGS)
bin_PROGRAMS = mbzoltan
mbzoltan_SOURCES = MBZoltan.cpp MBZoltan.hpp main.cpp
mbzoltan_LDADD = $(top_builddir)/libMOAB.la $(ZOLTAN_LIB_FLAGS) $(ZOLTAN_LIBS)
mbzoltan_DEPENDENCIES = $(top_builddir)/libMOAB.la


