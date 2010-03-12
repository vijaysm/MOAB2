LDADD = $(top_builddir)/src/libMOAB.la $(ZOLTAN_LIB_FLAGS) $(ZOLTAN_LIBS)

AM_CPPFLAGS += -DIS_BUILDING_MB\
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_builddir)/src \
               $(ZOLTAN_INC_FLAGS)
bin_PROGRAMS = mbzoltan
mbzoltan_SOURCES = MBZoltan.cpp MBZoltan.hpp main.cpp
mbzoltan_DEPENDENCIES = $(top_builddir)/src/libMOAB.la


