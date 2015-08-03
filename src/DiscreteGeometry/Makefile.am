# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

# Things to build
noinst_LTLIBRARIES = libDGM.la
         
# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB -DSRCDIR=$(srcdir) \
               -I$(srcdir)/..  -I.. -I$(srcdir)/../moab -I$(srcdir)/moab \
	       -I$(srcdir)/../parallel

# The directory in which to install headers
libDGM_la_includedir = $(includedir)

# Conditional sources
MOAB_LOCALDISCR_SRCS = \
	Solvers.cpp 

MOAB_LOCALDISCR_HDRS = \
	moab/Solvers.hpp 

# The list of source files, and any header files that do not need to be installed
libDGM_la_SOURCES = \
   $(MOAB_LOCALDISCR_SRCS)

# The list of header files which are to be installed
nobase_libDGM_la_include_HEADERS = \
  $(MOAB_LOCALDISCR_HDRS)
