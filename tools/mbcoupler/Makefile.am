AUTOMAKE_OPTIONS = foreign

@CGM_CONFIG_OPTIONS@

DEFS = $(DEFINES)
INCLUDES += -I$(top_srcdir) -I$(top_builddir) -I$(top_srcdir)/parallel

lib_LTLIBRARIES = libmbcoupler.la

libmbcoupler_la_SOURCES = \
   MBElemUtil.cpp \
   poly.c   \
   poly.h   \
   tensor.c \
   tensor.h \
   errmem.c \
   errmem.h \
   extrafindpt.h \
   findpt.c \
   findpt.h \
   minmax.h \
   types.h \
   MBCoupler.cpp

libmbcoupler_la_include_HEADERS = \
   MBElemUtil.hpp \
   MBCoupler.hpp

libmbcoupler_la_includedir = $(includedir)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

TESTS = findpt_test elem_util_test
check_PROGRAMS = $(TESTS)
findpt_test_SOURCES = findpt_test.cpp
findpt_test_LDADD = libmbcoupler.la $(top_builddir)/libMOAB.la
elem_util_test_SOURCES = ElemUtilTest.cpp
elem_util_test_LDADD = libmbcoupler.la $(top_builddir)/libMOAB.la

bin_PROGRAMS = 
if PARALLEL
#  TESTS += mbcoupler_test
#  check_PROGRAMS += mbcoupler_test
  bin_PROGRAMS += mbcoupler_test
  mbcoupler_test_SOURCES = mbcoupler_test.cpp
  mbcoupler_test_LDADD = libmbcoupler.la $(top_builddir)/libMOAB.la
  INCLUDES += -I$(srcdir)/parallel
endif

