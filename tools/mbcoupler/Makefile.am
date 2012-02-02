AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/test \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps \
               -I$(top_builddir)/itaps/imesh 

lib_LTLIBRARIES = libmbcoupler.la
libmbcoupler_la_LIBADD = $(top_builddir)/src/libMOAB.la $(top_builddir)/itaps/imesh/libiMesh.la

LDADD = libmbcoupler.la 

libmbcoupler_la_SOURCES = \
   ElemUtil.cpp \
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
   Coupler.cpp

libmbcoupler_la_include_HEADERS = \
   ElemUtil.hpp \
   Coupler.hpp

libmbcoupler_la_includedir = $(includedir)

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

TESTS = findpt_test elem_util_test element_test
check_PROGRAMS = $(TESTS)
findpt_test_SOURCES = findpt_test.cpp
elem_util_test_SOURCES = ElemUtilTest.cpp
element_test_SOURCES = ElementTest.cpp

noinst_PROGRAMS = 
if PARALLEL
#  TESTS += mbcoupler_test
#  check_PROGRAMS += mbcoupler_test
  noinst_PROGRAMS += mbcoupler_test ssn_test
  mbcoupler_test_SOURCES = mbcoupler_test.cpp
  ssn_test_SOURCES = ssn_test.cpp
  AM_CPPFLAGS += -I$(srcdir)/parallel
endif

