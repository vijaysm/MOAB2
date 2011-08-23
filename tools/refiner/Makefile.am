MESHDIR = $(top_srcdir)/tools/refiner/

AM_CPPFLAGS += -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_srcdir)/test \
               -DMESHDIR=$(MESHDIR)


lib_LTLIBRARIES = libMOABrefiner.la
libMOABrefiner_la_LIBADD = $(top_builddir)/src/libMOAB.la

LDADD = libMOABrefiner.la 

libMOABrefiner_la_SOURCES = \
  EdgeSizeEvaluator.cpp \
  EdgeSizeEvaluator.hpp \
  EdgeSizeSimpleImplicit.cpp \
  EdgeSizeSimpleImplicit.hpp \
  EntityRefiner.cpp \
  EntityRefiner.hpp \
  MeshOutputFunctor.cpp \
  MeshOutputFunctor.hpp \
  MeshRefiner.cpp \
  MeshRefiner.hpp \
  ProcessSet.cpp \
  ProcessSet.hpp \
  RefinerTagManager.cpp \
  RefinerTagManager.hpp \
  SimplexTemplateRefiner.cpp \
  SimplexTemplateRefiner.hpp \
  SimplexTemplateTagAssigner.cpp \
  SimplexTemplateTagAssigner.hpp \
  SplitVertices.hpp \
  SplitVertices.cpp

libMOABrefiner_la_include_HEADERS =

libMOABrefiner_la_includedir = $(includedir)

EXTRA_DIST = fourVolsBare.cub

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)

TESTS = test_mesh_refiner
check_PROGRAMS = $(TESTS)
test_mesh_refiner_SOURCES = test_mesh_refiner.cpp
test_mesh_refiner_LDADD = libMOABrefiner.la
test_mesh_refiner_CPPFLAGS = $(AM_CPPFLAGS)
