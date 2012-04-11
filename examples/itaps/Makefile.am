check_PROGRAMS = TagIterateC ListSetsNTagsCXX
if ENABLE_FORTRAN
check_PROGRAMS += FindConnectF TagIterateF ListSetsNTagsF90
endif

AM_CPPFLAGS += -I$(top_srcdir)/itaps/imesh \
	-I$(top_srcdir)/itaps/ \
	-I$(top_srcdir)/src/parallel/ \
	-I$(top_builddir)/itaps \
	-I$(top_builddir)/itaps/imesh

TESTS = $(check_PROGRAMS)
LDADD = $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}
TESTDEPS = $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la

FindConnectF_SOURCES = FindConnectF.F
FindConnectF_DEPENDENCIES = $(TESTDEPS)

TagIterateC_SOURCES = TagIterateC.cpp
TagIterateC_DEPENDENCIES = $(TESTDEPS)

TagIterateF_SOURCES = TagIterateF.F
TagIterateF_DEPENDENCIES = $(TESTDEPS)

ListSetsNTagsF90_SOURCES = ListSetsNTagsF90.F90
ListSetsNTagsF90_DEPENDENCIES = $(TESTDEPS)

ListSetsNTagsCXX_SOURCES = ListSetsNTagsCXX.cpp
ListSetsNTagsCXX_DEPENDENCIES = $(TESTDEPS)
