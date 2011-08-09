check_PROGRAMS = FindConnectF TagIterateC TagIterateF

AM_CPPFLAGS += -I$(top_srcdir)/itaps/imesh \
	-I$(top_srcdir)/itaps/ \
	-I$(top_builddir)/itaps/imesh

AM_FFLAGS = -fcray-pointer

TESTS = $(check_PROGRAMS)
LDADD = $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la ${MOAB_CXX_LINKFLAGS} ${MOAB_CXX_LIBS}
TESTDEPS = $(top_builddir)/itaps/imesh/libiMesh.la $(top_builddir)/src/libMOAB.la

FindConnectF_SOURCES = FindConnectF.F
FindConnectF_DEPENDENCIES = $(TESTDEPS)

TagIterateC_SOURCES = TagIterateC.c
TagIterateC_DEPENDENCIES = $(TESTDEPS)

TagIterateF_SOURCES = TagIterateF.F
TagIterateF_DEPENDENCIES = $(TESTDEPS)
