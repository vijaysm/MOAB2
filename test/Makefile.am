LDADD = $(top_builddir)/libiRel.la 

if ENABLE_IGEOM
 LDADD +=${IGEOM_DIR}/lib/libiGeom.la ${IGEOM_DIR}/lib/libcgm.la
endif

if ENABLE_IMESH
 LDADD += ${IMESH_DIR}/lib/libiMesh.la
endif

AM_CPPFLAGS = -DSRCDIR=$(srcdir) $(DEFINES) -I$(top_srcdir)/src -I$(top_srcdir)/src/lasso \
							-I$(top_builddir) -I$(top_srcdir)/itaps/irel -I$(top_builddir)/itaps/irel \
							-I${IGEOM_DIR}/include -I${IMESH_DIR}/include

TESTS = test_entset test_both test_status 
if ENABLE_FBIGEOM
  TESTS+= test_entsetFB
endif
check_PROGRAMS = ${TESTS}

# stuff for c test programs
test_entset_SOURCES = test_entset.cpp
test_entset_DEPENDENCIES = \
	$(LDADD)
test_entset_LDADD = $(LDADD)

test_both_SOURCES = test_both.cpp TestUtil.hpp
test_both_DEPENDENCIES = \
	$(LDADD)
test_both_LDADD = $(LDADD)

test_status_SOURCES = test_status.cpp TestUtil.hpp
test_status_DEPENDENCIES = \
	$(LDADD)
test_status_LDADD = $(LDADD)

if ENABLE_FBIGEOM
LDADD+=${FBIGEOM_DIR}/lib/libFBiGeomMOAB.la
test_entsetFB_SOURCES = test_entsetFB.cpp TestUtil.hpp
test_entsetFB_DEPENDENCIES = \
        $(LDADD)
endif
EXTRA_DIST = brick.cub
if ENABLE_FBIGEOM
   EXTRA_DIST+= shell.h5m shellQuad.h5m
endif 
