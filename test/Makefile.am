LDADD = $(top_builddir)/libiRel.la

AM_CPPFLAGS = -DSRCDIR=$(srcdir) $(DEFINES)

INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I. -I${IGEOM_DIR}/include -I${IMESH_DIR}/include

TESTS = test_entset test_both test_status 
if ENABLE_FBIGEOM
  TESTS+= test_entsetFB
endif
check_PROGRAMS = ${TESTS}

# stuff for c test programs
test_entset_SOURCES = test_entset.c
test_entset_DEPENDENCIES = \
	${IGEOM_DIR}/lib/libiGeom.la \
	${IMESH_DIR}/lib/libiMesh.la \
	$(LDADD)
test_entset_LDADD = $(LDADD)

test_both_SOURCES = test_both.cpp TestUtil.hpp
test_both_DEPENDENCIES = \
	${IGEOM_DIR}/lib/libiGeom.la \
	${IMESH_DIR}/lib/libiMesh.la \
	$(LDADD)
test_both_LDADD = $(LDADD)

test_status_SOURCES = test_status.cpp TestUtil.hpp
test_status_DEPENDENCIES = \
	${IGEOM_DIR}/lib/libiGeom.la \
	${IMESH_DIR}/lib/libiMesh.la \
	$(LDADD)
test_status_LDADD = $(LDADD)

if ENABLE_FBIGEOM
test_entsetFB_SOURCES = test_entsetFB.c TestUtil.hpp
test_entsetFB_DEPENDENCIES = \
        ${FBIGEOM_DIR}/lib/libFBiGeomMOAB.la \
        ${IMESH_DIR}/lib/libiMesh.la \
        $(LDADD)
endif
EXTRA_DIST = brick.cub
if ENABLE_FBIGEOM
   EXTRA_DIST+= shell.h5m shellQuad.h5m
endif 
