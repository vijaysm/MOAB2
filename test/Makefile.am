LDADD = $(top_builddir)/libiRel.la

INCLUDES = -I$(top_srcdir) -I$(top_builddir) -I. -I${IGEOM_DIR}/include -I${IMESH_DIR}/include

TESTS = test_entset
check_PROGRAMS = ${TESTS}

# stuff for c test program
test_entset_SOURCES = test_entset.c

test_entset_DEPENDENCIES = \
	${IGEOM_DIR}/lib/libiGeom.la \
	${IMESH_DIR}/lib/libiMesh.la \
	$(LDADD)

test_entset_LDADD =	$(LDADD)

EXTRA_DIST = brick.cub