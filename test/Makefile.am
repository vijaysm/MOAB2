SUBDIRS = . io dual obb perf
if HDF5_FILE
  SUBDIRS += h5file
endif
if PARALLEL
  SUBDIRS += parallel
endif
if OLD_HEADERS
  SUBDIRS += oldinc
endif

LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS += -DIS_BUILDING_MB \
               -DSRCDIR=$(srcdir) \
		-DMESHDIR=$(MESHDIR) \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src 

TESTS = moab_test \
        homxform_test \
        scdseq_test \
        test_adj \
        seq_man_test \
        geom_util_test \
        xform_test \
        obb_test \
	adaptive_kd_tree_tests \
	bsp_tree_test \
	file_options_test \
        kd_tree_test \
        var_len_test var_len_test_no_template \
	tag_test \
	mesh_set_test \
        mbcn_test \
        bsp_tree_poly_test
#                 merge_test \         # input files no longer exist?
#                 test_tag_server \    # fails

check_PROGRAMS = $(TESTS) \
		 kd_tree_tool \
		 kd_tree_time

MOSTLYCLEANFILES = mb_write.g \
                   mb_write2.g \
                   hex_mid_volume_nodes.g \
                   hex_mid_edge_face_vol_nodes.g \
                   hex_mid_edge_nodes.g \
                   tet_mid_edge_nodes.g \
                   tet_mid_face_nodes.g \
                   tet_mid_edge_face_nodes.g \
                   tri_mid_edge_face_nodes.g \
                   mb_stress_out.g \
                   merge_test.ncdf \
                   tree.h5m

# Tests and such

moab_test_SOURCES = MBTest.cpp
moab_test_DEPENDENCIES = \
	$(top_srcdir)/MeshFiles/unittest/mb_big_test.g \
	$(top_srcdir)/MeshFiles/unittest/cell1.gen \
	$(top_srcdir)/MeshFiles/unittest/cell2.gen \
	$(LDADD)

homxform_test_SOURCES = $(top_srcdir)/src/HomXform.cpp
homxform_test_CPPFLAGS = -DTEST $(AM_CPPFLAGS) $(CPPFLAGS)
homxform_test_DEPENDENCIES = $(LDADD)

scdseq_test_SOURCES = scdseq_test.cpp
scdseq_test_DEPENDENCIES = $(LDADD)

#scdseq_timing_SOURCES = scdseq_timing.cpp
#scdseq_timing_DEPENDENCIES = $(LDADD)

#merge_test_SOURCES = merge_test.cpp
#merge_test_DEPENDENCIES = $(LDADD)

#test_exo_SOURCES = test_exo.cpp
#test_exo_DEPENDENCIES = $(LDADD)

#test_tag_server_SOURCES = $(top_srcdir)/src/TagServer.cpp
#test_tag_server_CPPFLAGS = -DTEST $(CPPFLAGS)
#test_tag_server_DEPENDENCIES = $(LDADD)

test_adj_SOURCES = test_adj.cpp
test_adj_DEPENDENCIES = $(LDADD)

seq_man_test_SOURCES = TestUtil.hpp TestTypeSequenceManager.cpp 

mesh_set_test_SOURCES = Test_MBMeshSet.cpp TestUtil.hpp

geom_util_test_SOURCES = TestUtil.hpp GeomUtilTests.cpp 
geom_util_test_DEPENDENCIES = $(LDADD)

xform_test_SOURCES = $(top_srcdir)/src/AffineXform.cpp
xform_test_CPPFLAGS = -DTEST $(AM_CPPFLAGS) $(CPPFLAGS)
xform_test_DEPENDENCIES = $(LDADD)

obb_test_SOURCES = OBBTest.cpp
obb_test_DEPENDENCIES = $(LDADD)

adaptive_kd_tree_tests_SOURCES = adaptive_kd_tree_tests.cpp
adaptive_kd_tree_tests_DEPENDENCIES = $(LDADD)

kd_tree_tool_SOURCES = kd_tree_tool.cpp
kd_tree_tool_DEPENDENCIES = $(LDADD)

kd_tree_time_SOURCES = kd_tree_time.cpp
kd_tree_time_DEPENDENCIES = $(LDADD)

kd_tree_test_SOURCES = kd_tree_test.cpp

bsp_tree_test_SOURCES = bsp_tree_test.cpp

file_options_test_SOURCES = TestUtil.hpp $(top_srcdir)/src/FileOptions.cpp 
file_options_test_CPPFLAGS = -DTEST $(AM_CPPFLAGS) $(CPPFLAGS)
file_options_test_LDADD = 

var_len_test_SOURCES = TestUtil.hpp VarLenTagTest.cpp
var_len_test_no_template_SOURCES = $(var_len_test_SOURCES)
var_len_test_no_template_CPPFLAGS = -UTEMPLATE_SPECIALIZATION $(AM_CPPFLAGS) $(CPPFLAGS)

tag_test_SOURCES = TestUtil.hpp TagTest.cpp

mbcn_test_SOURCES = $(top_srcdir)/src/moab/CN.hpp \
                    $(top_srcdir)/src/CN.cpp \
                    mbcn_test.cc

bsp_tree_poly_test_SOURCES = bsp_tree_poly_test.cpp

$(MESHDIR)/mb_big_test.g: $(MESHDIR)/mb_big_test.g.gz
	$(ZCAT) $< > $@

$(MESHDIR)/cell1.gen: $(MESHDIR)/cell1.gen.gz
	$(ZCAT) $< > $@

$(MESHDIR)/cell2.gen: $(MESHDIR)/cell2.gen.gz
	$(ZCAT) $< > $@


# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
