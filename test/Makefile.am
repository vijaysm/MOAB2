# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects no-dependencies

SUBDIRS = . io dual obb perf
if HAVE_HDF5
  SUBDIRS += h5file
endif
if ENABLE_dagmc
  SUBDIRS += dagmc
endif
if PARALLEL
  SUBDIRS += parallel
endif
#if OLD_HEADERS
#  SUBDIRS += oldinc
#endif

LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS  = -DIS_BUILDING_MB \
               -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src  \
               -I$(top_srcdir)/src/LocalDiscretization \
               -I$(top_srcdir)/src/verdict \
               -I$(top_srcdir)/src/RefineMesh

if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

TESTS = range_test \
        moab_test \
        homxform_test \
        scdseq_test \
        scd_test_partn \
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
        bsp_tree_poly_test \
        reorder_test \
        test_prog_opt \
        coords_connect_iterate \
        elem_eval_test \
        spatial_locator_test \
        test_boundbox \
        adj_moab_test \
        uref_mesh_test \
        verdict_test \
        test_Matrix3

if HAVE_HDF5
  TESTS += mbfacet_test \
        gttool_test \
        cropvol_test \
        mergemesh_test \
        mbground_test \
        lloyd_smoother_test \
        mhdf_public_test
endif


#       merge_test \         # input files no longer exist?
#       test_tag_server \    # fails

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
                   merge_test.g \
                   mm_out.h5m \
                   tree.h5m 

# Tests and such

test_Matrix3_SOURCES = test_Matrix3.cpp
verdict_test_SOURCES = verdict_test.cpp
range_test_SOURCES = range_test.cpp
moab_test_SOURCES = MBTest.cpp
moab_test_DEPENDENCIES = \
	mb_big_test.g \
	cell1.gen \
	cell2.gen \
	$(LDADD)

homxform_test_SOURCES = homxform_test.cpp
scdseq_test_SOURCES = scdseq_test.cpp
scd_test_partn_SOURCES = scd_test_partn.cpp
#scdseq_timing_SOURCES = scdseq_timing.cpp
#merge_test_SOURCES = merge_test.cpp
#test_exo_SOURCES = test_exo.cpp
#test_tag_server_SOURCES = $(top_srcdir)/src/TagServer.cpp
#test_tag_server_CPPFLAGS = -DTEST $(CPPFLAGS)
test_adj_SOURCES = test_adj.cpp
seq_man_test_SOURCES = TestUtil.hpp TestTypeSequenceManager.cpp 
mesh_set_test_SOURCES = Test_MBMeshSet.cpp TestUtil.hpp
geom_util_test_SOURCES = TestUtil.hpp GeomUtilTests.cpp 
xform_test_SOURCES = affinexform_test.cpp
obb_test_SOURCES = OBBTest.cpp
adaptive_kd_tree_tests_SOURCES = adaptive_kd_tree_tests.cpp
kd_tree_tool_SOURCES = kd_tree_tool.cpp
kd_tree_time_SOURCES = kd_tree_time.cpp
kd_tree_test_SOURCES = kd_tree_test.cpp
bsp_tree_test_SOURCES = bsp_tree_test.cpp
reorder_test_SOURCES = reorder_test.cpp
elem_eval_test_SOURCES = elem_eval_test.cpp
file_options_test_SOURCES = TestUtil.hpp fileopts_test.cpp
var_len_test_SOURCES = TestUtil.hpp VarLenTagTest.cpp
var_len_test_no_template_SOURCES = $(var_len_test_SOURCES)
var_len_test_no_template_CPPFLAGS = -UTEMPLATE_SPECIALIZATION $(AM_CPPFLAGS) $(CPPFLAGS)
tag_test_SOURCES = TestUtil.hpp TagTest.cpp
spatial_locator_test_SOURCES = TestUtil.hpp spatial_locator_test.cpp
mergemesh_test_SOURCES = mergemesh_test.cpp

#mbcn_test_SOURCES = $(top_srcdir)/src/moab/CN.hpp \
#                    $(top_srcdir)/src/CN.cpp \
#                    mbcn_test.cc
mbcn_test_SOURCES = mbcn_test.cc

bsp_tree_poly_test_SOURCES = bsp_tree_poly_test.cpp
mbfacet_test_SOURCES = mbfacet_test.cpp
mbground_test_SOURCES = mbground_test.cpp
gttool_test_SOURCES = gttool_test.cpp
cropvol_test_SOURCES = crop_vol_test.cpp

test_prog_opt_SOURCES = test_prog_opt.cpp
test_prog_opt_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)
moab_test_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)
kd_tree_test_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)

coords_connect_iterate_SOURCES = coords_connect_iterate.cpp
coords_connect_iterate_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)

test_boundbox_SOURCES = test_boundbox.cpp
lloyd_smoother_test_SOURCES = lloyd_smoother_test.cpp

adj_moab_test_SOURCES = adj_moab_test.cpp
uref_mesh_test_SOURCES = urefine_mesh_test.cpp

if PARALLEL
moab_test_CPPFLAGS += -I$(top_srcdir)/src/parallel
kd_tree_test_CPPFLAGS += -I$(top_srcdir)/src/parallel
test_prog_opt_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

EXTRA_DIST = TestUtil.hpp TestRunner.hpp

mb_big_test.g: $(MESHDIR)/mb_big_test.g.gz
	$(AM_V_GEN)$(ZCAT) $< > $@

cell1.gen: $(MESHDIR)/cell1.gen.gz
	$(AM_V_GEN)$(ZCAT) $< > $@

cell2.gen: $(MESHDIR)/cell2.gen.gz
	$(AM_V_GEN)$(ZCAT) $< > $@


CLEANFILES = mb_big_test.g \
             cell1.gen \
             cell2.gen

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
