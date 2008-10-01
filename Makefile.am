# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign
# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I m4

# Subdirectories to build
SUBDIRS = 
if HDF5_FILE
  SUBDIRS += mhdf
endif
if PARALLEL
  SUBDIRS += parallel
endif
SUBDIRS += . test tools doc doxygen 


# Things to build
lib_LTLIBRARIES = libMOAB.la
TESTS = moab_test \
        homxform_test \
        scdseq_test \
        tqdcfr \
        test_adj \
        seq_man_test \
        geom_util_test \
        xform_test \
        obb_test \
	vtk_test \
	adaptive_kd_tree_tests \
	bsp_tree_test \
	file_options_test \
        kd_tree_test \
        var_len_test var_len_test_no_template \
	tag_test \
	mesh_set_test \
	cub_file_test
#                 merge_test \         # input files no longer exist?
#                 test_tag_server \    # fails

check_PROGRAMS = $(TESTS) \
		 kd_tree_tool \
		 kd_tree_time

#noinst_PROGRAMS = scdseq_timing


# Some variables
DEFS = $(DEFINES) -DIS_BUILDING_MB
INCLUDES += -I$(top_builddir) 

# The directory in which to install headers
libMOAB_la_includedir = $(includedir)



# Conditional sources
MOAB_EXTRA_SRCS =
MOAB_EXTRA_HDRS =
if NETCDF_FILE
  MOAB_EXTRA_SRCS +=  ReadNCDF.cpp  ReadNCDF.hpp \
                     WriteNCDF.cpp WriteNCDF.hpp \
                     WriteSLAC.cpp WriteSLAC.hpp 
endif
libMOAB_la_LIBADD = $(NETCDF_LIBS) $(HDF5_LIBS)
if HDF5_FILE
  libMOAB_la_LIBADD += $(top_builddir)/mhdf/libmhdf.la
  INCLUDES += -I$(srcdir)/mhdf/include
  MOAB_EXTRA_SRCS += ReadHDF5.cpp ReadHDF5.hpp WriteHDF5.cpp WriteHDF5.hpp
endif

if PARALLEL
  libMOAB_la_LIBADD += $(top_builddir)/parallel/libMOABpar.la
  INCLUDES += -I$(srcdir)/parallel
endif

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = moab.make

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/moab.make
	echo "MOAB_INCLUDEDIR=${includedir}" >> $(DESTDIR)$(cfgdir)/moab.make


# The list of source files, and any header files that do not need to be installed
libMOAB_la_SOURCES = \
  AEntityFactory.cpp \
  AEntityFactory.hpp \
  DualTool.cpp \
  ElementSequence.hpp \
  EntitySequence.cpp \
  EntitySequence.hpp \
  ExoIIUtil.cpp \
  ExoIIUtil.hpp \
  FileOptions.cpp \
  FileOptions.hpp \
  FileTokenizer.cpp \
  FileTokenizer.hpp \
  GeomTopoTool.cpp \
  HigherOrderFactory.cpp \
  HomXform.cpp \
  MBAdaptiveKDTree.cpp \
  MBAffineXform.cpp \
  MBAffineXform.hpp \
  MBAlloc.cpp \
  MBAlloc.hpp \
  MBAxisBox.cpp \
  MBAxisBox.hpp \
  MBBits.cpp \
  MBBits.hpp \
  MBBSPTree.cpp \
  MBBSPTree.hpp \
  MBCN.cpp \
  MBCNArrays.hpp \
  MBCartVect.cpp \
  MBMatrix3.cpp \
  MBMatrix3.hpp \
  MBCore.cpp \
  MBFactory.cpp \
  MBGeomUtil.cpp \
  MBInternals.hpp \
  MBMem.hpp \
  MBMeshSet.cpp \
  MBMeshSet.hpp \
  MBOrientedBox.cpp \
  MBOrientedBoxTreeTool.cpp \
  MBRange.cpp \
  MBRangeSeqIntersectIter.cpp \
  MBRangeSeqIntersectIter.hpp \
  MBReadUtil.cpp \
  MBReadUtil.hpp \
  MBReaderWriterSet.cpp \
  MBSkinner.cpp \
  MBSysUtil.cpp \
  MBSysUtil.hpp \
  MBUtil.cpp \
  MBWriteUtil.cpp \
  MBWriteUtil.hpp \
  MeshSetSequence.cpp \
  MeshSetSequence.hpp \
  MeshTopoUtil.cpp \
  PolyElementSeq.cpp \
  PolyElementSeq.hpp \
  RangeMap.hpp \
  ReadGmsh.cpp \
  ReadGmsh.hpp \
  ReadSTL.cpp \
  ReadSTL.hpp \
  ReadVtk.cpp \
  ReadVtk.hpp \
  ReadIDEAS.cpp \
  ReadIDEAS.hpp \
  ScdElementData.cpp \
  ScdElementData.hpp \
  ScdVertexData.cpp \
  ScdVertexData.hpp \
  SequenceData.hpp \
  SequenceData.cpp \
  SequenceManager.cpp \
  SequenceManager.hpp \
  SparseTagSuperCollection.cpp \
  SparseTagSuperCollection.hpp \
  SparseTagCollection.cpp \
  SparseTagCollection.hpp \
  StructuredElementSeq.cpp \
  StructuredElementSeq.hpp \
  TagCompare.hpp \
  TagInfo.cpp \
  TagInfo.hpp \
  TagServer.cpp \
  TagServer.hpp \
  Tqdcfr.cpp \
  Tqdcfr.hpp \
  TypeSequenceManager.cpp \
  TypeSequenceManager.hpp \
  UnstructuredElemSeq.cpp \
  UnstructuredElemSeq.hpp \
  VarLenTag.hpp \
  VertexSequence.hpp \
  VertexSequence.cpp \
  VtkUtil.cpp \
  VtkUtil.hpp \
  WriteAns.cpp \
  WriteAns.hpp \
  WriteGMV.cpp \
  WriteGMV.hpp \
  WriteGmsh.cpp \
  WriteGmsh.hpp \
  WriteSTL.cpp \
  WriteSTL.hpp \
  WriteTemplate.cpp \
  WriteVtk.cpp \
  WriteVtk.hpp \
  $(MOAB_EXTRA_SRCS)

#  MBMem.cpp

# The list of header files which are to be installed
libMOAB_la_include_HEADERS = \
  DualTool.hpp \
  ExoIIInterface.hpp \
  GeomTopoTool.hpp \
  HigherOrderFactory.hpp \
  HomXform.hpp \
  MBAdaptiveKDTree.hpp \
  MBCN.hpp \
  MBCN_FCDefs.h \
  MBCN_protos.h \
  MBCartVect.hpp \
  MBCore.hpp \
  MBEntityType.h \
  MBEntityHandle.h \
  MBError.hpp \
  MBForward.hpp \
  MBGeomUtil.hpp \
  MBInterface.hpp \
  MBOrientedBox.hpp \
  MBOrientedBoxTreeTool.hpp \
  MBParallelConventions.h \
  MBRange.hpp \
  MBReadUtilIface.hpp \
  MBReaderIface.hpp \
  MBReaderWriterSet.hpp \
  MBSkinner.hpp \
  MBTagConventions.hpp \
  MBTypes.h \
  MBUnknownInterface.h \
  MBUtil.hpp \
  MBWriteUtilIface.hpp \
  MBWriterIface.hpp \
  MeshTopoUtil.hpp \
  ReadWriteDefines.h \
  WriteTemplate.hpp \
  $(MOAB_EXTRA_HDRS)

CPPFLAGS += -DSRCDIR=$(srcdir)

# Tests and such

moab_test_SOURCES = MBTest.cpp
moab_test_LDADD = $(top_builddir)/libMOAB.la
moab_test_DEPENDENCIES = test/mb_big_test.g test/cell1.gen test/cell2.gen $(top_builddir)/libMOAB.la

homxform_test_SOURCES = HomXform.cpp
homxform_test_LDADD = $(top_builddir)/libMOAB.la
homxform_test_CXXFLAGS = -DTEST $(CXXFLAGS)
homxform_test_DEPENDENCIES = $(top_builddir)/libMOAB.la

scdseq_test_SOURCES = scdseq_test.cpp
scdseq_test_LDADD = $(top_builddir)/libMOAB.la
scdseq_test_DEPENDENCIES = $(top_builddir)/libMOAB.la

#scdseq_timing_SOURCES = scdseq_timing.cpp
#scdseq_timing_LDADD = $(top_builddir)/libMOAB.la
#scdseq_timing_DEPENDENCIES = $(top_builddir)/libMOAB.la

tqdcfr_SOURCES = Tqdcfr.cpp
tqdcfr_LDADD = $(top_builddir)/libMOAB.la
tqdcfr_CXXFLAGS = -DTEST_TQDCFR $(CXXFLAGS)
tqdcfr_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_rms_SOURCES = test_rms.cpp
#test_rms_LDADD = $(top_builddir)/libMOAB.la
#test_rms_DEPENDENCIES = $(top_builddir)/libMOAB.la

#merge_test_SOURCES = merge_test.cpp
#merge_test_LDADD = $(top_builddir)/libMOAB.la
#merge_test_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_exo_SOURCES = test_exo.cpp
#test_exo_LDADD = $(top_builddir)/libMOAB.la
#test_exo_DEPENDENCIES = $(top_builddir)/libMOAB.la

#test_tag_server_SOURCES = TagServer.cpp
#test_tag_server_LDADD = $(top_builddir)/libMOAB.la
#test_tag_server_CXXFLAGS = -DTEST $(CXXFLAGS)
#test_tag_server_DEPENDENCIES = $(top_builddir)/libMOAB.la

test_adj_SOURCES = test_adj.cpp
test_adj_LDADD = $(top_builddir)/libMOAB.la
test_adj_DEPENDENCIES = $(top_builddir)/libMOAB.la

seq_man_test_SOURCES = TestUtil.hpp TestTypeSequenceManager.cpp 
seq_man_test_LDADD = $(top_builddir)/libMOAB.la

mesh_set_test_SOURCES = Test_MBMeshSet.cpp TestUtil.hpp
mesh_set_test_LDADD = $(top_builddir)/libMOAB.la

geom_util_test_SOURCES = TestUtil.hpp GeomUtilTests.cpp 
geom_util_test_LDADD = $(top_builddir)/libMOAB.la
geom_util_test_DEPENDENCIES = $(geom_util_test_LDADD)

xform_test_SOURCES = MBAffineXform.cpp
xform_test_CXXFLAGS = -DTEST $(CXXFLAGS)
xform_test_LDADD = $(top_builddir)/libMOAB.la
xform_test_DEPENDENCIES = $(xform_test_LDADD)

obb_test_SOURCES = OBBTest.cpp
obb_test_LDADD = $(top_builddir)/libMOAB.la
obb_test_DEPENDENCIES = $(obb_test_LDADD)

vtk_test_SOURCES = TestUtil.hpp VtkTest.cpp
vtk_test_LDADD = $(top_builddir)/libMOAB.la
vtk_test_DEPENDENCIES = $(vtk_test_LDADD)

adaptive_kd_tree_tests_SOURCES = adaptive_kd_tree_tests.cpp
adaptive_kd_tree_tests_LDADD = $(top_builddir)/libMOAB.la
adaptive_kd_tree_tests_DEPENDENCIES = $(adaptive_kd_tree_tests_LDADD)

kd_tree_tool_SOURCES = kd_tree_tool.cpp
kd_tree_tool_LDADD = $(top_builddir)/libMOAB.la
kd_tree_tool_DEPENDENCIES = $(kd_tree_tool_LDADD)

kd_tree_time_SOURCES = kd_tree_time.cpp
kd_tree_time_LDADD = $(top_builddir)/libMOAB.la
kd_tree_time_DEPENDENCIES = $(kd_tree_time_LDADD)

kd_tree_test_SOURCES = kd_tree_test.cpp
kd_tree_test_LDADD = $(top_builddir)/libMOAB.la

bsp_tree_test_SOURCES = bsp_tree_test.cpp
bsp_tree_test_LDADD = $(top_builddir)/libMOAB.la

file_options_test_SOURCES = TestUtil.hpp FileOptions.cpp 
file_options_test_CPPFLAGS = -DTEST

var_len_test_SOURCES = TestUtil.hpp VarLenTagTest.cpp
var_len_test_no_template_SOURCES = $(var_len_test_SOURCES)
var_len_test_no_template_CPPFLAGS = -UTEMPLATE_SPECIALIZATION

tag_test_SOURCES = TestUtil.hpp TagTest.cpp
tag_test_LDADD = $(top_builddir)/libMOAB.la
tag_test_DEPENDENCIES = $(tag_test_LDADD)

cub_file_test_SOURCES = TestUtil.hpp cub_file_test.cc
cub_file_test_LDADD = $(top_builddir)/libMOAB.la
cub_file_test_DEPENDENCIES = $(tag_test_LDADD)


# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = dumped_acis.sat tree.h5m

test/mb_big_test.g: $(top_srcdir)/test/mb_big_test.g.gz
	$(ZCAT) $(top_srcdir)/test/mb_big_test.g.gz > test/mb_big_test.g

test/cell1.gen: $(top_srcdir)/test/cell1.gen.gz
	$(ZCAT) $(top_srcdir)/test/cell1.gen.gz > test/cell1.gen

test/cell2.gen: $(top_srcdir)/test/cell2.gen.gz
	$(ZCAT) $(top_srcdir)/test/cell2.gen.gz > test/cell2.gen

