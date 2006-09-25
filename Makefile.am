# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign



# Subdirectories to build
SUBDIRS = 
if HDF5_FILE
  SUBDIRS += mhdf
endif
SUBDIRS += . test tools doc doxygen


# Things to build
lib_LTLIBRARIES = libMOAB.la
check_PROGRAMS = moab_test \
                 homxform_test \
                 scdseq_test \
                 tqdcfr \
                 test_adj \
                 seq_man_test \
		 internals_test \
                 geom_util_test \
                 xform_test \
                 obb_test
#                 merge_test \         # input files no longer exist?
#                 test_tag_server \    # fails
TESTS = $(check_PROGRAMS)
         
#noinst_PROGRAMS = scdseq_timing


# Some variables
DEFS = $(DEFINES) -DIS_BUILDING_MB


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
if PARALLEL
  MOAB_PARALLEL_SRCS = WriteHDF5Parallel.cpp
  MOAB_PARALLEL_HDRS = WriteHDF5.hpp WriteHDF5Parallel.hpp
  PARALLEL_TESTS = test/h5file/parallel
else
  MOAB_PARALLEL_SRCS = WriteHDF5.hpp
  MOAB_PARALLEL_HDRS =
  PARALLEL_TESTS =
endif
if HDF5_FILE
  libMOAB_la_LIBADD = $(top_builddir)/mhdf/libmhdf.la
  INCLUDES += -I$(srcdir)/mhdf/include
  MOAB_EXTRA_SRCS += ReadHDF5.cpp ReadHDF5.hpp WriteHDF5.cpp $(MOAB_PARALLEL_SRCS)
  MOAB_EXTRA_HDRS += $(MOAB_PARALLEL_HDRS)
  TESTS += $(PARALLEL_TESTS)
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
  DenseTagCollections.cpp \
  DenseTagCollections.hpp \
  DualTool.cpp \
  EntitySequence.cpp \
  EntitySequence.hpp \
  EntitySequenceManager.cpp \
  EntitySequenceManager.hpp \
  ExoIIDefines.hpp \
  ExoIIUtil.cpp \
  ExoIIUtil.hpp \
  FileTokenizer.cpp \
  FileTokenizer.hpp \
  GeomTopoTool.cpp \
  HigherOrderFactory.cpp \
  HomXform.cpp \
  MBAffineXform.cpp \
  MBAffineXform.hpp \
  MBAlloc.cpp \
  MBAlloc.hpp \
  MBAxisBox.cpp \
  MBAxisBox.hpp \
  MBBits.cpp \
  MBBits.hpp \
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
  MBParallelComm.cpp \
  MBRange.cpp \
  MBRangeSeqIntersectIter.cpp \
  MBRangeSeqIntersectIter.hpp \
  MBReadUtil.cpp \
  MBReadUtil.hpp \
  MBReaderWriterSet.cpp \
  MBSkinner.cpp \
  MBUtil.cpp \
  MBWriteUtil.cpp \
  MBWriteUtil.hpp \
  MeshSetManager.cpp \
  MeshSetManager.hpp \
  MeshTopoUtil.cpp \
  PolyEntitySequence.cpp \
  PolyEntitySequence.hpp \
  ReadGmsh.cpp \
  ReadGmsh.hpp \
  ReadSTL.cpp \
  ReadSTL.hpp \
  ReadVtk.cpp \
  ReadVtk.hpp \
  ScdElementSeq.cpp \
  ScdElementSeq.hpp \
  ScdVertexSeq.cpp \
  ScdVertexSeq.hpp \
  SparseTagCollections.cpp \
  SparseTagCollections.hpp \
  TagServer.cpp \
  TagServer.hpp \
  Tqdcfr.cpp \
  Tqdcfr.hpp \
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
  MBCN.hpp \
  MBCartVect.hpp \
  MBCore.hpp \
  MBEntityType.h \
  MBError.hpp \
  MBGeomUtil.hpp \
  MBInterface.hpp \
  MBOrientedBox.hpp \
  MBOrientedBoxTreeTool.hpp \
  MBParallelComm.hpp \
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


# Tests and such

moab_test_SOURCES = MBTest.cpp
moab_test_LDADD = $(top_builddir)/libMOAB.la
moab_test_DEPENDENCIES = test/mb_big_test.g test/cell1.gen test/cell2.gen $(top_builddir)/libMOAB.la

internals_test_SOURCES = internals_test.cpp
internals_test_LDADD = $(top_builddir)/libMOAB.la

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

seq_man_test_SOURCES = EntitySequenceManager.cpp 
seq_man_test_LDADD = $(top_builddir)/libMOAB.la
seq_man_test_CXXFLAGS = -DTEST $(CXXFLAGS)

geom_util_test_SOURCES = GeomUtilTests.cpp 
geom_util_test_LDADD = $(top_builddir)/libMOAB.la
geom_util_test_DEPENDENCIES = $(geom_util_test_LDADD)

xform_test_SOURCES = MBAffineXform.cpp
xform_test_CXXFLAGS = -DTEST $(CXXFLAGS)
xform_test_LDADD = $(top_builddir)/libMOAB.la
xform_test_DEPENDENCIES = $(xform_test_LDADD)

obb_test_SOURCES = OBBTest.cpp
obb_test_LDADD = $(top_builddir)/libMOAB.la
obb_test_DEPENDENCIES = $(obb_test_LDADD)

# Other files to clean up (e.g. output from tests)
MOSTLYCLEANFILES = dumped_acis.sat

test/mb_big_test.g: $(top_srcdir)/test/mb_big_test.g.gz
	$(ZCAT) $(top_srcdir)/test/mb_big_test.g.gz > test/mb_big_test.g

test/cell1.gen: $(top_srcdir)/test/cell1.gen.gz
	$(ZCAT) $(top_srcdir)/test/cell1.gen.gz > test/cell1.gen

test/cell2.gen: $(top_srcdir)/test/cell2.gen.gz
	$(ZCAT) $(top_srcdir)/test/cell2.gen.gz > test/cell2.gen

