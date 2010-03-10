# Things to build
lib_LTLIBRARIES = libMOAB.la 

# Subdirectories to build
SUBDIRS = 
libMOAB_la_LIBADD = io/libmoabio.la
if PARALLEL
  SUBDIRS += parallel
  libMOAB_la_LIBADD += parallel/libMOABpar.la
  AM_CPPFLAGS += -I$(srcdir)/parallel
endif
SUBDIRS += io .

BUILT_SOURCES = MBCN_protos.h

# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I$(builddir) \
               -I$(srcdir)/io \
               -I$(srcdir)/io/moab \
               -I$(srcdir)/io/mhdf/include \
               -I$(builddir)/moab \
               -I$(srcdir)/moab \
               -I$(srcdir)/parallel/moab \
               -I$(builddir)/parallel/moab

# The directory in which to install headers
libMOAB_la_includedir = $(includedir)


# The list of source files, and any header files that do not need to be installed
libMOAB_la_SOURCES = \
  AEntityFactory.cpp \
  AEntityFactory.hpp \
  BitTagServer.cpp \
  BitTagServer.hpp \
  BSPTreePoly.cpp \
  DualTool.cpp \
  ElementSequence.hpp \
  EntitySequence.cpp \
  EntitySequence.hpp \
  FileOptions.cpp \
  FileOptions.hpp \
  GeomTopoTool.cpp \
  HigherOrderFactory.cpp \
  HomXform.cpp \
  MBAdaptiveKDTree.cpp \
  MBAffineXform.cpp \
  MBAffineXform.hpp \
  MBAxisBox.cpp \
  MBAxisBox.hpp \
  MBBSPTree.cpp \
  MBCN.cpp \
  MBCNArrays.hpp \
  MBCartVect.cpp \
  MBCartVect.hpp \
  MBMatrix3.cpp \
  MBMatrix3.hpp \
  MBCore.cpp \
  MBFactory.cpp \
  MBGeomUtil.cpp \
  MBInternals.hpp \
  MBMeshSet.cpp \
  MBMeshSet.hpp \
  MBOrientedBox.cpp \
  MBOrientedBox.hpp \
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
  SweptElementData.cpp \
  SweptElementData.hpp \
  SweptElementSeq.cpp \
  SweptElementSeq.hpp \
  SweptVertexData.cpp \
  SweptVertexData.hpp \
  TagCompare.hpp \
  TagInfo.cpp \
  TagInfo.hpp \
  TagServer.cpp \
  TagServer.hpp \
  TypeSequenceManager.cpp \
  TypeSequenceManager.hpp \
  UnstructuredElemSeq.cpp \
  UnstructuredElemSeq.hpp \
  VarLenTag.hpp \
  VertexSequence.hpp \
  VertexSequence.cpp 

# The list of header files which are to be installed
nobase_libMOAB_la_include_HEADERS = \
  moab/DualTool.hpp \
  moab/GeomTopoTool.hpp \
  moab/HigherOrderFactory.hpp \
  moab/HomXform.hpp \
  moab/MBAdaptiveKDTree.hpp \
  moab/MBBSPTree.hpp \
  moab/MBCN.hpp \
  moab/MBCore.hpp \
  moab/MBEntityType.h \
  moab/MBEntityHandle.h \
  moab/MBError.hpp \
  moab/MBForward.hpp \
  moab/MBGeomUtil.hpp \
  moab/MBInterface.hpp \
  moab/MBOrientedBoxTreeTool.hpp \
  moab/BSPTreePoly.hpp \
  moab/MBRange.hpp \
  moab/MBReadUtilIface.hpp \
  moab/MBReaderIface.hpp \
  moab/MBReaderWriterSet.hpp \
  moab/MBSkinner.hpp \
  moab/MBTagConventions.hpp \
  moab/MBTypes.h \
  moab/MBUnknownInterface.hpp \
  moab/MBUtil.hpp \
  moab/MBWriteUtilIface.hpp \
  moab/MBWriterIface.hpp \
  moab/MeshTopoUtil.hpp \
  moab/RangeMap.hpp \
  MBCN.h \
  MOAB_FCDefs.h \
  MBCN_protos.h 

MBCN_protos.h: MBCN.h $(top_srcdir)/itaps/mkprotos.sh
	$(top_srcdir)/itaps/mkprotos.sh MBCN MBCN $< $@ MOAB_FCDefs.h


MOAB_FCDefs.h: FCDefs.h
	cd .. && ./config.status src/MOAB_FCDefs.h

MOSTLYCLEANFILES = MOAB_FCDefs.h
