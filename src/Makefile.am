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
if OLD_HEADERS
  SUBDIRS += oldinc
endif

BUILT_SOURCES = MBCN_protos.h

# Some variables
AM_CPPFLAGS += -DIS_BUILDING_MB \
               -I. \
               -I$(srcdir)/io \
               -I$(srcdir)/io/mhdf/include \
               -I$(srcdir)/parallel

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
  Error.hpp \
  FileOptions.cpp \
  FileOptions.hpp \
  GeomTopoTool.cpp \
  HigherOrderFactory.cpp \
  HomXform.cpp \
  AdaptiveKDTree.cpp \
  AffineXform.cpp \
  AffineXform.hpp \
  AxisBox.cpp \
  AxisBox.hpp \
  BSPTree.cpp \
  CN.cpp \
  MBCNArrays.hpp \
  CartVect.cpp \
  Matrix3.cpp \
  Matrix3.hpp \
  Core.cpp \
  Factory.cpp \
  GeomUtil.cpp \
  Internals.hpp \
  MeshSet.cpp \
  MeshSet.hpp \
  OrientedBox.cpp \
  OrientedBox.hpp \
  OrientedBoxTreeTool.cpp \
  Range.cpp \
  RangeSeqIntersectIter.cpp \
  RangeSeqIntersectIter.hpp \
  ReadUtil.cpp \
  ReadUtil.hpp \
  ReaderWriterSet.cpp \
  Skinner.cpp \
  SysUtil.cpp \
  SysUtil.hpp \
  Util.cpp \
  WriteUtil.cpp \
  WriteUtil.hpp \
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
  moab/AdaptiveKDTree.hpp \
  moab/BSPTree.hpp \
  moab/CN.hpp \
  moab/CartVect.hpp \
  moab/Core.hpp \
  moab/EntityType.hpp \
  moab/EntityHandle.hpp \
  moab/Forward.hpp \
  moab/GeomUtil.hpp \
  moab/Interface.hpp \
  moab/OrientedBoxTreeTool.hpp \
  moab/BSPTreePoly.hpp \
  moab/Range.hpp \
  moab/ReadUtilIface.hpp \
  moab/ReaderIface.hpp \
  moab/ReaderWriterSet.hpp \
  moab/Skinner.hpp \
  moab/MBTagConventions.hpp \
  moab/Types.hpp \
  moab/UnknownInterface.hpp \
  moab/Util.hpp \
  moab/WriteUtilIface.hpp \
  moab/WriterIface.hpp \
  moab/MeshTopoUtil.hpp \
  moab/RangeMap.hpp \
  MBEntityType.h \
  MBCN.h \
  MOAB_FCDefs.h \
  MBCN_protos.h 

MBCN_protos.h: MBCN.h $(top_srcdir)/itaps/mkprotos.sh
	$(top_srcdir)/itaps/mkprotos.sh MBCN MBCN $< $@ MOAB_FCDefs.h


MOAB_FCDefs.h: FCDefs.h
	cd .. && ./config.status src/MOAB_FCDefs.h

MOSTLYCLEANFILES = MOAB_FCDefs.h
