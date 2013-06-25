
  # MOAB Library
  set ( MOAB_LIB_SRCS
    AdaptiveKDTree.cpp
    AEntityFactory.cpp
    AffineXform.cpp
    AxisBox.cpp
    BitPage.cpp
    BitTag.cpp
    BSPTree.cpp
    BSPTreePoly.cpp
    CN.cpp
    CartVect.cpp
    Core.cpp
    DebugOutput.cpp
    DenseTag.cpp
    DualTool.cpp
    EntitySequence.cpp
    Factory.cpp
    FileOptions.cpp
    GeomUtil.cpp
    GeomTopoTool.cpp
    HigherOrderFactory.cpp
    HomXform.cpp
    MeshSet.cpp
    MeshSetSequence.cpp
    MeshTag.cpp
    MeshTopoUtil.cpp
    OrientedBox.cpp
    OrientedBoxTreeTool.cpp
    PolyElementSeq.cpp
    Range.cpp
    RangeSeqIntersectIter.cpp
    ReadUtil.cpp
    ReaderWriterSet.cpp
    ScdElementData.cpp
    ScdInterface.cpp
    ScdVertexData.cpp
    SequenceData.cpp
    SequenceManager.cpp
    SetIterator.cpp
    Skinner.cpp
    SparseTag.cpp
    StructuredElementSeq.cpp
    SweptElementData.cpp
    SweptElementSeq.cpp
    SweptVertexData.cpp
    SysUtil.cpp
    TagInfo.cpp
    Types.cpp
    TypeSequenceManager.cpp
    UnstructuredElemSeq.cpp
    Util.cpp
    VarLenDenseTag.cpp
    VarLenSparseTag.cpp
    VertexSequence.cpp
    WriteUtil.cpp
  )
  include_directories(
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_BINARY_DIR}/src
    ${MOAB_SOURCE_DIR}/src/io
    ${MOAB_BINARY_DIR}/src/io
  )

  if ( MOAB_USE_MPI )
    find_package( MPI )
    # CMake FindMPI script is sorely lacking:
    if ( MPI_LIBRARY AND MPI_INCLUDE_PATH )
      set( MPI_FOUND 1 )
    endif ( MPI_LIBRARY AND MPI_INCLUDE_PATH )

    if ( MPI_FOUND )
      set ( MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
      include_directories(
        ${MPI_INCLUDE_PATH}
        ${MOAB_SOURCE_DIR}/parallel
      )
    endif ( MPI_FOUND )
  endif ( MOAB_USE_MPI )

  set_source_files_properties( ${MOAB_LIB_SRCS}
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}"
  )
  add_library( MOAB
    ${MOAB_LIB_SRCS}
  )
  
  add_subdirectory( io )
  target_link_libraries( MOAB moabio )

  if ( MOAB_USE_MPI AND MPI_FOUND )
    add_subdirectory( parallel )
    target_link_libraries( MOAB MOABpar )
  endif ( MOAB_USE_MPI AND MPI_FOUND )

  include( ${MOAB_SOURCE_DIR}/config/AutoconfHeader.cmake )
  # Define some more variables so they will be substituted properly in the autoconf files.
  set( MB_VERSION        "\"${MOAB_VERSION}\"" )
  set( MB_VERSION_MAJOR  ${MOAB_VERSION_MAJOR} )
  set( MB_VERSION_MINOR  ${MOAB_VERSION_MINOR} )
  if ( DEFINED MOAB_VERSION_PATCH )
    set( MB_VERSION_PATCH  ${MOAB_VERSION_PATCH} )
  endif ( DEFINED MOAB_VERSION_PATCH )
  set( MB_VERSION_STRING "\"${MOAB_VERSION_STRING}\"" )

  autoconf_header( ${MOAB_SOURCE_DIR}/src/moab/Version.h.in ${MOAB_BINARY_DIR}/src/moab/Version.h )
  autoconf_header( ${MOAB_SOURCE_DIR}/src/moab/EntityHandle.hpp.in ${MOAB_BINARY_DIR}/src/moab/EntityHandle.hpp )
  file( WRITE ${MOAB_BINARY_DIR}/src/MBCN_protos.h "" )
