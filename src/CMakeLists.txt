
  # MOAB Library
  set ( MOAB_LIB_SRCS
    AEntityFactory.cpp
    BitTagServer.cpp
    BSPTreePoly.cpp
    DualTool.cpp
    EntitySequence.cpp
    FileOptions.cpp
    GeomTopoTool.cpp
    HigherOrderFactory.cpp
    HomXform.cpp
    AdaptiveKDTree.cpp
    AffineXform.cpp
    AxisBox.cpp
    BSPTree.cpp
    CN.cpp
    CartVect.cpp
    Matrix3.cpp
    Core.cpp
    Factory.cpp
    GeomUtil.cpp
    MeshSet.cpp
    OrientedBox.cpp
    OrientedBoxTreeTool.cpp
    Range.cpp
    RangeSeqIntersectIter.cpp
    ReadUtil.cpp
    ReaderWriterSet.cpp
    Skinner.cpp
    SysUtil.cpp
    Util.cpp
    WriteUtil.cpp
    MeshSetSequence.cpp
    MeshTopoUtil.cpp
    PolyElementSeq.cpp
    ScdElementData.cpp
    ScdVertexData.cpp
    SequenceData.cpp
    SequenceManager.cpp
    SparseTagSuperCollection.cpp
    SparseTagCollection.cpp
    StructuredElementSeq.cpp
    SweptElementData.cpp
    SweptElementSeq.cpp
    SweptVertexData.cpp
    TagInfo.cpp
    TagServer.cpp
    TypeSequenceManager.cpp
    UnstructuredElemSeq.cpp
    VertexSequence.cpp
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

  autoconf_header( ${MOAB_SOURCE_DIR}/src/moab/MBVersion.h.in ${MOAB_BINARY_DIR}/src/moab/MBVersion.h )
  autoconf_header( ${MOAB_SOURCE_DIR}/src/moab/MBEntityHandle.h.in ${MOAB_BINARY_DIR}/src/moab/MBEntityHandle.h )
  file( WRITE ${MOAB_BINARY_DIR}/src/MBCN_protos.h "" )
