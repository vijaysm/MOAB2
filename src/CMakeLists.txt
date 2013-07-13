
  # MOAB Library
  set ( MOABCORE_LIB_SRCS
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

  set ( MOABPTLOC_LIB_SRCS
    lotte/findpt.c
    lotte/errmem.c
    lotte/poly.c
    lotte/tensor.c
  )

  set ( MOABIO_LIB_SRCS
    io/IODebugTrack.cpp
    io/ExoIIUtil.cpp
    io/FileTokenizer.cpp
    io/GmshUtil.cpp
    io/ReadABAQUS.cpp
    io/ReadGmsh.cpp
    io/ReadGCRM.cpp
    io/ReadIDEAS.cpp
    io/ReadMCNP5.cpp
    io/ReadNASTRAN.cpp
    io/ReadSmf.cpp
    io/ReadSms.cpp
    io/ReadSTL.cpp
    io/ReadTemplate.cpp
    io/ReadTetGen.cpp
    io/ReadTxt.cpp
    io/ReadVtk.cpp
    io/SMF_State.cpp
    io/Tqdcfr.cpp
    io/VtkUtil.cpp
    io/WriteAns.cpp
    io/WriteGMV.cpp
    io/WriteGmsh.cpp
    io/WriteSTL.cpp
    io/WriteSmf.cpp
    io/WriteTemplate.cpp
    io/WriteVtk.cpp
    ReaderWriterSet.cpp
  )
  include_directories(
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/src/io
    ${MOAB_SOURCE_DIR}/src/parallel
    ${MOAB_SOURCE_DIR}/src/moab/point_locater/lotte
  )

  if ( NetCDF_FOUND )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DNETCDF_FILE" )
    set ( MOABIO_LIB_SRCS
      ${MOABIO_LIB_SRCS}
      io/ReadNC.cpp
      io/ReadNCDF.cpp
      io/WriteNCDF.cpp
      io/WriteSLAC.cpp
      io/NCHelper.cpp
      io/NCHelperEuler.cpp
      io/NCHelperFV.cpp
      io/NCHelperHOMME.cpp
      SpectralMeshTool.cpp
    )
    include_directories(
      ${NetCDF_INCLUDE_DIRECTORIES}
    )
  endif ( NetCDF_FOUND )

  if ( HDF5_FOUND )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DHDF5_FILE" )
    check_function_exists( H5Pset_fapl_mpio MOAB_HDF_HAVE_PARALLEL )
    set ( MOABIO_LIB_SRCS
      ${MOABIO_LIB_SRCS}
      io/HDF5Common.cpp
      io/ReadHDF5.cpp
      io/ReadHDF5Dataset.cpp
      io/ReadHDF5VarLen.cpp
      io/WriteHDF5.cpp
    )

    include_directories(
      ${HDF5_INCLUDE_DIRECTORIES}
      io/mhdf/include
    )
    add_subdirectory( io/mhdf )
  endif ( HDF5_FOUND )


  SET(MOAB_LIB_SRCS ${MOABCORE_LIB_SRCS} ${MOABPTLOC_LIB_SRCS} ${MOABIO_LIB_SRCS})

  set_source_files_properties( ${MOAB_LIB_SRCS}
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}"
  )
  add_library( MOAB
    ${MOAB_LIB_SRCS}
  )

  if ( MOAB_USE_NETCDF AND NetCDF_FOUND )
    target_link_libraries( MOAB ${NetCDF_LIBRARIES} )
  endif ( MOAB_USE_NETCDF AND NetCDF_FOUND )

  if ( MOAB_USE_HDF AND HDF5_FOUND )
    target_link_libraries( MOAB mhdf ${HDF5_LIBRARIES} )
  endif ( MOAB_USE_HDF AND HDF5_FOUND )

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
