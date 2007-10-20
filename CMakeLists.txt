project( MOAB )

  set ( MOAB_VERSION_MAJOR 3 )
  set ( MOAB_VERSION_MINOR 0 )
  # Define this for release versions. Undefine it for beta versions.
  #set ( MOAB_VERSION_PATCH 0 )

  set ( MOAB_VERSION_STRING "${MOAB_VERSION_MAJOR}" )
  if ( DEFINED MOAB_VERSION_MINOR )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING}.${MOAB_VERSION_MINOR}" )
    if ( DEFINED MOAB_VERSION_PATCH )
      set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING}.${MOAB_VERSION_PATCH}" )
    else ( DEFINED MOAB_VERSION_PATCH )
      set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (Beta)" )
    endif ( DEFINED MOAB_VERSION_PATCH )
  endif ( DEFINED MOAB_VERSION_MINOR )

  set ( CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${MOAB_SOURCE_DIR}/cmake" )
  set ( EXECUTABLE_OUTPUT_PATH ${MOAB_BINARY_DIR}/bin CACHE PATH "Path to executables" FORCE )
  set ( LIBRARY_OUTPUT_PATH    ${MOAB_BINARY_DIR}/bin CACHE PATH "Path to libraries"   FORCE )
  mark_as_advanced(
    EXECUTABLE_OUTPUT_PATH
    LIBRARY_OUTPUT_PATH
  )

  include ( CheckIncludeFile )
  include ( CheckFunctionExists )
  include ( CheckTypeSize )

  # Shared libraries
  option ( BUILD_SHARED_LIBS "Should shared or static libraries be created?" ON )

  # HANDLE SIZE
  option ( MOAB_FORCE_64_BIT_HANDLES "Force MBEntityHandle to be 64 bits (uint64_t)" OFF )
  option ( MOAB_FORCE_32_BIT_HANDLES "Force MBEntityHandle to be 32 bits (uint32_t)" OFF )
  mark_as_advanced(
    MOAB_FORCE_64_BIT_HANDLES
    MOAB_FORCE_32_BIT_HANDLES
  )

  if ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )
    message( FATAL_ERROR
      "You may not turn both MOAB_FORCE_64_BIT_HANDLES and MOAB_FORCE_32_BIT_HANDLES on. Turn one off to continue."
    )
  endif ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )

  check_include_file( inttypes.h   MOAB_HAVE_INTTYPES_H )
  check_include_file( stdint.h     MOAB_HAVE_STDINT_H )
  check_include_file( stddef.h     MOAB_HAVE_STDDEF_H )
  check_include_file( stdlib.h     MOAB_HAVE_STDLIB_H )
  check_include_file( sys/types.h  MOAB_HAVE_SYS_TYPES_H )
  set( HAVE_SYS_TYPES_H ${MOAB_HAVE_SYS_TYPES_H} )
  set( HAVE_STDDEF_H    ${MOAB_HAVE_STDDEF_H} )
  set( HAVE_STDINT_H    ${MOAB_HAVE_STDINT_H} )
  set( HAVE_INTTYPES_H    ${MOAB_HAVE_INTTYPES_H} )
  set( HAVE_STDLIB_H    ${MOAB_HAVE_STDLIB_H} )
  check_include_file( memory.h     HAVE_MEMORY_H )

  if ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )
    if ( MOAB_HAVE_INTTYPES_H )
      set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};inttypes.h" )
    endif ( MOAB_HAVE_INTTYPES_H )
    if ( MOAB_HAVE_STDLIB_H )
      set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};stdlib.h" )
      #set ( CHECK_TYPE_SIZE_PREMAIN "${CHECK_TYPE_SIZE_PREMAIN}\n#include <stdlib.h>\n" )
    endif ( MOAB_HAVE_STDLIB_H )
    check_type_size(  size_t       HAVE_SIZE_T )
    check_type_size(  ptrdiff_t    HAVE_PTRDIFF_T )
    set ( MOAB_HAVE_SIZE_T ${HAVE_SIZE_T} )
    set ( MOAB_HAVE_PTRDIFF_T ${HAVE_PTRDIFF_T} )
  endif ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )

  # MPI
  option ( MOAB_USE_MPI "Should MOAB be compiled with MPI support?" OFF )
  if ( MOAB_USE_MPI )
    find_package( MPI )
  endif ( MOAB_USE_MPI )

  # HDF5
  find_package( HDF5 )
  if ( HDF5_FOUND )
    check_function_exists( H5Pset_fapl_mpio MOAB_HDF_HAVE_PARALLEL )
  endif ( HDF5_FOUND )
  # NetCDF
  # Documentation
  # Optional tools
  option ( MOAB_BUILD_MBCONVERT    "Build the MOAB mesh converter tool?"  ON )
  option ( MOAB_BUILD_HEXMODOPS    "Build the hexahedral modops tool?"    ON )
  option ( MOAB_BUILD_MBCHACO      "Build MOAB Chaco interface?"          OFF )
  option ( MOAB_BUILD_MBPERF       "Build the MOAB performance tool?"     ON )
  option ( MOAB_BUILD_QVDUAL       "Build the qvdual application?"        OFF )
  option ( MOAB_BUILD_MBSIZE       "Build the mbsize tool?"               ON )
  option ( MOAB_BUILD_MBSKIN       "Build the mbskin tool?"               ON )
  option ( MOAB_BUILD_MBTAGPROP    "Build the mbtagprop tool?"            ON )
  option ( MOAB_BUILD_SPHEREDECOMP "Build the sphere decomposition tool?" ON )
  option ( MOAB_BUILD_MBSURFPLOT   "Build the mbsurfplot application?"    ON )
  option ( MOAB_BUILD_MBZOLTAN     "Build the mbzoltan partitioner?"      OFF )
  option ( MOAB_BUILD_DAGMC        "Build the dagmc tool?"                ON )
  # iMesh
  option ( MOAB_BUILD_IMESH        "Build the iMesh interface?"           ON )

  include( ${MOAB_SOURCE_DIR}/cmake/AutoconfHeader.cmake )
  # Define some more variables so they will be substituted properly in the autoconf files.
  set( MB_VERSION_MAJOR  ${MOAB_VERSION_MAJOR} )
  set( MB_VERSION_MINOR  ${MOAB_VERSION_MINOR} )
  set( MB_VERSION_PATCH  ${MOAB_VERSION_PATCH} )
  set( MB_VERSION_STRING "\"${MOAB_VERSION_STRING}\"" )

  # Best we can do for now since there's no config.h.in in the repository...
  configure_file( ${MOAB_SOURCE_DIR}/config.h.cmake      ${MOAB_BINARY_DIR}/config.h )
  # ... but the line below is what should *really* happen
  #autoconf_header( ${MOAB_SOURCE_DIR}/config.h.in ${MOAB_BINARY_DIR}/config.h )

  autoconf_header( ${MOAB_SOURCE_DIR}/MBVersion.h.in ${MOAB_BINARY_DIR}/MBVersion.h )
  autoconf_header( ${MOAB_SOURCE_DIR}/MBEntityHandle.h.in ${MOAB_BINARY_DIR}/MBEntityHandle.h )

  # MOAB Library
  set ( MOAB_LIB_SRCS
    AEntityFactory.cpp
    DenseTagCollections.cpp
    DualTool.cpp
    EntitySequence.cpp
    EntitySequenceManager.cpp
    ExoIIUtil.cpp
    FileOptions.cpp
    FileTokenizer.cpp
    GeomTopoTool.cpp
    HigherOrderFactory.cpp
    HomXform.cpp
    MBAdaptiveKDTree.cpp
    MBAffineXform.cpp
    MBAlloc.cpp
    MBAxisBox.cpp
    MBBits.cpp
    MBCN.cpp
    MBCartVect.cpp
    MBHandleUtils.cpp
    MBMatrix3.cpp
    MBCore.cpp
    MBFactory.cpp
    MBGeomUtil.cpp
    MBMeshSet.cpp
    MBOrientedBox.cpp
    MBOrientedBoxTreeTool.cpp
    MBRange.cpp
    MBRangeSeqIntersectIter.cpp
    MBReadUtil.cpp
    MBReaderWriterSet.cpp
    MBSkinner.cpp
    MBUtil.cpp
    MBWriteUtil.cpp
    MeshSetSequence.cpp
    MeshTopoUtil.cpp
    PolyEntitySequence.cpp
    ReadGmsh.cpp
    ReadSTL.cpp
    ReadVtk.cpp
    ScdElementSeq.cpp
    ScdVertexSeq.cpp
    SparseTagCollections.cpp
    TagServer.cpp
    Tqdcfr.cpp
    VtkUtil.cpp
    WriteAns.cpp
    WriteGMV.cpp
    WriteGmsh.cpp
    WriteSTL.cpp
    WriteTemplate.cpp
    WriteVtk.cpp
  )
  set_source_files_properties( ${MOAB_LIB_SRCS}
    COMPILE_FLAGS "-DIS_BUILDING_MB"
  )
  include_directories(
    ${MOAB_SOURCE_DIR}
    ${MOAB_BINARY_DIR}
  )
  add_library( moab
    ${MOAB_LIB_SRCS}
  )

  if ( MOAB_USE_HDF )
    set ( MOAB_MHDF_LIB_SRCS
      ReadHDF5.cpp
      WriteHDF5.cpp
    )
    include_directories(
      ${HDF5_INCLUDE_DIRECTORIES}
      ${MOAB_SOURCE_DIR}/mhdf/include
    )
    add_library( mhdf
      ${MOAB_MHDF_LIB_SRCS}
    )
    target_link_libraries( mdf
      ${HDF5_LIBRARIES}
    )
  endif ( MOAB_USE_HDF )

  ## Tools/libraries dependent on MOAB:

  # Zoltan
  # Chaco

  # QVDual
  if ( MOAB_BUILD_QVDUAL )
    add_subdirectory( tools/qvdual )
  endif ( MOAB_BUILD_QVDUAL )

  # CGM
  if ( MOAB_BUILD_DAGMC )
    add_subdirectory( tools/dagmc )
  endif ( MOAB_BUILD_DAGMC )

  # iMesh Babel


