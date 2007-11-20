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

  # find Verdict
  find_package( verdict REQUIRED )

  # Build options
  option ( MOAB_USE_MPI    "Should MOAB be compiled with MPI support?"       OFF )
  option ( MOAB_USE_HDF    "Include HDF I/O in the build?"                   OFF )
  option ( MOAB_USE_NETCDF "Include NetCDF support (ExodusII) in the build?" OFF )
  
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

  # Compiler defines... this should really be in a config file.
  set( MOAB_DEFINES "" )

  # MOAB Library
  set ( MOAB_LIB_SRCS
    AEntityFactory.cpp
    DualTool.cpp
    EntitySequence.cpp
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
    MBEdgeSizeEvaluator.cpp
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
    PolyElementSeq.cpp
    ReadGmsh.cpp
    ReadSTL.cpp
    ReadVtk.cpp
    ScdElementData.cpp
    ScdVertexData.cpp
    SequenceData.cpp
    SequenceManager.cpp
    SparseTagCollections.cpp
    StructuredElementSeq.cpp
    TagServer.cpp
    Tqdcfr.cpp
    TypeSequenceManager.cpp
    UnstructuredElemSeq.cpp
    VertexSequence.cpp
    VtkUtil.cpp
    WriteAns.cpp
    WriteGMV.cpp
    WriteGmsh.cpp
    WriteSTL.cpp
    WriteTemplate.cpp
    WriteVtk.cpp
  )
  include_directories(
    ${MOAB_SOURCE_DIR}
    ${MOAB_BINARY_DIR}
  )

  if ( MOAB_USE_NETCDF )
    find_package( NetCDF )
    if ( NetCDF_FOUND )
      set ( MOAB_DEFINES "${MOAB_DEFINES} -DNETCDF_FILE" )
      set ( MOAB_LIB_SRCS
        ${MOAB_LIB_SRCS}
        ReadNCDF.cpp
        WriteNCDF.cpp
        WriteSLAC.cpp
      )
      include_directories(
        ${NetCDF_INCLUDE_DIRECTORIES}
      )
    endif ( NetCDF_FOUND )
  endif ( MOAB_USE_NETCDF )

  if ( MOAB_USE_HDF )
    # HDF5
    find_package( HDF5 )
    if ( HDF5_FOUND )
      set ( MOAB_DEFINES "${MOAB_DEFINES} -DHDF5_FILE" )
      check_function_exists( H5Pset_fapl_mpio MOAB_HDF_HAVE_PARALLEL )
      set ( MOAB_LIB_SRCS
        ${MOAB_LIB_SRCS}
        ReadHDF5.cpp
        WriteHDF5.cpp
      )
      include_directories(
        ${HDF5_INCLUDE_DIRECTORIES}
        ${MOAB_SOURCE_DIR}/mhdf/include
      )
      add_subdirectory( mhdf )
    endif ( HDF5_FOUND )
  endif ( MOAB_USE_HDF )

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
  if ( MOAB_USE_HDF AND HDF5_FOUND )
    target_link_libraries( MOAB
      mhdf
      ${HDF5_LIBRARIES}
    )
  endif ( MOAB_USE_HDF AND HDF5_FOUND )

  if ( MOAB_USE_MPI AND MPI_FOUND )
    add_subdirectory( parallel )
    target_link_libraries( MOAB
      MOABpar
    )
  endif ( MOAB_USE_MPI AND MPI_FOUND )

  if ( MOAB_USE_NETCDF AND NetCDF_FOUND )
    target_link_libraries( MOAB
      ${NetCDF_LIBRARIES}
    )
  endif ( MOAB_USE_NETCDF AND NetCDF_FOUND )
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

  include( ${MOAB_SOURCE_DIR}/cmake/AutoconfHeader.cmake )
  # Define some more variables so they will be substituted properly in the autoconf files.
  set( MB_VERSION_MAJOR  ${MOAB_VERSION_MAJOR} )
  set( MB_VERSION_MINOR  ${MOAB_VERSION_MINOR} )
  set( MB_VERSION_PATCH  ${MOAB_VERSION_PATCH} )
  set( MB_VERSION_STRING "\"${MOAB_VERSION_STRING}\"" )

  autoconf_header( ${MOAB_SOURCE_DIR}/MBVersion.h.in ${MOAB_BINARY_DIR}/MBVersion.h )
  autoconf_header( ${MOAB_SOURCE_DIR}/MBEntityHandle.h.in ${MOAB_BINARY_DIR}/MBEntityHandle.h )


# Testing
  # Always build tests
  enable_testing()
  # Point to the input decks for the tests
  configure_file( ${MOAB_SOURCE_DIR}/cmake/TestingConfig.h.in ${MOAB_BINARY_DIR}/testdir.h )

  if ( MOAB_USE_MPI AND MPI_FOUND )
     add_executable ( mbparallelcomm_test mbparallelcomm_test.cpp )
     target_link_libraries( mbparallelcomm_test MOAB )
     add_test( TestParallelCommunications ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test )
     set_source_files_properties( mbparallelcomm_test.cpp
       COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
  endif ( MOAB_USE_MPI AND MPI_FOUND )

  add_executable( moab_test MBTest.cpp )
  target_link_libraries( moab_test MOAB )
  add_test( TestMOAB ${EXECUTABLE_OUTPUT_PATH}/moab_test -nostress )
  # FIXME: moab_test depends on: test/mb_big_test.g test/cell1.gen test/cell2.gen

  add_executable( seq_man_test TestTypeSequenceManager.cpp )
  target_link_libraries( seq_man_test MOAB )
  set_source_files_properties( TestTypeSequenceManager.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
   add_test( TestTypeSequenceManager ${EXECUTABLE_OUTPUT_PATH}/seq_man_test )

  add_executable( internals_test internals_test.cpp )
  set_source_files_properties( internals_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
  target_link_libraries( internals_test MOAB )
  add_test( TestInternals ${EXECUTABLE_OUTPUT_PATH}/internals_test )

  add_executable( homxform_test HomXform.cpp )
  set_source_files_properties( HomXform.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES}" )
  target_link_libraries( homxform_test MOAB )
  add_test( TestHomXform ${EXECUTABLE_OUTPUT_PATH}/homxform_test )

  add_executable( scdseq_test scdseq_test.cpp )
  set_source_files_properties( scdseq_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
  target_link_libraries( scdseq_test MOAB )
  add_test( TestSCDSeq ${EXECUTABLE_OUTPUT_PATH}/scdseq_test )

  add_executable( tqdcfr Tqdcfr.cpp )
  set_source_files_properties( Tqdcfr.cpp
    COMPILE_FLAGS "-DTEST -DIS_BUILDING_MB ${MOAB_DEFINES}" )
  target_link_libraries( tqdcfr MOAB )
  add_test( TestTQDCFR ${EXECUTABLE_OUTPUT_PATH}/tqdcfr )


