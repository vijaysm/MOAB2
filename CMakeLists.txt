project( MOAB )

  cmake_minimum_required( VERSION 2.4 )
  if ( COMMAND cmake_policy )
    cmake_policy( SET CMP0003 NEW )
  endif ( COMMAND cmake_policy )
  # Always build tests
  enable_testing()

  list (APPEND CMAKE_MODULE_PATH "${MOAB_SOURCE_DIR}/config")
  include (config/GetAcInitVersion.cmake)

  # enable_language( Fortran )

  get_ac_init_version()
  set ( MOAB_VERSION_MAJOR  "${MAJOR_VERSION}"  )
  set ( MOAB_VERSION_MINOR  "${MINOR_VERSION}"  )
  set ( MOAB_VERSION        "${VERSION_STRING}" )
  set ( MOAB_VERSION_STRING "${VERSION_STRING}" )
  if ( DEFINED PATCH_VERSION )
    set ( MOAB_VERSION_PATCH "${PATCH_VERSION}" )
  else ( DEFINED PATCH_VERSION )
    if ( MOAB_VERSION_MINOR EQUAL 99 )
      set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (alpha)" )
    else ( MOAB_VERSION_MINOR EQUAL 99 )
      set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (beta)" )
    endif ( MOAB_VERSION_MINOR EQUAL 99 )
  endif ( DEFINED PATCH_VERSION )

  set ( EXECUTABLE_OUTPUT_PATH ${MOAB_BINARY_DIR}/bin CACHE PATH "Path to executables" FORCE )
  set ( LIBRARY_OUTPUT_PATH    ${MOAB_BINARY_DIR}/bin CACHE PATH "Path to libraries"   FORCE )
  mark_as_advanced(
    EXECUTABLE_OUTPUT_PATH
    LIBRARY_OUTPUT_PATH
  )

  # Compiler defines... this should really be in a config file.
  set( MOAB_DEFINES "" )

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

  try_compile( TEMPLATE_DEFS_INCLUDED ${MOAB_BINARY_DIR} ${MOAB_SOURCE_DIR}/config/TemplateSpecialization.cxx OUTPUT_VARIABLE BLORT )
  if ( TEMPLATE_DEFS_INCLUDED )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DTEMPLATE_DEFS_INCLUDED" )
  endif ( TEMPLATE_DEFS_INCLUDED )

  # find Verdict
  find_package( verdict REQUIRED )

  # Build options
  option ( MOAB_USE_MPI    "Should MOAB be compiled with MPI support?"       OFF )
  option ( MOAB_USE_HDF    "Include HDF I/O in the build?"                   OFF )
  option ( MOAB_USE_NETCDF "Include NetCDF support (ExodusII) in the build?" OFF )

  # iMesh
  option ( MOAB_BUILD_IMESH        "Build the iMesh interface?"           ON )

  # check for MPI package
  if ( MOAB_USE_MPI )
    find_package( MPI )
    # CMake FindMPI script is sorely lacking:
    if ( MPI_LIBRARY AND MPI_INCLUDE_PATH )
      set( MPI_FOUND 1 )
    endif ( MPI_LIBRARY AND MPI_INCLUDE_PATH )

    if ( MPI_FOUND )
      set ( MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
    endif ( MPI_FOUND )
  endif ( MOAB_USE_MPI )

  if ( MOAB_USE_NETCDF )
    find_package( NetCDF )
  endif ( MOAB_USE_NETCDF )

  if ( MOAB_USE_HDF )
    # HDF5
    find_package( HDF5 )
  endif ( MOAB_USE_HDF )

  add_subdirectory( src )
  add_subdirectory( itaps/imesh )
  add_subdirectory( tools )
  add_subdirectory( test )

