cmake_minimum_required( VERSION 2.8.12 )
if ( MPI_DIR )
  SET(CMAKE_USER_CXX_COMPILER_PATH ${MPI_DIR}/bin ${MPI_DIR})
  SET(CMAKE_GENERATOR_CC mpicc mpixlc_r mpixlc hcc mpxlc_r mpxlc sxmpicc mpifcc mpgcc mpcc cmpicc)
  SET(CMAKE_GENERATOR_CXX mpicxx mpic++ mpiCC sxmpic++ hcp mpxlC_r mpxlC mpixlcxx_r mpixlcxx mpg++ mpc++ mpCC cmpic++ mpiFCC)
  SET(CMAKE_GENERATOR_FC mpif95 mpxlf95_r mpxlf95 ftn mpif90 mpxlf90_r mpxlf90 mpf90 cmpif90c sxmpif90)
  SET(CMAKE_GENERATOR_F77 mpif77 hf77 mpxlf_r mpxlf mpifrt mpf77 cmpifc mpif95 mpxlf95_r mpxlf95 ftn mpif90 mpxlf90_r mpxlf90 mpf90 cmpif90c sxmpif90)
else ( MPI_DIR )
  SET(CMAKE_GENERATOR_CC gcc icc pgcc cc)
  SET(CMAKE_GENERATOR_CXX CCicpc pgCC pathCC sxc++ xlC_r xlC bgxlC_r bgxlC openCC sunCC crayCC icpc g++ c++ gpp aCC CC cxx cc++ cl.exe FCC KCC RCC)
  SET(CMAKE_GENERATOR_FC pathf95 ifort gfortran g95 f95 fort ifc efc openf95 sunf95 crayftn lf95 ftn xlf90 f90 pgf90 pghpf pathf90 epcf90 sxf90 openf90 sunf90)
  SET(CMAKE_GENERATOR_F77 ifort gfortran xlf f77 frt pgf77 pathf77 g77 cf77 fort77 fl32 af77)
endif ( MPI_DIR )
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
SET(CMAKE_CXX_FLAGS "-fPIC -DPIC")
SET(CMAKE_C_FLAGS "-fPIC -DPIC")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
SET(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

project( MOAB )

#Add our Cmake directory to the module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})

################################################################################
# Set up version info
################################################################################
include (config/GetAcInitVersion.cmake)
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

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0020 NEW)
  cmake_policy(SET CMP0042 OLD)
ENDIF ()

set ( abs_srcdir   ${CMAKE_SOURCE_DIR} )
set ( abs_builddir ${PROJECT_BINARY_DIR} )

include (config/CheckCompilerFlags.cmake)

################################################################################
# Install Related Settings
################################################################################

## Set the directory where the binaries will be stored
set( EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Directory where all executable will be stored"
)

## Set the directory where the libraries will be stored
set( LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Directory where all the libraries will be stored"
)
mark_as_advanced(
  EXECUTABLE_OUTPUT_PATH
  LIBRARY_OUTPUT_PATH)

include ( CheckIncludeFile )
include ( CheckFunctionExists )
include ( CheckTypeSize )

# Compiler defines... this should really be in a config file.
set( MOAB_DEFINES "" )
set( MOAB_LIBS )
set( MOAB_INSTALL_TARGETS )

################################################################################
# Options that the user controls
################################################################################
option ( BUILD_SHARED_LIBS   "Should shared or static libraries be created?"   ON  )
option ( MOAB_USE_SZIP       "Should build with szip support?"                 OFF )
option ( MOAB_USE_CGM        "Should build with CGM support?"                  OFF )
option ( MOAB_USE_CGNS       "Should build with CGNS support?"                 OFF )
if ( MPI_DIR )
  option ( MOAB_USE_MPI        "Should MOAB be compiled with MPI support?"       ON )
else ( MPI_DIR )
  option ( MOAB_USE_MPI        "Should MOAB be compiled with MPI support?"       OFF )
endif ( MPI_DIR )
option ( MOAB_USE_HDF        "Include HDF I/O in the build?"                   OFF )
option ( MOAB_USE_NETCDF     "Include NetCDF support (ExodusII) in the build?" OFF )
option ( MOAB_USE_PNETCDF    "Include parallel NetCDF support (ExodusII) in the build?" OFF )
option ( MOAB_USE_ZOLTAN     "Include Zoltan support for partitioning algorithms?" OFF )
option ( MOAB_ENABLE_TESTING "Enable Testing"                                  ON  )
option ( MOAB_FORCE_64_BIT_HANDLES "Force MBEntityHandle to be 64 bits (uint64_t)" OFF )
option ( MOAB_FORCE_32_BIT_HANDLES "Force MBEntityHandle to be 32 bits (uint32_t)" OFF )

option ( ENABLE_IMESH        "Should build IMESH?"       OFF )
option ( ENABLE_IGEOM        "Should build IGEOM?"       OFF )

mark_as_advanced(
    MOAB_FORCE_64_BIT_HANDLES
    MOAB_FORCE_32_BIT_HANDLES
  )

################################################################################
# Check for system include files
################################################################################
check_include_file( inttypes.h   MOAB_HAVE_INTTYPES_H )
check_include_file( stdint.h     MOAB_HAVE_STDINT_H )
check_include_file( stddef.h     MOAB_HAVE_STDDEF_H )
check_include_file( stdlib.h     MOAB_HAVE_STDLIB_H )
check_include_file( sys/types.h  MOAB_HAVE_SYS_TYPES_H )
check_type_size( size_t SIZE_T)
check_type_size( ptrdiff_t PTRDIFF_T )
set( MOAB_HAVE_SIZE_T ${HAVE_SIZE_T} )
set( MOAB_HAVE_PTRDIFF_T ${HAVE_PTRDIFF_T} )
set( HAVE_SYS_TYPES_H ${MOAB_HAVE_SYS_TYPES_H} )
set( HAVE_STDDEF_H    ${MOAB_HAVE_STDDEF_H} )
set( HAVE_STDINT_H    ${MOAB_HAVE_STDINT_H} )
set( HAVE_INTTYPES_H    ${MOAB_HAVE_INTTYPES_H} )
set( HAVE_STDLIB_H    ${MOAB_HAVE_STDLIB_H} )
check_include_file( memory.h     HAVE_MEMORY_H )

################################################################################
# Interger size Related Settings
################################################################################
if ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )
  message( FATAL_ERROR
      "You may not turn both MOAB_FORCE_64_BIT_HANDLES and MOAB_FORCE_32_BIT_HANDLES on. Turn one off to continue."
    )
endif ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )

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


################################################################################
# Find packages
################################################################################

if ( BUILD_SHARED_LIBS EQUAL OFF )
#Static start
set_target_properties(MOAB PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(MOAB PROPERTIES LINK_SEARCH_END_STATIC 1)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif ( BUILD_SHARED_LIBS EQUAL OFF )

find_package( verdict REQUIRED )

# check for MPI package
if ( MOAB_USE_MPI )
  find_package( MPI REQUIRED )
  # CMake FindMPI script is sorely lacking:
  if ( MPI_LIBRARY AND MPI_INCLUDE_PATH )
    set( MPI_FOUND 1 )
  endif ( MPI_LIBRARY AND MPI_INCLUDE_PATH )

  if ( MPI_FOUND )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    include_directories(MPI_INCLUDE_PATH)
  endif ( MPI_FOUND )
endif ( MOAB_USE_MPI )

if ( MOAB_USE_ZLIB )
  find_package( ZLIB REQUIRED )
endif ( )

if ( MOAB_USE_HDF )
  find_package( HDF5 REQUIRED )
  set( MOAB_DEFINES "-DHDF5_FILE ${MOAB_DEFINES}" )
  if ( HDF5_IS_PARALLEL AND MOAB_USE_MPI )
    set( MOAB_DEFINES "-DHDF5_PARALLEL ${MOAB_DEFINES}" )
    set (MOAB_HDF_HAVE_PARALLEL TRUE)
  endif ( HDF5_IS_PARALLEL AND MOAB_USE_MPI )
  set( MOAB_LIBS ${MOAB_LIBS} ${HDF5_LIBRARIES} )
  include_directories( ${HDF5_INCLUDE_DIR} src/io/mhdf/include )
endif ( MOAB_USE_HDF )

if ( MOAB_USE_NETCDF )
  find_package( NetCDF REQUIRED )
  set( MOAB_DEFINES "-DNETCDF_FILE ${MOAB_DEFINES}" )
  include_directories( ${NetCDF_INCLUDES} ${PNetCDF_INCLUDES} )
  set( MOAB_LIBS ${MOAB_LIBS} ${NetCDF_LIBRARIES} ${PNetCDF_LIBRARIES} )
endif ( MOAB_USE_NETCDF )

if ( MOAB_USE_ZOLTAN )
  find_package( Zoltan REQUIRED )
endif ( )

if ( MOAB_USE_CGM )
   find_package( CGM REQUIRED )
   set( MOAB_DEFINES "${CGM_DEFINES} -DCGM ${MOAB_DEFINES}" )
endif ()

if (MOAB_USE_CGNS)
  set( MOABIO_LIBS ${MOABIO_LIBS} ${CGNS_LIBRARIES} )
endif()

################################################################################
# Add Directories
################################################################################
add_subdirectory( src )
add_subdirectory( itaps )
add_subdirectory( tools )

################################################################################
# Testing Related Settings
################################################################################
#turn on ctest if we want testing
if ( MOAB_ENABLE_TESTING )
  enable_testing()
  add_subdirectory( test )
endif()

###############################################################################
#
###############################################################################
export(TARGETS ${MOAB_INSTALL_TARGETS}
              FILE "${PROJECT_BINARY_DIR}/MOABTargets.cmake")
install(EXPORT MOABTargets DESTINATION lib)

################################################################################
# Generate the MOABConfig.cmake file
################################################################################
set(CXX ${CMAKE_CXX_COMPILER})
set(CC ${CMAKE_C_COMPILER})
set(F77 "" )
set(FC "")
set(CXXFLAGS ${MOAB_DEFINES})
set(AM_CXXFLAGS)
set(CFLAGS ${MOAB_DEFINES})
set(AM_CFLAGS)
set(MOAB_FFLAGS)
set(AM_FFLAGS)

configure_file(MOABConfig.new.cmake.in "${PROJECT_BINARY_DIR}/MOABConfig.cmake" @ONLY)
configure_file(moab.make.in "${PROJECT_BINARY_DIR}/lib/moab.make" @ONLY)
if ( ENABLE_IGEOM )
  configure_file(itaps/igeom/iGeom-Defs.inc.in "${PROJECT_BINARY_DIR}/lib/iGeom-Defs.inc" @ONLY)
  install( FILES "${PROJECT_BINARY_DIR}/lib/iGeom-Defs.inc" DESTINATION lib )
endif ( ENABLE_IGEOM )

if ( ENABLE_IMESH )
  configure_file(itaps/imesh/iMesh-Defs.inc.in "${PROJECT_BINARY_DIR}/lib/iMesh-Defs.inc" @ONLY)
  install( FILES "${PROJECT_BINARY_DIR}/lib/iMesh-Defs.inc" DESTINATION lib )
endif ( ENABLE_IMESH )
install( FILES "${PROJECT_BINARY_DIR}/MOABConfig.cmake" "${PROJECT_BINARY_DIR}/lib/moab.make" DESTINATION lib )

# All done.
