set( moab_PARALLEL_SRCS
  ParallelComm.cpp
  ProcConfig.cpp
  ParallelData.cpp
  ReadParallel.cpp
  gs.cpp
  TupleList.cpp
  SharedSetData.cpp
)

include_directories(
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
  )

if ( MOAB_USE_HDF AND MOAB_HDF_HAVE_PARALLEL )
  set( moab_PARALLEL_SRCS
    ${moab_PARALLEL_SRCS}
    WriteHDF5Parallel.cpp
  )
endif ( MOAB_USE_HDF AND MOAB_HDF_HAVE_PARALLEL )

if ( MOAB_USE_HDF )
  include_directories(
    ${MOAB_SOURCE_DIR}/mhdf/include
  )
endif ( MOAB_USE_HDF )

# On Mac OS X, shared libraries may not have circular
# dependencies (e.g., MOAB depends on MOABpar depends on MOAB).
# Force a static MOABpar library in this case.
if ( APPLE AND BUILD_SHARED_LIBS )
  set ( MOABpar_LIB_TYPE "STATIC" )
else ( APPLE AND BUILD_SHARED_LIBS )
  if ( BUILD_SHARED_LIBS )
    set ( MOABpar_LIB_TYPE "SHARED" )
  else ( BUILD_SHARED_LIBS )
    set ( MOABpar_LIB_TYPE "STATIC" )
  endif ( BUILD_SHARED_LIBS )
endif ( APPLE AND BUILD_SHARED_LIBS )

set_source_files_properties( ${moab_PARALLEL_SRCS}
  COMPILE_FLAGS "-DIS_BUILDING_MB -DSRCDIR=\"${MOAB_SOURCE_DIR}\" ${MOAB_DEFINES}"
)
add_library( MOABpar ${MOABpar_LIB_TYPE}
  ${moab_PARALLEL_SRCS}
)

target_link_libraries( MOABpar
  ${MPI_LIBRARY}
  ${MPI_EXTRA_LIBRARY}
)

