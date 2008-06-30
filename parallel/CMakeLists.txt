set( moab_PARALLEL_SRCS
  MBParallelComm.cpp
  MBProcConfig.cpp
  MBParallelData.cpp
  ReadParallel.cpp
  crystal.c
  errmem.c
  transfer.c
  gs.c
  tuple_list.c
  sort.c
)

if ( MOAB_USE_HDF AND MOAB_HDF_HAVE_PARALLEL )
  set( moab_PARALLEL_SRCS
    ${moab_PARALLEL_SRCS}
    WriteHDF5Parallel.cpp
  )
endif ( MOAB_USE_HDF AND MOAB_HDF_HAVE_PARALLEL )

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
  COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}"
)
add_library( MOABpar ${MOABpar_LIB_TYPE}
  ${moab_PARALLEL_SRCS}
)

target_link_libraries( MOABpar
  ${MPI_LIBRARY}
  ${MPI_EXTRA_LIBRARY}
)

if ( MOAB_USE_MPI AND MPI_FOUND )
   add_executable ( mbparallelcomm_test mbparallelcomm_test.cpp )
   target_link_libraries( mbparallelcomm_test MOAB )
   add_test( TestParallelCommunications ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test )
   set_source_files_properties( mbparallelcomm_test.cpp
     COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
endif ( MOAB_USE_MPI AND MPI_FOUND )

