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

if ( MOAB_USE_MPI AND MPI_FOUND )
  add_executable ( mbparallelcomm_test mbparallelcomm_test.cpp )
  target_link_libraries( mbparallelcomm_test MOAB )
  set_source_files_properties( mbparallelcomm_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
  add_test( TestParallelComm-BcastDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} 0 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -1 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -2 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-Broadcast
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -3 ${MOAB_SOURCE_DIR}/parallel/ptest.cub )

  add_executable ( parallel_unit_tests parallel_unit_tests.cpp )
  target_link_libraries( parallel_unit_tests MOAB )
  add_test( TestParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/parallel_unit_tests ${MPIEXEC_POSTFLAGS} ${MOAB_SOURCE_DIR}/parallel/ptest.cub )
  set_source_files_properties( parallel_unit_tests.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )

  if ( MOAB_USE_HDF )
    add_executable( mhdf_parallel mhdf_parallel.c )
    target_link_libraries( mhdf_parallel MOAB MOABpar mhdf )
    add_test( TestMHDFParallel ${EXECUTABLE_OUTPUT_PATH}/mhdf_parallel )
    set_source_files_properties( mhdf_parallel.c
      COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES}" )
  endif ( MOAB_USE_HDF )
endif ( MOAB_USE_MPI AND MPI_FOUND )

