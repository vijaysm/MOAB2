set( moab_PARALLEL_SRCS
  ParallelComm.cpp
  ProcConfig.cpp
  ParallelData.cpp
  ParallelMergeMesh.cpp
  ReadParallel.hpp      ReadParallel.cpp
  SharedSetData.hpp     SharedSetData.cpp
  gs.cpp
)

#generate the moab_mpi_config.h file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/moab_mpi_config.h.in "${CMAKE_CURRENT_BINARY_DIR}/moab_mpi_config.h" @ONLY)

set(moab_PARALLEL_INSTALL_HDRS moab/gs.hpp
                               moab_mpi.h
                               moab_mpi_config.h
                               moab/ParallelComm.hpp
                               moab/ParallelMergeMesh.hpp
                               moab/ProcConfig.hpp
                               moab/ParallelData.hpp
                               MBParallelConventions.h )

include_directories(
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
    ${CMAKE_CURRENT_BINARY_DIR}
  )

if ( MOAB_USE_HDF AND MOAB_HDF_HAVE_PARALLEL )
  set( moab_PARALLEL_SRCS ${moab_PARALLEL_SRCS}
                          WriteHDF5Parallel.cpp )
  
  set( moab_PARALLEL_INSTALL_HDRS ${moab_PARALLEL_INSTALL_HDRS}
                                  WriteHDF5Parallel.hpp )
endif()

if ( MOAB_USE_HDF )
  include_directories(
    ${MOAB_SOURCE_DIR}/mhdf/include
  )
endif ( MOAB_USE_HDF )

set_source_files_properties( ${moab_PARALLEL_SRCS}
  COMPILE_FLAGS "-DIS_BUILDING_MB -DSRCDIR=\"${MOAB_SOURCE_DIR}\" ${MOAB_DEFINES}"
)
add_library( MOABpar STATIC ${moab_PARALLEL_INSTALL_HDRS} ${moab_PARALLEL_SRCS} )

target_link_libraries( MOABpar
  ${MPI_LIBRARY}
  ${MPI_EXTRA_LIBRARY}
)

include(AutoconfHeader)
moab_install_headers(${moab_PARALLEL_INSTALL_HDRS})

