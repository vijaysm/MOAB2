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

set_source_files_properties( ${moab_PARALLEL_SRCS}
  COMPILE_FLAGS "-DIS_BUILDING_MB"
)
add_library( MOABpar
  ${moab_PARALLEL_SRCS}
)

target_link_libraries( MOABpar
  MOAB
  ${MPI_LIBRARY}
  ${MPI_EXTRA_LIBRARY}
)
