include_directories(
  ${CMAKE_CURRENT_SOURCE_DIRECTORY}
)

set( MBCOUPLER_SRCS
  MBCoupler.cpp
  MBElemUtil.cpp
  errmem.c
  findpt.c
  poly.c
  tensor.c
)

add_library( mbcoupler
  ${MBCOUPLER_SRCS}
)
target_link_libraries( mbcoupler MOAB )
if ( MOAB_USE_MPI )
  target_link_libraries( mbcoupler MOABpar )
endif ( MOAB_USE_MPI )
set_source_files_properties( ${MBCOUPLER_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES}"
)

enable_testing()

# A test for the coupler library
add_executable( findpt_test findpt_test.cpp )
target_link_libraries( findpt_test mbcoupler MOAB )
if ( MOAB_USE_MPI )
  target_link_libraries( findpt_test MOABpar )
endif ( MOAB_USE_MPI )
add_test( TestMBCouplerFindPt ${EXECUTABLE_OUTPUT_PATH}/findpt_test )

if ( MOAB_USE_MPI )
  # Exercise the coupler and not just findpt:
  add_executable( mbcoupler_test mbcoupler_test.cpp )
  set_source_files_properties( mbcoupler_test.cpp
    COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbcoupler_test mbcoupler MOAB MOABpar )
  add_test( TestMBCoupler ${EXECUTABLE_OUTPUT_PATH}/mbcoupler_test )
endif ( MOAB_USE_MPI )
