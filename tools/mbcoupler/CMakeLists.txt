include_directories(
  ${CMAKE_CURRENT_SOURCE_DIRECTORY}
)

set( MBCOUPLER_SRCS
  MBElemUtil.cpp
  tensor.c
  poly.c
  findpt.c
  errmem.c
)

add_library( mbcoupler
  ${MBCOUPLER_SRCS}
)

enable_testing()

# A test for the coupler library
add_executable( findpt_test findpt_test.cpp )
target_link_libraries( findpt_test mbcoupler MOAB )
add_test( TestMBCouplerFindPt ${EXECUTABLE_OUTPUT_PATH}/findpt_test )
