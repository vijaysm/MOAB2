include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
    ${MOAB_SOURCE_DIR}/src/moab/point_locater/lotte
    ${MOAB_SOURCE_DIR}/itaps
    ${MOAB_SOURCE_DIR}/itaps/imesh
    ${MOAB_BINARY_DIR}
    ${MOAB_BINARY_DIR}/src
    ${MOAB_BINARY_DIR}/src/parallel
)

set( MBCOUPLER_SRCS
  Coupler.cpp
  ElemUtil.cpp
)

set_source_files_properties( ${MBCOUPLER_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES}"
)

add_library( mbcoupler
  ${MBCOUPLER_SRCS}
)

target_link_libraries( mbcoupler MOAB iMesh )
if ( MOAB_USE_MPI )
  target_link_libraries( mbcoupler MOABpar )
endif ( MOAB_USE_MPI )

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

  add_executable( ssn_test ssn_test.cpp )
  set_source_files_properties( ssn_test.cpp
    COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( ssn_test mbcoupler MOAB MOABpar )
  add_test( TestMBCoupler ${EXECUTABLE_OUTPUT_PATH}/ssn_test )
endif ( MOAB_USE_MPI )
