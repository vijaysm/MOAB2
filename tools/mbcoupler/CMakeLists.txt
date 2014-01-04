include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_SOURCE_DIR}/src/moab/point_locater/lotte
    ${CMAKE_SOURCE_DIR}/itaps
    ${CMAKE_SOURCE_DIR}/itaps/imesh
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/parallel )

set( MBCOUPLER_SRCS Coupler.cpp
                    ElemUtil.cpp )

set( MBCOUPLER_INCLUDE_SRCS Coupler.hpp
                            ElemUtil.hpp )

set_source_files_properties( ${MBCOUPLER_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES}" )

add_library( mbcoupler
  ${MBCOUPLER_SRCS}
  ${MBCOUPLER_INCLUDE_SRCS} )

target_link_libraries( mbcoupler MOAB iMesh )
if ( MOAB_USE_MPI )
  target_link_libraries( mbcoupler MOABpar )
endif ( )

if ( MOAB_USE_NETCDF )
  set(netcdf_progs update_coords)
endif()

install(TARGETS mbcoupler RUNTIME DESTINATION bin
                          LIBRARY DESTINATION lib
                          ARCHIVE DESTINATION lib )
install(FILES ${MBCOUPLER_INCLUDE_SRCS} DESTINATION includes)

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
