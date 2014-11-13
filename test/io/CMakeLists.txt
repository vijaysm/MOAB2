include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/test
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/parallel
)
  
set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

set( TESTS readutil_test.cpp VtkTest.cpp cub_file_test.cc gmsh_test.cc nastran_test.cc ideas_test.cc stl_test.cc smf_test.cc)
if( MOAB_USE_NETCDF )
  set( TESTS ${TESTS} exodus_test.cc )
  set( TESTS ${TESTS} read_nc.cpp )
  set( TESTS ${TESTS} read_ucd_nc.cpp )
  set( TESTS ${TESTS} read_mpas_nc.cpp )
endif()

if(MOAB_USE_CCMIO)
  set( TESTS ${TESTS} ccmio_test.cpp )
endif()

if( MOAB_USE_CGNS )
  set( TESTS ${TESTS} cgns_test.cpp )
endif()

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES}" )
  target_link_libraries( ${base} MOAB )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()


add_executable( tqdcfr tqdcfr_test.cpp )
set_source_files_properties( tqdcfr_test.cpp
    COMPILE_FLAGS "${MOAB_DEFINES}" )
target_link_libraries( tqdcfr MOAB )
add_test( TestTQDCFR ${EXECUTABLE_OUTPUT_PATH}/tqdcfr )
