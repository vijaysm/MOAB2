include_directories(  ${CMAKE_SOURCE_DIR}/src/io/mhdf/include )
set( TESTS h5file_test.cpp
           h5legacy.cpp
           h5varlen.cpp
           h5sets_test.cpp
           h5regression.cpp
           h5partial.cpp
           h5portable.cpp )

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR} ${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( ${base} MOAB )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  set_tests_properties ( ${base} PROPERTIES TIMEOUT 120)
endforeach()

add_executable( dump_sets dump_sets.c )
set_target_properties( dump_sets PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
target_link_libraries( dump_sets MOAB ${MOAB_LIBS} )

