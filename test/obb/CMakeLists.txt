set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest -DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR}")

set( TESTS obb_test.cpp )

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
  target_link_libraries( ${base} MOAB ${CGM_LIBRARIES} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

add_executable( obb_time obb_time.cpp)
set_target_properties( obb_time PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES}" )
target_link_libraries( obb_time MOAB ${CGM_LIBRARIES})

add_executable( obb_tree_tool obb_tree_tool.cpp)
set_target_properties( obb_tree_tool PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES}" )
target_link_libraries( obb_tree_tool MOAB ${CGM_LIBRARIES} )
