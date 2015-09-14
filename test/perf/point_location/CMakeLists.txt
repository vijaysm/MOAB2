set( LIBS MOAB ${CGM_LIBRARIES})
set( TESTS point_location.cpp
           tree_searching_perf.cpp
           sploc_searching_perf.cpp )

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES} -DIS_BUILDING_MB" )
  target_link_libraries( ${base} ${LIBS} )
endforeach()
