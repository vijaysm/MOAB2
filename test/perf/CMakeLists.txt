add_subdirectory(point_location)
set( LIBS MOAB )
set( TESTS perf.cpp
           seqperf.cpp
           adj_time.cpp
           perftool.cpp )
if(MOAB_BUILD_IMESH)
  set(TESTS ${TESTS} tstt_perf_binding.cpp)
  set(LIBS ${LIBS} iMesh)
endif()

set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" tmp ${fname} )
  string( REPLACE ".cc" "" base ${tmp} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS} -DIS_BUILDING_MB" )
  target_link_libraries( ${base} ${LIBS} )
endforeach()
