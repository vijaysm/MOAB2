set ( MOAB_REFINER_SRCS 
  MBEdgeSizeEvaluator.cpp
  MBEdgeSizeSimpleImplicit.cpp
  MBEntityRefiner.cpp
  MBMeshRefiner.cpp
  MBRefinerTagManager.cpp
  MBSimplexTemplateRefiner.cpp
  MBSimplexTemplateTagAssigner.cpp
  )

add_library( MOABrefiner
  ${MOAB_REFINER_SRCS}
  )

enable_testing()

add_executable( test_mesh_refiner test_mesh_refiner.cpp )
set_source_files_properties( test_mesh_refiner.cpp
  COMPILE_FLAGS "-DTEST ${MOAB_DEFINES}" )
target_link_libraries( test_mesh_refiner MOAB MOABrefiner )
add_test( TestMeshRefiner ${EXECUTABLE_OUTPUT_PATH}/test_mesh_refiner )

