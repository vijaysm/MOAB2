set ( MOAB_REFINER_SRCS 
  MBEdgeSizeEvaluator.cpp
  MBEdgeSizeSimpleImplicit.cpp
  MBEntityRefiner.cpp
  MBMeshOutputFunctor.cpp
  MBMeshRefiner.cpp
  MBProcessSet.cpp
  MBRefinerTagManager.cpp
  MBSimplexTemplateRefiner.cpp
  MBSimplexTemplateTagAssigner.cpp
  MBSplitVertices.cpp
  )

set_source_files_properties(
  ${MOAB_REFINER_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES}"
)

add_library( MOABrefiner
  ${MOAB_REFINER_SRCS}
  )

enable_testing()

add_executable( test_mesh_refiner test_mesh_refiner.cpp )
set_source_files_properties( test_mesh_refiner.cpp
  COMPILE_FLAGS "-DTEST ${MOAB_DEFINES}" )
target_link_libraries( test_mesh_refiner MOAB MOABrefiner )
add_test( TestMeshRefiner
  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
  ${EXECUTABLE_OUTPUT_PATH}/test_mesh_refiner ${MPIEXEC_POSTFLAGS}
  ${MOAB_SOURCE_DIR}/refiner/fourVolsBare.cub )

