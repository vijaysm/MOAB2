set ( MOAB_REFINER_SRCS 
  MBEdgeSizeEvaluator.cpp
  MBEdgeSizeSimpleImplicit.cpp
  MBEntityRefiner.cpp
  MBMeshRefiner.cpp
  )

add_library( MOABrefiner
  ${MOAB_REFINER_SRCS}
  )

