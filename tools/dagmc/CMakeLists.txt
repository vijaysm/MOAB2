# Find CGM directory
find_package( CGM )

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MOAB_BINARY_DIR}/src
  ${MOAB_SOURCE_DIR}/src
  ${MOAB_BINARY_DIR}/src/moab
  ${MOAB_SOURCE_DIR}/src/moab
)

set( MOAB_DAGMC_LIB_SRCS
  DagMC.cpp
)

set( MOAB_TEST_GEOM_SRCS
  test_geom.cc
)

if ( CGM_FOUND )
  set( MOAB_DAGMC_LIB_SRCS
    ${MOAB_DAGMC_LIB_SRCS}
  )
  set( MOAB_DAGMC_SRCS
    main.cc
    cubfile.c
  )
  include_directories(
    ${CGM_INCLUDE_DIRECTORIES}
  )
endif ( CGM_FOUND )

add_library( dagmc
  ${MOAB_DAGMC_LIB_SRCS}
)
target_link_libraries( dagmc
  MOAB
)

add_executable( test_geom
  ${MOAB_TEST_GEOM_SRCS}
)
target_link_libraries( test_geom
  dagmc
  MOAB
)

if ( CGM_FOUND )
  target_link_libraries( dagmc
    ${CGM_LIBRARIES}
  )
  target_link_libraries( test_geom
    ${CGM_LIBRARIES}
  )
  add_executable( cgm2moab
    ${MOAB_DAGMC_SRCS}
  )
  target_link_libraries( cgm2moab
    dagmc
    ${CGM_LIBRARIES}
  )
endif ( CGM_FOUND )
