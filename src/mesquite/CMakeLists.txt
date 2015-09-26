include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_BINARY_DIR}
                    parallel
                    LocalDiscretization
                    io io/mhdf
                    ${MOAB_BINARY_DIR}/src/parallel )

# Include the makefile to get all headers and sources
include (${CMAKE_CURRENT_SOURCE_DIR}/make.inc)
#add_subdirectory(Control)
#add_subdirectory(MappingFunction)
#add_subdirectory(Mesh)
#add_subdirectory(Misc)
#add_subdirectory(ObjectiveFunction)
#add_subdirectory(QualityAssessor)
#add_subdirectory(QualityImprover)
#add_subdirectory(QualityMetric)
#add_subdirectory(TargetCalculator)
#add_subdirectory(TargetMetric)
#add_subdirectory(Wrappers)

#generate the Version.h file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mesquite_version.h.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_version.h" @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mesquite_config.h.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_config.h" @ONLY)

set(MESQUITE_LIB_SOURCES ${MSQ_SRCS} )

set( MESQUITE_INSTALL_HEADERS ${MSQ_HDRS} )

set( MESQUITE_GENERATED_INSTALL_HEADERS
     ${CMAKE_CURRENT_BINARY_DIR}/include/mesquite_config.h
     ${CMAKE_CURRENT_BINARY_DIR}/mesquite_version.h
  )

add_library( moabmesquite OBJECT ${MESQUITE_LIB_SRCS} ${MESQUITE_INSTALL_HEADERS} ${MESQUITE_GENERATED_INSTALL_HEADERS} )
set_target_properties( moabmesquite PROPERTIES COMPILE_FLAGS "${MESQUITE_DEFINES} ${MOAB_DEFINES}")
include(AutoconfHeader)
moab_install_headers(${MESQUITE_INSTALL_HEADERS})

