include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_BINARY_DIR}/src/parallel
    ${Zoltan_INCLUDES}
)

if ( MOAB_USE_ZOLTAN )
  # Build the mbpart executable
  add_executable( mbpart mbpart.cpp )
  set_source_files_properties( mbpart.cpp
    COMPILE_FLAGS "${MOAB_DEFINES} -DMESHDIR=${MOAB_SOURCE_DIR}/MeshFiles/unittest" )
  target_link_libraries( mbpart MOAB ${Zoltan_LIBRARIES} )

endif ( MOAB_USE_ZOLTAN )
