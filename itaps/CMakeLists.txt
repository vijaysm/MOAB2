include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set( BASE_INCLUDE_HEADERS iBase.h iBase_f.h)
option ( ENABLE_IMESH        "Should build IMESH?"       OFF )
option ( ENABLE_IGEOM        "Should build IGEOM?"       OFF )

if(ENABLE_IMESH )
  add_subdirectory(imesh)
endif()

if(ENABLE_IGEOM )
  add_subdirectory(igeom)
endif()

install(FILES ${BASE_INCLUDE_HEADERS} DESTINATION include)
