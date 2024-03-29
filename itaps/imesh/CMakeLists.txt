
# Generate Fortran name mangling headers
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
  COMMAND ${CMAKE_COMMAND}
          -DDEF:STRING=MOAB
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMesh.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
          -Dprefix:STRING=iMesh
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh.h )
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
  COMMAND ${CMAKE_COMMAND}
          -DDEF:STRING=MOAB
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
          -Dprefix:STRING=iMesh
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h )

set ( MOAB_IMESH_SRCS
          iMesh_MOAB.hpp iMesh_MOAB.cpp
          MBIter.hpp )
set ( MOAB_IMESH_LIB_INCLUDES
          iMesh.h
          iMesh_extensions.h
          ${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
          iMesh_f.h
          ${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
          MBiMesh.hpp )

if ( MOAB_HAVE_MPI )
  include_directories( ${CMAKE_SOURCE_DIR}/src/parallel
                       ${CMAKE_BINARY_DIR}/src/parallel)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_extensions_protos.h
    COMMAND ${CMAKE_COMMAND}
            -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMeshP_extensions.h
            -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMeshP_extensions_protos.h
            -Dprefix:STRING=iMesh
            -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h
    COMMAND ${CMAKE_COMMAND}
            -DDEF:STRING=MOAB
            -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMeshP.h
            -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h
            -Dprefix:STRING=iMeshP
            -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMeshP.h )
  set(MOAB_IMESH_SRCS ${MOAB_IMESH_SRCS}
                      iMeshP_MOAB.cpp)
  set(MOAB_IMESH_LIB_INCLUDES
          ${MOAB_IMESH_LIB_INCLUDES}
          ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h
          ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_extensions_protos.h )
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/itaps
  ${CMAKE_BINARY_DIR}/itaps
  )

if ( MOAB_HAVE_HDF5 AND HDF5_FOUND )
  include_directories( ${HDF5_INCLUDE_DIRECTORIES}
                       ${CMAKE_SOURCE_DIR}/src/io/mhdf/include )
endif( )

set_source_files_properties(
  ${MOAB_IMESH_SRCS} COMPILE_FLAGS "${MOAB_DEFINES}"
  )

add_library(
  iMesh
  ${MOAB_IMESH_SRCS}
  ${MOAB_IMESH_LIB_INCLUDES}
  )

target_link_libraries( iMesh MOAB )

set_target_properties(
  iMesh
  PROPERTIES
  VERSION ${MOAB_VERSION}
  SOVERSION ${MOAB_VERSION_MAJOR}
  )

install(
  TARGETS iMesh
  EXPORT MOABTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
set( MOAB_INSTALL_TARGETS ${MOAB_INSTALL_TARGETS} iMesh)

install(
  FILES ${MOAB_IMESH_LIB_INCLUDES}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
