
# Generate Fortran name mangling headers
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
  COMMAND ${CMAKE_COMMAND}
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMesh.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
          -Dprefix:STRING=iMesh
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh.h
  )
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
  COMMAND ${CMAKE_COMMAND}
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
          -Dprefix:STRING=iMesh
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h
  )
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_extensions_protos.h
  COMMAND ${CMAKE_COMMAND}
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMeshP_extensions.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMeshP_extensions_protos.h
          -Dprefix:STRING=iMesh
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h
  )
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h
  COMMAND ${CMAKE_COMMAND}
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMeshP.h
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h
          -Dprefix:STRING=iMeshP
          -P ${CMAKE_SOURCE_DIR}/config/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMeshP.h
  )

set ( MOAB_IMESH_SRCS
      iMesh_MOAB.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
      ${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
      ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_extensions_protos.h
      ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h )

include_directories(
    ${MOAB_BINARY_DIR}
    ${MOAB_BINARY_DIR}/src
    ${MOAB_BINARY_DIR}/src/parallel
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/itaps
    ${MOAB_SOURCE_DIR}/itaps/imesh
    ${MOAB_BINARY_DIR}/itaps/imesh
  )

if ( MOAB_USE_HDF5 AND HDF5_FOUND )
  include_directories(
    ${HDF5_INCLUDE_DIRECTORIES}
    ${MOAB_SOURCE_DIR}/src/io/mhdf/include
  )
endif ( MOAB_USE_HDF5 AND HDF5_FOUND )

if ( MOAB_USE_MPI AND MPI_FOUND )
  LIST ( APPEND MOAB_IMESH_SRCS
         iMeshP_MOAB.cpp
         ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h )

  include_directories(
    ${MOAB_SOURCE_DIR}/src/parallel
  )
endif ( MOAB_USE_MPI AND MPI_FOUND )

set_source_files_properties( ${MOAB_IMESH_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES}"
)

add_library( iMesh
  ${MOAB_IMESH_SRCS}
  )

target_link_libraries( iMesh MOAB )
