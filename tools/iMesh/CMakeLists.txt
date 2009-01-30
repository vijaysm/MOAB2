
# Generate Fortran name mangling headers
ADD_CUSTOM_COMMAND( 
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
  COMMAND ${CMAKE_COMMAND} 
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMesh.h 
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h 
          -Dprefix:STRING=iMesh 
          -P ${CMAKE_SOURCE_DIR}/cmake/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh.h
  )
ADD_CUSTOM_COMMAND( 
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
  COMMAND ${CMAKE_COMMAND} 
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h 
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h 
          -Dprefix:STRING=iMesh 
          -P ${CMAKE_SOURCE_DIR}/cmake/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMesh_extensions.h
  )
ADD_CUSTOM_COMMAND( 
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h
  COMMAND ${CMAKE_COMMAND} 
          -Dinput_file:STRING=${CMAKE_CURRENT_SOURCE_DIR}/iMeshP.h 
          -Doutput_file:STRING=${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h 
          -Dprefix:STRING=iMeshP 
          -P ${CMAKE_SOURCE_DIR}/cmake/ITAPSFortranMangling.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/iMeshP.h
  )

set ( MOAB_IMESH_SRCS 
      iMesh_MOAB.cpp
      ${CMAKE_CURRENT_BINARY_DIR}/iMesh_protos.h
      ${CMAKE_CURRENT_BINARY_DIR}/iMesh_extensions_protos.h
      ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h )

if ( MOAB_USE_MPI AND MPI_FOUND )
  LIST ( APPEND MOAB_IMESH_SRCS 
         iMeshP_MOAB.cpp
         ${CMAKE_CURRENT_BINARY_DIR}/iMeshP_protos.h )
endif ( MOAB_USE_MPI AND MPI_FOUND )
         

include_directories( ${CMAKE_SOURCE_DIR}
                     ${CMAKE_BINARY_DIR}
                     ${CMAKE_SOURCE_DIR}/parallel
                     ${CMAKE_BINARY_DIR}/parallel
                     ${CMAKE_CURRENT_SOURCE_DIR}
                     ${CMAKE_CURRENT_BINARY_DIR}
                    )
                     

set_source_files_properties(
  ${MOAB_IMESH_SRCS}
  COMPILE_FLAGS "${MOAB_DEFINES}"
)

add_library( iMesh
  ${MOAB_IMESH_SRCS}
  )
