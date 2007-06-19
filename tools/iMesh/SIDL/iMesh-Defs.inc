include ${ITAPS_SIDL_DIR}/Itaps-Defs.inc

#=========================================
# specify implementation includes/libs here

iMesh_IMPL_DIR=${HOME}/MOAB_gcc4.2
include ${iMesh_IMPL_DIR}/lib/moab.make
include ${iMesh_IMPL_DIR}/lib/iMesh-Defs.inc
iMesh_IMPL_LIBS = ${IMESH_LIBS_LINK}
iMesh_IMPL_INCLUDES = ${IMESH_INCLUDES}
#=========================================

include ${iBase_DIR}/iBase-Defs.inc

iMesh_SIDL_SERVER_LIBS = \
      -L${iMesh_DIR} -liMeshserver -Wl,-rpath,${iMesh_DIR} \
      ${iBase_SIDL_SERVER_LIBS} \
      ${iMesh_IMPL_LIBS}

iMesh_SIDL_SERVER_FILES = \
      ${iMesh_DIR}/libiMeshserver.so ${iBase_SIDL_SERVER_FILES}

iMesh_SIDL_SERVER_INCLUDES = \
      -I${iMesh_DIR}/mserver ${iBase_SIDL_SERVER_INCLUDES} ${iMesh_IMPL_INCLUDES} 

iMesh_SIDL_CCLIENT_INCLUDES = -I${iMesh_DIR}/Cclient \
     ${iBase_SIDL_CCLIENT_INCLUDES}
iMesh_SIDL_CCLIENT_FILES = libiMeshCclient.so ${iBase_SIDL_CCLIENT_FILES}
iMesh_SIDL_CCLIENT_LIBS = -L${iMesh_DIR} -liMeshCclient
