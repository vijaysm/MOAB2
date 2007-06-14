include Itaps-Defs.inc

#=========================================
# specify implementation includes/libs here

MOAB_DIR=${HOME}/MOAB
include $(MOAB_DIR)/lib/moab.make
iMesh_IMPLEMENTATION_DIR = ${MOAB_LIBDIR}
iMesh_IMPLEMENTATION_LIBS = \
      ${MOAB_LIBS_LINK} -Wl,-rpath,${MOAB_LIBDIR}
#=========================================

include ${iBase_DIR}/iBase-Defs.inc

iMesh_CBIND_SERVER_INCLUDES = \
      -I${iMesh_DIR} ${iBase_CBIND_SERVER_INCLUDES}

iMesh_CBIND_SERVER_LIBS = \
      -L${iMesh_IMPLEMENTATION_DIR} -liMesh ${iMesh_IMPLEMENTATION_LIBS}

iMesh_CBIND_SERVER_FILES = \
      ${iMesh_IMPLEMENTATION_DIR}/libiMesh.a ${iMesh_DIR}/iMesh.h


iMesh_SERVER_INCLUDES = \
      -I${iMesh_DIR}/mserver ${iBase_SERVER_INCLUDES} ${iMesh_CBIND_SERVER_INCLUDES} 

iMesh_SERVER_LIBS = \
      -L${iMesh_DIR} -liMeshserver -Wl,-rpath,${iMesh_DIR} \
      ${iBase_SERVER_LIBS} ${iMesh_CBIND_SERVER_LIBS}

iMesh_SERVER_FILES = \
      ${iMesh_DIR}/libiMeshserver.so ${iMesh_CBIND_SERVER_FILES}

iMesh_CLIENT_C_INCLUDES = -I${iMesh_DIR}/Cclient 
iMesh_CLIENT_C_FILES = libiMeshCclient.so ${iMesh_SERVER_FILES}
iMesh_CLIENT_C_LIBS = -L${iMesh_DIR} -liMeshCclient ${iMesh_SERVER_LIBS}
