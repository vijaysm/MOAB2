INCLUDES += -I@BABEL_DIR@/include

BUILT_SOURCES =  $(libiMeshCclient_la_SOURCES) $(libiMeshCclient_la_include_HEADERS)

if USE_BABEL
  BABEL_CMD=@BABEL_DIR@/bin/babel -sC++ $(SIDL_FILES)
else
  BABEL_CMD=touch $(BUILT_SOURCES)
endif

BABEL=@BABEL_DIR@/bin/babel
SIDL_FILES = @srcdir@/../iBase.sidl \
             @srcdir@/../iMesh.sidl \
	     @srcdir@/../iMesh_SIDL.sidl

$(BUILT_SOURCES) : .timestamp

.timestamp: $(SIDL_FILES)
	$(BABEL_CMD)
	touch $@

lib_LTLIBRARIES = libiMeshCclient.la

libiMeshCclient_la_LIBADD = ../../libiMesh.la

libiMeshCclient_la_includedir = $(includedir)/SIDL

libiMeshCclient_la_SOURCES = \
  iBase_ArrTag_IOR.h iBase_CreationStatus_IOR.h iBase_EntSet_IOR.h    \
  iBase_EntTag_IOR.h iBase_EntityType_IOR.h iBase_ErrorActions_IOR.h          \
  iBase_ErrorType_IOR.h iBase_Error_IOR.h iBase_IOR.h iBase_SetBoolOps_IOR.h  \
  iBase_SetRelation_IOR.h iBase_SetTag_IOR.h iBase_StorageOrder_IOR.h         \
  iBase_TagValueType_IOR.h iBase_Tag_IOR.h iMesh_AdjacencyInfo_IOR.h          \
  iMesh_ArrMod_IOR.h iMesh_Arr_IOR.h iMesh_EntityTopology_IOR.h               \
  iMesh_Entity_IOR.h iMesh_IOR.h iMesh_Mesh_IOR.h iMesh_Modify_IOR.h          \
  iMesh_SIDL_IOR.h iMesh_SIDL_MeshSidl_IOR.h \
  iBase_ArrTag_Stub.c iBase_CreationStatus_Stub.c                    \
  iBase_EntSet_Stub.c iBase_EntTag_Stub.c iBase_EntityType_Stub.c             \
  iBase_ErrorActions_Stub.c iBase_ErrorType_Stub.c iBase_Error_Stub.c         \
  iBase_SetBoolOps_Stub.c iBase_SetRelation_Stub.c iBase_SetTag_Stub.c        \
  iBase_StorageOrder_Stub.c iBase_TagValueType_Stub.c iBase_Tag_Stub.c        \
  iMesh_AdjacencyInfo_Stub.c iMesh_ArrMod_Stub.c iMesh_Arr_Stub.c             \
  iMesh_EntityTopology_Stub.c iMesh_Entity_Stub.c iMesh_Mesh_Stub.c           \
  iMesh_Modify_Stub.c iMesh_SIDL_MeshSidl_Stub.c


libiMeshCclient_la_include_HEADERS = \
  iBase.h iBase_ArrTag.h iBase_CreationStatus.h iBase_EntSet.h       \
  iBase_EntTag.h iBase_EntityType.h iBase_Error.h iBase_ErrorActions.h        \
  iBase_ErrorType.h iBase_SetBoolOps.h iBase_SetRelation.h iBase_SetTag.h     \
  iBase_StorageOrder.h iBase_Tag.h iBase_TagValueType.h iMesh.h               \
  iMesh_AdjacencyInfo.h iMesh_Arr.h iMesh_ArrMod.h iMesh_Entity.h             \
  iMesh_EntityTopology.h iMesh_Mesh.h iMesh_Modify.h iMesh_SIDL.h             \
  iMesh_SIDL_MeshSidl.h


dist-hook:
	for i in $(BUILT_SOURCES); do $(RM) $(distdir)/$$i; done

