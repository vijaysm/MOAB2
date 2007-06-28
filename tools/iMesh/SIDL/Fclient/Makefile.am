INCLUDES += -I@BABEL_DIR@/include

BUILT_SOURCES = $(libiMeshF77client_la_SOURCES) $(libiMeshF77client_la_include_HEADERS)

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

lib_LTLIBRARIES = libiMeshF77client.la

libiMeshF77Client_la_LIBADD = ../../libiMesh.la

libiMeshF77client_la_includedir = $(includedir)/SIDL

libiMeshF77client_la_SOURCES = \
  iBase_ArrTag_IOR.h iBase_CreationStatus_IOR.h iBase_EntSet_IOR.h    \
  iBase_EntTag_IOR.h iBase_EntityType_IOR.h iBase_ErrorActions_IOR.h          \
  iBase_ErrorType_IOR.h iBase_Error_IOR.h iBase_IOR.h iBase_SetBoolOps_IOR.h  \
  iBase_SetRelation_IOR.h iBase_SetTag_IOR.h iBase_StorageOrder_IOR.h         \
  iBase_TagValueType_IOR.h iBase_Tag_IOR.h iMesh_AdjacencyInfo_IOR.h          \
  iMesh_ArrMod_IOR.h iMesh_Arr_IOR.h iMesh_EntityTopology_IOR.h               \
  iMesh_Entity_IOR.h iMesh_IOR.h iMesh_Mesh_IOR.h iMesh_Modify_IOR.h          \
  iMesh_SIDL_IOR.h iMesh_SIDL_MeshSidl_IOR.h \
  iBase_ArrTag_fStub.c iBase_CreationStatus_fStub.c                  \
  iBase_EntSet_fStub.c iBase_EntTag_fStub.c iBase_EntityType_fStub.c          \
  iBase_ErrorActions_fStub.c iBase_ErrorType_fStub.c iBase_Error_fStub.c      \
  iBase_SetBoolOps_fStub.c iBase_SetRelation_fStub.c iBase_SetTag_fStub.c     \
  iBase_StorageOrder_fStub.c iBase_TagValueType_fStub.c iBase_Tag_fStub.c     \
  iMesh_AdjacencyInfo_fStub.c iMesh_ArrMod_fStub.c iMesh_Arr_fStub.c          \
  iMesh_EntityTopology_fStub.c iMesh_Entity_fStub.c iMesh_Mesh_fStub.c        \
  iMesh_Modify_fStub.c iMesh_SIDL_MeshSidl_fStub.c sidl_BaseClass_fStub.c     \
  sidl_BaseException_fStub.c sidl_BaseInterface_fStub.c                       \
  sidl_ClassInfoI_fStub.c sidl_ClassInfo_fStub.c sidl_DFinder_fStub.c         \
  sidl_DLL_fStub.c sidl_Finder_fStub.c sidl_InvViolation_fStub.c              \
  sidl_Loader_fStub.c sidl_PostViolation_fStub.c sidl_PreViolation_fStub.c    \
  sidl_Resolve_fStub.c sidl_SIDLException_fStub.c sidl_Scope_fStub.c          \
  sidl_array_fStub.c sidl_bool_fStub.c sidl_char_fStub.c                      \
  sidl_dcomplex_fStub.c sidl_double_fStub.c sidl_fcomplex_fStub.c             \
  sidl_float_fStub.c sidl_int_fStub.c sidl_io_Deserializer_fStub.c            \
  sidl_io_IOException_fStub.c sidl_io_Serializeable_fStub.c                   \
  sidl_io_Serializer_fStub.c sidl_long_fStub.c sidl_opaque_fStub.c            \
  sidl_rmi_ConnectRegistry_fStub.c sidl_rmi_InstanceHandle_fStub.c            \
  sidl_rmi_InstanceRegistry_fStub.c sidl_rmi_Invocation_fStub.c               \
  sidl_rmi_NetworkException_fStub.c sidl_rmi_ProtocolFactory_fStub.c          \
  sidl_rmi_Response_fStub.c sidl_string_fStub.c


libiMeshF77client_la_include_HEADERS = \
  iBase_CreationStatus.inc iBase_EntityType.inc                \
  iBase_ErrorActions.inc iBase_ErrorType.inc iBase_StorageOrder.inc           \
  iBase_TagValueType.inc iMesh_AdjacencyInfo.inc iMesh_EntityTopology.inc     \
  sidl_Resolve.inc sidl_Scope.inc 

dist-hook:
	for i in $(BUILT_SOURCES); do $(RM) $(distdir)/$$i; done
