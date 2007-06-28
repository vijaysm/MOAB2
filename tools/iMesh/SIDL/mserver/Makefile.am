INCLUDES += -I@BABEL_DIR@/include

BUILT_SOURCES = $(libiMeshserver_la_SOURCES) $(libiMeshserver_la_include_HEADERS)

if USE_BABEL
  BABEL_CMD=@BABEL_DIR@/bin/babel -sC++ $(SIDL_FILES)
else
  BABEL_CMD=touch $(BUILT_SOURCES)
endif

SIDL_FILES = @srcdir@/../iBase.sidl \
             @srcdir@/../iMesh.sidl \
	     @srcdir@/../iMesh_SIDL.sidl

$(BUILT_SOURCES) : .timestamp

.timestamp: $(SIDL_FILES)
	$(BABEL_CMD)
	touch $@

lib_LTLIBRARIES = libiMeshserver.la

libiMeshserver_la_LIBADD = ../../libiMesh.la

libiMeshserver_la_includedir = $(includedir)/SIDL
 
libiMeshserver_la_SOURCES = \
  iBase_Error_Impl.hh iMesh_SIDL_MeshSidl_Impl.hh \
  iBase_Error_Impl.cc iMesh_SIDL_MeshSidl_Impl.cc \
  iBase_ArrTag_IOR.h iBase_CreationStatus_IOR.h iBase_EntSet_IOR.h    \
  iBase_EntTag_IOR.h iBase_EntityType_IOR.h iBase_ErrorActions_IOR.h          \
  iBase_ErrorType_IOR.h iBase_Error_IOR.h iBase_IOR.h iBase_SetBoolOps_IOR.h  \
  iBase_SetRelation_IOR.h iBase_SetTag_IOR.h iBase_StorageOrder_IOR.h         \
  iBase_TagValueType_IOR.h iBase_Tag_IOR.h iMesh_AdjacencyInfo_IOR.h          \
  iMesh_ArrMod_IOR.h iMesh_Arr_IOR.h iMesh_EntityTopology_IOR.h               \
  iMesh_Entity_IOR.h iMesh_IOR.h iMesh_Mesh_IOR.h iMesh_Modify_IOR.h          \
  iMesh_SIDL_IOR.h iMesh_SIDL_MeshSidl_IOR.h \
  iBase_Error_IOR.c iMesh_SIDL_MeshSidl_IOR.c \
  iBase_Error_Skel.cc iMesh_SIDL_MeshSidl_Skel.cc \
  iBase_ArrTag.cc iBase_EntSet.cc iBase_EntTag.cc iBase_Error.cc     \
  iBase_SetBoolOps.cc iBase_SetRelation.cc iBase_SetTag.cc iBase_Tag.cc       \
  iMesh_Arr.cc iMesh_ArrMod.cc iMesh_Entity.cc iMesh_Mesh.cc iMesh_Modify.cc  \
  iMesh_SIDL_MeshSidl.cc sidl_BaseClass.cc sidl_BaseException.cc              \
  sidl_BaseInterface.cc sidl_ClassInfo.cc sidl_ClassInfoI.cc sidl_DFinder.cc  \
  sidl_DLL.cc sidl_Finder.cc sidl_InvViolation.cc sidl_Loader.cc              \
  sidl_PostViolation.cc sidl_PreViolation.cc sidl_SIDLException.cc            \
  sidl_io_Deserializer.cc sidl_io_IOException.cc sidl_io_Serializeable.cc     \
  sidl_io_Serializer.cc sidl_rmi_ConnectRegistry.cc                           \
  sidl_rmi_InstanceHandle.cc sidl_rmi_InstanceRegistry.cc                     \
  sidl_rmi_Invocation.cc sidl_rmi_NetworkException.cc                         \
  sidl_rmi_ProtocolFactory.cc sidl_rmi_Response.cc

libiMeshserver_la_include_HEADERS = \
  iBase.hh iBase_ArrTag.hh iBase_CreationStatus.hh iBase_EntSet.hh   \
  iBase_EntTag.hh iBase_EntityType.hh iBase_Error.hh iBase_ErrorActions.hh    \
  iBase_ErrorType.hh iBase_SetBoolOps.hh iBase_SetRelation.hh iBase_SetTag.hh \
  iBase_StorageOrder.hh iBase_Tag.hh iBase_TagValueType.hh iMesh.hh           \
  iMesh_AdjacencyInfo.hh iMesh_Arr.hh iMesh_ArrMod.hh iMesh_Entity.hh         \
  iMesh_EntityTopology.hh iMesh_Mesh.hh iMesh_Modify.hh iMesh_SIDL.hh         \
  iMesh_SIDL_MeshSidl.hh sidl.hh sidl_BaseClass.hh sidl_BaseException.hh      \
  sidl_BaseInterface.hh sidl_ClassInfo.hh sidl_ClassInfoI.hh sidl_DFinder.hh  \
  sidl_DLL.hh sidl_Finder.hh sidl_InvViolation.hh sidl_Loader.hh              \
  sidl_PostViolation.hh sidl_PreViolation.hh sidl_Resolve.hh                  \
  sidl_SIDLException.hh sidl_Scope.hh sidl_io.hh sidl_io_Deserializer.hh      \
  sidl_io_IOException.hh sidl_io_Serializeable.hh sidl_io_Serializer.hh       \
  sidl_rmi.hh sidl_rmi_ConnectRegistry.hh sidl_rmi_InstanceHandle.hh          \
  sidl_rmi_InstanceRegistry.hh sidl_rmi_Invocation.hh                         \
  sidl_rmi_NetworkException.hh sidl_rmi_ProtocolFactory.hh                    \
  sidl_rmi_Response.hh

my_dist_sources = \
  iBase_Error_Impl.hh iMesh_SIDL_MeshSidl_Impl.hh \
  iBase_Error_Impl.cc iMesh_SIDL_MeshSidl_Impl.cc


dist-hook:
	for i in $(BUILT_SOURCES); do \
	  found=0; \
          for j in $(my_dist_sources); do \
	    if test "x$$i" = "x$$j"; then \
		found=1; \
	    fi; \
	  done; \
	  if test "$$found" = "0"; then \
	    $(RM) $(distdir)/$$i; \
	  fi; \
	done
