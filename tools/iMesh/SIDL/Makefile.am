include @abs_srcdir@/iMesh-SIDL-Defs.inc

vpath = mserver bserver 
VPATH = mserver bserver 

LD_FLAGS_PIC = -fpic -shared
SERVER_HDRS := iBase_ArrTag_IOR.h iBase_CreationStatus_IOR.h iBase_EntSet_IOR.h    \
  iBase_EntTag_IOR.h iBase_EntityType_IOR.h iBase_ErrorActions_IOR.h          \
  iBase_ErrorType_IOR.h iBase_IOR.h iBase_SetBoolOps_IOR.h  \
  iBase_SetRelation_IOR.h iBase_SetTag_IOR.h iBase_StorageOrder_IOR.h         \
  iBase_TagValueType_IOR.h iBase_Tag_IOR.h iMesh_AdjacencyInfo_IOR.h          \
  iMesh_ArrMod_IOR.h iMesh_Arr_IOR.h iMesh_EntityTopology_IOR.h               \
  iMesh_Entity_IOR.h iMesh_IOR.h iMesh_Mesh_IOR.h iMesh_Modify_IOR.h          \
  iMesh_SIDL_IOR.h iMesh_SIDL_MeshSidl_IOR.h                                  \
  iBase.hh iBase_ArrTag.hh iBase_CreationStatus.hh iBase_EntSet.hh   \
  iBase_EntTag.hh iBase_EntityType.hh iBase_ErrorActions.hh    \
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
  sidl_rmi_Response.hh iMesh_SIDL_MeshSidl_Impl.hh iBase_Error_Impl.hh        \
  iBase_Error_IOR.h iBase_Error.hh 
SERVER_SRCS = iMesh_SIDL_MeshSidl_IOR.c                                   \
  iBase_ArrTag.cc iBase_EntSet.cc iBase_EntTag.cc                             \
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
  sidl_rmi_ProtocolFactory.cc sidl_rmi_Response.cc iMesh_SIDL_MeshSidl_Skel.cc \
  iMesh_SIDL_MeshSidl_Impl.cc iBase_Error_Impl.cc iBase_Error_IOR.c           \
  iBase_Error_Skel.cc iBase_Error.cc
SERVER_OFILES = iMesh_SIDL_MeshSidl_IOR.o                                   \
  iBase_ArrTag.o iBase_EntSet.o iBase_EntTag.o                             \
  iBase_SetBoolOps.o iBase_SetRelation.o iBase_SetTag.o iBase_Tag.o       \
  iMesh_Arr.o iMesh_ArrMod.o iMesh_Entity.o iMesh_Mesh.o iMesh_Modify.o  \
  iMesh_SIDL_MeshSidl.o sidl_BaseClass.o sidl_BaseException.o              \
  sidl_BaseInterface.o sidl_ClassInfo.o sidl_ClassInfoI.o sidl_DFinder.o  \
  sidl_DLL.o sidl_Finder.o sidl_InvViolation.o sidl_Loader.o              \
  sidl_PostViolation.o sidl_PreViolation.o sidl_SIDLException.o            \
  sidl_io_Deserializer.o sidl_io_IOException.o sidl_io_Serializeable.o     \
  sidl_io_Serializer.o sidl_rmi_ConnectRegistry.o                           \
  sidl_rmi_InstanceHandle.o sidl_rmi_InstanceRegistry.o                     \
  sidl_rmi_Invocation.o sidl_rmi_NetworkException.o                         \
  sidl_rmi_ProtocolFactory.o sidl_rmi_Response.o iMesh_SIDL_MeshSidl_Skel.o \
  iMesh_SIDL_MeshSidl_Impl.o iBase_Error_Impl.o iBase_Error_IOR.o           \
  iBase_Error_Skel.o iBase_Error.o

#testcxx: ${iMesh_SIDL_FILES} testcxx.o 
#	$(CXX)  ${LD_FLAGS} -o $@ testcxx.o -Wl,-rpath,. -L. ${iMesh_SIDL_LIBS}

bin_PROGRAMS = testcxx
testcxx_SOURCES = testcxx.cpp

testcxx.o: bserver/babel.make mserver/babel.make

LDADD = @top_srcdir@/tools/iMesh/libiMesh.la $(top_builddir)/libMOAB.la @abs_srcdir@/libiMeshserver.so ${BABEL_DIR}/lib/libsidl.la

testcxx_DEPENDENCIES = bserver/babel.make @top_srcdir@/tools/iMesh/libiMesh.la $(top_builddir)/libMOAB.la @abs_srcdir@/libiMeshserver.so ${BABEL_DIR}/lib/libsidl.la ${SERVER_OFILES} bserver/babel.make mserver/babel.make

settings:
	@echo "FACTORY_OFILES = $(FACTORY_OFILES)"
	@echo "SERVER_OFILES  = $(SERVER_OFILES)"

#============= iBase =============
repo/.btimestamp: iBase.sidl
	-rm -rf repo
	@BABEL_DIR@/bin/babel -tXML -o repo iBase.sidl 
	touch $@ 

bserver/babel.make: repo/.btimestamp iBase.sidl
	-rm -f $@
	@BABEL_DIR@/bin/babel -sC++ -o bserver iBase.sidl

#============= iBase =============

repo/.mtimestamp: iMesh.sidl 
	@BABEL_DIR@/bin/babel -tXML -o repo -Rrepo iMesh.sidl
	touch $@ 

mserver/babel.make: repo/.mtimestamp repo/.btimestamp iMesh_SIDL.sidl iMesh.sidl 
	-rm -f $@
	@BABEL_DIR@/bin/babel -R"repo" -sC++ -o mserver iMesh_SIDL.sidl

Cclient: repo/.mtimestamp repo/.mtimestamp iMesh_SIDL.sidl
	-rm -rf Cclient
	@BABEL_DIR@/bin/babel -R"repo" -cC -o Cclient iMesh_SIDL.sidl

F77client: repo/.mtimestamp repo/.btimestamp iMesh_SIDL.sidl
	-rm -rf F77client
	@BABEL_DIR@/bin/babel -R"repo" -cF77 -o F77client iMesh_SIDL.sidl

${iMesh_SIDLLIBDIR}/libiMeshserver.so: bserver/babel.make mserver/babel.make ${SERVER_OFILES}
	$(CXX) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${SERVER_OFILES}

libiMeshCclient.so: Cclient ${CCLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${CCLIENT_OFILES}

libiMeshF77client.so: F77client ${F77CLIENT_OFILES}
	$(LD) $(LD_FLAGS) $(LD_FLAGS_PIC) -o $@ ${F77CLIENT_OFILES}

Cclient/babel.make: Cclient

F77client/babel.make: F77client

clean:
	-rm -f *.o *.so testc testcxx *~ core.*

clean_all: clean
	-rm -rf F77client Cclient repo
	cd mserver; rm -f ${SERVER_HDRS} ${SERVER_SRCS}

%.o : %.cxx
	$(CXX) $(CXX_FLAGS) ${iMesh_SIDL_INCLUDES} -o $@ -c $<

%.o : %.cpp
	$(CXX) $(CXX_FLAGS) ${iMesh_SIDL_INCLUDES} -o $@ -c $<

%.o : %.cc
	$(CXX) $(CXX_FLAGS) ${iMesh_SIDL_INCLUDES} -o $@ -c $<

%.o : %.c
	$(CC) $(CC_FLAGS) ${iMesh_SIDL_INCLUDES} -o $@ -c $<

%.o : %.f
	$(FC) $(FC_FLAGS) ${iMesh_SIDL_INCLUDES} -o $@ -c $<

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put 
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = iMesh-SIDL-Defs.inc

# By default, iMesh-SIDL-Defs.inc will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does 
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	echo "iMesh_DIR=${DESTDIR}${cfgdir}/.." >> $(DESTDIR)$(cfgdir)/iMesh-SIDL-Defs.inc
	echo "iMesh_SIDLLIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/iMesh-SIDL-Defs.inc
	echo "iMesh_SIDLINCLUDEDIR=${includedir}/SIDL" >> $(DESTDIR)$(cfgdir)/iMesh-SIDL-Defs.inc
	(cd mserver; install -d ${includedir}/SIDL; install -t ${includedir}/SIDL ${SERVER_HDRS})

#include make.dependencies
iMesh_SIDL_MeshSidl_IOR.o: mserver/iMesh_SIDL_MeshSidl_IOR.c
iBase_ArrTag.o: mserver/iBase_ArrTag.cc
iBase_EntSet.o: mserver/iBase_EntSet.cc
iBase_EntTag.o: mserver/iBase_EntTag.cc
iBase_SetBoolOps.o: mserver/iBase_SetBoolOps.cc
iBase_SetRelation.o: mserver/iBase_SetRelation.cc
iBase_SetTag.o: mserver/iBase_SetTag.cc
iBase_Tag.o: mserver/iBase_Tag.cc
iMesh_Arr.o: mserver/iMesh_Arr.cc
iMesh_ArrMod.o: mserver/iMesh_ArrMod.cc
iMesh_Entity.o: mserver/iMesh_Entity.cc
iMesh_Mesh.o: mserver/iMesh_Mesh.cc
iMesh_Modify.o: mserver/iMesh_Modify.cc
iMesh_SIDL_MeshSidl.o: mserver/iMesh_SIDL_MeshSidl.cc
sidl_BaseClass.o: mserver/sidl_BaseClass.cc
sidl_BaseException.o: mserver/sidl_BaseException.cc
sidl_BaseInterface.o: mserver/sidl_BaseInterface.cc
sidl_ClassInfo.o: mserver/sidl_ClassInfo.cc
sidl_ClassInfoI.o: mserver/sidl_ClassInfoI.cc
sidl_DFinder.o: mserver/sidl_DFinder.cc
sidl_DLL.o: mserver/sidl_DLL.cc
sidl_Finder.o: mserver/sidl_Finder.cc
sidl_InvViolation.o: mserver/sidl_InvViolation.cc
sidl_Loader.o: mserver/sidl_Loader.cc
sidl_PostViolation.o: mserver/sidl_PostViolation.cc
sidl_PreViolation.o: mserver/sidl_PreViolation.cc
sidl_SIDLException.o: mserver/sidl_SIDLException.cc
sidl_io_Deserializer.o: mserver/sidl_io_Deserializer.cc
sidl_io_IOException.o: mserver/sidl_io_IOException.cc
sidl_io_Serializeable.o: mserver/sidl_io_Serializeable.cc
sidl_io_Serializer.o: mserver/sidl_io_Serializer.cc
sidl_rmi_ConnectRegistry.o: mserver/sidl_rmi_ConnectRegistry.cc
sidl_rmi_InstanceHandle.o: mserver/sidl_rmi_InstanceHandle.cc
sidl_rmi_InstanceRegistry.o: mserver/sidl_rmi_InstanceRegistry.cc
sidl_rmi_Invocation.o: mserver/sidl_rmi_Invocation.cc
sidl_rmi_NetworkException.o: mserver/sidl_rmi_NetworkException.cc
sidl_rmi_ProtocolFactory.o: mserver/sidl_rmi_ProtocolFactory.cc
sidl_rmi_Response.o: mserver/sidl_rmi_Response.cc
iMesh_SIDL_MeshSidl_Skel.o: mserver/iMesh_SIDL_MeshSidl_Skel.cc
iMesh_SIDL_MeshSidl_Impl.o: mserver/iMesh_SIDL_MeshSidl_Impl.cc
iBase_Error_Impl.o: mserver/iBase_Error_Impl.cc
iBase_Error_IOR.o: mserver/iBase_Error_IOR.c
iBase_Error_Skel.o: mserver/iBase_Error_Skel.cc
iBase_Error.o: mserver/iBase_Error.cc
