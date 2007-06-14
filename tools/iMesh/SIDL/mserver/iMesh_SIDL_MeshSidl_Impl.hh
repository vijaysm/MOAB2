// 
// File:          iMesh_SIDL_MeshSidl_Impl.hh
// Symbol:        iMesh_SIDL.MeshSidl-v0.2
// Symbol Type:   class
// Babel Version: 0.10.12
// sidl Created:  20070614 17:39:48 CDT
// Generated:     20070614 17:39:54 CDT
// Description:   Server-side implementation for iMesh_SIDL.MeshSidl
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.12
// source-line   = 16
// source-url    = file:/home/tautges/MOAB/tools/iMesh/SIDL/iMesh_SIDL.sidl
// 

#ifndef included_iMesh_SIDL_MeshSidl_Impl_hh
#define included_iMesh_SIDL_MeshSidl_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_iMesh_SIDL_MeshSidl_IOR_h
#include "iMesh_SIDL_MeshSidl_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_iBase_Error_hh
#include "iBase_Error.hh"
#endif
#ifndef included_iBase_TagValueType_hh
#include "iBase_TagValueType.hh"
#endif
#ifndef included_iMesh_AdjacencyInfo_hh
#include "iMesh_AdjacencyInfo.hh"
#endif
#ifndef included_iMesh_CreationStatus_hh
#include "iMesh_CreationStatus.hh"
#endif
#ifndef included_iMesh_EntityTopology_hh
#include "iMesh_EntityTopology.hh"
#endif
#ifndef included_iMesh_EntityType_hh
#include "iMesh_EntityType.hh"
#endif
#ifndef included_iMesh_Mesh_hh
#include "iMesh_Mesh.hh"
#endif
#ifndef included_iMesh_StorageOrder_hh
#include "iMesh_StorageOrder.hh"
#endif
#ifndef included_iMesh_SIDL_MeshSidl_hh
#include "iMesh_SIDL_MeshSidl.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._includes)
// Insert-Code-Here {iMesh_SIDL.MeshSidl._includes} (includes or arbitrary code)
#include "iMesh.h"
// DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._includes)

namespace iMesh_SIDL { 

  /**
   * Symbol "iMesh_SIDL.MeshSidl" (version 0.2)
   */
  class MeshSidl_impl
  // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._inherits)
  // Insert-Code-Here {iMesh_SIDL.MeshSidl._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    MeshSidl self;

    // DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._implementation)
    // Insert-Code-Here {iMesh_SIDL.MeshSidl._implementation} (additional details)
    static iMesh_Instance imeshInstance;
    int imeshError;
    void processError() throw(::iBase::Error);
    
      // DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._implementation)

  private:
    // private default constructor (required)
    MeshSidl_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    MeshSidl_impl( struct iMesh_SIDL_MeshSidl__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~MeshSidl_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:
    /**
     * user defined static method.
     */
    static ::iMesh::Mesh
    newMesh (
      /* in */ const ::std::string& option
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * user defined non-static method.
     */
    void
    createTag (
      /* in */ const ::std::string& tag_name,
      /* in */ int32_t number_of_values,
      /* in */ ::iBase::TagValueType tag_type,
      /* out */ void*& tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    destroyTag (
      /* in */ void* tag_handle,
      /* in */ bool forced
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    ::std::string
    getTagName (
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getTagSizeValues (
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getTagSizeBytes (
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getTagHandle (
      /* in */ const ::std::string& tag_name
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    ::iBase::TagValueType
    getTagType (
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& tag_value,
      /* out */ int32_t& tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getIntData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    double
    getDblData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getEHData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<char> tag_value,
      /* in */ int32_t tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setIntData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ int32_t tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setDblData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ double tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEHData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ void* tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllTags (
      /* in */ void* entity_handle,
      /* inout */ ::sidl::array<void*>& tag_handles,
      /* out */ int32_t& tag_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvTag (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getIntArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<int32_t>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getDblArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<double>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEHArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<void*>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<char> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setIntArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<int32_t> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setDblArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<double> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEHArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<void*> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvArrTag (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& tag_value,
      /* in */ int32_t tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetIntData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* in */ int32_t tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetDblData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* in */ double tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetEHData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* in */ void* tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSetData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& tag_value,
      /* out */ int32_t& tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getEntSetIntData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    double
    getEntSetDblData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getEntSetEHData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllEntSetTags (
      /* in */ void* entity_set,
      /* inout */ ::sidl::array<void*>& tag_handles,
      /* out */ int32_t& tag_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntSetTag (
      /* in */ void* entity_set,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createEntSet (
      /* in */ bool isList,
      /* out */ void*& entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    destroyEntSet (
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isList (
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumEntSets (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSets (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops,
      /* inout */ ::sidl::array<void*>& contained_entset_handles,
      /* out */ int32_t& contained_entset_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntToSet (
      /* in */ void* entity_handle,
      /* inout */ void*& entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntFromSet (
      /* in */ void* entity_handle,
      /* inout */ void*& entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntArrToSet (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ void*& entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntArrFromSet (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ void*& entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntSet (
      /* in */ void* entity_set_to_add,
      /* inout */ void*& entity_set_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntSet (
      /* in */ void* entity_set_to_remove,
      /* inout */ void*& entity_set_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isEntContained (
      /* in */ void* containing_entity_set,
      /* in */ void* entity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isEntSetContained (
      /* in */ void* containing_entity_set,
      /* in */ void* contained_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addPrntChld (
      /* inout */ void*& parent_entity_set,
      /* inout */ void*& child_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvPrntChld (
      /* inout */ void*& parent_entity_set,
      /* inout */ void*& child_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isChildOf (
      /* in */ void* parent_entity_set,
      /* in */ void* child_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumChld (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumPrnt (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getChldn (
      /* in */ void* from_entity_set,
      /* in */ int32_t num_hops,
      /* inout */ ::sidl::array<void*>& child_handles,
      /* out */ int32_t& child_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getPrnts (
      /* in */ void* from_entity_set,
      /* in */ int32_t num_hops,
      /* inout */ ::sidl::array<void*>& parent_handles,
      /* out */ int32_t& parent_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    subtract (
      /* in */ void* entity_set_1,
      /* in */ void* entity_set_2,
      /* out */ void*& result_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    intersect (
      /* in */ void* entity_set_1,
      /* in */ void* entity_set_2,
      /* out */ void*& result_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    unite (
      /* in */ void* entity_set_1,
      /* in */ void* entity_set_2,
      /* out */ void*& result_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    load (
      /* in */ void* entity_set_handle,
      /* in */ const ::std::string& name,
      /* in */ const ::std::string& options
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    save (
      /* in */ void* entity_set_handle,
      /* in */ const ::std::string& name,
      /* in */ const ::std::string& options
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getRootSet() throw ( 
      ::iBase::Error
    );
    /**
     * user defined non-static method.
     */
    int32_t
    getGeometricDim() throw ( 
      ::iBase::Error
    );
    /**
     * user defined non-static method.
     */
    ::iMesh::StorageOrder
    getDfltStorage() throw ( 
      ::iBase::Error
    );
    /**
     * user defined non-static method.
     */
    void
    getAdjTable (
      /* inout */ ::sidl::array< ::iMesh::AdjacencyInfo>& adjacency_table,
      /* out */ int32_t& adjacency_table_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    areEHValid (
      /* in */ bool reset
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumOfType (
      /* in */ void* entity_set_handle,
      /* in */ ::iMesh::EntityType entity_type
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumOfTopo (
      /* in */ void* entity_set_handle,
      /* in */ ::iMesh::EntityTopology entity_topology
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllVtxCoords (
      /* in */ void* entity_set,
      /* inout */ ::sidl::array<double>& coords,
      /* out */ int32_t& coords_size,
      /* inout */ ::sidl::array<int32_t>& in_entity_set,
      /* out */ int32_t& in_entity_set_size,
      /* inout */ ::iMesh::StorageOrder& storage_order
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxCoordIndex (
      /* in */ void* entity_set,
      /* in */ ::iMesh::EntityType requested_entity_type,
      /* in */ ::iMesh::EntityTopology requested_entity_topology,
      /* in */ ::iMesh::EntityType entity_adjacency_type,
      /* inout */ ::sidl::array<int32_t>& offset,
      /* out */ int32_t& offset_size,
      /* inout */ ::sidl::array<int32_t>& index,
      /* out */ int32_t& index_size,
      /* inout */ ::sidl::array< ::iMesh::EntityTopology>& entity_topologies,
      /* out */ int32_t& entity_topologies_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntities (
      /* in */ void* entity_set,
      /* in */ ::iMesh::EntityType entity_type,
      /* in */ ::iMesh::EntityTopology entity_topology,
      /* inout */ ::sidl::array<void*>& entity_handles,
      /* out */ int32_t& entity_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxArrCoords (
      /* in */ ::sidl::array<void*> vertex_handles,
      /* in */ int32_t vertex_handles_size,
      /* inout */ ::iMesh::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& coords,
      /* out */ int32_t& coords_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAdjEntities (
      /* in */ void* entity_set,
      /* in */ ::iMesh::EntityType entity_type_requestor,
      /* in */ ::iMesh::EntityTopology entity_topology_requestor,
      /* in */ ::iMesh::EntityType entity_type_requested,
      /* inout */ ::sidl::array<void*>& adj_entity_handles,
      /* out */ int32_t& adj_entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& offset,
      /* out */ int32_t& offset_size,
      /* inout */ ::sidl::array<int32_t>& in_entity_set,
      /* out */ int32_t& in_entity_set_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    initEntIter (
      /* in */ void* entity_set_handle,
      /* in */ ::iMesh::EntityType requested_entity_type,
      /* in */ ::iMesh::EntityTopology requested_entity_topology,
      /* out */ void*& entity_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    getNextEntIter (
      /* in */ void* entity_iterator,
      /* out */ void*& entity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    resetEntIter (
      /* in */ void* entity_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    endEntIter (
      /* in */ void* entity_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    ::iMesh::EntityTopology
    getEntTopo (
      /* in */ void* entity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    ::iMesh::EntityType
    getEntType (
      /* in */ void* entity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxCoord (
      /* in */ void* vertex_handle,
      /* out */ double& x,
      /* out */ double& y,
      /* out */ double& z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntAdj (
      /* in */ void* entity_handle,
      /* in */ ::iMesh::EntityType entity_type_requested,
      /* inout */ ::sidl::array<void*>& adj_entity_handles,
      /* out */ int32_t& adj_entity_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    initEntArrIter (
      /* in */ void* entity_set_handle,
      /* in */ ::iMesh::EntityType requested_entity_type,
      /* in */ ::iMesh::EntityTopology requested_entity_topology,
      /* in */ int32_t requested_array_size,
      /* out */ void*& entArr_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    getNextEntArrIter (
      /* in */ void* entArr_iterator,
      /* inout */ ::sidl::array<void*>& entity_handles,
      /* out */ int32_t& entity_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    resetEntArrIter (
      /* in */ void* entArr_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    endEntArrIter (
      /* in */ void* entArr_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntArrTopo (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array< ::iMesh::EntityTopology>& topology,
      /* out */ int32_t& topology_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntArrType (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array< ::iMesh::EntityType>& type,
      /* out */ int32_t& type_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntArrAdj (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iMesh::EntityType entity_type_requested,
      /* inout */ ::sidl::array<void*>& adj_entity_handles,
      /* out */ int32_t& adj_entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& offset,
      /* out */ int32_t& offset_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setVtxCoords (
      /* in */ void* vertex_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createVtx (
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ void*& new_vertex_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createEnt (
      /* in */ ::iMesh::EntityTopology new_entity_topology,
      /* in */ ::sidl::array<void*> lower_order_entity_handles,
      /* in */ int32_t lower_order_entity_handles_size,
      /* out */ void*& new_entity_handle,
      /* out */ ::iMesh::CreationStatus& status
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    deleteEnt (
      /* in */ void* entity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setVtxArrCoords (
      /* in */ ::sidl::array<void*> vertex_handles,
      /* in */ int32_t vertex_handles_size,
      /* in */ ::iMesh::StorageOrder storage_order,
      /* in */ ::sidl::array<double> new_coords,
      /* in */ int32_t new_coords_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createVtxArr (
      /* in */ int32_t num_verts,
      /* in */ ::iMesh::StorageOrder storage_order,
      /* in */ ::sidl::array<double> new_coords,
      /* in */ int32_t new_coords_size,
      /* inout */ ::sidl::array<void*>& new_vertex_handles,
      /* out */ int32_t& new_vertex_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createEntArr (
      /* in */ ::iMesh::EntityTopology new_entity_topology,
      /* in */ ::sidl::array<void*> lower_order_entity_handles,
      /* in */ int32_t lower_order_entity_handles_size,
      /* inout */ ::sidl::array<void*>& new_entity_handles,
      /* out */ int32_t& new_entity_handles_size,
      /* inout */ ::sidl::array< ::iMesh::CreationStatus>& status,
      /* out */ int32_t& status_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    deleteEntArr (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size
    )
    throw ( 
      ::iBase::Error
    );

  };  // end class MeshSidl_impl

} // end namespace iMesh_SIDL

// DO-NOT-DELETE splicer.begin(iMesh_SIDL.MeshSidl._misc)
// Insert-Code-Here {iMesh_SIDL.MeshSidl._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(iMesh_SIDL.MeshSidl._misc)

#endif
