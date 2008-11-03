#ifndef IMESH_CBIND_H__
#define IMESH_CBIND_H__

  /** \mainpage The ITAPS Mesh Interface iMesh
   *
   * The ITAPS Mesh Interface iMesh provides a common interface for
   * accessing mesh and data associated with a mesh.  Applications written
   * to use this interface can use a variety of implementations, choosing
   * the one that best meets its needs.  They can also use tools written
   * to this interface, for example mesh smoothing, adaptive mesh refinement,
   * and parallel mesh support.
   *
   * \section ITAPS Data Model
   *
   * The ITAPS interfaces use a data model composed of four basic data types: \n
   * \em Entity: basic topological entities in a mesh, e.g. vertices, 
   * triangles, hexahedra. \n
   * \em Entity \em Set: arbitrary grouping of other entities and sets. 
   * Entity sets also support parent/child relations with other sets which
   * are distinct from entities contained in those sets.  Parent/child links
   * can be used to embed graph relationships between sets, e.g. to 
   * represent topological relationships between the sets. \n
   * \em Interface: the object with which mesh is associated and on which
   * functions in iMesh are called. \n
   * \em Tag: application data associated with objects of any of the other 
   * data types.  Each tag has a designated name, size, and data type.
   *
   * \section JTAPS Entity Type, Topology
   * Each entity has a specific Entity Type and Entity Topology.  The Entity 
   * Type is one of VERTEX, EDGE, FACE, and REGION, and is synonymous with
   * the topological dimension of the entity.  The Entity Topology denotes
   * the specific shape, for example TRIANGLE, QUADRILATERAL, TETRAHEDRON,
   * and HEXAHEDRON.  Entity Type and Entity Topology exist as enumerated
   * types, Entity Type in the iBase_EntityType enumeration, and
   * Entity Topology in the iMesh_EntityTopology enumeration.
   *
   * \section KTAPS Entity-, Array-, and Iterator-Based Access
   *
   * The iMesh interface provides functions for accessing entities
   * individually, as arrays of entities, or using iterators.  These access
   * methods have different memory versus execution time tradeoffs, 
   * depending on the implementation.
   *
   * \section LTAPS Lists Passed Through Interface
   *
   * Many of the functions in iMesh have arguments corresponding to lists of 
   * objects.  In-type arguments for lists consist of a pointer to an array and
   * a list size.  Lists returned from functions are passed in three arguments,
   * a pointer to the array representing the list, and pointers to the
   * allocated and occupied lengths of the array.  These three arguments are 
   * inout-type arguments, because they can be allocated by the application and
   * passed into the interface to hold the results of the function.  Lists
   * which are pre-allocated must be large enough to hold the results of the
   * function; if this is not the case, an error is generated.  Otherwise, the
   * occupied size is changed to the size output from the function.  If a list
   * argument is unallocated (the list pointer points to a NULL value) or if
   * the incoming value of the allocated size is zero, the list storage will be
   * allocated by the implementation.  IN ALL CASES, MEMORY ALLOCATED BY ITAPS
   * INTERFACE IMPLEMENTATIONS IS DONE USING THE C MALLOC FUNCTION, AND CAN BE
   * DE-ALLOCATED USING THE C FREE FUNCTION.
   *
   */

#ifndef ITAPS
#define ITAPS
#endif

#include "iBase.h"
#include "iMesh_protos.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**\brief  Type used to store iMesh interface handle
     *
     * Type used to store iMesh interface handle
     */
  typedef void* iMesh_Instance;

    /**\brief  Type used to store an iterator returned by iMesh
     *
     * Type used to store an iterator returned by iMesh
     */
  typedef void* iMesh_EntityIterator;

    /**\brief  Type used to store an array iterator returned by iMesh
     *
     * Type used to store an array iterator returned by iMesh
     */
  typedef void* iMesh_EntityArrIterator;

    /**\brief  Enumerator specifying entity topology
     *
     * Enumerator specifying entity topology.
     */
  enum iMesh_EntityTopology {
    iMesh_POINT = 0,              /**< a general zero-dimensional entity  */
    iMesh_LINE_SEGMENT,       /**< a general one-dimensional entity  */
    iMesh_POLYGON,            /**< a general two-dimensional element  */
    iMesh_TRIANGLE,           /**< a three-sided, two-dimensional element  */
    iMesh_QUADRILATERAL,      /**< a four-sided, two-dimensional element  */
    iMesh_POLYHEDRON,         /**< a general three-dimensional element */
    iMesh_TETRAHEDRON,        /**< a four-sided, three-dimensional element whose
			       *     faces are triangles */
    iMesh_HEXAHEDRON,         /**< a six-sided, three-dimensional element whose
			       *     faces are quadrilaterals */
    iMesh_PRISM,              /**< a five-sided, three-dimensional element which
			       *     has three quadrilateral faces and two
			       *     triangular faces  */
    iMesh_PYRAMID,            /**< a five-sided, three-dimensional element
			       *     which has one quadrilateral face and four
			       *     triangular faces */
    iMesh_SEPTAHEDRON,        /**< a hexahedral entity with one collapsed edge */
    iMesh_ALL_TOPOLOGIES      /**< allows the user to request information
			       *     about all the topology types */
  };

    /**\brief  Get the error type returned from the last iMesh function
     *
     * Get the error type returned from the last iMesh function.  Value
     * returned is a member of the iBase_ErrorType enumeration.
     * \param instance iMesh instance handle
     * \param *error_type Error type returned from last iMesh function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getErrorType(iMesh_Instance instance,
                          /*out*/ int *error_type, 
                          int *err);

    /**\brief  Get a description of the error returned from the last iMesh function
     *
     * Get a description of the error returned from the last iMesh function
     * \param instance iMesh instance handle
     * \param descr Pointer to a character string to be filled with a
     *        description of the error from the last iMesh function
     * \param *err Pointer to error type returned from function
     * \param descr_len Length of the character string pointed to by descr
     */
  void iMesh_getDescription(iMesh_Instance instance,
                            /*inout*/ char *descr, 
                            int *err, 
                            /*in*/ int descr_len);

    /**\brief  Construct a new iMesh instance
     *
     * Construct a new iMesh instance, using implementation-specific
     * options
     * \param options Pointer to implementation-specific options string
     * \param instance Pointer to iMesh instance handle returned from function
     * \param *err Pointer to error type returned from function
     * \param options_len Length of the character string pointed to by options
     */
  void iMesh_newMesh(const char *options,
                     /*out*/ iMesh_Instance *instance, 
                     /*out*/ int *err, 
                     /*in*/ int options_len);

    /**\brief  Destroy an iMesh instance
     *
     * Destroy an iMesh instance
     * \param instance iMesh instance to be destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_dtor(iMesh_Instance instance, 
                  /*out*/ int *err);

    /**\brief  Load a mesh from a file
     *
     * Load a mesh from a file.  If entity set is specified, loaded mesh
     * is added to that set; specify zero if that is not desired.
     * \param instance iMesh instance handle
     * \param entity_set_handle Set to which loaded mesh will be added, zero
     *        if not desired
     * \param name File name from which mesh is to be loaded
     * \param options Pointer to implementation-specific options string
     * \param *err Pointer to error type returned from function
     * \param name_len Length of the file name character string
     * \param options_len Length of the options character string
     */
  void iMesh_load(iMesh_Instance instance,
                  /*in*/ const iBase_EntitySetHandle entity_set_handle,
                  /*in*/ const char *name, 
                  /*in*/ const char *options,
                  /*out*/ int *err, 
                  /*in*/ int name_len, 
                  /*in*/ int options_len);

    /**\brief  Save a mesh to a file
     *
     * Save a mesh to a file.  If entity set is specified, save only the
     * mesh contained in that set.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being saved
     * \param name File name to which mesh is to be saved
     * \param options Pointer to implementation-specific options string
     * \param *err Pointer to error type returned from function
     * \param name_len Length of the file name character string
     * \param options_len Length of the options character string
     */
  void iMesh_save(iMesh_Instance instance,
                  /*in*/ const iBase_EntitySetHandle entity_set_handle,
                  /*in*/ const char *name, 
                  /*in*/ const char *options,
                  /*out*/ int *err, 
                  /*in*/ const int name_len, 
                  /*in*/ int options_len);

    /**\brief  Get handle of the root set for this instance
     *
     * Get handle of the root set for this instance.  All mesh in
     * this instance can be accessed from this set.
     * \param instance iMesh instance handle
     * \param root_set Pointer to set handle returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getRootSet(iMesh_Instance instance,
                        /*out*/ iBase_EntitySetHandle *root_set, 
                        /*out*/ int *err);

    /**\brief  Get the geometric dimension of mesh represented in this instance
     *
     * Get the geometric dimension of mesh represented in this instance
     * \param instance iMesh instance handle
     * \param geom_dim Pointer to dimension returned from this function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getGeometricDimension(iMesh_Instance instance,
                                   /*out*/ int *geom_dim, 
                                   /*out*/ int *err);

    /**\brief  Get the default storage order used by this implementation
     *
     * Get the default storage order used by this implementation.  Value
     * returned is a member of the iBase_StorageOrder enumeration.
     * \param instance iMesh instance handle
     * \param order Pointer to storage order returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getDfltStorage(iMesh_Instance instance,
                            /*out*/ int *order, 
                            /*out*/ int *err);

    /**\brief  Get the adjacency table for this implementation
     *
     * Get the adjacency table for this implementation.  This table 
     * is a 4x4 array, with indices 0-based, where A(i,j) (i=row, j=column) 
     * represents the relative cost of retrieving adjacencies between
     * entities of dimension i to entities of dimension j
     * \param instance iMesh instance handle
     * \param *adjacency_table Pointer to array representing adjacency table
     *        returned from function
     * \param adjacency_table_allocated Pointer to allocated size of 
     *        adjacency table
     * \param adjacency_table_size Pointer to occupied size of 
     *        adjacency table
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAdjTable (iMesh_Instance instance,
                          /*out*/ int** adjacency_table,
                          /*inout*/ int* adjacency_table_allocated,
                          /*out*/ int* adjacency_table_size, 
                          /*out*/ int *err);

    /**\brief  Set the adjacency table and interior entity information for this implementation
     *
     * Set the adjacency table and interior entity information for this 
     * implementation.  This table 
     * is a 4x4 array, with indices 0-based, where A(i,j) (i=row, j=column) 
     * non-zero requests that adjacencies be stored explicitly from type i
     * to type j.  Non-zero diagonal elements request that interior entities
     * of that dimension be represented explicitly, and created along with
     * higher-dimensional entities.
     * \param instance iMesh instance handle
     * \param *adjacency_table Array representing adjacency table
     *        passed to function
     * \param adjacency_table_size Size of adjacency table (should be 16)
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setAdjTable (iMesh_Instance instance,
                          /*in*/ int* adjacency_table,
                          /*in*/ int adjacency_table_size,
                          /*out*/ int *err);

    /**\brief  Get the number of entities with the specified type in the instance or set
     *
     * Get the number of entities with the specified type in the instance 
     * or set.  If entity set handle is zero, return information for instance,
     * otherwise for set.  Value of entity type must be from the
     * iBase_EntityType enumeration.  If iBase_ALL_TYPES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entity requested
     * \param num_type Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfType(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_type,
                          /*out*/ int *num_type, 
                          /*out*/ int *err);

    /**\brief  Get the number of entities with the specified topology in the instance or set
     *
     * Get the number of entities with the specified topology in the instance 
     * or set.  If entity set handle is zero, return information for instance,
     * otherwise for set.  Value of entity topology must be from the
     * iMesh_EntityTopology enumeration.  If iMesh_ALL_TOPOLOGIES is specified,
     * total number of entities (excluding entity sets) is returned.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_topology Topology of entity requested
     * \param num_topo Pointer to number of entities, returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumOfTopo(iMesh_Instance instance,
                          /*in*/ const iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const int entity_topology,
                          /*out*/ int *num_topo, 
                          /*out*/ int *err);

    /**\brief Return whether entity handles have changed since last reset or since
     *        instance construction
     *
     * Return whether entity handles have changed since last reset or since
     * instance construction.  If non-zero value is returned, it is not
     * guaranteed that a handle from before the last call to this function
     * represents the same entity as the same handle value does now.  If
     * doReset is non-zero, resets the starting point for this function.
     * \param instance iMesh instance handle
     * \param doReset Perform a reset on the starting point after which handles
     *        are invariant.
     * \param areHandlesInvariant Pointer to invariant flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_areEHValid(iMesh_Instance instance, 
                        /*in*/ int doReset,
                        /*out*/ int *areHandlesInvariant, 
                        /*out*/ int *err);

    /**\brief  Return coordinates of all vertices defined in this instance or set
     *
     * Return coordinates of all vertices defined in this instance or set.
     * If non-zero set is input, return coordinates for vertices either in
     * the set or contained by entities in the set.  If storage order is
     * a value other than iBase_UNDETERMINED, coordinates are returned with
     * that storage order; otherwise storage order is in native order with
     * respect to the implementation.  On return, storage order contains
     * order of the returned coordinates.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set for which vertex coordinates are
     *        requested
     * \param *coordinates Pointer to array holding coordinates
     * \param *coordinates_allocated Pointer to allocated size of coordinates
     *        array
     * \param *coordinates_size Pointer to occupied size of coordinates
     *        array
     * \param *in_entity_set Pointer to array of flags; if entity set was
     *        input, (*in_entity_set)[i]=1 indicates that vertex i was in
     *        the entity set, zero otherwise.
     * \param *in_entity_set_allocated Pointer to allocated size of 
     *        in_entity_set array
     * \param *in_entity_set_size Pointer to occupied size of 
     *        in_entity_set array
     * \param *storage_order Pointer to storage order requested/returned from
     *        function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAllVtxCoords (iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*inout*/ double** coordinates,
                              /*inout*/ int* coordinates_allocated,
                              /*out*/ int* coordinates_size,
                              /*inout*/ int** in_entity_set,
                              /*inout*/ int* in_entity_set_allocated,
                              /*out*/ int* in_entity_set_size,
                              /*inout*/ int* storage_order, /*out*/ int *err);

    /**\brief  Get adjacent entities as connectivity lists
     *
     * For a specified set handle, entity type and/or entity topology,
     * return connectivity of adjacent entities of specified dimension.  
     * Connectivity is expressed as index into array of vertices returned
     * by getEntities function.  \em offset[i] is the index of the first
     * vertex of the first entity adjacent to entity \em i in the entity set
     * or the mesh.  The topology of entities whose connectivity is returned
     * is given in the entity_topologies array.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param requested_entity_type Type of entity requested
     * \param requested_entity_topology Topology of entity requested
     * \param entity_adjacency_type Adjacency type
     * \param *offset Pointer to array of offsets returned from function
     * \param *offset_allocated Pointer to allocated size of offset array
     * \param *offset_size Pointer to occupied size of offset array
     * \param *index Pointer to array of indices returned from function
     * \param *index_allocated Pointer to allocated size of index array
     * \param *index_size Pointer to occupied size of index array
     * \param *entity_topologies Pointer to array of entity topologies 
     *        returned from function
     * \param *entity_topologies_allocated Pointer to allocated size of 
     *        entity_topologies array
     * \param *entity_topologies_size Pointer to occupied size of 
     *        entity_topologies array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getVtxCoordIndex (iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle entity_set_handle,
                               /*in*/ const int requested_entity_type,
                               /*in*/ const int requested_entity_topology,
                               /*in*/ const int entity_adjacency_type,
                               /*inout*/ int** offset,
                               /*inout*/ int* offset_allocated,
                               /*out*/ int* offset_size,
                               /*inout*/ int** index,
                               /*inout*/ int* index_allocated,
                               /*out*/ int* index_size,
                               /*inout*/  int** entity_topologies,
                               /*inout*/ int* entity_topologies_allocated,
                               /*out*/ int* entity_topologies_size, 
                               /*out*/ int *err);

    /**\brief  Get entities of specific type and/or topology in set or instance
     *
     * Get entities of specific type and/or topology in set or instance.  All 
     * entities of a given type or topology are requested by specifying
     * iBase_ALL_TOPOLOGIES or iBase_ALL_TYPES, respectively.  Specified type
     * or topology must be a value in the iBase_EntityType or iBase_EntityTopology
     * enumeration, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param entity_type Type of entities being requested
     * \param entity_topology Topology of entities being requested
     * \param *entity_handles Pointer to array of entity handles returned 
     *        from function
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntities(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int entity_type,
                         /*in*/ const int entity_topology,
                         /*out*/ iBase_EntityHandle** entity_handles,
                         /*out*/ int* entity_handles_allocated,
                         /*out*/ int* entity_handles_size,
                         /*out*/ int *err);

    /**\brief  Get coordinates of specified vertices
     *
     * Get coordinates of specified vertices.  If storage order is passed in
     * with a value other than iBase_UNDETERMINED, coordinates are returned
     * in the specified storage order, otherwise storage order is that native
     * to the implementation.  Storage order of returned coordinates is also
     * returned.
     * \param instance iMesh instance handle
     * \param vertex_handles Array of mesh vertex handles whose coordinates are
     *        being requested
     * \param vertex_handles_size Number of vertices in vertex_handles array
     * \param storage_order Pointer to storage order requested/returned from
     *        function
     * \param *coords Pointer to array of coordinates returned from function
     * \param *coords_allocated Pointer to allocated size of coords array
     * \param *coords_size Pointer to occupied size of coords array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getVtxArrCoords(iMesh_Instance instance,
                             /*in*/ const iBase_EntityHandle* vertex_handles,
                             /*in*/ const int vertex_handles_size,
                             /*inout*/ int* storage_order,
                             /*inout*/ double** coords,
                             /*inout*/ int* coords_allocated,
                             /*out*/ int* coords_size, 
                             /*out*/ int *err);


    /**\brief  Get entities adjacent to those of requested type/topology in set or instance
     *
     * Get entities adjacent to those of requested type/topology in set or 
     * instance.  If set is specified, also returns whether each adjacent entity 
     * is in the set ((*in_entity_set)[i]=1 indicates that adjacent entity
     * (*adj_entity_handles)[i] is in the specified set).  If adjacencies for 
     * all entities of specified type or
     * topology are desired, specify iMesh_ALL_TOPOLOGIES or iBase_ALL_TYPES,
     * respectively.  Types and topologies input to this function must be values
     * in the iBase_EntityType or iBase_EntityTopology enumerations, respectively.
     * \em offset[i] stores index of first adjacent entity to entity \em i from
     * mesh or set.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set whose adjacent entities are requested
     *        (zero indicates whole mesh)
     * \param entity_type_requestor Return entities adjacent to entities of this 
     *        type
     * \param entity_topology_requestor Return entities adjacent to entities of 
     *        this topology
     * \param entity_type_requested Return adjacent entities of this type
     * \param *adj_entity_handles Pointer to array of adjacent entity handles 
     *        returned from function
     * \param *adj_entity_handles_allocated Pointer to allocated size of 
     *        adj_entity_handles array
     * \param *adj_entity_handles_size Pointer to occupied size of 
     *        adj_entity_handles array
     * \param *offset Pointer to array of offsets returned from function
     * \param *offset_allocated Pointer to allocated size of offset array
     * \param *offset_size Pointer to occupied size of offset array
     * \param *in_entity_set Pointer to array of in_entity_set returned 
     *        from function
     * \param *in_entity_set_allocated Pointer to allocated size of 
     *        in_entity_set array
     * \param *in_entity_set_size Pointer to occupied size of in_entity_set 
     *        array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAdjEntities(iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle entity_set_handle,
                            /*in*/ const int entity_type_requestor,
                            /*in*/ const int entity_topology_requestor,
                            /*in*/ const int entity_type_requested,
                            /*inout*/ iBase_EntityHandle** adj_entity_handles,
                            /*inout*/ int* adj_entity_handles_allocated,
                            /*out*/ int* adj_entity_handles_size,
                            /*inout*/ int** offset,
                            /*inout*/ int* offset_allocated,
                            /*out*/ int* offset_size,
                            /*inout*/ int** in_entity_set,
                            /*inout*/ int* in_entity_set_allocated,
                            /*out*/ int* in_entity_set_size, 
                            /*out*/ int *err);

    /**\brief Initialize an array iterator over specified entity type, topology, and 
     *        size
     *
     * Initialize an array iterator over specified entity type, topology, and 
     * size, for a specified set or instance.  Iterator returned can be used 
     * as input to functions returning entities for the iterator.  If all 
     * entities of a specified type and/or topology are to be iterated, 
     * specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  
     * Specified type or topology must be a value in the iBase_EntityType or 
     * iMesh_EntityTopology enumerations, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being iterated
     * \param requested_entity_type Type of entity to iterate
     * \param requested_entity_topology Topology of entity to iterate
     * \param requested_array_size Size of chunks of handles returned for each
     *        value of the iterator
     * \param entArr_iterator Pointer to iterator returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_initEntArrIter(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle entity_set_handle,
                            /*in*/ const int requested_entity_type,
                            /*in*/ const int requested_entity_topology,
                            /*in*/ const int requested_array_size,
                            /*out*/ iMesh_EntityArrIterator* entArr_iterator,
                            /*out*/ int *err);

    /**\brief  Get entities contained in array iterator and increment iterator
     *
     * Get the entities contained in an array iterator, and increment the 
     * iterator.  Also return whether the next value of the iterator has
     * any entities (if non-zero, next iterator value is the end of the
     * iteration).
     * \param instance iMesh instance handle
     * \param entArr_iterator Iterator being queried
     * \param *entity_handles Pointer to array of entity handles contained in
     *        current value of iterator
     * \param *entity_handles_allocated Pointer to allocated size of 
     *        entity_handles array
     * \param *entity_handles_size Pointer to occupied size of entity_handles 
     *        array
     * \param has_data Pointer to flag; if returned non-zero, next iterator
     *        has a non-zero number of entities
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNextEntArrIter(iMesh_Instance instance,
                               /*in*/ iMesh_EntityArrIterator entArr_iterator,
                               /*inout*/ iBase_EntityHandle** entity_handles,
                               /*inout*/ int* entity_handles_allocated,
                               /*out*/ int* entity_handles_size,
                               /*out*/ int *has_data, 
                               /*out*/ int *err);


    /**\brief  Reset the array iterator
     *
     * Reset the array iterator
     * \param instance iMesh instance handle
     * \param entArr_iterator Iterator to reset
     * \param *err Pointer to error type returned from function
     */
  void iMesh_resetEntArrIter(iMesh_Instance instance,
                             /*in*/ iMesh_EntityArrIterator entArr_iterator, 
                             /*out*/ int *err);


    /**\brief  Destroy the specified array iterator
     *
     * Destroy the specified array iterator
     * \param instance iMesh instance handle
     * \param entArr_iterator Iterator which gets destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_endEntArrIter(iMesh_Instance instance,
                           /*in*/ iMesh_EntityArrIterator entArr_iterator, 
                           /*out*/ int *err);

    /**\brief  Get the entity topology for the specified entities
     *
     * Get the entity topology for the specified entities.  Topologies 
     * returned are values in the iMesh_EntityTopology enumeration.
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles being queried
     * \param entity_handles_size Number of entities in entity_handles array
     * \param *topology Pointer to array of entity topologies returned 
     *        from function
     * \param *topology_allocated Pointer to allocated size of topology array
     * \param *topology_size Pointer to occupied size of topology array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntArrTopo(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** topology,
                           /*inout*/ int* topology_allocated,
                           /*out*/ int* topology_size, 
                           /*out*/ int *err);


    /**\brief  Get the entity type for the specified entities
     *
     * Get the entity type for the specified entities.  Types
     * returned are values in the iBase_EntityType enumeration.
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles being queried
     * \param entity_handles_size Number of entities in entity_handles array
     * \param *type Pointer to array of types returned from function
     * \param *type_allocated Pointer to allocated size of type array
     * \param *type_size Pointer to occupied size of type array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntArrType(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*inout*/ int** type,
                           /*inout*/ int* type_allocated,
                           /*out*/ int* type_size, 
                           /*out*/ int *err);


    /**\brief  Get entities of specified type adjacent to entities
     *
     * Get entities of specified type adjacent to entities.  Specified type
     * must be value in the iBase_EntityType enumeration.  \em offset(i) is
     * index of first entity in adjacentEntityHandles array adjacent to 
     * entity_handles[i].
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles being queried
     * \param entity_handles_size Number of entities in entity_handles array
     * \param entity_type_requested Type of adjacent entities requested
     * \param *adjacentEntityHandles Pointer to array of adjacentEntityHandles 
     *        returned from function
     * \param *adjacentEntityHandles_allocated Pointer to allocated size of 
     *        adjacentEntityHandles array
     * \param *adj_entity_handles_size Pointer to occupied size of 
     *        adjacentEntityHandles array
     * \param *offset Pointer to array of offsets returned from function
     * \param *offset_allocated Pointer to allocated size of offset array
     * \param *offset_size Pointer to occupied size of offset array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntArrAdj(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const int entity_type_requested,
                          /*inout*/ iBase_EntityHandle** adjacentEntityHandles,
                          /*inout*/ int* adjacentEntityHandles_allocated,
                          /*out*/ int* adj_entity_handles_size,
                          /*inout*/ int** offset,
                          /*inout*/ int* offset_allocated,
                          /*out*/ int* offset_size,
                          /*out*/ int *err);

/**\brief Get "2nd order" adjacencies to an array of entities
 * Get "2nd order" adjacencies to an array of entities, that is, from each 
 * entity, through other entities of a specified "bridge" dimension, to other 
 * entities of another specified "to" dimension.
 *
 * \param instance iMesh instance for this call
 * \param entity_handles Entities from which adjacencies are requested
 * \param entity_handles_size Number of entities whose adjacencies are requested
 * \param order_adjacent_key Bridge dimension for 2nd order adjacencies
 * \param requested_entity_type Dimension of adjacent entities returned
 * \param *adj_entity_handles Adjacent entities
 * \param *adj_entity_handles_allocated Allocated size of returned array
 * \param *adj_entity_handles_size Occupied size of returned array
 * \param *offset Offset[i] is offset into adj_entity_handles of 2nd order 
 *        adjacencies of ith entity in entity_handles
 * \param *offset_allocated Allocated size of offset array
 * \param *offset_size Occupied size of offset array
 * \param *err Pointer to error type returned from function
 */
  void iMesh_getEntArr2ndAdj( iMesh_Instance instance,
                              iBase_EntityHandle const* entity_handles,
                              int entity_handles_size,
                              int order_adjacent_key,
                              int requested_entity_type,
                              iBase_EntityHandle** adj_entity_handles,
                              int* adj_entity_handles_allocated,
                              int* adj_entity_handles_size,
                              int** offset,
                              int* offset_allocated,
                              int* offset_size,
                              int* err );

    /**\brief  Create an entity set
     *
     * Create an entity set, either ordered (isList=1) or unordered 
     * (isList=0).  Unordered entity sets can contain a given entity or 
     * set only once.
     * \param instance iMesh instance handle
     * \param isList If non-zero, an ordered list is created, otherwise an
     *        unordered set is created.
     * \param entity_set_created Entity set created by function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createEntSet(iMesh_Instance instance,
                          /*in*/ const int isList,
                          /*out*/ iBase_EntitySetHandle* entity_set_created,
                          /*out*/ int *err);

    /**\brief  Destroy an entity set
     *
     * Destroy an entity set
     * \param instance iMesh instance handle
     * \param entity_set Entity set to be destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_destroyEntSet(iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set,
                           /*out*/ int *err);

    /**\brief  Return whether a specified set is ordered or unordered
     *
     * Return whether a specified set is ordered (*is_list=1) or 
     * unordered (*is_list=0)
     * \param instance iMesh instance handle
     * \param entity_set Entity set being queried
     * \param is_list Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isList(iMesh_Instance instance,
                    /*in*/ const iBase_EntitySetHandle entity_set,
                    /*out*/ int *is_list,
                    /*out*/ int *err);

    /**\brief  Get the number of entity sets contained in a set or interface
     *
     * Get the number of entity sets contained in a set or interface.  If
     * a set is input which is not the root set, num_hops indicates the 
     * maximum number of contained sets from entity_set_handle to one of the
     * contained sets, inclusive of the contained set.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to contained set,
     *        inclusive of the contained set
     * \param num_sets Pointer to the number of sets returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumEntSets(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const int num_hops,
                           /*out*/ int *num_sets,
                           /*out*/ int *err);


    /**\brief  Get the entity sets contained in a set or interface
     *
     * Get the entity sets contained in a set or interface.  If
     * a set is input which is not the root set, num_hops indicates the 
     * maximum number of contained sets from entity_set_handle to one of the
     * contained sets, inclusive of the contained set.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to contained set,
     *        inclusive of the contained set
     * \param *contained_set_handles Pointer to array of set handles returned
     *        from function
     * \param contained_set_handles_allocated Pointer to allocated length of
     *        contained_set_handles array
     * \param contained_set_handles_size Pointer to occupied length of
     *        contained_set_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSets(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set_handle,
                        /*in*/ const int num_hops,
                        /*out*/ iBase_EntitySetHandle** contained_set_handles,
                        /*out*/ int* contained_set_handles_allocated,
                        /*out*/ int* contained_set_handles_size,
                        /*out*/ int *err);

    /**\brief  Add an entity to a set
     *
     * Add an entity to a set
     * \param instance iMesh instance handle
     * \param entity_handle The entity being added
     * \param entity_set Pointer to the set being added to
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addEntToSet(iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle entity_handle,
                         /*inout*/ iBase_EntitySetHandle* entity_set,
                         /*out*/ int *err);

    /**\brief  Remove an entity from a set
     *
     * Remove an entity from a set
     *
     * \param instance iMesh instance handle
     * \param entity_handle The entity being removed
     * \param entity_set Pointer to the set being removed from
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntFromSet(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle entity_handle,
                           /*inout*/ iBase_EntitySetHandle* entity_set,
                           /*out*/ int *err);


    /**\brief  Add an array of entities to a set
     *
     * Add an array of entities to a set
     * \param instance iMesh instance handle
     * \param entity_handles Array of entities being added
     * \param entity_handles_size Number of entities in entity_handles array
     * \param entity_set Pointer to the set being added to
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addEntArrToSet(iMesh_Instance instance,
                            /*in*/ const iBase_EntityHandle* entity_handles,
                            /*in*/ const int entity_handles_size,
                            /*inout*/ iBase_EntitySetHandle* entity_set,
                            /*out*/ int *err);


    /**\brief  Remove an array of entities from a set
     *
     * Remove an array of entities from a set
     * \param instance iMesh instance handle
     * \param entity_handles Array of entities being remove
     * \param entity_handles_size Number of entities in entity_handles array
     * \param entity_set Pointer to the set being removed from
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntArrFromSet(iMesh_Instance instance,
                              /*in*/ const iBase_EntityHandle* entity_handles,
                              /*in*/ const int entity_handles_size,
                              /*inout*/ iBase_EntitySetHandle* entity_set,
                              /*out*/ int *err);


    /**\brief  Add an entity set to a set
     *
     * Add an entity set to a set
     * \param instance iMesh instance handle
     * \param entity_set_to_add The entity set being added
     * \param entity_set_handle Pointer to the set being added to
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addEntSet(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_set_to_add,
                       /*inout*/ iBase_EntitySetHandle* entity_set_handle,
                       /*out*/ int *err);


    /**\brief  Remove an entity set from a set
     *
     * Remove an entity set from a set
     * \param instance iMesh instance handle
     * \param entity_set_to_remove The entity set being removed
     * \param entity_set_handle Pointer to the set being removed from
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntSet(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle entity_set_to_remove,
                       /*inout*/ iBase_EntitySetHandle* entity_set_handle,
                       /*out*/ int *err);

    /**\brief  Return whether an entity is contained in another set
     *
     * Return whether an entity is contained (*is_contained=1) or not 
     * contained (*is_contained=0) in another set
     * \param instance iMesh instance handle
     * \param containing_entity_set Entity set being queried
     * \param contained_entity Entity potentially contained in 
     *        containing_entity_set
     * \param is_contained Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isEntContained(iMesh_Instance instance,
                            /*in*/ const iBase_EntitySetHandle containing_entity_set,
                            /*in*/ const iBase_EntitySetHandle contained_entity,
                            /*out*/ int *is_contained,
                            /*out*/ int *err);
                            
  void iMesh_isEntArrContained( iMesh_Instance instance,
                            /*in*/ iBase_EntitySetHandle containing_set,
                            /*in*/ const iBase_EntitySetHandle* entity_handles,
                            /*in*/ int num_entity_handles,
                         /*inout*/ int** is_contained,
                         /*inout*/ int* is_contained_allocated,
                           /*out*/ int* is_contained_size,
                           /*out*/ int* err );

    /**\brief  Return whether an entity set is contained in another set
     *
     * Return whether a set is contained (*is_contained=1) or not contained
     * (*is_contained=0) in another set
     * \param instance iMesh instance handle
     * \param containing_entity_set Entity set being queried
     * \param contained_entity_set Entity set potentially contained in 
     *        containing_entity_set
     * \param is_contained Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isEntSetContained(iMesh_Instance instance,
                               /*in*/ const iBase_EntitySetHandle containing_entity_set,
                               /*in*/ const iBase_EntitySetHandle contained_entity_set,
                               /*out*/ int *is_contained,
                               /*out*/ int *err);

    /**\brief  Add parent/child links between two sets
     *
     * Add parent/child links between two sets.  Makes parent point to child
     * and child point to parent.
     * \param instance iMesh instance handle
     * \param parent_entity_set Pointer to parent set
     * \param child_entity_set Pointer to child set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_addPrntChld(iMesh_Instance instance,
                         /*inout*/ iBase_EntitySetHandle* parent_entity_set,
                         /*inout*/ iBase_EntitySetHandle* child_entity_set,
                         /*out*/ int *err);

    /**\brief  Remove parent/child links between two sets
     *
     * Remove parent/child links between two sets.
     * \param instance iMesh instance handle
     * \param parent_entity_set Pointer to parent set
     * \param child_entity_set Pointer to child set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvPrntChld(iMesh_Instance instance,
                         /*inout*/ iBase_EntitySetHandle* parent_entity_set,
                         /*inout*/ iBase_EntitySetHandle* child_entity_set,
                         /*out*/ int *err);

    /**\brief  Return whether two sets are related by parent/child links
     *
     * Return whether two sets are related (*is_child=1) or not (*is_child=0)
     * by parent/child links
     * \param instance iMesh instance handle
     * \param parent_entity_set Pointer to parent set
     * \param child_entity_set Pointer to child set
     * \param is_child Pointer to flag returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_isChildOf(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle parent_entity_set,
                       /*in*/ const iBase_EntitySetHandle child_entity_set,
                       /*out*/ int *is_child,
                       /*out*/ int *err);

    /**\brief  Get the number of child sets linked from a specified set
     *
     * Get the number of child sets linked from a specified set.  If num_hops
     * is non-zero, this represents the maximum hops from entity_set to any
     * child in the count.
     * \param instance iMesh instance handle
     * \param entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to child set,
     *        inclusive of the child set
     * \param num_child Pointer to number of children returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumChld(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        /*out*/ int *num_child,
                        /*out*/ int *err);

    /**\brief  Get the number of parent sets linked from a specified set
     *
     * Get the number of parent sets linked from a specified set.  If num_hops
     * is non-zero, this represents the maximum hops from entity_set to any
     * parent in the count.
     * \param instance iMesh instance handle
     * \param entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to parent set,
     *        inclusive of the parent set
     * \param num_parent Pointer to number of parents returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNumPrnt(iMesh_Instance instance,
                        /*in*/ const iBase_EntitySetHandle entity_set,
                        /*in*/ const int num_hops,
                        /*out*/ int *num_parent,
                        /*out*/ int *err);

    /**\brief  Get the child sets linked from a specified set
     *
     * Get the child sets linked from a specified set.  If num_hops
     * is non-zero, this represents the maximum hops from entity_set to any
     * child.
     * \param instance iMesh instance handle
     * \param from_entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to child set,
     *        inclusive of the child set
     * \param *entity_set_handles Pointer to array of child sets
     *        returned from function
     * \param *entity_set_handles_allocated Pointer to allocated size of 
     *        entity_set_handles array
     * \param *entity_set_handles_size Pointer to occupied size of 
     *        entity_set_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getChldn(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*out*/ iBase_EntitySetHandle** entity_set_handles,
                      /*out*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size,
                      /*out*/ int *err);

    /**\brief  Get the parent sets linked from a specified set
     *
     * Get the parent sets linked from a specified set.  If num_hops
     * is non-zero, this represents the maximum hops from entity_set to any
     * parent.
     * \param instance iMesh instance handle
     * \param from_entity_set Entity set being queried
     * \param num_hops Maximum hops from entity_set_handle to parent set,
     *        inclusive of the parent set
     * \param *entity_set_handles Pointer to array of parent sets
     *        returned from function
     * \param *entity_set_handles_allocated Pointer to allocated size of 
     *        entity_set_handles array
     * \param *entity_set_handles_size Pointer to occupied size of 
     *        entity_set_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getPrnts(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle from_entity_set,
                      /*in*/ const int num_hops,
                      /*out*/ iBase_EntitySetHandle** entity_set_handles,
                      /*out*/ int* entity_set_handles_allocated,
                      /*out*/ int* entity_set_handles_size,
                      /*out*/ int *err);

    /**\brief  Set coordinates for an array of vertices
     *
     * Set coordinates for an array of vertices.  Specified storage 
     * order must be either iBase_INTERLEAVED or iBase_BLOCKED, and 
     * indicates order of x, y, and z coordinates in coordinate array.
     * \param instance iMesh instance handle
     * \param vertex_handles Array of vertex handles
     * \param vertex_handles_size Number of vertex handles in array
     * \param storage_order Storage order of coordinates in coordinate array
     * \param new_coords Coordinate array
     * \param new_coords_size Size of coordinate array; should be 
     *        3*vertex_handles_size
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setVtxArrCoords(iMesh_Instance instance,
                             /*in*/ iBase_EntityHandle* vertex_handles,
                             /*in*/ const int vertex_handles_size,
                             /*in*/ const int storage_order,
                             /*in*/ const double* new_coords,
                             /*in*/ const int new_coords_size,
                             /*out*/ int *err);


    /**\brief  Create an array of new vertices at specified coordinates
     *
     * Create an array of new vertices at specified coordinates.  Value of
     * storage_order must be either iBase_INTERLEAVED or iBase_BLOCKED.
     * \param instance iMesh instance handle
     * \param num_verts Number of new vertices to be created
     * \param storage_order Storage order of coordinates in new_coords array
     * \param new_coords Array of coordinates of new vertices
     * \param new_coords_size Number of coordinates in new_coords array, should
     *        be 3*num_verts
     * \param *new_vertex_handles Pointer to array of new vertex handles 
     *        returned from function
     * \param *new_vertex_handles_allocated Pointer to allocated size of 
     *        new_vertex_handles array
     * \param *new_vertex_handles_size Pointer to occupied size of 
     *        new_vertex_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createVtxArr(iMesh_Instance instance,
                          /*in*/ const int num_verts,
                          /*in*/ const int storage_order,
                          /*in*/ const double* new_coords,
                          /*in*/ const int new_coords_size,
                          /*inout*/ iBase_EntityHandle** new_vertex_handles,
                          /*inout*/ int* new_vertex_handles_allocated,
                          /*inout*/ int* new_vertex_handles_size,
                          /*out*/ int *err);


    /**\brief  Create an array of new entities with specified lower-order topology
     *
     * Create an array of new entities with specified lower-order topology.  
     * Specified new_entity_topology must be value in iMesh_EntityTopology
     * enumeration.  Values return in status array must be values in the
     * iBase_CreationStatus enumeration.
     * \param instance iMesh instance handle
     * \param new_entity_topology Topology of created entity
     * \param lower_order_entity_handles Array of lower order entity handles
     *        used to construct new entities
     * \param lower_order_entity_handles_size Number of entities in array of 
     *        lower order entity handles
     * \param *new_entity_handles Pointer to array of new_entity_handles 
     *        returned from function
     * \param *new_entity_handles_allocated Pointer to allocated size of 
     *        new_entity_handles array
     * \param *new_entity_handles_size Pointer to occupied size of 
     *        new_entity_handles array
     * \param *status Pointer to array of creation status returned from 
     *        function
     * \param *status_allocated Pointer to allocated size of status array
     * \param *status_size Pointer to occupied size of status array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createEntArr(iMesh_Instance instance,
                          /*in*/ const int new_entity_topology,
                          /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                          /*in*/ const int lower_order_entity_handles_size,
                          /*out*/ iBase_EntityHandle** new_entity_handles,
                          /*out*/ int* new_entity_handles_allocated,
                          /*out*/ int* new_entity_handles_size,
                          /*inout*/ int** status,
                          /*inout*/ int* status_allocated,
                          /*out*/ int* status_size,
                          /*out*/ int *err);


    /**\brief  Delete specified entities
     *
     * Delete specified entities
     * \param instance iMesh instance handle
     * \param entity_handles Array of entity handles to be deleted
     * \param entity_handles_size Number of entities in array to be deleted
     * \param *err Pointer to error type returned from function
     */
  void iMesh_deleteEntArr(iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*out*/ int *err);


    /**\brief  Create a tag with specified name, size, and type
     *
     * Create a tag with specified name, size, and type.  Tag size is in
     * units of size of tag_type data types.  Value input for tag type must be 
     * value in iBase_TagType enumeration.
     * \param instance iMesh instance handle
     * \param tag_name Character string indicating tag name
     * \param tag_size Size of each tag value, in units of number of tag_type 
     *        entities
     * \param tag_type Data type for data stored in this tag
     * \param tag_handle Pointer to tag handle returned from function
     * \param *err Pointer to error type returned from function
     * \param tag_name_len Length of tag name string
     */
  void iMesh_createTag(iMesh_Instance instance,
                       /*in*/ const char* tag_name,
                       /*in*/ const int tag_size,
                       /*in*/ const int tag_type,
                       /*out*/ iBase_TagHandle* tag_handle, 
                       /*out*/ int *err,
                       /*in*/ const int tag_name_len);


    /**\brief  Destroy a tag
     *
     * Destroy a tag.  If forced is non-zero and entities still have values
     * set for this tag, tag is deleted anyway and those values disappear,
     * otherwise tag is not deleted.
     * \param instance iMesh instance handle
     * \param tag_handle Handle of tag to be deleted
     * \param forced If non-zero, delete the tag even if entities have values
     *        set for that tag
     * \param *err Pointer to error type returned from function
     */
  void iMesh_destroyTag(iMesh_Instance instance,
                        /*in*/ iBase_TagHandle tag_handle,
                        /*in*/ const int forced,
                        /*out*/ int *err);

    /**\brief  Get the name for a given tag handle
     *
     * Get the name for a given tag handle
     * \param instance iMesh instance handle
     * \param tag_handle Tag handle being queried
     * \param name Pointer to character string to store name returned from 
     *        function
     * \param *err Pointer to error type returned from function
     * \param name_len Length of character string input to function
     */
  void iMesh_getTagName(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*inout*/ char *name,
                        /*out*/ int *err,
                        /*in*/ int name_len);

    /**\brief  Get size of a tag in units of numbers of tag data type
     *
     * Get size of a tag in units of numbers of tag data type
     * \param instance iMesh instance handle
     * \param tag_handle Handle of tag being queried
     * \param tag_size Pointer to tag size returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getTagSizeValues(iMesh_Instance instance,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ int *tag_size,
                              /*out*/ int *err);

    /**\brief  Get size of a tag in units of bytes
     *
     * Get size of a tag in units of bytes
     * \param instance iMesh instance handle
     * \param tag_handle Handle of tag being queried
     * \param tag_size Pointer to tag size returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getTagSizeBytes(iMesh_Instance instance,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ int *tag_size,
                             /*out*/ int *err);

    /**\brief  Get a the handle of an existing tag with the specified name
     *
     * Get a the handle of an existing tag with the specified name
     * \param instance iMesh instance handle
     * \param tag_name Name of tag being queried
     * \param tag_handle Pointer to tag handle returned from function
     * \param *err Pointer to error type returned from function
     * \param tag_name_len Length of tag name string
     */
  void iMesh_getTagHandle(iMesh_Instance instance,
                          /*in*/ const char* tag_name,
                          /*out*/ iBase_TagHandle *tag_handle, 
                          /*out*/ int *err,
                          int tag_name_len);

    /**\brief  Get the data type of the specified tag handle
     *
     * Get the data type of the specified tag handle.  Tag type is a value in
     * the iBase_TagType enumeration.
     * \param instance iMesh instance handle
     * \param tag_handle Handle for the tag being queried
     * \param tag_type Pointer to tag type returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getTagType(iMesh_Instance instance,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ int *tag_type,
                        /*out*/ int *err);

    /**\brief  Set a tag value of arbitrary type on an entity set
     *
     * Set a tag value of arbitrary type on an entity set.  Tag data is 
     * passed as char* type,
     * but really represents pointer to arbitrary data.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Pointer to tag data being set on entity set
     * \param tag_value_size Size in bytes of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetData(iMesh_Instance instance,
                           /*in*/ iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const char* tag_value,
                           /*in*/ const int tag_value_size,
                           /*out*/ int *err);


    /**\brief  Set a tag value of integer type on an entity set
     *
     * Set a tag value of integer type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetIntData(iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const int tag_value,
                              /*out*/ int *err);


    /**\brief  Set a tag value of double type on an entity set
     *
     * Set a tag value of double type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetDblData(iMesh_Instance instance,
                              /*in*/ iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*in*/ const double tag_value,
                              /*out*/ int *err);


    /**\brief  Set a tag value of entity handle type on an entity set
     *
     * Set a tag value of entity handle type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param tag_value Tag value being set on entity set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEntSetEHData(iMesh_Instance instance,
                             /*in*/ iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*in*/ const iBase_EntityHandle tag_value,
                             /*out*/ int *err);


    /**\brief  Get the value of a tag of arbitrary type on an entity set
     *
     * Get the value of a tag of arbitrary type on an entity set.  Tag data 
     * is passed back as char* type, but really represents arbitrary data.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *tag_value Pointer to tag data array being queried
     * \param *tag_value_allocated Pointer to tag data array allocated size
     * \param *tag_value_size Pointer to tag data array occupied size
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetData(iMesh_Instance instance,
                           /*in*/ const iBase_EntitySetHandle entity_set_handle,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ char** tag_value,
                           /*inout*/ int* tag_value_allocated,
                           /*inout*/ int* tag_value_size,
                           /*out*/ int *err);

    /**\brief  Get the value of a tag of integer type on an entity set
     *
     * Get the value of a tag of integer type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetIntData(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ int *out_data,
                              /*out*/ int *err);

    /**\brief  Get the value of a tag of double type on an entity set
     *
     * Get the value of a tag of double type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetDblData(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set,
                              /*in*/ const iBase_TagHandle tag_handle,
                              /*out*/ double *out_data,
                              /*out*/ int *err);

    /**\brief  Get the value of a tag of entity handle type on an entity set
     *
     * Get the value of a tag of entity handle type on an entity set.
     * \param instance iMesh instance handle
     * \param entity_set Entity set on which tag is being set
     * \param tag_handle Tag being set on an entity set
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntSetEHData(iMesh_Instance instance,
                             /*in*/ const iBase_EntitySetHandle entity_set,
                             /*in*/ const iBase_TagHandle tag_handle,
                             /*out*/ iBase_EntityHandle *out_data,
                             /*out*/ int *err);

    /**\brief  Get all the tags associated with a specified entity set
     *
     * Get all the tags associated with a specified entity set
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity being queried
     * \param *tag_handles Pointer to array of tag_handles returned from 
     *        function
     * \param *tag_handles_allocated Pointer to allocated size of tag_handles 
     *        array
     * \param *tag_handles_size Pointer to occupied size of tag_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAllEntSetTags(iMesh_Instance instance,
                              /*in*/ const iBase_EntitySetHandle entity_set_handle,
                              /*out*/ iBase_TagHandle** tag_handles,
                              /*out*/ int* tag_handles_allocated,
                              /*out*/ int* tag_handles_size,
                              /*out*/ int *err);

    /**\brief  Remove a tag value from an entity set
     *
     * Remove a tag value from an entity set
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set from which tag is being removed
     * \param tag_handle Tag handle of tag being removed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvEntSetTag(iMesh_Instance instance,
                          /*in*/ iBase_EntitySetHandle entity_set_handle,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*out*/ int *err);

    /**\brief  Set coordinates for a vertex
     *
     * Set coordinates for a vertex.
     * \param instance iMesh instance handle
     * \param vertex_handle vertex handle being set
     * \param x x coordinate being set
     * \param y y coordinate being set
     * \param z z coordinate being set
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setVtxCoords(iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle vertex_handle,
                          /*in*/ const double x, /*in*/ const double y,
                          /*in*/ const double z,
                          /*out*/ int *err);

    /**\brief  Create a new vertex at specified coordinates
     *
     * Create a new vertex at specified coordinates.
     * \param instance iMesh instance handle
     * \param x x coordinate of new vertex
     * \param y y coordinate of new vertex
     * \param z z coordinate of new vertex
     * \param new_vertex_handle Pointer to new vertex handles returned from 
     *        function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createVtx(iMesh_Instance instance,
                       /*in*/ const double x, /*in*/ const double y,
                       /*in*/ const double z,
                       /*out*/ iBase_EntityHandle* new_vertex_handle,
                       /*out*/ int *err);

    /**\brief  Create a new entity with specified lower-order topology
     *
     * Create a new entity with specified lower-order topology.  
     * Specified new_entity_topology must be value in iMesh_EntityTopology
     * enumeration.  Value returned as status must be a value in the
     * iBase_CreationStatus enumeration.
     * \param instance iMesh instance handle
     * \param new_entity_topology Topology of created entity
     * \param lower_order_entity_handles Array of lower order entity handles
     *        used to construct new entity
     * \param lower_order_entity_handles_size Number of entities in array of 
     *        lower order entity handles
     * \param new_entity_handle Pointer to new entity handle returned from 
     *        function
     * \param status Pointer to creation status returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_createEnt(iMesh_Instance instance,
                       /*in*/ const int new_entity_topology,
                       /*in*/ const iBase_EntityHandle* lower_order_entity_handles,
                       /*in*/ const int lower_order_entity_handles_size,
                       /*out*/ iBase_EntityHandle* new_entity_handle,
                       /*out*/ int* status,
                       /*out*/ int *err);

    /**\brief  Delete specified entity
     *
     * Delete specified entity
     * \param instance iMesh instance handle
     * \param entity_handle Entity to be deleted
     * \param *err Pointer to error type returned from function
     */
  void iMesh_deleteEnt(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*out*/ int *err);

    /**\brief  Get tag values of arbitrary type for an array of entities
     *
     * Get tag values of arbitrary type for an array of entities.  Tag data 
     * is returned as char* type, but really represents arbitrary data.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_values Pointer to tag data array being returned from 
     *        function
     * \param tag_values_allocated Pointer to allocated size of tag data array
     * \param tag_values_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getArrData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*inout*/ char** tag_values,
                        /*inout*/int* tag_values_allocated,
                        /*out*/ int* tag_values_size,
                        /*out*/ int *err);

    /**\brief  Get tag values of integer type for an array of entities
     *
     * Get tag values of integer type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_values Pointer to tag data array being returned from 
     *        function
     * \param tag_values_allocated Pointer to allocated size of tag data array
     * \param tag_values_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getIntArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ int** tag_values,
                           /*inout*/ int* tag_values_allocated,
                           /*out*/ int* tag_values_size,
                           /*out*/ int *err);

    /**\brief  Get tag values of double type for an array of entities
     *
     * Get tag values of double type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_values Pointer to tag data array being returned from 
     *        function
     * \param tag_values_allocated Pointer to allocated size of tag data array
     * \param tag_values_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getDblArrData(iMesh_Instance instance,
                           /*in*/ const iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*inout*/ double** tag_values,
                           /*inout*/ int* tag_values_allocated,
                           /*out*/ int* tag_values_size,
                           /*out*/ int *err);

    /**\brief  Get tag values of entity handle type for an array of entities
     *
     * Get tag values of entity handle type for an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param *tag_value Pointer to tag data array being returned from 
     *        function
     * \param tag_value_allocated Pointer to allocated size of tag data array
     * \param tag_value_size Pointer to occupied size of tag data array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEHArrData(iMesh_Instance instance,
                          /*in*/ const iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*inout*/ iBase_EntityHandle** tag_value,
                          /*inout*/ int* tag_value_allocated,
                          /*out*/ int* tag_value_size,
                          /*out*/ int *err);

    /**\brief  Set tag values of arbitrary type on an array of entities
     *
     * Set tag values of arbitrary type on an array of entities.  Tag data is 
     * passed as char* type, but really represents pointer to arbitrary data.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entity
     * \param tag_values_size Size in total bytes of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setArrData(iMesh_Instance instance,
                        /*in*/ iBase_EntityHandle* entity_handles,
                        /*in*/ const int entity_handles_size,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const char* tag_values,
                        /*in*/ const int tag_values_size,
                        /*out*/ int *err);

    /**\brief  Set tag values of integer type on an array of entities
     *
     * Set tag values of integer type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities
     * \param tag_values_size Size in total number of integers of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setIntArrData(iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const int* tag_values,
                           /*in*/ const int tag_values_size,
                           /*out*/ int *err);

    /**\brief  Set tag values of double type on an array of entities
     *
     * Set tag values of double type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities
     * \param tag_values_size Size in total number of doubles of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setDblArrData(iMesh_Instance instance,
                           /*in*/ iBase_EntityHandle* entity_handles,
                           /*in*/ const int entity_handles_size,
                           /*in*/ const iBase_TagHandle tag_handle,
                           /*in*/ const double* tag_values,
                           /*in*/ const int tag_values_size,
                           /*out*/ int *err);

    /**\brief  Set tag values of entity handle type on an array of entities
     *
     * Set tag values of entity handle type on an array of entities.
     * \param instance iMesh instance handle
     * \param entity_handles Entity array on which tag is being set
     * \param entity_handles_size Number of entities in array
     * \param tag_handle Tag being set on an entity
     * \param tag_values Pointer to tag data being set on entities
     * \param tag_values_size Size in total number of entity handles of tag 
     *        data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEHArrData(iMesh_Instance instance,
                          /*in*/ iBase_EntityHandle* entity_handles,
                          /*in*/ const int entity_handles_size,
                          /*in*/ const iBase_TagHandle tag_handle,
                          /*in*/ const iBase_EntityHandle* tag_values,
                          /*in*/ const int tag_values_size,
                          /*out*/ int *err);

    /**\brief  Remove a tag value from an array of entities
     *
     * Remove a tag value from an array of entities
     * \param instance iMesh instance handle
     * \param entity_handles Entity from which tag is being removed
     * \param entity_handles_size Number of entities in entity array
     * \param tag_handle Tag handle of tag being removed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvArrTag(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle* entity_handles,
                       /*in*/ const int entity_handles_size,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ int *err);

    /**\brief  Get the value of a tag of arbitrary type on an entity
     *
     * Get the value of a tag of arbitrary type on an entity.  Tag data 
     * is passed back as char* type, but really represents arbitrary data.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *tag_value Pointer to tag data array being queried
     * \param *tag_value_allocated Pointer to tag data array allocated size
     * \param *tag_value_size Pointer to tag data array occupied size
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getData(iMesh_Instance instance,
                     /*in*/ const iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle,
                     /*inout*/ char** tag_value,
                     /*inout*/ int *tag_value_allocated,
                     /*out*/ int *tag_value_size,
                     /*out*/ int *err);

    /**\brief  Get the value of a tag of integer type on an entity
     *
     * Get the value of a tag of integer type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getIntData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ int *out_data,
                        /*out*/ int *err);

    /**\brief  Get the value of a tag of double type on an entity
     *
     * Get the value of a tag of double type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getDblData(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*out*/ double *out_data,
                        /*out*/ int *err);

    /**\brief  Get the value of a tag of entity handle type on an entity
     *
     * Get the value of a tag of entity handle type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param *out_data Pointer to tag value returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEHData(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*out*/ iBase_EntityHandle *out_data,
                       /*out*/ int *err);

    /**\brief  Set a tag value of arbitrary type on an entity
     *
     * Set a tag value of arbitrary type on an entity.  Tag data is 
     * passed as char* type, but really represents pointer to arbitrary data.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Pointer to tag data being set on entity
     * \param tag_value_size Size in bytes of tag data
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setData(iMesh_Instance instance,
                     /*in*/ iBase_EntityHandle entity_handle,
                     /*in*/ const iBase_TagHandle tag_handle,
                     /*in*/ const char* tag_value,
                     /*in*/ const int tag_value_size,
                     /*out*/ int *err);

    /**\brief  Set a tag value of integer type on an entity
     *
     * Set a tag value of integer type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setIntData(iMesh_Instance instance,
                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const int tag_value,
                        /*out*/ int *err);

    /**\brief  Set a tag value of double type on an entity
     *
     * Set a tag value of double type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setDblData(iMesh_Instance instance,

                        /*in*/ iBase_EntityHandle entity_handle,
                        /*in*/ const iBase_TagHandle tag_handle,
                        /*in*/ const double tag_value,
                        /*out*/ int *err);

    /**\brief  Set a tag value of entity handle type on an entity
     *
     * Set a tag value of entity handle type on an entity.
     * \param instance iMesh instance handle
     * \param entity_handle Entity on which tag is being set
     * \param tag_handle Tag being set on an entity
     * \param tag_value Tag value being set on entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_setEHData(iMesh_Instance instance,
                       /*in*/ iBase_EntityHandle entity_handle,
                       /*in*/ const iBase_TagHandle tag_handle,
                       /*in*/ const iBase_EntityHandle tag_value,
                       /*out*/ int *err);

    /**\brief  Get all the tags associated with a specified entity handle
     *
     * Get all the tags associated with a specified entity handle
     * \param instance iMesh instance handle
     * \param entity_handle Entity being queried
     * \param *tag_handles Pointer to array of tag_handles returned from 
     *        function
     * \param *tag_handles_allocated Pointer to allocated size of tag_handles 
     *        array
     * \param *tag_handles_size Pointer to occupied size of tag_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getAllTags(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*inout*/ iBase_TagHandle** tag_handles,
                        /*inout*/ int* tag_handles_allocated,
                        /*out*/ int* tag_handles_size,
                        /*out*/ int *err);

    /**\brief  Remove a tag value from an entity
     *
     * Remove a tag value from an entity
     * \param instance iMesh instance handle
     * \param entity_handle Entity from which tag is being removed
     * \param tag_handle Tag handle of tag being removed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_rmvTag(iMesh_Instance instance,
                    /*in*/ iBase_EntityHandle entity_handle,
                    /*in*/ const iBase_TagHandle tag_handle,
                    /*out*/ int *err);

    /**\brief  Initialize an iterator over specified entity type, topology, and size
     *
     * Initialize an iterator over specified entity type, topology, and size,
     * for a specified set or instance.  Iterator returned can be used as input
     * to functions returning the entity for the iterator.  If all entities of 
     * a specified type and/or topology are to be iterated, specify 
     * iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  Specified type 
     * or topology must be a value in the iBase_EntityType or 
     * iMesh_EntityTopology enumerations, respectively.
     * \param instance iMesh instance handle
     * \param entity_set_handle Entity set being iterated
     * \param requested_entity_type Type of entity to iterate
     * \param requested_entity_topology Topology of entity to iterate
     * \param entity_iterator Pointer to iterator returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_initEntIter(iMesh_Instance instance,
                         /*in*/ const iBase_EntitySetHandle entity_set_handle,
                         /*in*/ const int requested_entity_type,
                         /*in*/ const int requested_entity_topology,
                         /*out*/ iMesh_EntityIterator* entity_iterator,
                         /*out*/ int *err);

    /**\brief  Get entity corresponding to an iterator and increment iterator
     *
     * Get the entity corresponding to an array iterator, and increment the 
     * iterator.  Also return whether the next value of the iterator has
     * an entity (if non-zero, next iterator value is the end of the
     * iteration).
     * \param instance iMesh instance handle
     * \param entity_iterator Iterator being queried
     * \param entity_handle Pointer to an entity handle corresponding to the
     *        current value of iterator
     * \param has_data Pointer to flag; if returned non-zero, next iterator
     *        has an entity
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getNextEntIter(iMesh_Instance instance,
                            /*in*/ iMesh_EntityIterator entity_iterator,
                            /*out*/ iBase_EntityHandle* entity_handle,
                            /*out*/ int *has_data,
                            /*out*/ int *err);

    /**\brief  Reset the iterator
     *
     * Reset the iterator
     * \param instance iMesh instance handle
     * \param entity_iterator Iterator to reset
     * \param *err Pointer to error type returned from function
     */
  void iMesh_resetEntIter(iMesh_Instance instance,
                          /*in*/ iMesh_EntityIterator entity_iterator,
                          /*out*/ int *err);

    /**\brief  Destroy the specified iterator
     *
     * Destroy the specified iterator
     * \param instance iMesh instance handle
     * \param entity_iterator Iterator which gets destroyed
     * \param *err Pointer to error type returned from function
     */
  void iMesh_endEntIter(iMesh_Instance instance,
                        /*in*/ iMesh_EntityIterator entity_iterator,
                        /*out*/ int *err);

    /**\brief  Get the entity topology for the specified entity
     *
     * Get the entity topology for the specified entity.  Topology
     * returned is a value in the iMesh_EntityTopology enumeration.
     * \param instance iMesh instance handle
     * \param entity_handle Entity handle being queried
     * \param *out_topo Pointer to entity topology returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntTopo(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*out*/ int *out_topo,
                        /*out*/ int *err);

    /**\brief  Get the entity type for the specified entity
     *
     * Get the entity type for the specified entity.  Type returned is a value
     * in the iBase_EntityType enumeration.
     * \param instance iMesh instance handle
     * \param entity_handle Entity handle being queried
     * \param *out_type Pointer to entity type returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntType(iMesh_Instance instance,
                        /*in*/ const iBase_EntityHandle entity_handle,
                        /*out*/ int *out_type,
                        /*out*/ int *err);

    /**\brief  Get coordinates of specified vertex
     *
     * Get coordinates of specified vertex.
     * \param instance iMesh instance handle
     * \param vertex_handle Mesh vertex being queried
     * \param *x Pointer to x coordinate returned from function
     * \param *y Pointer to y coordinate returned from function
     * \param *z Pointer to z coordinate returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getVtxCoord(iMesh_Instance instance,
                         /*in*/ const iBase_EntityHandle vertex_handle,
                         /*out*/ double *x, /*out*/ double *y, /*out*/ double *z,
                         /*out*/ int *err);

    /**\brief  Get entities of specified type adjacent to an entity
     *
     * Get entities of specified type adjacent to an entity.  Specified type
     * must be value in the iBase_EntityType enumeration.
     * \param instance iMesh instance handle
     * \param entity_handle Entity handle being queried
     * \param entity_type_requested Type of adjacent entities requested
     * \param *adj_entity_handles Pointer to array of adjacent entities
     *        returned from function
     * \param *adj_entity_handles_allocated Pointer to allocated size of 
     *        adj_entity_handles array
     * \param *adj_entity_handles_size Pointer to occupied size of 
     *        adj_entity_handles array
     * \param *err Pointer to error type returned from function
     */
  void iMesh_getEntAdj(iMesh_Instance instance,
                       /*in*/ const iBase_EntityHandle entity_handle,
                       /*in*/ const int entity_type_requested,
                       /*inout*/ iBase_EntityHandle** adj_entity_handles,
                       /*inout*/ int* adj_entity_handles_allocated,
                       /*out*/ int* adj_entity_handles_size,
                       /*out*/ int *err);

/**\brief Get "2nd order" adjacencies to an entity
 * Get "2nd order" adjacencies to an entity, that is, from an entity, through
 * other entities of a specified "bridge" dimension, to other entities of another 
 * specified "to" dimension.
 * \param instance iMesh instance for this call
 * \param entity_handle Entity from which adjacencies are requested
 * \param order_adjacent_key Bridge dimension for 2nd order adjacencies
 * \param requested_entity_type Dimension of adjacent entities returned
 * \param *adjacent_entities Pointer to array of adjacent entities
 * \param *adjacent_entities_allocated Pointer to allocated size of returned array
 * \param *adjacent_entities_size Pointer to occupied size of returned array
 * \param *err Pointer to error type returned from function
 */
  void iMesh_getEnt2ndAdj( iMesh_Instance instance,
                           iBase_EntityHandle entity_handle,
                           int order_adjacent_key,
                           int requested_entity_type,
                           iBase_EntityHandle** adjacent_entities,
                           int* adjacent_entities_allocated,
                           int* adjacent_entities_size,
                           int* err );

    /**\brief  Subtract contents of one entity set from another
     *
     * Subtract contents of one entity set from another
     * \param instance iMesh instance handle
     * \param entity_set_1 Entity set from which other set is being subtracted
     * \param entity_set_2 Entity set being subtracted from other set
     * \param result_entity_set Pointer to entity set returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_subtract(iMesh_Instance instance,
                      /*in*/ const iBase_EntitySetHandle entity_set_1,
                      /*in*/ const iBase_EntitySetHandle entity_set_2,
                      /*out*/ iBase_EntitySetHandle* result_entity_set,
                      /*out*/ int *err);

    /**\brief  Intersect contents of one entity set with another
     *
     * Intersect contents of one entity set with another
     * \param instance iMesh instance handle
     * \param entity_set_1 Entity set being intersected with another
     * \param entity_set_2 Entity set being intersected with another
     * \param result_entity_set Pointer to entity set returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_intersect(iMesh_Instance instance,
                       /*in*/ const iBase_EntitySetHandle entity_set_1,
                       /*in*/ const iBase_EntitySetHandle entity_set_2,
                       /*out*/ iBase_EntitySetHandle* result_entity_set,
                       /*out*/ int *err);

    /**\brief  Unite contents of one entity set with another
     *
     * Unite contents of one entity set with another
     * \param instance iMesh instance handle
     * \param entity_set_1 Entity set being united with another
     * \param entity_set_2 Entity set being united with another
     * \param result_entity_set Pointer to entity set returned from function
     * \param *err Pointer to error type returned from function
     */
  void iMesh_unite(iMesh_Instance instance,
                   /*in*/ const iBase_EntitySetHandle entity_set_1,
                   /*in*/ const iBase_EntitySetHandle entity_set_2,
                   /*out*/ iBase_EntitySetHandle* result_entity_set,
                   /*out*/ int *err);

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef IMESH_CBIND_H__ */
