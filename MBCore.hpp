/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifndef MB_IMPL_GENERAL_HPP
#define MB_IMPL_GENERAL_HPP

#define MOAB_IMPL_VERSION 1.01F
#define MOAB_IMPL_VERSION_STRING "1.01"

#include "MBInterface.hpp"
#include "MBProcConfig.hpp"
#include <map>

class MBWriteUtil;
class MBReadUtil;
class AEntityFactory;
class EntitySequenceManager;
class TagServer;
class MBError;
class MBReaderWriterSet;
class MeshSetManager;

#ifdef XPCOM_MB

#define MBCORE_CID \
{ 0x7cb5b7a0, 0x7d7, 0x11d3, { 0xba, 0xb2, 0x0, 0x0, 0x64, 0x65, 0x73, 0x74 } }

#define MBCORE_CONTRACTID "@sandia.gov/MB;1"

#endif


class MBCore : public MBInterface 
{

public:

  //!constructor
  MB_DLL_EXPORT MBCore( int rank = 0, int num_cpu = 1 );

  //!destructor
  MB_DLL_EXPORT ~MBCore();
  
  //! query an MB internal interface
  virtual MBErrorCode query_interface(const std::string& iface_name, void** iface);
 
  //! release an MB internal interface 
  virtual MBErrorCode release_interface(const std::string& iface_name, void* iface);

#if defined(XPCOM_MB)
  // this macro expands to all the nsISupports interface functions
  NS_DECL_ISUPPORTS
#endif

  virtual int QueryInterface (const MBuuid& uuid, MBUnknownInterface** iface );

    //! Returns the major.minor version number of the implementation
    /**
       \param iface_name If non-NULL, will be filled in with a string, possibly 
       containing implementation-specific information
    */
  virtual float impl_version(std::string *version_string = NULL);

  //! get the type from a handle, returns type
  virtual MBEntityType type_from_handle(const MBEntityHandle handle) const;
  
  //! get the id from a handle, returns id
  virtual unsigned int id_from_handle(const MBEntityHandle handle) const;
  
  //! get a handle from an id and type
  virtual MBErrorCode handle_from_id(const MBEntityType type, 
                                      const unsigned int id, 
                                      MBEntityHandle& handle) const;
  
  virtual int dimension_from_handle( const MBEntityHandle ) const;

  //! load mesh from data in file
  //! NOTE: if there is mesh already present, the new mesh will be added
  virtual MBErrorCode load_mesh(const char *file_name,
                                 const int *active_block_id_list = NULL,
                                 const int num_blocks = 0);
  
  virtual MBErrorCode write_mesh(const char *file_name,
                                  const MBEntityHandle *output_list = NULL,
                                  const int num_sets = 0);

  //! deletes all mesh entities from this datastore
  virtual MBErrorCode delete_mesh();

  //! get overall geometric dimension
  virtual MBErrorCode get_dimension(int &dim) const;

  //! set overall geometric dimension
  /** Returns error if setting to 3 dimensions, mesh has been created, and 
   *  there are only 2 dimensions on that mesh
   */
  virtual MBErrorCode set_dimension(const int dim);

  //! get blocked vertex coordinates for all vertices
  /** Blocked = all x, then all y, etc. 
   */
  virtual MBErrorCode get_vertex_coordinates(std::vector<double> &coords) const;

  //! get the coordinate information for this handle if it is of type Vertex
  //! otherwise, return an error
  virtual MBErrorCode  get_coords(const MBRange &entity_handles, 
                                   double *coords) const;
  
  virtual MBErrorCode  get_coords(const MBEntityHandle *entity_handles, 
                                   const int num_entities, 
                                   double *coords) const;
  
  virtual MBErrorCode  get_coords(const MBEntityHandle entity_handle, 
                                   const double *& x, const double *& y, const double *& z) const;
 
  //! set the coordinate information for this handle if it is of type Vertex
  //! otherwise, return an error
  virtual MBErrorCode  set_coords(MBEntityHandle *entity_handles, 
                                   const int num_entities,
                                   const double *coords);

      //! get global connectivity array for specified entity type
      /**  Assumes just vertices, no higher order nodes
       */
    virtual MBErrorCode get_connectivity_by_type(const MBEntityType type, 
                                                  std::vector<MBEntityHandle> &connect) const;

      //! Gets the connectivity for an element MBEntityHandle. 
      /** For non-element handles (ie, MBMeshSets), 
          returns an error. Connectivity data is copied from the database into the vector 
          <em>connectivity</em>. The nodes in <em>connectivity</em> are properly ordered.
          \param entity_handle MBEntityHandle to get connectivity of.
          \param connectivity Vector in which connectivity of <em>entity_handle</em> is returned.  
          Should contain MeshVertices.
          \param topological_connectivity If true, higher order nodes are ignored. 

          Example: \code 
          std::vector<MBEntityHandle> conn;
          get_connectivity( entity_handle, conn ); \endcode */
    virtual MBErrorCode  get_connectivity(const MBEntityHandle *entity_handles, 
                                           const int num_handles,
                                           std::vector<MBEntityHandle> &connectivity, 
                                           bool topological_connectivity = false) const;
 
    //! Gets the connectivity for a vector of elements
    /** Same as vector-based version except range is returned (unordered!)
     */
  virtual MBErrorCode  get_connectivity(const MBEntityHandle *entity_handles, 
                                        const int num_handles,
                                        MBRange &connectivity, 
                                        bool topological_connectivity = false) const;
 
    //! Gets a pointer to constant connectivity data of <em>entity_handle</em> 
      /** Sets <em>number_nodes</em> equal to the number of nodes of the <em> 
          entity_handle </em>.  Faster then the other <em>get_connectivity</em> function. 
          The nodes in 'connectivity' are properly ordered. 
          \param entity_handle MBEntityHandle to get connectivity of.
          \param connectivity Array in which connectivity of <em>entity_handle</em> is returned.
          Should contain MeshVertex's.
          \param num_nodes Number of MeshVertices in array <em>connectivity</em>. 

          Example: \code 
          const MBEntityHandle* conn;
          int number_nodes = 0;
          get_connectivity( entity_handle, conn, number_nodes ); \endcode */
    virtual MBErrorCode  get_connectivity(const MBEntityHandle entity_handle, 
                                           const MBEntityHandle *&connectivity, 
                                           int &num_nodes, 
                                           bool topological_connectivity = false) const;

      //! Sets the connectivity for an MBEntityHandle.  For non-element handles, return an error.
      /** Connectivity is stored exactly as it is ordered in vector <em>connectivity</em>. 
          \param entity_handle MBEntityHandle to set connectivity of.
          \param connect Vector containing new connectivity of <em>entity_handle</em>.
          \param num_connect Number of vertices in <em>connect</em>
   
          Example: \code 
          std::vector<MBEntityHandle> conn(3);
          conn[0] = node1;
          conn[1] = node2;
          conn[2] = node3;
          set_connectivity( entity_handle, conn, 3 ); \endcode */
    virtual MBErrorCode  set_connectivity(const MBEntityHandle entity_handle, 
                                          MBEntityHandle *connect,
                                          const int num_connect);

      //! get the adjacencies associated with a set of entities
      /** \param from_entities vector of MBEntityHandle to get adjacencies of.
          \param to_dimension Dimension of desired adjacency information.
          \param adj_entities Vector in which adjacent MBEntityHandles are returned. 
          \param operation_type enum of INTERSECT or UNION.  Defines whether to take
          the intersection or union of the set of adjacencies recovered for the from_entities.

          The adjacent entities in vector <em>adjacencies</em> are not in any particular 
          order. 

          Example: \code
            // get the set of edges that are adjacent to all entities in the from_entities list
            std::vector<MBEntityHandle> from_entities = {hex1, hex2};
            std::vector<MBEntityHandle> adjacencies;
            get_adjacencies( from_entities, MB_1D_ENTITY, adjacencies ); 
            \endcode */

    virtual MBErrorCode get_adjacencies(const MBEntityHandle *from_entities,
                                         const int num_entities,
                                         const int to_dimension,
                                         const bool create_if_missing,
                                         std::vector<MBEntityHandle>& adj_entities,
                                         const int operation_type = MBInterface::INTERSECT);

    virtual MBErrorCode get_adjacencies(const MBEntityHandle *from_entities,
                                        const int num_entities,
                                         const int to_dimension,
                                         const bool create_if_missing,
                                         MBRange &adj_entities,
                                         const int operation_type = MBInterface::INTERSECT);

    virtual MBErrorCode get_adjacencies(const MBRange &from_entities,
                                         const int to_dimension,
                                         const bool create_if_missing,
                                         MBRange &adj_entities,
                                         const int operation_type = MBInterface::INTERSECT);

      //! Adds adjacencies
      /** \param from_handle entities 
          \param both_ways add the adjacency information to both the
          to_handle and and the from_from :handle

          Example: \code
      */
    virtual MBErrorCode add_adjacencies(const MBEntityHandle from_handle, 
                                         const MBEntityHandle *to_handles,
                                         const int num_handles,
                                         bool both_ways);

      //! Adds adjacencies; same as vector-based, but with range instead
    virtual MBErrorCode add_adjacencies(const MBEntityHandle from_handle, 
                                        MBRange &adjacencies,
                                        bool both_ways);

      //! Removes adjacencies
      /** \param handle MBEntityHandle to get adjacencies of.

      Example: \code
      */
    virtual MBErrorCode remove_adjacencies(const MBEntityHandle from_handle, 
                                            const MBEntityHandle *to_handles,
                                            const int num_handles);

      //! Retrieves all entities in the database of given dimension.  
      /** \param dimension Dimension of entities desired.
          \param entities Range in which entities of dimension <em>dimension</em> are returned.

          Example: \code
          int dimension = 2;
          MBRange entities;
          get_entities_by_dimension( dimension, entities ); //get 2D MBEntityHandles in the database
          \endcode */
    virtual MBErrorCode get_entities_by_dimension(const MBEntityHandle meshset,
                                                   const int dimension, 
                                                   MBRange &entities,
                                                   const bool recursive = false) const;

      //! Retrieves all entities in the data base of given type.  
      /** \param type MBEntityType of entities desired (ie, MeshHex, MeshEdge, MeshTri, etc )
          \param entities Range in which entities of MBEntityType <em>type</em> are returned.

          Example: \code
          MBEntityType type = MeshTet;
          MBRange entities;
          get_entities_by_dimension( type, entities ); //get MeshTet type MBEntityHandles in the database
          \endcode */
    virtual MBErrorCode get_entities_by_type(const MBEntityHandle meshset,
                                              const MBEntityType type, 
                                              MBRange &entities,
                                              const bool recursive = false) const;

    virtual MBErrorCode get_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                      const MBEntityType type,
                                                      const MBTag *tag_handles,
                                                      const void* const* values,
                                                      const int num_tags,
                                                      MBRange &entities,
                                                      const int condition = MBInterface::INTERSECT,
                                                      const bool recursive = false) const;

      //! Retrieves all entities in the data base
      /** \param entities Range in which entities of MBEntityType <em>type</em> are returned.

      Example: \code
      MBRange entities;
      get_entities( entities ); //get MeshTet type MBEntityHandles in the database
      \endcode */
    virtual MBErrorCode get_entities_by_handle(const MBEntityHandle meshset,
                                      MBRange &entities,
                                      const bool recursive = false) const;

      //! Retrieves all entities in the data base
      /** \param entities Range in which entities of MBEntityType <em>type</em> are returned.

      Example: \code
      MBRange entities;
      get_entities( entities ); //get MeshTet type MBEntityHandles in the database
      \endcode */
    virtual MBErrorCode get_entities_by_handle(const MBEntityHandle meshset,
                                      std::vector<MBEntityHandle> &entities,
                                      const bool recursive = false) const;

      //! Retrieves all entities in the database of given dimension.  
      /** \param dimension Dimension of entities desired.
          \param entities Range in which entities of dimension <em>dimension</em> are returned.

          Example: \code
          int dimension = 2;
          MBRange entities;
          get_entities_by_dimension( dimension, entities ); //get 2D MBEntityHandles in the database
          \endcode */
    virtual MBErrorCode get_number_entities_by_dimension(const MBEntityHandle meshset,
                                                          const int dimension, 
                                                          int &num_entities,
                                                          const bool recursive = false) const;

      //! Retrieves all entities in the data base of given type.  
      /** \param type MBEntityType of entities desired (ie, MeshHex, MeshEdge, MeshTri, etc )
          \param entities Range in which entities of MBEntityType <em>type</em> are returned.

          Example: \code
          MBEntityType type = MeshTet;
          MBRange entities;
          get_entities_by_dimension( type, entities ); //get MeshTet type MBEntityHandles in the database
          \endcode */
    virtual MBErrorCode get_number_entities_by_type(const MBEntityHandle meshset,
                                                     const MBEntityType type, 
                                                     int &num_entities,
                                                     const bool recursive = false) const;

    virtual MBErrorCode get_number_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                             const MBEntityType type,
                                                             const MBTag *tag_handles,
                                                             const void** values,
                                                             const int num_tags,
                                                             int &num_entities,
                                                             const bool recursive = false) const;

      //! Retrieves all entities in the data base
      /** \param entities Range in which entities of MBEntityType <em>type</em> are returned.

      Example: \code
      MBRange entities;
      get_entities( entities ); //get MeshTet type MBEntityHandles in the database
      \endcode */
    virtual MBErrorCode get_number_entities_by_handle(const MBEntityHandle meshset,
                                             int &num_entities,
                                             const bool recursive = false) const;

      //! Creates an element based on the type and connectivity. 
      /** If connectivity vector is not correct for MBEntityType <em>type</em> (ie, a vector with 
          3 vertices is passed in to make an MeshQuad), the function returns MB_FAILURE. 
          \param type Type of element to create. (MeshTet, MeshTri, MeshKnife, etc.) 
          \param connectivity Vector containing connectivity of element to create.
          \param handle MBEntityHandle representing the newly created element in the database.

          Example: \code
          MBEntityType type = MeshQuad;
          std::vector<MBEntityHandle> connectivity(4);
          quad_conn[0] = vertex0;
          quad_conn[1] = vertex1;
          quad_conn[2] = vertex2;
          quad_conn[3] = vertex3;
          MBEntityHandle element_handle;
          create_element( type, connectivity, element_handle ); \endcode */
    virtual MBErrorCode create_element(const MBEntityType type, 
                                        const MBEntityHandle *connectivity,
                                        const int num_nodes, 
                                        MBEntityHandle &element_handle);

      /**\brief Create element given CPU ID and connectivity.
       *
       * Create an element with the specified processor ID
       *\param type The type of the element
       *\param processor_id The ID of the CPU on owning the element
       *\param connectivity The connectivity list for the element
       *\param num_nodes The length of the connectivity list
       *\param element_handle Output handle value.
       */
    virtual MBErrorCode create_element (const MBEntityType type, 
                                        const unsigned processor_id,
                                        const MBEntityHandle *connectivity,
                                        const int num_nodes, 
                                        MBEntityHandle &element_handle);

      //! Creates a vertex based on coordinates.  
      /**
         \param coordinates Array that has 3 doubles in it.
         \param entity_handle MBEntityHandle representing the newly created vertex in the database.

         Example: \code
         double *coordinates = double[3];
         coordinates[0] = 1.034;
         coordinates[1] = 23.23; 
         coordinates[2] = -0.432; 
         MBEntityHandle entity_handle = 0;
         create_vertex( coordinates, entity_handle ); \endcode */
    virtual MBErrorCode create_vertex(const double coordinates[3], 
                                       MBEntityHandle &entity_handle );

      /**\brief Create vertex given CPU ID and coordinates.
       *
       * Create a vertex with the specified processor ID
       *\param processor_id The ID of the CPU on owning the element
       *\param coordinates The vertex coordinates
       *\param entity_handle Output handle value.
       */
    virtual MBErrorCode create_vertex( const unsigned processor_id,
                                       const double coordinates[3], 
                                       MBEntityHandle &entity_handle );

      //! merges two entities
    virtual MBErrorCode merge_entities(MBEntityHandle entity_to_keep, 
                                        MBEntityHandle entity_to_remove,
                                        bool auto_merge,
                                        bool delete_removed_entity);

      //! Removes entities in a vector from the data base.  
      /** If any of the entities are contained in any meshsets, it is removed from those meshsets 
          which were created with MESHSET_TRACK_OWNER option bit set.  Tags for <em>entity<\em> are 
          removed as part of this function.
          \param entities 1d vector of entities to delete
          \param num_entities Number of entities in 1d vector
      */ 
    virtual MBErrorCode delete_entities(const MBEntityHandle *entities,
                                         const int num_entities);

      //! Removes entities in a range from the data base.  
      /** If any of the entities are contained in any meshsets, it is removed from those meshsets 
          which were created with MESHSET_TRACK_OWNER option bit set.  Tags for <em>entity<\em> are 
          removed as part of this function.
          \param entities Range of entities to delete
      */ 
    virtual MBErrorCode delete_entities(const MBRange &entities);

  virtual MBErrorCode list_entities(const MBRange &entities) const;
  
  virtual MBErrorCode list_entities(const MBEntityHandle *entities,
                                    const int num_entities) const;

  virtual MBErrorCode list_entity(const MBEntityHandle entity) const;

      //! function object for recieving events from MB of higher order nodes
      //! added to entities
    class HONodeAddedRemoved
    {
    public:
      HONodeAddedRemoved(){}
      virtual ~HONodeAddedRemoved(){}
        //! node_added called when a node was added to an element's connectivity array
        //! note: connectivity array of element may be incomplete (corner nodes will exist always)
      virtual void node_added(MBEntityHandle node, MBEntityHandle element);
      virtual void node_removed(MBEntityHandle node);
    };
  
    virtual MBErrorCode convert_entities(const MBEntityHandle meshset, 
                                          const bool mid_side,
                                          const bool mid_face, 
                                          const bool mid_volume, 
                                          MBInterface::HONodeAddedRemoved* function_object = 0);

      //! function to get the side number given two elements; returns
      //! MB_FAILURE if child not related to parent; does *not* create adjacencies
      //! between parent and child
    virtual MBErrorCode side_number(const MBEntityHandle parent,
                                     const MBEntityHandle child,
                                     int &side_number,
                                     int &sense,
                                     int &offset) const;

      //! given an entity and the connectivity and type of one of its subfacets, find the
      //! high order node on that subfacet, if any
    virtual MBErrorCode high_order_node(const MBEntityHandle parent_handle,
                                         const MBEntityHandle *subfacet_conn,
                                         const MBEntityType subfacet_type,
                                         MBEntityHandle &high_order_node) const;

      //! given an entity and a target dimension & side number, get that entity
    virtual MBErrorCode side_element(const MBEntityHandle source_entity,
                                      const int dim, 
                                      const int side_number,
                                      MBEntityHandle &target_entity) const;

      //-------------------------Tag Stuff-------------------------------------//

  //! Creates a dense tag with a name.
  /** Use to store data that is larger than 8 bits, on many 
      MBEntityHandles across all entity types.  Allows for storage of data of 
      <em>tag_size</em> bytes on any arbitrary entity.  
      \param tag_name String name of MBTag.
      \param tag_size Size of data to store on tag, in bytes.  For storing data 
      1 byte or less in size, use tag_create_bits(...)
      \param tag_handle MBTag to be created.
      \param default_value Default value tag data is set to when initially created.

      Example: \code
      std::string tag_name = "my_meshset_tag";
      int tag_size = sizeof(double); 
      MBTag tag_handle = 0;
      double value = 100.5;
      const void *default_value = &value;
      tag_create_dense( tag_name, 
              tag_size, default_value ); //Create a dense tag. 
                                         //The tag will hold data the size of a 
                                         //double and that data will initially be 
                                         //set to 100.5  \endcode */
  virtual MBErrorCode tag_create(const char *tag_name,
                                  const int tag_size, 
                                  const MBTagType type,
                                  MBTag &tag_handle, 
                                  const void *default_value);

    /** \brief Define a new tag.
     *
     * Define a new tag for storing application-defined data on MB entities.  
     *
     * \param name    The name of the tag.
     * \param size    The size of the tag data in bytes.
     * \param storage The tag storage type.
     * \param data    The tag data type.
     * \param handle  The tag handle (output)
     * \param def_val Optional default value for tag.
     * \param use_existing  If true, and a tag with the same name and
     *                same description exists, successfully return the
     *                handle of the tag.
     * \return - MB_ALREADY_ALLOCATED if a tag with name already exists.
     *         - MB_FAILURE if inconsistant arguments
     *         - MB_SUCCESS otherwise.
     */
  virtual MBErrorCode tag_create( const      char* name,
                                  const        int size,
                                  const  MBTagType storage,
                                  const MBDataType data,
                                            MBTag& handle,
                                  const      void* def_val,
                                              bool use_existing);

  //! Gets the tag name string of the tag_handle.
  /** \param tag_handle MBTag you want the name of.  
      \param tag_name Name string of <em>tag_handle</em>. 

      Example: \code
      MBTag tag_handle = 0;
      std::string tag_name = "my_special_tag";
      tag_get_name( tag_handle, tag_name );  //gets the MBTag from the tag's name string
      \endcode */
  virtual MBErrorCode  tag_get_name(const MBTag tag_handle, 
                                     std::string& tag_name) const;

  //! Gets tag handle from the tag's string name. 
  /**
      \param tag_name Name string of desired tag. 
      \param tag_handle MBTag to be retrieved.

      Example: \code
      MBTag tag_handle = 0;
      std::string tag_name = "quad_data_flag";
      tag_get_handle( tag_name, tag_handle ); \endcode */ 
  virtual MBErrorCode  tag_get_handle(const char *tag_name, 
                                       MBTag &tag_handle) const;

  //! Get handles for all tags defined on this entity
  virtual MBErrorCode tag_get_tags_on_entity(const MBEntityHandle entity,
                                             std::vector<MBTag> &tag_handles) const;

  //! get size of tag in bytes
  virtual MBErrorCode tag_get_size(const MBTag, int &tag_size) const;

    //! Get the default value of the specified tag
  virtual MBErrorCode tag_get_default_value(const MBTag tag, void *def_val) const;

  //! get type of tag (sparse, dense, etc.; 0 = dense, 1 = sparse, 2 = bit, 3 = mesh)
  virtual MBErrorCode tag_get_type(const MBTag, MBTagType &tag_type) const;

   /** \brief Get data type of tag.
    *
    * Get the type of the tag data.  The tag is data is assumed to
    * be a vector of this type.  If the tag data vetcor contains 
    * more than one value, then the tag size must be a multiple of
    * the size of this type.
    * \param tag  The tag 
    * \param type The type of the specified tag (output).
    */
  virtual MBErrorCode tag_get_data_type(const MBTag tag, MBDataType& type) const;

  //! get handles for all tags defined
  virtual MBErrorCode tag_get_tags(std::vector<MBTag> &tag_handles) const;

  virtual MBErrorCode  tag_get_data(const MBTag tag_handle, 
                                     const MBEntityHandle* entity_handles, 
                                     const int num_entities, 
                                     void *tag_data) const;

  virtual MBErrorCode  tag_get_data(const MBTag tag_handle, 
                                     const MBRange& entity_handles, 
                                     void *tag_data) const;

  //! Sets the data of a given EntityHandle and MBTag.  
  /** If the <em>tag_handle</em> and the entity type of <em>entity_handle</em> are not 
      compatible, data of <em>entity_handle</em> never existed and MB_FAILURE 
      is returned. 
      \param tag_handle MBTag indicating what data is to be set.
      \param entity_handle MBEntityHandle on which to set tag's data. 
      \param tag_data Data to set the <em>entity_handle</em>'s tag data to.

      Example: \code
      int tag_data = 1004;
      tag_set_data( tag_handle, entity_handle, &tag_data ); \endcode */
  virtual MBErrorCode  tag_set_data(const MBTag tag_handle, 
                                     const MBEntityHandle* entity_handles, 
                                     const int num_entities,
                                     const void *tag_data );
  
  virtual MBErrorCode  tag_set_data(const MBTag tag_handle, 
                                     const MBRange& entity_handles,
                                     const void *tag_data );

  //! Delete the data of a vector of entity handles and sparse tag
  /** Delete the data of a tag on a vector of entity handles.  Only sparse tag data are deleted with this
      function; dense tags are deleted by deleting the tag itself using tag_delete.
      \param tag_handle Handle of the (sparse) tag being deleted from entity
      \param entity_handles 1d vector of entity handles from which the tag is being deleted
      \param num_handles Number of entity handles in 1d vector
  */
  virtual MBErrorCode  tag_delete_data(const MBTag tag_handle, 
                                        const MBEntityHandle *entity_handles,
                                        const int num_handles);

  //! Delete the data of a range of entity handles and sparse tag
  /** Delete the data of a tag on a range of entity handles.  Only sparse tag data are deleted with this
      function; dense tags are deleted by deleting the tag itself using tag_delete.
      \param tag_handle Handle of the (sparse) tag being deleted from entity
      \param entity_range Range of entities from which the tag is being deleted
  */
  virtual MBErrorCode  tag_delete_data(const MBTag tag_handle, 
                                        const MBRange &entity_range);

  //! Removes the tag from the database and deletes all of its associated data.
  virtual MBErrorCode  tag_delete(MBTag tag_handle);

  /**a;dlfa;sfsdafasdfl; 
     a;dlfja;sljfl;sadfasd
     a;dlkfj;lsajdf */

  //! creates a mesh set
  virtual MBErrorCode create_meshset(const unsigned int options, 
                                     MBEntityHandle &ms_handle,
                                     int start_id = 0,
                                     int start_proc = -1);

  //! Empty a vector of mesh set
  /** Empty a mesh set.
      \param ms_handles 1d vector of handles of sets being emptied
      \param num_meshsets Number of entities in 1d vector
  */
  virtual MBErrorCode clear_meshset(MBEntityHandle *ms_handles, const int num_meshsets);

  //! Empty a range of mesh set
  /** Empty a mesh set.
      \param ms_handles Range of handles of sets being emptied
  */
  virtual MBErrorCode clear_meshset(MBRange &ms_handles);

  //! get the options of a mesh set
  virtual MBErrorCode get_meshset_options(const MBEntityHandle ms_handle, 
                                           unsigned int& options) const;

  //! subtracts meshset2 from meshset1 - modifies meshset1
  virtual MBErrorCode subtract_meshset(MBEntityHandle meshset1, 
                                        const MBEntityHandle meshset2);

  //! intersects meshset2 with meshset1 - modifies meshset1
  virtual MBErrorCode intersect_meshset(MBEntityHandle meshset1, 
                                         const MBEntityHandle meshset2);
    
  //! unites meshset2 with meshset1 - modifies meshset1
  virtual MBErrorCode unite_meshset(MBEntityHandle meshset1, 
                                     const MBEntityHandle meshset2);

  //! add entities to meshset
  virtual MBErrorCode add_entities(MBEntityHandle meshset, 
                                    const MBRange &entities);

  //! add entities to meshset
  virtual MBErrorCode add_entities(MBEntityHandle meshset, 
                                    const MBEntityHandle *entities,
                                    const int num_entities);
  
  //! remove entities from meshset
  virtual MBErrorCode remove_entities(MBEntityHandle meshset, 
                                       const MBRange &entities);

  //! remove entities from meshset
  virtual MBErrorCode remove_entities(MBEntityHandle meshset, 
                                       const MBEntityHandle *entities,
                                       const int num_entities);

  //------MeshSet Parent/Child functions------
  
  //! get parent meshsets
  virtual MBErrorCode get_parent_meshsets(const MBEntityHandle meshset,
                                           std::vector<MBEntityHandle> &parents, 
                                           const int num_hops = 1) const;

  //! get parent meshsets
  virtual MBErrorCode get_parent_meshsets(const MBEntityHandle meshset,
                                          MBRange &parents, 
                                          const int num_hops = 1) const;

  //! get child meshsets
  virtual MBErrorCode get_child_meshsets(const MBEntityHandle meshset, 
                                          std::vector<MBEntityHandle> &children, 
                                          const int num_hops = 1) const;

  //! get child meshsets
  virtual MBErrorCode get_child_meshsets(const MBEntityHandle meshset, 
                                         MBRange &children, 
                                          const int num_hops = 1) const;

  //! gets number of parent meshsets
  virtual MBErrorCode num_parent_meshsets(const MBEntityHandle meshset,  
                                          int *number,
                                          const int num_hops = 1) const;

  //! gets number of child meshsets
  virtual MBErrorCode num_child_meshsets(const MBEntityHandle meshset, 
                                         int *number, 
                                         const int num_hops = 1) const;

  //! add a parent meshset
  virtual MBErrorCode add_parent_meshset(MBEntityHandle meshset, 
                                          const MBEntityHandle parent_meshset);

  //! add a child meshset
  virtual MBErrorCode add_child_meshset(MBEntityHandle meshset, 
                                         const MBEntityHandle child_meshset);

  //! adds 'parent' to child's parent list and adds 'child' to parent's child list
  virtual MBErrorCode add_parent_child( MBEntityHandle parent, 
                                         MBEntityHandle child );

  //! removes 'parent' to child's parent list and removes 'child' to parent's child list
  virtual MBErrorCode remove_parent_child( MBEntityHandle parent, 
                                            MBEntityHandle child );

  //! remove parent meshset
  virtual MBErrorCode remove_parent_meshset(MBEntityHandle meshset, 
                                             const MBEntityHandle parent_meshset);
  
  //! remove child meshset
  virtual MBErrorCode remove_child_meshset(MBEntityHandle meshset, 
                                            const MBEntityHandle child_meshset);

  // ************************  error condition information *************** 

    //! return various specific tag handles
  MBTag material_tag();
  MBTag neumannBC_tag();
  MBTag dirichletBC_tag();
  MBTag globalId_tag();
  MBTag geom_dimension_tag();

    //! get/set the number of nodes
    //int total_num_nodes() const;
    //void total_num_nodes(const int val);
  
    //! get/set the number of elements
    //int total_num_elements() const;
    //void total_num_elements(const int val);

    //! return a reference to the tag server
  TagServer* tag_server() {return tagServer;}

    //! return a reference to the sequence manager
  EntitySequenceManager* sequence_manager() { return sequenceManager; }
  const EntitySequenceManager* sequence_manager() const { return sequenceManager; }

    //! return the a_entity_factory pointer
  AEntityFactory *a_entity_factory() { return aEntityFactory; }
  
    //! return set of registered IO tools
  MBReaderWriterSet* reader_writer_set() { return readerWriterSet; }


//-----------------MeshSet Interface Functions------------------//

  void print(const MBEntityHandle handle, const char *prefix,
             bool first_call = true) const;

  virtual MBErrorCode get_last_error(std::string& info) const;

  virtual std::string get_error_string(const MBErrorCode code) const;

    //! check all adjacencies for consistency
  MBErrorCode check_adjacencies();
  
    //! check some adjacencies for consistency
  MBErrorCode check_adjacencies(const MBEntityHandle *ents, int num_ents);
  
    //! return whether the input handle is valid or not
  bool is_valid(const MBEntityHandle this_ent);
  
  
  
  const MBProcConfig& proc_config() const 
    { return procInfo; }

private:

  const MBProcConfig procInfo;

    //! database init and de-init routines
  MBErrorCode initialize();
  void deinitialize();

    //! return the entity set representing the whole mesh
  MBEntityHandle get_root_set();
  
    // other interfaces for MB
  MBWriteUtil* mMBWriteUtil;
  MBReadUtil* mMBReadUtil;

    //! store the total number of elements defined in this interface
    //int totalNumElements;
  
    //! store the total number of nodes defined in this interface
    //int totalNumNodes;

    //! the overall geometric dimension of this mesh
  int geometricDimension;

  MBTag materialTag;
  MBTag neumannBCTag;
  MBTag dirichletBCTag;
  MBTag geomDimensionTag;
  MBTag globalIdTag;

    //! tag server for this interface
  TagServer* tagServer;

  EntitySequenceManager *sequenceManager;

  AEntityFactory *aEntityFactory;
  
  MBReaderWriterSet* readerWriterSet;

    //! a meshset for the overall mesh; used primarily to set tags on entire mesh
  MBEntityHandle myMeshSet;

  MBError* mError;

  static const char *errorStrings[];
  
//------------MeshSet Interface Private Functions & Data------------//

  MeshSetManager* mMeshSetManager;
};

inline float MBCore::impl_version(std::string *version_string) 
{
  if (NULL != version_string)
    *version_string = std::string("MOAB IMPLEMENTATION version ") + std::string(MOAB_IMPL_VERSION_STRING);
  return MOAB_IMPL_VERSION;
}
  
#endif   // MB_IMPL_GENERAL_HPP
