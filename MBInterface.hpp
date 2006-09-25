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

/**
 * \class MBInterface
 * \brief Main interface class to MOAB
 * \author Tim Tautges, Karl Merkley, Ray Meyers, Corey Ernst, Clinton Stimpson,
 * \author Hong-Jun Kim, Jason Kraftcheck
 * \version 1.00
 * \date April, 2004
 */
#ifndef MB_INTERFACE_HPP
#define MB_INTERFACE_HPP

#define MOAB_API_VERSION 1.01
#define MOAB_API_VERSION_STRING "1.01"

//! include files
#include <vector>
#include <string>
#include <functional>

#include "MBEntityType.h"

//! forward declarations
class MBRange;
class MBInterface;


//! component architecture definitions
#ifdef XPCOM_MB

#ifndef __gen_nsISupports_h__
#include "nsISupports.h"
#endif

#ifndef NS_NO_VTABLE
#define NS_NO_VTABLE
#endif

#define MBINTERFACE_IID_STR "f728830e-1dd1-11b2-9598-fb9f414f2465"

{0xf728830e, 0x1dd1, 0x11b2, \
  { 0x95, 0x98, 0xfb, 0x9f, 0x4
#define MBINTERFACE_IID \1, 0x4f, 0x24, 0x65 }}

#endif


#include "MBUnknownInterface.h"
#define MB_INTERFACE_VERSION "2.0.0"
static const MBuuid IDD_MBCore = MBuuid( 0x8956e0a, 0xc300, 0x4005,
                                         0xbd, 0xf6, 0xc3, 0x4e, 0xf7, 0x1f, 0x5a, 0x52 );


#ifdef WIN32
#ifdef MB_EXPORTS
#define MB_DLL_EXPORT __declspec(dllexport)
#else
#define MB_DLL_EXPORT
#endif
#else
#define MB_DLL_EXPORT
#endif

//! \name Types and names
//! Definitions of types used in the interface

  //@{

//! first, specific error codes returned from MB
enum MBErrorCode { MB_SUCCESS = 0,
                   MB_INDEX_OUT_OF_RANGE,
                   MB_TYPE_OUT_OF_RANGE,
                   MB_MEMORY_ALLOCATION_FAILED,
                   MB_ENTITY_NOT_FOUND,
                   MB_MULTIPLE_ENTITIES_FOUND,
                   MB_TAG_NOT_FOUND,
                   MB_FILE_DOES_NOT_EXIST,
                   MB_FILE_WRITE_ERROR,
                   MB_NOT_IMPLEMENTED,
                   MB_ALREADY_ALLOCATED,
                   MB_FAILURE};

//! tag types: for a more detailed description, see the MB User's Guide
//! MB_TAG_BIT: size measured in bits instead of bytes, otherwise identical to sparse
//! MB_TAG_SPARSE: tags stored in (entity handle, tag value) pairs
//! MB_TAG_DENSE: tags stored in vectors directly on entity sequences, cheaper for tags
//!     which go on lots of entities
enum MBTagType {
  MB_TAG_BIT = 0, 
  MB_TAG_SPARSE, 
  MB_TAG_DENSE, 
  MB_TAG_MESH, 
  MB_TAG_LAST};

enum MBDataType {
  MB_TYPE_OPAQUE  = 0,
  MB_TYPE_INTEGER = 1,
  MB_TYPE_DOUBLE  = 2,
  MB_TYPE_BIT     = 3,
  MB_TYPE_HANDLE  = 4 };

//! entity handle: data type used to access all entity types in MB; high-order 4 bits
//! stores the entity type (defined in Mesh), the rest is id
#ifdef USE_64_BIT_HANDLES
  // use 64 bit entity handles on a 64 bit machine
typedef unsigned long MBEntityHandle;
#else
  // use 32 bit integer entity handles on 32/64 bit machines.
typedef unsigned int MBEntityHandle;
#endif
 
//! Tag handle: used to reference tags; since they're so different from entities, we
//! use void** instead of a uint to prevent them from being confused as entity handles.
typedef void** MBTag;

//! Meshset options: used to pass property options through the interface
//! MESHSET_TRACK_OWNER: enable entity to meshset adjacencies
//! MESHSET_SET: stores an MBRange of handles
//! MESHSET_ORDERED: stores a vector of handles
#define MESHSET_TRACK_OWNER  0x1
#define MESHSET_SET          0x2
#define MESHSET_ORDERED      0x4

//! convenience items: not critical to the representation, but useful enough to appear here

//! global pointer for easy access to a single MB database, declared here but defined by
//! the application.  Applications are not required to define, since this variable is not used
//! inside MB.  Note that users of multiple simultaneous MB databases will have to manage 
//! their own MB instance pointers.
extern MBInterface* gMB;

//! typedef for handle vectors
typedef std::vector<MBEntityHandle> MBHandleVec;

  //@}

#if defined(XPCOM_MB)
class NS_NO_VTABLE MBInterface : public nsISupports {
#else
class MB_DLL_EXPORT MBInterface : public MBUnknownInterface {
#endif

public:

#ifdef XPCOM_MB
  NS_DEFINE_STATIC_IID_ACCESSOR(MBINTERFACE_IID)
#endif

      //! \name Interface-level functions

      //@{

      //! constructor
  MBInterface() {}

    //! destructor
  virtual ~MBInterface() {}

    //! return the entity set representing the whole mesh
  virtual MBEntityHandle get_root_set()=0;
  
    //! query an MB internal interface
  virtual MBErrorCode query_interface(const std::string& iface_name, void** iface)=0;
 
    //! release an MB internal interface 
  virtual MBErrorCode release_interface(const std::string& iface_name, void* iface)=0;

    //! Returns the major.minor version number of the interface
    /**
       \param version_string If non-NULL, will be filled in with a string, possibly 
       containing implementation-specific information
    */
  virtual float api_version(std::string *version_string = NULL);
    
    //! Returns the major.minor version number of the implementation
    /**
       \param version_string If non-NULL, will be filled in with a string, possibly 
       containing implementation-specific information
    */
  virtual float impl_version(std::string *version_string = NULL)=0;

    //@}    

    //! \name Type and id utility functions

    //@{

    //! Returns the entity type of an MBEntityHandle.
    /** Returns the MBEntityType (ie, MeshVertex, MeshQuad, MeshHex ) of <em>handle</em>.
        \param handle The MBEntityHandle you want to find the entity type of.
        \return type The entity type of <em>handle</em>. 

        Example: \code
        MBEntityType type = type_from_handle( handle); 
        if( type == MeshHex ) ...  \endcode 
    */
  virtual MBEntityType type_from_handle(const MBEntityHandle handle) const = 0;
 
    //! Returns the id from an MBEntityHandle.
    /** \param handle The MBEntityHandle you want to find the id of. 
        \return id Id of <em>handle</em>
     
        Example: \code
        int id = id_from_handle(handle); \endcode 
    */
  virtual unsigned int id_from_handle(const MBEntityHandle handle) const =0;

    //! Returns the topological dimension of an entity
    /** Returns the MBEntityType (ie, MeshVertex, MeshQuad, MeshHex ) of <em>handle</em>.
        \param handle The MBEntityHandle you want to find the dimension of.
        \return type The topological dimension of <em>handle</em>. 

        Example: \code
        int dim = dimension_from_handle( handle); 
        if( dim == 0 ) ...  \endcode 
    */
  virtual int dimension_from_handle(const MBEntityHandle handle) const = 0;

    //! Gets an entity handle from the data base, if it exists, according to type and id.
    /** Given an MBEntiyType and an id, this function gets the existent MBEntityHandle. 
        If no such MBEntityHandle exits, it returns MB_ENTITY_NOT_FOUND 
        and sets handle to zero.
        \param type The type of the MBEntityHandle to retrieve from the database.
        \param id The id of the MBEntityHandle to retrieve from the database.
        \param handle An MBEntityHandle of type <em>type</em> and <em>id</em>. 

        Example: \code
        MBEntityType handle;
        MBErrorCode error_code = handle_from_id(MeshTri, 204, handle );
        if( error_code == MB_ENTITY_NOT_FOUND ) ... \endcode
    */
  virtual MBErrorCode handle_from_id(const MBEntityType type, const unsigned int id, 
                                     MBEntityHandle& handle) const =0;

    //@}

    //! \name Mesh input/output functions

    //@{

    //! Loads a mesh file into the database.
    /** Loads the file 'file_name'; types of mesh which can be loaded depend on modules available
        at MB compile time.  If active_block_id_list is NULL, all material sets (blocks in the 
        ExodusII jargon) are loaded.  Individual material sets  can be loaded by specifying their 
        ids in 'active_block_id_list'.  All nodes are loaded on first call for a given file.  
        Subsequent calls for a file load any material sets not loaded in previous calls.
        \param file_name Name of file to load into database.
        \param active_block_id_list Material set/block ids to load.  If NULL, ALL blocks of 
        <em>file_name</em> are loaded.
        \param num_blocks Number of blocks in active_block_id_list

        Example: \code
        std::vector<int> active_block_id_list;
        int active_block_id_list[] = {1, 4, 10};
        load_mesh( "temp.gen", active_block_id_list, 3 );  //load blocks 1, 4, 10 \endcode 
    */
  virtual MBErrorCode load_mesh(const char *file_name,
                                const int *active_block_id_list = NULL,
                                const int num_blocks = 0)=0;

    //! Writes mesh to a file.
    /** Write mesh to file 'file_name'; if output_list is non-NULL, only material sets contained
        in that list will be written.
        \param file_name Name of file to write.
        \param output_list 1d array of material set handles to write; if NULL, all sets are written
        \param num_sets Number of sets in output_list array

        Example: \code
        MBEntityHandle output_list[] = {meshset1, meshset2, meshset3}; 
        write_mesh( "output_file.gen", output_list, 3 ); \endcode 
    */
  virtual MBErrorCode write_mesh(const char *file_name,
                                 const MBEntityHandle *output_list = NULL,
                                 const int num_sets = 0) = 0;

    //! Deletes all mesh entities from this MB instance
  virtual MBErrorCode delete_mesh()=0;

    //@}

    //! \name Geometric dimension functions

    //@{

    //! Get overall geometric dimension
  virtual MBErrorCode get_dimension(int &dim) const =0;

    //! Set overall geometric dimension
    /** Returns error if setting to 3 dimensions, mesh has been created, and 
     *  there are only 2 dimensions on that mesh
     */
  virtual MBErrorCode set_dimension(const int dim)=0;

    //@}

    //! \name Vertex coordinate functions

    //@{

    //! Get blocked vertex coordinates for all vertices
    /** Blocked = all x, then all y, etc. 
          
    Example: \code
    std::vector<double> coords;
    get_vertex_coordinates(coords);
    double xavg = 0;
    for (int i = 0; i < coords.size()/3; i++) xavg += coords[i]; \endcode
    */
  virtual MBErrorCode get_vertex_coordinates(std::vector<double> &coords) const =0;

    //! Gets xyz coordinate information for range of vertices
    /** Length of 'coords' should be at least 3*<em>entity_handles.size()</em> before making call.
        \param entity_handles Range of vertex handles (error if not of type MeshVertex)
        \param coords Array used to return x, y, and z coordinates.
   
        Example: \code 
        double coords[3];
        get_coords( vertex_handle, coords ); 
        std::cout<<"x = "<<coords[0]<<std::endl;
        std::cout<<"y = "<<coords[1]<<std::endl;
        std::cout<<"z = "<<coords[2]<<std::endl; \endcode 
    */
  virtual MBErrorCode  get_coords(const MBRange& entity_handles, 
                                  double *coords) const =0;
    
    //! Gets xyz coordinate information for vector of vertices
    /** Identical to range-based function, except entity handles are specified using a 1d vector
        and vector length.
    */
  virtual MBErrorCode  get_coords(const MBEntityHandle* entity_handles, 
                                  const int num_entities, 
                                  double *coords) const =0;
  
    //! Sets the xyz coordinates for a vector of vertices
    /** An error is returned if any entities in the vector are not vertices.
        \param entity_handles MBEntityHandle's to set coordinates of. (Must be of type MeshVertex)
        \param num_entities Number of entities in entity_handles
        \param coords Array containing new xyz coordinates.
 
        Example: \code
        double coords[3] = {0.234, -2.52, 12.023};
        set_coords( entity_handle, 1, coords ); \endcode 
    */
  virtual MBErrorCode  set_coords(MBEntityHandle *entity_handles, 
                                  const int num_entities,
                                  const double *coords)=0;

    //@}

    //! \name Connectivity functions

    //@{

    //! Get the connectivity array for all entities of the specified entity type
    /**  This function returns the connectivity of just the corner vertices, no higher order nodes
         \param type The entity type of elements whose connectivity is to be returned
         \param connect an STL vector used to return connectivity array (in the form of entity handles)
    */
  virtual MBErrorCode get_connectivity_by_type(const MBEntityType type, 
                                               std::vector<MBEntityHandle> &connect) const =0;

    //! Gets the connectivity for a vector of elements
    /** Same as vector-based version except range is returned (unordered!)
    */
  virtual MBErrorCode  get_connectivity(const MBEntityHandle *entity_handles, 
                                        const int num_handles,
                                        MBRange &connectivity, 
                                        bool topological_connectivity = false) const =0;
 
    //! Gets the connectivity for a vector of elements
    /** Corner vertices or all vertices (including higher-order nodes, if any) are returned.
        For non-element handles (ie, MB_MeshSets), returns an error. Connectivity data is copied 
        from the database into the vector.  Connectivity of a vertex is the same vertex.
        The nodes in <em>connectivity</em> are properly ordered according to that element's 
        canonical ordering.
        \param entity_handles Vector of element handles to get connectivity of.
        \param num_handles Number of entity handles in <em>entity_handles</em>
        \param connectivity Vector in which connectivity of <em>entity_handles</em> is returned.  
        \param topological_connectivity If true, higher order nodes are ignored. 
    */
  virtual MBErrorCode  get_connectivity(const MBEntityHandle *entity_handles, 
                                        const int num_handles,
                                        std::vector<MBEntityHandle> &connectivity, 
                                        bool topological_connectivity = false) const =0;
 
    //! Gets a pointer to constant connectivity data of <em>entity_handle</em> 
    /** Sets <em>number_nodes</em> equal to the number of nodes of the <em> 
        entity_handle </em>.  Faster then the other <em>get_connectivity</em> function because no
        data is copied.  The nodes in 'connectivity' are properly ordered according to the 
        element's canonical ordering.
        \param entity_handle MBEntityHandle to get connectivity of.
        \param connectivity Array in which connectivity of <em>entity_handle</em> is returned.
        \param num_nodes Number of MeshVertices in array <em>connectivity</em>. 
        \param topological_connectivity If true, num_nodes will be set to number of corner vertices
        for that element type.
    */
  virtual MBErrorCode  get_connectivity(const MBEntityHandle entity_handle, 
                                        const MBEntityHandle *&connectivity, 
                                        int &num_nodes, 
                                        bool topological_connectivity = false) const =0;

    //! Sets the connectivity for an MBEntityHandle.  For non-element handles, return an error.
    /** Connectivity is stored exactly as it is ordered in vector <em>connectivity</em>. 
        \param entity_handle MBEntityHandle to set connectivity of.
        \param connect Vector containing new connectivity of <em>entity_handle</em>.
        \param num_connect Number of vertices in <em>connect</em>
   
        Example: \code 
        MBEntityHandle conn[] = {node1, node2, node3};
        set_connectivity( tri_element, conn, 3 ); \endcode 
    */
  virtual MBErrorCode  set_connectivity(const MBEntityHandle entity_handle, 
                                        MBEntityHandle *connect,
                                        const int num_connect)=0;

    //! Sets the connectivity for an MBEntityHandle.  For non-element handles, return an error.
    //@}

    //! \name Adjacencies functions 

    //@{

    //! Get the adjacencies associated with a vector of entities to entities of a specfied dimension.
    /** \param from_entities Vector of MBEntityHandle to get adjacencies of.
        \param num_entities Number of entities in <em>from_entities</em>
        \param to_dimension Dimension of desired adjacencies
        \param create_if_missing If true, MB will create any entities of the specfied dimension
        which have not yet been created (only useful when <em>to_dimension < dim(*from_entities)</em>)
        \param adj_entities STL vector in which adjacent entities are returned. 
        \param operation_type Enum of INTERSECT or UNION.  Defines whether to take
        the intersection or union of the set of adjacencies recovered for the from_entities.

        The adjacent entities in vector <em>adjacencies</em> are not in any particular 
        order. 

        Example: \code
        std::vector<MBEntityHandle> adjacencies, from_entities = {hex1, hex2};
          // generate all edges for these two hexes
          get_adjacencies( from_entities, 2, 1, true, adjacencies, MBInterface::UNION); 
          adjacencies.clear();
            // now find the edges common to both hexes
            get_adjacencies( from_entities, 2, 1, false, adjacencies, MBInterface::INTERSECT); 
            \endcode 
    */
  virtual MBErrorCode get_adjacencies(const MBEntityHandle *from_entities,
                                      const int num_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      std::vector<MBEntityHandle>& adj_entities,
                                      const int operation_type = MBInterface::INTERSECT) = 0;

    //! Get the adjacencies associated with a vector of entities to entities of a specfied dimension.
    /** Identical to vector-based get_adjacencies function, except results are returned in a
        range instead of a vector.
    */
  virtual MBErrorCode get_adjacencies(const MBEntityHandle *from_entities,
                                      const int num_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      MBRange &adj_entities,
                                      const int operation_type = MBInterface::INTERSECT) = 0;

    //! Get the adjacencies associated with a range of entities to entities of a specfied dimension.
    /** Identical to vector-based get_adjacencies function, except "from" entities specified in a
        range instead of a vector.
    */
  virtual MBErrorCode get_adjacencies(const MBRange &from_entities,
                                      const int to_dimension,
                                      const bool create_if_missing,
                                      MBRange &adj_entities,
                                      const int operation_type = MBInterface::INTERSECT) = 0;

    //! Adds adjacencies between "from" and "to" entities.
    /** \param from_handle Entities on which the adjacencies are placed
        \param to_handles Vector of entities referenced by new adjacencies added to <em>from_handle</em>
        \param num_handles Number of entities in <em>to_handles</em>
        \param both_ways If true, add the adjacency information in both directions; if false,
        adjacencies are added only to <em>from_handle</em>
    */
  virtual MBErrorCode add_adjacencies(const MBEntityHandle from_handle, 
                                      const MBEntityHandle *to_handles,
                                      const int num_handles,
                                      bool both_ways) = 0;

    //! Adds adjacencies; same as vector-based, but with range instead
  virtual MBErrorCode add_adjacencies(const MBEntityHandle from_handle, 
                                      MBRange &adjacencies,
                                      bool both_ways) = 0;

    //! Removes adjacencies between handles.
    /** Adjacencies in both directions are removed.
        \param from_handle Entity from which adjacencies are being removed.
        \param to_handles Entities to which adjacencies are being removed.
        \param num_handles Number of handles in <em>to_handles</em>
    */
  virtual MBErrorCode remove_adjacencies(const MBEntityHandle from_handle, 
                                         const MBEntityHandle *to_handles,
                                         const int num_handles) = 0;

    //@}

    //! Enumerated type used in get_adjacencies() and other functions
  enum {INTERSECT, UNION};

    //! \name Functions for getting entities

    //@{

    //! Retrieves all entities of a given topological dimension in the database or meshset.
    /** \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param dimension Topological dimension of entities desired.
        \param entities Range in which entities of dimension <em>dimension</em> are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.

        Example: \code
          // get 1d (edge) elements in the entire mesh
          MBRange edges;
          get_entities_by_dimension( 0, 1, edges );
          \endcode 
    */
  virtual MBErrorCode get_entities_by_dimension(const MBEntityHandle meshset,
                                                const int dimension, 
                                                MBRange &entities,
                                                const bool recursive = false)  const = 0;

    //! Retrieve all entities of a given type in the database or meshset.
    /** \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param entities Range in which entities of type <em>type</em> are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBMESHSET is an error, as it would always 
                         result in an empty list.

        Example: \code
          // get the quadrilateral elements in meshset
          MBRange quads;
          get_entities_by_type( meshset, MeshQuad, quads );
          \endcode 
    */
  virtual MBErrorCode get_entities_by_type(const MBEntityHandle meshset,
                                           const MBEntityType type, 
                                           MBRange &entities,
                                           const bool recursive = false) const = 0;

    //! Retrieve entities in the database or meshset which have any or all of the tag(s) and (optionally)
    //! value(s) specified.
    /** \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param tag_handles Vector of tag handles entities must have
        \param values Vector of pointers to values of tags in <em>tag_handles</em>
        \param num_tags Number of tags and values in <em>tag_handles</em> and <em>values</em>
        \param entities Range in which entities are returned.
        \param condition Boolean condition, either MBInterface::UNION or MBInterface::INTERSECT
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBMESHSET is an error, as it would always 
                         result in an empty list.

        If MBInterface::UNION is specified as the condition, entities with <em>any</em> of the tags
        and values specified are returned.  If MBInterface::INTERSECT is specified, only entities with
        <em>all</em> of the tags/values are returned.

        If <em>values</em> is NULL, entities with the specified tags and any corresponding values are
        returned.  Note that if <em>values</em> is non-NULL, it is a vector of <em>pointers</em> to
        tag values.

        Example: \code
          // get the dirichlet sets in a mesh
          MBRange dir_sets;
          MBTag dir_tag;
          tag_get_handle(DIRICHLET_SET_TAG_NAME, dir_tag);
          get_entities_by_type_and_tag(0, MeshEntitySet, &dir_tag, NULL, 1, dir_sets, 
          MBInterface::UNION);
          \endcode 
    */
  virtual MBErrorCode get_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                   const MBEntityType type,
                                                   const MBTag *tag_handles,
                                                   const void* const* values,
                                                   const int num_tags,
                                                   MBRange &entities,
                                                   const int condition = MBInterface::INTERSECT,
                                                   const bool recursive = false) const = 0;

    //! Returns all entities in the data base or meshset, in a range (order not preserved)
    /** \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param entities Range in which entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.

        Example: \code
        MBRange entities;
          // get all non-meshset entities in meshset, including in contained meshsets
          get_entities_by_handle(meshset, entities, true);
          \endcode 
    */
  virtual MBErrorCode get_entities_by_handle(const MBEntityHandle meshset,
                                             MBRange &entities,
                                             const bool recursive = false) const = 0;

    //! Returns all entities in the data base or meshset, in a vector (order preserved)
    /** \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param entities STL vector in which entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.

        Example: \code
        std::vector<MBEntityHandle> entities;
          // get all non-meshset entities in meshset, including in contained meshsets
          get_entities_by_handle(meshset, entities, true);
          \endcode 
    */
  virtual MBErrorCode get_entities_by_handle(const MBEntityHandle meshset,
                                             std::vector<MBEntityHandle> &entities,
                                             const bool recursive = false) const = 0;

    //! Return the number of entities of given dimension in the database or meshset
    /** \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param dimension Dimension of entities desired.
        \param num_entities Number of entities of the given dimension
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual MBErrorCode get_number_entities_by_dimension(const MBEntityHandle meshset,
                                                       const int dimension, 
                                                       int &num_entities,
                                                       const bool recursive = false) const = 0;

    //! Retrieve the number of entities of a given type in the database or meshset.
    /** Identical to get_entities_by_dimension, except returns number instead of entities
        \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param num_entities Number of entities of type <em>type</em>
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBMESHSET is an error, as it would always 
                         result in an empty list.
    */
  virtual MBErrorCode get_number_entities_by_type(const MBEntityHandle meshset,
                                                  const MBEntityType type, 
                                                  int &num_entities,
                                                  const bool recursive = false) const = 0;

    //! Retrieve number of entities in the database or meshset which have any or all of the 
    //! tag(s) and (optionally) value(s) specified.
    /** Identical to get_entities_by_type_and_tag, except number instead of entities are returned
        \param meshset Meshset whose entities are being queried (zero if query is for entire mesh).
        \param type Type of entities to be returned
        \param tag_handles Vector of tag handles entities must have
        \param values Vector of pointers to values of tags in <em>tag_handles</em>
        \param num_tags Number of tags and values in <em>tag_handles</em> and <em>values</em>
        \param num_entities Range in which number of entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves.  Specifying 
                         both recursive=true and type=MBMESHSET is an error, as it would always 
                         result in an empty list.
    */
  virtual MBErrorCode get_number_entities_by_type_and_tag(const MBEntityHandle meshset,
                                                          const MBEntityType type,
                                                          const MBTag *tag_handles,
                                                          const void** values,
                                                          const int num_tags,
                                                          int &num_entities,
                                                          const bool recursive = false) const = 0;

    //! Returns number of entities in the data base or meshset
    /** Identical to get-entities_by_handle, except number instead of entities are returned
        \param meshset Meshset whose entities are being queried (zero if query is for the entire mesh).
        \param num_entities Range in which num_entities are returned.
        \param recursive If true, meshsets containing meshsets are queried recursively.  Returns
                         the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual MBErrorCode get_number_entities_by_handle(const MBEntityHandle meshset,
                                                    int &num_entities,
                                                    const bool recursive = false) const = 0;

    //@}

    //! \name Modifying the mesh

    //@{

    //! Create an element based on the type and connectivity. 
    /** Create a new element in the database.  Vertices composing this element must already exist,
        and connectivity must be specified in canonical order for the given element type.  If 
        connectivity vector is not correct for MBEntityType <em>type</em> (ie, a vector with 
        3 vertices is passed in to make an MeshQuad), the function returns MB_FAILURE. 
        \param type Type of element to create. (MeshTet, MeshTri, MeshKnife, etc.) 
        \param connectivity 1d vector containing connectivity of element to create.
        \param num_vertices Number of vertices in element
        \param element_handle Handle representing the newly created element in the database.

        Example: \code
        MBEntityHandle quad_conn[] = {vertex0, vertex1, vertex2, vertex3};
        MBEntityHandle quad_handle = 0;
        create_element( MeshQuad, quad_conn, 4, new_handle ); \endcode 
    */
  virtual MBErrorCode create_element(const MBEntityType type, 
                                     const MBEntityHandle *connectivity,
                                     const int num_vertices, 
                                     MBEntityHandle &element_handle) = 0;

    //! Creates a vertex with the specified coordinates.  
    /**
       \param coordinates Array that has 3 doubles in it.
       \param entity_handle MBEntityHandle representing the newly created vertex in the database.

       Example: \code
       double coordinates[] = {1.034, 23.23, -0.432};
       MBEntityHandle new_handle = 0;
       create_vertex( coordinates, entity_handle ); \endcode 
    */
  virtual MBErrorCode create_vertex(const double coordinates[3], 
                                    MBEntityHandle &entity_handle ) = 0;

    //! Merge two entities into a single entity
    /** Merge two entities into a single entities, with <em>entity_to_keep</em> receiving
        adjacencies that were on <em>entity_to_remove</em>.
        \param entity_to_keep Entity to be kept after merge
        \param entity_to_remove Entity to be merged into <em>entity_to_keep</em>
        \param auto_merge If false, <em>entity_to_keep</em> and <em>entity_to_remove</em> must share
        the same lower-dimensional entities; if true, MB tries to merge those entities automatically
        \param delete_removed_entity If true, <em>entity_to_remove</em> is deleted after merge is complete
    */
  virtual MBErrorCode merge_entities(MBEntityHandle entity_to_keep, 
                                     MBEntityHandle entity_to_remove,
                                     bool auto_merge,
                                     bool delete_removed_entity) = 0;

    //! Removes entities in a vector from the data base.  
    /** If any of the entities are contained in any meshsets, it is removed from those meshsets 
        which were created with MESHSET_TRACK_OWNER option bit set.  Tags for <em>entity</em> are 
        removed as part of this function.
        \param entities 1d vector of entities to delete
        \param num_entities Number of entities in 1d vector
    */ 
  virtual MBErrorCode delete_entities(const MBEntityHandle *entities,
                                      const int num_entities) = 0;

    //! Removes entities in a range from the data base.  
    /** If any of the entities are contained in any meshsets, it is removed from those meshsets 
        which were created with MESHSET_TRACK_OWNER option bit set.  Tags for <em>entity</em> are 
        removed as part of this function.
        \param entities Range of entities to delete
    */ 
  virtual MBErrorCode delete_entities(const MBRange &entities) = 0;

    //@}

    //! \name Listing entities

    //@{

    //! List entities to standard output
    /** Lists all data pertaining to entities (i.e. vertex coordinates if vertices, connectivity if
        elements, set membership if set).  Useful for debugging, but output can become quite long
        for large databases.
    */
  virtual MBErrorCode list_entities(const MBRange &entities) const = 0;
  
    //! List entities, or number of entities in database, to standard output
    /** Lists data pertaining to entities to standard output.  If <em>entities</em> is NULL and
        <em>num_entities</em> is zero, lists only the number of entities of each type in the 
        database.  If <em>entities</em> is NULL and <em>num_entities</em> is non-zero, lists all
        information for all entities in the database.
        \param entities 1d vector of entities to list
        \param num_entities Number of entities in 1d vector
    */
  virtual MBErrorCode list_entities(const MBEntityHandle *entities,
                                    const int num_entities) const = 0;

    //! List a single entity; no header printed
    /** Lists a single entity, including its connectivity and its adjacencies.
     *  No header is printed, because calling function might print information between header
     *  and information printed by this function.
     *  \param entity The entity to be listed.
     */
  virtual MBErrorCode list_entity(const MBEntityHandle entity) const = 0;

    //@}

    //! \name Functions for higher-order elements

    //@{

    //! function object for recieving events from MB of higher order nodes added to entities
  class HONodeAddedRemoved
  {
  public:
      //! Constructor
    HONodeAddedRemoved(){}
 
      //! Destructor
    virtual ~HONodeAddedRemoved(){}

      //! node_added called when a node was added to an element's connectivity array
      //! note: connectivity array of element may be incomplete (corner nodes will exist always)
      /** 
       * \param node Node being added
       * \param element Element node is being added to
       */
    virtual void node_added(MBEntityHandle node, MBEntityHandle element) = 0;

      //! node_added called when a node was added to an element's connectivity array
      //! note: connectivity array of element may be incomplete (corner nodes will exist always)
      /**
       * \param node Node being removed.
       */
    virtual void node_removed(MBEntityHandle node) = 0;
  };
  
    //! Convert entities to higher-order elements by adding mid nodes
    /** This function causes MB to create mid-nodes on all edges, faces, and element interiors 
        for all entities in <em>meshset</em>.  Higher order nodes appear in an element's connectivity
        array according to the algorithm described in the documentation for Mesh.  If 
        <em>HONodeAddedRemoved</em> function is input, this function is called to notify the application
        of nodes being added/removed from the mesh.
        \param meshset The set of entities being converted
        \param mid_edge If true, mid-edge nodes are created 
        \param mid_face If true, mid-face nodes are created 
        \param mid_region If true, mid-element nodes are created 
        \param function_object If non-NULL, the node_added or node_removed functions on this object 
        are called when nodes are added or removed from an entity, respectively
    */
  virtual MBErrorCode convert_entities(const MBEntityHandle meshset, 
                                       const bool mid_edge,
                                       const bool mid_face, 
                                       const bool mid_region, 
                                       HONodeAddedRemoved* function_object = 0) = 0;

    //! Returns the side number, in canonical ordering, of <em>child</em> with respect to <em>parent</em>
    /** Given a parent and child entity, returns the canonical ordering information side number, sense, 
        and offset of <em>child</em> with respect to <em>parent</em>.  This function returns
        MB_FAILURE if <em>child</em> is not related to <em>parent</em>.  This function does *not* 
        create adjacencies between <em>parent</em> and <em>child</em>.
        \param parent Parent entity to be compared
        \param child Child entity to be compared
        \param side_number Side number in canonical ordering of <em>child</em> with respect to 
        <em>parent</em>
        \param sense Sense of <em>child</em> with respect to <em>parent</em>, assuming ordering of 
        <em>child</em> as given by get_connectivity called on <em>child</em>; sense is 1, -1
        for forward/reverse sense, resp.
        \param offset Offset between first vertex of <em>child</em> and first vertex of side 
        <em>side_number</em> on <em>parent</em>
    */
  virtual MBErrorCode side_number(const MBEntityHandle parent,
                                  const MBEntityHandle child,
                                  int &side_number,
                                  int &sense,
                                  int &offset) const = 0;

    //! Find the higher-order node on a subfacet of an entity
    /** Given an entity and the connectivity and type of one of its subfacets, find the
        high order node on that subfacet, if any.  The number of vertices in <em>subfacet_conn</em>
        is derived from <em>subfacet_type</em> and the canonical numbering for that type.
        \param parent_handle The element whose subfacet is being queried
        \param subfacet_conn The connectivity of the subfacet being queried
        \param subfacet_type The type of subfacet being queried
        \param high_order_node If the subfacet has a high-order node defined on <em>parent_handle</em>,
        the handle for that node.
    */
  virtual MBErrorCode high_order_node(const MBEntityHandle parent_handle,
                                      const MBEntityHandle *subfacet_conn,
                                      const MBEntityType subfacet_type,
                                      MBEntityHandle &high_order_node) const = 0;

    //! Return the handle of the side element of a given dimension and index
    /** Given a parent entity and a target dimension and side number, return the handle of
        the entity corresponding to that side.  If an entity has not been created to represent
        that side, one is not created by this function, and zero is returned in <em>target_entity</em>.
        \param source_entity The entity whose side is being queried.
        \param dim The topological dimension of the side being queried.
        \param side_number The canonical index of the side being queried.
        \param target_entity The handle of the entity representing this side, if any.
    */
  virtual MBErrorCode side_element(const MBEntityHandle source_entity,
                                   const int dim, 
                                   const int side_number,
                                   MBEntityHandle &target_entity) const = 0;

    //@}

    //! \name Tag functions

    //@{

    //! Create a tag with the specified name, type and length
    /** Create a "tag", used to store application-defined data on MB entities.  If MB_ALREADY_ALLOCATED
        is returned, a tag with this name has already been created.  Tags created with this function
        are assigned to entities using the tag_set_data function described below.
        \param tag_name Name of this tag
        \param tag_size Size of data to store on tag, in bytes (MB_TAG_DENSE, MB_TAG_SPARSE) or
        bits (MB_TAG_BITS).
        \param type Type of tag to create (MB_TAG_BIT, MB_TAG_SPARSE, MB_TAG_DENSE, MB_TAG_MESH)
        \param tag_handle Tag handle created
        \param default_value Default value tag data is set to when initially created

        Example: \code
        MBTag tag_handle;
        double value = 100.0;
          // create a dense tag with default value of 100
          tag_create( "my_tag", sizeof(double), MB_TAG_DENSE, tag_handle, &value );\endcode
    */
  virtual MBErrorCode tag_create(const char *tag_name,
                                 const int tag_size, 
                                 const MBTagType type,
                                 MBTag &tag_handle, 
                                 const void *default_value) = 0;

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
                                              bool use_existing = false ) = 0;

    //! Get the name of a tag corresponding to a handle
    /** \param tag_handle Tag you want the name of.  
        \param tag_name Name string for <em>tag_handle</em>. 
    */
  virtual MBErrorCode  tag_get_name(const MBTag tag_handle, 
                                    std::string& tag_name) const = 0;

    //! Gets the tag handle corresponding to a name
    /** If a tag of that name does not exist, returns MB_TAG_NOT_FOUND
        \param tag_name Name of the desired tag. 
        \param tag_handle Tag handle corresponding to <em>tag_name</em>
    */ 
  virtual MBErrorCode  tag_get_handle(const char *tag_name, 
                                      MBTag &tag_handle) const = 0;

    //! Get the size of the specified tag
    /** Get the size of the specified tag, in bytes (MB_TAG_SPARSE, MB_TAG_DENSE, MB_TAG_MESH) or
        bits (MB_TAG_BIT).
        \param tag Handle of the desired tag. 
        \param tag_size Size of the specified tag
    */ 
  virtual MBErrorCode tag_get_size(const MBTag tag, int &tag_size) const = 0;

    //! Get the type of the specified tag
    /** Get the type of the specified tag
        \param tag Handle of the desired tag. 
        \param tag_type Type of the specified tag
    */ 
  virtual MBErrorCode tag_get_type(const MBTag tag, MBTagType &tag_type) const = 0;

    /** \brief Get data type of tag.
     *
     * Get the type of the tag data.  The tag is data is assumed to
     * be a vector of this type.  If the tag data vetcor contains 
     * more than one value, then the tag size must be a multiple of
     * the size of this type.
     * \param tag  The tag 
     * \param type The type of the specified tag (output).
     */
   virtual MBErrorCode tag_get_data_type(const MBTag tag, MBDataType& type) const = 0;

    //! Get the default value of the specified tag
    /** Get the default value of the specified tag
        \param tag Handle of the desired tag. 
        \param def_value Pointer to memory where default value of the specified tag is written
        \return - MB_ENTITY_NOT_FOUND If no default value is set for tag.
                - MB_SUCCESS          If success.
                - MB_FAILURE          If <code>def_val</code> is NULL.
                - MB_TAG_NOT_FOUND    If <code>tag_handle</code> is invalid.
    */ 
  virtual MBErrorCode tag_get_default_value(const MBTag tag, void *def_val) const = 0;

    //! Get handles for all tags defined in the mesh instance
    /** Get handles for all tags defined on the mesh instance.
        \param tag_handles STL vector of all tags
    */
  virtual MBErrorCode tag_get_tags(std::vector<MBTag> &tag_handles) const = 0;

    //! Get handles for all tags defined on this entity
    /** Get handles for all tags defined on this entity; if zero, get all tags defined 
        on mesh instance
        \param entity Entity for which you want tags
        \param tag_handles STL vector of all tags defined on <em>entity</em>
    */
  virtual MBErrorCode tag_get_tags_on_entity(const MBEntityHandle entity,
                                             std::vector<MBTag> &tag_handles) const = 0;

    //! Get the value of the indicated tag on the specified entities in the specified vector
    /** Get the value of the indicated tag on the specified entities; <em>tag_data</em> must contain
        enough space (i.e. tag_size*num_entities bytes or bits) to hold all tag data.  MB does <em>not</em>
        check whether this space is available before writing to it.
        \param tag_handle Tag whose values are being queried.
        \param entity_handles 1d vector of entity handles whose tag values are being queried
        \param num_entities Number of entities in 1d vector of entity handles
        \param tag_data Pointer to memory into which tag data will be written
    */
  virtual MBErrorCode  tag_get_data(const MBTag tag_handle, 
                                    const MBEntityHandle* entity_handles, 
                                    const int num_entities, 
                                    void *tag_data) const = 0;

    //! Get the value of the indicated tag on the specified entities in the specified range
    /** Identical to previous function, except entities are specified using a range instead of a 1d vector.
        \param tag_handle Tag whose values are being queried.
        \param entity_handles Range of entity handles whose tag values are being queried
        \param tag_data Pointer to memory into which tag data will be written
    */
  virtual MBErrorCode  tag_get_data(const MBTag tag_handle, 
                                    const MBRange& entity_handles, 
                                    void *tag_data) const = 0;

    //! Set the value of the indicated tag on the specified entities in the specified vector
    /** Set the value of the indicated tag on the specified entities; <em>tag_data</em> contains the
        values, <em>one value per entity in <em>entity_handles</em></em>.
        \param tag_handle Tag whose values are being set
        \param entity_handles 1d vector of entity handles whose tag values are being set
        \param num_entities Number of entities in 1d vector of entity handles
        \param tag_data Pointer to memory holding tag values to be set, <em>one entry per entity handle</em>
    */
  virtual MBErrorCode  tag_set_data(const MBTag tag_handle, 
                                    const MBEntityHandle* entity_handles, 
                                    const int num_entities,
                                    const void *tag_data ) = 0;
  
    //! Set the value of the indicated tag on the specified entities in the specified range
    /** Identical to previous function, except entities are specified using a range instead of a 1d vector.
        \param tag_handle Tag whose values are being set
        \param entity_handles Range of entity handles whose tag values are being set
        \param tag_data Pointer to memory holding tag values to be set, <em>one entry per entity handle</em>
    */
  virtual MBErrorCode  tag_set_data(const MBTag tag_handle, 
                                    const MBRange& entity_handles,
                                    const void *tag_data ) = 0;

    //! Delete the data of a vector of entity handles and sparse tag
    /** Delete the data of a tag on a vector of entity handles.  Only sparse tag data are deleted with this
        function; dense tags are deleted by deleting the tag itself using tag_delete.
        \param tag_handle Handle of the (sparse) tag being deleted from entity
        \param entity_handles 1d vector of entity handles from which the tag is being deleted
        \param num_handles Number of entity handles in 1d vector
    */
  virtual MBErrorCode  tag_delete_data(const MBTag tag_handle, 
                                       const MBEntityHandle *entity_handles,
                                       const int num_handles) = 0;

    //! Delete the data of a range of entity handles and sparse tag
    /** Delete the data of a tag on a range of entity handles.  Only sparse tag data are deleted with this
        function; dense tags are deleted by deleting the tag itself using tag_delete.
        \param tag_handle Handle of the (sparse) tag being deleted from entity
        \param entity_range Range of entities from which the tag is being deleted
    */
  virtual MBErrorCode  tag_delete_data(const MBTag tag_handle, 
                                       const MBRange &entity_range) = 0;

    //! Remove a tag from the database and delete all of its associated data
    /** Deletes a tag and all associated data.
     */
  virtual MBErrorCode  tag_delete(MBTag tag_handle) = 0;

    //@}

    //! \name Meshset functions

    //@{

    //! Create a new mesh set
    /** Create a new mesh set.  Meshsets can store entities ordered or unordered. A set can include entities
        at most once (MESHSET_SET) or more than once.  Meshsets can optionally track its members using
        adjacencies (MESHSET_TRACK_OWNER); if set, entities are deleted from tracking meshsets before
        being deleted.  This adds data to mesh entities, which can be expensive.
        \param options Options bitmask for the new meshset, possible values defined above
        \param ms_handle Handle for the meshset created
    */
  virtual MBErrorCode create_meshset(const unsigned int options, 
                                     MBEntityHandle &ms_handle,
                                     int start_id = 0,
                                     int start_proc = -1) = 0;

    //! Empty a vector of mesh set
    /** Empty a mesh set.
        \param ms_handles 1d vector of handles of sets being emptied
        \param num_meshsets Number of entities in 1d vector
    */
  virtual MBErrorCode clear_meshset(MBEntityHandle *ms_handles, const int num_meshsets) = 0;

    //! Empty a range of mesh set
    /** Empty a mesh set.
        \param ms_handles Range of handles of sets being emptied
    */
  virtual MBErrorCode clear_meshset(MBRange &ms_handles) = 0;

    //! Get the options of a mesh set
    /** Get the options of a mesh set.
        \param ms_handle Handle for mesh set being queried
        \param options Bit mask in which mesh set options are returned
    */
  virtual MBErrorCode get_meshset_options(const MBEntityHandle ms_handle, 
                                          unsigned int& options) const = 0;

    //! Subtract meshsets
    /** Subtract <em>meshset2</em> from <em>meshset1</em>, placing the results in meshset1.
        \param meshset1 Mesh set being subtracted from, also used to pass back result
        \param meshset2 Mesh set being subtracted from <em>meshset1</em>
    */
  virtual MBErrorCode subtract_meshset(MBEntityHandle meshset1, 
                                       const MBEntityHandle meshset2) = 0;

    //! Intersect meshsets
    /** Intersect <em>meshset1</em> with <em>meshset2</em>, placing the results in meshset1.
        \param meshset1 Mesh set being intersected, also used to pass back result
        \param meshset2 Mesh set being intersected with <em>meshset1</em>
    */
  virtual MBErrorCode intersect_meshset(MBEntityHandle meshset1, 
                                        const MBEntityHandle meshset2) = 0;
    
    //! Unite meshsets
    /** Unite <em>meshset1</em> with <em>meshset2</em>, placing the results in meshset1.
        \param meshset1 Mesh set being united, also used to pass back result
        \param meshset2 Mesh set being united with <em>meshset1</em>
    */
  virtual MBErrorCode unite_meshset(MBEntityHandle meshset1, 
                                    const MBEntityHandle meshset2) = 0;

    //! Add to a meshset entities in specified range
    /** Add to a meshset entities in specified range.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies are also added to entities in <em>entities</em>.
        \param meshset Mesh set being added to
        \param entities Range of entities being added to meshset
    */
  virtual MBErrorCode add_entities(MBEntityHandle meshset, 
                                   const MBRange &entities) = 0;

    //! Add to a meshset entities in specified vector
    /** Add to a meshset entities in specified vector.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies are also added to entities in <em>entities</em>.
        \param meshset Mesh set being added to
        \param entities 1d vector of entities being added to meshset
        \param num_entities Number of entities in 1d vector
    */
  virtual MBErrorCode add_entities(MBEntityHandle meshset, 
                                   const MBEntityHandle *entities,
                                   const int num_entities) = 0;
  
    //! Remove from a meshset entities in specified range
    /** Remove from a meshset entities in specified range.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies in entities in <em>entities</em> are updated.
        \param meshset Mesh set being removed from
        \param entities Range of entities being removed from meshset
    */
  virtual MBErrorCode remove_entities(MBEntityHandle meshset, 
                                      const MBRange &entities) = 0;

    //! Remove from a meshset entities in specified vector
    /** Remove from a meshset entities in specified vector.  If <em>meshset</em> has MESHSET_TRACK_OWNER
        option set, adjacencies in entities in <em>entities</em> are updated.
        \param meshset Mesh set being removed from
        \param entities 1d vector of entities being removed from meshset
        \param num_entities Number of entities in 1d vector
    */
  virtual MBErrorCode remove_entities(MBEntityHandle meshset, 
                                      const MBEntityHandle *entities,
                                      const int num_entities) = 0;

    //@}

    //! \name MeshSet parent/child functions

    //@{
  
    //! Get parent mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate parents are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose parents are being queried
        \param parents STL vector holding the parents returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual MBErrorCode get_parent_meshsets(const MBEntityHandle meshset,
                                          std::vector<MBEntityHandle> &parents, 
                                          const int num_hops = 1) const = 0;

    //! Get parent mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate parents are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose parents are being queried
        \param parents MBRange holding the parents returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual MBErrorCode get_parent_meshsets(const MBEntityHandle meshset,
                                          MBRange &parents,
                                          const int num_hops = 1) const = 0;

    //! Get child mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate children are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose children are being queried
        \param children STL vector holding the children returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual MBErrorCode get_child_meshsets(const MBEntityHandle meshset, 
                                         std::vector<MBEntityHandle> &children, 
                                         const int num_hops = 1) const = 0;

    //! Get child mesh sets of a mesh set
    /** If <em>num_hops</em> is 1, only immediate children are returned.  If <em>num_hops</em> is zero,
        all ancenstors are returned.  Otherwise, <em>num_hops</em> specifies the maximum number of 
        generations to traverse.
        \param meshset The mesh set whose children are being queried
        \param children MBRange holding the children returned by this function
        \param num_hops Number of generations to traverse (0 = all)
    */
  virtual MBErrorCode get_child_meshsets(const MBEntityHandle meshset, 
                                         MBRange &children, 
                                         const int num_hops = 1) const = 0;

    //! Get the number of parent mesh sets of a mesh set
    /** Identical to get_parent_meshsets, only number is returned instead of actual parents.
        \param meshset The mesh set whose parents are being queried
        \param number Number of parents
    */
  virtual MBErrorCode num_parent_meshsets(const MBEntityHandle meshset,  
                                          int *number,
                                          const int num_hops = 1) const = 0;

    //! Get the number of child mesh sets of a mesh set
    /** Identical to get_child_meshsets, only number is returned instead of actual children.
        \param meshset The mesh set whose children are being queried
        \param number Number of children
    */
  virtual MBErrorCode num_child_meshsets(const MBEntityHandle meshset, 
                                         int *number,
                                         const int num_hops = 1) const = 0;

    //! Add a parent mesh set to a mesh set
    /** Make <em>parent_meshset</em> a new parent of <em>child_meshset</em>.  This function does 
        <em>not</em> add a corresponding child link to <em>parent_meshset</em>.
        \param child_meshset The child mesh set being given a new parent.
        \param parent_meshset The parent being added to <em>child_meshset</em>
    */
  virtual MBErrorCode add_parent_meshset(MBEntityHandle child_meshset, 
                                         const MBEntityHandle parent_meshset) = 0;

    //! Add a child mesh set to a mesh set
    /** Make <em>child_meshset</em> a new child of <em>parent_meshset</em>.  This function does 
        <em>not</em> add a corresponding parent link to <em>child_meshset</em>.
        \param parent_meshset The parent mesh set being given a new child.
        \param child_meshset The child being added to <em>parent_meshset</em>
    */
  virtual MBErrorCode add_child_meshset(MBEntityHandle parent_meshset, 
                                        const MBEntityHandle child_meshset) = 0;

    //! Add parent and child links between mesh sets
    /** Makes <em>child_meshset</em> a new child of <em>parent_meshset</em>, and vica versa.
        \param parent The parent mesh set being given a new child, and the new parent
        \param child The child being given a new parent, and the new child
    */
  virtual MBErrorCode add_parent_child( MBEntityHandle parent, 
                                        MBEntityHandle child ) = 0;

    //! Remove parent and child links between mesh sets
    /** Removes parent/child links between <em>child_meshset</em> and <em>parent_meshset</em>.
        \param parent The parent mesh set being removed from <em>child</em>
        \param child The child mesh set being removed from <em>parent</em>
    */
  virtual MBErrorCode remove_parent_child( MBEntityHandle parent, 
                                           MBEntityHandle child ) = 0;

    //! Remove a parent mesh set from a mesh set
    /** Removes <em>parent_meshset</em> from the parents of <em>child_meshset</em>.  This function does 
        <em>not</em> remove a corresponding child link from <em>parent_meshset</em>.
        \param child_meshset The child mesh whose parent is being removed
        \param parent_meshset The parent being removed from <em>meshset</em>
    */
  virtual MBErrorCode remove_parent_meshset(MBEntityHandle child_meshset, 
                                            const MBEntityHandle parent_meshset) = 0;
  
    //! Remove a child mesh set from a mesh set
    /** Removes <em>child_meshset</em> from the children of <em>parent_meshset</em>.  This function does 
        <em>not</em> remove a corresponding parent link from <em>child_meshset</em>.
        \param parent_meshset The parent mesh set whose child is being removed
        \param child_meshset The child being removed from <em>parent_meshset</em>
    */
  virtual MBErrorCode remove_child_meshset(MBEntityHandle parent_meshset, 
                                           const MBEntityHandle child_meshset) = 0;

    //@}

    //! \name Error condition information 

    //@{

    //! Return information about the last error
    /** \param info std::string into which information on the last error is written.
     */
  virtual MBErrorCode get_last_error(std::string& info) const = 0;

    //! Return string representation of given error code
    /** \param code Error code for which string is wanted
     */
  virtual std::string get_error_string(const MBErrorCode code) const = 0;

    //@}
  
};

//! predicate for STL algorithms.  Returns true if the entity handle is
//! of the specified type.  For example, to remove all the tris out of a list
//! of 2D entities retrieved using get_adjacencies you could do
//! \example std::remove_if(list.begin(), list.end(), type_equals(gMB, MeshTri));
class type_equals : public std::unary_function<MBEntityHandle, bool>
{
public:
    //! interface object
  MBInterface* meshDB;

    //! type corresponding to this predicate
  const MBEntityType test_type;

    //! Constructor
  type_equals(MBInterface* mdb, const MBEntityType type) : meshDB(mdb), test_type(type){}

    //! operator predicate
  bool operator()(MBEntityHandle handle) const
    { 
      return (meshDB->type_from_handle(handle) ==  test_type); 
    } 
};

//! predicate for STL algorithms.  Returns true if the entity handle is not
//! of the specified type.  For example, to remove all but the tris out of a list
//! of 2D entities retrieved using get_adjacencies you could do
//! \example std::remove_if(list.begin(), list.end(), type_not_equals(gMB, MeshTri));
class type_not_equals : public std::unary_function<MBEntityHandle, bool>
{
public:

    //! interface object
  MBInterface* meshDB;

    //! type corresponding to this predicate
  const MBEntityType test_type;

    //! Constructor
  type_not_equals(MBInterface* mdb, const MBEntityType type) : meshDB(mdb), test_type(type){}

    //! operator predicate
  bool operator()(MBEntityHandle handle) const
    { 
      return (meshDB->type_from_handle(handle) !=  test_type); 
    } 
};

inline float MBInterface::api_version(std::string *version_string) 
{
  if (NULL != version_string)
    *version_string = std::string("MOAB API version ") + std::string(MOAB_API_VERSION_STRING);
  return MOAB_API_VERSION;
}


#endif   // MB_INTERFACE_HPP

  
