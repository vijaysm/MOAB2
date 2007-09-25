/*
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

/**\file MBFileReadTool.hpp
 *\author Jason Kraftcheck <kraftche@cae.wisc.edu>
 *\date 2007-04-25
 */

#ifndef MB_FILE_READ_TOOL_HPP
#define MB_FILE_READ_TOOL_HPP

#include "RangeMap.hpp"
#include "MBCore.hpp"

#include <iosfwd>

class MBError;
class MBEntitySequence;

/**\brief Non-template portions of MBFileReadTool */
class MBFileReadToolBase {
public:

  MBFileReadToolBase( MBCore* moab );

  /**\brief Set last error message */
  void report_error( std::string error );

  /**\brief Set last error message */
  void report_error( const char* error, ... )
#ifdef __GNUC__
  __attribute__((format(printf,2,3)))
#endif
  ;
  
  inline MBInterface* moab() const 
    { return mMOAB; }
    
  inline unsigned parallel_rank() const
    { return parallelRank; }
  
  /**\brief delete all read entities
   *
   * Delete all entities created through the methods
   * provided by this class.
   */
  MBErrorCode delete_all_entities();
  
  /**\brief Get entity set containing all read entities 
   *
   * Get a new entity set containing all the entities
   * allocated via methods of this object.
   */
  inline MBEntityHandle get_file_set() const 
    { return fileSet; }
  
  /**\brief Delete entities with specified tag
   *
   * Delete all entities read from file that do not have
   * the specified tag or tag values.
   *\param tagname  The name of the tag -- must have single integer value
   *\param tag_values If not null, a list of integer tag values
   *\param num_tag_values The length of tag_values
   */
  MBErrorCode delete_non_tagged_entities( const char* tagname,
                                          const int* tag_values = 0,
                                          int num_tag_values = 0 );

protected:

  /**\brief Allocate vertex sequence
   *
   * Allocate a vertex sequence and update any bookkeeping in
   * this object.
   *\param preferred_start_id First MBEntityID, if available
   *\param count              Number of vertices to allocate
   *\param start_handle       Output: first vertex handle
   *\param x                  Output: MOAB internal storage for vertex X-coordinates.
   *\param y                  Output: MOAB internal storage for vertex Y-coordinates.
   *\param z                  Output: MOAB internal storage for vertex Z-coordinates.
   *\param proc_id            NULL or requested parallel processor ID
   */
  MBErrorCode create_node_arrays( MBEntityID preferred_start_id,
                                  MBEntityID count,
                                  MBEntityHandle& start_handle,
                                  double *& x, double *& y, double *& z,
                                  const int* proc_id );
  
  /**\brief Allocate element sequence
   *
   * Allocate a element sequence and update any bookkeeping in
   * this object.
   *\param preferred_start_id First MBEntityID, if available
   *\param count              Number of elements to allocate
   *\param type               Element type
   *\param conn_length        Vertices per element, zero for MBPOLYGON and MBPOLYHEDRON
   *\param start_handle       Output: first vertex handle
   *\param conn_array         Output: MOAB internal storage for element connectivity
   *\param last_index_array   Cannot be NULL for MBPOLYGON or MBPOLYHEDRON.  Must 
   *                          be NULL for other types.  Output: Pointer to MOAB 
   *                          internal storage of last index in connectivity
   *                          array for each element.
   *\param proc_id            NULL or requested parallel processor ID
   */
  MBErrorCode create_element_array( MBEntityID preferred_start_id,
                                    MBEntityID count,
                                    MBEntityType type,
                                    MBEntityID conn_length,
                                    MBEntityHandle& start_handle,
                                    MBEntityHandle*& conn_array,
                                    int** last_index_array = 0,
                                    const int* proc_id = 0);
  
  /**\brief Allocate a meshset sequence
   *
   * Allocate a meshset sequence and update any bookkeeping in
   * this object.
   *\param preferred_start_id First MBEntityID, if available
   *\param count              Number of entity sets to allocate
   *\param start_handle       Output: first entity set handle
   *\param flags              entity set creation flags (one val for all entity sets)
   *\param proc_id            NULL or requested parallel processor ID
   */
  MBErrorCode create_meshset_block( MBEntityID preferred_start_id, 
                                    MBEntityID count,
                                    MBEntityHandle& start_handle_out,
                                    unsigned flags,
                                    const int* proc_id );
  
  /**\brief Allocate a meshset sequence
   *
   * Allocate a meshset sequence and update any bookkeeping in
   * this object.
   *\param preferred_start_id First MBEntityID, if available
   *\param count              Number of entity sets to allocate
   *\param start_handle       Output: first entity set handle
   *\param flags              entity set creation flags (one val per entity sets)
   *\param proc_id            NULL or requested parallel processor ID
   */
  MBErrorCode create_meshset_block( MBEntityID preferred_start_id, 
                                    MBEntityID count,
                                    MBEntityHandle& start_handle_out,
                                    const unsigned* flags,
                                    const int* proc_id );
  
  /**\brief Set tag data for a block of entity handles. */
  MBErrorCode set_tag_data( MBTag tag,
                            MBEntityHandle beg,
                            MBEntityHandle count,
                            const void* data );

private:

  MBErrorCode create_sequence( MBEntityID start_id,
                               MBEntityID count,
                               MBEntityType type,
                               MBEntityID conn_len,
                               const int* proc_id,
                               MBEntitySequence*& seq );

  MBCore* mMOAB;
  MBEntityHandle fileSet;
  unsigned parallelRank;
};

/**\brief Helper class for file readers.
 *
 * This class provides methods for creation of entities
 * in MOAB and maintains a map from file IDs to MBEntityHandles.
 * All methods accept values as file IDs and map them to MBEntityHandles
 * internally.  
 *
 * The type of a file ID is a template parameter, allowing a 
 * particular reader to work directly with the identifier/index
 * type used by an underlying IO library.
 */
template <typename FileIDType>
class MBFileReadTool : public MBFileReadToolBase
{
public:

  /**\brief Create reader tool 
   * 
   *\param vertex_unique_ids File ID space for vertices is separate from 
   *                         other entities.  If this value is true, then
   *                         the file may contain both a vertex with ID 1
   *                         and an element with the same ID value of 1.
   *\param element_type_unique_ids File ID space is separate for different
   *                         element types.  If this value is true, then
   *                         the file may contain an edge and a 
   *                         triangle element that have the same ID.
   *\param set_unique_ids File ID space for entity sets is separate from 
   *                         other entities.
   */
  inline
  MBFileReadTool( MBCore* moab_instance,
                  bool vertex_unique_ids,
                  bool element_type_unique_ids,
                  bool set_unique_ids );
  
  /**\brief Allocate a block of nodes.
   * 
   * Allocate a block of nodes and return pointers to MOAB's internal
   * storage of the vertex coordinaets for those nodes, such that the
   * file reader may read the coordinate values directly into the 
   * vertex storage allocated by this method.
   *
   *\param start_id  The file ID (or index) of the first vertex.  The
   *                 allocated vertices will be assigned file IDs in 
   *                 the range [start_id,start_id+count).
   *\param count     Number of vertices to allocate.
   *\param x_array_out  Output: location at which to write X-coordinates.
   *\param y_array_out  Output: location at which to write Y-coordinates.
   *\param z_array_out  Output: location at which to write Z-coordinates.
   *\param proc_id   For parallel, ID of processor to associate with entities.
   *                 If null, current processor rank is assumed.
   */
  inline
  MBErrorCode create_nodes( FileIDType start_id,
                            FileIDType count,
                            double *& x_array_out,
                            double *& y_array_out,
                            double *& z_array_out,
                            const int* proc_id = 0 );
  
  /**\brief Create a block of elements.
   * 
   * Create mesh elements.
   *
   *\param start_id  The file ID (or index) of the first vertex.  The
   *                 allocated vertices will be assigned file IDs in 
   *                 the range [start_id,start_id+count).
   *\param count     Number of vertices to allocate.
   *\param type      The element type.
   *\param connectivity_length Nodes per element.  Zero for MBPOLY* types.
   *\param interleaved_connectivity Element connectivity list, specified
   *                 as vertex IDs as assigned by the create_nodes method
   *                 of this class (face IDs if type == POLYHEDRON).
   *\param poly_end_indices For MBPOLYGON and MBPOLYHEDRON types, the
   *                 *last* index in interleaved_connectivity for each element.
   *\param proc_id   For parallel, ID of processor to associate with entities.
   *                 If null, current processor rank is assumed.
   */
  inline
  MBErrorCode create_elements( FileIDType start_id,
                               FileIDType count,
                               MBEntityType type,
                               FileIDType connectivity_length,
                               const FileIDType* interleaved_connectivity,
                               const FileIDType* poly_end_indices = 0,
                               const int* proc_id = 0 );

  /**\brief Specifiy explicit adjacencies
   *
   * Create explicit adjacencies between entities.
   *\param from_id   Source entity for adjacency
   *\param to_ids    Target entities for adjacencies
   *\param to_count  Length of to_ids
   *\param from_type Type of source entity.  The caller
   *                 must specify this value if it cannot
   *                 be inferred from the source entity ID
   *                 (see constructor arguments.)
   *\param to_type   Type of target entities.  The caller
   *                 must specify this value if it cannot
   *                 be inferred from the target entity IDs
   *                 (see constructor arguments.)
   */
  inline
  MBErrorCode add_adjacencies( FileIDType from_id,
                               const FileIDType* to_ids,
                               FileIDType to_count,
                               MBEntityType from_type = MBMAXTYPE,
                               MBEntityType to_type = MBMAXTYPE );
  
  /**\brief Create tag for storing entity handles.
   *
   *\param name         Tag name
   *\param num_handles  Number of handles in each tag value.  If
   *                    this value is greater than one, then each
   *                    entity has an array of handles for this tag.
   *\param storage_type Tag storage type.
   *\param handle_out   Output: tag handle
   *\param mesh_value   Global, or mesh value for tag.
   *\param default_value Default value for tag.
   *\param null_id      If non-null, the value pointed to by this
   *                    argument is treated as special: a "NULL" file ID.
   */
  inline 
  MBErrorCode create_entity_handle_tag( const char* name,
                                        unsigned num_handles,
                                        MBTagType storage_type,
                                        MBTag& handle_out,
                                        const FileIDType* mesh_value = 0,
                                        const FileIDType* default_value = 0,
                                        const FileIDType* null_id = 0 );
  
  /**\brief Set tag data for a block of entities
   *
   * Set tag data for non-handle type tags.  This method will
   * fail if the tag has a data type MB_TAG_HANDLE.
   *
   *\param handle       Tag handle
   *\param start_id     File ID of first entity
   *\param id_count     Number of entities.
   *\param tag_data     Tag values
   *\param entity_type  The type of the entities for which the tag
   *                    data is being set.  The caller must specify
   *                    this if it cannot be inferred from the file
   *                    IDs (see constructor arguments).
   *\param value_type   The type of the entities in the tag values.  
   *                    The caller must specify
   *                    this if it cannot be inferred from the file
   *                    IDs (see constructor arguments).
   */
  inline
  MBErrorCode set_tag_data( MBTag handle,
                            FileIDType start_id,
                            FileIDType id_count,
                            const void* tag_data,
                            MBEntityType entity_type = MBMAXTYPE,
                            MBEntityType val_type = MBMAXTYPE );
  
  /**\brief Set tag data for entity handle tag for a block of entities
   *
   * This method will fail for any tag that does not have a data type
   * of MB_TYPE_HANDLE.
   *
   *\param handle       Tag handle
   *\param start_id     File ID of first entity
   *\param id_count     Number of entities.
   *\param tag_data     Tag values as file IDs
   *\param entity_type  The type of the entities for which the tag
   *                    data is being set.  The caller must specify
   *                    this if it cannot be inferred from the file
   *                    IDs (see constructor arguments).
   *\param value_type   The type of the entities in the tag values.  
   *                    The caller must specify
   *                    this if it cannot be inferred from the file
   *                    IDs (see constructor arguments).
   *\param null_id      If non-null, the value pointed to by this
   *                    argument is treated as special: a "NULL" file ID.
   */
  inline
  MBErrorCode set_tag_data( MBTag handle,
                            const FileIDType* id_array,
                            FileIDType id_count,
                            const FileIDType* tag_data,
                            MBEntityType entity_type = MBMAXTYPE,
                            MBEntityType value_type = MBMAXTYPE,
                            const FileIDType* null_id = 0 );

  /**\brief Create a block of entity sets
   *
   * Create count entity sets with consecutive file IDs, beginning
   * with the specified start_id.
   *\param start_id   File ID for first entity set
   *\param count      Number of entity sets to allocate
   *\param options    Set creation flags for each entity set.
   *\param proc_id    For parallel, ID of processor to associate with entities.
   *                  If null, current processor rank is assumed.
   */
  inline
  MBErrorCode create_meshsets( FileIDType start_id,
                               MBEntityID count,
                               const unsigned* options,
                               const int* proc_id = 0 );
  
  /**\brief Create a block of entity sets
   *
   * Create count entity sets with consecutive file IDs, beginning
   * with the specified start_id.
   *\param start_id   File ID for first entity set
   *\param count      Number of entity sets to allocate
   *\param options    Set creation flags for each entity set.
   *\param proc_id    For parallel, ID of processor to associate with entities.
   *                  If null, current processor rank is assumed.
   */
  inline
  MBErrorCode create_meshsets( FileIDType start_id,
                               FileIDType count,
                               unsigned options,
                               const int* proc_id = 0 );
  
  /**\brief Specify set contents for a block of entity sets
   *
   *\param start_id  File id of first entity set
   *\param count     Number of entity sets
   *\param concatenated_content_lists The concatenation of the list of set
   *                 contents.
   *\param end_indices For each entity set, the index of the last entry
   *                 for that set in concatenated_content_lists
   *\param type      The type of the entities in concatenated_content_list.
   *                 The caller must specify this value if it cannot be
   *                 inferred from the file IDs (see constructor args.)
   */
  inline
  MBErrorCode add_meshset_contents( FileIDType start_id,
                                    FileIDType count,
                                    const FileIDType* concatenated_content_lists,
                                    const FileIDType* end_indices ,
                                    MBEntityType type = MBMAXTYPE );
                                    
  /**\brief Specify set contents for a block of entity sets
   *
   * Specify contents for a block of entity sets, where the set
   * contents are specified as ranges of file IDs.
   *
   *\param start_id  File id of first entity set
   *\param count     Number of entity sets
   *\param concatenated_range_list The concatenation of the list of set
   *                 contents.  The contents for each set are specified
   *                 as pairs of file IDS, where each pair corresponds
   *                 to the first and last file ID of a range of entities.
   *                 To specify a single entity, both values in the range
   *                 pair should be the same.
   *\param end_indices For each entity set, the index of the last entry
   *                 for that set in concatenated_range_list
   *\param type      The type of the entities in concatenated_range_list.
   *                 The caller must specify this value if it cannot be
   *                 inferred from the file IDs (see constructor args.)
   */
  inline
  MBErrorCode add_meshset_range_contents( FileIDType start_id,
                                          FileIDType count,
                                          const FileIDType* concatenated_range_list,
                                          const FileIDType* end_indices,
                                          MBEntityType type = MBMAXTYPE );

  /**\brief Specify set children for a block of entity sets
   *
   *\param start_parent_id  File id of first entity set
   *\param parent_count     Number of entity sets
   *\param concatenated_children_lists The concatenation of the list of set
   *                        children.
   *\param end_indices      For each entity set, the index of the last entry
   *                        for that set in concatenated_children_lists
   *\param null_id          If non-null, the value pointed to by this
   *                        argument is treated as special: a "NULL" file ID.
   *\param bidirectional    If true, also create link from child up to parent.
   *                        If false, links are single-directional.
   */
  inline
  MBErrorCode add_meshset_children( FileIDType start_parent_id,
                                    FileIDType parent_count,
                                    const FileIDType* concatenated_children_lists,
                                    const FileIDType* end_indices,
                                    bool bidirectional = true );

  /**\brief Specify set parents for a block of entity sets
   *
   *\param start_child_id  File id of first entity set
   *\param child_count     Number of entity sets
   *\param concatenated_parents_lists The concatenation of the list of set
   *                       parents.
   *\param end_indices     For each entity set, the index of the last entry
   *                       for that set in concatenated_parents_lists
   *\param null_id         If non-null, the value pointed to by this
   *                       argument is treated as special: a "NULL" file ID.
   */
  inline
  MBErrorCode add_meshset_parents(  FileIDType start_child_id,
                                    FileIDType child_count,
                                    const FileIDType* concatenated_parents_lists,
                                    const FileIDType* end_indices );
  
  /**\brief Get number of entities of each type */
  inline void get_entity_counts( FileIDType counts[MBMAXTYPE] ) const;


  /**\brief Convert from file IDs to MBEntityHandles
   *
   * Convert from file IDs to MBEntityHandles
   *
   *\param type    The type of the entities, or MBMAXTYPE if
   *               the entity type can be inferred from the file ID.
   *               The type can be inferred if all arguments to the
   *               constructor that control the file ID space mapping
   *               are false.  I.e. there is to potential ID conflict.
   *\param ids     The input array of IDs for which to get the 
   *               corresponding MBEntityHandles.
   *\param handles The array in which to store the result MBEntityHandles
   *\param count   The length of ids (and handles).
   *\param null_id If this argument is not NULL, then the value it points
   *               to will be assumed to be a 'NULL' file ID value, such
   *               that any input FileIDType with this value will be
   *               assigned a zero (NULL) MBEntityHandle.
   *\return        true if success, false if failure.
   */
  inline bool ids_to_handles( MBEntityType type,
                              const FileIDType* ids, 
                              MBEntityHandle* handles,
                              FileIDType count,
                              const FileIDType* null_id = 0 );
  
private:
  RangeMap<FileIDType,MBEntityHandle> typeMap[MBMAXTYPE+1];
  const bool vertexUniqueIds, elemTypeUniqueIds, setUniqueIds;
  std::vector<MBEntityHandle> handleArray;
};

template <typename FileIDType> inline
bool MBFileReadTool<FileIDType>::ids_to_handles( MBEntityType type,
                                                 const FileIDType* ids, 
                                                 MBEntityHandle* handles,
                                                 FileIDType count,
                                                 const FileIDType* null_id  )
{
  bool ok = true;
  const RangeMap<FileIDType,MBEntityHandle>& map = typeMap[type];
  for (FileIDType i = 0; i < count; ++i) {
    handles[i] = map.find( ids[i] );
    if (!handles[i] && (!null_id || *null_id != ids[i]))
      ok = false;
  }
  return ok;
} 

template <typename FileIDType> inline
MBFileReadTool<FileIDType>::MBFileReadTool( MBCore* moab_iface,
                                            bool vertex_unique_ids,
                                            bool element_type_unique_ids,
                                            bool set_unique_ids )
  : MBFileReadToolBase( moab_iface ),
    vertexUniqueIds( vertex_unique_ids ),
    elemTypeUniqueIds( element_type_unique_ids ),
    setUniqueIds( set_unique_ids )
  {}
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::create_nodes( FileIDType start_id,
                                                      FileIDType count,
                                                      double *& x,
                                                      double *& y,
                                                      double *& z,
                                                      const int* proc_id )
{
  MBEntityHandle handle;
  MBErrorCode rval = create_node_arrays( start_id, count, handle, x, y, z, proc_id );
  if (MB_SUCCESS != rval) 
    return rval;
  
  typeMap[MBVERTEX].insert( handle, count, start_id );
  if (!vertexUniqueIds)
    typeMap[MBMAXTYPE].insert( handle, count, start_id );
}
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::create_elements( FileIDType start_id,
                                                         FileIDType count,
                                                         MBEntityType type,
                                                         FileIDType len,
                                                         const FileIDType* conn_ids,
                                                         const FileIDType* poly,
                                                         const int* proc_id )
{
  if (type == MBPOLYGON || type == MBPOLYHEDRON)  
    if (!poly || poly[count-1] != len)
      return MB_FAILURE;
  
  int* index_array;
  MBErrorCode rval;
  MBEntityHandle handle, *conn;
  
  rval = create_element_array( start_id, count, type, len, handle, conn, &index_array, proc_id );
  if (MB_SUCCESS != rval)
    return rval;
    
  typeMap[type].insert( handle, count, start_id );
  if (!elemTypeUniqueIds)
    typeMap[MBMAXTYPE].insert( handle, count, start_id );
  
  bool valid;
  if (type != MBPOLYHEDRON)
    valid = convert_ids_to_handles( MBVERTEX, conn_ids, conn, count );
  else if (elemTypeUniqueIds)
    valid = convert_ids_to_handles( MBPOLYGON, conn_ids, conn, count );
  else
    valid = convert_ids_to_handles( MBMAXTYPE, conn_ids, conn, count );
  
  if (type == MBPOLYGON || type == MBPOLYHEDRON)  
    for (FileIDType i = 0; i < count; ++i)
      index_array[i] = poly[i];
  
  return valid ? MB_SUCCESS : MB_FAILURE;;
}
    
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::add_adjacencies( FileIDType from_id,
                                              const FileIDType* to_ids,
                                              FileIDType to_count,
                                              MBEntityType from_type,
                                              MBEntityType to_type )
{
  MBEntityHandle from_handle;
  if (!convert_ids_to_handles( from_type, &from_id, &from_handle, 1 ))
    return MB_FAILURE;
  
  handleArray.resize( to_count );
  if (!convert_ids_to_handles( to_type, to_ids, &handleArray[0], to_count ))
    return MB_FAILURE;
  
  return moab()->add_adjacencies( from_handle, &handleArray[0], to_count, false );
}
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::create_entity_handle_tag( 
                                        const char* name,
                                        unsigned num_handles,
                                        MBTagType storage_type,
                                        MBTag& handle_out,
                                        const FileIDType* mesh_value,
                                        const FileIDType* default_value,
                                        const FileIDType* null_id )
{
  handleArray.resize( 2*num_handles );
  MBEntityHandle *mesh_handles = 0, *default_handles = 0;
  if (mesh_value) {
    mesh_handles = &handleArray[0];
    if (!ids_to_handles( MBMAXTYPE, mesh_value, mesh_handles, num_handles, null_id ))
      return MB_FAILURE;
  }
  if (default_value) {
    default_handles = &handleArray[num_handles];
    if (!ids_to_handles( MBMAXTYPE, default_value, default_handles, num_handles, null_id ))
      return MB_FAILURE;
  }
  
  MBErrorCode rval =  moab()->tag_create( name, 
                             num_handles * sizeof(MBEntityHandle),
                             storage_type,
                             MB_TYPE_HANDLE,
                             handle_out,
                             default_handles,
                             true );
  if (MB_SUCCESS != rval)
    return rval;
  
  return set_tag_data( handle_out, 0, 0, mesh_handles );
}
  
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::set_tag_data( MBTag handle,
                                                    FileIDType start_id,
                                                    FileIDType id_count,
                                                    const void* tag_data,
                                                    MBEntityType type,
                                                    MBEntityType val_type )
{
  MBErrorCode rval, result = MB_SUCCESS;
  std::vector<MBEntityHandle> ranges, hdata;
  if (!typeMap[type].find( start_id, id_count, ranges ))
    return MB_FAILURE;
  
  MBDataType tag_type;
  rval = moab()->tag_get_data_type( handle, tag_type );
  if (MB_SUCCESS != rval)
    return rval;
  
  int tag_size;
  rval = moab()->tag_get_size( handle, tag_size );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (MB_TYPE_HANDLE == tag_type)
    return MB_TYPE_OUT_OF_RANGE;
                      
  const char* data = static_cast<const char*>(tag_data);
  for (MBEntityID i = 0; i < ranges.size(); i += 2) {
    MBEntityHandle start_handle = ranges[i];
    MBEntityID count = ranges[i+1];
    rval = set_tag_data( handle, start_handle, count, data );
    if (MB_SUCCESS != rval)
      result = rval;
    data += (count + tag_size);
  }
 
  return result;
}
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::set_tag_data( MBTag handle,
                                                    const FileIDType* id_array,
                                                    FileIDType id_count,
                                                    const FileIDType* tag_data,
                                                    MBEntityType type,
                                                    MBEntityType val_type,
                                                    const FileIDType* null_id )
{
  MBErrorCode rval;
  MBDataType tag_type;
  rval = moab()->tag_get_data_type( handle, tag_type );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (tag_type != MB_TYPE_HANDLE) 
    return MB_TYPE_OUT_OF_RANGE;

  int tag_size;
  rval = moab()->tag_get_size( handle, tag_size );
  if (MB_SUCCESS != rval)
    return rval;

  tag_size /= sizeof(MBEntityHandle);
  handleArray.resize( id_count * (tag_size + 1) );
  if (!ids_to_handles( val_type, 
                       tag_data,
                       &handleArray[id_count],
                       null_id ))
    return MB_FAILURE;

  tag_data = &handleArray[id_count];

  if (!ids_to_handles( type, id_array, &handleArray[0], id_count ))
    return MB_FAILURE;
  
  return set_tag_data( handle, &handleArray[0], id_count, tag_data );
}


template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::create_meshsets( FileIDType start_id,
                               MBEntityID count,
                               const unsigned* options,
                               const int* proc_id )
{
  MBEntityHandle h;
  MBErrorCode rval;
  
  rval = craete_meshset_block( start_id, count, h, options, proc_id );
  if (MB_SUCCESS != rval)
    return rval;
  
  typeMap[MBENTITYSET].insert( h, count, start_id );
  if (!setUniqueIds)
    typeMap[MBMAXTYPE].insert( h, count, start_id );
  
  return MB_SUCCESS;
}
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::create_meshsets( FileIDType start_id,
                               FileIDType count,
                               unsigned options,
                               const int* proc_id )
{
  MBEntityHandle h;
  MBErrorCode rval;
  
  rval = craete_meshset_block( start_id, count, h, options, proc_id );
  if (MB_SUCCESS != rval)
    return rval;
  
  typeMap[MBENTITYSET].insert( h, count, start_id );
  if (!setUniqueIds)
    typeMap[MBMAXTYPE].insert( h, count, start_id );
  
  return MB_SUCCESS;
}
  
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::add_meshset_contents( 
                                    FileIDType start_id,
                                    FileIDType count,
                                    const FileIDType* list,
                                    const FileIDType* end_indices,
                                    MBEntityType type )
{
  MBEntityHandle handle;
  FileIDType id, prev = 0, n;
  MBErrorCode result = MB_SUCCESS;
  
  for (FileIDType i = 0; i < count; prev = end_indices[i++]) {
    id = start_id + i;
    if (!ids_to_handles( MBENTITYSET, &id, &handle, 1 )) {
      result = MB_FAILURE;
      continue;
    }
    
    n = end_indices[i] - prev;
    handleArray.resize( n );
    if (!ids_to_handles( type, list+prev, &handleArray[0], n))
      result = MB_FAILURE;
    else if (MB_SUCCESS != moab()->add_entities( handle, &handleArray[0], n ))
      result = MB_FAILURE;
  }
  return result;
}

                                    
template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::add_meshset_range_contents(
                                          FileIDType start_id,
                                          FileIDType count,
                                          const FileIDType* list,
                                          const FileIDType* end_indices,
                                          MBEntityType type )
{
  MBEntityHandle handle;
  FileIDType id, prev = 0, n;
  MBErrorCode rval, result = MB_SUCCESS;
  MBRange range;
  MBRange::iterator h;
  
  for (FileIDType i = 0; i < count; ++i) {
    id = start_id + i;
    if (!ids_to_handles( MBENTITYSET, &id, &handle, 1 )) {
      result = MB_FAILURE;
      prev = end_indices[i];
      continue;
    }
    
    range.clear();
    h = range.begin();
    for ( ; prev < end_indices[i]; prev += 2) 
      typeMap[type].find( list[i], list[i+i], range, h );
    
    if (range.size() != count)
      result = MB_FAILURE;
    
    rval = moab()->add_entities( handle, range );
    if (MB_SUCCESS != rval)
      result = rval;
  }
    
  return result;
}


template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::add_meshset_children( 
                                    FileIDType start_id,
                                    FileIDType count,
                                    const FileIDType* list,
                                    const FileIDType* end_indices,
                                    bool bidirectional )
{
  MBEntityHandle handle;
  FileIDType id, prev = 0, n;
  MBErrorCode result = MB_SUCCESS;
  
  for (FileIDType i = 0; i < count; prev = end_indices[i++]) {
    id = start_id + i;
    if (!ids_to_handles( MBENTITYSET, &id, &handle, 1 )) {
      result = MB_FAILURE;
      continue;
    }
    
    n = end_indices[i] - prev;
    handleArray.resize( n );
    if (!ids_to_handles( MBENTITYSET, list+prev, &handleArray[0], n))
      result = MB_FAILURE;
    else {
      if (bidirectional) {
        for (FileIDType i = 0; i < n; ++i) 
          if (MB_SUCCESS != moab()->add_parent_child( handle, handleArray[i] ))
            result = MB_FAILURE;
      }
      else if (MB_SUCCESS != moab()->add_child_meshsets( handle, &handleArray[0], n ))
        result = MB_FAILURE;
    }
  }
  return result;
}

template <typename FileIDType> inline
MBErrorCode MBFileReadTool<FileIDType>::add_meshset_parents(  
                                    FileIDType start_id,
                                    FileIDType count,
                                    const FileIDType* list,
                                    const FileIDType* end_indices )
{
  MBEntityHandle handle;
  FileIDType id, prev = 0, n;
  MBErrorCode result = MB_SUCCESS;
  
  for (FileIDType i = 0; i < count; prev = end_indices[i++]) {
    id = start_id + i;
    if (!ids_to_handles( MBENTITYSET, &id, &handle, 1 )) {
      result = MB_FAILURE;
      continue;
    }
    
    n = end_indices[i] - prev;
    handleArray.resize( n );
    if (!ids_to_handles( MBENTITYSET, list+prev, &handleArray[0], n))
      result = MB_FAILURE;
    else if (MB_SUCCESS != moab()->add_parent_meshsets( handle, &handleArray[0], n ))
      result = MB_FAILURE;
  }
  return result;
}

template <typename FileIDType> inline
void MBFileReadTool<FileIDType>::get_entity_counts( FileIDType counts[MBMAXTYPE] ) const
{
  for (MBEntityType t = MBVERTEX; t < MBMAXTYPE; ++t)
    counts[t] = typeMap[t].size();
}

  
#endif
