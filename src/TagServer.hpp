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


/**********************************************
 * Filename   :     TagServer.hpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/



#ifndef TAG_SERVER_HPP
#define TAG_SERVER_HPP

#ifndef IS_BUILDING_MB
#error "TagServer.hpp isn't supposed to be included into an application"
#endif

#include <string>
#include <vector>

#include "moab/Types.hpp"
#include "moab/Range.hpp"
#include "Internals.hpp"
#include "TagInfo.hpp"

namespace moab {

class SparseTagSuperCollection;
class SequenceManager;
class BitTagServer;

//! SparseTagServer class which associates tag data with entities
class TagServer
{
public:
  //! constructor
  TagServer( SequenceManager* seqman );
  //! destructor
  virtual ~TagServer();

  //! add a tag
  ErrorCode add_tag(const char *tag_name, 
                       const int data_size,
                       const TagType storage,
                       const DataType data_type,
                       Tag &tag_handle,
                       const void *default_data = NULL,
                       int default_value_size = 0);

  // tag name is the name of the tag
  // entity_type is the type of entity (in implementation - it means
  // tag data is stored in a different namespace)
  // entity_handle is the entity to which data is tagged to
  // tag_usage defines which storage mechanism to use
  // data is the actual data 
  // data_size is the size of the data in bytes

  //! remove a tag from this tag server
  ErrorCode remove_tag(const Tag tag_handle);

  //! resets all tag data associated with this handle back to default data or NULL
  ErrorCode reset_data(EntityHandle entity_handle);

  //! cleans out all data tagged on all entities
  ErrorCode reset_all_data();

  //! set global/mesh value of tag
  ErrorCode set_mesh_data( const Tag tag_handle, const void* data, int size = 0);

  /**\brief Set value for {Tag,EntityHandle) tuple.
   * 
   * Set tag value.
   *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH for variable-length tags.
   *\param tag_handle    The tag.
   *\param entity_handle The entity.
   *\param data          Pointer to tag value.
   */
  ErrorCode set_data(const Tag tag_handle, const EntityHandle entity_handle, const void* data )
    { return set_data( tag_handle, &entity_handle, 1, data ); }
  
  /**\brief Set tag values for an array of entity handles.
   * 
   * Set tag values.
   *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH for variable-length tags.
   *\param tag_handle     The tag.
   *\param entity_handles Array of entity handles.
   *\param num_entities   Length of entity_handles array.
   *\param data           Pointer to memory containing concatenation of
   *                      tag values for all entities, in the order the
   *                      entities are specified in entity_handles.
   */
  ErrorCode set_data(const Tag tag_handle, const EntityHandle* entity_handles, const int num_entities, const void* data );
  
  /**\brief Set tag values for an Range of entity handles.
   * 
   * Set tag values.
   *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH for variable-length tags.
   *\param tag_handle     The tag.
   *\param entity_handles Range of entity handles.
   *\param data           Pointer to memory containing concatenation of
   *                      tag values for all entities, in the order the
   *                      entities are specified in entity_handles.
   */
  ErrorCode set_data(const Tag tag_handle, const Range& entity_handles, const void* data );
  
  /**\brief Set tag values for an array of entity handles.
   * 
   * Set tag values.
   *\param tag_handle     The tag.
   *\param entity_handles Array of entity handles.
   *\param num_entities   Length of entity_handles array.
   *\param data           Array of pointers to per-entity tag values.
   *\param lengths        Length of each entity's tag value.  Ignored
   *                      for fixed-length tags.
   */
  ErrorCode set_data( const Tag tag_handle, 
                        const EntityHandle* entity_handles, 
                        const int num_entities, 
                        void const* const* data,
                        const int* lengths = 0 );
  
  /**\brief Set tag values for an Range of entity handles.
   * 
   * Set tag values.
   *\param tag_handle     The tag.
   *\param entity_handles entity handles.
   *\param data           Array of pointers to per-entity tag values.
   *\param lengths        Length of each entity's tag value.  Ignored
   *                      for fixed-length tags.
   */
  ErrorCode set_data( const Tag tag_handle, 
                        const Range& entity_handles, 
                        void const* const* data,
                        const int* lengths = 0 );

  /**\brief Set tags on all entities to same value
   *\param tag The tag to set
   *\param entities Entities on which to set tag
   *\param value    Pointer to tag value.
   *\param length   For variable-length tags, the length of the
   *                tag value.  Ignored otherwise.
   */
  ErrorCode clear_data( Tag tag, 
                        const Range& entities, 
                        const void* value, 
                        int length = 0 );
                        
  /**\brief Set tags on all entities to same value
   *\param tag The tag to set
   *\param entities Entities on which to set tag
   *\param value    Pointer to tag value.
   *\param length   For variable-length tags, the length of the
   *                tag value.  Ignored otherwise.
   */
  ErrorCode clear_data( Tag tag, 
                        const EntityHandle* entities, 
                        int num_entities, 
                        const void* value, 
                        int length = 0 );

  //! get global/mesh value of tag
  ErrorCode get_mesh_data( const Tag tag_handle, void* data ) const;

  //! get pointer/reference to mesh data
  ErrorCode get_mesh_data( Tag tag_handle, const void*& data_ptr, int& size ) const;

  /**\Brief Get tag value
   *
   * Get the value for a {Tag,EntityHandle} tuple.
   *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH for variable-length tags.
   *\param tag_handle    The tag
   *\param entity_handle The entity
   *\param data          Pointer to memory location to which to copy tag value.
   */
  ErrorCode get_data(const Tag tag_handle, const EntityHandle entity_handle, void* data )
    { return get_data( tag_handle, &entity_handle, 1, data ); }
  
  /**\Brief Get tag values
   *
   * For a single tag, get the tag value for an array of entities.
   *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH for variable-length tags.
   *\param tag_handle     The tag
   *\param entity_handles Entity handle array
   *\param num_entiites   Length of entity_handles array
   *\param data           Pointer to memory location to which to copy tag values.
   *                      Writes the concatenation of tag values, in the order
   *                      of the entity handles in the input array.
   */
  ErrorCode get_data(const Tag tag_handle, const EntityHandle* entity_handles, const int num_ents, void* data );
  
  /**\Brief Get tag values
   *
   * For a single tag, get the tag value for an Range of entities.
   *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH for variable-length tags.
   *\param tag_handle     The tag
   *\param entity_handles Entity handles
   *\param data           Pointer to memory location to which to copy tag values.
   *                      Writes the concatenation of tag values, in the order
   *                      of the entity handles in the input array.
   */
  ErrorCode get_data(const Tag tag_handle, const Range& entity_handles, void* data );
  
  /**\brief Get pointers to tag values for an array of entity handles.
   * 
   * Get pointers to tag values.
   *\param tag_handle     The tag.
   *\param entity_handles Array of entity handles.
   *\param num_entities   Length of entity_handles array.
   *\param data           Output: Array of pointers to per-entity tag values.
   *\param lengths        Output: Length of each entity's tag value.  
   *                      Optional for fixed-length tags.
   */
  ErrorCode get_data( const Tag tag_handle, 
                        const EntityHandle* entity_handles, 
                        const int num_entities, 
                        const void** data,
                        int* lengths = 0 );
  
  /**\brief Get pointers to tag values for an Range of entity handles.
   * 
   * Get pointers to tag values.
   *\param tag_handle     The tag.
   *\param entity_handles entity handles.
   *\param data           Output: Array of pointers to per-entity tag values.
   *\param lengths        Output: Length of each entity's tag value.  
   *                      Optional for fixed-length tags.
   */
  ErrorCode get_data( const Tag tag_handle, 
                        const Range& entity_handles, 
                        const void** data,
                        int* lengths = 0 );

  /**\brief Access tag data via direct pointer into contiguous blocks
   *
   * Iteratively obtain direct access to contiguous blocks of tag
   * storage.  This function cannot be used with bit tags because
   * of the compressed bit storage.  This function cannot be used
   * with variable length tags because it does not provide a mechanism
   * to determine the length of the value for each entity.  This
   * function may be used with sparse tags, but if it is used, it
   * will return data for a single entity at a time.  
   *
   *\param tag_handle  The handle of the tag for which to access data
   *\param iter        As input, the first entity for which to return
   *                   data.  As output, one past the last entity for
   *                   which data was returned.
   *\param end         One past the last entity for which data is desired
   *\param data_ptr    Output: pointer to tag storage.
   *  
   *\Note If this function is called for entities for which no tag value
   *      has been set, but for which a default value exists, it will 
   *      force the allocation of explicit storage for each such entity
   *      even though MOAB would normally not explicitly store tag values
   *      for such entities.
   */
  ErrorCode tag_iterate( Tag tag_handle,
                         Range::iterator& iter,
                         const Range::iterator& end,
                         void*& data_ptr );


  //! remove global/mesh value of tag
  ErrorCode remove_mesh_data( const Tag tag_handle );

  //! remove the tag data on an entity
  ErrorCode remove_data( const Tag tag_handle, const EntityHandle entity_handle );

  //! gets all entity handles that match a type and tag
  ErrorCode get_entities( const Tag tag_handle, 
                             const EntityType type,
                             Range &entities);

  //! gets all entity handles that match a tag
  ErrorCode get_entities( const Tag tag_handle, 
                             Range &entities);

  //! For the set of entities in the input range, return those
  //! that match the specified type and have a value for the 
  //! specified tag.
  ErrorCode get_entities( const Range &input_range,
                             const Tag tag_handle, 
                             const EntityType type,
                             Range &entities);

  //! gets all entity handles that match a type, tag and tag value 
  ErrorCode get_entities_with_tag_value( const EntityType type,
                                            const Tag tag_handle,
                                            const void* value,
                                            Range &entities,
                                            int value_size = 0 );
  
  //! For the set of entities in the input range, return those
  //! that match the specified type and have the specified tag value.
  ErrorCode get_entities_with_tag_value( const Range &input_range,
                                            const EntityType type,
                                            const Tag tag_handle,
                                            const void* value,
                                            Range &entities,
                                            int value_size = 0 );
  
  ErrorCode get_entities_with_tag_values( const Range &input_range,
                                             const EntityType type,
                                             const Tag *tags,
                                             const void* const* values,
                                             const int num_tags,
                                             Range &entities,
                                             const int condition);
  
  //! gets number of entities that match a type and tag
  ErrorCode get_number_entities( const Tag tag_handle, const EntityType type,
                             int& num_entities);

  //! get number of entities with tag set
  ErrorCode get_number_entities( const Tag tag_handle, unsigned long& num_ents );

  //! gets number of entities that match a type and tag
  ErrorCode get_number_entities(const Range &input_range,
                                   const Tag tag_handle, const EntityType type,
                                    int& num_entities);

  //! gets a tag handle by name and entity handle
  Tag get_handle(const char *tag_name) const;

  //! get all the tags which have been defined for this entity
  ErrorCode get_tags(const EntityHandle entity, std::vector<Tag> &all_tags);
 
  //! get all the tags which have a global/mesh value
  ErrorCode get_mesh_tags(std::vector<Tag> &all_tags) const;
 
  //! get all the tags which have been defined
  ErrorCode get_tags(std::vector<Tag> &all_tags);
  
    //! get the default value for a given tag
  ErrorCode get_default_data(const Tag tag_handle, void *data, int& size);

    //! get the default value for a given tag
  ErrorCode get_default_data_ref(const Tag tag_handle, const void *& data, int& size);

  //! get information about a tag
  inline const TagInfo* get_tag_info(const char *tag_name ) const;
  inline const TagInfo* get_tag_info( Tag tag_handle ) const;
  inline TagInfo* get_tag_info( Tag tag_handle );
  inline const TagInfo* get_tag_info( TagId id, TagType storage ) const;
  inline TagInfo* get_tag_info( TagId id, TagType storage );
  
  unsigned long get_memory_use( Tag tag_handle ) const;
  
  ErrorCode get_memory_use( Tag tag_handle,
                              unsigned long& total,
                              unsigned long& per_entity ) const;

private:

  //! Table of tag ids and tag information
  //! Primary (array) index is tag type.  
  //! Secondary (std::vector) index is tag id less one (tag ids begin with 1).
  std::vector<TagInfo> mTagTable[MB_TAG_LAST+1];

  //! container for storing the sparse data and tag ids
  SparseTagSuperCollection* mSparseData;
  
  //! SequenceManager required to access dense tag data
  SequenceManager* sequenceManager;

  //! manager for the bit data
  BitTagServer* mBitServer;

};

inline const TagInfo* TagServer::get_tag_info( const char *tag_name ) const
{
  const Tag handle = get_handle( tag_name );
  return handle ? get_tag_info( handle ) : 0;
}

inline const TagInfo* TagServer::get_tag_info( Tag tag ) const
{
  return get_tag_info( ID_FROM_TAG_HANDLE(tag), PROP_FROM_TAG_HANDLE(tag) );
}

inline TagInfo* TagServer::get_tag_info( Tag tag )
{
  return get_tag_info( ID_FROM_TAG_HANDLE(tag), PROP_FROM_TAG_HANDLE(tag) );
}

inline const TagInfo* TagServer::get_tag_info( TagId id, TagType type ) const
{
  if (id > 0 && id <= mTagTable[type].size() && mTagTable[type][id-1].is_valid())
    return &mTagTable[type][id-1];
  else
    return NULL;
}

inline TagInfo* TagServer::get_tag_info( TagId id, TagType type )
{
  if (id > 0 && id <= mTagTable[type].size() && mTagTable[type][id-1].is_valid())
    return &mTagTable[type][id-1];
  else
    return NULL;
}

} // namespace moab

#endif //TAG_SERVER_HPP




