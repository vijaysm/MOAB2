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
 * Filename   :     SparseTagSuperCollection.hpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/



#ifndef SPARSE_TAG_SUPER_COLLECTION_HPP
#define SPARSE_TAG_SUPER_COLLECTION_HPP

#ifndef IS_BUILDING_MB
#error "SparseTagSuperCollection.hpp isn't supposed to be included into an application"
#endif

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include <map>
#include <vector>

#include "moab/Types.hpp"
#include "Internals.hpp"
#include "moab/Range.hpp"
#include "TagInfo.hpp"

#define get_collection( A ) ((A) < mDataTags.size() ? mDataTags[(A)] : 0)

namespace moab {

class SparseTagCollection;

//! a collection of SparseTagCollections
class SparseTagSuperCollection
{
public:
  //! constructor
  SparseTagSuperCollection(){}
  
  //! destructor
  virtual ~SparseTagSuperCollection();

  void reset_data();

  //! allocate new tag id
  ErrorCode reserve_tag_id(int data_size, TagId tag_id);

  //! releases an TagId
  ErrorCode release_tag_id(TagId tag_id);

  //! size of data tag
  int tag_size(const TagId tag_id) const;

    /** Set fixed-length tag values.
     *\NOTE Will fail for variable-length tag data.
     */
  ErrorCode set_data( TagId tag_handle,
                        const EntityHandle* handles,
                        int num_handles,
                        const void* data );

    /** Set tag values for array of entity handles
     *\param tag_id      The tag.
     *\param handles     Array of entity handles.
     *\param num_handles Length of 'handles' array.
     *\param values      Array of pointers to tag values, one pointer for each handle
     *\param lengths     Length of each tag value.  Ignored for fixed-length tags.
     */
  ErrorCode set_data( TagId tag_handle,
                        const EntityHandle* handles,
                        int num_handles,
                        void const* const* data_pointers,
                        const int* lengths = 0 );

    /** Set fixed-length tag value for an Range of entities
     *\NOTE Will fail for variable-length tag data
     */
  ErrorCode set_data( TagId tag_handle,
                        const Range& handles,
                        const void* data );

    /** Set tag data for an Range of entities.
     *
     *\param tag_id  The tag
     *\param handles The entities
     *\param values  An array of pointers, one per entity, pointing to the
     *               tag value for the corresponding entity.
     *\param lengths An array of integers, one per entity, indicating the
     *               length of the tag value for each entity.  Ingored
     *               for fixed-length tags.
     */
  ErrorCode set_data( TagId tag_handle,
                        const Range& handles,
                        void const* const* data_pointers,
                        const int* lengths = 0 );
  
    /** Get fixed-length tag values for array of entities
     *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called
     *      for variable-length tag.
     */
  ErrorCode get_data( TagId tag_handle,
                        const EntityHandle* handles,
                        int num_handles,
                        void* data,
                        const void* default_value ) const;
  
    /** Get pointers to tag data for array of entities
     *\param tag_id      The Tag.
     *\param handles     Array of entity handles.
     *\param num_handles Length of 'handles' array.
     *\param tag_ptrs    Pointers to tag values, one pointer for each input handle.
     *\param lengths     Length of each tag value.  Ignored for fixed-length tags.
     *\param default_value Pointer to default value for tag, or NULL if none.
     *\param default_value_length  Length of default tag value.  Ingored for
     *                   fixed-length tags.
     */
  ErrorCode get_data( TagId tag_handle,
                        const EntityHandle* handles,
                        int num_handles,
                        const void** data,
                        int* lengths,
                        const void* default_value,
                        int default_value_length ) const;
  
    /** Get fixed-length tag value for an Range of entities
     *\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called
     *      for variable-length tag.
     */
  ErrorCode get_data( TagId tag_handle,
                        const Range& handles,
                        void* data,
                        const void* default_value ) const;
  
    /** Get pointers to tag data for an Range of entities.
     *
     *\param tag_id   The tag.
     *\param handles  The entities.
     *\param values   Array of pointers of type 'const void*'.  Array
     *                must be the same length as the size of 'entities'.
     *                The array will be populated with pointers to the
     *                internal storage of the tag data for each entity.
     *\param lengths  Array of integers.  Will be populated with the 
     *                length of the tag value for each entity.  Argument
     *                is optional for fixed-length tags.
     *\param default_value The default value for the tag.
     *\param default_value_length The length of the default tag value.
     */
  ErrorCode get_data( TagId tag_handle,
                        const Range& handles,
                        const void** data,
                        int* lengths,
                        const void* default_value,
                        int default_value_length ) const;

  //! removes data
  ErrorCode remove_data(const TagId tag_handle, const EntityHandle entity_handle);

  //! gets all entity handles that match a type and tag
  ErrorCode get_entities(const TagId tag_handle, Range &entities);

  //! gets all entity handles that match a tag
  ErrorCode get_entities(const TagId tag_handle, const EntityType type,
                            Range &entities);

  //! gets all entity handles that match a tag
  ErrorCode get_entities(const Range &range,
                            const TagId tag_handle, const EntityType type,
                            Range &entities);

  //! gets all tags on a given entity handle
  ErrorCode get_tags(const EntityHandle entity,
                       std::vector<Tag> &all_tags);
  
  //! gets all entity handles that match a tag
  ErrorCode get_entities_with_tag_value( const TagId tag_handle, 
                                           const TagInfo& tag_info,
                                           const EntityType type,
                                           Range &entities,
                                           const void* tag_value,
                                           int value_size);

  //! gets all entity handles that match a tag
  ErrorCode get_entities_with_tag_value( const Range &range,
                                           const TagId tag_handle, 
                                           const TagInfo& tag_info,
                                           const EntityType type,
                                           Range &entities,
                                           const void* tag_value,
                                           int value_size);

  //! gets the number of entities that match a tag
  ErrorCode get_number_entities(const TagId tag_handle, const EntityType type, int& num_ent);


  //! gets the number of entities that match a tag
  ErrorCode get_number_entities(const Range &range,
                                   const TagId tag_handle, const EntityType type, int& num_ent);

  ErrorCode get_memory_use( TagId tag_id, 
                              unsigned long& total,
                              unsigned long& per_entity );

private:

  std::vector<SparseTagCollection*> mDataTags;
};

} // namespace moab

#endif //SPARSE_TAG_SUPER_COLLECTION_HPP




