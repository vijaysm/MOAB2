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
 * Filename   :     SparseTagCollections.hpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/



#ifndef SPARSE_TAG_COLLECTIONS_HPP
#define SPARSE_TAG_COLLECTIONS_HPP

#ifndef IS_BUILDING_MB
#error "SparseTagCollections.hpp isn't supposed to be included into an application"
#endif

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include <map>
#include <vector>

#include "MBTypes.h"
#include "MBInternals.hpp"
class MBRange;

//! allocator for tag data
class SparseTagDataAllocator
{
public:
  //! constructor
  SparseTagDataAllocator(){}
  //! destructor
  ~SparseTagDataAllocator(){}
  //! allocates memory of size and returns pointer
  void* allocate(size_t data_size) { return malloc(data_size); }
  //! frees the memory
  void destroy(void* p){ free(p); }
};


//! collection of tag data associated with entity ids
class SparseTagCollection
{
public:
  
  //! constructor takes tag data size
  SparseTagCollection(int data_size);
  
  //! destructor
  ~SparseTagCollection();
  
  //! set the tag data for an entity id
  MBErrorCode set_data(const MBEntityHandle entity_handle, const void* data);

  //! get the tag data for an entity id
  MBErrorCode get_data(const MBEntityHandle entity_handle, void* data);

  //! removes the data
  MBErrorCode remove_data(const MBEntityHandle entity_handle);
  
  //! finds the first entity handle that matches this data
  MBEntityHandle find_entity( const void* data );

  //! get number of entities of type
  MBErrorCode get_number_entities(MBEntityType type, int& num_entities);
  
  //! get number of entities
  unsigned long get_number_entities()
    { return mData.size(); }

  //! gets all entity handles that match a type and tag
  MBErrorCode get_entities(MBEntityType type, MBRange &entities);

  //! gets all entity handles that match a type, tag, tag_value
  MBErrorCode get_entities_with_tag_value(MBEntityType type, 
                                           MBRange &entities, 
                                           const void* tag_value);

  //! if this collection contains this entity, return true, otherwise false
  bool contains(const MBEntityHandle entity) const;
  
  int tag_size() const { return mDataSize; }

protected:
  
  //! hidden constructor
  SparseTagCollection(){}
  
  //! size of the data
  int mDataSize;

  //! allocator for this collection
  SparseTagDataAllocator mAllocator;

  //! map of entity id and tag data
  std::map<MBEntityHandle /*entity_handle*/ , void* /*data*/ > mData;

};

inline bool SparseTagCollection::contains(const MBEntityHandle entity) const
{
  return (mData.find(entity) == mData.end() ? false : true);
}

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
  MBErrorCode reserve_tag_id(int data_size, MBTagId tag_id);

  //! releases an MBTagId
  MBErrorCode release_tag_id(MBTagId tag_id);

  //! size of data tag
  int tag_size(const MBTagId tag_id) const;

  //! set the data of a tag
  MBErrorCode set_data(const MBTagId tag_handle, const MBEntityHandle entity_handle, const void* data);

  //! get the data of a tag
  MBErrorCode get_data(const MBTagId tag_handle, const MBEntityHandle entity_handle, void* data);

  //! removes data
  MBErrorCode remove_data(const MBTagId tag_handle, const MBEntityHandle entity_handle);
  
  //! finds the entity handle with this data
  MBEntityHandle find_entity( const MBTagId tag_handle, const void* data );

  //! gets all entity handles that match a tag
  MBErrorCode get_entities(const MBTagId tag_handle, const MBEntityType type,
                            MBRange &entities);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities(const MBRange &range,
                            const MBTagId tag_handle, const MBEntityType type,
                            MBRange &entities);

  //! gets all tags on a given entity handle
  MBErrorCode get_tags(const MBEntityHandle entity,
                       std::vector<MBTag> &all_tags);
  
  //! gets all entity handles that match a tag
  MBErrorCode get_entities_with_tag_value(const MBTagId tag_handle, const MBEntityType type,
                                           MBRange &entities,
                                           const void* tag_value);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities_with_tag_value(const MBRange &range,
                                           const MBTagId tag_handle, const MBEntityType type,
                                           MBRange &entities,
                                           const void* tag_value);

  //! gets the number of entities that match a tag
  MBErrorCode get_number_entities(const MBTagId tag_handle, const MBEntityType type, int& num_ent);


  //! gets the number of entities that match a tag
  MBErrorCode get_number_entities(const MBRange &range,
                                   const MBTagId tag_handle, const MBEntityType type, int& num_ent);

  MBErrorCode get_memory_use( MBTagId tag_id, 
                              unsigned long& total,
                              unsigned long& per_entity );

private:

  std::vector<SparseTagCollection*> mDataTags;
};


#endif //SPARSE_TAG_COLLECTIONS_HPP




