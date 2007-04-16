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
 * Filename   :     SparseTagCollections.cpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/


#include <memory.h>

#include "SparseTagCollections.hpp"
#include "MBRange.hpp"

/*
  SparseTagSuperCollection functions -----------------------------
*/

#define get_collection( A ) ((A) < mDataTags.size() ? mDataTags[(A)] : 0)

SparseTagSuperCollection::~SparseTagSuperCollection()
{
  std::vector<SparseTagCollection*>::iterator tag_iterator;
  for(tag_iterator = mDataTags.begin(); tag_iterator != mDataTags.end(); ++tag_iterator)
    delete *tag_iterator;
  mDataTags.clear();
}

void SparseTagSuperCollection::reset_data()
{
  std::vector<SparseTagCollection*>::iterator tag_iterator;
  for(tag_iterator = mDataTags.begin(); tag_iterator != mDataTags.end(); ++tag_iterator)
  {
    if (*tag_iterator) {
      int data_size = (*tag_iterator)->tag_size();
      delete *tag_iterator;
      *tag_iterator = new SparseTagCollection(data_size);
    }
  }
  
}

MBErrorCode SparseTagSuperCollection::reserve_tag_id(int data_size, MBTagId& tag_id)
{
  if(data_size<=0)
    return MB_FAILURE;

  // start at 1
  tag_id = 1;
 
  while (tag_id < mDataTags.size() && mDataTags[tag_id])
    ++tag_id;
  
  if (tag_id >= mDataTags.size())
    mDataTags.resize( tag_id+1, 0 );
    
  mDataTags[tag_id] = new SparseTagCollection(data_size);
  return MB_SUCCESS;
}

MBErrorCode SparseTagSuperCollection::release_tag_id(MBTagId tag_id)
{
  if (tag_id >= mDataTags.size() || !mDataTags[tag_id])
    return MB_TAG_NOT_FOUND;
  
  delete mDataTags[tag_id];
  mDataTags[tag_id] = 0;
  return MB_SUCCESS;
}

int SparseTagSuperCollection::tag_size(const MBTagId tag_id) const
{
  SparseTagCollection* coll = get_collection(tag_id);
  return coll ? coll->tag_size() : 0;
}

MBErrorCode SparseTagSuperCollection::set_data(const MBTagId tag_handle, 
    const MBEntityHandle entity_handle, const void* data)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  return coll ? coll->set_data( entity_handle, data ) : MB_TAG_NOT_FOUND;
}

MBErrorCode SparseTagSuperCollection::get_data(const MBTagId tag_handle, 
    const MBEntityHandle entity_handle, void* data)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  return coll ? coll->get_data( entity_handle, data ) : MB_TAG_NOT_FOUND;
}


MBErrorCode SparseTagSuperCollection::remove_data( const MBTagId tag_handle, 
    const MBEntityHandle entity_handle )
{
  SparseTagCollection* coll = get_collection(tag_handle);
  return coll ? coll->remove_data( entity_handle ) : MB_TAG_NOT_FOUND;
}


MBEntityHandle SparseTagSuperCollection::find_entity( const MBTagId tag_handle, const void* data )
{
  SparseTagCollection* coll = get_collection(tag_handle);
  return coll ? coll->find_entity( data ) : MB_TAG_NOT_FOUND;
}


//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities(const MBTagId tag_handle, const MBEntityType type,
                                                    MBRange &entities)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  return coll ? coll->get_entities(type, entities) : MB_TAG_NOT_FOUND;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities(const MBRange &range,
                                                    const MBTagId tag_handle, const MBEntityType type,
                                                    MBRange &entities)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  if (!coll)
    return MB_TAG_NOT_FOUND;
  
  MBRange dum_range;
  MBErrorCode result = coll->get_entities(type, dum_range);

  std::set_intersection(dum_range.begin(), dum_range.end(),
                        range.begin(), range.end(),
                        mb_range_inserter(entities));
  
  return result;
}

MBErrorCode SparseTagSuperCollection::get_tags(const MBEntityHandle entity,
                                                std::vector<MBTag> &all_tags)
{
  for (MBTagId id = 0; id < mDataTags.size(); ++id)
    if (mDataTags[id] && mDataTags[id]->contains(entity))
      all_tags.push_back( TAG_HANDLE_FROM_ID( id, MB_TAG_SPARSE ) );
  
  return MB_SUCCESS;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities_with_tag_value(const MBTagId tag_handle, const MBEntityType type,
                                                                   MBRange &entities, const void* tag_value)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  if (!coll)
    return MB_TAG_NOT_FOUND;
  
  return coll->get_entities_with_tag_value(type, entities, tag_value);
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities_with_tag_value(const MBRange &range,
                                                                   const MBTagId tag_handle, const MBEntityType type,
                                                                   MBRange &entities, const void* tag_value)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  if (!coll)
    return MB_TAG_NOT_FOUND;

  MBRange dum_range;
  MBErrorCode result = coll->get_entities_with_tag_value(type, dum_range, tag_value);

    // do this the hard way to preserve order in the vector
  std::set_intersection(range.begin(), range.end(),
                        dum_range.begin(), dum_range.end(),
                        mb_range_inserter(entities));
  
  return result;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_number_entities(const MBTagId tag_handle, const MBEntityType type,
                                                           int& num_entities)
{
  SparseTagCollection* coll = get_collection(tag_handle);
  return coll ? coll->get_number_entities( type, num_entities ) : MB_TAG_NOT_FOUND;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_number_entities(const MBRange &range,
                                                           const MBTagId tag_handle, const MBEntityType type,
                                                           int& num_entities)
{
  MBRange dum_range;
  MBErrorCode result = get_entities(range, tag_handle, type, dum_range);
  num_entities = dum_range.size();
  return result;
}


MBErrorCode SparseTagSuperCollection::get_memory_use( MBTagId tag_id,
                                              unsigned long& total,
                                              unsigned long& per_entity )
{
  SparseTagCollection* coll = get_collection(tag_id);
  if (!coll)
    return MB_TAG_NOT_FOUND;

    // 3*sizeof(void*)                      - std::map RB tree node
    // sizeof(void*)*sizeof(MBEntityHandle) - data in std::map node
    // coll->tag_size()                     - the actual tag data
  per_entity = 4*sizeof(void*)+sizeof(MBEntityHandle)+coll->tag_size();
    
    // Count number of occupied slots in mDataTags vector
  unsigned num_coll =0;
  for (unsigned i = 0; i < mDataTags.size(); ++i)
    if (mDataTags[i])
      ++num_coll;

    // amortized space in mDataTags vector
  total = sizeof(SparseTagCollection*) * mDataTags.capacity() / num_coll;
    // SparseTagCollection object for this tag
  total += sizeof(SparseTagCollection);
    // Per-entity data in SparseTagCollection
  total += per_entity * coll->get_number_entities();
  return MB_SUCCESS;
}

/*
  SparseTagCollection functions ----------------------------------
*/

SparseTagCollection::SparseTagCollection(int data_size)
{
  mDataSize = data_size;
}

SparseTagCollection::~SparseTagCollection()
{
  void* tag_data = NULL;

  std::map<MBEntityHandle, void*>::iterator tag_iterator;
  for(tag_iterator = mData.begin(); tag_iterator != mData.end(); ++tag_iterator)
  {
    tag_data = tag_iterator->second;
    if(tag_data != NULL)
      mAllocator.destroy(tag_data);
  }
  mData.clear();
}

MBErrorCode SparseTagCollection::set_data(const MBEntityHandle entity_handle, const void* data)
{
  MBErrorCode ret_val = MB_TAG_NOT_FOUND;

  std::map<MBEntityHandle, void*>::iterator iterator =
    mData.lower_bound(entity_handle);
  
  // data space already exists
  if (iterator!= mData.end() && iterator->first == entity_handle)
  {
    memcpy( iterator->second, data, mDataSize);
    ret_val = MB_SUCCESS;
  }
  // we need to make some data space
  else 
  {
    void* new_data = mAllocator.allocate(mDataSize);
    memcpy(new_data, data, mDataSize);
    mData.insert(iterator, std::pair<const MBEntityHandle,void*>(entity_handle, new_data));
    ret_val = MB_SUCCESS;
  }

  return ret_val;
}

MBErrorCode SparseTagCollection::get_data(const MBEntityHandle entity_handle, void* data)
{
  std::map<MBEntityHandle, void*>::iterator iter =
    mData.find(entity_handle);

  if(iter == mData.end())
    return MB_TAG_NOT_FOUND;
  
  memcpy(data, iter->second, mDataSize);
  return MB_SUCCESS;
}



MBErrorCode SparseTagCollection::remove_data( const MBEntityHandle entity_handle )
{
  std::map<MBEntityHandle, void*>::iterator iterator =
    mData.find(entity_handle);

  if(iterator != mData.end())
  {
    mAllocator.destroy(iterator->second);
    mData.erase(iterator);
    return MB_SUCCESS;
  }
  return MB_ENTITY_NOT_FOUND;
}



MBEntityHandle SparseTagCollection::find_entity( const void* data )
{
  std::map<MBEntityHandle, void*>::iterator iterator;
  for(iterator = mData.begin();
      iterator != mData.end();
      iterator++)
  {
    if(memcmp(data, iterator->second, mDataSize) == 0)
      return iterator->first;
  }

  return 0;
}

//! get number of entities of type
MBErrorCode SparseTagCollection::get_number_entities(MBEntityType type, int& num_entities)
{
  std::map<MBEntityHandle, void*>::iterator iter;
  for(iter = mData.begin(); iter != mData.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(iter->first) == type)
      num_entities++;
  }
  return MB_SUCCESS;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagCollection::get_entities(MBEntityType type, MBRange &entities)
{
  std::map<MBEntityHandle, void*>::iterator iter;
  for(iter = mData.begin(); iter != mData.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(iter->first) == type)
      entities.insert(iter->first);    
  }
  return MB_SUCCESS;
}


//! gets all entity handles that match a type, tag and tag value
MBErrorCode SparseTagCollection::get_entities_with_tag_value(MBEntityType type, 
                                                              MBRange &entities, const void* tag_value)
{
  std::map<MBEntityHandle, void*>::iterator iter;

  for(iter = mData.begin(); iter != mData.end(); ++iter)
  {
    if(TYPE_FROM_HANDLE(iter->first) == type) {
#ifndef NDEBUG
      MBEntityHandle this_ent = iter->first;
      void *this_tag = iter->second;
#endif
      if( memcmp( iter->second, tag_value, mDataSize ) == 0) 
        entities.insert(iter->first);    
    }
  }

  return MB_SUCCESS;
}




