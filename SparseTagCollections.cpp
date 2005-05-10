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


SparseTagSuperCollection::~SparseTagSuperCollection()
{
  SparseTagCollection* tag_collection = NULL;

  std::map<MBTagId, SparseTagCollection*>::iterator tag_iterator;
  for(tag_iterator = mDataTags.begin(); tag_iterator != mDataTags.end(); ++tag_iterator)
  {
    tag_collection = (*tag_iterator).second;
    if(tag_collection != NULL)
      delete tag_collection;
  }
  mDataTags.clear();
}

void SparseTagSuperCollection::reset_data()
{
  std::map<MBTagId, SparseTagCollection*>::iterator tag_iterator;
  for(tag_iterator = mDataTags.begin(); tag_iterator != mDataTags.end(); ++tag_iterator)
  {
    int data_size = tag_iterator->second->tag_size();
    delete tag_iterator->second;
    tag_iterator->second = new SparseTagCollection(data_size);
  }
  
}

MBErrorCode SparseTagSuperCollection::reserve_tag_id(int data_size, MBTagId& tag_id)
{
  if(data_size<=0)
    return MB_FAILURE;

  // start at 1
  tag_id = 1;
 
  // go until we find one that isn't used 
  for( std::map<MBTagId, SparseTagCollection*>::iterator iter = mDataTags.begin();
       iter != mDataTags.end(); ++iter)
  {
    if(tag_id == iter->first)
      tag_id++;
    else
      break;
  }

  SparseTagCollection* new_tag_collection = new SparseTagCollection(data_size);
  std::pair<MBTagId,SparseTagCollection*> tmp(tag_id, new_tag_collection);
  mDataTags.insert( tmp );

  return MB_SUCCESS;
}

MBErrorCode SparseTagSuperCollection::release_tag_id(MBTagId tag_id)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_id);

  if(iterator != mDataTags.end())
  {
    delete iterator->second;
    mDataTags.erase(iterator);
    return MB_SUCCESS;
  }
  return MB_TAG_NOT_FOUND;
}

int SparseTagSuperCollection::tag_size(const MBTagId tag_id) const
{
  std::map<MBTagId, SparseTagCollection*>::const_iterator iterator =
    mDataTags.find(tag_id);
  if(iterator != mDataTags.end())
  {
    return iterator->second->tag_size();
  }
  return 0;
}

MBErrorCode SparseTagSuperCollection::set_data(const MBTagId tag_handle, 
    const MBEntityHandle entity_handle, const void* data)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator != mDataTags.end())
    return iterator->second->set_data(entity_handle, data);
  else
    return MB_TAG_NOT_FOUND;
}

MBErrorCode SparseTagSuperCollection::get_data(const MBTagId tag_handle, const MBEntityHandle entity_handle, void* data)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator != mDataTags.end())
    return iterator->second->get_data(entity_handle, data);
  else
    return MB_TAG_NOT_FOUND;
}


MBErrorCode SparseTagSuperCollection::remove_data( const MBTagId tag_handle, const MBEntityHandle entity_handle )
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator != mDataTags.end())
    return iterator->second->remove_data(entity_handle);

  return MB_TAG_NOT_FOUND;
}


MBEntityHandle SparseTagSuperCollection::find_entity( const MBTagId tag_handle, const void* data )
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator != mDataTags.end())
    return iterator->second->find_entity(data);
  else
    return 0;

}


//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities(const MBTagId tag_handle, const MBEntityType type,
                                                    MBRange &entities)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator != mDataTags.end())
    return iterator->second->get_entities(type, entities);

  return MB_TAG_NOT_FOUND;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities(const MBRange &range,
                                                    const MBTagId tag_handle, const MBEntityType type,
                                                    MBRange &entities)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  MBErrorCode result = MB_TAG_NOT_FOUND;
  MBRange dum_range;
  
  if(iterator != mDataTags.end())
    result = iterator->second->get_entities(type, dum_range);

  std::set_intersection(dum_range.begin(), dum_range.end(),
                        range.begin(), range.end(),
                        mb_range_inserter(entities));
  
  return result;
}

MBErrorCode SparseTagSuperCollection::get_tags(const MBEntityHandle entity,
                                                std::vector<MBTag> &all_tags)
{
  std::map<MBTagId, SparseTagCollection*>::iterator tag_it;
  for (tag_it = mDataTags.begin(); tag_it != mDataTags.end(); tag_it++) {
    if (tag_it->second->contains(entity))
      all_tags.push_back(TAG_HANDLE_FROM_ID(tag_it->first, MB_TAG_SPARSE));
  }
  
  return MB_SUCCESS;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities_with_tag_value(const MBTagId tag_handle, const MBEntityType type,
                                                                   MBRange &entities, const void* tag_value)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator != mDataTags.end())
    return iterator->second->get_entities_with_tag_value(type, entities, tag_value);

  return MB_TAG_NOT_FOUND;
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_entities_with_tag_value(const MBRange &range,
                                                                   const MBTagId tag_handle, const MBEntityType type,
                                                                   MBRange &entities, const void* tag_value)
{
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  MBRange dum_range;
  MBErrorCode result = MB_TAG_NOT_FOUND;

  if(iterator != mDataTags.end())
    result = iterator->second->get_entities_with_tag_value(type, dum_range, tag_value);

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
  std::map<MBTagId, SparseTagCollection*>::iterator iterator =
    mDataTags.find(tag_handle);

  if(iterator == mDataTags.end())
    return MB_TAG_NOT_FOUND;

  return iterator->second->get_number_entities(type, num_entities);
}

//! gets all entity handles that match a type and tag
MBErrorCode SparseTagSuperCollection::get_number_entities(const MBRange &range,
                                                           const MBTagId tag_handle, const MBEntityType type,
                                                           int& num_entities)
{
  MBRange dum_range;
  MBErrorCode result = get_entities(range, tag_handle, type, dum_range);
  if (MB_SUCCESS != result) return result;
  
  num_entities = dum_range.size();
  return result;
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
  if(iterator->first == entity_handle)
  {
    memcpy( iterator->second, data, mDataSize);
    ret_val = MB_SUCCESS;
  }
  // we need to make some data space
  else 
  {
    void* new_data = mAllocator.allocate(mDataSize);
    memcpy(new_data, data, mDataSize);
    mData.insert( iterator, std::make_pair(entity_handle, new_data));
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
    if(TYPE_FROM_HANDLE(iter->first) == type)
      if( memcmp( iter->second, tag_value, mDataSize ) == 0) 
        entities.insert(iter->first);    
  }

  return MB_SUCCESS;
}




