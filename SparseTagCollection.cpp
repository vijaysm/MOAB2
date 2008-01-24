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
 * Filename   :     SparseTagCollection.cpp
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
#include <algorithm>

#include "SparseTagCollection.hpp"
#include "MBRange.hpp"
#include "TagCompare.hpp"

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
    if (mDataSize == MB_VARIABLE_LENGTH)
      reinterpret_cast<VarLenTag*>(tag_data)->clear();
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
    if (mDataSize == MB_VARIABLE_LENGTH)
      reinterpret_cast<VarLenTag*>(iterator->second)->clear();
    mAllocator.destroy(iterator->second);
    mData.erase(iterator);
    return MB_SUCCESS;
  }
  return MB_ENTITY_NOT_FOUND;
}


//! get number of entities of type
MBErrorCode SparseTagCollection::get_number_entities(MBEntityType type, int& num_entities)
{
  num_entities = 0;
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
MBErrorCode SparseTagCollection::get_entities_with_tag_value(
                                                    const TagInfo& tag_info,
                                                    MBEntityType type, 
                                                    MBRange &entities, 
                                                    const void* tag_value,
                                                    int value_size)
{
  std::map<MBEntityHandle, void*>::iterator iter, end;
  int junk;
  iter = mData.lower_bound( CREATE_HANDLE( type, MB_START_ID, junk ) );
  end = mData.upper_bound( CREATE_HANDLE( type, MB_END_ID, junk ) );
  find_tag_values_equal( tag_info, tag_value, value_size, iter, end, entities );
  return MB_SUCCESS;
}





