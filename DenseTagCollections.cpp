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



/*  Dense tag for MB
 *
 *  File   :      DenseTagCollections.cpp
 *  Creator:      Clinton Stimpson
 *  Date   :      10-28-2002
 */

#include "DenseTagCollections.hpp"

#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

/*! page size of 1024 bytes for storage */
const int DensePage::mPageSize = 1024;

/*! returns an available tag id to use when getting and setting data */
MBErrorCode DenseTagSuperCollection::reserve_tag_id(int num_bytes, const void* default_data, MBTagId& tag_id)
{
  // make sure we get a good number of bytes
  if(num_bytes <= 0 )
  {
    tag_id = 0;
    return MB_FAILURE;
  }

  // see if we have any page groups that aren't being used
  int which_one = -1;
  for(std::vector<DensePageGroup*>::iterator iter = (*mDensePageGroups).begin();
      (iter != (*mDensePageGroups).end()) && (which_one == -1);
      ++iter)
  {
    if(*iter == NULL)
    {
      which_one = iter - (*mDensePageGroups).begin();
    }
  }
  
  // if they are all being used, make a new one
  if(which_one == -1)
  {
    for(int i=0; i<(int)MBMAXTYPE; i++)
      mDensePageGroups[i].push_back(NULL);
    mDefaultData.push_back(NULL);
    which_one = (*mDensePageGroups).size() - 1;
  }

  for(int i=0; i<(int)MBMAXTYPE; i++)
    mDensePageGroups[i][which_one] = new DensePageGroup(num_bytes);

  if(default_data)
  {
    mDefaultData[which_one] = malloc(num_bytes);
    memcpy(mDefaultData[which_one], default_data, num_bytes);
  }

  tag_id = which_one;

  mDensePageGroupsSize = (*mDensePageGroups).size();

  // can we really fail?
  return MB_SUCCESS;
  
}


/*! give back a tag id that was used to set and get data */
MBErrorCode DenseTagSuperCollection::release_tag_id(MBTagId tag_id)
{
  // make sure tag_id is good
  if(tag_id >= mDensePageGroupsSize)
    return MB_TAG_NOT_FOUND;

  // mark it unused
  for(int i=0; i<(int)MBMAXTYPE; i++)
  {
    delete mDensePageGroups[i][tag_id];
    mDensePageGroups[i][tag_id] = NULL;
    if(mDefaultData[tag_id])
    {
      free(mDefaultData[tag_id]);
      mDefaultData[tag_id] = NULL;
    }
  }

  // clean up a bit if this is the last one
  if((tag_id+1) == (unsigned short)(*mDensePageGroups).size())
  {
    for(int i=0; i<(int)MBMAXTYPE; i++)
    {
      mDensePageGroups[i].resize(mDensePageGroups[i].size()-1);
    }
    mDefaultData.resize(mDefaultData.size()-1);
  }
  
  mDensePageGroupsSize = (*mDensePageGroups).size();
  
  return MB_SUCCESS;
}
  
void DenseTagSuperCollection::reset_data()
{
  for(int i = 0; i<(int)MBMAXTYPE; i++)
  {
    for (std::vector<DensePageGroup*>::iterator iter = mDensePageGroups[i].begin();
        iter != mDensePageGroups[i].end();
        ++iter)
    {
      if(*iter == NULL)
        continue;
      int tag_size = (*iter)->tag_size();
      delete *iter;
      *iter = new DensePageGroup(tag_size);
    }
  }
}

// get the entities with a tag
MBErrorCode DenseTagSuperCollection::get_entities(const MBTagId tag_id, const MBEntityType type, MBRange& entities)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  return mDensePageGroups[type][tag_id]->get_entities(type, entities);
}

// get the entities with a tag
MBErrorCode DenseTagSuperCollection::get_entities(const MBRange &range,
                                                   const MBTagId tag_id, const MBEntityType type, MBRange& entities)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  MBRange dum_range;
  MBErrorCode result = mDensePageGroups[type][tag_id]->get_entities(type, dum_range);
  if (MB_SUCCESS != result) return result;

  std::set_intersection(dum_range.begin(), dum_range.end(),
                        range.begin(), range.end(),
                        mb_range_inserter(entities));
  
  return result;
}

// get number of entities with a tag
MBErrorCode DenseTagSuperCollection::get_number_entities(const MBTagId tag_id, const MBEntityType type, int& entities)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  return mDensePageGroups[type][tag_id]->get_number_entities(type, entities);
}

// get number of entities with a tag
MBErrorCode DenseTagSuperCollection::get_number_entities(const MBRange &range,
                                                          const MBTagId tag_id, const MBEntityType type, int& entities)
{
  MBRange dum_range;
  MBErrorCode result = get_entities(range, tag_id, type, dum_range);
  if (MB_SUCCESS != result) return result;
  entities = dum_range.size();

  return result;
}

//! get the entities with a value
MBErrorCode DenseTagSuperCollection::get_entities_with_tag_value(const MBTagId tag_id, const MBEntityType type, 
                                                                  MBRange &entities, const void* value)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_TAG_NOT_FOUND;
  
  bool equals_default = (NULL != mDefaultData[tag_id] ? 
    (memcmp(mDefaultData[tag_id], value, (*mDensePageGroups)[tag_id]->tag_size()) == 0) :
                         false);

    // for now, return if default value is requested
  if (equals_default)
    return MB_FAILURE;

  return mDensePageGroups[type][tag_id]->get_entities_with_tag_value(type, value, entities, 
                                                                     equals_default);
}


//! get the entities with a value
MBErrorCode DenseTagSuperCollection::get_entities_with_tag_value(const MBRange &range,
                                                                  const MBTagId tag_id,
                                                                  const MBEntityType type,
                                                                  MBRange &entities, 
                                                                  const void* value)
{
  MBRange dum_range;
  MBErrorCode result = get_entities_with_tag_value(tag_id, type, dum_range, value);
  if (MB_SUCCESS != result) return result;
  
  std::set_intersection(range.begin(), range.end(),
                        dum_range.begin(), dum_range.end(),
                        mb_range_inserter(entities));
  return result;
}

MBErrorCode DensePageGroup::get_entities_with_tag_value(const MBEntityType type,
                                                        const void* value, 
                                                        MBRange &entities,
                                                        const bool equals_default)
{

  void* test_data = malloc(mBytesPerFlag);

  int dum = 0;
  MBEntityHandle handle = CREATE_HANDLE(type, MB_START_ID, 0, dum);
  MBEntityHandle end_handle = handle + mDensePages.size() * DensePage::mPageSize;
  MBRange::iterator insert_iter = entities.begin();
  for(; handle < end_handle; handle++)
  {
    MBErrorCode result = get_data(handle, test_data);
    if((result == MB_SUCCESS && !memcmp(test_data, value, mBytesPerFlag))
     ||(result == MB_TAG_NOT_FOUND && equals_default))
        insert_iter = entities.insert(insert_iter, handle, handle);
  }

  free(test_data);

  return MB_SUCCESS;
}

