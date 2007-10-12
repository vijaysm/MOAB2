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
MBErrorCode DenseTagSuperCollection::reserve_tag_id(int num_bytes, const void* default_data, MBTagId tag_id)
{
  // make sure we get a good number of bytes
  if(num_bytes <= 0 )
    return MB_FAILURE;
  
  // make sure we have storage for tag id
  if (mDensePageGroups.size() <= tag_id)
    mDensePageGroups.resize( tag_id+1, 0 );
  
  // make sure tag_id isn't already in use
  if (mDensePageGroups[tag_id])
    return MB_FAILURE;
    
  // allocate tag data
  mDensePageGroups[tag_id] = new DensePageGroup( num_bytes, default_data ) ;

  return MB_SUCCESS;
  
}


/*! give back a tag id that was used to set and get data */
MBErrorCode DenseTagSuperCollection::release_tag_id(MBTagId tag_id)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  delete *group;
  *group = 0;
  
  // clean up a bit if this is the last one
  while (!mDensePageGroups.back())
    mDensePageGroups.pop_back();
  
  return MB_SUCCESS;
}
  
void DenseTagSuperCollection::reset_data()
{
  for (std::vector<DensePageGroup*>::iterator iter = mDensePageGroups.begin();
      iter != mDensePageGroups.end();
      ++iter)
  {
    if(*iter == NULL)
      continue;
    DensePageGroup* newgroup = new DensePageGroup( (*iter)->tag_size(), (*iter)->get_default_value() );
    delete *iter;
    *iter = newgroup;
  }
}

// get the entities with a tag
MBErrorCode DenseTagSuperCollection::get_entities(const MBTagId tag_id, const MBEntityType type, MBRange& entities)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->get_entities(type, entities);
}

// get the entities with a tag
MBErrorCode DenseTagSuperCollection::get_entities(const MBRange &range,
                                                   const MBTagId tag_id, const MBEntityType type, MBRange& entities)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  MBRange dum_range;
  MBErrorCode result = (*group)->get_entities(type, dum_range);
  if (MB_SUCCESS != result) return result;

  std::set_intersection(dum_range.begin(), dum_range.end(),
                        range.begin(), range.end(),
                        mb_range_inserter(entities));
  
  return result;
}

// get number of entities with a tag
MBErrorCode DenseTagSuperCollection::get_number_entities(const MBTagId tag_id, const MBEntityType type, int& entities)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->get_number_entities(type, entities);
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
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->get_entities_with_tag_value(type, value, entities);
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

MBErrorCode DenseTagSuperCollection::get_memory_use( MBTagId tag_id,
                                             unsigned long& total,
                                             unsigned long& per_entity )
{
  std::vector<DensePageGroup*>::iterator group;
  
    // get memory use by dense page group
  group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;
  (*group)->get_memory_use( total, per_entity );
  
    // count number of occupied slots in mDensePageGroups
  unsigned num_used = 0;
  for (group = mDensePageGroups.begin(); group != mDensePageGroups.end(); ++group)
    if (*group)
      ++num_used;
  
    // add in amortized storage in mDensePageGroups vector
  total += sizeof(DensePageGroup*) * mDensePageGroups.capacity() / num_used;

  return MB_SUCCESS;
}

MBErrorCode DensePageGroup::get_entities_with_tag_value(const MBEntityType type,
                                                        const void* value, 
                                                        MBRange &entities)
{
    // for now, return if default value is requested
  if (mDefaultValue && !memcmp(value,mDefaultValue,mBytesPerFlag))
    return MB_FAILURE;

    // iterate over dense pages
  std::vector<DensePage>::iterator page_it;
  const std::vector<DensePage>::iterator end = mDensePages[type].end();
  int dum =0;
  MBEntityHandle handle = CREATE_HANDLE(type, MB_START_ID, dum);
  MBRange::iterator insert_iter = entities.begin();
  for(page_it = mDensePages[type].begin(); page_it != end; 
      ++page_it, handle += DensePage::mPageSize)
  {
    if (page_it->has_data()) {
      for (int i = 1; i <= DensePage::mPageSize; i++) {
        if (!page_it->memcmp(i, mBytesPerFlag, value))
          entities.insert(handle+i-1);
      }
    }
  }

  return MB_SUCCESS;
}

MBErrorCode DensePageGroup::get_memory_use( unsigned long& total,
                                            unsigned long& per_entity )
{
  per_entity = tag_size();
  
  total = sizeof(*this);
  if (mDefaultValue)
    total += tag_size();
  
  for (unsigned i = 0; i < MBMAXTYPE; ++i) {
    total += mDensePages[i].capacity() * sizeof(DensePage);
    for (unsigned long j = 0; j < mDensePages[i].size(); ++j)
      if (mDensePages[i][j].has_data())
        total += DensePage::mPageSize * tag_size();
  }
  
  return MB_SUCCESS;
}
