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



#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif

#include <algorithm>
#include <assert.h>
#include "TagServer.hpp"
#include "MBRange.hpp"

using namespace std;

int TagInfo::size_from_data_type( MBDataType t )
{
  static const int sizes[] = { 1, 
                               sizeof(int), 
                               sizeof(double), 
                               1, 
                               sizeof(MBEntityHandle),
                               0 };  
   return sizes[t];
}

/*
  TagServer functions ----------------------------------
*/

TagServer::TagServer()
{
  //sanity checks
  assert((TAG_PROP_MASK | TAG_BIT_PROPERTIES[MB_TAG_BIT]) == TAG_PROP_MASK);
  assert((TAG_PROP_MASK | TAG_BIT_PROPERTIES[MB_TAG_SPARSE]) == TAG_PROP_MASK);
  assert((TAG_PROP_MASK | TAG_BIT_PROPERTIES[MB_TAG_DENSE]) == TAG_PROP_MASK);
  assert((TAG_PROP_MASK | TAG_BIT_PROPERTIES[MB_TAG_MESH]) == TAG_PROP_MASK);
  //assert((TAG_PROP_MASK | TAG_BIT_PROPERTIES[MB_TAG_STATIC) == TAG_PROP_MASK);

  // we need these tag properties to be in order.
  // if this order is changed, then change reset_data() as well
  assert(TAG_BIT_PROPERTIES[MB_TAG_BIT] < TAG_BIT_PROPERTIES[MB_TAG_SPARSE]);
  assert(TAG_BIT_PROPERTIES[MB_TAG_SPARSE] < TAG_BIT_PROPERTIES[MB_TAG_DENSE]);
  assert(TAG_BIT_PROPERTIES[MB_TAG_DENSE] < TAG_BIT_PROPERTIES[MB_TAG_MESH]);
  assert(TAG_BIT_PROPERTIES[MB_TAG_MESH] < TAG_BIT_PROPERTIES[MB_TAG_LAST]);
  
  mSparseData = new SparseTagSuperCollection;
  mDenseData = new DenseTagSuperCollection;
  mBitServer = new MBBitServer;
}

TagServer::~TagServer()
{
  if(mSparseData)
    delete mSparseData;
  if(mDenseData)
    delete mDenseData;
  if(mBitServer)
    delete mBitServer;
}


MBErrorCode TagServer::reset_all_data()
{
  mSparseData->reset_data();
  mDenseData->reset_data();
  mBitServer->reset_data();

  return MB_SUCCESS;
}


MBErrorCode TagServer::add_tag( const char *tag_name, 
                                const int data_size,
                                const MBTagType storage,
                                const MBDataType data_type,
                                MBTag &tag_handle,
                                const void *default_value)
{

  if(NULL != tag_name && strcmp(tag_name, "") != 0)
  {
    // verify that the name doesn't already exist for this entity type
    for(std::map<MBTag, TagInfo>::iterator tag_iterator = mTagTable.begin();
        tag_iterator != mTagTable.end();
        ++tag_iterator)
    {
      // if the type and name matches, another tag, return a "null" handle
      if(strcmp(tag_name, tag_iterator->second.get_name().c_str()) == 0)
      {
        tag_handle = tag_iterator->first;
        return MB_ALREADY_ALLOCATED;
      }
    }
  }
    
  MBErrorCode result;
  MBTagId id;

    // Input size must be a multiple of the size of the data type.
  int typesize = TagInfo::size_from_data_type( data_type );
  if (data_size % typesize)
    return MB_FAILURE;

  if(storage == MB_TAG_BIT)
  {
    if (data_type != MB_TYPE_BIT)
      return MB_FAILURE;
    result = mBitServer->reserve_tag_id(data_size, id);
  }
  else if(storage == MB_TAG_SPARSE || storage == MB_TAG_MESH)
  {
    result = mSparseData->reserve_tag_id(data_size, id);
  }
  else if(storage == MB_TAG_DENSE)
  {
    result = mDenseData->reserve_tag_id(data_size, default_value, id);
  }
  else
  {
    return MB_FAILURE;
  }
  
  if(result != MB_SUCCESS)
    return result;

  unsigned long tmp_handle = id;
  //tmp_handle |= (entity_type << 16);
  tmp_handle |= TAG_BIT_PROPERTIES[storage];
  tag_handle = reinterpret_cast<MBTag>(tmp_handle);

  // we have a valid id, lets register it
  if(tag_handle > 0)
  {
    TagInfo tag_info(tag_name, data_size, data_type, default_value);
    mTagTable.insert( std::pair<MBTag, TagInfo>( tag_handle, tag_info ) );
  }

  return MB_SUCCESS;
}


MBErrorCode TagServer::remove_tag(const MBTag tag_handle)
{

  const std::map<MBTag, TagInfo>::iterator iterator = mTagTable.find(tag_handle);
  
  if(iterator == mTagTable.end())
    return MB_TAG_NOT_FOUND;

  mTagTable.erase(iterator);
  MBErrorCode status = MB_TAG_NOT_FOUND;

  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    status = mDenseData->release_tag_id(ID_FROM_TAG_HANDLE(tag_handle));
  }
  else if((reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE]) ||
          (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_MESH]))
  {
    status = mSparseData->release_tag_id(ID_FROM_TAG_HANDLE(tag_handle));
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    status = mBitServer->release_tag_id(ID_FROM_TAG_HANDLE(tag_handle));
  }

  return status;

}


//! resets all data tagged to this entity handle back to default data or NULL.
//! this is used to clean out stale data that might be referenced again
MBErrorCode TagServer::reset_data(MBEntityHandle entity_handle)
{
  // note: this algorithm assumes that tag properties are
  // BITS < SPARSE < DENSE < STATIC.
  // if that is not the case anymore, then rearrange this algorithm

  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;

  std::map<MBTag, TagInfo>::iterator iter;

  // go through and clean out the bits
  MBTag max_tag = reinterpret_cast<MBTag>(TAG_BIT_PROPERTIES[MB_TAG_SPARSE]);
  for(iter = mTagTable.begin(); iter != mTagTable.end() && iter->first < max_tag; ++iter)
  {
    // default data for bits is zero
    mBitServer->weak_set_bits(ID_FROM_TAG_HANDLE(iter->first), entity_handle, 0);
  }

  // now clean out the sparse data
  max_tag = reinterpret_cast<MBTag>(TAG_BIT_PROPERTIES[MB_TAG_DENSE]);
  for(iter = mTagTable.begin(); iter != mTagTable.end() && iter->first < max_tag; ++iter)
  {
    mSparseData->remove_data(ID_FROM_TAG_HANDLE(iter->first), entity_handle);
  }
  
  // now clean out the dense data
  max_tag = reinterpret_cast<MBTag>(TAG_BIT_PROPERTIES[MB_TAG_LAST]);
  for(iter = mTagTable.begin(); iter != mTagTable.end() && iter->first < max_tag; ++iter)
  {
    mDenseData->remove_data(ID_FROM_TAG_HANDLE(iter->first), entity_handle,
        iter->second.default_value());
  }

  return MB_SUCCESS;
}


//! set the value of a tag
MBErrorCode TagServer::set_bits(const MBTag tag_handle, const MBEntityHandle entity_handle, unsigned char data )
{
  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;
  return mBitServer->set_bits(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
}

//! get the value of a tag
MBErrorCode TagServer::get_bits(const MBTag tag_handle, const MBEntityHandle entity_handle, unsigned char& data )
{
  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;
  return mBitServer->get_bits(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
}

MBErrorCode TagServer::set_data(const MBTag tag_handle, 
                                 const MBEntityHandle entity_handle, 
                                 const void* data)
{

  // this assumes that the entity_handle is valid, 
  // the check for valid handles are done one level up
  
  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;
  
  if( reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    return mDenseData->set_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    return mSparseData->set_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    return set_bits(tag_handle, entity_handle, *((unsigned char *)data));
  }
  else {
      // if we get here, we didn't find the right tag to set
    return MB_TAG_NOT_FOUND;
  }
  //else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_STATIC)
  //{
    //return mStaticSparseData.set_data(ID_FROM_TAG_HANDLE(tag_handle), data);
  //}
  
  return MB_SUCCESS;
}


MBErrorCode TagServer::set_data(const MBTag tag_handle, 
                       const MBEntityHandle* entity_handles, 
                       const int num_entities,
                       const void* data)
{

  // this assumes that the entity_handle is valid, 
  // the check for valid handles are done one level up
  
  MBErrorCode result = MB_SUCCESS;
  
  
  const MBTagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const unsigned char* mydata = static_cast<const unsigned char*>(data);
  const MBEntityHandle* end = entity_handles+num_entities;
  
  if( reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    const int data_size = mDenseData->tag_size(tag_id);
    for(const MBEntityHandle* iter = entity_handles; iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mDenseData->set_data(tag_id, *iter, mydata);
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    const int data_size = mSparseData->tag_size(tag_id);
    for(const MBEntityHandle* iter = entity_handles; iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mSparseData->set_data(tag_id, *iter, mydata);
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    if (num_entities == 1)
      return set_data(tag_handle, entity_handles[0], data);
    else
      // don't support this right now - not sure how to pass in multiple bit tags
      return MB_FAILURE;
  }
  else {
      // if we get here, we didn't find the right tag to set
    return MB_TAG_NOT_FOUND;
  }

  
  return MB_SUCCESS;
}

MBErrorCode TagServer::set_data(const MBTag tag_handle, 
                       const MBRange& entity_handles, 
                       const void* data)
{

  // this assumes that the entity_handle is valid, 
  // the check for valid handles are done one level up
  
  MBErrorCode result = MB_SUCCESS;
  
  const MBTagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const unsigned char* mydata = static_cast<const unsigned char*>(data);
  const MBRange::const_iterator end = entity_handles.end();
  
  if( reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    const int data_size = mDenseData->tag_size(tag_id);
    for(MBRange::const_iterator iter = entity_handles.begin(); iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mDenseData->set_data(tag_id, *iter, mydata);
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    const int data_size = mSparseData->tag_size(tag_id);
    for(MBRange::const_iterator iter = entity_handles.begin(); iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mSparseData->set_data(tag_id, *iter, mydata);
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    if (entity_handles.size() == 1)
      return set_data(tag_handle, *entity_handles.begin(), data);
    else
      // don't support this right now - not sure how to pass in multiple bit tags
      return MB_FAILURE;
  }
  else {
      // if we get here, we didn't find the right tag to set
    return MB_TAG_NOT_FOUND;
  }

  
  return MB_SUCCESS;
}


MBErrorCode TagServer::get_data(const MBTag tag_handle,
                       const MBEntityHandle entity_handle,
                       void* data)
{

  MBErrorCode result = MB_TAG_NOT_FOUND;

  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;

  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    result = mDenseData->get_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    result = mSparseData->get_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    result = get_bits(tag_handle, entity_handle, *((unsigned char *)data));
  }
  //else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_STATIC)
  //{
    //result = mStaticSparseData.get_data(ID_FROM_TAG_HANDLE(tag_handle), data);
  //}

  // if we couldn't get a value
  // try to get a default value
  if(result == MB_TAG_NOT_FOUND)
  {
    result = get_default_data(tag_handle, data);
      // if failure is returned, change it back to tag not found, since it's
      // ok to look for data and not find any
    if (result == MB_FAILURE) result = MB_TAG_NOT_FOUND;
  }

  return result;
}

MBErrorCode TagServer::get_data(const MBTag tag_handle,
                       const MBEntityHandle* entity_handles,
                       const int num_entities,
                       void* data)
{

  MBErrorCode result = MB_SUCCESS;

  const MBTagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  unsigned char* mydata = static_cast<unsigned char*>(data);
  const MBEntityHandle* end = entity_handles+num_entities;

  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    const int data_size = mDenseData->tag_size(tag_id);
    for(const MBEntityHandle* iter = entity_handles; iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mDenseData->get_data(tag_id, *iter, mydata);
      if(result == MB_TAG_NOT_FOUND)
      {
        result = get_default_data(tag_handle, mydata);
          // if failure is returned, change it back to tag not found, since it's
          // ok to look for data and not find any
        if (result == MB_FAILURE) result = MB_TAG_NOT_FOUND;
      }
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    const int data_size = mSparseData->tag_size(tag_id);
    for(const MBEntityHandle* iter = entity_handles; iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mSparseData->get_data(tag_id, *iter, mydata);
      if(result == MB_TAG_NOT_FOUND)
      {
        result = get_default_data(tag_handle, mydata);
          // if failure is returned, change it back to tag not found, since it's
          // ok to look for data and not find any
        if (result == MB_FAILURE) result = MB_TAG_NOT_FOUND;
      }
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    if (num_entities == 1)
      return get_data(tag_handle, entity_handles[0], data);
    else
      // don't support this right now - not sure how to pass in multiple bit tags
      return MB_FAILURE;
  }
  else {
      // if we get here, we didn't find the right tag to set
    return MB_TAG_NOT_FOUND;
  }


  return MB_SUCCESS;
}

MBErrorCode TagServer::get_data(const MBTag tag_handle,
                       const MBRange& entity_handles,
                       void* data)
{

  MBErrorCode result = MB_SUCCESS;

  const MBTagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  unsigned char* mydata = static_cast<unsigned char*>(data);
  const MBRange::const_iterator end = entity_handles.end();

  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    const int data_size = mDenseData->tag_size(tag_id);
    for(MBRange::const_iterator iter = entity_handles.begin(); iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mDenseData->get_data(tag_id, *iter, mydata);
      if(result == MB_TAG_NOT_FOUND)
      {
        result = get_default_data(tag_handle, mydata);
          // if failure is returned, change it back to tag not found, since it's
          // ok to look for data and not find any
        if (result == MB_FAILURE) result = MB_TAG_NOT_FOUND;
      }
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    const int data_size = mSparseData->tag_size(tag_id);
    for(MBRange::const_iterator iter = entity_handles.begin(); iter != end; ++iter)
    {
      if(TYPE_FROM_HANDLE(*iter) >= MBMAXTYPE)
        return MB_TYPE_OUT_OF_RANGE;
      result = mSparseData->get_data(tag_id, *iter, mydata);
      if(result == MB_TAG_NOT_FOUND)
      {
        result = get_default_data(tag_handle, mydata);
          // if failure is returned, change it back to tag not found, since it's
          // ok to look for data and not find any
        if (result == MB_FAILURE) result = MB_TAG_NOT_FOUND;
      }
      if(result != MB_SUCCESS)
        return result;
      mydata += data_size;
    }
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    if (entity_handles.size() == 1)
      return get_data(tag_handle, *entity_handles.begin(), data);
    else
      // don't support this right now - not sure how to pass in multiple bit tags
      return MB_FAILURE;
  }
  else {
      // if we get here, we didn't find the right tag to set
    return MB_TAG_NOT_FOUND;
  }


  return MB_SUCCESS;
}

MBTag TagServer::get_handle(const char *tag_name)
{

  // perhaps speed this up since tag handles are sorted by tag properties
  // then sorted by entity type
  std::map<MBTag, TagInfo>::iterator iterator;
  for(iterator = mTagTable.begin(); iterator != mTagTable.end(); ++iterator)
  {
    if (strcmp(tag_name, iterator->second.get_name().c_str()) == 0)
    {
        return iterator->first;
    }
  }

  return 0;
}

MBErrorCode TagServer::get_tags(std::vector<MBTag> &all_tags)
{
  std::map<MBTag, TagInfo>::iterator iterator;
  for(iterator = mTagTable.begin(); iterator != mTagTable.end(); ++iterator)
  {
    all_tags.push_back(iterator->first);
  }

  return MB_SUCCESS;
}

MBErrorCode TagServer::get_tags(const MBEntityHandle entity, std::vector<MBTag> &all_tags)
{
  MBErrorCode result = MB_SUCCESS;
  
  MBErrorCode tmp_result = mSparseData->get_tags(entity, all_tags);
  if (MB_SUCCESS != tmp_result) result = tmp_result;
  
  tmp_result = mDenseData->get_tags(entity, all_tags);
  if (MB_SUCCESS != tmp_result) result = tmp_result;
  
  tmp_result = mBitServer->get_tags(entity, all_tags);
  if (MB_SUCCESS != tmp_result) result = tmp_result;
  
  return result;
}

MBErrorCode TagServer::get_default_data_ref(const MBTag tag_handle, const void *& data) 
{

  const TagInfo* tag_info = get_tag_info(tag_handle);
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  // get the default value
  const void *def_data = tag_info->default_value();
  // if we have a default value, copy it
  if(def_data != NULL)
  {
    data = def_data; 
    return MB_SUCCESS;
  }

  return MB_FAILURE;

}


MBErrorCode TagServer::get_default_data(const MBTag tag_handle, void *data) 
{

  const TagInfo* tag_info = get_tag_info(tag_handle);
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  // get the default value
  const void *def_data = tag_info->default_value();
  // if we have a default value, copy it
  if(def_data != NULL)
  {
    memcpy(data, def_data, tag_info->get_size());
    return MB_SUCCESS;
  }
  return MB_FAILURE;

}


MBErrorCode TagServer::remove_data( const MBTag tag_handle, const MBEntityHandle entity_handle )
{


  // there is no remove_data equivalent for Dense tags
  if( reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    return MB_SUCCESS;
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    return mSparseData->remove_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle);
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
      // can't take bit tags off entities currently
    return MB_FAILURE;
  }


  // use the delete_tag function instead of this
  //else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_STATIC)
  //{
    //return mStaticSparseData.remove_data(ID_FROM_TAG_HANDLE(tag_handle));
  //}

  return MB_TAG_NOT_FOUND;

}


MBEntityHandle TagServer::find_entity( const MBTag tag_handle, const void* data )
{
  
  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    return mSparseData->find_entity(ID_FROM_TAG_HANDLE(tag_handle), data);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    // todo:  call RMS find entity
    return 0;
  }
  else 
    return 0;

}

//! gets all entity handles that match a type and tag
MBErrorCode TagServer::get_entities(const MBTag tag_handle, const MBEntityType type,
                                     MBRange &entities)
{
  MBErrorCode result = MB_TAG_NOT_FOUND;
  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    result = mSparseData->get_entities(ID_FROM_TAG_HANDLE(tag_handle), type, entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    result = mDenseData->get_entities(ID_FROM_TAG_HANDLE(tag_handle), type, entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    result = mBitServer->get_entities(ID_FROM_TAG_HANDLE(tag_handle), type, entities);
  }
  
  return result;
}

//! gets all entity handles that match a type and tag
MBErrorCode TagServer::get_entities(const MBRange &range,
                                     const MBTag tag_handle, const MBEntityType type,
                                     MBRange &entities)
{
  MBErrorCode result = MB_TAG_NOT_FOUND;
  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    result = mSparseData->get_entities(range, ID_FROM_TAG_HANDLE(tag_handle), type, entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    result = mDenseData->get_entities(range, ID_FROM_TAG_HANDLE(tag_handle), type, entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    result = mBitServer->get_entities(range, ID_FROM_TAG_HANDLE(tag_handle), type, entities);
  }
  
  return result;
}

MBErrorCode TagServer::get_entities_with_tag_value( const MBEntityType type,
                                                     const MBTag tag_handle,
                                                     const void* value,
                                                     MBRange &entities ) 
{

  MBErrorCode result = MB_TAG_NOT_FOUND;

  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    result = mSparseData->get_entities_with_tag_value(ID_FROM_TAG_HANDLE(tag_handle), type, entities, value);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    result = mDenseData->get_entities_with_tag_value(ID_FROM_TAG_HANDLE(tag_handle), type, entities, value);
  }
  else if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    result = mBitServer->get_entities_with_tag_value(ID_FROM_TAG_HANDLE(tag_handle), type, 
                                                     entities, *((const unsigned char*)value));
  }

  return result;
  
}

MBErrorCode TagServer::get_entities_with_tag_value( const MBRange &range,
                                                     const MBEntityType type,
                                                     const MBTag tag_handle,
                                                     const void* value,
                                                     MBRange &entities ) 
{

  MBErrorCode result = MB_TAG_NOT_FOUND;

  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    result = mSparseData->get_entities_with_tag_value(range, ID_FROM_TAG_HANDLE(tag_handle), type, entities, value);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    result = mDenseData->get_entities_with_tag_value(range, ID_FROM_TAG_HANDLE(tag_handle), type, entities, value);
  }

  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    result = mBitServer->get_entities_with_tag_value(range, ID_FROM_TAG_HANDLE(tag_handle), type, 
                                                     entities, *((const unsigned char*)value));
  }

  return result;
  
}

MBErrorCode TagServer::get_entities_with_tag_values( MBEntityType type,
                                                      const MBTag *tags,
                                                      const void** values,
                                                      const int num_tags,
                                                      MBRange &entities,
                                                      const int condition) 
{
  MBErrorCode result;
  MBRange temp1, temp2;

  if (condition != MBInterface::INTERSECT &&
      condition != MBInterface::UNION)
    return MB_FAILURE;

    // if there aren't any values we're looking for, it has to be union
  //This doesn't make sense to me so I removed it -- J.Kraftcheck
  //int temp_condition = (NULL == values ? MBInterface::UNION : condition);
  
  for (unsigned int it = 0; it < (unsigned int) num_tags; it++) {
      // get all entities with this tag/value combo

      // running result is in entities; temp1 and temp2 are working lists
    temp1.clear();
    temp2.clear();

      // get the sets with this tag/value combo in temp1
    if (NULL == values || NULL == values[it]) 
      result = get_entities(tags[it], type, temp1);
    else
      result = get_entities_with_tag_value(type, tags[it], values[it], temp1);

      // if we're doing a running intersection and we're just starting and
      // the list comes in empty, the 1st result is the start
    if (0 == it && condition == MBInterface::INTERSECT && entities.empty()) {
      entities = temp1;
      if (entities.empty()) return MB_SUCCESS;
    }

      // else if we're doing a running intersection, intersect this result (temp1)
      // with the running result (entities) into temp2, then move that to the running
      // result (entities)
    else if (condition == MBInterface::INTERSECT) {
      std::set_intersection(entities.begin(), entities.end(),
                            temp1.begin(), temp1.end(),
                            mb_range_inserter(temp2));
      entities = temp2;
      if (entities.empty()) return MB_SUCCESS;
    }

      // else if we're doing a union, put these results (temp1) into the running result (entities)
      // and re-sort the running result
    else if (condition == MBInterface::UNION) {
      entities.merge(temp1);
    }
  }

    // running result is in entities, where it should be
  return MB_SUCCESS;
}

MBErrorCode TagServer::get_entities_with_tag_values( const MBRange &range,
                                                      const MBEntityType type,
                                                      const MBTag *tags,
                                                      const void** values,
                                                      const int num_tags,
                                                      MBRange &entities,
                                                      const int condition) 
{
  MBErrorCode result;
  MBRange temp1, temp2;

  if (condition != MBInterface::INTERSECT &&
      condition != MBInterface::UNION)
    return MB_FAILURE;

    // if there aren't any values we're looking for, it has to be union
  //This doesn't make sense to me so I removed it -- J.Kraftcheck
  //int temp_condition = (NULL == values ? MBInterface::UNION : condition);
  
  for (unsigned int it = 0; it < (unsigned int) num_tags; it++) {
      // get all entities with this tag/value combo

      // running result is in entities; temp1 and temp2 are working lists
    temp1.clear();
    temp2.clear();

      // get the sets with this tag/value combo in temp1
    if (NULL == values || NULL == values[it]) 
      result = get_entities(tags[it], type, temp1);
    else
      result = get_entities_with_tag_value(type, tags[it], values[it], temp1);

      // if we're doing a running intersection and we're just starting and
      // the list comes in empty, the 1st result is the start
    if (0 == it && condition == MBInterface::INTERSECT && entities.empty()) {
      temp1 = entities;
    }

      // else if we're doing a running intersection, intersect this result (temp1)
      // with the running result (entities) into temp2, then move that to the running
      // result (entities)
    else if (condition == MBInterface::INTERSECT) {
      std::set_intersection(entities.begin(), entities.end(),
                            temp1.begin(), temp1.end(),
                            mb_range_inserter(temp2));
      entities = temp2;
      if (entities.empty()) return MB_SUCCESS;
    }

      // else if we're doing a union, put these results (temp1) into the running result (entities)
      // and re-sort the running result
    else if (condition == MBInterface::UNION) {
      entities.merge(temp1);
    }

    if (!range.empty()) {
        // need to intersect results with what's in range
      temp1.clear();
      std::set_intersection(entities.begin(), entities.end(),
                            range.begin(), range.end(),
                            mb_range_inserter(temp1));
    }
  }


    // running result is in entities, where it should be
  return MB_SUCCESS;
}

MBErrorCode TagServer::get_number_entities( const MBTag tag_handle, 
                                             const MBEntityType type,
                                             int& num_entities)
{
  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    return mSparseData->get_number_entities(ID_FROM_TAG_HANDLE(tag_handle), type, num_entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    return mDenseData->get_number_entities(ID_FROM_TAG_HANDLE(tag_handle), type, num_entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    return mBitServer->get_number_entities(ID_FROM_TAG_HANDLE(tag_handle), type, num_entities);
  }
  return MB_TAG_NOT_FOUND;
}

MBErrorCode TagServer::get_number_entities( const MBRange &range,
                                             const MBTag tag_handle, 
                                             const MBEntityType type,
                                             int& num_entities)
{
  if(reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_SPARSE])
  {
    return mSparseData->get_number_entities(range, ID_FROM_TAG_HANDLE(tag_handle), type, num_entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_DENSE])
  {
    return mDenseData->get_number_entities(range, ID_FROM_TAG_HANDLE(tag_handle), type, num_entities);
  }
  else if (reinterpret_cast<long>(tag_handle) & TAG_BIT_PROPERTIES[MB_TAG_BIT])
  {
    return mBitServer->get_number_entities(range, ID_FROM_TAG_HANDLE(tag_handle), type, num_entities);
  }
  return MB_TAG_NOT_FOUND;
}

#ifdef TEST

#include <iostream>
#include <time.h>
#include <assert.h>
#include <memory.h>

#include "TagServer.hpp"


const int MAX_ID = 10000;
const int SET_TAG_LOOPS = 0xFFF;
const int TEST_LOOPS = 3000;
const int TAG_NAME_SIZE = 3;
const int MAX_TAG_DATA_SIZE = 10;

int main()
{
  TagServer tag_server;

  std::string tag_name;
  tag_name.reserve(TAG_NAME_SIZE+1);
  //char tag_name[TAG_NAME_SIZE+1] = {0};
  MBTag my_tag_handle;
  unsigned char tag_data[MAX_TAG_DATA_SIZE] = {0};
  unsigned char tag_size;
  unsigned int entity_type;
  MBTag tag_handle = (MBTag)3245;


  tag_server.add_tag("densex", 4, TAG_BIT_PROPERTIES[MB_TAG_DENSE]);

  // test for robustness
  tag_server.get_data(tag_handle, 324, tag_data);

  for(int test_loops = 0; test_loops < TEST_LOOPS; test_loops++)
  {
    tag_name.resize(0);
    srand(clock());
    for(int stn=0; stn<TAG_NAME_SIZE; stn++)
    {
      // get a letter between A-D
      tag_name += static_cast<char>((rand()/(RAND_MAX+1.0))*('D' - 'A') + 'A');
    }

    // tag size between 1-MAX
    tag_size = static_cast<char>((rand()/(RAND_MAX+1.0))*MAX_TAG_DATA_SIZE + 1);

    // random entity type
    entity_type = (0x3 & rand()) + 1;

    if((my_tag_handle = tag_server.add_tag(tag_name, tag_size, TAG_BIT_PROPERTIES[MB_TAG_SPARSE])) != 0)
    {
      std::cout << "adding tag - " << tag_name << " " << my_tag_handle << " with entity type " << entity_type << " of size " << static_cast<int>(tag_size) << std::endl;

      MBTag tmp_tag_handle = tag_server.get_handle(tag_name);

      if(tmp_tag_handle != my_tag_handle)
      {
        std::cout << " error ---- " << tag_name << " is " << my_tag_handle << " but get_id returned " << tmp_tag_handle << std::endl;
        assert(0);
      }

      // test for robustness
      tag_server.get_data(my_tag_handle, rand(), tag_data);

      // test the find_entity function.  This test succeeds on failure
      tag_data[0] = rand();
      tag_server.set_data(my_tag_handle, 200, tag_data);
      assert( 200 != tag_server.find_entity(my_tag_handle, tag_data));

      for( int i=0; i<SET_TAG_LOOPS; i++ )
      {
        MBEntityHandle id = (MAX_ID*rand())/(RAND_MAX+1);
        for(int j=0; j<tag_size; j++)
          tag_data[j] = rand();

        unsigned char tag_data_save[MAX_TAG_DATA_SIZE];
        memcpy(tag_data_save, tag_data, sizeof(tag_data));
        if(MB_SUCCESS == tag_server.set_data(my_tag_handle, id, &tag_data))
        {
          for(int j=0; j<tag_size; j++)
            tag_data[j] = rand();

          if(MB_SUCCESS == tag_server.get_data(my_tag_handle, id, &tag_data))
            assert(memcmp(tag_data, tag_data_save, tag_size) == 0);
        }
        tag_server.remove_data(my_tag_handle, id);
      }

    }
    else
    {
      unsigned int entity_type_tmp = (0x3 & rand()) + 1;
      // try to get the handle and remove it
      MBTag tmp = tag_server.get_handle(tag_name); 
      if(tmp)
      {
        std::cout << " removing tag - " << tag_name << " with entity type " << entity_type_tmp << std::endl;
        tag_server.remove_tag(tmp);
      }
    }
  }
 
  std::cout << std::endl << "---SUCCESS---" << std::endl << std::endl; 
  return 0;
}


#endif
