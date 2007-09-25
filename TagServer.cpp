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
#include "SparseTagCollections.hpp"
#include "MBBits.hpp"
#include "DenseTagCollections.hpp"
#include "MBInterface.hpp"

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


void TagInfo::invalidate()
{
  mTagName.clear();
  isValid = false;
  delete [] mDefaultValue;
  mDefaultValue = 0;
  delete [] mMeshValue;
  mMeshValue = 0;
}
  

/*
  TagServer functions ----------------------------------
*/

TagServer::TagServer()
{
  mSparseData = new SparseTagSuperCollection;
  mDenseData = new DenseTagSuperCollection;
  mBitServer = new MBBitServer;
}

TagServer::~TagServer()
{
  delete mSparseData;
  delete mDenseData;
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
    // Check if name is already in use
    // if so, pass back the existing tag handle.
    // NOTE: If name is NULL or empty (tag is unnamed),
    // get_handle will return zero, so no explicit 
    // check is required here.
  tag_handle = get_handle( tag_name );
  if (tag_handle)
    return MB_ALREADY_ALLOCATED;

    // Input size must be a multiple of the size of the data type.
  int typesize = TagInfo::size_from_data_type( data_type );
  if (data_size % typesize)
    return MB_FAILURE;
  
    // data type must be BIT of tag storage type is BIT
  if (storage == MB_TAG_BIT && data_type != MB_TYPE_BIT)
    return MB_FAILURE;
  
    // find an unused tag id
  std::vector<TagInfo>& list = mTagTable[storage];
  std::vector<TagInfo>::iterator i;
  for (i = list.begin(); i != list.end(); ++i)
    if (!i->is_valid())
      break;

    // add TagInfo entry for new tag
  if (i == list.end()) 
    i = list.insert( i, TagInfo( tag_name, data_size, data_type, default_value ) );
  else
    *i = TagInfo( tag_name, data_size, data_type, default_value );


  MBTagId tag_id = i - list.begin() + 1;
  tag_handle = TAG_HANDLE_FROM_ID( tag_id, storage );

  MBErrorCode result = MB_FAILURE;
  switch (storage) {
    case MB_TAG_BIT:
      result = mBitServer->reserve_tag_id( data_size, tag_id );
      break;
    case MB_TAG_SPARSE:
      result = mSparseData->reserve_tag_id(data_size, tag_id);
      break;
    case MB_TAG_DENSE:
      result = mDenseData->reserve_tag_id(data_size, default_value, tag_id);
      break;
    case MB_TAG_MESH:
      result = MB_SUCCESS;
      break;
  }
  
  if (MB_SUCCESS != result)
    i->invalidate();
  return result;
}


MBErrorCode TagServer::remove_tag(const MBTag tag_handle)
{
  const MBTagId tag_id = ID_FROM_TAG_HANDLE( tag_handle );
  const MBTagType tag_type = PROP_FROM_TAG_HANDLE( tag_handle );
  const unsigned tag_idx = tag_id - 1;
  if (tag_idx >= mTagTable[tag_type].size() ||
      !mTagTable[tag_type][tag_idx].is_valid())
    return MB_TAG_NOT_FOUND;
  
  MBErrorCode status = MB_FAILURE;
  switch (tag_type) {
    case MB_TAG_BIT:
      status = mBitServer->release_tag_id(tag_id);
      break;
    case MB_TAG_SPARSE:
      status = mSparseData->release_tag_id(tag_id);
      break;
    case MB_TAG_DENSE:
      status = mDenseData->release_tag_id(tag_id);
      break;
    case MB_TAG_MESH:
      status = MB_SUCCESS;
      break;
  }
  
  mTagTable[tag_type][tag_idx].invalidate();
  return status;
}


//! resets all data tagged to this entity handle back to default data or NULL.
//! this is used to clean out stale data that might be referenced again
MBErrorCode TagServer::reset_data(MBEntityHandle entity_handle)
{
  std::vector<TagInfo>::iterator i;

  MBTagId tag_id;

  for (tag_id = 1; tag_id <= mTagTable[MB_TAG_BIT].size(); ++tag_id) 
    if (mTagTable[MB_TAG_BIT][tag_id-1].is_valid())
        // default data for bits is zero
      mBitServer->weak_set_bits( tag_id, entity_handle, 0 );

  for (tag_id = 1; tag_id <= mTagTable[MB_TAG_SPARSE].size(); ++tag_id) 
    if (mTagTable[MB_TAG_SPARSE][tag_id-1].is_valid())
      mSparseData->remove_data( tag_id, entity_handle );

  for (tag_id = 1; tag_id <= mTagTable[MB_TAG_DENSE].size(); ++tag_id) 
    if (mTagTable[MB_TAG_DENSE][tag_id-1].is_valid())
      mDenseData->remove_data( tag_id, entity_handle );

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

MBErrorCode TagServer::set_mesh_data( const MBTag tag_handle,
                                      const void* data )
{
  TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  
  info->set_mesh_value( data );
  return MB_SUCCESS;
}

MBErrorCode TagServer::set_data(const MBTag tag_handle, 
                                 const MBEntityHandle entity_handle, 
                                 const void* data)
{

  // this assumes that the entity_handle is valid, 
  // the check for valid handles are done one level up
  
  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;
  
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_BIT:
      return set_bits(tag_handle, entity_handle, *((unsigned char *)data));
    case MB_TAG_SPARSE:
      return mSparseData->set_data(id, entity_handle, data);
    case MB_TAG_DENSE:
      return mDenseData->set_data(id, entity_handle, data);
    case MB_TAG_MESH:
      return MB_FAILURE;
  }
  return MB_FAILURE;
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
  
  if( PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_SPARSE)
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
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
  
  if( PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_SPARSE)
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
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


MBErrorCode TagServer::get_mesh_data( const MBTag tag_handle,
                                      void* data ) const
{
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info || !info->get_mesh_value())
    return MB_TAG_NOT_FOUND;
  
  memcpy( data, info->get_mesh_value(), info->get_size() );
  return MB_SUCCESS;
}

MBErrorCode TagServer::get_data(const MBTag tag_handle,
                       const MBEntityHandle entity_handle,
                       void* data)
{

  MBErrorCode result = MB_TAG_NOT_FOUND;

  if(TYPE_FROM_HANDLE(entity_handle) >= MBMAXTYPE)
    return MB_TYPE_OUT_OF_RANGE;

  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_DENSE:
      result = mDenseData->get_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
      break;
    case MB_TAG_SPARSE:
      result = mSparseData->get_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle, data);
      break;
    case MB_TAG_BIT:
      result = get_bits(tag_handle, entity_handle, *((unsigned char *)data));
      break;
    case MB_TAG_MESH:
      result = MB_FAILURE;
      break;
  }

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

  if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
  {
    const int data_size = mDenseData->tag_size(tag_id);
    for(const MBEntityHandle* iter = entity_handles; iter != end; ++iter)
    {
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_SPARSE)
  {
    const int data_size = mSparseData->tag_size(tag_id);
    for(const MBEntityHandle* iter = entity_handles; iter != end; ++iter)
    {
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
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

  if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_SPARSE)
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
  else if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
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

MBTag TagServer::get_handle(const char *tag_name) const
{
  if (tag_name && *tag_name)
    for (int i = 0; i < MB_TAG_LAST+1; ++i) 
      for (MBTagId j = 0; j < mTagTable[i].size(); ++j) 
        if (mTagTable[i][j].is_valid() && mTagTable[i][j].get_name() == tag_name)
          return TAG_HANDLE_FROM_ID( j + 1, (MBTagType)i );
  
  return 0;
}

MBErrorCode TagServer::get_tags(std::vector<MBTag> &all_tags)
{
  for (int i = 0; i < MB_TAG_LAST+1; ++i) 
    for (MBTagId j = 0; j < mTagTable[i].size(); ++j) 
      if (mTagTable[i][j].is_valid())
        all_tags.push_back( TAG_HANDLE_FROM_ID( j + 1, (MBTagType)i ) );

  return MB_SUCCESS;
}

MBErrorCode TagServer::get_tags(MBDataType type, std::vector<MBTag> &all_tags)
{
  for (int i = 0; i < MB_TAG_LAST+1; ++i) 
    for (MBTagId j = 0; j < mTagTable[i].size(); ++j) 
      if (mTagTable[i][j].is_valid() && mTagTable[i][j].get_data_type() == type)
        all_tags.push_back( TAG_HANDLE_FROM_ID( j + 1, (MBTagType)i ) );

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

MBErrorCode TagServer::get_mesh_tags( std::vector<MBTag>& all_tags ) const
{
  for (int i = 0; i < MB_TAG_LAST+1; ++i) 
    for (MBTagId j = 0; j < mTagTable[i].size(); ++j) 
      if (mTagTable[i][j].is_valid() && mTagTable[i][j].get_mesh_value())
        all_tags.push_back( TAG_HANDLE_FROM_ID( j + 1, (MBTagType)i ) );
  
  return MB_SUCCESS;
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


MBErrorCode TagServer::remove_mesh_data( const MBTag tag_handle )
{
  TagInfo* info = get_tag_info( tag_handle );
  if (!info || !info->get_mesh_value())
    return MB_TAG_NOT_FOUND;
  info->remove_mesh_value();
  return MB_SUCCESS;
}


MBErrorCode TagServer::remove_data( const MBTag tag_handle, const MBEntityHandle entity_handle )
{
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_DENSE:
      return MB_SUCCESS;
    case MB_TAG_SPARSE:
      return mSparseData->remove_data(ID_FROM_TAG_HANDLE(tag_handle), entity_handle);
    case MB_TAG_BIT:
      return MB_FAILURE;
    case MB_TAG_MESH:
      return MB_FAILURE;
  }
  
  return MB_TAG_NOT_FOUND;
}


MBEntityHandle TagServer::find_entity( const MBTag tag_handle, const void* data )
{
  
  if(PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_SPARSE)
  {
    return mSparseData->find_entity(ID_FROM_TAG_HANDLE(tag_handle), data);
  }
  else if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_DENSE)
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
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities(id, type, entities);
      break;
    case MB_TAG_DENSE:
      result = mDenseData->get_entities(id, type, entities);
      break;
    case MB_TAG_BIT:
      result = mBitServer->get_entities(id, type, entities);
      break;
    case MB_TAG_MESH:
      result = MB_TYPE_OUT_OF_RANGE;
      break;
  }
  
  return result;
}

//! gets all entity handles that match a type and tag
MBErrorCode TagServer::get_entities(const MBRange &range,
                                     const MBTag tag_handle, const MBEntityType type,
                                     MBRange &entities)
{
  MBErrorCode result = MB_TAG_NOT_FOUND;
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities(range, id, type, entities);
      break;
    case MB_TAG_DENSE:
      result = mDenseData->get_entities(range, id, type, entities);
      break;
    case MB_TAG_BIT:
      result = mBitServer->get_entities(range, id, type, entities);
      break;
    case MB_TAG_MESH:
      result = MB_TYPE_OUT_OF_RANGE;
      break;
  }
  
  return result;
}

MBErrorCode TagServer::get_entities_with_tag_value( const MBEntityType type,
                                                     const MBTag tag_handle,
                                                     const void* value,
                                                     MBRange &entities ) 
{

  MBErrorCode result = MB_TAG_NOT_FOUND;
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities_with_tag_value(id, type, entities, value);
      break;
    case MB_TAG_DENSE:
      result = mDenseData->get_entities_with_tag_value(id, type, entities, value);
      break;
    case MB_TAG_BIT:
      result = mBitServer->get_entities_with_tag_value(id, type, 
                                entities, *((const unsigned char*)value));
      break;
    case MB_TAG_MESH:
      result = MB_TYPE_OUT_OF_RANGE;
      break;
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
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities_with_tag_value(range, id, type, entities, value);
      break;
    case MB_TAG_DENSE:
      result = mDenseData->get_entities_with_tag_value(range, id, type, entities, value);
      break;
    case MB_TAG_BIT:
      result = mBitServer->get_entities_with_tag_value(range, id, type, 
                                   entities, *((const unsigned char*)value));
      break;
    case MB_TAG_MESH:
      result = MB_TYPE_OUT_OF_RANGE;
      break;
  }

  return result;
  
}

MBErrorCode TagServer::get_entities_with_tag_values( const MBRange &input_range,
                                                      const MBEntityType type,
                                                      const MBTag *tags,
                                                      const void* const* values,
                                                      const int num_tags,
                                                      MBRange &entities,
                                                      const int condition) 
{
    // range should never come in empty
  assert(!input_range.empty());
  if (input_range.empty()) return MB_FAILURE;
  MBRange range = input_range;
  
  MBErrorCode result;
  MBRange temp1;

  if (condition != MBInterface::INTERSECT &&
      condition != MBInterface::UNION)
    return MB_FAILURE;

  for (unsigned int it = 0; it < (unsigned int) num_tags; it++) {
      // get all entities with this tag/value combo

      // running result is in entities; temp1 and temp2 are working lists
    temp1.clear();

      // get the sets with this tag/value combo in temp1
    if (NULL == values || NULL == values[it]) 
      result = get_entities(range, tags[it], type, temp1);
    else
      result = get_entities_with_tag_value(range, type, tags[it], values[it], temp1);

      // if we're doing a running intersection and we're just starting and
      // the list comes in empty, the 1st result is the start
    if (0 == it && condition == MBInterface::INTERSECT && entities.empty()) {
      entities = temp1.intersect(range);
    }

      // else if we're doing a running intersection, intersect this result (temp1)
      // with the running result (entities) into temp2, then move that to the running
      // result (entities)
    else if (condition == MBInterface::INTERSECT) {
      entities = entities.intersect(temp1);
      if (entities.empty()) return MB_SUCCESS;

        // also restrict the range at which we look; entities has already been 
        // intersected with range (through input to get_entities above) so just assign
      range = entities;
    }

      // else if we're doing a union, put these results (temp1) into the running 
      // result (entities)
    else if (condition == MBInterface::UNION) {
      entities.merge(temp1);
    }
  }

    // running result is in entities, where it should be
  return MB_SUCCESS;
}

MBErrorCode TagServer::get_number_entities( const MBTag tag_handle,
                                            unsigned long& num_entities )
{
  num_entities = 0;
  for (MBEntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
    int type_num_entities;
    MBErrorCode rval = get_number_entities( tag_handle, t, type_num_entities );
    if (MB_SUCCESS != rval) 
      return rval;
    num_entities += type_num_entities;
  }
  return MB_SUCCESS;
}

MBErrorCode TagServer::get_number_entities( const MBTag tag_handle, 
                                             const MBEntityType type,
                                             int& num_entities)
{
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      return mSparseData->get_number_entities(id, type, num_entities);
    case MB_TAG_DENSE:
      return mDenseData->get_number_entities(id, type, num_entities);
    case MB_TAG_BIT:
      return mBitServer->get_number_entities(id, type, num_entities);
    case MB_TAG_MESH:
      return MB_TYPE_OUT_OF_RANGE;
  }
  return MB_TAG_NOT_FOUND;
}

MBErrorCode TagServer::get_number_entities( const MBRange &range,
                                             const MBTag tag_handle, 
                                             const MBEntityType type,
                                             int& num_entities)
{
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      return mSparseData->get_number_entities(range, id, type, num_entities);
    case MB_TAG_DENSE:
      return mDenseData->get_number_entities(range, id, type, num_entities);
    case MB_TAG_BIT:
      return mBitServer->get_number_entities(range, id, type, num_entities);
    case MB_TAG_MESH:
      return MB_TYPE_OUT_OF_RANGE;
  }
  return MB_TAG_NOT_FOUND;
}

unsigned long TagServer::get_memory_use( MBTag tag_handle ) const
{
  if (!get_tag_info(tag_handle))
    return 0;

  unsigned long result = 0, tmp;
  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      mSparseData->get_memory_use( id, result, tmp );
      break;
    case MB_TAG_DENSE:
      mDenseData->get_memory_use( id, result, tmp );
      break;
    case MB_TAG_BIT:
      mBitServer->get_memory_use( id, result, tmp );
      break;
    case MB_TAG_MESH:
      break;
  }

    // add in size of entry in mTagTable
  return result + sizeof(MBTag) + sizeof(TagInfo) + 3*sizeof(void*);
}

MBErrorCode TagServer::get_memory_use( MBTag tag_handle,
                                       unsigned long& total,
                                       unsigned long& per_entity ) const
{
  const TagInfo* tag_info = get_tag_info(tag_handle);
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  MBTagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      mSparseData->get_memory_use( id, total, per_entity );
      break;
    case MB_TAG_DENSE:
      mDenseData->get_memory_use( id, total, per_entity );
      break;
    case MB_TAG_BIT:
      mBitServer->get_memory_use( id, total, per_entity );
      break;
    case MB_TAG_MESH:
      break;
  }
  
    // size of entry in mTagTable map
  total += sizeof(MBTag) + sizeof(TagInfo) + 3*sizeof(void*);
  if (tag_info->default_value())
    total += tag_info->get_size();
  if (tag_info->get_mesh_value())
    total += tag_info->get_size();
  total += tag_info->get_name().size();
  
  return MB_SUCCESS;
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


  tag_server.add_tag("densex", 4, MB_TAG_DENSE, MB_TYPE_OPAQUE, my_tag_handle);

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

    if(MB_SUCCESS == tag_server.add_tag(tag_name.c_str(), tag_size, MB_TAG_SPARSE, MB_TYPE_OPAQUE, my_tag_handle))
    {
      std::cout << "adding tag - " << tag_name << " " << my_tag_handle << " with entity type " << entity_type << " of size " << static_cast<int>(tag_size) << std::endl;

      MBTag tmp_tag_handle = tag_server.get_handle(tag_name.c_str());

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
      MBTag tmp = tag_server.get_handle(tag_name.c_str()); 
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
