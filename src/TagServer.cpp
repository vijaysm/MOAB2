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
#include "moab/Range.hpp"
#include "SparseTagSuperCollection.hpp"
#include "BitTagServer.hpp"
#include "moab/Interface.hpp"
#include "SequenceManager.hpp"
#include "TagCompare.hpp"

namespace moab {

using namespace std;
  
static inline bool check_lengths( const int* length_array, int expected_value, int count )
{
  if (length_array)
    for (int i = 0; i < count; ++i)
      if (length_array[i] != expected_value)
        return false;
  return true;
}

/*
  TagServer functions ----------------------------------
*/

TagServer::TagServer( SequenceManager* seqman )
  : sequenceManager(seqman)
{
  mSparseData = new SparseTagSuperCollection;
  mBitServer = new BitTagServer;
}

TagServer::~TagServer()
{
  sequenceManager->reset_tag_data();
  delete mSparseData;
  delete mBitServer;
}


ErrorCode TagServer::reset_all_data()
{
  mSparseData->reset_data();
  mBitServer->reset_data();
  sequenceManager->reset_tag_data();

  return MB_SUCCESS;
}


ErrorCode TagServer::add_tag( const char *tag_name, 
                                const int data_size,
                                const TagType storage,
                                const DataType data_type,
                                Tag &tag_handle,
                                const void *default_value,
                                int default_value_size)
{
  if (!default_value) {
    default_value_size = 0;
  }
  else if (storage == MB_TAG_BIT) {
      // bit tags cannot be zero-length or variable-length
    if (data_size < 1)
      return MB_INVALID_SIZE;
    default_value_size = (data_size+7)/8; // convert from bits to bytes
  }
  else if (data_size != MB_VARIABLE_LENGTH) 
    default_value_size = data_size;
  else if (default_value_size % TagInfo::size_from_data_type(data_type))
    return MB_INVALID_SIZE;

    // Check if name is already in use
    // if so, pass back the existing tag handle.
    // NOTE: If name is NULL or empty (tag is unnamed),
    // get_handle will return zero, so no explicit 
    // check is required here.
  tag_handle = get_handle( tag_name );
  if (tag_handle)
    return MB_ALREADY_ALLOCATED;

    // Input size must be a multiple of the size of the data type.
  if (data_size != MB_VARIABLE_LENGTH) {
    int typesize = TagInfo::size_from_data_type( data_type );
    if (data_size % typesize)
      return MB_INVALID_SIZE;
  }
  
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
    i = list.insert( i, TagInfo( tag_name, data_size, data_type, default_value, default_value_size ) );
  else
    *i = TagInfo( tag_name, data_size, data_type, default_value, default_value_size );


  TagId tag_id = i - list.begin() + 1;
  tag_handle = TAG_HANDLE_FROM_ID( tag_id, storage );

  ErrorCode result = MB_FAILURE;
  switch (storage) {
    case MB_TAG_BIT:
      result = mBitServer->reserve_tag_id( data_size, tag_id );
      break;
    case MB_TAG_SPARSE:
      result = mSparseData->reserve_tag_id(data_size, tag_id);
      break;
    case MB_TAG_DENSE:
      result = sequenceManager->reserve_tag_id( data_size, tag_id );
      break;
    case MB_TAG_MESH:
      result = MB_SUCCESS;
      break;
  }
  
  if (MB_SUCCESS != result)
    i->invalidate();
  return result;
}


ErrorCode TagServer::remove_tag(const Tag tag_handle)
{
  const TagId tag_id = ID_FROM_TAG_HANDLE( tag_handle );
  const TagType tag_type = PROP_FROM_TAG_HANDLE( tag_handle );
  const unsigned tag_idx = tag_id - 1;
  if (tag_idx >= mTagTable[tag_type].size() ||
      !mTagTable[tag_type][tag_idx].is_valid())
    return MB_TAG_NOT_FOUND;
  
  ErrorCode status = MB_FAILURE;
  switch (tag_type) {
    case MB_TAG_BIT:
      status = mBitServer->release_tag_id(tag_id);
      break;
    case MB_TAG_SPARSE:
      status = mSparseData->release_tag_id(tag_id);
      break;
    case MB_TAG_DENSE:
      status = sequenceManager->release_tag(tag_id);
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
ErrorCode TagServer::reset_data(EntityHandle entity_handle)
{
  std::vector<TagInfo>::iterator i;

  TagId tag_id;

  for (tag_id = 1; tag_id <= mTagTable[MB_TAG_BIT].size(); ++tag_id) 
    if (mTagTable[MB_TAG_BIT][tag_id-1].is_valid())
      mBitServer->clear_bits( tag_id, &entity_handle, 1,
                              reinterpret_cast<const unsigned char*>
                              (mTagTable[MB_TAG_DENSE][tag_id-1].default_value()) );

  for (tag_id = 1; tag_id <= mTagTable[MB_TAG_SPARSE].size(); ++tag_id) 
    if (mTagTable[MB_TAG_SPARSE][tag_id-1].is_valid())
      mSparseData->remove_data( tag_id, entity_handle );

  for (tag_id = 1; tag_id <= mTagTable[MB_TAG_DENSE].size(); ++tag_id) 
    if (mTagTable[MB_TAG_DENSE][tag_id-1].is_valid())
      sequenceManager->remove_tag_data( tag_id, entity_handle,
                                        mTagTable[MB_TAG_DENSE][tag_id-1].default_value() );

  return MB_SUCCESS;
}


ErrorCode TagServer::set_mesh_data( const Tag tag_handle,
                                      const void* data,
                                      int size )
{
  TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  
  if (info->get_size() == MB_VARIABLE_LENGTH) {
    if (!size || !info->check_valid_sizes( &size, 1 ))
      return MB_INVALID_SIZE;
  }
  else if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
    size = 1; // one byte max for bit tag values
  else
    size = info->get_size();
  
  info->set_mesh_value( data, size );
  return MB_SUCCESS;
}


ErrorCode TagServer::set_data( const Tag tag_handle, 
                                 const EntityHandle* entity_handles, 
                                 const int num_entities,
                                 const void* data )
{
  ErrorCode rval;
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType tag_type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* tag_info;
  switch (tag_type) {
    case MB_TAG_DENSE:
      if (!(tag_info = get_tag_info( tag_id, tag_type )))
        rval = MB_TAG_NOT_FOUND;
      else
        rval = sequenceManager->set_tag_data( tag_id, entity_handles, num_entities, 
                                            data, tag_info->default_value() );
      break;
        
    case MB_TAG_SPARSE:
      rval = sequenceManager->check_valid_entities( entity_handles, num_entities );
      if (MB_SUCCESS == rval)
        rval = mSparseData->set_data( tag_id, entity_handles, num_entities, data );
      break;
      
    case MB_TAG_BIT:
      if (!(tag_info = get_tag_info( tag_id, tag_type )))
        rval = MB_TAG_NOT_FOUND;
      else
        rval = sequenceManager->check_valid_entities( entity_handles, num_entities );
      if (MB_SUCCESS == rval)
        rval= mBitServer->set_bits( tag_id, entity_handles, num_entities, 
                                    reinterpret_cast<const unsigned char*>(data),
                                    reinterpret_cast<const unsigned char*>(tag_info->default_value()));
      break;
    
    default:
      rval = MB_TAG_NOT_FOUND;
      break;
  }
  return rval;
}

ErrorCode TagServer::set_data( const Tag tag_handle, 
                                 const Range& entity_handles, 
                                 const void* data )
{
  ErrorCode rval;
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType tag_type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* tag_info;
  switch (tag_type) {
    case MB_TAG_DENSE:
      if (!(tag_info = get_tag_info( tag_id, tag_type )))
        rval = MB_TAG_NOT_FOUND;
      else
        rval = sequenceManager->set_tag_data( tag_id, entity_handles, 
                                            data, tag_info->default_value() );
      break;
    
    case MB_TAG_SPARSE:
      rval = sequenceManager->check_valid_entities( entity_handles );
      if (MB_SUCCESS == rval)
        rval = mSparseData->set_data( tag_id, entity_handles, data );
      break;
    
    case MB_TAG_BIT:
      if (!(tag_info = get_tag_info( tag_id, tag_type )))
        rval = MB_TAG_NOT_FOUND;
      else
        rval = sequenceManager->check_valid_entities( entity_handles );
      if (MB_SUCCESS == rval)
        rval= mBitServer->set_bits( tag_id, entity_handles, 
                                    reinterpret_cast<const unsigned char*>(data),
                                    reinterpret_cast<const unsigned char*>(tag_info->default_value()));
      break;
    
    default:
      rval = MB_TAG_NOT_FOUND;
      break;
  }
  
  return rval;
}

ErrorCode TagServer::set_data( const Tag tag_handle, 
                                 const EntityHandle* entity_handles, 
                                 const int num_entities,
                                 void const* const* data,
                                 const int* lengths )
{
  ErrorCode rval;
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType tag_type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* tag_info = get_tag_info( tag_id, tag_type );
  if (!tag_info)
    return MB_TAG_NOT_FOUND;
    
  if (tag_info->get_size() == MB_VARIABLE_LENGTH &&
      !tag_info->check_valid_sizes( lengths, num_entities ))
    return MB_INVALID_SIZE;
    
  switch (tag_type) {
    case MB_TAG_DENSE:
      rval = sequenceManager->set_tag_data( tag_id, entity_handles, num_entities, 
                   data, lengths, tag_info->default_value() );
      break;
      
    case MB_TAG_SPARSE:
      rval = sequenceManager->check_valid_entities( entity_handles, num_entities );
      if (MB_SUCCESS == rval)
        rval = mSparseData->set_data( tag_id, entity_handles, num_entities, data, lengths );
      break;
    
    case MB_TAG_BIT:
      rval = MB_FAILURE;
      break;
    
    default:
      rval = MB_TAG_NOT_FOUND;
      break;
  }
  
  return rval;
}

ErrorCode TagServer::set_data( const Tag tag_handle, 
                                 const Range& entity_handles, 
                                 void const* const* data,
                                 const int* lengths )
{
  ErrorCode rval;
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType tag_type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* tag_info = get_tag_info( tag_id, tag_type );
  if (!tag_info)
    return MB_TAG_NOT_FOUND;
    
  if (tag_info->get_size() == MB_VARIABLE_LENGTH &&
      !tag_info->check_valid_sizes( lengths, entity_handles.size() ))
    return MB_INVALID_SIZE;
    
  switch (tag_type) {
    case MB_TAG_DENSE:
      rval = sequenceManager->set_tag_data( tag_id, entity_handles, 
                      data, lengths, tag_info->default_value() );
      break;
      
    case MB_TAG_SPARSE:
      rval = sequenceManager->check_valid_entities( entity_handles );
      if (MB_SUCCESS == rval)
        rval = mSparseData->set_data( tag_id, entity_handles, data, lengths );
      break;
    
    case MB_TAG_BIT:
      rval = MB_FAILURE;
      break;
    
    default:
      rval = MB_TAG_NOT_FOUND;
      break;
  }
  
  return rval;
}



ErrorCode TagServer::get_mesh_data( const Tag tag_handle,
                                      void* data ) const
{
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  
  if (info->get_size() == MB_VARIABLE_LENGTH)
    return MB_VARIABLE_DATA_LENGTH;
  
  if (info->get_mesh_value())
    memcpy( data, info->get_mesh_value(), info->get_mesh_value_size() );
  else if (info->default_value())
    memcpy( data, info->default_value(), info->default_value_size() );
  else
    return MB_TAG_NOT_FOUND;
  
  return MB_SUCCESS;
}

ErrorCode TagServer::get_mesh_data( Tag tag_handle,
                                      const void*& data,
                                      int& size ) const
{
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  
  if (info->get_mesh_value()) {
    size = info->get_mesh_value_size();
    data = info->get_mesh_value();
  }
  else if (info->default_value()) {
    size = info->default_value_size();
    data = info->default_value();
  }
  else {
    data = 0;
    size = 0;
    return MB_TAG_NOT_FOUND;
  }
    
    // Mesh value and default value sizes are the number of bytes.
    // For bit tags, the number of bytes is always 1.  We want the
    // number of bits instead.
  if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
    size = info->get_size();
  
  return MB_SUCCESS;
}

ErrorCode TagServer::get_data( const Tag tag_handle,
                                 const EntityHandle* entity_handles,
                                 const int num_entities,
                                 void* data )
{
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  const void* const default_val = info->default_value();
  switch (type) {
    case MB_TAG_DENSE:
      return sequenceManager->get_tag_data( tag_id, entity_handles, num_entities, data, default_val );
    case MB_TAG_SPARSE:
      return mSparseData->get_data( tag_id, entity_handles, num_entities, data, default_val );
    case MB_TAG_BIT:
      return mBitServer->get_bits( tag_id, entity_handles, num_entities, 
                                   reinterpret_cast<unsigned char*>(data),
                                   reinterpret_cast<const unsigned char*>(default_val));
    default:
      return MB_TAG_NOT_FOUND;
  }
}

ErrorCode TagServer::get_data( const Tag tag_handle,
                                 const Range& entity_handles,
                                 void* data )
{
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  const void* const default_val = info->default_value();
  switch (type) {
    case MB_TAG_DENSE:
      return sequenceManager->get_tag_data( tag_id, entity_handles, data, default_val );
    case MB_TAG_SPARSE:
      return mSparseData->get_data( tag_id, entity_handles, data, default_val );
    case MB_TAG_BIT:
      return mBitServer->get_bits( tag_id, entity_handles,
                                   reinterpret_cast<unsigned char*>(data),
                                   reinterpret_cast<const unsigned char*>(default_val));
    default:
      return MB_TAG_NOT_FOUND;
  }
}

ErrorCode TagServer::get_data( const Tag tag_handle,
                                 const EntityHandle* entity_handles,
                                 const int num_entities,
                                 const void** data,
                                 int* lengths )
{
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  const void* const default_val = info->default_value();
  const int def_val_len = info->default_value_size();
  switch (type) {
    case MB_TAG_DENSE:
      return sequenceManager->get_tag_data( tag_id, entity_handles, num_entities, data, lengths, default_val, def_val_len );
    case MB_TAG_SPARSE:
      return mSparseData->get_data( tag_id, entity_handles, num_entities, data, lengths, default_val, def_val_len );
    case MB_TAG_BIT:
      return MB_FAILURE;
    default:
      return MB_TAG_NOT_FOUND;
  }
}

ErrorCode TagServer::get_data( const Tag tag_handle,
                                 const Range& entity_handles,
                                 const void** data,
                                 int* lengths )
{
  const TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  const TagType type = PROP_FROM_TAG_HANDLE(tag_handle);
  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  const void* const default_val = info->default_value();
  const int def_val_len = info->default_value_size();
  switch (type) {
    case MB_TAG_DENSE:
      return sequenceManager->get_tag_data( tag_id, entity_handles, data, lengths, default_val, def_val_len );
    case MB_TAG_SPARSE:
      return mSparseData->get_data( tag_id, entity_handles, data, lengths, default_val, def_val_len );
    case MB_TAG_BIT:
      return MB_FAILURE;
    default:
      return MB_TAG_NOT_FOUND;
  }
}


Tag TagServer::get_handle(const char *tag_name) const
{
  if (tag_name && *tag_name)
    for (int i = 0; i < MB_TAG_LAST+1; ++i) 
      for (TagId j = 0; j < mTagTable[i].size(); ++j) 
        if (mTagTable[i][j].is_valid() && mTagTable[i][j].get_name() == tag_name)
          return TAG_HANDLE_FROM_ID( j + 1, (TagType)i );
  
  return 0;
}

ErrorCode TagServer::get_tags(std::vector<Tag> &all_tags)
{
  for (int i = 0; i < MB_TAG_LAST+1; ++i) 
    for (TagId j = 0; j < mTagTable[i].size(); ++j) 
      if (mTagTable[i][j].is_valid())
        all_tags.push_back( TAG_HANDLE_FROM_ID( j + 1, (TagType)i ) );

  return MB_SUCCESS;
}

ErrorCode TagServer::get_tags(const EntityHandle entity, std::vector<Tag> &all_tags)
{
  ErrorCode result = MB_SUCCESS;
  
  ErrorCode tmp_result = mSparseData->get_tags(entity, all_tags);
  if (MB_SUCCESS != tmp_result) result = tmp_result;
  
  tmp_result = sequenceManager->get_entity_tags(entity, all_tags);
  if (MB_SUCCESS != tmp_result) result = tmp_result;
  
  tmp_result = mBitServer->get_tags(entity, all_tags);
  if (MB_SUCCESS != tmp_result) result = tmp_result;
  
  return result;
}

ErrorCode TagServer::get_mesh_tags( std::vector<Tag>& all_tags ) const
{
  for (int i = 0; i < MB_TAG_LAST+1; ++i) 
    for (TagId j = 0; j < mTagTable[i].size(); ++j) 
      if (mTagTable[i][j].is_valid() && mTagTable[i][j].get_mesh_value())
        all_tags.push_back( TAG_HANDLE_FROM_ID( j + 1, (TagType)i ) );
  
  return MB_SUCCESS;
}

ErrorCode TagServer::get_default_data_ref( const Tag tag_handle, 
                                             const void *& data,
                                             int& size) 
{

  const TagInfo* tag_info = get_tag_info(tag_handle);
  if(!tag_info)
    return MB_TAG_NOT_FOUND;
  
  size = tag_info->default_value_size();
  if (size) {
    data = tag_info->default_value();
    if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
      size = tag_info->get_size(); // get number of bits rather than bytes
    return MB_SUCCESS;
  }
  else {
    data = 0;
    return MB_ENTITY_NOT_FOUND;
  }
}


ErrorCode TagServer::get_default_data( const Tag tag_handle, 
                                         void *data,
                                         int& size) 
{

  const TagInfo* tag_info = get_tag_info(tag_handle);
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  size = tag_info->default_value_size();
  if (!size)
    return MB_ENTITY_NOT_FOUND;
  if (data)
    memcpy( data, tag_info->default_value(), size );
  if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
    size = tag_info->get_size(); // get bits rather than bytes
  return MB_SUCCESS;
}


ErrorCode TagServer::remove_mesh_data( const Tag tag_handle )
{
  TagInfo* info = get_tag_info( tag_handle );
  if (!info || !info->get_mesh_value_size())
    return MB_TAG_NOT_FOUND;
  info->remove_mesh_value();
  return MB_SUCCESS;
}


ErrorCode TagServer::remove_data( const Tag tag_handle, const EntityHandle entity_handle )
{
  TagId tag_id = ID_FROM_TAG_HANDLE(tag_handle);
  TagType storage = PROP_FROM_TAG_HANDLE(tag_handle);
  if (tag_id > mTagTable[storage].size())
    return MB_TAG_NOT_FOUND;
  const void* defval = mTagTable[storage][tag_id-1].default_value();
  
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_DENSE:
      sequenceManager->remove_tag_data( tag_id, entity_handle, defval );
      return MB_SUCCESS;
    case MB_TAG_SPARSE:
      return mSparseData->remove_data(tag_id, entity_handle);
    case MB_TAG_BIT:
      return mBitServer->clear_bits( tag_id, &entity_handle, 1,
                                static_cast<const unsigned char*>(defval) );
    case MB_TAG_MESH:
      return MB_FAILURE;
  }
  
  return MB_TAG_NOT_FOUND;
}


//! gets all entity handles that match a type and tag
ErrorCode TagServer::get_entities(const Tag tag_handle, const EntityType type,
                                     Range &entities)
{
  ErrorCode result = MB_TAG_NOT_FOUND;
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities(id, type, entities);
      break;
    case MB_TAG_DENSE:
      result = sequenceManager->get_tagged_entities(id, type, entities);
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

//! gets all entity handles that match a tag
ErrorCode TagServer::get_entities(const Tag tag_handle, Range &entities)
{
  ErrorCode result = MB_TAG_NOT_FOUND;
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities(id, entities);
      break;
    case MB_TAG_DENSE:
      for (EntityType t = MBENTITYSET; t >= MBVERTEX; --t)
        result = sequenceManager->get_tagged_entities(id, t, entities);
      break;
    case MB_TAG_BIT:
      result = mBitServer->get_entities(id, entities);
      break;
    case MB_TAG_MESH:
      result = MB_TYPE_OUT_OF_RANGE;
      break;
  }
  
  return result;
}

//! gets all entity handles that match a type and tag
ErrorCode TagServer::get_entities(const Range &range,
                                     const Tag tag_handle, const EntityType type,
                                     Range &entities)
{
  ErrorCode result = MB_TAG_NOT_FOUND;
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities(range, id, type, entities);
      break;
    case MB_TAG_DENSE: {
      Range temp;
      result = sequenceManager->get_tagged_entities(id, type, temp);
      entities.merge( intersect( range, temp) );
      }
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

ErrorCode TagServer::get_entities_with_tag_value( const EntityType type,
                                                    const Tag tag_handle,
                                                    const void* value,
                                                    Range &entities,
                                                    int value_size ) 
{

  ErrorCode result = MB_TAG_NOT_FOUND;
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);

  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
  if (!value_size && info->get_size() != MB_VARIABLE_LENGTH)
    value_size = info->get_size();
  
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities_with_tag_value(id, *info, type, entities, value, value_size);
      break;
    case MB_TAG_DENSE:
      result = sequenceManager->get_entities_with_tag_value(id, *info, type, entities, value, value_size);
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

ErrorCode TagServer::get_entities_with_tag_value( const Range &range,
                                                    const EntityType type,
                                                    const Tag tag_handle,
                                                    const void* value,
                                                    Range &entities,
                                                    int value_size ) 
{

  ErrorCode result = MB_TAG_NOT_FOUND;
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);

  const TagInfo* info = get_tag_info( tag_handle );
  if (!info)
    return MB_TAG_NOT_FOUND;
    
  if (!value_size) {
    if (info->get_size() == MB_VARIABLE_LENGTH)
      return MB_VARIABLE_DATA_LENGTH;
    value_size = info->get_size();
  }
  
    // If tag value is default value, then we want every entity
    // in 'range' of the correct type, except those with a different tag value.
  bool equals_default = false;
  if (info->default_value()) {
    if (PROP_FROM_TAG_HANDLE(tag_handle) == MB_TAG_BIT)
      equals_default = (*(char*)value == *(char*)info->default_value());
    else
      equals_default = !memcmp( value, info->default_value(), info->default_value_size() );
  }
  Range tmp_ents;
  if (equals_default)
    tmp_ents.swap( entities );
  
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      result = mSparseData->get_entities_with_tag_value(range, id, *info, type, entities, value, value_size);
      break;
    case MB_TAG_DENSE:
      result = sequenceManager->get_entities_with_tag_value(range, id, *info, type, entities, value, value_size);
      break;
    case MB_TAG_BIT:
      result = mBitServer->get_entities_with_tag_value(range, id, type, 
                                   entities, *((const unsigned char*)value));
      break;
    case MB_TAG_MESH:
      result = MB_TYPE_OUT_OF_RANGE;
      break;
  }

  if (MB_SUCCESS != result || !equals_default)
    return result;

    // If tag value is default value, then we want every entity
    // in 'range' of the correct type, except those with a different tag value.
  
  Range all_tagged;
  result = get_entities( range, tag_handle, type, all_tagged );
  if (MB_SUCCESS != result) 
    return result;
  
    // get entities with a different tag value
  entities = subtract( all_tagged, entities );
    // get everything that does not have a different value
  entities = subtract( range, entities );
    // remove entities of the incorrect type, if any
  std::pair<Range::iterator,Range::iterator> p = entities.equal_range(type);
  entities.erase( entities.begin(), p.first );
  entities.erase( p.second, entities.end() );
    // merge with entities initially in range
  entities.merge( tmp_ents );
    
  return MB_SUCCESS;
  
}

ErrorCode TagServer::get_entities_with_tag_values( const Range &input_range,
                                                      const EntityType type,
                                                      const Tag *tags,
                                                      const void* const* values,
                                                      const int num_tags,
                                                      Range &entities,
                                                      const int condition) 
{
    // range should never come in empty
  assert(!input_range.empty());
  if (input_range.empty()) return MB_FAILURE;
  Range range = input_range;
  
  ErrorCode result;
  Range temp1;

  if (condition != Interface::INTERSECT &&
      condition != Interface::UNION)
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
    if (0 == it && condition == Interface::INTERSECT && entities.empty()) {
      entities = intersect( temp1, range);
    }

      // else if we're doing a running intersection, intersect this result (temp1)
      // with the running result (entities) into temp2, then move that to the running
      // result (entities)
    else if (condition == Interface::INTERSECT) {
      entities = intersect( entities, temp1);
      if (entities.empty()) return MB_SUCCESS;

        // also restrict the range at which we look; entities has already been 
        // intersected with range (through input to get_entities above) so just assign
      range = entities;
    }

      // else if we're doing a union, put these results (temp1) into the running 
      // result (entities)
    else if (condition == Interface::UNION) {
      entities.merge(temp1);
    }
  }

    // running result is in entities, where it should be
  return MB_SUCCESS;
}

ErrorCode TagServer::get_number_entities( const Tag tag_handle,
                                            unsigned long& num_entities )
{
  num_entities = 0;
  for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
    int type_num_entities;
    ErrorCode rval = get_number_entities( tag_handle, t, type_num_entities );
    if (MB_SUCCESS != rval) 
      return rval;
    num_entities += type_num_entities;
  }
  return MB_SUCCESS;
}

ErrorCode TagServer::get_number_entities( const Tag tag_handle, 
                                             const EntityType type,
                                             int& num_entities)
{
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      return mSparseData->get_number_entities(id, type, num_entities);
    case MB_TAG_DENSE:
      return sequenceManager->count_tagged_entities(id, type, num_entities);
    case MB_TAG_BIT:
      return mBitServer->get_number_entities(id, type, num_entities);
    case MB_TAG_MESH:
      return MB_TYPE_OUT_OF_RANGE;
  }
  return MB_TAG_NOT_FOUND;
}

ErrorCode TagServer::get_number_entities( const Range &range,
                                             const Tag tag_handle, 
                                             const EntityType type,
                                             int& num_entities)
{
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      return mSparseData->get_number_entities(range, id, type, num_entities);
    case MB_TAG_DENSE: {
      Range temp;
      ErrorCode rval = get_entities( range, tag_handle, type, temp );
      num_entities = temp.size();
      return rval;
    }
    case MB_TAG_BIT:
      return mBitServer->get_number_entities(range, id, type, num_entities);
    case MB_TAG_MESH:
      return MB_TYPE_OUT_OF_RANGE;
  }
  return MB_TAG_NOT_FOUND;
}

unsigned long TagServer::get_memory_use( Tag tag_handle ) const
{
  const TagInfo* tag_info = get_tag_info(tag_handle);
  if (!tag_info)
    return 0;

  unsigned long result = 0, tmp;
  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      mSparseData->get_memory_use( id, result, tmp );
      break;
    case MB_TAG_DENSE:
      sequenceManager->get_tag_memory_use( id, result, tmp );
      break;
    case MB_TAG_BIT:
      mBitServer->get_memory_use( id, result, tmp );
      break;
    case MB_TAG_MESH:
      break;
  }

    // add in size of entry in mTagTable
  result += sizeof(Tag) + sizeof(TagInfo) + 3*sizeof(void*);
  result += tag_info->default_value_size();
  result += tag_info->get_mesh_value_size();
  result += tag_info->get_name().size();
  return result;
}

ErrorCode TagServer::get_memory_use( Tag tag_handle,
                                       unsigned long& total,
                                       unsigned long& per_entity ) const
{
  const TagInfo* tag_info = get_tag_info(tag_handle);
  if(!tag_info)
    return MB_TAG_NOT_FOUND;

  TagId id = ID_FROM_TAG_HANDLE(tag_handle);
  switch (PROP_FROM_TAG_HANDLE(tag_handle)) {
    case MB_TAG_SPARSE:
      mSparseData->get_memory_use( id, total, per_entity );
      break;
    case MB_TAG_DENSE:
      sequenceManager->get_tag_memory_use( id, total, per_entity );
      break;
    case MB_TAG_BIT:
      mBitServer->get_memory_use( id, total, per_entity );
      break;
    case MB_TAG_MESH:
      break;
  }
  
    // size of entry in mTagTable map
  total += sizeof(Tag) + sizeof(TagInfo) + 3*sizeof(void*);
  total += tag_info->default_value_size();
  total += tag_info->get_mesh_value_size();
  total += tag_info->get_name().size();
  
  return MB_SUCCESS;
}
    

} // namespace moab

#ifdef TEST

#include <iostream>
#include <time.h>
#include <assert.h>
#include <memory.h>

#include "TagServer.hpp"

using namespace moab;

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
  Tag my_tag_handle;
  unsigned char tag_data[MAX_TAG_DATA_SIZE] = {0};
  unsigned char tag_size;
  unsigned int entity_type;
  Tag tag_handle = (Tag)3245;


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

      Tag tmp_tag_handle = tag_server.get_handle(tag_name.c_str());

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
      Range tmp_range;
      tag_server.get_entities_with_tag_value( TYPE_FROM_HANDLE(200), my_tag_handle, tag_data, tmp_range );
      assert( tmp_range.size() == 1 && tmp_range.begin() == 200 );

      for( int i=0; i<SET_TAG_LOOPS; i++ )
      {
        EntityHandle id = (MAX_ID*rand())/(RAND_MAX+1);
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
      Tag tmp = tag_server.get_handle(tag_name.c_str()); 
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
