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


/*  Dense tag storage for MB
 *
 *  File   :      DenseTagCollections.hpp
 *  Creator:      Clinton Stimpson
 *  Date   :      10-28-2002
 *
 * DenseTagSuperCollection: contains a vector of DensePageGroup objects
 * DensePageGroup: vector of vectors of DensePage objects, one vector per 
 *   entity type, blocking possible handles in the type
 * DensePage: block of values corresponding to a block of handles of a given
 *   type; may or may not have data allocated
 */

#ifndef DENSE_TAG_COLLECTIONS_HPP
#define DENSE_TAG_COLLECTIONS_HPP

#ifndef IS_BUILDING_MB
#error "DenseTagCollections.hpp isn't supposed to be included into an application"
#endif

#include "MBInternals.hpp"
#include "MBTypes.h"
#include "MBRange.hpp"

#include <vector>
#include <assert.h>



//! dense page class
/*! This class stores a page of memory for storing dense data
*/
class DensePage
{

public:

  //! the page size
  static const int mPageSize;
  
  // default constructor
  DensePage() : mByteArray(0) {}
  
  // default destructor
  ~DensePage() { delete [] mByteArray; }
  
  //! get the bytes from a page
  MBErrorCode get_bytes(int offset, int num_bytes_per_flag, 
      void* data) const;
  
  //! set the bytes in a page
  MBErrorCode set_bytes(int offset, int num_bytes_per_flag, 
      const void* default_data, const void* data);
  
  //! remove the bytes in page only if space has been allocated
  MBErrorCode remove_data(int offset, int num_bytes_per_flag, 
      const void* default_data);

  //! return whether has data or not
  bool has_data() const { return mByteArray != 0; }

    //! do a memcmp on one of the values in this page
  bool memcmp(int offset, int num_bytes_per_flag, const void *data) const;

  //!std::auto_ptr-style behavior
  DensePage( const DensePage& other )
    : mByteArray(other.mByteArray) 
    { const_cast<DensePage&>(other).mByteArray = 0; }
  DensePage& operator=( const DensePage& other ) 
    {
      assert(!mByteArray);
      delete [] mByteArray;
      mByteArray = other.mByteArray; 
      const_cast<DensePage&>(other).mByteArray = 0;
      return *this; 
    }

private:

  //! byte array,  uses lazy allocation
  unsigned char* mByteArray;
};

/*! get the bytes from a page
    takes offset into the page
    takes how many bytes to get
    return data
*/
inline MBErrorCode DensePage::get_bytes(int offset, int num_bytes_per_flag, void* data) const
{
  // if no memory has been allocated, get the default value
  if(!mByteArray)
  {
    return MB_TAG_NOT_FOUND;
  }

  unsigned char* data_to_copy = mByteArray + offset*num_bytes_per_flag;
  memcpy(data, data_to_copy, num_bytes_per_flag);

  return MB_SUCCESS;
}

//! do a memcmp on one of the values in this page
inline bool DensePage::memcmp(int offset, int num_bytes_per_flag, const void *data) const
{
  // if no memory has been allocated, get the default value
  assert(mByteArray);

  return ::memcmp(mByteArray + offset*num_bytes_per_flag, data, num_bytes_per_flag);
}

/*! set the bytes in a page
    takes byte offset into the page
    takes how many bytes to set
    takes the data to set
*/
inline MBErrorCode DensePage::set_bytes(int offset, 
    int num_bytes_per_flag, const void* default_data, const void* data)
{
  // if memory hasn't been allocated, allocate it and zero the memory
  if(!mByteArray)
  {
    mByteArray = new unsigned char [mPageSize*num_bytes_per_flag];
    if (!mByteArray) 
      return MB_MEMORY_ALLOCATION_FAILED;
    
    if(!default_data)
    {
      memset(mByteArray, 0, mPageSize*num_bytes_per_flag);
    }
    else
    {
      unsigned char* byte_array = mByteArray;
      unsigned char* byte_array_end = byte_array + mPageSize*num_bytes_per_flag;
      for(; byte_array < byte_array_end; byte_array += num_bytes_per_flag)
      {
        memcpy(byte_array, default_data, num_bytes_per_flag);
      }
    }
  }

  unsigned char* data_to_copy = mByteArray + offset*num_bytes_per_flag;
  memcpy(data_to_copy, data, num_bytes_per_flag);
  return MB_SUCCESS;
  
}

/*! remove data if memory has been allocated
    takes byte offset
    takes number of bytes to set
    takes the data to set
*/
inline MBErrorCode DensePage::remove_data(int offset, 
    int num_bytes_per_flag, const void* default_data)
{
  if(mByteArray)
  {
    unsigned char* data_to_copy = mByteArray + offset*num_bytes_per_flag;
    if(default_data)
      memcpy(data_to_copy, default_data, num_bytes_per_flag);
    else
      memset(data_to_copy, 0, num_bytes_per_flag);
  }
  return MB_SUCCESS;
}

//! class which is a collection of byte pages
class DensePageGroup
{
public:
  // default constructor
  DensePageGroup(int bytes_per_flag, const void* default_value)
      : mBytesPerFlag(bytes_per_flag),
        mDefaultValue(0)
  {
    if (default_value) {
      mDefaultValue = malloc(bytes_per_flag);
      memcpy( mDefaultValue, default_value, bytes_per_flag );
    }
  }

  // default destructor
  ~DensePageGroup() { free(mDefaultValue); }

  //! get data from byte pages
  MBErrorCode get_data(MBEntityHandle handle, void* data) const;

  //! set data in byte pages
  MBErrorCode set_data(MBEntityHandle handle, const void* data);

  //! remove data associated to an entity
  MBErrorCode remove_data(MBEntityHandle handle);

  //! get number of entities of type
  MBErrorCode get_number_entities(MBEntityType type, int& num_entities) const;

  //! get the entities
  MBErrorCode get_entities(MBEntityType type, MBRange& entities) const;
  
  //! get the entities
  MBErrorCode get_entities(MBRange& entities) const;
  
  //! get the entities with a value
  MBErrorCode get_entities_with_tag_value(const MBEntityType type, const void* value, 
                                          MBRange &entities) const;

    //! return true if this page group contains this entity, false otherwise
  bool contains(const MBEntityHandle entity) const;
  

  int tag_size() const { return mBytesPerFlag; }
  
  const void* get_default_value() const { return mDefaultValue; }
  
  MBErrorCode get_memory_use( unsigned long& total,
                              unsigned long& per_entity ) const;

private:
  
  //!  number of bytes for each entity  
  unsigned short mBytesPerFlag;
  
  //! default value for tag
  void* mDefaultValue;

  //! vectors of dense byte pages
  //! allocate more than MBMAXTYPE vectors - instead allocate for all possible
  //! values of the type bits within a handle so we don't need to check the
  //! size of the type all the time - just let in index an empty list and just
  //! return not-found if it is out of bounds.
  std::vector<DensePage> mDensePages[1<<MB_TYPE_WIDTH];
  
  //! don't allow copy of this class
  DensePageGroup& operator=(const DensePageGroup&);

  //! don't allow copy of this class
  DensePageGroup(const DensePageGroup&);
};

/*! get the data from pages
    takes entity handle
    return the data
*/
inline MBErrorCode DensePageGroup::get_data(MBEntityHandle handle, void* data) const
{
  // strip off the entity type
  const MBEntityType type = TYPE_FROM_HANDLE( handle );
  const MBEntityID entity_id = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  const MBEntityID which_page = entity_id / DensePage::mPageSize;
  const unsigned int offset = entity_id % DensePage::mPageSize;
  
  std::vector<DensePage>::const_iterator page = mDensePages[type].begin() + which_page;
  if (page >= mDensePages[type].end())
    return MB_TAG_NOT_FOUND;
  
  return page->get_bytes( offset, mBytesPerFlag, data);
}

inline bool DensePageGroup::contains(const MBEntityHandle handle) const
{
  // strip off the entity type
  const MBEntityType type = TYPE_FROM_HANDLE( handle );
  const MBEntityID entity_id = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  const MBEntityID which_page = entity_id / DensePage::mPageSize;

  std::vector<DensePage>::const_iterator page = mDensePages[type].begin() + which_page;
  return page < mDensePages[type].end() && page->has_data();
}

/*! set the data in pages
    takes entity handle
    takes the data to set
*/
inline MBErrorCode DensePageGroup::set_data(MBEntityHandle handle, const void* data)
{
  // strip off the entity type
  const MBEntityType type = TYPE_FROM_HANDLE( handle );
  const MBEntityID entity_id = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  const MBEntityID which_page = entity_id / DensePage::mPageSize;
  const unsigned int offset = entity_id  % DensePage::mPageSize;

  std::vector<DensePage>::iterator page = mDensePages[type].begin() + which_page;
  if (page >= mDensePages[type].end()) {
    mDensePages[type].resize(mDensePages[type].size() + which_page + 1);
    page = mDensePages[type].begin() + which_page;
  }

  // return data in page
  return page->set_bytes( offset, mBytesPerFlag, mDefaultValue, data);
}


/*! set the data in pages
    takes entity handle
    takes the default data to set
*/
inline MBErrorCode DensePageGroup::remove_data(MBEntityHandle handle)
{
  // strip off the entity type
  const MBEntityType type = TYPE_FROM_HANDLE( handle );
  const MBEntityID entity_id = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  const MBEntityID which_page = entity_id / DensePage::mPageSize;
  const unsigned int offset = entity_id  % DensePage::mPageSize;
  
  std::vector<DensePage>::iterator page = mDensePages[type].begin() + which_page;
  // if the page doesn't exist, return
  // Return value changed from MB_SUCCESS to MB_FAILURE - j.k. 2006-8-23
  if (page >= mDensePages[type].end())
    return MB_FAILURE;
 
  // try to clean out data
  return page->remove_data( offset, mBytesPerFlag, mDefaultValue );
}

//! get the entities
inline MBErrorCode DensePageGroup::get_entities(MBEntityType type, MBRange& entities) const
{
  std::vector<DensePage>::const_iterator iter;
  const std::vector<DensePage>::const_iterator end = mDensePages[type].end();
  int dum =0;
  MBEntityHandle handle = CREATE_HANDLE(type, 0, dum);
  MBEntityID first_time = MB_START_ID; // Don't want zero-ID handle at start of range.
  MBRange::iterator insert_pos = entities.begin();
  for(iter = mDensePages[type].begin(); iter != end; ++iter, handle += DensePage::mPageSize)
  {
    if (iter->has_data())
      insert_pos = entities.insert( insert_pos, handle + first_time, handle + DensePage::mPageSize - 1 );
    first_time = 0;
  }
  return MB_SUCCESS;
}

//! get the entities
inline MBErrorCode DensePageGroup::get_entities(MBRange& entities) const
{
  std::vector<DensePage>::const_iterator iter;
  int dum =0;
  for (MBEntityType type = MBENTITYSET; type >= MBVERTEX; type--) {
    const std::vector<DensePage>::const_iterator end = mDensePages[type].end();
    MBEntityHandle handle = CREATE_HANDLE(type, 0, dum);
    MBEntityID first_time = MB_START_ID; // Don't want zero-ID handle at start of range.
    MBRange::iterator insert_pos = entities.begin();
    for(iter = mDensePages[type].begin(); iter != end; ++iter, handle += DensePage::mPageSize)
    {
      if (iter->has_data())
        insert_pos = entities.insert( insert_pos, handle + first_time, handle + DensePage::mPageSize - 1 );
      first_time = 0;
    }
  }
  
  return MB_SUCCESS;
}

//! get number of entities of type
inline MBErrorCode DensePageGroup::get_number_entities(MBEntityType type, int& entities) const
{
  entities = 0;
  std::vector<DensePage>::const_iterator iter;
  const std::vector<DensePage>::const_iterator end = mDensePages[type].end();
  MBEntityID first_time = MB_START_ID;
  for(iter = mDensePages[type].begin(); iter != end; ++iter) {
    if(iter->has_data())
      entities += DensePage::mPageSize - first_time;
    first_time = 0;
  }
  return MB_SUCCESS;
}

//! DenseTagSuperCollection class provides constant time 
//! lookup for data tagged on entities
class DenseTagSuperCollection
{
public:
  //! default constructor
  DenseTagSuperCollection() {}
  
  //! default destructor
  ~DenseTagSuperCollection()
  {
    // clean things out
    clear();
  }

  void reset_data();

  //! allocate new tag id
  MBErrorCode reserve_tag_id(int data_size, const void* default_data, MBTagId tag_id);
  //! release a tag id for reuse
  MBErrorCode release_tag_id(MBTagId tag_id);
  //! get the tag size
  int tag_size(const MBTagId tag_id) const;
  //! set the data associated with an entity handle
  MBErrorCode set_data(MBTagId tag_id, MBEntityHandle handle, const void* data);
  //! get the data associated with an entity handle
  MBErrorCode get_data(const MBTagId tag_id, const MBEntityHandle handle, void* data);
  //! remove/clean out data associated with an entity handle, only if memory has been allocated
  MBErrorCode remove_data(MBTagId tag_id, MBEntityHandle handle);

  //! get the entities with a tag
  MBErrorCode get_number_entities(const MBTagId tag_id, const MBEntityType type, int& num_entities);

  //! get the entities with a tag
  MBErrorCode get_number_entities(const MBRange &range,
                                   const MBTagId tag_id, const MBEntityType type, int& num_entities);

  //! get the entities with a tag
  MBErrorCode get_entities(const MBTagId tag_id, const MBEntityType type, MBRange& entities);
  
  //! get the entities with a tag
  MBErrorCode get_entities(const MBTagId tag_id, MBRange& entities);
  
  //! get the entities with a tag
  MBErrorCode get_entities(const MBRange &range,
                            const MBTagId tag_id, const MBEntityType type, MBRange& entities);
  
  //! get the entities with a value
  MBErrorCode get_entities_with_tag_value(const MBTagId tag_id, const MBEntityType type, 
                                           MBRange &entities, const void* value);
  
  //! get the entities with a value
  MBErrorCode get_entities_with_tag_value(const MBRange &range,
                                           const MBTagId tag_id, const MBEntityType type, 
                                           MBRange &entities, const void* value);

    //! get all tags defined on an entity
  MBErrorCode get_tags(const MBEntityHandle entity,
                       std::vector<MBTag> &tags);
  
  MBErrorCode get_memory_use( MBTagId tag_id,
                              unsigned long& total,
                              unsigned long& per_entity );
private:

  //! clean things out
  void clear() 
  {
    for (std::vector<DensePageGroup*>::iterator i = mDensePageGroups.begin();
         i != mDensePageGroups.end(); ++i) 
      delete *i;
    mDensePageGroups.clear();
  }


  //! dense pages are indexed by tag id
  std::vector< DensePageGroup* > mDensePageGroups;
};

/*! get some data based on a tag id and handle */
inline MBErrorCode DenseTagSuperCollection::get_data(const MBTagId tag_id, 
                                                      const MBEntityHandle handle, void* data)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;
  
  return (*group)->get_data( handle, data );
}

/*! set some data based on a tag id and handle */
inline MBErrorCode DenseTagSuperCollection::set_data(MBTagId tag_id, 
    MBEntityHandle handle, const void* data)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->set_data(handle, data);
}

/*! set some data based on a tag id and handle only if memory has been allocated*/
inline MBErrorCode DenseTagSuperCollection::remove_data(MBTagId tag_id, MBEntityHandle handle )
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->remove_data(handle);
}

inline int DenseTagSuperCollection::tag_size(const MBTagId tag_id) const 
{
  std::vector<DensePageGroup*>::const_iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->tag_size();
}

    //! get all tags defined on an entity
inline MBErrorCode DenseTagSuperCollection::get_tags(const MBEntityHandle entity,
                                                      std::vector<MBTag> &tags) 
{
  for (unsigned tagid = 0; tagid < mDensePageGroups.size(); ++tagid) 
    if (mDensePageGroups[tagid] && mDensePageGroups[tagid]->contains(entity))
      tags.push_back( TAG_HANDLE_FROM_ID( tagid, MB_TAG_DENSE ) );
  return MB_SUCCESS;
}

// get the entities with a tag
inline MBErrorCode DenseTagSuperCollection::get_entities(const MBTagId tag_id, 
                                                         MBRange& entities)
{
  std::vector<DensePageGroup*>::iterator group = mDensePageGroups.begin() + tag_id;
  if (group >= mDensePageGroups.end() || !*group)
    return MB_TAG_NOT_FOUND;

  return (*group)->get_entities(entities);
}

#endif


