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
 */

#ifndef DENSE_TAG_COLLECTIONS_HPP
#define DENSE_TAG_COLLECTIONS_HPP

#ifndef IS_BUILDING_MB
#error "DenseTagCollections.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBRange.hpp"

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
  DensePage() : mByteArray(NULL) {}
  
  // default destructor
  ~DensePage()
  {  
    if(mByteArray)
      delete [] mByteArray;
  }
  
  //! get the bytes from a page
  MBErrorCode get_bytes(int offset, int num_bytes_per_flag, 
      void* data);
  
  //! set the bytes in a page
  MBErrorCode set_bytes(int offset, int num_bytes_per_flag, 
      const void* default_data, const void* data);
  
  //! remove the bytes in page only if space has been allocated
  MBErrorCode remove_data(int offset, int num_bytes_per_flag, 
      const void* default_data);

  //! return whether has data or not
  bool has_data() const { return mByteArray == NULL ? false : true; }


private:
  //! byte array,  uses lazy allocation
  unsigned char* mByteArray;
  
  //! don't allow copying of these pages
  DensePage(const DensePage&) 
    : mByteArray(NULL) 
  { 
    // not even by self
    assert(0);
  }
  
  //! don't allow copying of these pages
  DensePage& operator=(const DensePage&) 
  { 
    // not even by self
    assert(0);
    return *this;
  }
  
};

/*! get the bytes from a page
    takes offset into the page
    takes how many bytes to get
    return data
*/
inline MBErrorCode DensePage::get_bytes(int offset, int num_bytes_per_flag, void* data)
{
  // if no memory has been allocated, get the default value
  if(!mByteArray)
  {
    return MB_FAILURE;
  }

  unsigned char* data_to_copy = mByteArray + offset*num_bytes_per_flag;
  memcpy(data, data_to_copy, num_bytes_per_flag);

  return MB_SUCCESS;
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
    mByteArray = new unsigned char [mPageSize];
    if(!default_data)
    {
      memset(mByteArray, 0, mPageSize);
    }
    else
    {
      unsigned char* byte_array = mByteArray;
      unsigned char* byte_array_end = byte_array + mPageSize;
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
  DensePageGroup(int bytes_per_flag)
      : mBytesPerFlag(bytes_per_flag)
  {
    //compute the offset factor based on the number of bytes for each entity
    mOffsetFactor = compute_offset_factor(bytes_per_flag);
    mDensePagesSize=0;
  }

  // default destructor
  ~DensePageGroup() 
  { 
    //clear(); 
    // delete each page
    for (std::vector<DensePage*>::iterator iter = mDensePages.begin(); iter != mDensePages.end(); ++iter)
      delete *iter;

    // clean out the vector of pointers
    mDensePages.clear(); 
  }

  //! get data from byte pages
  MBErrorCode get_data(MBEntityHandle handle, void* data);

  //! set data in byte pages
  MBErrorCode set_data(MBEntityHandle handle, const void* default_data, const void* data);

  //! remove data associated to an entity
  MBErrorCode remove_data(MBEntityHandle handle, const void* default_value);

  //! get number of entities of type
  MBErrorCode get_number_entities(MBEntityType type, int& num_entities);

  //! get the entities
  MBErrorCode get_entities(MBEntityType type, MBRange& entities);
  
  //! get the entities with a value
  MBErrorCode get_entities_with_tag_value(const MBEntityType type, const void* value, 
                                          MBRange &entities,
                                          const bool equals_default);

    //! return true if this page group contains this entity, false otherwise
  bool contains(const MBEntityHandle entity) const;
  

  int tag_size() const { return mBytesPerFlag; }

private:

  // compute offset factor to use when computing which page to index into
  int compute_offset_factor(int num_bytes)
  {
    return (DensePage::mPageSize) / num_bytes;
  }
  
  //!  number of bytes for each entity  
  unsigned short mBytesPerFlag;
  
  //! offset factor used when computing which page to jump to
  //! can think of this as the number of entities per page
  unsigned short mOffsetFactor;

  //! cached size of mDensePages
  unsigned int mDensePagesSize;

  //! vector of dense byte pages
  std::vector<DensePage*> mDensePages;
  
  //! don't allow copy of this class
  DensePageGroup& operator=(const DensePageGroup&)
  {
    assert(0);
    return *this;
  }

  //! don't allow copy of this class
  DensePageGroup(const DensePageGroup&)
  {
    assert(0);
  }

  //! don't allow bare constructor
  DensePageGroup()
  {
    assert(0);
  }

};

/*! get the data from pages
    takes entity handle
    return the data
*/
inline MBErrorCode DensePageGroup::get_data(MBEntityHandle handle, void* data)
{
  // strip off the entity type
  unsigned int tmp_handle = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  unsigned int which_page = tmp_handle / mOffsetFactor;
  
  // if the page isn't there, just return failure
  if(which_page >= mDensePagesSize)
    return MB_TAG_NOT_FOUND;

  // return data from page
  return mDensePages[which_page]->get_bytes( (tmp_handle - ( which_page * mOffsetFactor )), mBytesPerFlag, data);
}

inline bool DensePageGroup::contains(const MBEntityHandle handle) const
{
  // strip off the entity type
  unsigned int entity_id = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  unsigned int which_page = entity_id / mOffsetFactor;
  
  // if the page isn't there, the entity isn't assigned this tag
  return (which_page >= mDensePagesSize) ? false : true;
}

/*! set the data in pages
    takes entity handle
    takes the data to set
*/
inline MBErrorCode DensePageGroup::set_data(MBEntityHandle handle, const void* default_data, const void* data)
{
  // strip off the entity type
  handle = ID_FROM_HANDLE(handle);
  // figure out which page to jump to
  unsigned int which_page = handle / mOffsetFactor;

  // if the page doesn't exist, make one
  if(which_page >= mDensePagesSize)
  {
    for(int j= which_page - mDensePagesSize +1; j--;)
      mDensePages.push_back(new DensePage());
    mDensePagesSize = mDensePages.size();
  }

  // return data in page
  return mDensePages[which_page]->set_bytes( (handle - ( which_page * mOffsetFactor )), 
      mBytesPerFlag, default_data, data);
}


/*! set the data in pages
    takes entity handle
    takes the default data to set
*/
inline MBErrorCode DensePageGroup::remove_data(MBEntityHandle handle, const void* default_data)
{
  // strip off entity type
  handle = ID_FROM_HANDLE(handle);
  // find out which page to jump to
  unsigned int which_page = handle / mOffsetFactor;
  
  // if the page doesn't exist, return
  if(which_page >= mDensePagesSize)
    return MB_SUCCESS;
 
  // try to clean out data
  return mDensePages[which_page]->remove_data( (handle - ( which_page * mOffsetFactor )), 
      mBytesPerFlag, default_data);
}

//! get the entities
inline MBErrorCode DensePageGroup::get_entities(MBEntityType type, MBRange& entities)
{
  std::vector<DensePage*>::iterator iter;
  int dum =0;
  MBEntityHandle handle = CREATE_HANDLE(type, 0, dum);
  int first_time = 1; // Don't want zero-ID handle at start of range.
  for(iter = mDensePages.begin(); iter < mDensePages.end(); ++iter)
  {
    if (*iter && (*iter)->has_data())
    {
      entities.insert( handle + first_time, handle + mOffsetFactor - 1 );
    }
    first_time = 0;
    handle += mOffsetFactor;
  }
  return MB_SUCCESS;
}

//! get number of entities of type
inline MBErrorCode DensePageGroup::get_number_entities(MBEntityType , int& entities)
{
  std::vector<DensePage*>::iterator iter;
  for(iter = mDensePages.begin(); iter < mDensePages.end(); ++iter)
  {
    if(*iter)
    {
      if((*iter)->has_data())
      {
        entities += mOffsetFactor;
      }        
    }
  }
  return MB_SUCCESS;
}

//! DenseTagSuperCollection class provides constant time 
//! lookup for data tagged on entities
class DenseTagSuperCollection
{
public:
  //! default constructor
  DenseTagSuperCollection()
  {
    mDensePageGroupsSize = 0;
  }
  //! default destructor
  ~DenseTagSuperCollection()
  {
    // clean things out
    clear();
  }

  void reset_data();

  //! return an available tag id for use
  MBErrorCode reserve_tag_id(int data_size, const void* default_data, MBTagId& tag_id);
  //! release a tag id for reuse
  MBErrorCode release_tag_id(MBTagId tag_id);
  //! get the tag size
  int tag_size(const MBTagId tag_id) const;
  //! set the data associated with an entity handle
  MBErrorCode set_data(MBTagId tag_id, MBEntityHandle handle, const void* data);
  //! get the data associated with an entity handle
  MBErrorCode get_data(const MBTagId tag_id, const MBEntityHandle handle, void* data);
  //! remove/clean out data associated with an entity handle, only if memory has been allocated
  MBErrorCode remove_data(MBTagId tag_id, MBEntityHandle handle, const void* default_data);

  //! get the entities with a tag
  MBErrorCode get_number_entities(const MBTagId tag_id, const MBEntityType type, int& num_entities);

  //! get the entities with a tag
  MBErrorCode get_number_entities(const MBRange &range,
                                   const MBTagId tag_id, const MBEntityType type, int& num_entities);

  //! get the entities with a tag
  MBErrorCode get_entities(const MBTagId tag_id, const MBEntityType type, MBRange& entities);
  
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
  
private:

  //! clean things out
  void clear() 
  {
    for(int i = 0; i<(int)MBMAXTYPE; i++)
    {
      for (std::vector<DensePageGroup*>::iterator iter = mDensePageGroups[i].begin();
          iter != mDensePageGroups[i].end();
          ++iter)
      {
        if(*iter == NULL)
          continue;
        delete *iter;
      }

      mDensePageGroups[i].clear();
    }
  }


  //! dense pages are indexed by tag id and entity type
  std::vector< DensePageGroup* > mDensePageGroups[MBMAXTYPE];
  unsigned long mDensePageGroupsSize;
  std::vector< void* > mDefaultData;

};

/*! get some data based on a tag id and handle */
inline MBErrorCode DenseTagSuperCollection::get_data(const MBTagId tag_id, 
                                                      const MBEntityHandle handle, void* data)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  MBErrorCode result = 
    mDensePageGroups[TYPE_FROM_HANDLE(handle)][tag_id]->get_data(handle, data);
  
  if(result == MB_FAILURE || result == MB_TAG_NOT_FOUND)
  {
    if(mDefaultData[tag_id])
    {
      memcpy(data, mDefaultData[tag_id], (*mDensePageGroups)[tag_id]->tag_size());
      return MB_SUCCESS;
    }
    else
      return result;
  }
  
  return MB_SUCCESS;
}

/*! set some data based on a tag id and handle */
inline MBErrorCode DenseTagSuperCollection::set_data(MBTagId tag_id, 
    MBEntityHandle handle, const void* data)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_TAG_NOT_FOUND;

  return mDensePageGroups[TYPE_FROM_HANDLE(handle)][tag_id]->set_data(handle, mDefaultData[tag_id], data);
}

/*! set some data based on a tag id and handle only if memory has been allocated*/
inline MBErrorCode DenseTagSuperCollection::remove_data(MBTagId tag_id, 
    MBEntityHandle handle, const void* default_data)
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return MB_SUCCESS;

  return mDensePageGroups[TYPE_FROM_HANDLE(handle)][tag_id]->remove_data(handle, default_data);
}

inline int DenseTagSuperCollection::tag_size(const MBTagId tag_id) const 
{
  if(tag_id >= mDensePageGroupsSize || (*mDensePageGroups)[tag_id] == NULL)
    return 0;

  return (*mDensePageGroups)[tag_id]->tag_size();
}

    //! get all tags defined on an entity
inline MBErrorCode DenseTagSuperCollection::get_tags(const MBEntityHandle entity,
                                                      std::vector<MBTag> &tags) 
{
    // get the tags defined for this type
  MBEntityType this_type = TYPE_FROM_HANDLE(entity);

  for (long i = 0; i < (long) mDensePageGroups[this_type].size(); i++) {
    if (mDensePageGroups[this_type][i] != NULL &&
        mDensePageGroups[this_type][i]->contains(entity))
      tags.push_back(TAG_HANDLE_FROM_ID(i, MB_TAG_DENSE));
  }
  return MB_SUCCESS;
}


#endif


