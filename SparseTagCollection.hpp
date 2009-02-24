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
 * Filename   :     SparseTagCollection.hpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/



#ifndef SPARSE_TAG_COLLECTION_HPP
#define SPARSE_TAG_COLLECTION_HPP

#ifndef IS_BUILDING_MB
#error "SparseTagCollection.hpp isn't supposed to be included into an application"
#endif

#ifdef WIN32
#pragma warning(disable : 4786)
#endif


#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#ifdef HAVE_UNORDERED_MAP
# include STRINGIFY(HAVE_UNORDERED_MAP)
#else
# include <map>
#endif
#include <vector>

#include "MBTypes.h"
#include "MBInternals.hpp"
#include "MBRange.hpp"
#include "TagInfo.hpp"

//! allocator for tag data
class SparseTagDataAllocator
{
public:
  //! constructor
  SparseTagDataAllocator(){}
  //! destructor
  ~SparseTagDataAllocator(){}
  //! allocates memory of size and returns pointer
  void* allocate(size_t data_size) { return malloc(data_size); }
  //! frees the memory
  void destroy(void* p){ free(p); }
};


//! collection of tag data associated with entity ids
class SparseTagCollection
{
public:
  
  //! constructor takes tag data size
  SparseTagCollection(int data_size);
  
  //! destructor
  ~SparseTagCollection();
  
  //! set the tag data for an entity id
  //!\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called for 
  //!      variable-length tag.
  MBErrorCode set_data(const MBEntityHandle entity_handle, const void* data);

  //! get the tag data for an entity id
  //!\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called for 
  //!      variable-length tag.
  MBErrorCode get_data(const MBEntityHandle entity_handle, void* data);
  
  //! set variable-length tag data for an entity id
  //!
  //!\NOTE If called for fixed-length tag, size must be either zero or the tag size.
  //!
  //!\NOTE If called with zero size for a variable-length tag, is equivalent
  //!      to remove_data().
  MBErrorCode set_data(const MBEntityHandle entity_handle, const void* data, int length);

  //! get the variable-length data for an entity id
  MBErrorCode get_data(const MBEntityHandle entity_handle, const void*& data, int& length);

  //! removes the data
  MBErrorCode remove_data(const MBEntityHandle entity_handle);

  //! get number of entities of type
  MBErrorCode get_number_entities(MBEntityType type, int& num_entities);
  
  //! get number of entities
  unsigned long get_number_entities()
    { return mData.size(); }

  //! gets all entity handles that match a type and tag
  MBErrorCode get_entities(MBEntityType type, MBRange &entities);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities(MBRange &entities) const;

  //! gets all entity handles that match a type, tag, tag_value
  MBErrorCode get_entities_with_tag_value( const TagInfo& info,
                                           MBEntityType type, 
                                           MBRange &entities, 
                                           const void* tag_value,
                                           int value_size);

  //! if this collection contains this entity, return true, otherwise false
  bool contains(const MBEntityHandle entity) const;
  
  int tag_size() const { return mDataSize; }

protected:
  
  //! hidden constructor
  SparseTagCollection(){}
  
  //! size of the data
  int mDataSize;

  //! allocator for this collection
  SparseTagDataAllocator mAllocator;

  //! map of entity id and tag data
#ifdef HAVE_UNORDERED_MAP
  typedef UNORDERED_MAP_NS::unordered_map<MBEntityHandle,void*> myMapType;
#else
  typedef std::map<MBEntityHandle /*entity_handle*/ , void* /*data*/ > myMapType;
#endif

  myMapType mData;
};

inline bool SparseTagCollection::contains(const MBEntityHandle entity) const
{
  return (mData.find(entity) == mData.end() ? false : true);
}

inline MBErrorCode SparseTagCollection::get_entities(MBRange &entities) const 
{
  for (myMapType::const_iterator mit = mData.begin(); mit != mData.end(); mit++) 
    entities.insert((*mit).first);

  return MB_SUCCESS;
}


#endif //SPARSE_TAG_COLLECTION_HPP




