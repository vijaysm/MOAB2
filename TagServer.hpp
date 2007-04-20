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



#ifndef TAG_SERVER_HPP
#define TAG_SERVER_HPP

#ifndef IS_BUILDING_MB
#error "TagServer.hpp isn't supposed to be included into an application"
#endif

#include <string>
#include <vector>
#include <map>

#include "MBTypes.h"

class MBRange;
class DenseTagSuperCollection;
class SparseTagSuperCollection;
class MBBitServer;

// ! stores information about a tag
class TagInfo
{
public:

  //! constructor
  TagInfo() : mTagName(""), mDataSize(0), mDefaultValue(NULL) {}
  
  //! destructor
  ~TagInfo();
  
  //! copy constructor
  TagInfo(const TagInfo&);

  //! constructor that takes all parameters
  TagInfo(const char *, int, MBDataType type, const void *);

  //! assignment operator
  TagInfo &operator=(const TagInfo &rhs);
  
  //! set the name of the tag
  void set_name( const std::string& name) { mTagName = name; }

  //! set the size of the data
  void set_size( const int size ) { mDataSize = size; }

  //! get the name of the tag
  const std::string& get_name() const { return mTagName; }

  //! get the size of the data
  int get_size() const { return mDataSize; }

    //! get the default data
  const void *default_value() const  { return mDefaultValue;}
  
  inline MBDataType get_data_type() const     { return dataType; }
  
  inline void set_data_type( MBDataType t )   { dataType = t; }

  static int size_from_data_type( MBDataType t );

private:    

  //! stores the tag name
  std::string mTagName;

  //! stores the size of the data for this tag
  unsigned short mDataSize;

  //! stores the default data, if any
  unsigned char *mDefaultValue;
  
  //! type of tag data
  MBDataType dataType;

};


//! SparseTagServer class which associates tag data with entities
class TagServer
{
public:
  //! constructor
  TagServer();
  //! destructor
  virtual ~TagServer();

  //! add a tag
  MBErrorCode add_tag(const char *tag_name, 
                       const int data_size,
                       const MBTagType storage,
                       const MBDataType data_type,
                       MBTag &tag_handle,
                       const void *default_data = NULL);

  // tag name is the name of the tag
  // entity_type is the type of entity (in implementation - it means
  // tag data is stored in a different namespace)
  // entity_handle is the entity to which data is tagged to
  // tag_usage defines which storage mechanism to use
  // data is the actual data 
  // data_size is the size of the data in bytes

  //! remove a tag from this tag server
  MBErrorCode remove_tag(const MBTag tag_handle);

  //! resets all tag data associated with this handle back to default data or NULL
  MBErrorCode reset_data(MBEntityHandle entity_handle);

  //! cleans out all data tagged on all entities
  MBErrorCode reset_all_data();

  //! set the value of a tag
  MBErrorCode set_data(const MBTag tag_handle, const MBEntityHandle entity_handle, const void* data );
  
  MBErrorCode set_data(const MBTag tag_handle, const MBEntityHandle* entity_handles, const int num_entities, const void* data );
  
  MBErrorCode set_data(const MBTag tag_handle, const MBRange& entity_handles, const void* data );

  //! get the value of a tag
  MBErrorCode get_data(const MBTag tag_handle, const MBEntityHandle entity_handle, void* data );
  
  MBErrorCode get_data(const MBTag tag_handle, const MBEntityHandle* entity_handles, const int num_ents, void* data );
  
  MBErrorCode get_data(const MBTag tag_handle, const MBRange& entity_handles, void* data );
  
  //! set the value of a tag
  MBErrorCode set_bits(const MBTag tag_handle, const MBEntityHandle entity_handle, unsigned char data );

  //! get the value of a tag
  MBErrorCode get_bits(const MBTag tag_handle, const MBEntityHandle entity_handle, unsigned char& data );

  //! remove the tag data on an entity
  MBErrorCode remove_data( const MBTag tag_handle, const MBEntityHandle entity_handle );

  //! gets all entity handles that match a type and tag
  MBErrorCode get_entities( const MBTag tag_handle, 
                             const MBEntityType type,
                             MBRange &entities);

  //! gets all entity handles that match a type and tag
  MBErrorCode get_entities( const MBRange &input_range,
                             const MBTag tag_handle, 
                             const MBEntityType type,
                             MBRange &entities);

  //! gets all entity handles that match a type, tag and tag value 
  MBErrorCode get_entities_with_tag_value( const MBEntityType type,
                                            const MBTag tag_handle,
                                            const void* value,
                                            MBRange &entities );
  
  //! gets all entity handles that match a type, tag and tag value 
  MBErrorCode get_entities_with_tag_value( const MBRange &input_range,
                                            const MBEntityType type,
                                            const MBTag tag_handle,
                                            const void* value,
                                            MBRange &entities );
  
  MBErrorCode get_entities_with_tag_values( MBEntityType type,
                                             const MBTag *tags,
                                             const void* const* values,
                                             const int num_tags,
                                             MBRange &entities,
                                             const int condition);
  
  MBErrorCode get_entities_with_tag_values( const MBRange &input_range,
                                             const MBEntityType type,
                                             const MBTag *tags,
                                             const void* const* values,
                                             const int num_tags,
                                             MBRange &entities,
                                             const int condition);
  
  //! gets number of entities that match a type and tag
  MBErrorCode get_number_entities( const MBTag tag_handle, const MBEntityType type,
                             int& num_entities);

  //! get number of entities with tag set
  MBErrorCode get_number_entities( const MBTag tag_handle, unsigned long& num_ents );

  //! gets number of entities that match a type and tag
  MBErrorCode get_number_entities(const MBRange &input_range,
                                   const MBTag tag_handle, const MBEntityType type,
                                    int& num_entities);

  //! finds the first entity handle that matches this data
  MBEntityHandle find_entity( const MBTag tag_handle, const void* data );

  //! gets a tag handle by name and entity handle
  MBTag get_handle(const char *tag_name);

  //! get all the tags which have been defined for this entity
  MBErrorCode get_tags(const MBEntityHandle entity, std::vector<MBTag> &all_tags);
  
  //! get all the tags which have been defined
  MBErrorCode get_tags(std::vector<MBTag> &all_tags);
  
    //! get the default value for a given tag
  MBErrorCode get_default_data(const MBTag tag_handle, void *data);

    //! get the default value for a given tag
  MBErrorCode get_default_data_ref(const MBTag tag_handle, const void *& data);

  //! get information about a tag
  const TagInfo* get_tag_info(const char *tag_name ) const;
  const TagInfo* get_tag_info( MBTag tag_handle ) const;
  
  unsigned long get_memory_use( MBTag tag_handle ) const;
  
  MBErrorCode get_memory_use( MBTag tag_handle,
                              unsigned long& total,
                              unsigned long& per_entity ) const;

private:

  //! table of tag ids and tag information
  //! do we really need to do it this way?
  //! we at least need a table between names and tag ids
  //! and we need information between tag ids and tag usage
  std::map< MBTag, TagInfo >  mTagTable;

  //! container for storing the sparse data and tag ids
  SparseTagSuperCollection* mSparseData;

  //! container for storing the static sparse data and tag ids
  //StaticSparseTagCollection* mStaticSparseData;
  // static tags currently don't fit in OOA

  //! manager for dense data
  DenseTagSuperCollection* mDenseData;

  //! manager for the bit data
  MBBitServer* mBitServer;

};



inline TagInfo::TagInfo(const TagInfo& copy)
  : mTagName( copy.mTagName ),
    mDataSize( copy.mDataSize ),
    dataType( copy.dataType )
{
  if (NULL != copy.mDefaultValue)
  {
    mDefaultValue = new unsigned char[mDataSize];
    memcpy(mDefaultValue, copy.mDefaultValue, mDataSize);
  }
  else 
    mDefaultValue = NULL;
}

inline TagInfo::TagInfo( const char* name, 
                         int size, 
                         MBDataType type,
                         const void* default_value)
 : mTagName( name ),
   mDataSize( size ),
   dataType( type )
{
  if (NULL != default_value) 
  {
    mDefaultValue = new unsigned char[size];
    memcpy(mDefaultValue, default_value, size);
  }
  else 
    mDefaultValue = NULL;
}

inline TagInfo &TagInfo::operator=(const TagInfo &rhs)
{
  mTagName = rhs.mTagName;
  mDataSize = rhs.mDataSize;
  if (NULL != rhs.mDefaultValue) 
  {
    // delete the old data and make a new one (could be different size)
    if(mDefaultValue != NULL)
      delete [] mDefaultValue;
    mDefaultValue = new unsigned char[mDataSize];
    memcpy(mDefaultValue, rhs.mDefaultValue, mDataSize);
  }
  else
  {
    // delete old data
    if(mDefaultValue != NULL)
      delete [] mDefaultValue;
    mDefaultValue = NULL;
  }
  
  return *this;
}

inline TagInfo::~TagInfo() 
{
  // clean up default value
  if (NULL != mDefaultValue)
    delete [] mDefaultValue;
}


inline const TagInfo* TagServer::get_tag_info( const char *tag_name ) const
{
  if(NULL == tag_name || strcmp(tag_name, "") == 0)
    return NULL;

  std::map<MBTag, TagInfo>::const_iterator iterator;
  const std::string temp_name = tag_name;
  for(iterator = mTagTable.begin(); iterator != mTagTable.end(); ++iterator)
  {
      if ( temp_name.size()  == iterator->second.get_name().size() &&
           std::equal(temp_name.begin(), temp_name.end(), iterator->second.get_name().begin() ) )
      {
        return &(iterator->second);
      }
  }
  return NULL;
}

inline const TagInfo* TagServer::get_tag_info( MBTag tag_handle ) const
{
  std::map<MBTag, TagInfo>::const_iterator iterator = mTagTable.find(tag_handle);

  if ( iterator != mTagTable.end() )
  {
    return &(iterator->second);
  }
  
  return NULL;
}


#endif //TAG_SERVER_HPP




