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

#include "MBTypes.h"
#include "MBInternals.hpp"
#include "TagInfo.hpp"

class MBRange;
class SparseTagSuperCollection;
class SequenceManager;
class MBBitServer;

//! SparseTagServer class which associates tag data with entities
class TagServer
{
public:
  //! constructor
  TagServer( SequenceManager* seqman );
  //! destructor
  virtual ~TagServer();

  //! add a tag
  MBErrorCode add_tag(const char *tag_name, 
                       const int data_size,
                       const MBTagType storage,
                       const MBDataType data_type,
                       MBTag &tag_handle,
                       const void *default_data = NULL,
                       int default_value_size = 0);

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

  //! set global/mesh value of tag
  MBErrorCode set_mesh_data( const MBTag tag_handle, const void* data, int size = 0);

  //! set the value of a tag
  MBErrorCode set_data(const MBTag tag_handle, const MBEntityHandle entity_handle, const void* data );
  
  MBErrorCode set_data(const MBTag tag_handle, const MBEntityHandle* entity_handles, const int num_entities, const void* data );
  
  MBErrorCode set_data(const MBTag tag_handle, const MBRange& entity_handles, const void* data );

  //! get global/mesh value of tag
  MBErrorCode get_mesh_data( const MBTag tag_handle, void* data, int& size ) const;

  //! get the value of a tag
  MBErrorCode get_data(const MBTag tag_handle, const MBEntityHandle entity_handle, void* data );
  
  MBErrorCode get_data(const MBTag tag_handle, const MBEntityHandle* entity_handles, const int num_ents, void* data );
  
  MBErrorCode get_data(const MBTag tag_handle, const MBRange& entity_handles, void* data );
  
  //! set the value of a tag
  MBErrorCode set_bits(const MBTag tag_handle, const MBEntityHandle entity_handle, unsigned char data );

  //! get the value of a tag
  MBErrorCode get_bits(const MBTag tag_handle, const MBEntityHandle entity_handle, unsigned char& data );

  //! remove global/mesh value of tag
  MBErrorCode remove_mesh_data( const MBTag tag_handle );

  //! remove the tag data on an entity
  MBErrorCode remove_data( const MBTag tag_handle, const MBEntityHandle entity_handle );

  //! gets all entity handles that match a type and tag
  MBErrorCode get_entities( const MBTag tag_handle, 
                             const MBEntityType type,
                             MBRange &entities);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities( const MBTag tag_handle, 
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
                                            MBRange &entities,
                                            int value_size = 0 );
  
  //! gets all entity handles that match a type, tag and tag value 
  MBErrorCode get_entities_with_tag_value( const MBRange &input_range,
                                            const MBEntityType type,
                                            const MBTag tag_handle,
                                            const void* value,
                                            MBRange &entities,
                                            int value_size = 0 );
  
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

  //! gets a tag handle by name and entity handle
  MBTag get_handle(const char *tag_name) const;

  //! get all the tags which have been defined for this entity
  MBErrorCode get_tags(const MBEntityHandle entity, std::vector<MBTag> &all_tags);
 
  //! get all the tags which have a global/mesh value
  MBErrorCode get_mesh_tags(std::vector<MBTag> &all_tags) const;
 
  //! get all the tags which have been defined
  MBErrorCode get_tags(std::vector<MBTag> &all_tags);
  
    //! get the default value for a given tag
  MBErrorCode get_default_data(const MBTag tag_handle, void *data, int& size);

    //! get the default value for a given tag
  MBErrorCode get_default_data_ref(const MBTag tag_handle, const void *& data, int& size);

  //! get information about a tag
  const TagInfo* get_tag_info(const char *tag_name ) const;
  const TagInfo* get_tag_info( MBTag tag_handle ) const;
  TagInfo* get_tag_info( MBTag tag_handle );
  
  unsigned long get_memory_use( MBTag tag_handle ) const;
  
  MBErrorCode get_memory_use( MBTag tag_handle,
                              unsigned long& total,
                              unsigned long& per_entity ) const;

private:

  //! Table of tag ids and tag information
  //! Primary (array) index is tag type.  
  //! Secondary (std::vector) index is tag id less one (tag ids begin with 1).
  std::vector<TagInfo> mTagTable[MB_TAG_LAST+1];

  //! container for storing the sparse data and tag ids
  SparseTagSuperCollection* mSparseData;
  
  //! SequenceManager required to access dense tag data
  SequenceManager* sequenceManager;

  //! manager for the bit data
  MBBitServer* mBitServer;

};

inline const TagInfo* TagServer::get_tag_info( const char *tag_name ) const
{
  const MBTag handle = get_handle( tag_name );
  return handle ? get_tag_info( handle ) : 0;
}

inline const TagInfo* TagServer::get_tag_info( MBTag tag_handle ) const
{
  const MBTagId id = ID_FROM_TAG_HANDLE( tag_handle );
  const MBTagType type = PROP_FROM_TAG_HANDLE( tag_handle );
  if (id <= mTagTable[type].size() && mTagTable[type][id-1].is_valid())
    return &mTagTable[type][id-1];
  else
    return NULL;
}

inline TagInfo* TagServer::get_tag_info( MBTag tag_handle )
{
  const MBTagId id = ID_FROM_TAG_HANDLE( tag_handle );
  const MBTagType type = PROP_FROM_TAG_HANDLE( tag_handle );
  if (id <= mTagTable[type].size() && mTagTable[type][id-1].is_valid())
    return &mTagTable[type][id-1];
  else
    return NULL;
}

#endif //TAG_SERVER_HPP




