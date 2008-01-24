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
 * Filename   :     SparseTagSuperCollection.hpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/



#ifndef SPARSE_TAG_SUPER_COLLECTION_HPP
#define SPARSE_TAG_SUPER_COLLECTION_HPP

#ifndef IS_BUILDING_MB
#error "SparseTagSuperCollection.hpp isn't supposed to be included into an application"
#endif

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include <map>
#include <vector>

#include "MBTypes.h"
#include "MBInternals.hpp"
#include "MBRange.hpp"
#include "TagInfo.hpp"

#define get_collection( A ) ((A) < mDataTags.size() ? mDataTags[(A)] : 0)

class SparseTagCollection;

//! a collection of SparseTagCollections
class SparseTagSuperCollection
{
public:
  //! constructor
  SparseTagSuperCollection(){}
  
  //! destructor
  virtual ~SparseTagSuperCollection();

  void reset_data();

  //! allocate new tag id
  MBErrorCode reserve_tag_id(int data_size, MBTagId tag_id);

  //! releases an MBTagId
  MBErrorCode release_tag_id(MBTagId tag_id);

  //! size of data tag
  int tag_size(const MBTagId tag_id) const;

  //! set the data of a tag
  MBErrorCode set_data(const MBTagId tag_handle, const MBEntityHandle entity_handle, const void* data);

  //! get the data of a tag
  MBErrorCode get_data(const MBTagId tag_handle, const MBEntityHandle entity_handle, void* data);

  //! removes data
  MBErrorCode remove_data(const MBTagId tag_handle, const MBEntityHandle entity_handle);

  //! gets all entity handles that match a type and tag
  MBErrorCode get_entities(const MBTagId tag_handle, MBRange &entities);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities(const MBTagId tag_handle, const MBEntityType type,
                            MBRange &entities);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities(const MBRange &range,
                            const MBTagId tag_handle, const MBEntityType type,
                            MBRange &entities);

  //! gets all tags on a given entity handle
  MBErrorCode get_tags(const MBEntityHandle entity,
                       std::vector<MBTag> &all_tags);
  
  //! gets all entity handles that match a tag
  MBErrorCode get_entities_with_tag_value( const MBTagId tag_handle, 
                                           const TagInfo& tag_info,
                                           const MBEntityType type,
                                           MBRange &entities,
                                           const void* tag_value,
                                           int value_size);

  //! gets all entity handles that match a tag
  MBErrorCode get_entities_with_tag_value( const MBRange &range,
                                           const MBTagId tag_handle, 
                                           const TagInfo& tag_info,
                                           const MBEntityType type,
                                           MBRange &entities,
                                           const void* tag_value,
                                           int value_size);

  //! gets the number of entities that match a tag
  MBErrorCode get_number_entities(const MBTagId tag_handle, const MBEntityType type, int& num_ent);


  //! gets the number of entities that match a tag
  MBErrorCode get_number_entities(const MBRange &range,
                                   const MBTagId tag_handle, const MBEntityType type, int& num_ent);

  MBErrorCode get_memory_use( MBTagId tag_id, 
                              unsigned long& total,
                              unsigned long& per_entity );

private:

  std::vector<SparseTagCollection*> mDataTags;
};


#endif //SPARSE_TAG_SUPER_COLLECTION_HPP




