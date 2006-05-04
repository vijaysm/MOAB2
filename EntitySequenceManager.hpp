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

/*!
 *  \class   MBEntitySequenceManager
 *  \authors Karl Merkley & Corey Ernst
 *  \date    3/27/02
 *  \brief   The MBEntitySequence represents a contiguous range of 
 *           mesh entities of a single type.  MBEntitySequence manages
 *           the internal ids of those entities using a starting id and 
 *           the number of entities in the range.  All MBEntitySequences 
 *           for a given MBEntity type are stored in an stl Map.  When a 
 *           MBHandle references an entity in mdb, mdb finds the Entity-
 *           Sequence associated with the MBHandle from the stl Map.  
 *           All enquires for any data associated with the MBEntity are
 *           directed through the MBEntitySequence object to which the
 *           MBEntity's MBHandle refers.  A MBEntitySequence knows 
 *           about all DenseTags defined for that MBEntitySequence, and
 *           contains one TagRange for each of these DenseTags.  The 
 *           MBEntitySequence can query for and return any of the DenseTag
 *           data from the TagRanges it manages. 
 *          
 */ 

#ifndef ENTITY_SEQUENCE_MANAGER_HPP
#define ENTITY_SEQUENCE_MANAGER_HPP

#ifndef IS_BUILDING_MB
#error "EntitySequenceManager.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "EntitySequence.hpp"
#include "MBInternals.hpp"
#include "MBRange.hpp"
#include "HomXform.hpp"
#include <map>

//! class for managing entity sequences
class EntitySequenceManager
{
public:

  //! constructor
  EntitySequenceManager();

  //! destructor
  ~EntitySequenceManager();

  MBErrorCode create_scd_sequence(const int imin, const int jmin, const int kmin,
                                   const int imax, const int jmax, const int kmax,
                                   const MBEntityType type,
                                   const int hint_start_id,
                                   MBEntityHandle &start,
                                   MBEntitySequence *&seq);
  
  MBErrorCode create_scd_sequence(const HomCoord &coord_min,
                                   const HomCoord &coord_max,
                                   const MBEntityType type,
                                   const int hint_start_id,
                                   MBEntityHandle &start,
                                   MBEntitySequence *&seq);
  
  //! creates an entity sequence that will fit number of entities.  uses hint_start
  //! for the start handle if possible
  //! returns the start handle and a pointer to the entity sequence.
  MBErrorCode create_entity_sequence(MBEntityType type, int num_ent, int num_nodes_per, 
                                     int hint_start_id, int hint_start_proc,
                                     MBEntityHandle& start, MBEntitySequence *&);

  //! finds the specific MBEntitySequence in which MBEntityHandle resides
  MBErrorCode find( MBEntityHandle     entity_handle,
                     MBEntitySequence*& sequence ) const;

  //! get entities from the entity sequences according to type 
  MBErrorCode get_entities(MBEntityType,
                            MBRange &entities) const;
  
  //! get entities from the entity sequences according to type 
  MBErrorCode get_number_entities(MBEntityType, int& entities);

  //! deletes an entity from the database
  MBErrorCode delete_entity( MBEntityHandle entity );

  //! creates a vertex in the database
  MBErrorCode create_vertex(const double coords[3], MBEntityHandle& vertex);

  //! creates an element in the database
  MBErrorCode create_element(MBEntityType type, 
                              const MBEntityHandle *conn_array,
                              const int num_vertices,
                              MBEntityHandle& element);

  //! return a const reference to the map of sequences
  const std::map<MBEntityHandle, MBEntitySequence*>*
    entity_map( MBEntityType type ) const { return &(mSequenceMap[type]); }

  void entity_sequence_created(MBEntitySequence* seq);
  void entity_sequence_deleted(MBEntitySequence* seq);
  void notify_full(MBEntitySequence* seq);
  void notify_not_full(MBEntitySequence* seq);

private:
  
  //! creates an entity sequence with a start handle and number of entities
  MBErrorCode private_create_entity_sequence(MBEntityHandle start,
                                      int num_ent, int num_nodes,
                                      bool full, MBEntitySequence *&);

  void delete_all();

  //! last instance of entity sequence that was accessed
  //! each entity type keeps its own cache
  mutable MBEntitySequence* mLastAccessed[MBMAXTYPE];
  
  //! map of all the EntitySequences that are created 
  //! each entity type has its own map
  std::map< MBEntityHandle, MBEntitySequence* >  mSequenceMap[MBMAXTYPE];

  std::map< MBEntityHandle, MBEntitySequence* >  mPartlyFullSequenceMap[MBMAXTYPE];

    //! get a valid start handle for this type and a hinted-at start id
  MBEntityHandle get_start_handle(int hint_start, int hint_proc, MBEntityType type, int num_ent);
  
};

  //! create a structured sequence of vertices or elements
inline MBErrorCode EntitySequenceManager::create_scd_sequence(const HomCoord &coord_min,
                                                               const HomCoord &coord_max,
                                                               const MBEntityType type,
                                                               const int hint_start_id,
                                                               MBEntityHandle &start_handle,
                                                               MBEntitySequence *&seq) 
{
  return create_scd_sequence(coord_min.i(), coord_min.j(), coord_min.k(), 
                             coord_max.i(), coord_max.j(), coord_max.k(), 
                             type, hint_start_id,
                             start_handle, seq);
}

#endif


