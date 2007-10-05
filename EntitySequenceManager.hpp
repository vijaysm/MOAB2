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

#include "MBForward.hpp"
#include "MBHandleUtils.hpp"
#include <map>

class MBEntitySequence;
class HomCoord;

//! class for managing entity sequences
class EntitySequenceManager
{
public:

  //! constructor
  EntitySequenceManager( const MBHandleUtils &handle_utils);

  //! destructor
  ~EntitySequenceManager();

  MBErrorCode create_scd_sequence(const int imin, const int jmin, const int kmin,
                                   const int imax, const int jmax, const int kmax,
                                   const MBEntityType type,
                                   const MBEntityID hint_start_id,
                                   MBEntityHandle &start,
                                   MBEntitySequence *&seq);
  
  MBErrorCode create_scd_sequence(const HomCoord &coord_min,
                                   const HomCoord &coord_max,
                                   const MBEntityType type,
                                   const MBEntityID hint_start_id,
                                   MBEntityHandle &start,
                                   MBEntitySequence *&seq);
  
  //! creates an entity sequence that will fit number of entities.  uses hint_start
  //! for the start handle if possible
  //! returns the start handle and a pointer to the entity sequence.
  MBErrorCode create_entity_sequence(MBEntityType type, MBEntityID num_ent, int num_nodes_per, 
                                     MBEntityID hint_start_id, int hint_start_proc,
                                     MBEntityHandle& start, MBEntitySequence *&);

  //! Create entity sequence for mesh sets.  If flags is NULL, 
  //! sequence will be created with all handles unused.  If
  //! flags is not NULL, it should be an array of num_ent values.
  MBErrorCode create_meshset_sequence( MBEntityID num_ent,
                                       MBEntityID hint_start_id,
                                       int hint_start_proc,
                                       const unsigned* flags,
                                       MBEntitySequence *&seq );
                                       
  //! finds the specific MBEntitySequence in which MBEntityHandle resides
  MBErrorCode find( MBEntityHandle     entity_handle,
                     MBEntitySequence*& sequence ) const;

  //! get entities from the entity sequences according to type 
  MBErrorCode get_entities(MBEntityType,
                            MBRange &entities) const;
  
  //! get entities from the entity sequences according to type 
  MBErrorCode get_number_entities(MBEntityType, MBEntityID& entities) const;

  //! deletes an entity from the database
  MBErrorCode delete_entity( MBEntityHandle entity );

  //! creates a vertex in the database
  MBErrorCode create_vertex( const unsigned processor_id,
                             const double coords[3], 
                             MBEntityHandle& vertex );

  //! creates an element in the database
  MBErrorCode create_element( MBEntityType type, 
                              const unsigned processor_id,
                              const MBEntityHandle *conn_array,
                              const unsigned num_vertices,
                              MBEntityHandle& element);

  MBErrorCode create_mesh_set( unsigned proc_id,
                               unsigned flags,
                               MBEntityHandle& h );
                               
  MBErrorCode allocate_mesh_set( MBEntityHandle handle, unsigned flags );

  //! return a const reference to the map of sequences
  const std::map<MBEntityHandle, MBEntitySequence*>*
    entity_map( MBEntityType type ) const { return &(mSequenceMap[type]); }

  void entity_sequence_created(MBEntitySequence* seq);
  void entity_sequence_deleted(MBEntitySequence* seq);
  void notify_full(MBEntitySequence* seq);
  void notify_not_full(MBEntitySequence* seq);
  
  void get_memory_use( unsigned long& total_entity_storage,
                       unsigned long& total_storage ) const;
  void get_memory_use( MBEntityType type,
                       unsigned long& total_entity_storage,
                       unsigned long& total_storage ) const;
  MBErrorCode get_memory_use( const MBRange& entities,
                              unsigned long& total_entity_storage,
                              unsigned long& total_amortized_storage ) const;

private:
  
  //! creates an entity sequence with a start handle and number of entities
  MBErrorCode private_create_entity_sequence(MBEntityHandle start,
                                      MBEntityID num_ent, int num_nodes,
                                      bool full, MBEntitySequence *&);

  void delete_all();

  //! last instance of entity sequence that was accessed
  //! each entity type keeps its own cache
  mutable MBEntitySequence* mLastAccessed[MBMAXTYPE];
  
  //! map of all the EntitySequences that are created 
  //! each entity type has its own map, and each processor
  //! ID has it's own map.
  typedef std::map< MBEntityHandle, MBEntitySequence* > SeqMap;
  SeqMap mSequenceMap[MBMAXTYPE];
  SeqMap mPartlyFullSequenceMap[MBMAXTYPE];

    //! get a valid start handle for this type and a hinted-at start id
  MBErrorCode get_start_handle( MBEntityID hint_start, 
                                int hint_proc, 
                                MBEntityType type, 
                                MBEntityID num_ent,
                                MBEntityHandle& start_handle );
                                
    //! helper function for get_start_handle, if the requested entity ID
    //! is not available, find any available range of handles
  MBErrorCode find_free_handles( int proc, 
                                 MBEntityType type, 
                                 MBEntityID num_ent,
                                 MBEntityHandle& start_handle );
  
  const MBHandleUtils handleUtils;
};

#endif


