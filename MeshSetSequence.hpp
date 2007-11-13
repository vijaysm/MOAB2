/*
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

/**\file MeshSetSequence.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-04-30
 */

#ifndef MESH_SET_SEQUENCE_HPP
#define MESH_SET_SEQUENCE_HPP

#include "EntitySequence.hpp"
#include "MBMeshSet.hpp"
#include "SequenceData.hpp"

class SequenceManager;

class MeshSetSequence : public EntitySequence
{
public:

  MeshSetSequence( MBEntityHandle start,
                   MBEntityID count,
                   const unsigned* flags,
                   SequenceData* data );
  
  MeshSetSequence( MBEntityHandle start,
                   MBEntityID count,
                   unsigned flags,
                   SequenceData* data );

  MeshSetSequence( MBEntityHandle start,
                   MBEntityID count,
                   const unsigned* flags,
                   MBEntityID sequence_size );
  
  MeshSetSequence( MBEntityHandle start,
                   MBEntityID count,
                   unsigned flags,
                   MBEntityID sequence_size );

  virtual ~MeshSetSequence();

  EntitySequence* split( MBEntityHandle here );
  
  SequenceData* create_data_subset( MBEntityHandle, MBEntityHandle ) const
    { return 0; }
  
  MBErrorCode pop_back( MBEntityID count );
  MBErrorCode pop_front( MBEntityID count );
  MBErrorCode push_back( MBEntityID count, const unsigned* flags );
  MBErrorCode push_front( MBEntityID count, const unsigned* flags );
  
  void get_const_memory_use( unsigned long& bytes_per_entity,
                             unsigned long& size_of_sequence ) const;
  unsigned long get_per_entity_memory_use( MBEntityHandle first,
                                           MBEntityHandle last ) const;


  inline MBMeshSet* get_set( MBEntityHandle h );
  inline const MBMeshSet* get_set( MBEntityHandle h ) const;
  
  MBErrorCode get_entities( MBEntityHandle set, std::vector<MBEntityHandle>& entities ) const;
  MBErrorCode get_entities(  SequenceManager const* seqman, MBEntityHandle set,                    MBRange& entities, bool recursive ) const;
  MBErrorCode get_dimension( SequenceManager const* seqman, MBEntityHandle set, int dim,           MBRange& entities, bool recursive ) const;
  MBErrorCode get_type(      SequenceManager const* seqman, MBEntityHandle set, MBEntityType type, MBRange& entities, bool recursive ) const;
  
  MBErrorCode num_entities(  SequenceManager const* seqman, MBEntityHandle set,                    int& count, bool recursive ) const;
  MBErrorCode num_dimension( SequenceManager const* seqman, MBEntityHandle set, int dim,           int& count, bool recursive ) const;
  MBErrorCode num_type(      SequenceManager const* seqman, MBEntityHandle set, MBEntityType type, int& count, bool recursive ) const;

  MBErrorCode get_parents ( SequenceManager const* seqman, MBEntityHandle of, std::vector<MBEntityHandle>& parents,  int num_hops ) const;
  MBErrorCode get_children( SequenceManager const* seqman, MBEntityHandle of, std::vector<MBEntityHandle>& children, int num_hops ) const;
  MBErrorCode num_parents ( SequenceManager const* seqman, MBEntityHandle of, int& number, int num_hops ) const;
  MBErrorCode num_children( SequenceManager const* seqman, MBEntityHandle of, int& number, int num_hops ) const;
  
private:

  MeshSetSequence( MeshSetSequence& split_from, MBEntityHandle split_at )
    : EntitySequence( split_from, split_at )
    {}

  void initialize( const unsigned* set_flags );
  
  MBErrorCode get_parent_child_meshsets( MBEntityHandle meshset,
                                    SequenceManager const* set_sequences,
                                    std::vector<MBEntityHandle>& results,
                                    int num_hops, bool parents ) const;
                                    
  static MBErrorCode recursive_get_sets( MBEntityHandle start_set,
                            SequenceManager const* set_sequences,
                            std::vector<const MBMeshSet*>& sets_out );
  static MBErrorCode recursive_get_sets( MBEntityHandle start_set,
                            SequenceManager* set_sequences,
                            std::vector<MBMeshSet*>& sets_out );
  
  enum {
    SET_SIZE = (sizeof(MBMeshSet_MBRange) > sizeof(MBMeshSet_Vector)) ?
                sizeof(MBMeshSet_MBRange) : sizeof(MBMeshSet_Vector)
  };

  inline const unsigned char* array() const
    { return reinterpret_cast<const unsigned char*>(data()->get_sequence_data(0)); }

  inline unsigned char* array()
    { return reinterpret_cast<unsigned char*>(data()->get_sequence_data(0)); }
    
  inline void allocate_set( unsigned flags, MBEntityID index )
  {
    const bool tracking = (0 != (flags&MESHSET_TRACK_OWNER));
    unsigned char* const ptr = array() + index * SET_SIZE;
    if (flags & MESHSET_ORDERED)
      new (ptr) MBMeshSet_Vector(tracking);
    else
      new (ptr) MBMeshSet_MBRange(tracking);
  }
    
  inline void deallocate_set( MBEntityID index ) 
  {
    MBMeshSet* set = reinterpret_cast<MBMeshSet*>(array() + SET_SIZE * index );
    if (set->vector_based())
      reinterpret_cast<MBMeshSet_Vector*>(set)->~MBMeshSet_Vector();
    else
      reinterpret_cast<MBMeshSet_MBRange*>(set)->~MBMeshSet_MBRange();
  }
};

inline MBMeshSet* MeshSetSequence::get_set( MBEntityHandle h )
{
  return reinterpret_cast<MBMeshSet*>(array() + SET_SIZE*(h - data()->start_handle()));
}
inline const MBMeshSet* MeshSetSequence::get_set( MBEntityHandle h ) const
{
  return reinterpret_cast<const MBMeshSet*>(array() + SET_SIZE*(h - data()->start_handle()));
}

#endif
