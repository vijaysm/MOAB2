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

class MeshSetSequence : public MBEntitySequence
{
public:

  MeshSetSequence( EntitySequenceManager* seq_man,
                   MBEntityHandle start_handle,
                   MBEntityID num_entities,
                   unsigned set_flags );

  MeshSetSequence( EntitySequenceManager* seq_man,
                   MBEntityHandle start_handle,
                   MBEntityID num_entities,
                   const unsigned* set_flags = 0 );

  virtual ~MeshSetSequence();
  virtual MBEntityHandle get_unused_handle();
  MBEntityHandle add_meshset( unsigned flags );
  virtual void free_handle( MBEntityHandle handle );
  virtual void get_entities( MBRange& entities ) const;
  virtual MBEntityID get_next_free_index( MBEntityID prev_free_index ) const;
  virtual void get_memory_use( unsigned long& ,unsigned long& ) const;
  virtual unsigned long get_memory_use( MBEntityHandle ) const;
  
  MBMeshSet* get_set( MBEntityHandle h );
  const MBMeshSet* get_set( MBEntityHandle h ) const;
  
  MBErrorCode get_entities( MBEntityHandle set, MBRange& entities, bool recursive ) const;
  MBErrorCode get_entities( MBEntityHandle set, std::vector<MBEntityHandle>& entities ) const;
  MBErrorCode get_dimension( MBEntityHandle set, int dim, MBRange& entities, bool recursive ) const;
  MBErrorCode get_type( MBEntityHandle set, MBEntityType type, MBRange& entities, bool recursive ) const;
  
  MBErrorCode num_entities( MBEntityHandle set, int& count, bool recursive ) const;
  MBErrorCode num_dimension( MBEntityHandle set, int dim, int& count, bool recursive ) const;
  MBErrorCode num_type( MBEntityHandle set, MBEntityType type, int& count, bool recursive ) const;

  MBErrorCode get_parents ( MBEntityHandle of, std::vector<MBEntityHandle>& parents, int num_hops ) const;
  MBErrorCode get_children( MBEntityHandle of, std::vector<MBEntityHandle>& children, int num_hops ) const;
  MBErrorCode num_parents ( MBEntityHandle of, int& number, int num_hops ) const;
  MBErrorCode num_children( MBEntityHandle of, int& number, int num_hops ) const;
  
private:

  void initialize( EntitySequenceManager* seq_man,
                   MBEntityHandle start_handle,
                   MBEntityID num_entities,
                   const unsigned* set_flags );
  
  MBErrorCode get_parent_child_meshsets( MBEntityHandle meshset,
                                    std::vector<MBEntityHandle>& results,
                                    int num_hops, bool parents ) const;
                                    
  MBErrorCode recursive_get_sets( MBEntityHandle start_set,
                            std::vector<MBMeshSet*>& sets_out ) const ;
  
  enum {
    SET_SIZE = (sizeof(MBMeshSet_MBRange) > sizeof(MBMeshSet_Vector)) ?
                sizeof(MBMeshSet_MBRange) : sizeof(MBMeshSet_Vector)
  };
  
  unsigned char* mSets;
  
  inline MBEntityID& next_free( MBEntityID index )
    { return *reinterpret_cast<MBEntityID*>(mSets + SET_SIZE * index ); }
  inline MBEntityID next_free( MBEntityID index ) const
    { return *reinterpret_cast<MBEntityID*>(mSets + SET_SIZE * index ); }
    
  inline void allocate_set( unsigned flags, MBEntityID index )
  {
    const bool tracking = (0 != (flags&MESHSET_TRACK_OWNER));
    unsigned char* const ptr = mSets + index * SET_SIZE;
    if (flags & MESHSET_ORDERED)
      new (ptr) MBMeshSet_Vector(tracking);
    else
      new (ptr) MBMeshSet_MBRange(tracking);
  }
    
  inline void deallocate_set( MBEntityID index ) 
  {
    MBMeshSet* set = reinterpret_cast<MBMeshSet*>(mSets + SET_SIZE * index );
    if (set->vector_based())
      reinterpret_cast<MBMeshSet_Vector*>(set)->~MBMeshSet_Vector();
    else
      reinterpret_cast<MBMeshSet_MBRange*>(set)->~MBMeshSet_MBRange();
  }
};

inline MBMeshSet* MeshSetSequence::get_set( MBEntityHandle h )
{
  assert(is_valid_entity(h));
  return reinterpret_cast<MBMeshSet*>(mSets + SET_SIZE*(h - get_start_handle()));
}
inline const MBMeshSet* MeshSetSequence::get_set( MBEntityHandle h ) const
{
  assert(is_valid_entity(h));
  return reinterpret_cast<MBMeshSet*>(mSets + SET_SIZE*(h - get_start_handle()));
}

#endif
