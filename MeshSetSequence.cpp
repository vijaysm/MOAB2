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

/**\file MeshSetSequence.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-04-30
 */

#include "MeshSetSequence.hpp"
#include "SequenceManager.hpp"

MeshSetSequence::MeshSetSequence( MBEntityHandle start,
                                  MBEntityID count,
                                  const unsigned* flags,
                                  SequenceData* data )
  : EntitySequence( start, count, data )
{
  initialize( flags );
}

MeshSetSequence::MeshSetSequence( MBEntityHandle start,
                                  MBEntityID count,
                                  unsigned flags,
                                  SequenceData* data )
  : EntitySequence( start, count, data )
{
  std::vector<unsigned> vect( count, flags );
  initialize( &vect[0] );
}

MeshSetSequence::MeshSetSequence( MBEntityHandle start,
                                  MBEntityID count,
                                  const unsigned* flags,
                                  MBEntityID data_size )
  : EntitySequence( start, count, new SequenceData( 1, start, start+data_size-1) )
{
  initialize( flags );
}

MeshSetSequence::MeshSetSequence( MBEntityHandle start,
                                  MBEntityID count,
                                  unsigned flags,
                                  MBEntityID data_size )
  : EntitySequence( start, count, new SequenceData( 1, start, start+data_size-1) )
{
  std::vector<unsigned> vect( count, flags );
  initialize( &vect[0] );
}

void MeshSetSequence::initialize( const unsigned* flags )
{
  if (!data()->get_sequence_data(0))
    data()->create_sequence_data( 0, SET_SIZE );
 
  MBEntityID offset = start_handle() - data()->start_handle();
  for (MBEntityID i = 0; i < size(); ++i)
    allocate_set( flags[i], i+offset );
}

MeshSetSequence::~MeshSetSequence()
{
  MBEntityID offset = start_handle() - data()->start_handle();
  MBEntityID count = size();
  for (MBEntityID i = 0; i < count; ++i)
    deallocate_set( i + offset );
}

EntitySequence* MeshSetSequence::split( MBEntityHandle here )
{
  return new MeshSetSequence( *this, here );
}

MBErrorCode MeshSetSequence::pop_back( MBEntityID count )
{
  MBEntityID offset = end_handle() + 1 - count - data()->start_handle();
  MBErrorCode rval = EntitySequence::pop_back(count);
  if (MB_SUCCESS == rval)
    for (MBEntityID i = 0; i < count; ++i)
      deallocate_set( i + offset );
  return rval;
}

MBErrorCode MeshSetSequence::pop_front( MBEntityID count )
{
  MBEntityID offset = start_handle() - data()->start_handle();
  MBErrorCode rval = EntitySequence::pop_front(count);
  if (MB_SUCCESS == rval)
    for (MBEntityID i = 0; i < count; ++i)
      deallocate_set( i + offset );
  return rval;
}

MBErrorCode MeshSetSequence::push_back( MBEntityID count, const unsigned* flags )
{
  MBEntityID offset = end_handle() + 1 - data()->start_handle();
  MBErrorCode rval = EntitySequence::append_entities( count );
  if (MB_SUCCESS == rval)
    for (MBEntityID i = 0; i < count; ++i)
      allocate_set( flags[i], i + offset );
  return rval;
}

MBErrorCode MeshSetSequence::push_front( MBEntityID count, const unsigned* flags )
{
  MBEntityID offset = start_handle() - data()->start_handle() - count;
  MBErrorCode rval = EntitySequence::prepend_entities( count );
  if (MB_SUCCESS == rval)
    for (MBEntityID i = 0; i < count; ++i)
      allocate_set( flags[i], i + offset );
  return rval;
}

void MeshSetSequence::get_const_memory_use( unsigned long& per_ent,
                                            unsigned long& seq_size ) const
{
  per_ent = SET_SIZE;
  seq_size = sizeof(*this);
}

unsigned long MeshSetSequence::get_per_entity_memory_use( MBEntityHandle first,
                                                          MBEntityHandle last
                                                        ) const
{
  if (first < start_handle())
    first = start_handle();
  if (last > end_handle())
    last = end_handle();
  
  unsigned long sum = 0;
  for (MBEntityHandle h = first; h <= last; ++h) 
    sum += get_set(h)->get_memory_use();
  return sum;
}

MBErrorCode MeshSetSequence::get_entities( const SequenceManager* seqman,
                                           MBEntityHandle handle,
                                           MBRange& entities,
                                           bool recursive ) const
{
  if (!recursive) {
    get_set(handle)->get_entities( entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<const MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, seqman, &list );
    for (std::vector<const MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_non_set_entities( entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::get_entities( MBEntityHandle handle,
                                std::vector<MBEntityHandle>& entities ) const
{
  get_set(handle)->get_entities( entities );
  return MB_SUCCESS;
}

MBErrorCode MeshSetSequence::get_dimension( const SequenceManager* seqman,
                                            MBEntityHandle handle,
                                            int dimension,
                                            std::vector<MBEntityHandle>& entities,
                                            bool recursive ) const
{
  if (!recursive) {
    get_set(handle)->get_entities_by_dimension( dimension, entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<const MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, seqman, &list );
    for (std::vector<const MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_dimension( dimension, entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::get_dimension( const SequenceManager* seqman,
                                            MBEntityHandle handle,
                                            int dimension,
                                            MBRange& entities,
                                            bool recursive ) const
{
  if (!recursive) {
    get_set(handle)->get_entities_by_dimension( dimension, entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<const MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, seqman, &list );
    for (std::vector<const MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_dimension( dimension, entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::get_type( const SequenceManager* seqman,
                                       MBEntityHandle handle,
                                       MBEntityType type,
                                       std::vector<MBEntityHandle>& entities,
                                       bool recursive ) const
{
  if (!recursive) {
    get_set(handle)->get_entities_by_type( type, entities );
    return MB_SUCCESS;
  }
  else if (type == MBENTITYSET) {
    return recursive_get_sets( handle, seqman, 0, 0, &entities );
  }
  else {
    std::vector<const MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, seqman, &list );
    for (std::vector<const MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_type( type, entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::get_type( const SequenceManager* seqman,
                                       MBEntityHandle handle,
                                       MBEntityType type,
                                       MBRange& entities,
                                       bool recursive ) const
{
  if (!recursive) {
    get_set(handle)->get_entities_by_type( type, entities );
    return MB_SUCCESS;
  }
  else if (type == MBENTITYSET) {
    return recursive_get_sets( handle, seqman, 0, &entities );
  }
  else {
    std::vector<const MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, seqman, &list );
    for (std::vector<const MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_type( type, entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::num_entities( const SequenceManager* seqman,
                                           MBEntityHandle handle,
                                           int& number,
                                           bool recursive ) const
{
  if (!recursive) {
    number = get_set(handle)->num_entities();
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_entities( seqman, handle, range, true );
    number = range.size();
    return result;
  }
}

MBErrorCode MeshSetSequence::num_dimension( const SequenceManager* seqman,
                                            MBEntityHandle handle,
                                            int dimension,
                                            int& number,
                                            bool recursive ) const
{
  if (!recursive) {
    number = get_set(handle)->num_entities_by_dimension(dimension);
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_dimension( seqman, handle, dimension, range, true );
    number = range.size();
    return result;
  }
}
 
MBErrorCode MeshSetSequence::num_type( const SequenceManager* seqman,
                                       MBEntityHandle handle,
                                       MBEntityType type,
                                       int& number,
                                       bool recursive ) const
{
  if (!recursive) {
    number = get_set(handle)->num_entities_by_type(type);
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_type( seqman, handle, type, range, true );
    number = range.size();
    return result;
  }
}

MBErrorCode MeshSetSequence::recursive_get_sets( MBEntityHandle start_set,
                              const SequenceManager* seq_sets,
                              std::vector<const MBMeshSet*>* sets,
                              MBRange* set_handles,
                              std::vector<MBEntityHandle>* set_vector )
{
  std::set<MBEntityHandle> visited;
  std::vector<MBEntityHandle> stack;
  stack.push_back( start_set );
  bool remove_start_set = true;
  while (!stack.empty()) {
    MBEntityHandle handle = stack.back();
    stack.pop_back();
    
    if (!visited.insert(handle).second) {
      if (handle == start_set)
        remove_start_set = false;
      continue;
    }
    
    const EntitySequence* seq;
    MBErrorCode rval = seq_sets->find( handle, seq );
    if (MB_SUCCESS != rval)
      return rval;
    
    const MeshSetSequence* mseq = reinterpret_cast<const MeshSetSequence*>(seq);
    const MBMeshSet *ms_ptr = mseq->get_set( handle );
    if (sets)
      sets->push_back( ms_ptr );
    
    MBRange tmp_range;
    ms_ptr->get_entities_by_type( MBENTITYSET, tmp_range );
    std::copy( tmp_range.begin(), tmp_range.end(), std::back_inserter( stack ) );
  }
  
  if (set_handles) {
    if (remove_start_set)
      visited.erase( start_set );
    MBRange::iterator hint = set_handles->begin();
    std::set<MBEntityHandle>::iterator it;
    for (it = visited.begin(); it != visited.end(); ++it)
      hint = set_handles->insert( hint, *it, *it );
  }
  
  if (set_vector) {
    if (remove_start_set)
      visited.erase( start_set );
    std::copy( visited.begin(), visited.end(), std::back_inserter(*set_vector) );
  }
    
  return MB_SUCCESS;
}

MBErrorCode MeshSetSequence::recursive_get_sets( MBEntityHandle start_set,
                              SequenceManager* seq_sets,
                              std::vector<MBMeshSet*>& sets )
{
  std::set<MBEntityHandle> visited;
  std::vector<MBEntityHandle> stack;
  stack.push_back( start_set );
  while (!stack.empty()) {
    MBEntityHandle handle = stack.back();
    stack.pop_back();
    
    if (!visited.insert(handle).second)
      continue;
    
    EntitySequence* seq;
    MBErrorCode rval = seq_sets->find( handle, seq );
    if (MB_SUCCESS != rval)
      return rval;
    
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    MBMeshSet *ms_ptr = mseq->get_set( handle );
    sets.push_back( ms_ptr );
    
    MBRange tmp_range;
    ms_ptr->get_entities_by_type( MBENTITYSET, tmp_range );
    std::copy( tmp_range.begin(), tmp_range.end(), std::back_inserter( stack ) );
  }
    
  return MB_SUCCESS;
}

MBErrorCode MeshSetSequence::get_parent_child_meshsets( MBEntityHandle meshset,
                                    const SequenceManager* seq_sets,
                                    std::vector<MBEntityHandle>& results,
                                    int num_hops, bool parents ) const
{
  MBErrorCode result = MB_SUCCESS;
  std::vector<MBEntityHandle>::iterator i;
  const MBEntityHandle* array;
  int k, count;
  
    // Skip any meshsets already in input vector (yes, don't
    // get their children either even if num_hops would indicate
    // that we should.)  There is an exception to that if the
    // input meshset is in the list, which is handled by the order
    // of checks in the main loop below.
  std::set<MBEntityHandle> visited;
  for (i = results.begin(); i != results.end(); ++i)
    visited.insert( *i );
    
    // Two lists for breadth-first search
  std::vector<MBEntityHandle> lists[2];
  int index = 0;  // which list to read from (write to lists[1-index])
  lists[index].push_back( meshset ); // begin with input set
    // loop for num_hops (or until no more sets)
  for( ; num_hops && !lists[index].empty(); --num_hops) {
      // for each set at the current num_hops
    for (i = lists[index].begin(); i != lists[index].end(); ++i) {
        // get meshset from handle
      const EntitySequence* seq;
      MBErrorCode rval = seq_sets->find( *i, seq );
      if (MB_SUCCESS != rval)
        return rval;
      const MBMeshSet *ms_ptr = reinterpret_cast<const MeshSetSequence*>(seq)->get_set( *i );
      
        // querying for parents or children?
      array = parents ? ms_ptr->get_parents(count) : ms_ptr->get_children(count);
        // copy any parents/children we haven't visited yet into list
      for (k = 0; k < count; ++k) 
        if (visited.insert(array[k]).second) 
          lists[1-index].push_back(array[k]);
    }
    
      // iterate
    lists[index].clear();
    index = 1 - index;
      // append each level of sets to the output list.
      // note: to make a more useful search (e.g. get entities 3 hops away, 
      // rather than entities up to and including 3 hops) move this outside
      // the loop, but then need to handle the get all (num_hops < 0) case
      // specially.
    std::copy( lists[index].begin(), lists[index].end(), std::back_inserter(results) );
  }
  
  return result;
}

MBErrorCode MeshSetSequence::get_parents( const SequenceManager* seqman,
                                          MBEntityHandle handle,
                                          std::vector<MBEntityHandle>& parents,
                                          int num_hops ) const
{
  if (num_hops == 1) {
    int count;
    const MBEntityHandle* array = get_set( handle )->get_parents(count);  
    if (parents.empty()) {
      parents.resize(count);
      std::copy( array, array + count, parents.begin() );
      return MB_SUCCESS;
    }
    else if (!count) {
      return MB_SUCCESS;
    }
  }
  
  if (num_hops > 0)
    return get_parent_child_meshsets( handle, seqman, parents, num_hops, true );
  else
    return get_parent_child_meshsets( handle, seqman, parents, -1, true );
}

MBErrorCode MeshSetSequence::get_children( const SequenceManager* seqman,
                                           MBEntityHandle handle,
                                           std::vector<MBEntityHandle>& children,
                                           int num_hops ) const
{
  if (num_hops == 1) {
    int count;
    const MBEntityHandle* array = get_set( handle )->get_children(count);  
    if (children.empty()) {
      children.resize(count);
      std::copy( array, array + count, children.begin() );
      return MB_SUCCESS;
    }
    else if (!count) {
      return MB_SUCCESS;
    }
  }

  if (num_hops > 0) 
    return get_parent_child_meshsets( handle, seqman, children, num_hops, false );
  else 
    return get_parent_child_meshsets( handle, seqman, children, -1, false );
}

MBErrorCode MeshSetSequence::num_parents( const SequenceManager* seqman,
                                          MBEntityHandle handle,
                                          int& number,
                                          int num_hops ) const
{
  if (num_hops == 1) {
    number = get_set( handle )->num_parents();
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle> parents;
  MBErrorCode result = get_parents( seqman, handle, parents, num_hops );
  number = parents.size();
  return result;
}


MBErrorCode MeshSetSequence::num_children( const SequenceManager* seqman,
                                           MBEntityHandle handle,
                                           int& number,
                                           int num_hops ) const
{
  if (num_hops == 1) {
    number = get_set( handle )->num_children();
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle> children;
  MBErrorCode result = get_children( seqman, handle, children, num_hops );
  number = children.size();
  return result;
}
