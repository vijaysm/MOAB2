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
#include "EntitySequenceManager.hpp"

MeshSetSequence::MeshSetSequence( EntitySequenceManager* man,
                                  MBEntityHandle start_handle,
                                  MBEntityID num_entities,
                                  unsigned flags )
  : MBEntitySequence( man, start_handle, num_entities )
{
  std::vector<unsigned> flag_vect(num_entities, flags);
  initialize( man, start_handle, num_entities, &flag_vect[0] );
}

MeshSetSequence::MeshSetSequence( EntitySequenceManager* man,
                                  MBEntityHandle start_handle,
                                  MBEntityID num_entities,
                                  const unsigned* flags )
  : MBEntitySequence( man, start_handle, num_entities )
{
  initialize( man, start_handle, num_entities, flags );
}

void MeshSetSequence::initialize( EntitySequenceManager* man,
                                  MBEntityHandle start_handle,
                                  MBEntityID num_entities,
                                  const unsigned* flags )
{
  man->entity_sequence_created( this );
  
    // allocate storage
  mSets = new unsigned char[SET_SIZE * num_entities];
  mFreeEntities.clear();
  if (flags) {
    mNumEntities = num_entities;
    mFirstFreeIndex = -1;
    for (MBEntityID i = 0; i < num_entities; ++i) 
      allocate_set( flags[i], i );
  }
  else {
    man->notify_not_full( this );
    mNumEntities = 0;
    mFirstFreeIndex = 0;
    mFreeEntities.resize( mNumAllocated, true );
    for (MBEntityID i = 0; i < num_entities; ++i)
      next_free(i) = i + 1;
    next_free(num_entities-1) = -1; 
  }
}

MeshSetSequence::~MeshSetSequence()
{
  mSequenceManager->entity_sequence_deleted( this );
  if (mSets) {
    for (MBEntityID i = 0; i < mNumAllocated; ++i) 
      if (is_valid_entity(get_start_handle()+i)) 
        deallocate_set( i );
    delete [] mSets;
  }
}

MBEntityHandle MeshSetSequence::get_unused_handle()
{
  return 0;
}

MBEntityHandle MeshSetSequence::add_meshset( unsigned flags )
{
  if (mFirstFreeIndex == -1)
    return 0;
  const MBEntityID index = mFirstFreeIndex;
  
  mFreeEntities[index] = false;
  mFirstFreeIndex = next_free(index);
  if (mLastDeletedIndex == index) 
    mLastDeletedIndex = -1;
  
  allocate_set( flags, index );
  mNumEntities++;
  if (mNumEntities == mNumAllocated) {
    mSequenceManager->notify_full(this);
    std::vector<bool> empty;
    mFreeEntities.swap( empty);
  }
  
  return get_start_handle() + index;
}

void MeshSetSequence::free_handle( MBEntityHandle handle )
{
  if (!is_valid_entity(handle))
    return;
  const MBEntityID index = handle - get_start_handle();
  
    // free any memory allocated by the MBMeshSet
  deallocate_set( index );

    // decerement count of valid entities
  if(mNumEntities == mNumAllocated) {
    mSequenceManager->notify_not_full(this);
    mFreeEntities.resize( mNumAllocated, false );
  }
  --mNumEntities;

    // mark this entity as invalid
  mFreeEntities[index] = true;
  
    // Add this entity to the free list.
    // Free list is maintained in sorted order.

    // insert at beginning?
  if (mFirstFreeIndex == -1 || mFirstFreeIndex > index) {
    next_free(index) = mFirstFreeIndex;
    mFirstFreeIndex = mLastDeletedIndex = index;
    return;
  }
  
    // Find entry to insert after.
    
  MBEntityID prev_index = mFirstFreeIndex;
    // mLastDeletedIndex is a cache of the last deleted
    // entity.  Used to speed up sequential deletion.
  if (mLastDeletedIndex != -1 && mLastDeletedIndex < index)
     prev_index = mLastDeletedIndex;
    // Search list for location to insert at
  MBEntityID next = next_free(prev_index);
  while (next != -1 && next < index) {
    prev_index = next;
    next = next_free(next);
  }
    // insert in free list
  mLastDeletedIndex = index;
  next_free(index) = next_free(prev_index);
  next_free(prev_index) = index;
}  

void MeshSetSequence::get_entities( MBRange& entities ) const
{
  MBRange::iterator iter = entities.insert(get_start_handle(), get_end_handle());
  for(MBEntityID index = mFirstFreeIndex; index != -1; index = next_free(index))
  {
    iter += get_start_handle() + index - *iter;
    iter = entities.erase(iter);
  }
}

MBEntityID MeshSetSequence::get_next_free_index( MBEntityID prev_free_index ) const
{
  return prev_free_index < 0 ? mFirstFreeIndex : next_free(prev_free_index);
}

void MeshSetSequence::get_memory_use( unsigned long& used,
                                      unsigned long& allocated ) const
{
  used = 0;
  allocated = sizeof(*this) + mFreeEntities.capacity()/8;
  allocated += mNumAllocated * SET_SIZE;
  for (MBEntityHandle h = get_start_handle(); h <= get_end_handle(); ++h) {
    if (is_valid_entity(h)) {
      unsigned long m = get_set(h)->get_memory_use();
      used += SET_SIZE + m;
      allocated += m;
    }
  }
}

unsigned long MeshSetSequence::get_memory_use( MBEntityHandle h ) const
{
  return get_set(h)->get_memory_use();
}

MBErrorCode MeshSetSequence::get_entities( MBEntityHandle handle,
                                           MBRange& entities,
                                           bool recursive ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  if (!recursive) {
    get_set(handle)->get_entities( entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, list );
    for (std::vector<MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_non_set_entities( entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::get_entities( MBEntityHandle handle,
                                std::vector<MBEntityHandle>& entities ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  get_set(handle)->get_entities( entities );
  return MB_SUCCESS;
}

MBErrorCode MeshSetSequence::get_dimension( MBEntityHandle handle,
                                            int dimension,
                                            MBRange& entities,
                                            bool recursive ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  if (!recursive) {
    get_set(handle)->get_entities_by_dimension( dimension, entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, list );
    for (std::vector<MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_dimension( dimension, entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::get_type(      MBEntityHandle handle,
                                            MBEntityType type,
                                            MBRange& entities,
                                            bool recursive ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  if (!recursive) {
    get_set(handle)->get_entities_by_type( type, entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<MBMeshSet*> list;
    MBErrorCode rval = recursive_get_sets( handle, list );
    for (std::vector<MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_type( type, entities );
    return rval;
  }
}

MBErrorCode MeshSetSequence::num_entities( MBEntityHandle handle,
                                           int& number,
                                           bool recursive ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  if (!recursive) {
    number = get_set(handle)->num_entities();
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_entities( handle, range, true );
    number = range.size();
    return result;
  }
}

MBErrorCode MeshSetSequence::num_dimension( MBEntityHandle handle,
                                            int dimension,
                                            int& number,
                                            bool recursive ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  if (!recursive) {
    number = get_set(handle)->num_entities_by_dimension(dimension);
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_dimension( handle, dimension, range, true );
    number = range.size();
    return result;
  }
}
 
MBErrorCode MeshSetSequence::num_type( MBEntityHandle handle,
                                       MBEntityType type,
                                       int& number,
                                       bool recursive ) const
{
  if (!is_valid_entity(handle))
    return MB_ENTITY_NOT_FOUND;
  
  if (!recursive) {
    number = get_set(handle)->num_entities_by_type(type);
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_type( handle, type, range, true );
    number = range.size();
    return result;
  }
}

MBErrorCode MeshSetSequence::recursive_get_sets( MBEntityHandle start_set,
                              std::vector<MBMeshSet*>& sets ) const
{
  std::set<MBEntityHandle> visited;
  std::vector<MBEntityHandle> stack;
  stack.push_back( start_set );
  while (!stack.empty()) {
    MBEntityHandle handle = stack.back();
    stack.pop_back();
    
    if (!visited.insert(handle).second)
      continue;
    
    MBEntitySequence* seq;
    MBErrorCode rval = mSequenceManager->find( handle, seq );
    if (MB_SUCCESS != rval)
      return rval;
    
    MeshSetSequence* mseq = reinterpret_cast<MeshSetSequence*>(seq);
    if (!mseq->is_valid_entity(handle))
      return MB_ENTITY_NOT_FOUND;
    
    
    MBMeshSet *ms_ptr = mseq->get_set( handle );
    sets.push_back( ms_ptr );
    
    MBRange tmp_range;
    ms_ptr->get_entities_by_type( MBENTITYSET, tmp_range );
    std::copy( tmp_range.begin(), tmp_range.end(), std::back_inserter( stack ) );
  }
    
  return MB_SUCCESS;
}

MBErrorCode MeshSetSequence::get_parent_child_meshsets( MBEntityHandle meshset,
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
      MBEntitySequence* seq;
      MBErrorCode rval = mSequenceManager->find( *i, seq );
      if (MB_SUCCESS != rval)
        return rval;
      if (!seq->is_valid_entity(*i))
        return MB_ENTITY_NOT_FOUND;
      MBMeshSet *ms_ptr = reinterpret_cast<MeshSetSequence*>(seq)->get_set( *i );
      
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

MBErrorCode MeshSetSequence::get_parents( MBEntityHandle handle,
                            std::vector<MBEntityHandle>& parents,
                            int num_hops ) const
{
  if (num_hops == 1) {
   if (!is_valid_entity(handle))
      return MB_ENTITY_NOT_FOUND;
    
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
    return get_parent_child_meshsets( handle, parents, num_hops, true );
  else
    return get_parent_child_meshsets( handle, parents, -1, true );
}

MBErrorCode MeshSetSequence::get_children( MBEntityHandle handle,
                            std::vector<MBEntityHandle>& children,
                            int num_hops ) const
{
  if (num_hops == 1) {
    if (!is_valid_entity(handle))
      return MB_ENTITY_NOT_FOUND;
    
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
    return get_parent_child_meshsets( handle, children, num_hops, false );
  else 
    return get_parent_child_meshsets( handle, children, -1, false );
}

MBErrorCode MeshSetSequence::num_parents( MBEntityHandle handle,
                                         int& number,
                                         int num_hops ) const
{
  if (num_hops == 1) {
    if (!is_valid_entity(handle))
      return MB_ENTITY_NOT_FOUND;
    
    number = get_set( handle )->num_parents();
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle> parents;
  MBErrorCode result = get_parents( handle, parents, num_hops );
  number = parents.size();
  return result;
}


MBErrorCode MeshSetSequence::num_children( MBEntityHandle handle,
                                         int& number,
                                         int num_hops ) const
{
  if (num_hops == 1) {
    if (!is_valid_entity(handle))
      return MB_ENTITY_NOT_FOUND;
    
    number = get_set( handle )->num_children();
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle> children;
  MBErrorCode result = get_children( handle, children, num_hops );
  number = children.size();
  return result;
}
