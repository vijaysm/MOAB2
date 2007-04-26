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

/**\file MeshSetManager.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-08-10
 */

#include "MeshSetManager.hpp"
#include "MBMeshSet.hpp"
#include "MBInternals.hpp"

#include <set>

static inline MBEntityType split_handle( MBEntityHandle handle, 
                                         unsigned& i, 
                                         unsigned& j, 
                                         unsigned& k )
{
  k = handle & (MESHSET_MANAGER_LEVEL_THREE_COUNT-1);
  handle = handle >> MESHSET_MANAGER_LEVEL_THREE_BITS;
  j = handle & (MESHSET_MANAGER_LEVEL_TWO_COUNT-1);
  handle = handle >> MESHSET_MANAGER_LEVEL_TWO_BITS;
  i = handle & (MESHSET_MANAGER_LEVEL_ONE_COUNT-1);
  handle = handle >> MESHSET_MANAGER_LEVEL_ONE_BITS;
  return (MBEntityType)handle;
}

static inline MBEntityHandle make_handle( unsigned i, unsigned j, unsigned k )

{
  MBEntityHandle result = (MBEntityHandle)MBENTITYSET;
  result = (result << MESHSET_MANAGER_LEVEL_ONE_BITS) | i;
  result = (result << MESHSET_MANAGER_LEVEL_TWO_BITS) | j;
  result = (result << MESHSET_MANAGER_LEVEL_THREE_BITS) | k;
  return result;
}

MeshSetManager::MeshSetManager( AEntityFactory* a_ent_fact, 
                                const MBProcConfig& proc_info )
  : aEntityFactory( a_ent_fact ), procInfo( proc_info )
{
  memset( setArrays, 0, sizeof(setArrays) );
}
  
MeshSetManager::~MeshSetManager()
{
  for (unsigned i = 0; i < MESHSET_MANAGER_LEVEL_ONE_COUNT; ++i) {
    if (!setArrays[i])
      continue;
    for (unsigned j = 0; j < MESHSET_MANAGER_LEVEL_TWO_COUNT; ++j) {
      if (!setArrays[i][j])
        continue;
#ifndef HAVE_BOOST_POOL_OBJECT_POOL_HPP
      for (unsigned k = 0; k < MESHSET_MANAGER_LEVEL_THREE_COUNT; ++k) {
        delete setArrays[i][j][k];
      }
#endif
      free( setArrays[i][j] );
    }
    free( setArrays[i] );
  }
  
#ifdef HAVE_BOOST_POOL_OBJECT_POOL_HPP
  rangePool.~object_pool();
  vectorPool.~object_pool();
  new (&rangePool) boost::object_pool<MBMeshSet_MBRange>;
  new (&vectorPool) boost::object_pool<MBMeshSet_Vector>;
#endif
}  

MBEntityHandle MeshSetManager::find_next_free_handle(unsigned proc_id)
{
  if (proc_id >= lastID.size())
    lastID.resize( proc_id + 1, MB_START_ID - 1 );
  
  int junk;
  MBEntityID start_id = lastID[proc_id];
  MBEntityID id = start_id;
  MBEntityHandle handle;

  do {
    ++id;
    if (id > procInfo.max_id())
      id = MB_START_ID;
    if (id == start_id)
      return 0; // exhausted ID space??
    
    handle = CREATE_HANDLE( MBENTITYSET, procInfo.id(id, proc_id), junk );  
  } while (get_mesh_set( handle ));
  
  lastID[proc_id] = id;
  return handle;
}


MBErrorCode MeshSetManager::create_mesh_set( unsigned options, 
                                             MBEntityID start_id,
                                             unsigned int proc_id,
                                             MBEntityHandle& handle )
{
    // get handle
  if (start_id == 0) {
    handle = find_next_free_handle( proc_id );
    if (!handle) // no free handles
      return MB_FAILURE;
  }
  else {
    int junk;
    handle = CREATE_HANDLE( MBENTITYSET, procInfo.id(start_id, proc_id), junk );
    if (get_mesh_set( handle ))
      return MB_ALREADY_ALLOCATED;
  }
  
    // Make sure slot is allocated in array
  unsigned i, j, k;
  split_handle( handle, i, j, k );
  MBMeshSet*** two = setArrays[i];
  if (!two) {
    size_t size = sizeof(MBMeshSet**) * MESHSET_MANAGER_LEVEL_TWO_COUNT;
    two = setArrays[i] = (MBMeshSet***)malloc( size );
    memset( two, 0, size );
  }
  MBMeshSet** three = two[j];
  if (!three) {
    size_t size = sizeof(MBMeshSet*) * MESHSET_MANAGER_LEVEL_THREE_COUNT;
    three = two[j] = (MBMeshSet**)malloc( size );
    memset( three, 0, size );
  }
  MBMeshSet *& ms_ptr = three[k];

    // Create mesh set instance
  bool track = options & MESHSET_TRACK_OWNER ? true : false;
#ifdef HAVE_BOOST_POOL_OBJECT_POOL_HPP
  if (options & MESHSET_ORDERED) 
    ms_ptr = vectorPool.construct( track );
  else 
    ms_ptr = rangePool.construct( track );
#else
  if (options & MESHSET_ORDERED) 
    ms_ptr = new MBMeshSet_Vector( track );
  else 
    ms_ptr = new MBMeshSet_MBRange( track );
#endif
  return MB_SUCCESS;
}

MBMeshSet* MeshSetManager::get_mesh_set( MBEntityHandle handle ) const
{
  unsigned i, j, k;
  if (split_handle(handle,i,j,k) != MBENTITYSET)
    return 0;
  
  MBMeshSet*** two = setArrays[i];
  if (!two)
    return 0;
  MBMeshSet** three = two[j];
  if (!three)
    return 0;
  return three[k];
}

MBErrorCode MeshSetManager::get_options( MBEntityHandle handle, 
                                         unsigned& options ) const
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  if (dynamic_cast<MBMeshSet_MBRange*>(ms_ptr))
    options = MESHSET_SET;
  else if (dynamic_cast<MBMeshSet_Vector*>(ms_ptr))
    options = MESHSET_ORDERED;
  else
    options = 0;
  
  if (ms_ptr->tracking())
    options |= MESHSET_TRACK_OWNER;
  
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::destroy_mesh_set( MBEntityHandle handle )
{
  unsigned i, j, k;
  if (split_handle(handle,i,j,k) != MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
    // Get MeshSet pointer
  MBMeshSet*** two = setArrays[i];
  if (!two)
    return MB_ENTITY_NOT_FOUND;
  MBMeshSet** three = two[j];
  if (!three || !three[k])
    return MB_ENTITY_NOT_FOUND;
  MBMeshSet* ms_ptr = three[k];
  
    // Remove from global array
  three[k] = 0;
    // Check entire block is empty
  MBMeshSet **iter3, **end3 = three + MESHSET_MANAGER_LEVEL_THREE_COUNT;
  for (iter3 = three; iter3 != end3 && *iter3; ++iter3);
  if (iter3 == end3) {
      // free block
    free( three );
    two[j] = 0;
      // check if entire superblock is empty
    MBMeshSet ***iter2, ***end2 = two + MESHSET_MANAGER_LEVEL_TWO_COUNT;
    for (iter2 = two; iter2 != end2 && *iter2; ++iter2);
    if (iter2 == end2) {
        // free superblock
      free( two );
      setArrays[i] = 0;
    }
  }
  
    // update parents
  int count;
  const MBEntityHandle* parents = ms_ptr->get_parents(count);
  for (int i = 0; i < count; ++i) 
    if (MBMeshSet* parent_ptr = get_mesh_set( parents[i] ))
      parent_ptr->remove_child( handle );

    // update children
  const MBEntityHandle* children = ms_ptr->get_children( count );
  for (int i = 0; i < count; ++i) 
    if (MBMeshSet* child_ptr = get_mesh_set( children[i] ))
      child_ptr->remove_parent( handle );

    // remove contents before deleting so adjacency info is removed
  MBErrorCode rval = ms_ptr->clear( handle, aEntityFactory );
#ifdef HAVE_BOOST_POOL_OBJECT_POOL_HPP
  if (dynamic_cast<MBMeshSet_MBRange*>(ms_ptr))
    rangePool.destroy( (MBMeshSet_MBRange*)ms_ptr );
  else
    vectorPool.destroy( (MBMeshSet_Vector*)ms_ptr );
#else
  delete ms_ptr;
#endif
  return rval;
}

MBErrorCode MeshSetManager::clear_mesh_set( MBEntityHandle handle )
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  return ms_ptr->clear( handle, aEntityFactory );
}

void MeshSetManager::get_all_mesh_sets( MBRange& range ) const
{
  MBRange::iterator ins_pos = range.lower_bound( MBENTITYSET );
  
  for (unsigned i = 0; i < MESHSET_MANAGER_LEVEL_ONE_COUNT; ++i) {
    MBMeshSet*** two = setArrays[i];
    if (!two)
      continue;
    
    for (unsigned j = 0; j < MESHSET_MANAGER_LEVEL_TWO_COUNT; ++j) {
      MBMeshSet** three = two[j];
      if (!three) 
        continue;
      
      unsigned k = 0;
      for (;;) {
        for (; k < MESHSET_MANAGER_LEVEL_THREE_COUNT && !three[k]; ++k);
        if (k == MESHSET_MANAGER_LEVEL_THREE_COUNT)
          break;
        unsigned k0 = k;
        for (; k < MESHSET_MANAGER_LEVEL_THREE_COUNT && three[k]; ++k);
        ins_pos = range.insert( ins_pos, make_handle( i, j, k0 ), make_handle( i, j, k-1 ) );
      }
    } // for j
  } // for k
}

MBErrorCode MeshSetManager::get_entities( MBEntityHandle handle,
                                          MBRange& entities,
                                          bool recursive ) const
{
  if (!recursive) {
    MBMeshSet* ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    ms_ptr->get_entities( entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<MBMeshSet*> list;
    MBErrorCode result = recursive_get_sets( handle, list );
    for (std::vector<MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_non_set_entities( entities );
    return result;
  }
}  

MBErrorCode MeshSetManager::get_entities( MBEntityHandle handle,
                                 std::vector<MBEntityHandle>& entities ) const
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  ms_ptr->get_entities( entities );
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::get_dimension( MBEntityHandle handle,
                                           int dimension,
                                           MBRange& entities,
                                           bool recursive ) const
{
  if (!recursive) {
    MBMeshSet* ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    ms_ptr->get_entities_by_dimension( dimension, entities );
    return MB_SUCCESS;
  }
  else {
    std::vector<MBMeshSet*> list;
    MBErrorCode result = recursive_get_sets( handle, list );
    for (std::vector<MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_dimension( dimension, entities );
    return result;
  }
}

MBErrorCode MeshSetManager::get_type( MBEntityHandle handle,
                                      MBEntityType type,
                                      MBRange& entities,
                                      bool recursive ) const
{
  if (!recursive) {
    MBMeshSet* ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    ms_ptr->get_entities_by_type( type, entities );
    return MB_SUCCESS;
  }
  else {
    if (type == MBENTITYSET) // will never return anything
      return MB_TYPE_OUT_OF_RANGE;
  
    std::vector<MBMeshSet*> list;
    MBErrorCode result = recursive_get_sets( handle, list );
    for (std::vector<MBMeshSet*>::iterator i = list.begin(); i != list.end(); ++i)
      (*i)->get_entities_by_type( type, entities );
    return result;
  }
}

MBErrorCode MeshSetManager::num_entities( MBEntityHandle handle,
                                          int& number,
                                          bool recursive ) const
{
  if (!recursive) {
    MBMeshSet *ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    number = ms_ptr->num_entities( );
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_entities( handle, range, true );
    number = range.size();
    return result;
  }
}
 
MBErrorCode MeshSetManager::num_dimension( MBEntityHandle handle,
                                           int dim,
                                           int& number,
                                           bool recursive ) const
{
  if (!recursive) {
    MBMeshSet *ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    number = ms_ptr->num_entities_by_dimension( dim );
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_dimension( handle, dim, range, true );
    number = range.size();
    return result;
  }
}
 
MBErrorCode MeshSetManager::num_type( MBEntityHandle handle,
                                      MBEntityType type,
                                      int& number,
                                      bool recursive ) const
{
  if (!recursive) {
    MBMeshSet *ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    number = ms_ptr->num_entities_by_type( type );
    return MB_SUCCESS;
  }
  else {
    MBRange range;
    MBErrorCode result = get_type( handle, type, range, true );
    number = range.size();
    return result;
  }
}

MBErrorCode MeshSetManager::recursive_get_sets( MBEntityHandle start_set, 
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
    
    MBMeshSet *ms_ptr = get_mesh_set( handle );
    if (!ms_ptr) 
      return MB_ENTITY_NOT_FOUND;
    sets.push_back( ms_ptr );
    
    MBRange tmp_range;
    ms_ptr->get_entities_by_type( MBENTITYSET, tmp_range );
    std::copy( tmp_range.begin(), tmp_range.end(), std::back_inserter( stack ) );
  }
    
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::subtract( MBEntityHandle h1, MBEntityHandle h2 )
{
  MBMeshSet *s1, *s2;
  s1 = get_mesh_set( h1 );
  s2 = get_mesh_set( h2 );
  if (!s1 || !s2)
    return MB_ENTITY_NOT_FOUND;
  
  return s1->subtract( s2, h1, aEntityFactory );
}

MBErrorCode MeshSetManager::intersect( MBEntityHandle h1, MBEntityHandle h2 )
{
  MBMeshSet *s1, *s2;
  s1 = get_mesh_set( h1 );
  s2 = get_mesh_set( h2 );
  if (!s1 || !s2)
    return MB_ENTITY_NOT_FOUND;
  
  return s1->intersect( s2, h1, aEntityFactory );
}

MBErrorCode MeshSetManager::unite( MBEntityHandle h1, MBEntityHandle h2 )
{
  MBMeshSet *s1, *s2;
  s1 = get_mesh_set( h1 );
  s2 = get_mesh_set( h2 );
  if (!s1 || !s2)
    return MB_ENTITY_NOT_FOUND;
  
  return s1->unite( s2, h1, aEntityFactory );
}

MBErrorCode MeshSetManager::add_entities( MBEntityHandle handle,
                                          const MBRange& entities )
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  return ms_ptr->add_entities( entities, handle, aEntityFactory );
}

MBErrorCode MeshSetManager::add_entities( MBEntityHandle handle,
                                          const MBEntityHandle* entities,
                                          int num_entities )
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  return ms_ptr->add_entities( entities, num_entities, handle, aEntityFactory );
}

MBErrorCode MeshSetManager::remove_entities( MBEntityHandle handle,
                                             const MBRange& entities )
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  return ms_ptr->remove_entities( entities, handle, aEntityFactory );
}

MBErrorCode MeshSetManager::remove_entities( MBEntityHandle handle,
                                             const MBEntityHandle* entities,
                                             int num_entities )
{
  MBMeshSet* ms_ptr = get_mesh_set( handle );
  if (!ms_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  return ms_ptr->remove_entities( entities, num_entities, handle, aEntityFactory );
}


MBErrorCode MeshSetManager::get_parent_child_meshsets( MBEntityHandle meshset,
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
      MBMeshSet* ms_ptr = get_mesh_set( *i );
      if (!ms_ptr) {
        result = MB_ENTITY_NOT_FOUND;
        continue;
      }
      
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

MBErrorCode MeshSetManager::get_parents( MBEntityHandle handle,
                            std::vector<MBEntityHandle>& parents,
                            int num_hops ) const
{
  if (num_hops == 1) {
    MBMeshSet *ms_ptr = get_mesh_set( handle );
    if( !ms_ptr )
      return MB_ENTITY_NOT_FOUND;
    
    int count;
    const MBEntityHandle* array = ms_ptr->get_parents(count);  
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

MBErrorCode MeshSetManager::get_children( MBEntityHandle handle,
                            std::vector<MBEntityHandle>& children,
                            int num_hops ) const
{
  if (num_hops == 1) {
    MBMeshSet *ms_ptr = get_mesh_set( handle );
    if( !ms_ptr )
      return MB_ENTITY_NOT_FOUND;
      
    int count;
    const MBEntityHandle* array = ms_ptr->get_children(count);  
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

MBErrorCode MeshSetManager::num_parents( MBEntityHandle handle,
                                         int& number,
                                         int num_hops ) const
{
  if (num_hops == 1) {
    MBMeshSet* ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    
    number = ms_ptr->num_parents();
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle> parents;
  MBErrorCode result = get_parents( handle, parents, num_hops );
  number = parents.size();
  return result;
}


MBErrorCode MeshSetManager::num_children( MBEntityHandle handle,
                                         int& number,
                                         int num_hops ) const
{
  if (num_hops == 1) {
    MBMeshSet* ms_ptr = get_mesh_set( handle );
    if (!ms_ptr)
      return MB_ENTITY_NOT_FOUND;
    
    number = ms_ptr->num_children();
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle> children;
  MBErrorCode result = get_children( handle, children, num_hops );
  number = children.size();
  return result;
}

MBErrorCode MeshSetManager::add_parents( MBEntityHandle to_handle,
                                         const MBEntityHandle* parents,
                                         int count )
{
    // get target meshset
  MBMeshSet* to_set_ptr = get_mesh_set( to_handle );
  if (!to_set_ptr)
    return MB_ENTITY_NOT_FOUND;
  
    // make sure all parent handles are valid
  for (int i = 0; i < count; ++i)
    if (!get_mesh_set( parents[i] ))
      return MB_ENTITY_NOT_FOUND;
  
    // add parent handles to child list on target meshset
  for (int i = 0; i < count; ++i)
    to_set_ptr->add_parent( parents[i] );
    
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::add_children( MBEntityHandle to_handle,
                                          const MBEntityHandle* children,
                                          int count )
{
    // get target meshset
  MBMeshSet* to_set_ptr = get_mesh_set( to_handle );
  if (!to_set_ptr)
    return MB_ENTITY_NOT_FOUND;
  
    // make sure all child handles are valid
  for (int i = 0; i < count; ++i)
    if (!get_mesh_set( children[i] ))
      return MB_ENTITY_NOT_FOUND;
  
    // add child handles to parent list on target meshset
  for (int i = 0; i < count; ++i)
    to_set_ptr->add_child( children[i] );
    
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::add_parent_child( MBEntityHandle parent_handle,
                                              MBEntityHandle child_handle )
{
  MBMeshSet* parent_ptr = get_mesh_set( parent_handle );
  MBMeshSet* child_ptr = get_mesh_set( child_handle );
  if (!parent_ptr || !child_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  parent_ptr->add_child( child_handle );
  child_ptr->add_parent( parent_handle );
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::remove_parent_child( MBEntityHandle parent_handle,
                                                 MBEntityHandle child_handle )
{
  MBMeshSet* parent_ptr = get_mesh_set( parent_handle );
  MBMeshSet* child_ptr = get_mesh_set( child_handle );
  if (!parent_ptr || !child_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  if (!parent_ptr->remove_child( child_handle ))
    return MB_FAILURE;
  if (child_ptr->remove_parent( child_handle ))
    return MB_SUCCESS;
  parent_ptr->add_child( child_handle );
  return MB_FAILURE;
}

MBErrorCode MeshSetManager::remove_parent( MBEntityHandle from_handle,
                                           MBEntityHandle parent_handle )
{
  MBMeshSet* from_ptr = get_mesh_set( from_handle );
  if (!from_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  from_ptr->remove_parent( parent_handle );
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::remove_child( MBEntityHandle from_handle,
                                           MBEntityHandle child_handle )
{
  MBMeshSet* from_ptr = get_mesh_set( from_handle );
  if (!from_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  from_ptr->remove_child( child_handle );
  return MB_SUCCESS;
}

void MeshSetManager::get_memory_use( unsigned long& entity_total,
                                     unsigned long& memory_total) const
{
  memory_total = sizeof(*this);
  memory_total += sizeof(MBEntityID) * lastID.capacity();
  entity_total = 0;

  for (MBEntityID i = 0; i < MESHSET_MANAGER_LEVEL_ONE_COUNT; ++i) {
    if (!setArrays[i])
      continue;
    memory_total += sizeof(MBMeshSet**) * MESHSET_MANAGER_LEVEL_TWO_COUNT;
    for (MBEntityID j = 0; j < MESHSET_MANAGER_LEVEL_TWO_COUNT; ++j) {
      if (!setArrays[i][j])
        continue;
      memory_total += sizeof(MBMeshSet*) * MESHSET_MANAGER_LEVEL_THREE_COUNT;
      for (MBEntityID k = 0; k < MESHSET_MANAGER_LEVEL_THREE_COUNT; ++k) {
        if (setArrays[i][j][k])
          entity_total += setArrays[i][j][k]->get_memory_use();
      }
    }
  }
  memory_total += entity_total;
}

MBErrorCode MeshSetManager::get_memory_use( const MBRange& range,
                                            unsigned long& per_entity,
                                            unsigned long& amortized )
{
  per_entity =0;
  amortized = 0;
  MBRange::iterator r = range.lower_bound( MBENTITYSET );

  unsigned long ent_count = 0, total_count = 0;
  for (MBEntityID i = 0; i < MESHSET_MANAGER_LEVEL_ONE_COUNT; ++i) {
    if (!setArrays[i]) 
      continue;
    amortized += sizeof(MBMeshSet**) * MESHSET_MANAGER_LEVEL_TWO_COUNT;
    for (MBEntityID j = 0; j < MESHSET_MANAGER_LEVEL_TWO_COUNT; ++j) {
      if (!setArrays[i][j]) 
        continue;
      amortized += sizeof(MBMeshSet*) * MESHSET_MANAGER_LEVEL_THREE_COUNT;
      for (MBEntityID k = 0; k < MESHSET_MANAGER_LEVEL_THREE_COUNT; ++k) {
        if (!setArrays[i][j][k])
          continue;
        
        ++total_count;
        MBEntityHandle h = make_handle(i,j,k);
        r = range.lower_bound( r, range.end(), h );
        if (r != range.end() && *r == h) {
          ++ent_count;
          per_entity += setArrays[i][j][k]->get_memory_use();
        }
      }
    }
  }
  
  amortized += sizeof(*this) + sizeof(MBEntityID) * lastID.capacity();
  amortized = (unsigned long)( (double)amortized * ent_count / total_count );
  amortized += per_entity;
  return MB_SUCCESS;
}

