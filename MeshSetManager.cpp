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

MeshSetManager::MeshSetManager( AEntityFactory* a_ent_fact )
  : procSet( MB_PROC_COUNT ),
    aEntityFactory( a_ent_fact )
  {}
  
MeshSetManager::~MeshSetManager()
{
  std::vector<ProcData>::iterator p;
  std::vector<MBMeshSet*>::iterator i;
  for (p = procSet.begin(); p != procSet.end(); ++p) 
    for (i = p->list.begin(); i != p->list.end(); ++i) 
      delete *i;
}  

MBErrorCode MeshSetManager::create_mesh_set( unsigned options, 
                                             MBEntityHandle start_id,
                                             unsigned int proc_id,
                                             MBEntityHandle& handle )
{
    // Get data for specified processor
  if (proc_id >= procSet.size())
    return MB_INDEX_OUT_OF_RANGE;
  ProcData& data = procSet[proc_id];
  
    // Don't put first index into data.free because an
    // ID of zero isn't valid.  Allocate it w/out putting
    // it in the free list so it is never used.
  if (data.list.empty()) 
    data.list.resize(MB_START_ID,0);
  
    // Get id for mesh set
  if (start_id) {
    if (start_id <= MB_START_ID)
      return MB_INDEX_OUT_OF_RANGE;
    if (start_id < data.list.size() && data.list[start_id])
      return MB_ALREADY_ALLOCATED;
  }
  else if (!data.free.empty()) {
    start_id = data.free.front(); data.free.pop_front();
  }
  else {
    start_id = data.list.size();
  }
  
    // Create handle
  int err;
  handle = CREATE_HANDLE( MBENTITYSET, start_id, proc_id, err );
  if (err)
    return MB_INDEX_OUT_OF_RANGE;
  
    // Create hole in the ID space necessary to accomodate
    // the requested start_id.
  while (start_id > data.list.size()) {
      // normally, deleted sets go on the back because we want to delay
      // reusing them as long as possible.  These were never used, so
      // put on the front.
    data.free.push_front( data.list.size() );
      // set pointer for unused spots to NULL
    data.list.push_back( 0 );
  }
  
    // Create mesh set instance
  MBMeshSet* ms_ptr = 0;
  bool track = options & MESHSET_TRACK_OWNER ? true : false;
  if (options & MESHSET_ORDERED) 
    ms_ptr = new MBMeshSet_Vector( track );
  else 
    ms_ptr = new MBMeshSet_MBRange( track );
  
    // Put a pointer to the MBMeshSet in the list of sets
  if (start_id == data.list.size()) 
    data.list.push_back( ms_ptr );
  else
    data.list[start_id] = ms_ptr;

  return MB_SUCCESS;
}

MBMeshSet* MeshSetManager::get_mesh_set( MBEntityHandle handle ) const
{
  if (TYPE_FROM_HANDLE(handle) != MBENTITYSET)
    return 0;
  unsigned proc = PROC_FROM_HANDLE(handle);
  MBEntityHandle id = ID_FROM_HANDLE(handle);
  if (proc >= procSet.size() || id >= procSet[proc].list.size())
    return 0;
  return procSet[proc].list[id];
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
  if (TYPE_FROM_HANDLE(handle) != MBENTITYSET)
    return MB_TYPE_OUT_OF_RANGE;
  
    // get location
  unsigned proc = PROC_FROM_HANDLE(handle);
  MBEntityHandle id = ID_FROM_HANDLE(handle);
  if (proc >= procSet.size() || 
      id >= procSet[proc].list.size() ||
      !procSet[proc].list[id])
    return MB_ENTITY_NOT_FOUND;
  
    // remove from global list
  MBMeshSet* ms_ptr = procSet[proc].list[id];
  procSet[proc].list[id] = 0;
  procSet[proc].free.push_back( id );
  
    // update parents
  const std::vector<MBEntityHandle>& parents = ms_ptr->get_parents();
  for (size_t i = 0; i < parents.size(); ++i) 
    if ((ms_ptr = get_mesh_set( parents[i] )))
      ms_ptr->remove_child( handle );

    // update children
  const std::vector<MBEntityHandle>& children = ms_ptr->get_children();
  for (size_t i = 0; i < children.size(); ++i) 
    if ((ms_ptr = get_mesh_set( children[i] )))
      ms_ptr->remove_parent( handle );

    // remove contents before deleting so adjacency info is removed
  MBErrorCode rval = ms_ptr->clear( handle, aEntityFactory );
  delete ms_ptr;
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
  int junk;
  MBRange::iterator ins_pos = range.lower_bound( MBENTITYSET );
  
  for (unsigned p = 0; p < procSet.size(); ++p) {
    const ProcData& data = procSet[p];
    const MBEntityHandle count = data.list.size();
    MBEntityHandle i = MB_START_ID;
    for (;;) {
        // skip unoccupied slots
      for (; i < count && !data.list[i]; ++i); //<- Note ending ';'
        // check loop termination condition
      if (i >= count)
        break;
        // store first used entity handle
      const MBEntityHandle beg = CREATE_HANDLE( MBENTITYSET, i, p, junk );
        // loop over all occupied slots
      for (++i; i < count && data.list[i]; ++i); //<- Note ending ';'
        // insert into range
      const MBEntityHandle end = CREATE_HANDLE( MBENTITYSET, i-1, p, junk );
      ins_pos = range.insert( ins_pos, beg, end );
        // we know the current one is NULL, so don't bother checking
        // it at the start of the next iteration.
      ++i;
    }
  }
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
  std::vector<MBEntityHandle>::const_iterator k;
  
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
      
        // querying for parents or childre?
      const std::vector<MBEntityHandle>& setlist = parents ?
                                     ms_ptr->get_parents() : 
                                     ms_ptr->get_children();
        // copy any parents/children we haven't visited yet into list
      for (k = setlist.begin(); k != setlist.end(); ++k) 
        if (visited.insert(*k).second) 
          lists[1-index].push_back(*k);
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
      
    const std::vector<MBEntityHandle>& sp = ms_ptr->get_parents();
    if (parents.empty()) {
      parents = sp;
      return MB_SUCCESS;
    }
    else if (sp.empty()) {
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
      
    const std::vector<MBEntityHandle>& sp = ms_ptr->get_children();
    if (children.empty()) {
      children = sp;
      return MB_SUCCESS;
    }
    else if (sp.empty()) {
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

MBErrorCode MeshSetManager::add_parent( MBEntityHandle to_handle,
                                        MBEntityHandle parent_handle )
{
  MBMeshSet* to_set_ptr = get_mesh_set( to_handle );
  MBMeshSet* parent_ptr = get_mesh_set( parent_handle );
  if (!to_set_ptr || !parent_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  to_set_ptr->add_parent( parent_handle );
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::add_child( MBEntityHandle to_handle,
                                        MBEntityHandle child_handle )
{
  MBMeshSet* to_set_ptr = get_mesh_set( to_handle );
  MBMeshSet* child_ptr = get_mesh_set( child_handle );
  if (!to_set_ptr || !child_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  to_set_ptr->add_child( child_handle );
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
  
  if (!from_ptr->remove_parent( parent_handle ))
    return MB_FAILURE;
  
  return MB_SUCCESS;
}

MBErrorCode MeshSetManager::remove_child( MBEntityHandle from_handle,
                                           MBEntityHandle child_handle )
{
  MBMeshSet* from_ptr = get_mesh_set( from_handle );
  if (!from_ptr)
    return MB_ENTITY_NOT_FOUND;
  
  if (!from_ptr->remove_child( child_handle ))
    return MB_FAILURE;
  
  return MB_SUCCESS;
}

  

