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

#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "PolyEntitySequence.hpp"
#include "EntitySequenceManager.hpp"
#include "MBRange.hpp"

#define CHUNK_SIZE 4096

PolyEntitySequence::PolyEntitySequence(EntitySequenceManager* seq_manager,
                                       MBEntityHandle start_handle, 
                                       MBEntityID num_entities,
                                       int nodes_per_element, bool all_handles_used,
                                       bool allocate_connect)
: ElementEntitySequence(seq_manager, start_handle, num_entities, 0, all_handles_used, false)
{
  if (allocate_connect) {
      // for a poly sequence, the nodes_per_element is really the total # vertices/faces
    polyConn.reserve(nodes_per_element);
    mElements = &polyConn[0];
    mLastIndex.reserve(num_entities);
  }
  
  mNodesPerElement = nodes_per_element;
}

PolyEntitySequence::~PolyEntitySequence() 
{
    // set the mElement pointer to NULL to avoid deleting in ~ElementEntitySequence
  mElements = NULL;
}

bool PolyEntitySequence::is_valid_entity(MBEntityHandle entity) const
{
  return ( (entity-mStartEntityHandle) < (mNumEntities + mDeadEntities.size()) && 
           entity >= mStartEntityHandle &&
           std::find(mDeadEntities.begin(), mDeadEntities.end(), entity) == 
           mDeadEntities.end());
}

MBErrorCode PolyEntitySequence::get_connectivity(MBEntityHandle entity,
                                                 const MBEntityHandle*& conn,
                                                 int &num_vertices,
                                                 const bool /*topological_connectivity*/) const
{
  MBEntityID index = entity - mStartEntityHandle;
  if (!is_valid_entity(entity))
    return MB_ENTITY_NOT_FOUND;
  
  if (0 == index) {
    conn = &polyConn[0];
    num_vertices = mLastIndex[0]+1;
  }
  else {
    conn = &polyConn[mLastIndex[index-1]+1];
    num_vertices = mLastIndex[index] - mLastIndex[index-1];
  }
  
  return MB_SUCCESS;
}

MBErrorCode PolyEntitySequence::get_connectivity(MBEntityHandle entity,
                                                 std::vector<MBEntityHandle> &conn,
                                                 const bool) const
{
  if (!is_valid_entity(entity)) return MB_ENTITY_NOT_FOUND;

  const MBEntityHandle *this_conn;
  int numv;
  MBErrorCode result = get_connectivity(entity, this_conn, numv);
  if (MB_SUCCESS != result) return result;
  
  conn.reserve(numv);
  conn.insert(conn.end(), this_conn, this_conn+numv);
  return MB_SUCCESS;
}

MBErrorCode PolyEntitySequence::get_connectivity_array(MBEntityHandle*& conn_array)
{
  conn_array = mElements;
  return MB_SUCCESS;
}

MBErrorCode PolyEntitySequence::get_index_array(int*& index_array)
{
  index_array = &mLastIndex[0];
  return MB_SUCCESS;
}

MBErrorCode PolyEntitySequence::add_entity(const MBEntityHandle *conn,
                                           const int num_conn, MBEntityHandle &handle)
{
    // make sure input connectivity is the right type
#ifndef NDEBUG
  const int conn_dim = get_type() == MBPOLYGON ? 0 : 2;
  for (const MBEntityHandle *it = conn; it < conn+num_conn; it++) {
    if (MBCN::Dimension(TYPE_FROM_HANDLE(*it)) != conn_dim)
      return MB_FAILURE;
  }
#endif

  if (mNumEntities+mDeadEntities.size() == 0 || 
      mLastIndex[mNumEntities+mDeadEntities.size()-1]+num_conn+1 > mNumAllocated) {
      // reserve more space
    if (mNumEntities+mDeadEntities.size() == 0) {
      polyConn.reserve(CHUNK_SIZE);
      mNumAllocated = CHUNK_SIZE;
    }
    else {
      mNumAllocated *= 2;
      polyConn.reserve(mNumAllocated);
    }
    
      // need to re-assign mElements
    mElements = &polyConn[0];
  }

    // check a few things: num entities agrees with index array size...
  assert(mNumEntities+mDeadEntities.size() == mLastIndex.size());

    // ... last index agrees with connectivity array size...
  assert(0 == mNumEntities+mDeadEntities.size() || 
         mLastIndex[mNumEntities+mDeadEntities.size()-1] == (int) polyConn.size()-1);
  
  handle = mStartEntityHandle + mNumEntities + mDeadEntities.size();

    // add the last index to the index array
  if (0 == mNumEntities+mDeadEntities.size()) 
    mLastIndex.push_back(num_conn-1);
  else
    mLastIndex.push_back(mLastIndex[mNumEntities+mDeadEntities.size()-1]+num_conn);

  mNumEntities++;
  
    // put the connectivity on the end of that array
  polyConn.insert(polyConn.end(), conn, conn+num_conn);

  return MB_SUCCESS;
}

MBErrorCode PolyEntitySequence::set_connectivity(MBEntityHandle entity, 
                                                 const MBEntityHandle *conn,
                                                 const int num_vertices)
{
    // only works if entity is next handle to be allocated
  if (entity != mStartEntityHandle+mNumEntities+mDeadEntities.size())
    return MB_NOT_IMPLEMENTED;

  assert(TYPE_FROM_HANDLE(entity) == TYPE_FROM_HANDLE(mStartEntityHandle));
  
    // ok, it's the right one; just add an entity
  MBEntityHandle dum_handle;
  MBErrorCode result = add_entity(conn, num_vertices, dum_handle);
  assert(entity == dum_handle);
  return result;
}

void PolyEntitySequence::get_entities(MBRange& entities) const
{
  entities.insert(mStartEntityHandle, mStartEntityHandle+mNumEntities+mDeadEntities.size()-1);
  for (std::vector<MBEntityHandle>::const_iterator vit = mDeadEntities.begin();
       vit != mDeadEntities.end(); vit++)
    entities.erase(*vit);
}

void PolyEntitySequence::free_handle(MBEntityHandle entity) 
{
    // add this handle to the dead list, and zero out its connectivity
  MBEntityID index = entity - mStartEntityHandle;
  if (!is_valid_entity(entity)) return;
  
  int start_index;
  if (0 == index) start_index = 0;
  else start_index = mLastIndex[index-1]+1;
  for (int i = start_index; i <= mLastIndex[index]; i++)
    mElements[i] = 0;

    // now add it to the dead list
  mDeadEntities.push_back(entity);

    // decrement number of entities
  mNumEntities--;
}

MBEntityHandle PolyEntitySequence::get_unused_handle() 
{
  return mNumEntities + mDeadEntities.size() + mStartEntityHandle;
}


void PolyEntitySequence::get_memory_use( unsigned long& used, 
                                         unsigned long& allocated) const
{
  allocated = sizeof(*this)
       + mFreeEntities.size() / 8
       + polyConn.capacity() * sizeof(MBEntityHandle)
       + mLastIndex.capacity() * sizeof(int)
       + mDeadEntities.capacity() * sizeof(MBEntityHandle)
       ;
  used = 0;
  for (MBEntityHandle h = get_start_handle(); h <= get_end_handle(); ++h)
    if (is_valid_entity(h))
      used += get_memory_use( h );
}

unsigned long PolyEntitySequence::get_memory_use( MBEntityHandle h ) const
{
  MBEntityID id = h - get_start_handle();
  unsigned long result = 0;
  if (!h)
    result = mLastIndex[h];
  else if (h < mLastIndex.size())
    result = mLastIndex[h] - mLastIndex[h-1];
  return result * sizeof(MBEntityHandle) + sizeof(int);
}
