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

#define CHUNK_SIZE 4096

PolyEntitySequence::PolyEntitySequence(EntitySequenceManager* seq_manager,
                                       MBEntityHandle start_handle, int num_entities,
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
/*
  if(all_handles_used)
  {
    mNumEntities = num_entities;
    mFirstFreeIndex = -1;
  }
  else
  {
    seq_manager->notify_not_full(this);
    std::vector<bool>(mNumAllocated, true).swap(mFreeEntities);

      // do the minimum to support this, for now
    mNumEntities = 0;
    mFirstFreeIndex = 0;
  }
*/
}

PolyEntitySequence::~PolyEntitySequence() 
{
    // set the mElement pointer to NULL to avoid deleting in ~ElementEntitySequence
  mElements = NULL;
}

bool PolyEntitySequence::is_valid_entity(MBEntityHandle entity) const
{
  return ((int) (entity-mStartEntityHandle) < mNumEntities && entity >= mStartEntityHandle);
}

MBErrorCode PolyEntitySequence::get_connectivity(MBEntityHandle entity,
                                                 const MBEntityHandle*& conn,
                                                 int &num_vertices,
                                                 const bool /*topological_connectivity*/) const
{
  int index = entity - mStartEntityHandle;
  assert (0 <= index);
  
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
  if (!is_valid_entity(entity)) return MB_FAILURE;

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
#ifdef NDEBUG
  MBEntityType target_type = (TYPE_FROM_HANDLE(mStartEntityHandle) == MBPOLYGON ? 
                              MBVERTEX : MBPOLYGON);
  for (const MBEntityHandle *it = conn; it < conn+num_conn; it++) {
    if (TYPE_FROM_HANDLE(*it) != target_type)
      return MB_FAILURE;
  }
#endif

  if (mNumEntities == 0 || mLastIndex[mNumEntities-1]+num_conn+1 > mNumAllocated) {
      // reserve more space
    if (mNumEntities == 0) {
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
  assert(mNumEntities == (int) mLastIndex.size());

    // ... last index agrees with connectivity array size...
  assert(0 == mNumEntities || mLastIndex[mNumEntities-1] == (int) polyConn.size()-1);
  
  handle = mStartEntityHandle + mNumEntities;

    // add the last index to the index array
  if (0 == mNumEntities) 
    mLastIndex.push_back(num_conn-1);
  else
    mLastIndex.push_back(mLastIndex[mNumEntities-1]+num_conn);

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
  if (entity != mStartEntityHandle+mNumEntities)
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
  entities.insert(mStartEntityHandle, mStartEntityHandle+mNumEntities-1);
}

