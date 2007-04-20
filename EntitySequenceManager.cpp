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
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "PolyEntitySequence.hpp"
#include "EntitySequence.hpp"
#include "EntitySequenceManager.hpp"
#include "MBRange.hpp"
#include "ScdElementSeq.hpp"
#include "ScdVertexSeq.hpp"
#include <assert.h>
#include <algorithm>


EntitySequenceManager::EntitySequenceManager( const MBProcConfig& proc_info )
  : procInfo( proc_info )
{
  memset(mLastAccessed, 0, MBMAXTYPE*sizeof(void*));
}

EntitySequenceManager::~EntitySequenceManager()
{
  // delete the entity sequences
  delete_all();
}

void EntitySequenceManager::entity_sequence_created(MBEntitySequence* seq)
{
  mSequenceMap[seq->get_type()].insert(
    std::pair<MBEntityHandle, MBEntitySequence*> (seq->get_start_handle(), seq));
}

void EntitySequenceManager::entity_sequence_deleted(MBEntitySequence* seq)
{
  mSequenceMap[seq->get_type()].erase(seq->get_start_handle());
  mPartlyFullSequenceMap[seq->get_type()].erase(seq->get_start_handle());  
}

void EntitySequenceManager::notify_full(MBEntitySequence* seq)
{
  mPartlyFullSequenceMap[seq->get_type()].erase(seq->get_start_handle());
}

void EntitySequenceManager::notify_not_full(MBEntitySequence* seq)
{
  mPartlyFullSequenceMap[seq->get_type()].insert(
    std::pair<MBEntityHandle, MBEntitySequence*>(seq->get_start_handle(), seq));
}

  //! create a structured sequence of vertices or elements
MBErrorCode EntitySequenceManager::create_scd_sequence(const int imin, const int jmin, const int kmin,
                                                        const int imax, const int jmax, const int kmax,
                                                        const MBEntityType type,
                                                        const MBEntityID hint_start_id,
                                                        MBEntityHandle &start_handle,
                                                        MBEntitySequence *&seq) 
{
  int this_dim = MBCN::Dimension(type);

    // use > instead of != in the following assert to also catch cases where imin > imax, etc.
  assert((this_dim < 3 || kmax > kmin) &&
         (this_dim < 2 || jmax > jmin) &&
         (this_dim < 1 || imax > imin));

    // compute # entities; not as easy as it would appear...
  MBEntityID num_ent;
  if (MBVERTEX == type)
    num_ent = (MBEntityID)(imax-imin+1)*(MBEntityID)(jmax-jmin+1)*(MBEntityID)(kmax-kmin+1);
  else {
    num_ent = (imax-imin) *
      (this_dim >= 2 ? (jmax-jmin) : 1) *
      (this_dim >= 3 ? (kmax-kmin) : 1);
  }
  
    // get a start handle
  MBErrorCode rval = get_start_handle(hint_start_id, procInfo.rank(), type, num_ent, start_handle);
  if (MB_SUCCESS != rval) return rval;
  
  if (MBVERTEX == type)
      // create a new structured vertex sequence
    seq = new ScdVertexSeq(this, start_handle, 
                           imin, jmin, kmin, imax, jmax, kmax);
  else
    seq = new ScdElementSeq(this, start_handle, 
                           imin, jmin, kmin, imax, jmax, kmax);
    
  if (NULL == seq) return MB_FAILURE;

  else return MB_SUCCESS;
}


  //! create a structured sequence of vertices or elements
MBErrorCode EntitySequenceManager::create_scd_sequence(const HomCoord &coord_min,
                                                       const HomCoord &coord_max,
                                                       const MBEntityType type,
                                                       const MBEntityID hint_start_id,
                                                       MBEntityHandle &start_handle,
                                                       MBEntitySequence *&seq) 
{
  return create_scd_sequence(coord_min.i(), coord_min.j(), coord_min.k(), 
                             coord_max.i(), coord_max.j(), coord_max.k(), 
                             type, hint_start_id,
                             start_handle, seq);
}

/*!
  creates an entity sequence based on number of entities and type.
  uses the hint_start as the start id for the entity handles if possible
  returns the actual start handle and the entity sequence pointer
*/
MBErrorCode EntitySequenceManager::create_entity_sequence( MBEntityType type, MBEntityID num_ent, int num_nodes,
                                                           MBEntityID hint_start, int hint_start_proc, 
                                                           MBEntityHandle& start_handle, 
                                                            MBEntitySequence*& seq)
{
  MBErrorCode rval = get_start_handle(hint_start, hint_start_proc, type, num_ent, start_handle);
  if (MB_SUCCESS != rval)
    return rval;
  
  // actually create the sequence
  return private_create_entity_sequence( start_handle, num_ent, num_nodes, true, seq);
}

MBErrorCode EntitySequenceManager::get_start_handle( MBEntityID hint_start, 
                                                     int proc, 
                                                     MBEntityType type,
                                                     MBEntityID num_ent,
                                                     MBEntityHandle& start_hint_handle) 
{
  if (hint_start < MB_START_ID)
    return find_free_handles( proc, type, num_ent, start_hint_handle );

  // Create handles from input parameters
  int dum = 0;
  start_hint_handle = CREATE_HANDLE(type, procInfo.id(hint_start, proc), dum);
  MBEntityHandle end_hint_handle = start_hint_handle + num_ent - 1;
  MBEntityHandle last_handle = CREATE_HANDLE( type, procInfo.last_id(proc), dum );
  
  // Check if the handle type can accomodate the requested number of handles
  if (end_hint_handle > last_handle)
    return find_free_handles( proc, type, num_ent, start_hint_handle );

  // Find the first entity sequence with a handle greater than requested handle.
  SeqMap::iterator iter = mSequenceMap[type].upper_bound(start_hint_handle);
  
  // Check that the requested handle range is before the beginning
  // of the next sequence.
  if (iter != mSequenceMap[type].end() && iter->first <= end_hint_handle)
    return find_free_handles( proc, type, num_ent, start_hint_handle );
  
  // If there is no previous sequence, then done.
  if (iter == mSequenceMap[type].begin())
    return MB_SUCCESS;
  // Otherwise make sure the start of the requested range is
  // after the previous sequence.
  --iter;
  if (iter->second->get_end_handle() < start_hint_handle)
    return MB_SUCCESS;
  else
    return find_free_handles( proc, type, num_ent, start_hint_handle );
}

MBErrorCode EntitySequenceManager::find_free_handles( int proc,
                                                      MBEntityType type,
                                                      MBEntityID num_ent,
                                                      MBEntityHandle& handle_out )
{
  int dum = 0;

    // get first and largest possible handle for specified proc and type
  MBEntityHandle last_handle = CREATE_HANDLE( type, procInfo.last_id(proc), dum );
  handle_out = CREATE_HANDLE( type, procInfo.first_id(proc), dum );

    // check that handle space is large enough to accomodate requested
    // number of entities
  if (last_handle - handle_out + 1 < (MBEntityHandle)num_ent)
    return MB_MEMORY_ALLOCATION_FAILED;

    // get the first entity sequence for the NEXT rank (processor id)
  SeqMap::iterator iter = mSequenceMap[type].upper_bound( last_handle );
  
    // If map is empty or all sequences have larger processor ID
  if (iter == mSequenceMap[type].begin()) 
    return MB_SUCCESS;
  
    // decrement to get first sequence with same (or smaller) 
    // rank (processor ID)
  --iter;
  
    // If previous sequence is for previous processor, then there
    // are currently no handles allocated for this type and proc.
    // Return the first handle.
  if (procInfo.rank(iter->second->get_end_handle()) < (unsigned)proc)
    return MB_SUCCESS;

    // Otherwise try the handle after those currently allocated
    // for this type and processor.
  handle_out = iter->second->get_end_handle() + 1;
    // Check if enough IDs at end of range
  if (last_handle - handle_out + 1 >= (MBEntityHandle)num_ent)
    return MB_SUCCESS;

    // Not enough available handles at the end of the range of handles
    // for the specified entity type and processor.  Scan the entire
    // range of handles looking for a large enough hole before we give
    // up entirely.
  for (;;) {
    MBEntityHandle last_start = iter->second->get_start_handle();
      // If this is the first sequence for this processor
    if (iter == mSequenceMap[type].begin() ||
        procInfo.rank((--iter)->second->get_end_handle()) < (unsigned)proc) {
      MBEntityHandle first_handle = CREATE_HANDLE( type, procInfo.first_id(proc), dum );
      if (first_handle - last_start > (MBEntityHandle)num_ent) {
        handle_out = last_start - num_ent;
        return MB_SUCCESS;
      }
      break;
    }

      // check if we fit in the current range
    if (last_start - iter->second->get_end_handle() - 1 > (MBEntityHandle)num_ent) {
      handle_out = iter->second->get_end_handle() + 1;
      return MB_SUCCESS;
    }
  }
  return MB_MEMORY_ALLOCATION_FAILED;
}

MBErrorCode EntitySequenceManager::private_create_entity_sequence(MBEntityHandle start,
                                                           MBEntityID num_ent, int num_nodes,
                                                           bool full,
                                                           MBEntitySequence *& seq)
{
  MBEntityType type = TYPE_FROM_HANDLE(start);
  if(type == MBVERTEX)
    seq = new VertexEntitySequence(this, start, num_ent, full);
  else if(type == MBPOLYGON || type == MBPOLYHEDRON)
    seq = new PolyEntitySequence(this, start, num_ent, num_nodes, full);
  else
    seq = new ElementEntitySequence(this, start, num_ent, num_nodes, full);
  
  return MB_SUCCESS;
}


MBErrorCode EntitySequenceManager::create_vertex( const unsigned processor_id,
                                                  const double coords[3], 
                                                  MBEntityHandle& handle )
{
  VertexEntitySequence* seq = 0;

  // see if there is an existing sequence that can take this new vertex
  SeqMap& seq_map = mPartlyFullSequenceMap[MBVERTEX];
  for (SeqMap::iterator i = seq_map.begin(); i != seq_map.end(); ++i) {
    if (procInfo.rank(i->second->get_start_handle()) == processor_id) {
      seq = static_cast<VertexEntitySequence*>(i->second);
      break;
    }
  }

  if (!seq) {
    MBErrorCode rval = get_start_handle( MB_START_ID, processor_id, MBVERTEX, 4096, handle );
    if (MB_SUCCESS != rval)
      return rval;
    seq = new VertexEntitySequence( this, handle, 4096, false );
  }

  handle = seq->get_unused_handle();
  seq->set_coordinates(handle, coords[0], coords[1], coords[2]);
  
  return MB_SUCCESS;
}


MBErrorCode EntitySequenceManager::create_element( MBEntityType type, 
                                                   const unsigned processor_id,
                                                   const MBEntityHandle *conn, 
                                                   const unsigned num_vertices,
                                                   MBEntityHandle& handle )
{
  MBErrorCode rval;
  const bool poly = (MBPOLYGON == type || MBPOLYHEDRON == type);
  const unsigned connlen = poly ? 0 : num_vertices;
  
  ElementEntitySequence* seq = 0;
  SeqMap& seq_map = mPartlyFullSequenceMap[type];
  for (SeqMap::iterator i = seq_map.begin(); i != seq_map.end(); ++i) {
    ElementEntitySequence* tseq = reinterpret_cast<ElementEntitySequence*>(i->second);
    if (tseq->nodes_per_element() == connlen && 
        procInfo.rank( tseq->get_start_handle() ) == processor_id) {
      seq = tseq;
      break;
    }
  }
  
  if (poly) {
    PolyEntitySequence* pseq = reinterpret_cast<PolyEntitySequence*>(seq);
    if (!pseq) {
      rval = get_start_handle( MB_START_ID, processor_id, type, 4096, handle );
      if (MB_SUCCESS != rval)
        return rval;
      pseq = new PolyEntitySequence(this, handle, 4096, 0, false);
    }
    handle = pseq->get_unused_handle();
    return pseq->add_entity(conn, num_vertices, handle);
  }
  
  if (!seq) {
    rval = get_start_handle( MB_START_ID, processor_id, type, 4096, handle );
    if (MB_SUCCESS != rval)
      return rval;
    seq = new ElementEntitySequence(this, handle, 4096, num_vertices, false);
  } 

  handle = seq->get_unused_handle();
  return seq->set_connectivity(handle, conn, num_vertices);
}


MBErrorCode EntitySequenceManager::get_entities(MBEntityType type, MBRange &entities) const
{
  
  //index into the static sequence map to get the sequences according to type only 
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator beg_seq, end_seq;
  beg_seq = mSequenceMap[type].begin();
  end_seq = mSequenceMap[type].end();

  //for each sequence, get all the entity handles it contains

  for (; beg_seq != end_seq; ++beg_seq)
  {
    const MBEntitySequence* tmp_seq = beg_seq->second;
    tmp_seq->get_entities(entities);
  }

  return MB_SUCCESS;
}

MBErrorCode EntitySequenceManager::get_number_entities(MBEntityType type, MBEntityID& num_entities) const 
{
  num_entities = 0;
  //index into the static sequence map to get the sequences according to type only 
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator beg_seq, end_seq;
  beg_seq = mSequenceMap[type].begin();
  end_seq = mSequenceMap[type].end();

  //for each sequence, get all the entity handles it contains

  for (; beg_seq != end_seq; ++beg_seq)
  {
    num_entities += beg_seq->second->number_entities();
  }
  return MB_SUCCESS;
}


void EntitySequenceManager::delete_all()
{
  for(MBEntityType i = MBVERTEX; i<MBMAXTYPE; i++)
  {
    while(!mSequenceMap[i].empty())
      delete mSequenceMap[i].begin()->second;
    mSequenceMap[i].clear();
  } 
}


MBErrorCode EntitySequenceManager::find( MBEntityHandle entity_handle,
                                      MBEntitySequence*& sequence ) const
{
  MBEntityType ent_type = TYPE_FROM_HANDLE(entity_handle);

  // check to see if the sequence is cached
  sequence = mLastAccessed[ent_type];
  if (sequence && sequence->get_start_handle() <= entity_handle &&
     sequence->get_end_handle() >= entity_handle)
  {
    return MB_SUCCESS;
  }
  
  sequence = 0;

  if ( !mSequenceMap[ent_type].empty() )
  {
    if (ID_FROM_HANDLE(entity_handle) < MB_START_ID) 
      return MB_INDEX_OUT_OF_RANGE;

    // create an iterator to look for which sequence entity handle will be found in
    // using lower bounds function of map
    std::map<MBEntityHandle, MBEntitySequence*>::const_iterator iter =
      mSequenceMap[ent_type].upper_bound(entity_handle);
    if (iter == mSequenceMap[ent_type].begin())
      return MB_ENTITY_NOT_FOUND;
    --iter;

    // ensure that the entity is bounded by this sequence.  
    // upper_bound will indicate that this sequence can hold this entity
    // but doesn't guarantee that it does!
    MBEntitySequence* seq = iter->second; 
    if ( (entity_handle >= seq->get_start_handle()) &&
         (entity_handle <= seq->get_end_handle()))
    {
      sequence = seq;
      mLastAccessed[ent_type] = seq;
      return MB_SUCCESS; 
    }
    else
      return MB_ENTITY_NOT_FOUND;
  }

  return MB_ENTITY_NOT_FOUND;
}



MBErrorCode EntitySequenceManager::delete_entity( MBEntityHandle entity )
{
  MBEntitySequence* seq;
  find(entity, seq);
  if(seq != NULL)
  {
    seq->free_handle(entity);
    // leave the sequences around for a while if it is empty
    
    return MB_SUCCESS;
  }
  return MB_FAILURE;
}

void EntitySequenceManager::get_memory_use( unsigned long& entity,
                                            unsigned long& total) const
{
  total = entity = 0;
  unsigned long used, allocated;
  for (unsigned type = 0; type < MBMAXTYPE; ++type)
    for (SeqMap::const_iterator i = mSequenceMap[type].begin(); i != mSequenceMap[type].end(); ++i) {
      i->second->get_memory_use( used, allocated );
      entity += used;
      total += allocated;
    }
}

MBErrorCode EntitySequenceManager::get_memory_use( 
                                     const MBRange& entities,
                                     unsigned long& min_per_entity,
                                     unsigned long& amortized ) const
{
  min_per_entity = 0;
  amortized =0;
  if (entities.empty())
    return MB_SUCCESS;
  
  MBRange::iterator e, i = entities.begin();
  for (unsigned type = TYPE_FROM_HANDLE( *i ); type < MBENTITYSET; ++type) {
    if (i == entities.end())
      break;
    
    const SeqMap& map = mSequenceMap[type];
    if (map.empty())
      continue;
      
    for (SeqMap::const_iterator j = map.begin(); j != map.end(); ++j)
    {
      MBEntitySequence* seq = j->second;
      if (seq->get_end_handle() < *i)
        continue;
      
      if (*i < seq->get_start_handle())
        return MB_ENTITY_NOT_FOUND;
      
      e = entities.upper_bound( i, entities.end(), seq->get_end_handle() );
      MBEntityID count = 0;
      for (; i != e; ++i) {
        min_per_entity += seq->get_memory_use( *i );
        ++count;
      }
      
      unsigned long used, allocated;
      seq->get_memory_use( used, allocated );
      amortized += (unsigned long)( (double)count * allocated / seq->number_entities() );
    }
  }
        
  return MB_SUCCESS;
}


#ifdef TEST

#include <iostream>

#define check( A ) \
  do {  \
    if ((A) == MB_SUCCESS) break; \
    std::cerr << "Error at " << __FILE__ << ":" << __LINE__ << std::endl; \
    return 1; \
  } while(false)
  
#define ensure( A ) \
  do { \
    if ((A)) break; \
    std::cerr << "Failed test at " << __FILE__ << ":" << __LINE__ << std::endl; \
    std::cerr << "  " #A << std::endl; \
    return 1; \
  } while(false)
    

int main()
{
  EntitySequenceManager manager( MBProcConfig(0,1) );
  MBEntitySequence* seq;
  
  // create some sequences
  const unsigned NV[] = { 100, 5, 1000 };
  MBEntityHandle vh[3];
  MBEntitySequence* vs[3];
  MBEntityHandle th;
  MBEntitySequence *ts;
  const MBEntityID TRI_START_ID = 3;
  const unsigned NT = 640;
  check(manager.create_entity_sequence( MBVERTEX, NV[0], 0, 35, 0, vh[0], vs[0] ));
  check(manager.create_entity_sequence( MBVERTEX, NV[1], 0, 8, 0, vh[1], vs[1] ));
  check(manager.create_entity_sequence( MBVERTEX, NV[2], 0, 1, 0, vh[2], vs[2] ));
  check(manager.create_entity_sequence( MBTRI, NT, 3, TRI_START_ID, 0, th, ts ));
  
  // check that returned entity sequences are valid
  ensure( ID_FROM_HANDLE(vh[0]) > 0 && TYPE_FROM_HANDLE(vh[0]) == MBVERTEX );
  ensure( ID_FROM_HANDLE(vh[0]) > 0 && TYPE_FROM_HANDLE(vh[1]) == MBVERTEX );
  ensure( ID_FROM_HANDLE(vh[0]) > 0 && TYPE_FROM_HANDLE(vh[2]) == MBVERTEX );
  ensure( ID_FROM_HANDLE( th  ) > 0 && TYPE_FROM_HANDLE( th  ) == MBTRI    );
  ensure( ID_FROM_HANDLE(th) == TRI_START_ID );
  ensure( vs[0] );
  ensure( vs[1] );  
  ensure( vs[2] );  
  ensure( ts );  
  ensure( vs[0]->get_type() == MBVERTEX );
  ensure( vs[1]->get_type() == MBVERTEX );
  ensure( vs[2]->get_type() == MBVERTEX );
  ensure( ts   ->get_type() == MBTRI    );
  ensure( vs[0]->get_start_handle() == vh[0] );
  ensure( vs[1]->get_start_handle() == vh[1] );
  ensure( vs[2]->get_start_handle() == vh[2] );
  ensure( ts   ->get_start_handle() == th    );
  ensure( vs[0]->get_end_handle() - vs[0]->get_start_handle() == NV[0] - 1 );
  ensure( vs[1]->get_end_handle() - vs[1]->get_start_handle() == NV[1] - 1 );
  ensure( vs[2]->get_end_handle() - vs[2]->get_start_handle() == NV[2] - 1 );
  ensure( ts   ->get_end_handle() - ts   ->get_start_handle() == NT    - 1 );
  
  // construct ranges
  MBRange vertices;
  vertices.insert( vs[0]->get_start_handle(), vs[0]->get_end_handle() );
  vertices.insert( vs[1]->get_start_handle(), vs[1]->get_end_handle() );
  vertices.insert( vs[2]->get_start_handle(), vs[2]->get_end_handle() );
  MBRange tris( ts->get_start_handle(), ts->get_end_handle() );
  
  // check find sequence given first
  check( manager.find( vs[0]->get_start_handle(), seq ) );
  ensure( seq == vs[0] );
  check( manager.find( vs[1]->get_start_handle(), seq ) );
  ensure( seq == vs[1] );
  check( manager.find( vs[2]->get_start_handle(), seq ) );
  ensure( seq == vs[2] );
  check( manager.find( ts   ->get_start_handle(), seq ) );
  ensure( seq == ts    );
  
  // check find sequence given last
  check( manager.find( vs[0]->get_end_handle(), seq ) );
  ensure( seq == vs[0] );
  check( manager.find( vs[1]->get_end_handle(), seq ) );
  ensure( seq == vs[1] );
  check( manager.find( vs[2]->get_end_handle(), seq ) );
  ensure( seq == vs[2] );
  check( manager.find( ts   ->get_end_handle(), seq ) );
  ensure( seq == ts    );
  
  // check find of invalid handle
  MBEntityHandle badhandle = (MBEntityHandle)(ts->get_end_handle() + 1);
  ensure( manager.find( badhandle, seq ) == MB_ENTITY_NOT_FOUND );
  
  // check find of invalid type
  int junk;
  badhandle = CREATE_HANDLE( MBTET, 1, junk );
  ensure( manager.find( badhandle, seq ) == MB_ENTITY_NOT_FOUND );
  
  // check get_entities
  MBRange chkverts, chktris, chkhexes;
  check( manager.get_entities( MBVERTEX, chkverts ) );
  check( manager.get_entities( MBTRI   , chktris  ) );
  check( manager.get_entities( MBHEX   , chkhexes ) );
  ensure( chkverts.size() == vertices.size() && 
          std::equal(chkverts.begin(), chkverts.end(), vertices.begin() ) );
  ensure( chktris.size() == tris.size() && 
          std::equal(chktris.begin(), chktris.end(), tris.begin() ) );
  ensure( chkhexes.empty()     );
  
  // check get_number_entities
  MBEntityID num_vtx, num_tri, num_hex;
  check( manager.get_number_entities( MBVERTEX, num_vtx ) );
  check( manager.get_number_entities( MBTRI,    num_tri ) );
  check( manager.get_number_entities( MBHEX,    num_hex ) );
  ensure( num_vtx >= 0 &&  (unsigned)num_vtx == chkverts.size() );
  ensure( num_tri >= 0 &&  (unsigned)num_tri == chktris .size() );
  ensure( num_hex == 0 );
  
  return 0;
}

#endif


