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


EntitySequenceManager::EntitySequenceManager()
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
                                                        const int hint_start_id,
                                                        MBEntityHandle &start_handle,
                                                        MBEntitySequence *&seq) 
{
  int this_dim = MBCN::Dimension(type);

    // use > instead of != in the following assert to also catch cases where imin > imax, etc.
  assert((this_dim < 3 || kmax > kmin) &&
         (this_dim < 2 || jmax > jmin) &&
         (this_dim < 1 || imax > imin));

    // compute # entities; not as easy as it would appear...
  int num_ent;
  if (MBVERTEX == type)
    num_ent = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
  else {
    num_ent = (imax-imin) *
      (this_dim >= 2 ? (jmax-jmin) : 1) *
      (this_dim >= 3 ? (kmax-kmin) : 1);
  }
  
    // get a start handle
  start_handle = get_start_handle(hint_start_id, MB_PROC_RANK, type, num_ent);
  assert(0 != start_handle);
  if (0 == start_handle) return MB_FAILURE;
  
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

  
/*!
  creates an entity sequence based on number of entities and type.
  uses the hint_start as the start id for the entity handles if possible
  returns the actual start handle and the entity sequence pointer
*/
MBErrorCode EntitySequenceManager::create_entity_sequence( MBEntityType type, int num_ent, int num_nodes,
                                                           int hint_start, int hint_start_proc, 
                                                           MBEntityHandle& start_handle, 
                                                            MBEntitySequence*& seq)
{
  start_handle = get_start_handle(hint_start, hint_start_proc, type, num_ent);
  
  // actually create the sequence
  return private_create_entity_sequence( start_handle, num_ent, num_nodes, true, seq);
}

MBEntityHandle EntitySequenceManager::get_start_handle(int hint_start, int hint_start_proc, MBEntityType type,
                                                        int num_ent) 
{
  // need to find unused space in the MBEntityHandle ID space
  int dum = 0;
  MBEntityHandle start_hint_handle = CREATE_HANDLE(type, hint_start, hint_start_proc, dum);
 
  // this is the first of this type we are making 
  if(mSequenceMap[type].empty()) {
    if (hint_start < MB_START_ID)
      return CREATE_HANDLE( type, MB_START_ID, hint_start_proc, dum );
    else
      return start_hint_handle;
  }

  // see if we can use the start hint
  std::map<MBEntityHandle, MBEntitySequence*>::iterator iter =
    mSequenceMap[type].upper_bound(start_hint_handle);

  MBEntityHandle last_of_type;
  if (iter == mSequenceMap[type].end())
    last_of_type = CREATE_HANDLE( type, MB_END_ID, MB_PROC_COUNT - 1, dum );
  else
    last_of_type = iter->second->get_start_handle() - 1;
  
  --iter;
  if (start_hint_handle > iter->second->get_end_handle()
   && start_hint_handle + num_ent - 1 <= last_of_type) 
    return start_hint_handle;
  
  return mSequenceMap[type].rbegin()->second->get_end_handle() + 1;
}

MBErrorCode EntitySequenceManager::private_create_entity_sequence(MBEntityHandle start,
                                                           int num_ent, int num_nodes,
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


MBErrorCode EntitySequenceManager::create_vertex(const double coords[3], MBEntityHandle& handle)
{
    
  // see if there is an existing sequence that can take this new vertex
  if(!mPartlyFullSequenceMap[MBVERTEX].empty())
  {
    MBEntitySequence* seq = mPartlyFullSequenceMap[MBVERTEX].begin()->second;
    handle = seq->get_unused_handle();
    
    static_cast<VertexEntitySequence*>(seq)->
      set_coordinates(handle, coords[0], coords[1], coords[2]);

    return MB_SUCCESS;
  }

  // we need to make a new entity sequence
  if(!mSequenceMap[MBVERTEX].empty())
  {
    handle = mSequenceMap[MBVERTEX].rbegin()->second->get_end_handle() + 1;
  }
  else
  {
    int err=0;
    handle = CREATE_HANDLE(MBVERTEX, MB_START_ID, err);
  }

  VertexEntitySequence* seq = new VertexEntitySequence(this, handle, 4096, false);
  handle = seq->get_unused_handle();
  seq->set_coordinates(handle, coords[0], coords[1], coords[2]);
  
  return MB_SUCCESS;

}


MBErrorCode EntitySequenceManager::create_element(MBEntityType type, 
                                                   const MBEntityHandle *conn, 
                                                   const int num_vertices,
                                                   MBEntityHandle& handle)
{
  // see if there is an existing sequence that can take this new element
  std::map<MBEntityHandle, MBEntitySequence*>::iterator iter;
  for(iter = mPartlyFullSequenceMap[type].begin();
      iter != mPartlyFullSequenceMap[type].end();
      ++iter)
  {
    ElementEntitySequence* seq = dynamic_cast<ElementEntitySequence*>(iter->second);
    if(seq->nodes_per_element() == (unsigned int) num_vertices ||
       seq->nodes_per_element() == 0)
    {
      if (MBPOLYGON == type || MBPOLYHEDRON == type) {
        return dynamic_cast<PolyEntitySequence*>(iter->second)->add_entity(conn, num_vertices, handle);
      }
      else {
        handle = seq->get_unused_handle();
        return seq->set_connectivity(handle, conn, num_vertices);
      }
    }

  }

  // we need to make a new entity sequence
  if(!mSequenceMap[type].empty())
  {
    handle = mSequenceMap[type].rbegin()->second->get_end_handle() + 1;
  }
  else
  {
    int err=0;
    handle = CREATE_HANDLE(type, MB_START_ID, err);
  }

  ElementEntitySequence* seq;
  if (MBPOLYGON == type || MBPOLYHEDRON == type)
    seq = new PolyEntitySequence(this, handle, 0, 0, false);
  else
    seq = new ElementEntitySequence(this, handle, 4096, num_vertices, false);

  handle = seq->get_unused_handle();
  MBErrorCode result = seq->set_connectivity(handle, conn, num_vertices);
  
  return result;
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

MBErrorCode EntitySequenceManager::get_number_entities(MBEntityType type, int& num_entities) const 
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
  EntitySequenceManager manager;
  MBEntitySequence* seq;
  
  // create some sequences
  const unsigned NV[] = { 100, 5, 1000 };
  MBEntityHandle vh[3];
  MBEntitySequence* vs[3];
  MBEntityHandle th;
  MBEntitySequence *ts;
  const unsigned TRI_START_ID = 3;
  const unsigned NT = 640;
  check(manager.create_entity_sequence( MBVERTEX, NV[0], 0, 1, 0, vh[0], vs[0] ));
  check(manager.create_entity_sequence( MBVERTEX, NV[1], 0, 1, 0, vh[1], vs[1] ));
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
  int num_vtx, num_tri, num_hex;
  check( manager.get_number_entities( MBVERTEX, num_vtx ) );
  check( manager.get_number_entities( MBTRI,    num_tri ) );
  check( manager.get_number_entities( MBHEX,    num_hex ) );
  ensure( num_vtx >= 0 &&  (unsigned)num_vtx == chkverts.size() );
  ensure( num_tri >= 0 &&  (unsigned)num_tri == chktris .size() );
  ensure( num_hex == 0 );
  
  return 0;
}

#endif


