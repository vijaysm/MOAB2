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

/**\file MBFileReadTool.cpp
 *\author Jason Kraftcheck <kraftche@cae.wisc.edu>
 *\date 2007-08-21
 */

#include "MBFileReadTool.hpp"
#include "MBError.hpp"
#include <set>
#include "MBInternals.hpp"
#include "EntitySequence.hpp"
#include "PolyEntitySequence.hpp"
#include "EntitySequenceManager.hpp"

MBFileReadToolBase::MBFileReadToolBase( MBCore* moab )
  : mMOAB( moab ), fileSet( 0 ),
    parallelRank( moab->proc_config().rank() )
{ }

void MBFileReadToolBase::report_error( std::string err )
{
  mMOAB->get_error_handler()->set_last_error( err );
}

void MBFileReadToolBase::report_error( const char* fmt, ... )
{
  va_list args;
  va_start(args, fmt);
  mMOAB->get_error_handler()->set_last_error( fmt, args );
  va_end( args );
}

MBErrorCode MBFileReadToolBase::delete_all_entities()
{
  MBErrorCode rval;
  MBRange list;
  
  if (fileSet == 0)  // all entities already deleted
    return MB_SUCCESS;
  
  rval = moab()->get_entities_by_handle( fileSet, list );
  if (MB_SUCCESS != rval)
    return rval;
    
  list.insert( fileSet );
  rval = moab()->delete_entities( list );
  if (MB_SUCCESS != rval)
    return rval;
  
  fileSet = 0;
  return MB_SUCCESS;
}

MBErrorCode MBFileReadToolBase::delete_non_tagged_entities(
                                              const char* tagname, 
                                              const int* tag_values, 
                                              int num_tag_values )
{
  MBErrorCode rval;
  
    // no entities?
  if (fileSet == 0) {
    report_error( "No entities with tag \"%s\" in file.\n", tagname );
    return MB_FAILURE;
  }
  
    // get tag
  MBTag tag;
  rval = moab()->tag_get_handle( tagname, tag );
  if (MB_SUCCESS != rval) {
    report_error( "Tag not found: \"%s\"\n", tagname );
    return rval;
  }
  
    // check type of tag data
  MBDataType type;
  int size = 0;
  moab()->tag_get_size( tag, size );
  moab()->tag_get_data_type( tag, type );
  if (size != sizeof(int) || (type != MB_TYPE_INTEGER && type != MB_TYPE_OPAQUE)) {
    report_error( "Invalid type for tag \"%s\"\n", tagname );
    return MB_TYPE_OUT_OF_RANGE;
  }

    // construct argument arrays
  std::vector<MBTag> taghandles( num_tag_values ? num_tag_values : 1, tag );
  std::vector<const void*> tagvals( num_tag_values );
  for (int i = 0; i < num_tag_values; ++i)
    tagvals[i] = tag_values + i;
  
    // get entity sets with specified tag
  MBRange sets;  
  rval = moab()->get_entities_by_type_and_tag( fileSet,
                                               MBENTITYSET,
                                               &taghandles[0],
                                               num_tag_values ? &tagvals[0] : 0,
                                               taghandles.size(),
                                               sets,
                                               MBInterface::UNION );

    // get closure of entity sets
  MBRange entities, temp;
  MBRange::iterator iter;
    // recursively descend all containd and child entity sets
  std::set<MBEntityHandle> visited;
  std::vector<MBEntityHandle> children;
  while (!sets.empty()) {
    MBEntityHandle entset = sets.pop_front();
    if (!(visited.insert(entset).second))
      continue;
    entities.insert( entset );
    
    temp.clear();
    moab()->get_entities_by_handle( entset, temp );
    iter = temp.lower_bound( MBENTITYSET );
    entities.merge( temp.begin(), iter );
    std::copy( iter, temp.end(), mb_range_inserter( sets ) );
    
    children.clear();
    rval = moab()->get_child_meshsets( entset, children );
    if (MB_SUCCESS != rval)
      return rval;
    std::copy( children.begin(), children.end(), mb_range_inserter( sets ) );
  }
    // get faces for all polyhedra
  const MBEntityHandle* conn;
  int len; 
  iter = entities.lower_bound( MBPOLYHEDRON );
  for (; iter != entities.end() && TYPE_FROM_HANDLE(*iter) < MBENTITYSET; ++iter) {
    rval = moab()->get_connectivity( *iter, conn, len );
    if (MB_SUCCESS != rval)
      return rval;
    std::copy( conn, conn+len, mb_range_inserter( entities ) );
  }
    // get vertices for all elements
  iter = entities.lower_bound( MBEDGE );
  for (; iter != entities.end() && TYPE_FROM_HANDLE(*iter) < MBPOLYHEDRON; ++iter) {
    rval = moab()->get_connectivity( *iter, conn, len );
    if (MB_SUCCESS != rval)
      return rval;
    std::copy( conn, conn+len, mb_range_inserter( entities ) );
  }
  
    // We now have the list of entities to keep.
    // Get the list of entities to remove.
  MBRange all, dead;
  moab()->get_entities_by_handle( fileSet, all );
  dead = all.subtract( entities );
  
    // Destroy the entities
  rval = moab()->remove_entities( fileSet, dead );
  if (MB_SUCCESS != rval)
    return rval;
  rval = moab()->delete_entities( dead );
  if (MB_SUCCESS != rval)
    return rval;
  
  return MB_SUCCESS;
}

MBErrorCode MBFileReadToolBase::create_sequence( MBEntityID start_id,
                                                 MBEntityID count,
                                                 MBEntityType type,
                                                 MBEntityID conn_len,
                                                 const int* proc_id,
                                                 MBEntitySequence*& seq )
{
  MBErrorCode rval;
  
    // if we haven't created a set for all file entities yet,
    // do it now
  if (!fileSet) {
    rval = moab()->create_meshset( MESHSET_SET, fileSet );
    if (MB_SUCCESS != rval) {
      fileSet =0;
      return rval;
    }
  }
  
  MBEntityHandle start_handle;
  const int proc = proc_id ? *proc_id : parallel_rank();
  rval = mMOAB->sequence_manager()->create_entity_sequence(
                                                    type, 
                                                    count,
                                                    conn_len,
                                                    start_id,
                                                    proc,
                                                    start_handle,
                                                    seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  MBRange new_handles( seq->get_start_handle(), seq->get_end_handle() );
  assert( new_handles.size() == (size_t)count );
  rval = moab()->add_entities( fileSet, new_handles );
  if (MB_SUCCESS != rval) {
    delete seq;
    return rval;
  }
  
  return MB_SUCCESS;
}


MBErrorCode MBFileReadToolBase::create_node_arrays( 
                               MBEntityID preferred_start_id,
                               MBEntityID count,
                               MBEntityHandle& start_handle,
                               double *& x, double *& y, double *& z,
                               const int* proc_id )
{
  MBEntitySequence* seq = 0;
  MBErrorCode rval = create_sequence( preferred_start_id,
                                      count,
                                      MBVERTEX,
                                      0,
                                      proc_id, 
                                      seq );
  if (MB_SUCCESS != rval)
    return rval;
    
  VertexEntitySequence* vseq = static_cast<VertexEntitySequence*>(seq);
  vseq->get_coordinate_arrays( x, y, z );
  
  start_handle = seq->get_start_handle();
  return MB_SUCCESS;
}

MBErrorCode MBFileReadToolBase::create_element_array( 
                                              MBEntityID preferred_start_id,
                                              MBEntityID count,
                                              MBEntityType type,
                                              MBEntityID conn_length,
                                              MBEntityHandle& start_handle,
                                              MBEntityHandle*& conn_array,
                                              int** last_index_array,
                                              const int* proc_id )
{
    // caller needs last_index_array for poly-type elements
  if (type == MBPOLYGON || type == MBPOLYHEDRON) {
    if (!last_index_array)
      return MB_FAILURE;
  }
    // no last_index_array for other element types
  else if (last_index_array)
    *last_index_array = 0;
    
  MBEntitySequence* seq = 0;
  MBErrorCode rval = create_sequence( preferred_start_id,
                                      count,
                                      type,
                                      conn_length,
                                      proc_id, 
                                      seq );
  if (MB_SUCCESS != rval)
    return rval;
    
  if (type == MBPOLYGON || type == MBPOLYHEDRON) {
    PolyEntitySequence* pseq = static_cast<PolyEntitySequence*>(seq);
    pseq->get_connectivity_array( conn_array );
    pseq->get_index_array( *last_index_array );
  }
  else {
    ElementEntitySequence* eseq = static_cast<ElementEntitySequence*>(seq);
    eseq->get_connectivity_array( conn_array );
  }
  
  start_handle = seq->get_start_handle();
  return MB_SUCCESS;
}  
  
MBErrorCode MBFileReadToolBase::create_meshset_block( MBEntityID start_id,
                                                      MBEntityID count,
                                                      MBEntityHandle& start_handle,
                                                      unsigned flags,
                                                      const int* proc_id )
{
  std::vector<unsigned> flag_array( count, flags );
  return create_meshset_block( start_id, count, start_handle, &flag_array[0], proc_id );
}
  
MBErrorCode MBFileReadToolBase::create_meshset_block( MBEntityID start_id,
                                                      MBEntityID count,
                                                      MBEntityHandle& start_handle,
                                                      const unsigned* flags,
                                                      const int* proc_id )
{
  MBErrorCode rval;
  
  const int proc = proc_id ? *proc_id : parallel_rank();
  MBEntitySequence* seq;
  rval = mMOAB->sequence_manager()->create_meshset_sequence( count, 
                                                             start_id,
                                                             proc,
                                                             flags,
                                                             seq );
  if (MB_SUCCESS != rval)
    return rval;
  
  start_handle = seq->get_start_handle();
  return MB_SUCCESS;
}

MBErrorCode MBFileReadToolBase::set_tag_data( MBTag tag,
                                              MBEntityHandle beg,
                                              MBEntityHandle num,
                                              const void* data )
{
  if (num == 0)
    return MB_SUCCESS;
    
  MBRange range( beg, beg + num - 1 );
  return moab()->tag_set_data( tag, range, data );
}


