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

#include "MBReadUtil.hpp"
#include "MBCore.hpp"
#include "AEntityFactory.hpp"
#include "MBError.hpp"
#include "EntitySequenceManager.hpp"
#include "PolyEntitySequence.hpp"


MBReadUtil::MBReadUtil(MBCore* mdb, MBError* error_handler) 
    : MBReadUtilIface(), mMB(mdb), mError(error_handler)
{
}

unsigned  MBReadUtil::parallel_rank() const
  { return mMB->proc_config().rank(); }
   
  /** Update reference counts for connectivity links.  Does nothing if 
   *  not compiled with reference counting enabled.   
   */
MBErrorCode MBReadUtil::increment_reference_count( 
                               const MBEntityHandle* ent_array,
                               size_t num_ent )
{
#ifdef MOAB_WITH_REFCOUNT
  return mMB->increment_reference_count( ent_array, num_ent );
#else
  return MB_SUCCESS;
#endif
}

MBErrorCode MBReadUtil::get_node_arrays(
    const int /*num_arrays*/,
    const int num_nodes, 
    const int preferred_start_id,
    const int preferred_start_proc,
    MBEntityHandle& actual_start_handle, 
    std::vector<double*>& arrays)
{

  MBErrorCode error;
  MBEntitySequence* seq = 0;

  MBEntityHandle preferred_start_handle;
  static int err;
  preferred_start_handle = CREATE_HANDLE(MBVERTEX, mMB->proc_config().id(preferred_start_id, 
                                         preferred_start_proc), err);
 
  // create an entity sequence for these nodes 
  error = mMB->sequence_manager()->create_entity_sequence(
      MBVERTEX, num_nodes, 0, preferred_start_handle, preferred_start_proc, actual_start_handle,
      seq);

  if(error != MB_SUCCESS)
    return error;

  arrays.resize(3);

  error = static_cast<VertexEntitySequence*>(seq)->get_coordinate_arrays(arrays[0], arrays[1], arrays[2]);
  
  return error;
}

MBErrorCode MBReadUtil::get_element_array(
    const int num_elements, 
    const int verts_per_element,
    const MBEntityType mdb_type,
    const int preferred_start_id, 
    const int proc, 
    MBEntityHandle& actual_start_handle, 
    MBEntityHandle*& array)
{

  MBErrorCode error;
  MBEntitySequence* seq;
  
  if (mdb_type <= MBVERTEX || mdb_type >= MBPOLYHEDRON || mdb_type == MBPOLYGON)
    return MB_TYPE_OUT_OF_RANGE;

  // make an entity sequence to hold these elements
  unsigned proc_id = proc < 0 ? parallel_rank() : proc;
  error = mMB->sequence_manager()->create_entity_sequence(
      mdb_type, num_elements, verts_per_element, preferred_start_id, 
      proc_id, actual_start_handle, seq);
  if (MB_SUCCESS != error)
    return error;

  // get an array for the connectivity
  error = static_cast<ElementEntitySequence*>(seq)->get_connectivity_array(array);

  return error;
  
}

MBErrorCode MBReadUtil::get_poly_element_array(
      const int num_poly, 
      const int conn_list_length,
      const MBEntityType mdb_type,
      const int preferred_start_id, 
      const int preferred_start_proc,
      MBEntityHandle& actual_start_handle, 
      int*& last_index_array,
      MBEntityHandle*& connectivity_array )
{

  MBErrorCode error;
  MBEntitySequence* seq;
  
  if (mdb_type != MBPOLYGON && mdb_type != MBPOLYHEDRON)
    return MB_TYPE_OUT_OF_RANGE;

  error = mMB->sequence_manager()->create_entity_sequence(
      mdb_type, num_poly, conn_list_length, preferred_start_id, 
      preferred_start_proc, actual_start_handle, seq);
  if (MB_SUCCESS != error || NULL == seq)
    return error;

  PolyEntitySequence *pseq = dynamic_cast<PolyEntitySequence*>(seq);
  assert(NULL != pseq);
  
  pseq->get_connectivity_array( connectivity_array );
  pseq->get_index_array( last_index_array );
  return MB_SUCCESS;
}

MBErrorCode MBReadUtil::create_entity_sets( MBEntityID num_sets,
                                            const unsigned* flags ,
                                            MBEntityID start_id,
                                            int proc,
                                            MBEntityHandle& start_handle )
{
  MBErrorCode error;
  MBEntitySequence* seq;
  error = mMB->sequence_manager()->create_meshset_sequence( num_sets,
                                                            start_id,
                                                            proc,
                                                            flags,
                                                            seq );
  if (seq)
    start_handle = seq->get_start_handle();
  return error;
}


MBErrorCode MBReadUtil::update_adjacencies(
      const MBEntityHandle start_handle,
      const int number_elements,
      const int number_vertices_per_element,
      const MBEntityHandle* conn_array)
{

  MBEntityHandle tmp_hndl = start_handle;
  AEntityFactory* adj_fact = mMB->a_entity_factory();

  // iterator over the elements and update adjacency information
  if(adj_fact != NULL && adj_fact->vert_elem_adjacencies())
  {
    int j=0;
    for(int i=0; i<number_elements; i++)
    {
      adj_fact->notify_create_entity( tmp_hndl, (conn_array+j), number_vertices_per_element);
      tmp_hndl++;
      j+=number_vertices_per_element;
    }
  }
  return MB_SUCCESS;
}



MBErrorCode MBReadUtil::report_error( const std::string& error )
{
  if(mError)
    return mError->set_last_error(error);
  else
    return MB_FAILURE;
}


MBErrorCode MBReadUtil::report_error( const char* error, ... )
{
  va_list args;
  va_start(args, error);
  MBErrorCode result = mError->set_last_error(error, args);
  va_end(args);
  return result;
}


