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
#include "SequenceManager.hpp"
#include "VertexSequence.hpp"
#include "ElementSequence.hpp"

#define RR if (MB_SUCCESS != result) return result

MBReadUtil::MBReadUtil(MBCore* mdb, MBError* error_handler) 
    : MBReadUtilIface(), mMB(mdb), mError(error_handler)
{
}

unsigned  MBReadUtil::parallel_rank() const
  { return mMB->proc_rank(); }

MBErrorCode MBReadUtil::get_node_arrays(
    const int /*num_arrays*/,
    const int num_nodes, 
    const int preferred_start_id,
    const int preferred_start_proc,
    MBEntityHandle& actual_start_handle, 
    std::vector<double*>& arrays)
{

  MBErrorCode error;
  EntitySequence* seq = 0;

  // create an entity sequence for these nodes 
  error = mMB->sequence_manager()->create_entity_sequence(
    MBVERTEX, num_nodes, 0, preferred_start_id, 
    preferred_start_proc, actual_start_handle,
    seq);

  if(error != MB_SUCCESS)
    return error;

  arrays.resize(3);

  error = static_cast<VertexSequence*>(seq)->get_coordinate_arrays(arrays[0], arrays[1], arrays[2]);
  for (unsigned i = 0; i< arrays.size(); ++i)
    if (arrays[i])
      arrays[i] += (actual_start_handle - seq->start_handle());
  
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
  EntitySequence* seq;
  
//  if (mdb_type <= MBVERTEX || mdb_type >= MBPOLYHEDRON || mdb_type == MBPOLYGON)
//    return MB_TYPE_OUT_OF_RANGE;

  // make an entity sequence to hold these elements
  unsigned proc_id = proc < 0 ? parallel_rank() : proc;
  error = mMB->sequence_manager()->create_entity_sequence(
      mdb_type, num_elements, verts_per_element, preferred_start_id, 
      proc_id, actual_start_handle, seq);
  if (MB_SUCCESS != error)
    return error;

  // get an array for the connectivity
  array = static_cast<ElementSequence*>(seq)->get_connectivity_array();
  if (!array)
    return MB_FAILURE;
  array += (actual_start_handle - seq->start_handle()) 
         * static_cast<ElementSequence*>(seq)->nodes_per_element();

  return error;
  
}

MBErrorCode MBReadUtil::create_entity_sets( MBEntityID num_sets,
                                            const unsigned* flags ,
                                            MBEntityID start_id,
                                            int proc,
                                            MBEntityHandle& start_handle )
{
  MBErrorCode error;
  EntitySequence* seq;
  error = mMB->sequence_manager()->create_meshset_sequence( num_sets,
                                                            start_id,
                                                            proc,
                                                            flags,
                                                            start_handle,
                                                            seq );
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

MBErrorCode MBReadUtil::gather_related_ents(MBRange &partition,
                                            MBRange &related_ents,
                                            MBRange *all_sets) 
{
    // first, related ents includes the partition itself
  related_ents.merge(partition);
  
    // loop over any sets, getting contained ents
  std::pair<MBRange::const_iterator, MBRange::const_iterator> pair_it =
    partition.equal_range(MBENTITYSET);

  MBErrorCode result = MB_SUCCESS;
  for (MBRange::const_iterator rit = pair_it.first; 
       rit != pair_it.second; rit++) {
    MBErrorCode tmp_result = 
      mMB->get_entities_by_handle(*rit, related_ents, 
                                  MBInterface::UNION);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  RR;

    // gather adjacent ents of lower dimension
  MBRange tmp_ents;
  for (int dim = 2; dim >= 0; dim--) {
    MBEntityType lower_type = MBCN::TypeDimensionMap[dim+1].first,
      upper_type = MBCN::TypeDimensionMap[3].second;
    
    MBRange::const_iterator bit = related_ents.lower_bound(lower_type),
      eit = related_ents.upper_bound(upper_type);
    MBRange from_ents;
    from_ents.merge(bit, eit);
    tmp_ents.clear();
    MBErrorCode tmp_result = mMB->get_adjacencies(from_ents, dim, false, 
                                                  tmp_ents, 
                                                  MBInterface::UNION);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
    else related_ents.merge(tmp_ents);
  }
  RR;
  
    // get related sets
  MBRange tmp_ents3;
  if (!all_sets) all_sets = &tmp_ents3;
  result = mMB->get_entities_by_type(0, MBENTITYSET, *all_sets);
  for (MBRange::iterator rit = all_sets->begin(); 
       rit != all_sets->end(); rit++) {
    tmp_ents.clear();
    result = mMB->get_entities_by_handle(*rit, tmp_ents, true); RR;
    MBRange tmp_ents2 = tmp_ents.intersect(related_ents);
    
      // if the intersection is not empty, set is related
    if (!tmp_ents2.empty()) related_ents.insert(*rit);
  }

  return MB_SUCCESS;
}
