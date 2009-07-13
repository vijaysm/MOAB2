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

MBErrorCode MBReadUtil::get_node_arrays(
    const int /*num_arrays*/,
    const int num_nodes, 
    const int preferred_start_id,
    MBEntityHandle& actual_start_handle, 
    std::vector<double*>& arrays)
{

  MBErrorCode error;
  EntitySequence* seq = 0;
  
  if (num_nodes < 1) {
    actual_start_handle = 0;
    arrays.clear();
    return MB_INDEX_OUT_OF_RANGE;
  }

  // create an entity sequence for these nodes 
  error = mMB->sequence_manager()->create_entity_sequence(
    MBVERTEX, num_nodes, 0, preferred_start_id, 
    actual_start_handle,
    seq);

  if(error != MB_SUCCESS)
    return error;

  if (seq->start_handle() > actual_start_handle ||
      seq->end_handle() < actual_start_handle ||
      seq->end_handle() - actual_start_handle + 1 < (unsigned)num_nodes)
    return MB_FAILURE;

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
    MBEntityHandle& actual_start_handle, 
    MBEntityHandle*& array)
{

  MBErrorCode error;
  EntitySequence* seq;
  
  if (num_elements < 1) {
    actual_start_handle = 0;
    array = 0;
    return MB_INDEX_OUT_OF_RANGE;
  }
  
//  if (mdb_type <= MBVERTEX || mdb_type >= MBPOLYHEDRON || mdb_type == MBPOLYGON)
//    return MB_TYPE_OUT_OF_RANGE;

  // make an entity sequence to hold these elements
  error = mMB->sequence_manager()->create_entity_sequence(
      mdb_type, num_elements, verts_per_element, preferred_start_id, 
      actual_start_handle, seq);
  if (MB_SUCCESS != error)
    return error;

  if (seq->start_handle() > actual_start_handle ||
      seq->end_handle() < actual_start_handle ||
      seq->end_handle() - actual_start_handle + 1 < (unsigned)num_elements)
    return MB_FAILURE;

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
                                            MBEntityHandle& start_handle )
{
  if (num_sets < 1) {
    start_handle = 0;
    return MB_INDEX_OUT_OF_RANGE;
  }

  MBErrorCode error;
  EntitySequence* seq;
  error = mMB->sequence_manager()->create_meshset_sequence( num_sets,
                                                            start_id,
                                                            flags,
                                                            start_handle,
                                                            seq );
  if (MB_SUCCESS != error)
    return error;

  if (seq->start_handle() > start_handle ||
      seq->end_handle() < start_handle ||
      seq->end_handle() - start_handle + 1 < (MBEntityHandle)num_sets)
    return MB_FAILURE;
    
  return MB_SUCCESS;
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

    // gather adjacent ents of other dimensions
  MBRange tmp_ents;
  for (int dim = 3; dim >= 0; dim--) {
    tmp_ents.clear();
    MBErrorCode tmp_result = mMB->get_adjacencies(related_ents, dim, false, 
                                                  tmp_ents, 
                                                  MBInterface::UNION);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
    else related_ents.merge(tmp_ents);
  }
  RR;
  
    // related ents includes the partition itself
  related_ents.merge(partition);
  
    // get contains-related sets
  MBRange tmp_ents3, last_related;
  if (!all_sets) all_sets = &tmp_ents3;
  result = mMB->get_entities_by_type(0, MBENTITYSET, *all_sets); RR;
  while (related_ents.size() != last_related.size()) {
    last_related = related_ents;
    for (MBRange::iterator rit = all_sets->begin(); 
         rit != all_sets->end(); rit++) {
      if (related_ents.find(*rit) != related_ents.end()) continue;
      
      tmp_ents.clear();
      result = mMB->get_entities_by_handle(*rit, tmp_ents, true); RR;
      MBRange tmp_ents2 = intersect( tmp_ents, related_ents);
    
        // if the intersection is not empty, set is related
      if (!tmp_ents2.empty()) related_ents.insert(*rit);
    }
  }
  
    // get parent/child-related sets
  last_related.clear();
  while (related_ents.size() != last_related.size()) {
    last_related = related_ents;
    std::pair<MBRange::const_iterator, MBRange::const_iterator> it_pair = 
      last_related.equal_range(MBENTITYSET);

    for (MBRange::const_iterator rit = it_pair.first;
         rit != it_pair.second; rit++) {
        // get all parents/children and add to related ents
      tmp_ents.clear();
      result = mMB->get_parent_meshsets(*rit, tmp_ents, 0); RR;
      related_ents.merge(tmp_ents);
      
      tmp_ents.clear();
      result = mMB->get_child_meshsets(*rit, tmp_ents, 0); RR;
      related_ents.merge(tmp_ents);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBReadUtil::get_ordered_vertices(MBEntityHandle *bound_ents,
                                             int *sense, 
                                             int bound_size,
                                             int dim,
                                             MBEntityHandle *bound_verts, 
                                             MBEntityType &etype) 
{
    // get dimension of bounding entities
  int bound_dim = MBCN::Dimension(TYPE_FROM_HANDLE(bound_ents[0]));
  int indices[MB_MAX_SUB_ENTITY_VERTICES];
  const MBEntityHandle *connect;
  std::vector<MBEntityHandle> tmp_connect;
  
    // find the right entity type based on # bounding ents
  int numv = 0, num_connect;
  MBErrorCode result;
  for (MBEntityType t = MBEDGE; t < MBENTITYSET; t++) {
    int nindex = MBCN::NumSubEntities(t, bound_dim);
    if (MBCN::Dimension(t) != dim || nindex != bound_size) 
      continue;

      // fill in vertices from bounding entity vertices
    int nverts = MBCN::VerticesPerEntity(t);
    std::fill(bound_verts, bound_verts+nverts, 0);
    for (int index = 0; index < nindex; index++) {
      result = mMB->get_connectivity(bound_ents[index], connect, num_connect,
                                     false, &tmp_connect);
      if (MB_SUCCESS != result) return result;
      
      MBCN::SubEntityVertexIndices(t, bound_dim, index, indices);

      for (int c = 0; c < num_connect; c++) {
        if (!bound_verts[indices[c]]) {
          bound_verts[indices[c]] = (sense[index] > 0) ?
              connect[c] : connect[num_connect - c - 1];
          numv++;
        }
      }
      if (numv == nverts) {
        etype = t;
        return MB_SUCCESS;
      }
    }
  }
  
    // if we get here, we didn't get full connectivity
  etype = MBMAXTYPE;
  return MB_FAILURE;
}
