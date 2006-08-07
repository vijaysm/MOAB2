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

#include "HigherOrderFactory.hpp"
#include "EntitySequenceManager.hpp"
#include "AEntityFactory.hpp"
#include "MBCore.hpp"
#include "MBCN.hpp"
#include <assert.h>
#include <algorithm>

using namespace std;

HigherOrderFactory::HigherOrderFactory(MBCore* MB, MBInterface::HONodeAddedRemoved* function_object) 
  : mMB(MB), mHONodeAddedRemoved(function_object)
{
  initialize_map();
}
HigherOrderFactory::~HigherOrderFactory() {}

//bool HigherOrderFactory::mMapInitialized = false;

void HigherOrderFactory::initialize_map()
{
 // if(mMapInitialized)
  //  return;

  for(MBEntityType i=MBVERTEX; i<MBMAXTYPE; i++)
  {
    const MBCN::ConnMap& canon_map = MBCN::mConnectivityMap[i][0];
    unsigned char (&this_map)[8][8] = mNodeMap[i];
    int num_node = MBCN::VerticesPerEntity(i);
    for(int j=0; j<canon_map.num_sub_elements; j++)
    {
      unsigned char x = canon_map.conn[j][0];
      unsigned char y = canon_map.conn[j][1];
      this_map[x][y] = num_node;
      this_map[y][x] = num_node;
      num_node++;
    }
  }

  //mMapInitialized = true;
}

MBErrorCode HigherOrderFactory::convert(const MBEntityHandle meshset, const bool mid_edge_nodes, 
                                         const bool mid_face_nodes, const bool mid_volume_nodes)
{

  // TODO --  add some more code to prevent from splitting of entity sequences when we don't need to.
  // Say we have all hex8's in our mesh and 3 falses are passed in.  In the end, no conversion will
  // happen, but the sequences could still be split up.

  MBRange entities;
  mMB->get_entities_by_handle(meshset, entities, true);
  
  // find out what entity sequences we need to convert 
  // and where, if necessary, to split them

  EntitySequenceManager* seq_manager = mMB->sequence_manager();

  std::vector< std::pair<ElementEntitySequence*,MBRange> > sequence_range_pairs;

  // get all the sequences that have these entities
  // could speed this up a bit by using a MBRange::pair_iterator
  for(MBRange::iterator p_iter = entities.begin(); p_iter != entities.end(); ++p_iter)
  {
    if(TYPE_FROM_HANDLE(*p_iter) == MBVERTEX || TYPE_FROM_HANDLE(*p_iter) >= MBENTITYSET)
      continue;

    MBEntitySequence* seq = NULL;
    MBErrorCode rval = seq_manager->find(*p_iter, seq);
    if(MB_SUCCESS == rval)
    {
      if(sequence_range_pairs.empty() || seq != sequence_range_pairs.rbegin()->first)
      {
        sequence_range_pairs.push_back(
            std::pair<ElementEntitySequence*,MBRange>( static_cast<ElementEntitySequence*>(seq), MBRange() )
            );
      }
      sequence_range_pairs.rbegin()->second.insert(*p_iter);
    }
  }

  std::vector<ElementEntitySequence*> sequences_to_process;
  

  // iterate through each entity sequence and split up entity sequences where we need to
  // put results into sequences list to process
  std::vector< std::pair<ElementEntitySequence*, MBRange> >::iterator s_iter;
  for(s_iter = sequence_range_pairs.begin(); s_iter != sequence_range_pairs.end(); ++s_iter)
  {
    // find all the entities we weren't told to convert
    MBRange holes;

    // we have a hole at the beginning of the sequence
    if(*(s_iter->second.begin()) != s_iter->first->get_start_handle())
    {
      holes.insert(s_iter->first->get_start_handle(), *(s_iter->second.begin())-1);
    }

    //we have a hole at the end of the sequence
    if(*(s_iter->second.rbegin()) != s_iter->first->get_end_handle())
    {
      holes.insert(*(s_iter->second.rbegin())+1, s_iter->first->get_end_handle());
    }
    // any holes in the middle
    MBRange::pair_iterator r_iter = s_iter->second.begin();
    MBRange::pair_iterator prev_r_iter = r_iter;
    r_iter++;
    for(; r_iter != s_iter->second.end(); )
    {
      holes.insert(prev_r_iter->second+1, r_iter->first-1);
      ++r_iter;
      ++prev_r_iter;
    }

    //we have a hole at the end of the sequence
    if(*(s_iter->second.rbegin()) != s_iter->first->get_end_handle())
    {
      holes.insert(*(s_iter->second.rbegin())+1, s_iter->first->get_end_handle());
    }

    // remove entities from the holes list that are unallocated entities
    for(MBRange::iterator rm_iter = holes.begin(); rm_iter != holes.end(); )
    {
      if(!s_iter->first->is_valid_entity(*rm_iter))
        rm_iter = holes.erase(rm_iter);
      else
        ++rm_iter;
    }

    // no holes, convert the whole sequence
    if(holes.empty())
    {
      sequences_to_process.push_back(s_iter->first);
    }
    else
    {
      MBEntitySequence* new_seq=NULL;
      // now we know where we need to split

      r_iter = holes.begin();

      //split at beginning
      if(r_iter->first == s_iter->first->get_start_handle())
      {
        s_iter->first->split(r_iter->second+1, new_seq);
        s_iter->first = static_cast<ElementEntitySequence*>(new_seq);
        sequences_to_process.push_back(static_cast<ElementEntitySequence*>(new_seq));
      }
      else
      {
        sequences_to_process.push_back(s_iter->first);
        s_iter->first->split(r_iter->first, new_seq);
        s_iter->first = static_cast<ElementEntitySequence*>(new_seq);
        if(r_iter->second != new_seq->get_end_handle())
        {
          s_iter->first->split(r_iter->second+1, new_seq);
          s_iter->first = static_cast<ElementEntitySequence*>(new_seq);
          sequences_to_process.push_back(static_cast<ElementEntitySequence*>(new_seq));
        }
      }

      // split in middle
      ++r_iter;
      MBRange::pair_iterator next_r_iter = r_iter;
      ++next_r_iter;
      for( ; r_iter != holes.end(); )
      {
        s_iter->first->split(r_iter->first, new_seq);
        s_iter->first = static_cast<ElementEntitySequence*>(new_seq);
        if( (next_r_iter != holes.end()) || (next_r_iter == holes.end() && r_iter->second != new_seq->get_end_handle()) )
        {
          s_iter->first->split(r_iter->second+1, new_seq);
          s_iter->first = static_cast<ElementEntitySequence*>(new_seq);
          sequences_to_process.push_back(static_cast<ElementEntitySequence*>(new_seq));
        }
        ++r_iter;
        ++next_r_iter;
      }

    }
  }

  // now convert the entity sequences
  for(std::vector<ElementEntitySequence*>::iterator ziter = sequences_to_process.begin(); 
      ziter != sequences_to_process.end(); ++ziter)
  {
    convert_sequence(*ziter, mid_edge_nodes, mid_face_nodes, mid_volume_nodes);
  }


  return MB_SUCCESS;
}


MBErrorCode HigherOrderFactory::convert_sequence(ElementEntitySequence* seq, const bool mid_edge_nodes, 
                                                  const bool mid_face_nodes, const bool mid_volume_nodes)
{

  MBErrorCode status = MB_SUCCESS;

  bool temp_mid_edge = mid_edge_nodes;
  bool temp_mid_face = mid_face_nodes;
  bool temp_mid_volume = mid_volume_nodes;
  
  // lets make sure parameters are ok before we continue
  MBEntityType this_type = seq->get_type();
  if(this_type == MBEDGE)
    temp_mid_face = temp_mid_volume = false;
  if(this_type == MBTRI || this_type == MBQUAD)
    temp_mid_volume = false;

  MBTag deletable_nodes;
  mMB->tag_create("", 1, MB_TAG_BIT, deletable_nodes, NULL);

  // will return if we need to create them
  seq->convert_realloc(temp_mid_edge, temp_mid_face, temp_mid_volume, mMB, deletable_nodes);

  // gather nodes that were marked
  MBRange nodes;
  mMB->get_entities_by_type_and_tag(0, MBVERTEX, &deletable_nodes, NULL, 1, nodes);

  MBEntityHandle low_meshset;
  int dum;
  low_meshset = CREATE_HANDLE(MBENTITYSET, 0, dum);
  
  for(MBRange::iterator iter = nodes.begin(); iter != nodes.end(); ++iter)
  {
    unsigned char marked = 0;
    mMB->tag_get_data(deletable_nodes, &(*iter), 1, &marked);
    if(marked)
    {
      // we can delete it
      if(mHONodeAddedRemoved)
        mHONodeAddedRemoved->node_removed( *iter );
      mMB->delete_entities(&(*iter), 1);
    }
  }
  
  mMB->tag_delete(deletable_nodes);

  // create mid edge nodes if necessary
  if(temp_mid_edge)
    status = add_mid_edge_nodes(seq);
  // create mid face nodes if necessary
  if(temp_mid_face)
    status = add_mid_face_nodes(seq);
  // create mid volume nodes if necessary
  if(temp_mid_volume)
    status = add_mid_volume_nodes(seq);

  return status;

}


MBErrorCode HigherOrderFactory::add_mid_volume_nodes(ElementEntitySequence* seq)
{
  MBEntityType this_type = seq->get_type();
  EntitySequenceManager* seq_manager = mMB->sequence_manager();

  // find out where in the connectivity list to add these new mid volume nodes
  int edge_factor = seq->has_mid_edge_nodes() ? 1 : 0;
  int face_factor = seq->has_mid_face_nodes() ? 1 : 0;
  // offset by number of higher order nodes on edges if they exist
  int num_corner_nodes = MBCN::VerticesPerEntity(this_type);
  int new_node_index = num_corner_nodes;
  new_node_index += edge_factor * MBCN::mConnectivityMap[this_type][0].num_sub_elements;
  new_node_index += face_factor * MBCN::mConnectivityMap[this_type][1].num_sub_elements;

  MBEntityHandle* element = NULL;
  seq->get_connectivity_array(element);
  MBEntityHandle curr_handle = seq->get_start_handle();
  int nodes_per_element = seq->nodes_per_element();
  MBEntityHandle* end_element = element + nodes_per_element * (seq->number_allocated());

  // iterate over the elements
  for(; element < end_element; element+=nodes_per_element)
  {
    // this element isn't being used, skip it
    if(!seq->is_valid_entity(curr_handle))
    {
      curr_handle++;
      continue;
    }

    // find the centroid of this element
    double tmp_coords[3], sum_coords[3] = {0,0,0};
    MBEntitySequence* seq=NULL;
    for(int i=0; i<num_corner_nodes; i++)
    {
      seq_manager->find(element[i], seq);
      static_cast<VertexEntitySequence*>(seq)->get_coordinates(
          element[i], tmp_coords[0], tmp_coords[1], tmp_coords[2]
          );
      sum_coords[0] += tmp_coords[0];
      sum_coords[1] += tmp_coords[1];
      sum_coords[2] += tmp_coords[2];
    }
    sum_coords[0] /= num_corner_nodes;
    sum_coords[1] /= num_corner_nodes;
    sum_coords[2] /= num_corner_nodes;

    // create a new vertex at the centroid
    mMB->create_vertex(sum_coords, element[new_node_index]);
    
    if(mHONodeAddedRemoved)
      mHONodeAddedRemoved->node_added(element[new_node_index], curr_handle);

    curr_handle++;
  }

  return MB_SUCCESS;
}


MBErrorCode HigherOrderFactory::add_mid_face_nodes(ElementEntitySequence* seq)
{
  MBEntityType this_type = seq->get_type();
  EntitySequenceManager* seq_manager = mMB->sequence_manager();
  int num_vertices = MBCN::VerticesPerEntity(this_type);
  int num_edges = MBCN::mConnectivityMap[this_type][0].num_sub_elements;
  num_edges = seq->has_mid_edge_nodes() ? num_edges : 0;
  int num_faces = MBCN::mConnectivityMap[this_type][1].num_sub_elements;

  const MBCN::ConnMap& entity_faces = MBCN::mConnectivityMap[this_type][1];

  MBEntityHandle* element = NULL;
  seq->get_connectivity_array(element);
  MBEntityHandle curr_handle = seq->get_start_handle();
  int nodes_per_element = seq->nodes_per_element();
  MBEntityHandle* end_element = element + nodes_per_element * (seq->number_allocated());

  MBEntityHandle tmp_face_conn[4];  // max face nodes = 4
  std::vector<MBEntityHandle> adjacent_entities(4);

  double tmp_coords[3];
  
  // iterate over the elements
  for(; element < end_element; element+=nodes_per_element)
  {
    // this element isn't being used, skip it
    if(!seq->is_valid_entity(curr_handle))
    {
      curr_handle++;
      continue;
    }

    // for each edge in this entity
    for(int i=0; i<num_faces; i++)
    {
      // a node was already assigned
      if(element[i+num_edges+num_vertices] != 0)
        continue;

      tmp_face_conn[0] = element[entity_faces.conn[i][0]];
      tmp_face_conn[1] = element[entity_faces.conn[i][1]];
      tmp_face_conn[2] = element[entity_faces.conn[i][2]];
      if(entity_faces.num_nodes_per_sub_element[i] == 4)
        tmp_face_conn[3] = element[entity_faces.conn[i][3]];
      else
        tmp_face_conn[3] = 0;

      MBEntityHandle already_made_node = center_node_exist(tmp_face_conn, adjacent_entities);
      
      if(already_made_node)
      {
        element[i+num_edges+num_vertices] = already_made_node;
      }
      // create a node
      else
      {
        MBEntitySequence* tmp_sequence = NULL;
        double sum_coords[3] = {0,0,0};
        int max_nodes = entity_faces.num_nodes_per_sub_element[i];
        for(int k=0; k<max_nodes; k++)
        {
          seq_manager->find(tmp_face_conn[k], tmp_sequence);
          static_cast<VertexEntitySequence*>(tmp_sequence)->get_coordinates( 
              tmp_face_conn[k], tmp_coords[0], tmp_coords[1], tmp_coords[2]);
          sum_coords[0] += tmp_coords[0];
          sum_coords[1] += tmp_coords[1];
          sum_coords[2] += tmp_coords[2];
        }

        sum_coords[0] /= max_nodes;
        sum_coords[1] /= max_nodes;
        sum_coords[2] /= max_nodes;

        mMB->create_vertex(sum_coords, element[i+num_edges+num_vertices]);
      }

      if(mHONodeAddedRemoved)
        mHONodeAddedRemoved->node_added(element[i+num_edges+num_vertices], curr_handle);
      
    }

    curr_handle++;

  }

  return MB_SUCCESS;
}


MBErrorCode HigherOrderFactory::add_mid_edge_nodes(ElementEntitySequence* seq)
{
  // for each node, need to see if it was already created.
  MBEntityType this_type = seq->get_type();
  EntitySequenceManager* seq_manager = mMB->sequence_manager();

  // offset by number of corner nodes
  int num_vertices = MBCN::VerticesPerEntity(this_type);
  int num_edges = MBCN::mConnectivityMap[this_type][0].num_sub_elements;

  const MBCN::ConnMap& entity_edges = MBCN::mConnectivityMap[this_type][0];
  

  MBEntityHandle* element = NULL;
  seq->get_connectivity_array(element);
  MBEntityHandle curr_handle = seq->get_start_handle();
  int nodes_per_element = seq->nodes_per_element();
  MBEntityHandle* end_element = element + nodes_per_element * (seq->number_allocated());

  MBEntityHandle tmp_edge_conn[2];
  std::vector<MBEntityHandle> adjacent_entities(32);

  double tmp_coords[3];

  // iterate over the elements
  for(; element < end_element; element+=nodes_per_element)
  {
    // this element isn't being used, skip it
    if(!seq->is_valid_entity(curr_handle))
    {
      curr_handle++;
      continue;
    }

    // for each edge in this entity
    for(int i=0; i<num_edges; i++)
    {
      // a node was already assigned
      if(element[i+num_vertices] != 0)
        continue;

      tmp_edge_conn[0] = element[entity_edges.conn[i][0]];
      tmp_edge_conn[1] = element[entity_edges.conn[i][1]];

      MBEntityHandle already_made_node = center_node_exist(tmp_edge_conn[0], tmp_edge_conn[1],
          adjacent_entities);
      
      if(already_made_node)
      {
        element[i+num_vertices] = already_made_node;
      }
      // create a node
      else
      {
        MBEntitySequence* tmp_sequence = NULL;
        double sum_coords[3] = {0,0,0};
        seq_manager->find(tmp_edge_conn[0], tmp_sequence);
        static_cast<VertexEntitySequence*>(tmp_sequence)->get_coordinates( 
            tmp_edge_conn[0], tmp_coords[0], tmp_coords[1], tmp_coords[2]);
        sum_coords[0] += tmp_coords[0];
        sum_coords[1] += tmp_coords[1];
        sum_coords[2] += tmp_coords[2];
        seq_manager->find(tmp_edge_conn[1], tmp_sequence);
        static_cast<VertexEntitySequence*>(tmp_sequence)->get_coordinates( 
            tmp_edge_conn[1], tmp_coords[0], tmp_coords[1], tmp_coords[2]);
        sum_coords[0] = (sum_coords[0] + tmp_coords[0]) /2;
        sum_coords[1] = (sum_coords[1] + tmp_coords[1]) /2;
        sum_coords[2] = (sum_coords[2] + tmp_coords[2]) /2;

        mMB->create_vertex(sum_coords, element[i+num_vertices]);
      }

      if(mHONodeAddedRemoved)
        mHONodeAddedRemoved->node_added(element[i+num_vertices], curr_handle);
      
    }

    curr_handle++;

  }

  return MB_SUCCESS;
}

MBEntityHandle HigherOrderFactory::center_node_exist( MBEntityHandle corner1, 
    MBEntityHandle corner2, std::vector<MBEntityHandle>& adj_entities)
{
  AEntityFactory* a_fact = mMB->a_entity_factory();
  std::vector<MBEntityHandle> adj_corner1(32);
  std::vector<MBEntityHandle> adj_corner2(32);

  // create needed vertex adjacencies
  if (!a_fact->vert_elem_adjacencies())
    a_fact->create_vert_elem_adjacencies();

  // vectors are returned sorted
  
  a_fact->get_adjacencies(corner1, adj_corner1);
  a_fact->get_adjacencies(corner2, adj_corner2);

  // these are the entities adjacent to both nodes
  adj_entities.clear();
  std::set_intersection(adj_corner1.begin(), adj_corner1.end(), adj_corner2.begin(),
      adj_corner2.end(), std::back_inserter<std::vector<MBEntityHandle> >(adj_entities));


  // iterate of the entities to find a mid node
  const MBEntityHandle* conn;
  int conn_size = 0;
  for(std::vector<MBEntityHandle>::iterator iter = adj_entities.begin();
      iter != adj_entities.end(); )
  {
    MBEntityType this_type = TYPE_FROM_HANDLE(*iter);
    mMB->get_connectivity(*iter, conn, conn_size);
    // if this entity has mid edge nodes
    if(MBCN::HasMidEdgeNodes(this_type, conn_size))
    {
      // find out at which index the mid node should be at
      int first_node = std::find(conn, conn+conn_size, corner1) - conn;
      int second_node = std::find(conn, conn+conn_size, corner2) - conn;
      if(first_node == conn_size || second_node == conn_size)
        assert("We should always find our nodes no matter what" == NULL);
      int high_node_index = mNodeMap[this_type][first_node][second_node];
      if(conn[high_node_index] != 0)
        return conn[high_node_index];
      ++iter;
    }
    else
    {
      iter = adj_entities.erase(iter);
    }
  }

  return 0;

}

MBEntityHandle HigherOrderFactory::center_node_exist( MBEntityHandle corners[4], 
    std::vector<MBEntityHandle>& adj_entities)
{
  AEntityFactory* a_fact = mMB->a_entity_factory();
  std::vector<MBEntityHandle> adj_corner[4];
  int num_nodes = corners[3] == 0 ? 3 : 4;
  int i = 0;

  // create needed vertex adjacencies
  if (!a_fact->vert_elem_adjacencies())
    a_fact->create_vert_elem_adjacencies();

  // vectors are returned sorted
  for(i=0; i<num_nodes; i++)
    a_fact->get_adjacencies(corners[i], adj_corner[i]);

  // these are the entities adjacent to both nodes
  for(i=1; i<num_nodes; i++)
  {
    adj_entities.clear();
    std::set_intersection(adj_corner[i-1].begin(), adj_corner[i-1].end(), adj_corner[i].begin(),
      adj_corner[i].end(), std::back_inserter<std::vector<MBEntityHandle> >(adj_entities));
    adj_corner[i].swap(adj_entities);
  }
  adj_entities.swap(adj_corner[i-1]);
  

  // iterate of the entities to find a mid node
  const MBEntityHandle* conn;
  int conn_size = 0;
  for(std::vector<MBEntityHandle>::iterator iter = adj_entities.begin();
      iter != adj_entities.end(); )
  {
    MBEntityType this_type = TYPE_FROM_HANDLE(*iter);
    const MBCN::ConnMap& entity_faces = MBCN::mConnectivityMap[this_type][1];
    mMB->get_connectivity(*iter, conn, conn_size);
    int offset = MBCN::VerticesPerEntity(this_type);
    if(MBCN::HasMidEdgeNodes(this_type, conn_size))
      offset += MBCN::mConnectivityMap[this_type][0].num_sub_elements;

    // if this entity has mid face nodes
    if(MBCN::HasMidFaceNodes(this_type, conn_size))
    {
      int k;
      int indexes[4];
      for(k=0; k<num_nodes; k++)
        indexes[k] = std::find(conn, conn+conn_size, corners[k]) - conn;
      
      // find out at which index the mid node should be at
      for(k=0; k<entity_faces.num_sub_elements; k++)
      {
        if(MBCN::VerticesPerEntity(entity_faces.target_type[k]) != num_nodes)
          continue;

        int* pivot = std::find(indexes, indexes+num_nodes, entity_faces.conn[k][0]);
        if(pivot == indexes+num_nodes)
          continue;

        if(pivot != indexes)
          std::rotate(indexes, pivot, indexes+num_nodes);

        if(std::equal(indexes, indexes+num_nodes, entity_faces.conn[k]))
        {
          if(conn[k+offset] != 0)
            return conn[k+offset];
          k=entity_faces.num_sub_elements;
        }
        else
        {
          int temp = indexes[1];
          indexes[1] = indexes[num_nodes-1];
          indexes[num_nodes-1] = temp;
          if(std::equal(indexes, indexes+num_nodes, entity_faces.conn[k]))
          {
            if(conn[k+offset] != 0)
              return conn[k+offset];
            k=entity_faces.num_sub_elements;
          }
        }
      }
      ++iter;
    }
    else
    {
      iter = adj_entities.erase(iter);
    }
  }

  return 0;

}

bool HigherOrderFactory::add_center_node(MBEntityType this_type, MBEntityHandle* element_conn, 
    int conn_size, MBEntityHandle corner_node1, MBEntityHandle corner_node2, 
    MBEntityHandle center_node)
{
  int first_node = std::find(element_conn, element_conn+conn_size, corner_node1) - element_conn;
  int second_node = std::find(element_conn, element_conn+conn_size, corner_node2) - element_conn;
  if(first_node == conn_size || second_node == conn_size)
    assert("We should always find our nodes no matter what" == NULL);
  int high_node_index = mNodeMap[this_type][first_node][second_node];
  element_conn[high_node_index] = center_node;  
  return true;
}

