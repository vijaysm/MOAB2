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

#include "MBSkinner.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <iostream>
#include "MBUtil.hpp"
#include "MBTagConventions.hpp"

#define SKINNER_PI 3.1415926535897932384626




MBSkinner::~MBSkinner()
{
  // delete the adjacency tag
}


void MBSkinner::initialize()
{
  // go through and mark all the target dimension entities
  // that already exist as not deleteable
  // also get the connectivity tags for each type
  // also populate adjacency information
  MBEntityType type;
  MBDimensionPair target_ent_types = MBCN::TypeDimensionMap[mTargetDim];

  void* null_ptr = NULL;
  MBErrorCode result;

  if(mAdjTag == 0)
  {
    result = thisMB->tag_create("skinner adj", sizeof(void*), MB_TAG_DENSE, mAdjTag, &null_ptr);
    assert(MB_SUCCESS == result);
  }

  result = thisMB->tag_create("skinner deletable", 1, MB_TAG_BIT, mDeletableMBTag, NULL);
  assert(MB_SUCCESS == result);
  
  MBRange entities;
 
  // go through each type at this dimension 
  for(type = target_ent_types.first; type <= target_ent_types.second; ++type)
  {
    // get the entities of this type in the MB
    thisMB->get_entities_by_type(0, type, entities);

    // go through each entity of this type in the MB
    // and set its deletable tag to NO
    MBRange::iterator iter, end_iter;
    end_iter = entities.end();
    for(iter = entities.begin(); iter != end_iter; ++iter)
    {
      unsigned char bit = 0x1;
      result = thisMB->tag_set_data(mDeletableMBTag, &(*iter), 1, &bit);
      assert(MB_SUCCESS == result);
      // add adjacency information too
      add_adjacency(*iter);
    }
  }
}

void MBSkinner::deinitialize()
{
  MBErrorCode result = thisMB->tag_delete( mDeletableMBTag);
  assert(MB_SUCCESS == result);

  // remove the adjaceny tag
  result = thisMB->tag_delete(mAdjTag);
  assert(MB_SUCCESS == result);
  mAdjTag = NULL;
}


void MBSkinner::add_adjacency(MBEntityHandle entity)
{
  std::vector<MBEntityHandle> *adj = NULL;
  const MBEntityHandle *nodes;
  int num_nodes;
  MBErrorCode result = thisMB->get_connectivity(entity, nodes, num_nodes);
  assert(MB_SUCCESS == result);
  const MBEntityHandle *iter =
    std::min_element(nodes, nodes+num_nodes);

  if(iter == nodes+num_nodes)
    return;

  // add this entity to the node
  if(thisMB->tag_get_data(mAdjTag, iter, 1, &adj) == MB_SUCCESS && adj != NULL)
  {
    adj->push_back(entity);    
  }
  // create a new vector and add it
  else
  {
    adj = new std::vector<MBEntityHandle>;
    adj->push_back(entity);
    result = thisMB->tag_set_data(mAdjTag, iter, 1, &adj);
    assert(MB_SUCCESS == result);
  }
}

void MBSkinner::add_adjacency(MBEntityHandle entity, 
                               const MBEntityHandle *nodes,
                               const int num_nodes)
{
  std::vector<MBEntityHandle> *adj = NULL;
  const MBEntityHandle *iter = 
    std::min_element(nodes, nodes+num_nodes);

  if(iter == nodes+num_nodes)
    return;

  // add this entity to the node
  if(thisMB->tag_get_data(mAdjTag, iter, 1, &adj) == MB_SUCCESS && adj != NULL)
  {
    adj->push_back(entity);    
  }
  // create a new vector and add it
  else
  {
    adj = new std::vector<MBEntityHandle>;
    adj->push_back(entity);
    MBErrorCode result = thisMB->tag_set_data(mAdjTag, iter, 1, &adj);
    assert(MB_SUCCESS == result);
  }
}

MBErrorCode MBSkinner::find_geometric_skin(MBRange &forward_target_entities) 
{
    // attempts to find whole model skin, using geom topo sets first then
    // normal find_skin function
  bool debug = true;

    // look for geom topo sets
  MBTag geom_tag;
  MBErrorCode result = thisMB->tag_create(GEOM_DIMENSION_TAG_NAME, 4, 
                                            MB_TAG_SPARSE, geom_tag, NULL);

  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
    return result;
  
    // get face sets (dimension = 2)
  MBRange face_sets;
  int two = 2;
  const void *two_ptr = &two;
  result = thisMB->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag, &two_ptr, 1,
                                                 face_sets);

  MBRange::iterator it;
  if (MB_SUCCESS != result)
    return result;
  else if (face_sets.empty())
    return MB_ENTITY_NOT_FOUND;

    // ok, we have face sets; use those to determine skin
  MBRange skin_sets;
  if (debug) std::cout << "Found " << face_sets.size() << " face sets total..." << std::endl;
  
  for (it = face_sets.begin(); it != face_sets.end(); it++) {
    int num_parents;
    result = thisMB->num_parent_meshsets(*it, &num_parents);
    if (MB_SUCCESS != result)
      return result;
    else if (num_parents == 1)
      skin_sets.insert(*it);
  }

  if (debug) std::cout << "Found " << skin_sets.size() << " 1-parent face sets..." << std::endl;

  if (skin_sets.empty())
    return MB_FAILURE;
      
    // ok, we have the shell; gather up the elements, putting them all in forward for now
  for (it = skin_sets.begin(); it != skin_sets.end(); it++) {
    result = thisMB->get_entities_by_handle(*it, forward_target_entities, true);
    if (MB_SUCCESS != result) 
      return result;
  }
        
  return result;
}

MBErrorCode MBSkinner::find_skin(const MBRange &source_entities,
                                   MBRange &forward_target_entities,
                                   MBRange &reverse_target_entities)
{
  if(source_entities.empty())
    return MB_FAILURE;
  
  // get our working dimensions
  MBEntityType type = thisMB->type_from_handle(*(source_entities.begin()));
  const int source_dim = MBCN::Dimension(type);
  mTargetDim = source_dim - 1;

  // make sure we can handle the working dimensions
  if(mTargetDim < 0 || source_dim > 3)
    return MB_FAILURE;

  // initialize
  initialize();

  MBRange::const_iterator iter, end_iter;
  end_iter = source_entities.end();
  const MBEntityHandle *conn;
  MBEntityHandle match;

  direction direct;
  MBErrorCode result;
    // assume we'll never have more than 32 vertices on a facet (checked
    // with assert later)
  static MBEntityHandle sub_conn[32];
  int num_nodes;

  // for each source entity
  for(iter = source_entities.begin(); iter != end_iter; ++iter)
  {
    // get the connectivity of this entity
    result = thisMB->get_connectivity(*iter, conn, num_nodes);
    assert(MB_SUCCESS == result);
    
    type = thisMB->type_from_handle(*iter);
    MBRange::iterator seek_iter;
    
    // get connectivity of each n-1 dimension entity
    const struct MBCN::ConnMap* conn_map = &(MBCN::mConnectivityMap[type][mTargetDim-1]);
    for(int i=0; i<conn_map->num_sub_elements; i++)
    {
      int num_sub_nodes = conn_map->num_nodes_per_sub_element[i];
      assert(num_sub_nodes <= 32);
      for(int j=0; j<num_sub_nodes; j++)
        sub_conn[j] = conn[conn_map->conn[i][j]];
      
      // see if we can match this connectivity with
      // an existing entity
      find_match( conn_map->target_type[i], sub_conn, num_sub_nodes, match, direct );
  
      // if there is no match, create a new entity
      if(match == 0)
      {
        MBEntityHandle tmphndl=0;
        result = thisMB->create_element(conn_map->target_type[i], sub_conn, num_sub_nodes,
                                tmphndl);
        assert(MB_SUCCESS == result);
        add_adjacency(tmphndl, sub_conn, num_sub_nodes);
        forward_target_entities.insert(tmphndl);
      }
      // if there is a match, delete the matching entity
      // if we can. 
      else
      {
        if ( (seek_iter = forward_target_entities.find(match)) != forward_target_entities.end())
        {
          forward_target_entities.erase(seek_iter);
          remove_adjacency(match);
          if(entity_deletable(match))
          {
            result = thisMB->delete_entities(&match, 1);
            assert(MB_SUCCESS == result);
          }
        }
        else if ( (seek_iter = reverse_target_entities.find(match)) != reverse_target_entities.end())
        {
          reverse_target_entities.erase(seek_iter);
          remove_adjacency(match);
          if(entity_deletable(match))
          {
            result = thisMB->delete_entities(&match, 1);
            assert(MB_SUCCESS == result);
          }
        }
        else
        {
          if(direct == FORWARD)
          {
            forward_target_entities.insert(match);
          }
          else
          {
            reverse_target_entities.insert(match);
          }
        }
        
      }
    }
  }

  deinitialize();

  return MB_SUCCESS;
}


void MBSkinner::find_match( MBEntityType type, 
                             const MBEntityHandle *conn,
                             const int num_nodes,
                             MBEntityHandle& match,
                             MBSkinner::direction &direct)
{
  match = 0;

  const MBEntityHandle *iter = std::min_element(conn, conn+num_nodes);

  std::vector<MBEntityHandle> *adj = NULL;

  MBErrorCode result = thisMB->tag_get_data(mAdjTag, iter, 1, &adj);
  if(result == MB_FAILURE || adj == NULL)
  {
    return;
  }

  std::vector<MBEntityHandle>::iterator jter, end_jter;
  end_jter = adj->end();

  const MBEntityHandle *tmp;
  int num_verts;

  for(jter = adj->begin(); jter != end_jter; ++jter)
  {
    MBEntityType tmp_type;
    tmp_type = thisMB->type_from_handle(*jter);

    if( type != tmp_type )
      continue;

    result = thisMB->get_connectivity(*jter, tmp, num_verts);
    assert(MB_SUCCESS == result && num_verts >= num_nodes);
    if(connectivity_match(conn, tmp, num_verts, direct))
    {
      match = *jter;
      break;
    }        
  }
}

bool MBSkinner::connectivity_match( const MBEntityHandle *conn1,
                                     const MBEntityHandle *conn2,
                                     const int num_verts,
                                     MBSkinner::direction &direct)
{
  const MBEntityHandle *iter =
    std::find(conn2, conn2+num_verts, conn1[0]);
  if(iter == conn2+num_verts)
    return false;

  bool they_match = true;

  int i;
  unsigned int j = iter - conn2;
    
  // first compare forward
  for(i = 1; i<num_verts; ++i)
  {
    if(conn1[i] != conn2[(j+i)%num_verts])
    {
      they_match = false;
      break;
    }
  }
  
  if(they_match == true)
  {
    direct = FORWARD;
    return true;
  }
  
  they_match = true;
  
  // then compare reverse
  j += num_verts;
  for(i = 1; i < num_verts; )
  {
    if(conn1[i] != conn2[(j-i)%num_verts])
    {
      they_match = false;
      break;
    }
    ++i;
  }
  if (they_match)
  {
    direct = REVERSE;
  }
  return they_match;
}

  
void MBSkinner::remove_adjacency(MBEntityHandle entity)
{
  std::vector<MBEntityHandle> nodes, *adj = NULL;
  MBErrorCode result = thisMB->get_connectivity(&entity, 1, nodes);
  assert(MB_SUCCESS == result);
  std::vector<MBEntityHandle>::iterator iter = 
    std::min_element(nodes.begin(), nodes.end());

  if(iter == nodes.end())
    return;

  // remove this entity from the node
  if(thisMB->tag_get_data(mAdjTag, &(*iter), 1, &adj) == MB_SUCCESS && adj != NULL)
  {
    iter = std::find(adj->begin(), adj->end(), entity);
    if(iter != adj->end())
      adj->erase(iter);
  }
}

bool MBSkinner::entity_deletable(MBEntityHandle entity)
{
  unsigned char deletable=0;
  MBErrorCode result = thisMB->tag_get_data(mDeletableMBTag, &entity, 1, &deletable);
  assert(MB_SUCCESS == result);
  if(deletable == 1)
    return false;
  return true;
}

MBErrorCode MBSkinner::classify_2d_boundary( const MBRange &boundary,
                                               const MBRange &bar_elements,
                                               MBEntityHandle boundary_edges,
                                               MBEntityHandle inferred_edges,
                                               MBEntityHandle non_manifold_edges,
                                               MBEntityHandle other_edges,
                                               int &number_boundary_nodes)
{
  MBRange bedges, iedges, nmedges, oedges;
  MBErrorCode result = classify_2d_boundary(boundary, bar_elements,
                                             bedges, iedges, nmedges, oedges,
                                             number_boundary_nodes);
  if (MB_SUCCESS != result) return result;
  
    // now set the input meshsets to the output ranges
  result = thisMB->clear_meshset(&boundary_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(boundary_edges, bedges);
  if (MB_SUCCESS != result) return result;

  result = thisMB->clear_meshset(&inferred_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(inferred_edges, iedges);
  if (MB_SUCCESS != result) return result;

  result = thisMB->clear_meshset(&non_manifold_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(non_manifold_edges, nmedges);
  if (MB_SUCCESS != result) return result;

  result = thisMB->clear_meshset(&other_edges, 1);
  if (MB_SUCCESS != result) return result;
  result = thisMB->add_entities(other_edges, oedges);
  if (MB_SUCCESS != result) return result;

  return MB_SUCCESS;
}

MBErrorCode MBSkinner::classify_2d_boundary( const MBRange &boundary,
                                               const MBRange &bar_elements,
                                               MBRange &boundary_edges,
                                               MBRange &inferred_edges,
                                               MBRange &non_manifold_edges,
                                               MBRange &other_edges,
                                               int &number_boundary_nodes)
{

  // clear out the edge lists

  boundary_edges.clear();
  inferred_edges.clear();
  non_manifold_edges.clear();
  other_edges.clear();

  number_boundary_nodes = 0;

  // make sure we have something to work with
  if(boundary.empty())
  {
    return MB_FAILURE;
  }
  
  // get our working dimensions
  MBEntityType type = thisMB->type_from_handle(*(boundary.begin()));
  const int source_dim = MBCN::Dimension(type);

  // make sure we can handle the working dimensions
  if(source_dim != 2)
  {
    return MB_FAILURE;
  }
  mTargetDim = source_dim - 1;

  // initialize
  initialize();

  // additional initialization for this routine
  // define a tag for MBEDGE which counts the occurances of the edge below
  // default should be 0 for existing edges, if any

  MBTag count_tag;
  int default_count = 0;
  MBErrorCode result = thisMB->tag_create("mdbskinner count edges", sizeof(int),
                                            MB_TAG_DENSE, count_tag, &default_count);
  assert(MB_SUCCESS == result);

 
  MBRange::const_iterator iter, end_iter;
  end_iter = boundary.end();

  std::vector<MBEntityHandle> conn;
  static MBEntityHandle sub_conn[32];
  MBEntityHandle match;

  MBRange edge_list;
  MBRange boundary_nodes;
  MBSkinner::direction direct;
  
  // now, process each entity in the boundary

  for(iter = boundary.begin(); iter != end_iter; ++iter)
  {
    // get the connectivity of this entity
    conn.clear();
    result = thisMB->get_connectivity(&(*iter), 1, conn, true);
    assert(MB_SUCCESS == result);

    // add node handles to boundary_node range
    std::copy(conn.begin(), conn.end(), mb_range_inserter(boundary_nodes));

    type = thisMB->type_from_handle(*iter);
    
    // get connectivity of each n-1 dimension entity (edge in this case)
    const struct MBCN::ConnMap* conn_map = &(MBCN::mConnectivityMap[type][0]);
    for(int i=0; i<conn_map->num_sub_elements; i++)
    {
      int num_sub_nodes = conn_map->num_nodes_per_sub_element[i];
      assert(num_sub_nodes <= 32);
      for(int j=0; j<num_sub_nodes; j++)
        sub_conn[j] = conn[conn_map->conn[i][j]];
      
      // see if we can match this connectivity with
      // an existing entity
      find_match( conn_map->target_type[i], sub_conn, num_sub_nodes, match, direct );
  
      // if there is no match, create a new entity
      if(match == 0)
      {
        MBEntityHandle tmphndl=0;
        result = thisMB->create_element(conn_map->target_type[i], sub_conn, 
                                         num_sub_nodes, tmphndl);
        assert(MB_SUCCESS == result);
        add_adjacency(tmphndl, sub_conn, num_sub_nodes);
        //target_entities.insert(tmphndl);
        edge_list.insert(tmphndl);
        int count;
        result = thisMB->tag_get_data(count_tag, &tmphndl, 1, &count);
        assert(MB_SUCCESS == result);
        count++;
        result = thisMB->tag_set_data(count_tag, &tmphndl, 1, &count);
        assert(MB_SUCCESS == result);

      }
      else
      {
        // We found a match, we must increment the count on the match
        int count;
        result = thisMB->tag_get_data(count_tag, &match, 1, &count);
        assert(MB_SUCCESS == result);
        count++;
        result = thisMB->tag_set_data(count_tag, &match, 1, &count);
        assert(MB_SUCCESS == result);

        // if the entity is not deletable, it was pre-existing in
        // the database.  We therefore may need to add it to the
        // edge_list.  Since it will not hurt the range, we add
        // whether it was added before or not
        if(!entity_deletable(match))
        {
          edge_list.insert(match);
        }
      }
    }
  }

  // Any bar elements in the model should be classified separately
  // If the element is in the skin edge_list, then it should be put in
  // the non-manifold edge list.  Edges not in the edge_list are stand-alone
  // bars, and we make them simply boundary elements

  if (!bar_elements.empty())
  {
    MBRange::iterator bar_iter;
    for(iter = bar_elements.begin(); iter != bar_elements.end(); ++iter)
    {
      MBEntityHandle handle = *iter;
      bar_iter = edge_list.find(handle);
      if (bar_iter != edge_list.end())
      {
        // it is in the list, erase it and put in non-manifold list
        edge_list.erase(bar_iter);
        non_manifold_edges.insert(handle);
      }
      else
      {
        // not in the edge list, make it a boundary edge
        boundary_edges.insert(handle);
      }
    }
  }

  // now all edges should be classified.  Go through the edge_list,
  // and put all in the appropriate lists

  MBRange::iterator edge_iter, edge_end_iter;
  edge_end_iter = edge_list.end();
  int count;
  for(edge_iter = edge_list.begin(); edge_iter != edge_end_iter; edge_iter++)
  {
    // check the count_tag
    result = thisMB->tag_get_data(count_tag, &(*edge_iter), 1, &count);
    assert(MB_SUCCESS == result);
    if (count == 1)
    {
      boundary_edges.insert(*edge_iter);
   }
    else if (count == 2)
    {
      other_edges.insert(*edge_iter);
    }
    else
    {
      non_manifold_edges.insert(*edge_iter);
    }
  }

  // find the inferred edges from the other_edge_list

  double min_angle_degrees = 20.0;
  find_inferred_edges(const_cast<MBRange&> (boundary), other_edges, inferred_edges, min_angle_degrees);

  // we now want to remove the inferred_edges from the other_edges

  MBRange temp_range;
 
  std::set_difference(other_edges.begin(), other_edges.end(),
                      inferred_edges.begin(), inferred_edges.end(),
                      mb_range_inserter(temp_range),
                      std::less<MBEntityHandle>() );

  other_edges = temp_range;

  // get rid of count tag and deinitialize

  result = thisMB->tag_delete(count_tag);
  assert(MB_SUCCESS == result);
  deinitialize();

  // set the node count
  number_boundary_nodes = boundary_nodes.size();

  return MB_SUCCESS;
} 

void MBSkinner::find_inferred_edges(MBRange &skin_boundary,
                                     MBRange &candidate_edges,
                                     MBRange &inferred_edges,
                                     double reference_angle_degrees)
{

  // mark all the entities in the skin boundary
  MBTag mark_tag;
  MBErrorCode result = thisMB->tag_create("find inferred edges mark", 1, MB_TAG_BIT, mark_tag, NULL);
  assert(MB_SUCCESS == result);
  for(MBRange::iterator mark_iter = skin_boundary.begin();
      mark_iter != skin_boundary.end(); ++mark_iter)
  {
    unsigned char bit = true;
    result = thisMB->tag_set_data(mark_tag, &(*mark_iter), 1, &bit);
    assert(MB_SUCCESS == result);
  }

  // find the cosine of the reference angle

  double reference_cosine = cos(reference_angle_degrees*SKINNER_PI/180.0);
  
  // check all candidate edges for an angle greater than the minimum

  MBRange::iterator iter, end_iter = candidate_edges.end();
  std::vector<MBEntityHandle> adjacencies;
  std::vector<MBEntityHandle>::iterator adj_iter;
  MBEntityHandle face[2];

  for(iter = candidate_edges.begin(); iter != end_iter; ++iter)
  {

    // get the 2D elements connected to this edge
    adjacencies.clear();
    result = thisMB->get_adjacencies(&(*iter), 1, 2, false, adjacencies);
    if (MB_SUCCESS != result) 
      continue;

    // there should be exactly two, that is why the edge is classified as nonBoundary
    // and manifold

    int faces_found = 0;
    for(adj_iter = adjacencies.begin(); adj_iter != adjacencies.end() && faces_found < 2; ++adj_iter)
    {
      // we need to find two of these which are in the skin
      unsigned char is_marked = 0;
      result = thisMB->tag_get_data(mark_tag, &(*adj_iter), 1, &is_marked);
      assert(MB_SUCCESS == result);
      if(is_marked)
      {
        face[faces_found] = *adj_iter;
        faces_found++;
      } 
    }

//    assert(faces_found == 2 || faces_found == 0);
    if (2 != faces_found) 
      continue;

    // see if the two entities have a sufficient angle

    if ( has_larger_angle(face[0], face[1], reference_cosine) )
    {
       inferred_edges.insert(*iter);
    }
  }
  
  result = thisMB->tag_delete(mark_tag);
  assert(MB_SUCCESS == result);
}

bool MBSkinner::has_larger_angle(MBEntityHandle &entity1,
                                 MBEntityHandle &entity2,
                                 double reference_angle_cosine)
{
  // compare normals to get angle.  We assume that the surface quads
  // which we test here will be approximately planar

  double norm[2][3];
  MBUtil::normal(thisMB, entity1, norm[0][0], norm[0][1], norm[0][2]);
  MBUtil::normal(thisMB, entity2, norm[1][0], norm[1][1], norm[1][2]);

  double cosine = norm[0][0] * norm[1][0] + norm[0][1] * norm[1][1] + norm[0][2] * norm[1][2];

  if (cosine < reference_angle_cosine)
  {
    return true;
  }


  return false;
}
