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
#pragma warning (disable : 4786)
#endif

#include "MeshTopoUtil.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"

    //! generate all the AEntities bounding the vertices
MBErrorCode MeshTopoUtil::construct_aentities(const MBRange &vertices) 
{
  MBRange out_range;
  MBErrorCode result;
  result = mbImpl->get_adjacencies(vertices, 1, true, out_range, MBInterface::UNION);
  if (MB_SUCCESS != result) return result;
  out_range.clear();
  result = mbImpl->get_adjacencies(vertices, 2, true, out_range, MBInterface::UNION);
  if (MB_SUCCESS != result) return result;
  out_range.clear();
  result = mbImpl->get_adjacencies(vertices, 3, true, out_range, MBInterface::UNION);

  return result;
}

    //! given an entity, get its average position (avg vertex locations)
MBErrorCode MeshTopoUtil::get_average_position(const MBEntityHandle entity,
                                               double *avg_position) 
{
  const MBEntityHandle *connect;
  int num_connect;
  if (MBVERTEX == mbImpl->type_from_handle(entity))
    return mbImpl->get_coords(&entity, 1, avg_position);
    
  MBErrorCode result = mbImpl->get_connectivity(entity, connect, num_connect);
  if (MB_SUCCESS != result) return result;
  double dum_pos[3];
  avg_position[0] = avg_position[1] = avg_position[2] = 0.0;
  for (int i = 0; i < num_connect; i++) {
    result = mbImpl->get_coords(connect+i, 1, dum_pos);
    if (MB_SUCCESS != result) return result;
    avg_position[0] += dum_pos[0]; 
    avg_position[1] += dum_pos[1]; 
    avg_position[2] += dum_pos[2]; 
  }
  avg_position[0] /= (double) num_connect;
  avg_position[1] /= (double) num_connect;
  avg_position[2] /= (double) num_connect;

  return MB_SUCCESS;
}

  // given an entity, find the entities of next higher dimension around
  // that entity, ordered by connection through next higher dimension entities; 
  // if any of the star entities is in only entity of next higher dimension, 
  // on_boundary is returned true
MBErrorCode MeshTopoUtil::star_entities(const MBEntityHandle star_center,
                                        std::vector<MBEntityHandle> &star_entities,
                                        bool &bdy_entity,
                                        const MBEntityHandle starting_star_entity,
                                        std::vector<MBEntityHandle> *star_entities_dp1,
                                        MBRange *star_candidates_dp1)
{
    // now start the traversal
  bdy_entity = false;
  MBEntityHandle last_entity = starting_star_entity, last_dp1 = 0, next_entity, next_dp1;
  std::vector<MBEntityHandle> star_dp1;

  do {
      // get the next star entity
    MBErrorCode result = star_next_entity(star_center, last_entity, last_dp1,
                                          star_candidates_dp1,
                                          next_entity, next_dp1);
    if (MB_SUCCESS != result) return result;
    
      // if we're at a bdy and bdy_entity hasn't been set yet, we're at the
      // first bdy; reverse the lists and start traversing in the other direction; but,
      // pop the last star entity off the list and find it again, so that we properly
      // check for next_dp1
    if (0 == next_dp1 && !bdy_entity) {
      star_entities.push_back(next_entity);
      bdy_entity = true;
      std::reverse(star_entities.begin(), star_entities.end());
      std::reverse(star_dp1.begin(), star_dp1.end());
      star_entities.pop_back();
      last_entity = star_entities.back();
      last_dp1 = star_dp1.back();
    }
      // else if we're not on the bdy and next_entity is already in star, that means
      // we've come all the way around; don't put next_entity on list again, and
      // zero out last_dp1 to terminate while loop
    else if (!bdy_entity && 
             std::find(star_entities.begin(), star_entities.end(), next_entity) != 
             star_entities.end())
    {
      last_dp1 = 0;
    }

      // else, just assign last entities seen and go on to next iteration
    else {
      star_entities.push_back(next_entity);
      if (0 != next_dp1) star_dp1.push_back(next_dp1);
      last_entity = next_entity;
      last_dp1 = next_dp1;
    }
  }
  while (0 != last_dp1);
  
    // copy over the star_dp1 list, if requested
  if (NULL != star_entities_dp1) 
    (*star_entities_dp1).swap(star_dp1);
  
  return MB_SUCCESS;
}

MBErrorCode MeshTopoUtil::star_next_entity(const MBEntityHandle star_center,
                                           const MBEntityHandle last_entity,
                                           const MBEntityHandle last_dp1,
                                           MBRange *star_candidates_dp1,
                                           MBEntityHandle &next_entity,
                                           MBEntityHandle &next_dp1) 
{
    // given a star_center, a last_entity (whose dimension should be 1 greater than center)
    // and last_dp1 (dimension 2 higher than center), returns the next star entity across
    // last_dp1, and the next dp1 entity sharing next_entity; if star_candidates is non-empty,
    // star must come from those
  MBRange from_ents, to_ents;
  from_ents.insert(star_center);
  if (0 != last_dp1) from_ents.insert(last_dp1);
  int dim = mbImpl->dimension_from_handle(star_center);
  
  MBErrorCode result = mbImpl->get_adjacencies(from_ents, dim+1, false, to_ents);
  if (MB_SUCCESS != result) return result;
  
    // remove last_entity from result, and should only have 1 left, if any
  if (0 != last_entity) to_ents.erase(last_entity);
  if (!to_ents.empty()) next_entity = *to_ents.begin();
  else {
    next_entity = 0;
    next_dp1 = 0;
    return MB_SUCCESS;
  }
  
    // get next_dp1
  if (0 != star_candidates_dp1) to_ents = *star_candidates_dp1;
  else to_ents.clear();
  
  result = mbImpl->get_adjacencies(&next_entity, 1, dim+2, false, to_ents);
  if (MB_SUCCESS != result) return result;

    // can't be last one
  if (0 != last_dp1) to_ents.erase(last_dp1);
  
  if (!to_ents.empty()) next_dp1 = *to_ents.begin();

    // could be zero, means we're at bdy
  else next_dp1 = 0;

  return MB_SUCCESS;
}

    //! get "bridge" or "2nd order" adjacencies, going through dimension bridge_dim
MBErrorCode MeshTopoUtil::get_bridge_adjacencies(const MBEntityHandle from_entity,
                                                 const int bridge_dim,
                                                 const int to_dim,
                                                 MBRange &to_adjs) 
{
    // get pointer to connectivity for this entity
  const MBEntityHandle *connect;
  int num_connect;
  MBErrorCode result = mbImpl->get_connectivity(from_entity, connect, num_connect);
  if (MB_SUCCESS != result) return result;
  
  MBEntityType from_type = TYPE_FROM_HANDLE(from_entity);
  if (from_type >= MBENTITYSET) return MB_FAILURE;

  int from_dim = MBCN::Dimension(from_type);
  
  MBRange to_ents;

  if (MB_SUCCESS != result) return result;

  if (bridge_dim < from_dim) {
      // looping over each sub-entity of dimension bridge_dim...
    static MBEntityHandle bridge_verts[MB_MAX_SUB_ENTITIES];
    for (int i = 0; i < MBCN::NumSubEntities(from_type, bridge_dim); i++) {

        // get the vertices making up this sub-entity
      int num_bridge_verts;
      MBCN::SubEntityConn(connect, from_type, bridge_dim, i, &bridge_verts[0], num_bridge_verts);
    
        // get the to_dim entities adjacent
      to_ents.clear();
      MBErrorCode tmp_result = mbImpl->get_adjacencies(bridge_verts, num_bridge_verts,
                                                       to_dim, false, to_ents, MBInterface::INTERSECT);
      if (MB_SUCCESS != tmp_result) result = tmp_result;
    
      to_adjs.merge(to_ents);
    }

  }

    // now get the direct ones too, or only in the case where we're 
    // going to higher dimension for bridge
  MBRange bridge_ents, tmp_ents;
  tmp_ents.insert(from_entity);
  MBErrorCode tmp_result = mbImpl->get_adjacencies(tmp_ents, bridge_dim,
                                                   false, bridge_ents, 
                                                   MBInterface::UNION);
  if (MB_SUCCESS != tmp_result) return tmp_result;
    
  tmp_result = mbImpl->get_adjacencies(bridge_ents, to_dim, false, to_adjs, 
                                       MBInterface::UNION);
  if (MB_SUCCESS != tmp_result) return tmp_result;
  
    // if to_dimension is same as that of from_entity, make sure from_entity isn't
    // in list
  if (to_dim == from_dim) to_adjs.erase(from_entity);
  
  return result;
}

    //! return a common entity of the specified dimension, or 0 if there isn't one
MBEntityHandle MeshTopoUtil::common_entity(const MBEntityHandle ent1,
                                           const MBEntityHandle ent2,
                                           const int dim) 
{
  MBRange tmp_range, tmp_range2;
  tmp_range.insert(ent1); tmp_range.insert(ent2);
  MBErrorCode result = mbImpl->get_adjacencies(tmp_range, dim, false, tmp_range2);
  if (MB_SUCCESS == result || tmp_range2.empty()) return 0;
  else return *tmp_range2.begin();
}

