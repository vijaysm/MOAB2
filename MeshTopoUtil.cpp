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
                                        MBRange *star_entities_candidates_dp1)
{
  MBEntityHandle last_entity = 0, last_dp1 = 0;
  MBErrorCode result;
  MBRange from_ents, adj_ents;
  unsigned int star_dim = MBCN::Dimension(mbImpl->type_from_handle(star_center))+1;
  std::vector<MBEntityHandle> star_dp1;
  
  if (0 != starting_star_entity) {
    last_entity = starting_star_entity;
    star_entities.push_back(last_entity);
  }
  
    // now start the traversal
  bdy_entity = false;
  bool first = true;
  
  while (first || 0 != last_entity) {
    first = false;
    from_ents.clear();
    from_ents.insert(star_center);
    
      // get dp1 entities which are adjacent to star_center and last_entity
    if (0 != last_entity) from_ents.insert(last_entity);
      // if candidates were specified, also needs to be among those
    if (NULL != star_entities_candidates_dp1) 
      adj_ents = *star_entities_candidates_dp1;
    else adj_ents.clear();
    
    MBErrorCode tmp_result = mbImpl->get_adjacencies(from_ents, star_dim+1, false, adj_ents);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      break;
    }
    if (adj_ents.empty()) {
      bdy_entity = true;
      return MB_SUCCESS;
    }

      // get a dp1 entity which isn't last_dp1
      // if both last_entity and last_dp1 are zero, we're just starting, so do nothing
    if (0 == last_entity && 0 == last_dp1);
    else if (*adj_ents.begin() != last_dp1) last_dp1 = *adj_ents.begin();
    else if (adj_ents.size() > 1) last_dp1 = *adj_ents.rbegin();
    else {
        // if we're here, we're at an open bdy; try reversing the list and going
        // in the other direction; also jumble the from_ents list so it appears we
        // last checked the new last_entity
      bdy_entity = true;
      std::reverse(star_entities.begin(), star_entities.end());
      std::reverse(star_dp1.begin(), star_dp1.end());
      if (0 != last_entity) from_ents.erase(last_entity);
      last_entity = *star_entities.rbegin();
      last_dp1 = *star_dp1.rbegin();
      from_ents.insert(last_entity);
    }

      // get the other star-dimension entity adj to star_center and last_dp1
    if (0 != last_entity) from_ents.erase(last_entity);
    if (0 != last_dp1) from_ents.insert(last_dp1);
    adj_ents.clear();
    tmp_result = mbImpl->get_adjacencies(from_ents, star_dim, false, adj_ents);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      break;
    }
    
      // next star entity is the one not equal to last_entity; must be one at this point
    if (*adj_ents.begin() != last_entity) last_entity = *adj_ents.begin();
    else last_entity = *adj_ents.rbegin();

    if (std::find(star_entities.begin(), star_entities.end(), last_entity) != 
        star_entities.end()) {
        // either we're back where we started or one after; either way, we're done
      last_entity = 0;
        // if we're not a bdy entity, the last dp1 entity wasn't put on list
      if (0 != last_dp1) star_dp1.push_back(last_dp1);
    }
    else {
        // add star and last_dp1 to the list
      star_entities.push_back(last_entity);
      if (0 != last_dp1) star_dp1.push_back(last_dp1);
    }
  } // end while

    // copy over the star_dp1 list, if requested
  if (NULL != star_entities_dp1) 
    (*star_entities_dp1).swap(star_dp1);
  
  return MB_SUCCESS;
}

