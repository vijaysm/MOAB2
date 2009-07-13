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
#include "MBInterface.hpp"
#include "MBCN.hpp"

#include <assert.h>

#define RR {if (MB_SUCCESS != result) return result;}

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
MBErrorCode MeshTopoUtil::get_average_position(MBRange &entities,
                                               double *avg_position) 
{
  std::vector<MBEntityHandle> ent_vec;
  std::copy(entities.begin(), entities.end(), std::back_inserter(ent_vec));
  return get_average_position(&ent_vec[0], ent_vec.size(), avg_position);
}

//! given an entity, get its average position (avg vertex locations)
MBErrorCode MeshTopoUtil::get_average_position(const MBEntityHandle *entities,
                                               const int num_entities,
                                               double *avg_position) 
{
  double dum_pos[3];
  avg_position[0] = avg_position[1] = avg_position[2] = 0.0;

  MBRange connect;
  MBErrorCode result = mbImpl->get_adjacencies(entities, num_entities, 0, false,
                                               connect, MBInterface::UNION);
  if (MB_SUCCESS != result) return result;

  if (connect.empty()) return MB_FAILURE;
  
  for (MBRange::iterator rit = connect.begin(); rit != connect.end(); rit++) {
    result = mbImpl->get_coords(&(*rit), 1, dum_pos);
    if (MB_SUCCESS != result) return result;
    avg_position[0] += dum_pos[0]; 
    avg_position[1] += dum_pos[1]; 
    avg_position[2] += dum_pos[2]; 
  }
  avg_position[0] /= (double) connect.size();
  avg_position[1] /= (double) connect.size();
  avg_position[2] /= (double) connect.size();

  return MB_SUCCESS;
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

  return get_average_position(connect, num_connect, avg_position);
}

  // given an entity, find the entities of next higher dimension around
  // that entity, ordered by connection through next higher dimension entities; 
  // if any of the star entities is in only one entity of next higher dimension, 
  // on_boundary is returned true
MBErrorCode MeshTopoUtil::star_entities(const MBEntityHandle star_center,
                                        std::vector<MBEntityHandle> &star_entities,
                                        bool &bdy_entity,
                                        const MBEntityHandle starting_star_entity,
                                        std::vector<MBEntityHandle> *star_entities_dp2,
                                        MBRange *star_candidates_dp2)
{
    // now start the traversal
  bdy_entity = false;
  MBEntityHandle last_entity = starting_star_entity, last_dp2 = 0, next_entity, next_dp2;
  std::vector<MBEntityHandle> star_dp2;
  MBErrorCode result;
  int center_dim = mbImpl->dimension_from_handle(star_center);
  
  MBRange tmp_candidates_dp2;
  if (NULL != star_candidates_dp2) tmp_candidates_dp2 = *star_candidates_dp2;
  else {
    result = mbImpl->get_adjacencies(&star_center, 1, 
                                     center_dim+2,
                                     false, tmp_candidates_dp2);
    if (MB_SUCCESS != result) return result;
  }

  do {
      // get the next star entity
    result = star_next_entity(star_center, last_entity, last_dp2,
                              &tmp_candidates_dp2,
                              next_entity, next_dp2);
    if (MB_SUCCESS != result) return result;

      // special case: if starting_star_entity isn't connected to any entities of next
      // higher dimension, it's the only entity in the star; put it on the list and return
    if (star_entities.empty() && next_entity == 0 && next_dp2 == 0) {
      star_entities.push_back(last_entity);
      bdy_entity = true;
      return MB_SUCCESS;
    }

      // if we're at a bdy and bdy_entity hasn't been set yet, we're at the
      // first bdy; reverse the lists and start traversing in the other direction; but,
      // pop the last star entity off the list and find it again, so that we properly
      // check for next_dp2
    if (0 == next_dp2 && !bdy_entity) {
      star_entities.push_back(next_entity);
      bdy_entity = true;
      std::reverse(star_entities.begin(), star_entities.end());
      star_entities.pop_back();
      last_entity = star_entities.back();
      if (!star_dp2.empty()) {
        std::reverse(star_dp2.begin(), star_dp2.end());
        last_dp2 = star_dp2.back();
      }
    }
      // else if we're not on the bdy and next_entity is already in star, that means
      // we've come all the way around; don't put next_entity on list again, and
      // zero out last_dp2 to terminate while loop
    else if (!bdy_entity && 
             std::find(star_entities.begin(), star_entities.end(), next_entity) != 
             star_entities.end() &&
             (std::find(star_dp2.begin(), star_dp2.end(), next_dp2) != 
             star_dp2.end() || !next_dp2))
    {
      last_dp2 = 0;
    }

      // else, just assign last entities seen and go on to next iteration
    else {
      if (std::find(star_entities.begin(), star_entities.end(), next_entity) == 
          star_entities.end())
        star_entities.push_back(next_entity);
      if (0 != next_dp2) {
        star_dp2.push_back(next_dp2);
        tmp_candidates_dp2.erase(next_dp2);
      }
      last_entity = next_entity;
      last_dp2 = next_dp2;
    }
  }
  while (0 != last_dp2);
  
    // copy over the star_dp2 list, if requested
  if (NULL != star_entities_dp2) 
    (*star_entities_dp2).swap(star_dp2);
  
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

    // if no last_dp1, contents of to_ents should share dp1-dimensional entity with last_entity
  if (0 != last_entity && 0 == last_dp1) {
    MBRange tmp_to_ents;
    for (MBRange::iterator rit = to_ents.begin(); rit != to_ents.end(); rit++) {
      if (0 != common_entity(last_entity, *rit, dim+2))
        tmp_to_ents.insert(*rit);
    }
    to_ents = tmp_to_ents;
  }

  if (0 == last_dp1 && to_ents.size() > 1 && NULL != star_candidates_dp1 && 
      !star_candidates_dp1->empty()) {
      // if we have a choice of to_ents and no previous dp1 and dp1 candidates, 
      // the one we choose needs to be adjacent to one of the candidates
    result = mbImpl->get_adjacencies(*star_candidates_dp1, dim+1, false,
                                     from_ents, MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    to_ents = intersect( to_ents, from_ents);
  }
  
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

MBErrorCode MeshTopoUtil::star_entities_nonmanifold(const MBEntityHandle star_entity,
                                                    std::vector<std::vector<MBEntityHandle> > &stars,
                                                    std::vector<bool> *bdy_flags,
                                                    std::vector<std::vector<MBEntityHandle> > *dp2_stars) 
{
    // Get a series of (d+1)-dimensional stars around a d-dimensional entity, such that
    // each star is on a (d+2)-manifold containing the d-dimensional entity; each star
    // is either open or closed, and also defines a (d+2)-star whose entities are bounded by
    // (d+1)-entities on the star and on the (d+2)-manifold
    //
    // Algorithm:
    // get the (d+2)-manifold entities; for d=1 / d+2=3, just assume all connected elements, since
    //   we don't do 4d yet
    // get intersection of (d+1)-entities adjacent to star entity and union of (d+1)-entities 
    //   adjacent to (d+2)-manifold entities; these will be the entities in the star
    // while (d+1)-entities
    //   remove (d+1)-entity from (d+1)-entities
    //   get the (d+1)-star and (d+2)-star around that (d+1)-entity (using star_entities)
    //   save that star to the star list, and the bdy flag and (d+2)-star if requested
    //   remove (d+2)-entities from the (d+2)-manifold entities
    //   remove (d+1)-entities from the (d+1)-entities
    // (end while)

  int this_dim = mbImpl->dimension_from_handle(star_entity);
  if (3 <= this_dim || 0 > this_dim) return MB_FAILURE;
  
    // get the (d+2)-manifold entities; for d=1 / d+2=3, just assume all connected elements, since
    //   we don't do 4d yet
  MBRange dp2_manifold;
  MBErrorCode result = get_manifold(star_entity, this_dim+2, dp2_manifold);
  if (MB_SUCCESS != result) return result;
  
    // get intersection of (d+1)-entities adjacent to star and union of (d+1)-entities 
    //   adjacent to (d+2)-manifold entities; also add manifold (d+1)-entities, to catch
    //   any not connected to (d+2)-entities
  MBRange dp1_manifold;
  result = mbImpl->get_adjacencies(dp2_manifold, this_dim+1, false, dp1_manifold,
                                   MBInterface::UNION);
  if (MB_SUCCESS != result) return result;

  result = mbImpl->get_adjacencies(&star_entity, 1, this_dim+1, false, dp1_manifold);
  if (MB_SUCCESS != result) return result;

  result = get_manifold(star_entity, this_dim+1, dp1_manifold);
  if (MB_SUCCESS != result) return result;
  
    // while (d+1)-entities
  while (!dp1_manifold.empty()) {
    
      //   get (d+1)-entity from (d+1)-entities (don't remove it until after star,
      //     since the star entities must come from dp1_manifold)
    MBEntityHandle this_ent = *dp1_manifold.begin();
    
      //   get the (d+1)-star and (d+2)-star around that (d+1)-entity (using star_entities)
    std::vector<MBEntityHandle> this_star_dp1, this_star_dp2;
    bool on_bdy;
    result = star_entities(star_entity, this_star_dp1, on_bdy, this_ent, &this_star_dp2, 
                           &dp2_manifold);
    if (MB_SUCCESS != result) return result;

      // if there's no star entities, it must mean this_ent isn't bounded by any dp2 
      // entities (wasn't put into star in star_entities 'cuz we're passing in non-null
      // dp2_manifold above); put it in
    if (this_star_dp1.empty()) {
      MBRange dum_range;
      result = mbImpl->get_adjacencies(&this_ent, 1, this_dim+2, false, dum_range);
      if (MB_SUCCESS != result) return result;
      if (dum_range.empty()) this_star_dp1.push_back(this_ent);
    }
    
      // now we can remove it
    dp1_manifold.erase(dp1_manifold.begin());

      //   save that star to the star list, and the bdy flag and (d+2)-star if requested
    if (!this_star_dp1.empty()) {
      stars.push_back(this_star_dp1);
      if (NULL != bdy_flags) bdy_flags->push_back(on_bdy);
      if (NULL != dp2_stars) dp2_stars->push_back(this_star_dp2);
    }
    
      //   remove (d+2)-entities from the (d+2)-manifold entities
    for (std::vector<MBEntityHandle>::iterator vit = this_star_dp2.begin(); 
         vit != this_star_dp2.end(); vit++)
      dp2_manifold.erase(*vit);
    
      //   remove (d+1)-entities from the (d+1)-entities
    for (std::vector<MBEntityHandle>::iterator vit = this_star_dp1.begin(); 
         vit != this_star_dp1.end(); vit++)
      dp1_manifold.erase(*vit);
    
      // (end while)
  }

    // check for leftover dp2 manifold entities, these should be in one of the 
    // stars
  if (!dp2_manifold.empty()) {
    for (MBRange::iterator rit = dp2_manifold.begin(); rit != dp2_manifold.end(); rit++) {
    }
  }
    
  return MB_SUCCESS;
}

//! get (target_dim)-dimensional manifold entities connected to star_entity; that is,
//! the entities with <= 1 connected (target_dim+2)-dimensional adjacent entities;
//! for target_dim=3, just return all of them
//! just insert into the list, w/o clearing manifold list first
MBErrorCode MeshTopoUtil::get_manifold(const MBEntityHandle star_entity,
                                       const int target_dim,
                                       MBRange &manifold) 
{
    // get all the entities of target dimension connected to star
  MBRange tmp_range;
  MBErrorCode result = mbImpl->get_adjacencies(&star_entity, 1, target_dim, false, tmp_range);
  if (MB_SUCCESS != result) return result;
  
    // now save the ones which are (target_dim+1)-dimensional manifold;
    // for target_dim=3, just return whole range, since we don't do 4d
  if (target_dim == 3) {
    manifold.merge(tmp_range);
    return MB_SUCCESS;
  }
  
  for (MBRange::iterator rit = tmp_range.begin(); rit != tmp_range.end(); rit++) {
    MBRange dum_range;
      // get (target_dim+1)-dimensional entities
    result = mbImpl->get_adjacencies(&(*rit), 1, target_dim+1, false, dum_range);
    if (MB_SUCCESS != result) return result;
    
      // if there are only 1 or zero, add to manifold list
    if (1 >= dum_range.size()) manifold.insert(*rit);
  }
  
  return MB_SUCCESS;
}

    //! get "bridge" or "2nd order" adjacencies, going through dimension bridge_dim
MBErrorCode MeshTopoUtil::get_bridge_adjacencies(MBRange &from_entities,
                                                 int bridge_dim,
                                                 int to_dim, 
                                                 MBRange &to_ents,
                                                 int num_layers)
{
  MBRange bridge_ents, last_toents, new_toents(from_entities);
  MBErrorCode result;
  if (0 == num_layers || from_entities.empty()) return MB_FAILURE;
  
    // for each layer, get bridge-adj entities and accumulate
  for (int nl = 0; nl < num_layers; nl++) {
    MBRange new_bridges;
      // get bridge ents
    result = mbImpl->get_adjacencies(new_toents, bridge_dim, true, new_bridges,
                                     MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    
      // get to_dim adjacencies, merge into to_ents
    last_toents =  to_ents;
    if (-1 == to_dim) {
      result = mbImpl->get_adjacencies(new_bridges, 3, false, to_ents,
				       MBInterface::UNION);
      if (MB_SUCCESS != result) return result;
      for (int d = 2; d >= 1; d--) {
	result = mbImpl->get_adjacencies(to_ents, d, true, to_ents,
					 MBInterface::UNION);
	if (MB_SUCCESS != result) return result;
      }
    }
    else {
      result = mbImpl->get_adjacencies(new_bridges, to_dim, false, to_ents,
				       MBInterface::UNION);
      if (MB_SUCCESS != result) return result;
    }

      // subtract last_toents to get new_toents
    if (nl < num_layers-1)
      new_toents = subtract( to_ents, last_toents);
  }

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
  MBErrorCode result = MB_SUCCESS;
  MBEntityType from_type = TYPE_FROM_HANDLE(from_entity);
  if (from_type == MBVERTEX) {
    connect = &from_entity;
    num_connect = 1;
  }
  else {
    result = mbImpl->get_connectivity(from_entity, connect, num_connect);
    if (MB_SUCCESS != result) return result;
  }
  
  if (from_type >= MBENTITYSET) return MB_FAILURE;

  int from_dim = MBCN::Dimension(from_type);
  
  MBRange to_ents;

  if (MB_SUCCESS != result) return result;

  if (bridge_dim < from_dim) {
      // looping over each sub-entity of dimension bridge_dim...
    static MBEntityHandle bridge_verts[MB_MAX_SUB_ENTITIES];
    int bridge_indices[MB_MAX_SUB_ENTITIES];
    for (int i = 0; i < MBCN::NumSubEntities(from_type, bridge_dim); i++) {

        // get the vertices making up this sub-entity
      int num_bridge_verts = MBCN::VerticesPerEntity( MBCN::SubEntityType( from_type, bridge_dim, i ) );
      MBCN::SubEntityVertexIndices( from_type, bridge_dim, i, bridge_indices );
      for (int j = 0; j < num_bridge_verts; ++j)
        bridge_verts[j]= connect[bridge_indices[j]];
      //MBCN::SubEntityConn(connect, from_type, bridge_dim, i, &bridge_verts[0], num_bridge_verts);
    
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
  if (MB_SUCCESS != result || tmp_range2.empty()) return 0;
  else return *tmp_range2.begin();
}

  //! return the opposite side entity given a parent and bounding entity.
  //! This function is only defined for certain types of parent/child types;
  //! See MBCN.hpp::OppositeSide for details.
  //!
  //! \param parent The parent element
  //! \param child The child element
  //! \param opposite_element The index of the opposite element
MBErrorCode MeshTopoUtil::opposite_entity(const MBEntityHandle parent,
                                          const MBEntityHandle child,
                                          MBEntityHandle &opposite_element) 
{
    // get the side no.
  int side_no, offset, sense;
  MBErrorCode result = mbImpl->side_number(parent, child, side_no, 
                                           offset, sense);
  if (MB_SUCCESS != result) return result;
  
    // get the child index from MBCN
  int opposite_index, opposite_dim;
  int status = MBCN::OppositeSide(mbImpl->type_from_handle(parent),
                                  side_no, mbImpl->dimension_from_handle(child),
                                  opposite_index, opposite_dim);
  if (0 != status) return MB_FAILURE;
  
    // now get the side element from MOAB
  result = mbImpl->side_element(parent, opposite_dim, opposite_index, 
                                opposite_element);
  if (MB_SUCCESS != result) return result;
  
  return MB_SUCCESS;
}

MBErrorCode MeshTopoUtil::split_entities_manifold(MBRange &entities,
                                                  MBRange &new_entities,
                                                  MBRange *fill_entities) 
{
  MBRange tmp_range, *tmp_ptr_fill_entity;
  if (NULL != fill_entities) tmp_ptr_fill_entity = &tmp_range;
  else tmp_ptr_fill_entity = NULL;
  
  for (MBRange::iterator rit = entities.begin(); rit != entities.end(); rit++) {
    MBEntityHandle new_entity;
    if (NULL != tmp_ptr_fill_entity) tmp_ptr_fill_entity->clear();

    MBEntityHandle this_ent = *rit;
    MBErrorCode result = split_entities_manifold(&this_ent, 1, &new_entity, 
                                                 tmp_ptr_fill_entity);
    if (MB_SUCCESS != result) return result;

    new_entities.insert(new_entity);
    if (NULL != fill_entities) fill_entities->merge(*tmp_ptr_fill_entity);
  }
  
  return MB_SUCCESS;
}


MBErrorCode MeshTopoUtil::split_entities_manifold(MBEntityHandle *entities,
                                                  const int num_entities,
                                                  MBEntityHandle *new_entities,
                                                  MBRange *fill_entities,
                                                  MBEntityHandle *gowith_ents)
{
    // split entities by duplicating them; splitting manifold means that there is at
    // most two higher-dimension entities bounded by a given entity; after split, the
    // new entity bounds one and the original entity bounds the other

#define ITERATE_RANGE(range, it) for (MBRange::iterator it = range.begin(); it != range.end(); it++)
#define GET_CONNECT_DECL(ent, connect, num_connect) \
  const MBEntityHandle *connect; int num_connect; \
  {MBErrorCode connect_result = mbImpl->get_connectivity(ent, connect, num_connect); \
   if (MB_SUCCESS != connect_result) return connect_result;}
#define GET_CONNECT(ent, connect, num_connect) \
  {MBErrorCode connect_result = mbImpl->get_connectivity(ent, connect, num_connect);\
   if (MB_SUCCESS != connect_result) return connect_result;}
#define TC if (MB_SUCCESS != tmp_result) {result = tmp_result; continue;}

  MBErrorCode result = MB_SUCCESS;
  for (int i = 0; i < num_entities; i++) {
    MBErrorCode tmp_result;
      
      // get original higher-dimensional bounding entities
    MBRange up_adjs[4];
      // can only do a split_manifold if there are at most 2 entities of each
      // higher dimension; otherwise it's a split non-manifold
    bool valid_up_adjs = true;
    for (int dim = 1; dim <= 3; dim++) {
      tmp_result = mbImpl->get_adjacencies(entities+i, 1, dim, false, up_adjs[dim]); TC;
      if (dim > MBCN::Dimension(TYPE_FROM_HANDLE(entities[i])) && up_adjs[dim].size() > 2) {
        valid_up_adjs = false;
        break;
      }
    }
    if (!valid_up_adjs) return MB_FAILURE;
      
      // ok to split; create the new entity, with connectivity of the original
    GET_CONNECT_DECL(entities[i], connect, num_connect);
    MBEntityHandle new_entity;
    result = mbImpl->create_element(mbImpl->type_from_handle(entities[i]), connect, num_connect, 
                                    new_entity); TC;
      
      // by definition, new entity and original will be equivalent; need to add explicit
      // adjs to distinguish them; don't need to check if there's already one there,
      // 'cuz add_adjacency does that for us
    for (int dim = 1; dim <= 3; dim++) {
      if (up_adjs[dim].empty() || dim == MBCN::Dimension(TYPE_FROM_HANDLE(entities[i]))) continue;

      if (dim < MBCN::Dimension(TYPE_FROM_HANDLE(entities[i]))) {
          // adjacencies from other entities to this one; if any of those are equivalent entities,
          // need to make explicit adjacency to new entity too
        for (MBRange::iterator rit = up_adjs[dim].begin(); rit != up_adjs[dim].end(); rit++) {
          if (equivalent_entities(*rit))
            result = mbImpl->add_adjacencies(*rit, &new_entity, 1, false);
        }
      }
      else {
        
          // get the two up-elements
        MBEntityHandle up_elem1 = *(up_adjs[dim].begin()),
            up_elem2 = (up_adjs[dim].size() > 1 ? *(up_adjs[dim].rbegin()) : 0);
        
          // if two, and a gowith entity was input, make sure the new entity goes with
          // that one
        if (gowith_ents && up_elem2 && 
            gowith_ents[i] != up_elem1 && gowith_ents[i] == up_elem2) {
          MBEntityHandle tmp_elem = up_elem1;
          up_elem1 = up_elem2;
          up_elem2 = tmp_elem;
        }
        
        tmp_result = mbImpl->remove_adjacencies(entities[i], &up_elem1, 1);
          // (ok if there's an error, that just means there wasn't an explicit adj)

        tmp_result = mbImpl->add_adjacencies(new_entity, &up_elem1, 1, false); TC;
        if (!up_elem2) continue;

          // add adj to other up_adj
        tmp_result = mbImpl->add_adjacencies(entities[i], &up_elem2, 1, false); TC;
      }
    }

      // if we're asked to build a next-higher-dimension object, do so
    MBEntityHandle fill_entity = 0;
    MBEntityHandle tmp_ents[2];
    if (NULL != fill_entities) {
        // how to do this depends on dimension
      switch (MBCN::Dimension(TYPE_FROM_HANDLE(entities[i]))) {
        case 0:
          tmp_ents[0] = entities[i];
          tmp_ents[1] = new_entity;
          tmp_result = mbImpl->create_element(MBEDGE, tmp_ents, 2, fill_entity); TC;
          break;
        case 1:
          tmp_result = mbImpl->create_element(MBPOLYGON, connect, 2, fill_entity); TC;
            // need to create explicit adj in this case
          tmp_result = mbImpl->add_adjacencies(entities[i], &fill_entity, 1, false); TC;
          tmp_result = mbImpl->add_adjacencies(new_entity, &fill_entity, 1, false); TC;
          break;
        case 2:
          tmp_ents[0] = entities[i];
          tmp_ents[1] = new_entity;
          tmp_result = mbImpl->create_element(MBPOLYHEDRON, tmp_ents, 2, fill_entity); TC;
          break;
      }
      if (0 == fill_entity) {
        result = MB_FAILURE;
        continue;
      }
      fill_entities->insert(fill_entity);
    }

    new_entities[i] = new_entity;
    
  } // end for over input entities

  return result;
}

MBErrorCode MeshTopoUtil::split_entity_nonmanifold(MBEntityHandle split_ent,
                                                   MBRange &old_adjs,
                                                   MBRange &new_adjs,
                                                   MBEntityHandle &new_entity) 
{
    // split an entity into two entities; new entity gets explicit adj to new_adjs,
    // old to old_adjs

    // make new entities and add adjacencies
    // create the new entity
  MBEntityType split_type = mbImpl->type_from_handle(split_ent);
  
  MBErrorCode result;
  if (MBVERTEX == split_type) {
    double coords[3];
    result = mbImpl->get_coords(&split_ent, 1, coords); RR;
    result = mbImpl->create_vertex(coords, new_entity); RR;
  }
  else {
    const MBEntityHandle *connect;
    int num_connect;
    result = mbImpl->get_connectivity(split_ent, connect, num_connect); RR;
    result = mbImpl->create_element(split_type, connect, num_connect, new_entity); RR;

      // remove any explicit adjacencies between new_adjs and split entity
    for (MBRange::iterator rit = new_adjs.begin(); rit != new_adjs.end(); rit++)
      mbImpl->remove_adjacencies(split_ent, &(*rit), 1);
  }
      
  if (MBVERTEX != split_type) {
        //  add adj's between new_adjs & new entity, old_adjs & split_entity
    for (MBRange::iterator rit = new_adjs.begin(); rit != new_adjs.end(); rit++)
      mbImpl->add_adjacencies(new_entity, &(*rit), 1, true);
    for (MBRange::iterator rit = old_adjs.begin(); rit != old_adjs.end(); rit++)
      mbImpl->add_adjacencies(split_ent, &(*rit), 1, true);
  }
  else if (split_ent != new_entity) {
      // in addition to explicit adjs, need to check if vertex is part of any
      // other entities, and check those entities against ents in old and new adjs
    MBRange other_adjs;
    for (int i = 1; i < 4; i++) {
      result = mbImpl->get_adjacencies(&split_ent, 1, i, false, other_adjs, 
                                       MBInterface::UNION); RR;
    }
    other_adjs = subtract( other_adjs, old_adjs);
    other_adjs = subtract( other_adjs, new_adjs);
    for (MBRange::iterator rit1 = other_adjs.begin(); rit1 != other_adjs.end(); rit1++) {
        // find an adjacent lower-dimensional entity in old_ or new_ adjs
      bool found = false;
      for (MBRange::iterator rit2 = old_adjs.begin(); rit2 != old_adjs.end(); rit2++) {
        if (mbImpl->dimension_from_handle(*rit1) != mbImpl->dimension_from_handle(*rit2) &&
            common_entity(*rit1, *rit2, mbImpl->dimension_from_handle(*rit1))) {
          found = true; old_adjs.insert(*rit1); break;
        }
      }
      if (found) continue;
      for (MBRange::iterator rit2 = new_adjs.begin(); rit2 != new_adjs.end(); rit2++) {
        if (mbImpl->dimension_from_handle(*rit1) != mbImpl->dimension_from_handle(*rit2) &&
            common_entity(*rit1, *rit2, mbImpl->dimension_from_handle(*rit1))) {
          found = true; new_adjs.insert(*rit1); break;
        }
      }
      if (!found) return MB_FAILURE;
    }
          
      // instead of adjs replace in connectivity
    std::vector<MBEntityHandle> connect;
    for (MBRange::iterator rit = new_adjs.begin(); rit != new_adjs.end(); rit++) {
      connect.clear();
      result = mbImpl->get_connectivity(&(*rit), 1, connect); RR;
      std::replace(connect.begin(), connect.end(), split_ent, new_entity);
      result = mbImpl->set_connectivity(*rit, &connect[0], connect.size()); RR;
    }
  }
  
  return result;

/*

Commented out for now, because I decided to do a different implementation
for the sake of brevity.  However, I still think this function is the right
way to do it, if I ever get the time.  Sigh.

    // split entity d, producing entity nd; generates various new entities,
    // see algorithm description in notes from 2/25/05
  const MBEntityHandle split_types = {MBEDGE, MBPOLYGON, MBPOLYHEDRON};
  MBErrorCode result = MB_SUCCESS;
  const int dim = MBCN::Dimension(TYPE_FROM_HANDLE(d));
  MeshTopoUtil mtu(this);

    // get all (d+2)-, (d+1)-cells connected to d
  MBRange dp2s, dp1s, dp1s_manif, dp2s_manif;
  result = get_adjacencies(&d, 1, dim+2, false, dp2s); RR;
  result = get_adjacencies(&d, 1, dim+1, false, dp1s); RR;

    // also get (d+1)-cells connected to d which are manifold
  get_manifold_dp1s(d, dp1s_manif);
  get_manifold_dp2s(d, dp2s_manif);
  
    // make new cell nd, then ndp1
  result = copy_entity(d, nd); RR;
  MBEntityHandle tmp_connect[] = {d, nd};
  MBEntityHandle ndp1;
  result = create_element(split_types[dim],
                          tmp_connect, 2, ndp1); RR;
  
    // modify (d+2)-cells, depending on what type they are
  ITERATE_RANGE(dp2s, dp2) {
      // first, get number of connected manifold (d+1)-entities
    MBRange tmp_range, tmp_range2(dp1s_manif);
    tmp_range.insert(*dp2);
    tmp_range.insert(d);
    tmp_result = get_adjacencies(tmp_range, 1, false, tmp_range2); TC;
    MBEntityHandle ndp2;
    
      // a. manif (d+1)-cells is zero
    if (tmp_range2.empty()) {
        // construct new (d+1)-cell
      MBEntityHandle ndp1a;
      MBEntityHandle tmp_result = create_element(split_types[dim],
                                                 tmp_connect, 2, ndp1a); TC;
        // now make new (d+2)-cell
      MBEntityHandle tmp_connect2[] = {ndp1, ndp1a};
      tmp_result = create_element(split_types[dim+1],
                                  tmp_connect2, 2, ndp2); TC;
        // need to add explicit adjacencies, since by definition ndp1, ndp1a will be equivalent
      tmp_result = add_adjacencies(ndp1a, &dp2, 1, false); TC;
      tmp_result = add_adjacencies(ndp1a, &ndp2, 1, false); TC;
      tmp_result = add_adjacencies(ndp1, &ndp2, 1, false); TC;

        // now insert nd into connectivity of dp2, right after d if dim < 1
      std::vector<MBEntityHandle> connect;
      tmp_result = get_connectivity(&dp2, 1, connect); TC;
      if (dim < 1) {
        std::vector<MBEntityHandle>::iterator vit = std::find(connect.begin(), connect.end(), d);
        if (vit == connect.end()) {
          result = MB_FAILURE;
          continue;
        }
        connect.insert(vit, nd);
      }
      else
        connect.push_back(nd);
      tmp_result = set_connectivity(dp2, connect); TC;

        // if dim < 1, need to add explicit adj from ndp2 to higher-dim ents, since it'll
        // be equiv to other dp2 entities
      if (dim < 1) {
        MBRange tmp_dp3s;
        tmp_result = get_adjacencies(&dp2, 1, dim+3, false, tmp_dp3s); TC;
        tmp_result = add_adjacencies(ndp2, tmp_dp3s, false); TC;
      }
    } // end if (tmp_range2.empty())

      // b. single manifold (d+1)-cell, which isn't adjacent to manifold (d+2)-cell
    else if (tmp_range2.size() == 1) {
        // b1. check validity, and skip if not valid

        // only change if not dp1-adjacent to manifold dp2cell; check that...
      MBRange tmp_adjs(dp2s_manif);
      tmp_result = get_adjacencies(&(*tmp_range2.begin()), 1, dim+2, false, tmp_adjs); TC;
      if (!tmp_adjs.empty()) continue;

      MBEntityHandle dp1 = *tmp_range2.begin();

        // b2. make new (d+1)- and (d+2)-cell next to dp2

        // get the (d+2)-cell on the other side of dp1
      tmp_result = get_adjacencies(&dp1, 1, dim+2, false, tmp_adjs); TC;
      MBEntityHandle odp2 = *tmp_adjs.begin();
      if (odp2 == dp2) odp2 = *tmp_adjs.rbegin();

        // get od, the d-cell on dp1_manif which isn't d
      tmp_result = get_adjacencies(&dp1_manif, 1, dim, false, tmp_adjs); TC;
      tmp_adjs.erase(d);
      if (tmp_adjs.size() != 1) {
        result = MB_FAILURE;
        continue;
      }
      MBEntityHandle od = *tmp_adjs.begin();

        // make a new (d+1)-cell from od and nd
      tmp_adjs.insert(nd);
      tmp_result = create_element(split_types[1], tmp_adjs, ndp1a); TC;
      
        // construct new (d+2)-cell from dp1, ndp1, ndp1a
      tmp_adjs.clear();
      tmp_adjs.insert(dp1); tmp_adjs.insert(ndp1); tmp_adjs.insert(ndp1a);
      tmp_result = create_element(split_types[2], tmp_adjs, ndp2); TC;

        // b3. replace d, dp1 in connect/adjs of odp2
      std::vector<MBEntityHandle> connect;
      tmp_result = get_connectivity(&odp2, 1, connect); TC;
      if (dim == 0) {
        *(std::find(connect.begin(), connect.end(), d)) = nd;
        remove_adjacency(dp1, odp2);
        

      
        // if dp1 was explicitly adj to odp2, remove it
      remove_adjacency

...

*/
}

    //! return whether entity is equivalent to any other of same type and same vertices;
    //! if equivalent entity is found, it's returned in equiv_ents and return value is true,
    //! false otherwise.
bool MeshTopoUtil::equivalent_entities(const MBEntityHandle entity,
                                       MBRange *equiv_ents) 
{
  const MBEntityHandle *connect;
  int num_connect;
  MBErrorCode result = mbImpl->get_connectivity(entity, connect, num_connect);
  if (MB_SUCCESS != result) return false;

  MBRange dum;
  result = mbImpl->get_adjacencies(connect, num_connect, 
                                   mbImpl->dimension_from_handle(entity),
                                   false, dum);
  dum.erase(entity);

  if (NULL != equiv_ents) {
    equiv_ents->swap(dum);
  }
  
  if (!dum.empty()) return true;
  else return false;
}

  
