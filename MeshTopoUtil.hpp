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

/*!
 *  \class   MeshTopoUtil
 *  \authors Tim Tautges
 *  \date    2/04
 *  \brief   MeshTopoUtil contains general mesh utility functions 
 *          
 */ 

#ifndef MESH_TOPO_UTIL_HPP
#define MESH_TOPO_UTIL_HPP

#include "MBForward.hpp"

class MeshTopoUtil
{
public:
  MeshTopoUtil(MBInterface *impl) : mbImpl(impl) {}
  
  ~MeshTopoUtil() {}

    //! generate all the AEntities bounding the vertices
  MBErrorCode construct_aentities(const MBRange &vertices);

    //! given an entity, get its average position (avg vertex locations)
  MBErrorCode get_average_position(MBRange &entities,
                                   double *avg_position);

    //! given an entity, get its average position (avg vertex locations)
  MBErrorCode get_average_position(const MBEntityHandle entity,
                                   double *avg_position);

    //! given a set of entities, get their average position (avg vertex locations)
  MBErrorCode get_average_position(const MBEntityHandle *entities,
                                   const int num_entities,
                                   double *avg_position);

    //! get (target_dim)-dimensional manifold entities connected to star_entity; that is,
    //! the entities with <= 1 connected (target_dim+2)-dimensional adjacent entities;
    //! for target_dim=3, just return all of them
    //! just insert into the list, w/o clearing manifold list first
  MBErrorCode get_manifold(const MBEntityHandle star_entity,
                           const int target_dim,
                           MBRange &manifold);
  
  //! given an entity, find the entities of next higher dimension around
  //! that entity, ordered by connection through next higher dimension entities; 
  //! if any of the star entities is in only entity of next higher dimension, 
  //! on_boundary is returned true
  MBErrorCode star_entities(const MBEntityHandle star_center,
                            std::vector<MBEntityHandle> &star_entities,
                            bool &bdy_entity,
                            const MBEntityHandle starting_star_entity = 0,
                            std::vector<MBEntityHandle> *star_entities_dp1 = NULL,
                            MBRange *star_entities_candidates_dp1 = NULL);

    //! Get a series of (d+1)-dimensional stars around a d-dimensional entity, such that
    //! each star is on a (d+2)-manifold containing the d-dimensional entity; each star
    //! is either open or closed, and also defines a (d+2)-star whose entities are bounded by
    //! (d+1)-entities on the star and on the (d+2)-manifold
  MBErrorCode star_entities_nonmanifold(const MBEntityHandle star_entity,
                                        std::vector<std::vector<MBEntityHandle> > &stars,
                                        std::vector<bool> *bdy_flags = NULL,
                                        std::vector<std::vector<MBEntityHandle> > *dp2_stars = NULL);

    //! given a star_center, a last_entity (whose dimension should be 1 greater than center)
    //! and last_dp1 (dimension 2 higher than center), returns the next star entity across
    //! last_dp1, and the next dp1 entity sharing next_entity; if star_candidates is non-empty,
    //! star must come from those
  MBErrorCode star_next_entity(const MBEntityHandle star_center,
                               const MBEntityHandle last_entity,
                               const MBEntityHandle last_dp1,
                               MBRange *star_candidates_dp1,
                               MBEntityHandle &next_entity,
                               MBEntityHandle &next_dp1);
  
    //! get "bridge" or "2nd order" adjacencies, going through dimension bridge_dim
  MBErrorCode get_bridge_adjacencies(const MBEntityHandle from_entity,
                                     const int bridge_dim,
                                     const int to_dim,
                                     MBRange &to_adjs);

    //! return a common entity of the specified dimension, or 0 if there isn't one
  MBEntityHandle common_entity(const MBEntityHandle ent1,
                               const MBEntityHandle ent2,
                               const int dim);
  
    //! split entity which is non-manifold, that is, which has > 2 connected entities
    //! of next higher dimension; assumes that there are >= 2 connected regions of
    //! (d+2)-dimensional entities; a new d-entity is created for each region after the
    //! first, and it's made explicitly-adjacent to the region to which it corresponds
  MBErrorCode split_entity_nonmanifold(MBEntityHandle split_ent,
                                       MBRange &new_ents);
  
    //! split entities that are manifold (shared by two or less entities of each higher dimension),
    //! optionally creating an entity of next higher dimension to fill the gap
    /**
       \param entities The entities to be split
       \param new_entities New entities, in order of correspondence to that of entities
       \param fill_entities If non-NULL, create an entity of next higher dimension to fill the gap,
                       passing it back in *fill_entities
    */
  MBErrorCode split_entities_manifold(MBRange &entities,
                                      MBRange &new_entities,
                                      MBRange *fill_entities);
  
    //! split entities that are manifold (shared by two or less entities of each higher dimension),
    //! optionally creating an entity of next higher dimension to fill the gap
    /**
       \param entities The entities to be split
       \param new_entities New entities, in order of correspondence to that of entities
       \param fill_entities If non-NULL, create an entity of next higher dimension to fill the gap,
                       passing it back in *fill_entities
    */
  MBErrorCode split_entities_manifold(MBEntityHandle *entities,
                                      const int num_entities,
                                      MBEntityHandle *new_entities,
                                      MBRange *fill_entities);

    //! return whether entity is equivalent to any other of same type and same vertices;
    //! if equivalent entity is found, it's returned in equiv_ents and return value is true,
    //! false otherwise.
  bool equivalent_entities(const MBEntityHandle entity,
                           MBRange *equiv_ents = NULL);
  
                                  
private:
  MBInterface *mbImpl;
  
};


#endif

