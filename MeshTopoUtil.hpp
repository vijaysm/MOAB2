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

#include "MBInterface.hpp"

class MeshTopoUtil
{
public:
  MeshTopoUtil(MBInterface *impl) : mbImpl(impl) {}
  
  ~MeshTopoUtil() {}

    //! generate all the AEntities bounding the vertices
  MBErrorCode construct_aentities(const MBRange &vertices);

    //! given an entity, get its average position (avg vertex locations)
  MBErrorCode get_average_position(const MBEntityHandle entity,
                                   double *avg_position);

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
  
private:
  MBInterface *mbImpl;
  
};


#endif

