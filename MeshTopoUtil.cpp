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

    //! given a mesh edge, find the ordered faces around the edge; if any
    //! of the faces is in only one region, on_boundary is returned true
MBErrorCode MeshTopoUtil::star_faces(const MBEntityHandle edge,
                                     std::vector<MBEntityHandle> &rad_faces,
                                     bool &bdy_edge,
                                     std::vector<MBEntityHandle> *rad_regions)
{
  MBRange all_faces;
  MBRange in_range, out_range;
  in_range.insert(edge);
  MBErrorCode result = mbImpl->get_adjacencies(in_range, 2, false, all_faces);
  if (MB_SUCCESS != result) return result;

    // if any of the faces have a single connected region, choose that,
    // otherwise choose any
  MBRange::iterator rit;
  MBEntityHandle last_face = 0, first_face = 0;
  for (rit = all_faces.begin(); rit != all_faces.end(); rit++) {
    in_range.clear();
    in_range.insert(*rit);
    out_range.clear();
    result = mbImpl->get_adjacencies(in_range, 3, false, out_range);
    if (MB_SUCCESS != result) return result;
    if (out_range.size() == 1) {
        // if we have a single-region face, take it off the all_faces list, since 
        // we'll never get back to it going around the edge
      all_faces.erase(*rit);
      last_face = *rit;
      first_face = last_face;
      rad_faces.push_back(first_face);
      break;
    }
  }

    // if no single-region faces, just pick the last one; don't take it off
    // the list, though, so that we get back to it
  if (0 == last_face)
    last_face = *all_faces.rbegin();
  
  MBRange regions;
  MBEntityHandle last_region;

  while (!all_faces.empty()) {
      // during each iteration:
      // - start with a last_face & edge
      // - find a region that's not on the list, and the other face in the region
      //   sharing that edge
      // - remove that face from all_faces & put on face list & assign to last_face
      // - add the region to the region list
      // - proceed to next iteration

      // get 3d elements common to face and edge
    in_range.clear();
    in_range.insert(edge);
    in_range.insert(last_face);
    out_range.clear();
    result = mbImpl->get_adjacencies(in_range, 3, false, out_range, 
                                     MBInterface::INTERSECT);
    if (MB_SUCCESS != result) return result;

      // find one which hasn't been treated yet
    last_region = 0;
    for (rit = out_range.begin(); rit != out_range.end(); rit++) {
      if (regions.find(*rit) == regions.end()) {
        last_region = *rit;
        break;
      }
    }

      // if we got here and we didn't find an untreated region
    if (0 == last_region) return MB_FAILURE;
    
      // get the other face sharing the edge
    in_range.clear(); out_range.clear();
    in_range.insert(edge);
    in_range.insert(last_region);
    result = mbImpl->get_adjacencies(in_range, 2, false, out_range);
    if (MB_SUCCESS != result) return result;
    else if (out_range.size() != 2) return MB_FAILURE;

    rit = out_range.begin();
    if (last_face != *rit)
      last_face = *rit;
    else if (last_face != *(++rit))
      last_face = *rit;
    else return MB_FAILURE;

      // remove the face from all_faces and add region to regions
    all_faces.erase(last_face);
    regions.insert(last_region);

      // get dual vertex for the region & put on list
    rad_faces.push_back(last_face);
    if (NULL != rad_regions) rad_regions->push_back(last_region);
  }

    // if it's a closed loop, we got all the vertices; if not, we need 2 more;
    // closed is indicated by first_face being zero
  if (0 != first_face) {
    bdy_edge = true;
  }
  else {
    bdy_edge = false;
  }
  
  return MB_SUCCESS;
}

