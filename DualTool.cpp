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

#include "DualTool.hpp"
#include "MeshTopoUtil.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBSkinner.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <assert.h>

#define RR if (MB_SUCCESS != result) return result

bool debug = false;

  //! tag name for dual surfaces
const char *DualTool::DUAL_SURFACE_TAG_NAME = "DUAL_SURFACE";

  //! tag name for dual curves
const char *DualTool::DUAL_CURVE_TAG_NAME = "DUAL_CURVE";

  //! tag name for dual cells
const char *DualTool::IS_DUAL_CELL_TAG_NAME = "IS_DUAL_CELL";

  //! tag name for dual entities
const char *DualTool::DUAL_ENTITY_TAG_NAME = "DUAL_ENTITY";

  //! tag name for extra dual entities
const char *DualTool::EXTRA_DUAL_ENTITY_TAG_NAME = "__EXTRA_DUAL_ENTITY";

  //! tag name for graphics point
const char *DualTool::DUAL_GRAPHICS_POINT_TAG_NAME = "__DUAL_GRAPHICS_POINT";

//const int DualTool::GP_SIZE = 20;

DualTool::DualTool(MBInterface *impl) 
    : mbImpl(impl)
{
  int dummy = 0;

  MBErrorCode result = mbImpl->tag_create(DUAL_SURFACE_TAG_NAME, sizeof(MBEntityHandle), 
                                          MB_TAG_SPARSE, 
                                          dualSurfaceTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
  result = mbImpl->tag_create(DUAL_CURVE_TAG_NAME, sizeof(MBEntityHandle), MB_TAG_SPARSE, 
                              dualCurveTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
  result = mbImpl->tag_create(IS_DUAL_CELL_TAG_NAME, 1, MB_TAG_BIT, 
                              isDualCellTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);

  result = mbImpl->tag_create(DUAL_ENTITY_TAG_NAME, sizeof(MBEntityHandle), MB_TAG_DENSE, 
                              dualEntityTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);

  result = mbImpl->tag_create(EXTRA_DUAL_ENTITY_TAG_NAME, sizeof(MBEntityHandle), MB_TAG_SPARSE, 
                              extraDualEntityTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
  static const char dum_name[CATEGORY_TAG_NAME_LENGTH] = "\0";
  result = mbImpl->tag_create(CATEGORY_TAG_NAME, CATEGORY_TAG_NAME_LENGTH, MB_TAG_SPARSE, 
                              categoryTag, dum_name);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
  DualTool::GraphicsPoint dum_pt(0.0, 0.0, 0.0, -1);
  result = mbImpl->tag_create(DUAL_GRAPHICS_POINT_TAG_NAME, 
                              sizeof(DualTool::GraphicsPoint), MB_TAG_DENSE, 
                              dualGraphicsPointTag, &dum_pt);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
}

DualTool::~DualTool() 
{}

  //! construct the dual entities for the entire mesh
MBErrorCode DualTool::construct_dual() 
{
    // allocate a dual entity for each primal entity in the mesh, starting
    // with highest dimension and working downward; do each dimension in a separate code
    // block, since they're all handled slightly differently
  
  MBRange all_regions, all_faces, all_edges, all_vertices;

    // first, construct all the aentities, since they're currently needed to 
    // compute the dual
  MBErrorCode result = mbImpl->get_entities_by_dimension(0, 0, all_vertices);
  if (MB_SUCCESS != result) return result;

  result = MeshTopoUtil(mbImpl).construct_aentities(all_vertices);
  if (MB_SUCCESS != result) return result;

    // get all edges, faces and regions now, so we don't need to filter out dual
    // entities later
  
  result = mbImpl->get_entities_by_dimension(0, 1, all_edges);
  if (MB_SUCCESS != result) return result;
  result = mbImpl->get_entities_by_dimension(0, 2, all_faces);
  if (MB_SUCCESS != result) return result;
  result = mbImpl->get_entities_by_dimension(0, 3, all_regions);
  if (MB_SUCCESS != result) return result;

  MBRange dual_verts;
  result = construct_dual_vertices(all_regions, dual_verts);
  if (MB_SUCCESS != result || dual_verts.size() != all_regions.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_verts.size() << " dual vertices." << std::endl;

    // don't really need dual edges, but construct 'em anyway
  MBRange dual_edges;
  result = construct_dual_edges(all_faces, dual_edges);
  if (MB_SUCCESS != result || dual_edges.size() != all_faces.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_edges.size() << " dual edges." << std::endl;

    // construct dual faces
  MBRange dual_faces;
  result = construct_dual_faces(all_edges, dual_faces);
  if (MB_SUCCESS != result || dual_faces.size() != all_edges.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_faces.size() << " dual faces." << std::endl;

    // construct dual cells
  MBRange dual_cells;
  result = construct_dual_cells(all_vertices, dual_cells);
  if (MB_SUCCESS != result || dual_cells.size() != all_vertices.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_cells.size() << " dual cells." << std::endl;

  return MB_SUCCESS;
}


MBErrorCode DualTool::construct_dual_vertices(const MBRange &all_regions,
                                              MBRange &dual_ents) 
{
  if (all_regions.empty()) return MB_SUCCESS;
  
    // make sure they're all regions
  assert(3 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_regions.begin())) &&
         3 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_regions.rbegin())));
  
  MBRange::const_iterator rit;
  MBEntityHandle dual_ent;
  MBErrorCode tmp_result = MB_SUCCESS;
  MBErrorCode result = MB_SUCCESS;
  
  for (rit = all_regions.begin(); rit != all_regions.end(); rit++) {
    if (tmp_result != MB_SUCCESS) result = tmp_result;
    
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS == tmp_result && 0 != dual_ent) {
      dual_ents.insert(dual_ent);
      continue;
    }
    else if (MB_SUCCESS != tmp_result) continue;

    tmp_result = construct_dual_vertex(*rit, dual_ent, false, true);
    if (MB_SUCCESS != tmp_result) continue;

      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
  }

  return result;
}

MBErrorCode DualTool::construct_dual_vertex(MBEntityHandle entity, 
                                            MBEntityHandle &dual_ent, 
                                            const bool extra,
                                            const bool add_graphics_pt)
{
    // no dual entity; construct one; first need the avg coordinates
  unsigned int is_dual = 0x1;
  double avg_pos[3];
  MBErrorCode result = MeshTopoUtil(mbImpl).get_average_position(entity, avg_pos);
  if (MB_SUCCESS != result) return result;
    
    // now construct the new dual entity
  result = mbImpl->create_vertex(avg_pos, dual_ent);
  if (MB_SUCCESS != result) return result;
    
    // tag it indicating it's a dual entity
  result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
  if (MB_SUCCESS != result) return result;
    
    // tag the primal entity with its dual entity and vica versa
  if (extra) 
    result = mbImpl->tag_set_data(extraDualEntity_tag(), &(entity), 1, &dual_ent);
  else
    result = mbImpl->tag_set_data(dualEntity_tag(), &(entity), 1, &dual_ent);
  if (MB_SUCCESS != result) return result;
    
  result = mbImpl->tag_set_data(dualEntity_tag(), &dual_ent, 1, &(entity));
  if (MB_SUCCESS != result) return result;

  if (add_graphics_pt)
        // put a graphics point on that vertex too
    result = add_graphics_point(dual_ent, avg_pos);

  return result;
}

MBErrorCode DualTool::add_graphics_point(MBEntityHandle entity,
                                         double *avg_pos) 
{
    // add a graphics pt, placed at the same position as the vertex
  double my_pos[3];
  MBErrorCode result;
  
  if (NULL == avg_pos) {
    result = MeshTopoUtil(mbImpl).get_average_position(entity, my_pos);
    if (MB_SUCCESS != result) return result;
  }
  else
    for (int i = 0; i < 3; i++) my_pos[i] = avg_pos[i];
  
  DualTool::GraphicsPoint dum_pt(my_pos, -1);
  result = mbImpl->tag_set_data(dualGraphicsPoint_tag(), &entity, 1, 
                                &dum_pt);
  return result;
}

MBErrorCode DualTool::construct_dual_edges(const MBRange &all_faces,
                                           MBRange &dual_ents) 
{
  if (all_faces.empty()) return MB_SUCCESS;
  
    // make sure they're all faces
  assert(2 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_faces.begin())) &&
         2 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_faces.rbegin())));
  
  MBRange::const_iterator rit;
  MBEntityHandle dual_ent;
  unsigned int is_dual = 0x1;
  MBErrorCode tmp_result = MB_SUCCESS;
  MBErrorCode result = MB_SUCCESS;
  
  for (rit = all_faces.begin(); rit != all_faces.end(); rit++) {
    if (tmp_result != MB_SUCCESS) result = tmp_result;
    
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS == tmp_result && 0 != dual_ent) {
      dual_ents.insert(dual_ent);
      continue;
    }
    
      // no dual entity; construct one; get the bounding regions
    std::vector<MBEntityHandle> in_ents, out_ents;
    tmp_result = mbImpl->get_adjacencies(&(*rit), 1, 3, false, out_ents);
    if (MB_SUCCESS != tmp_result) continue;

      // get the dual vertices
    std::vector<MBEntityHandle> dual_verts(out_ents.size());
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &out_ents[0], out_ents.size(), 
                                      &dual_verts[0]);
    if (MB_SUCCESS != tmp_result) continue;
    assert(dual_verts.size() <= 2);
    
    double avg_pos[3];      
    bool bdy_face = (dual_verts.size() == 1 ? true : false);
    if (bdy_face) {
        // boundary face - make a dual vertex at the face center and put in list
      tmp_result = construct_dual_vertex(*rit, dual_ent, true, true);
      
        // put it on vertex list
      dual_verts.push_back(dual_ent);
    }
    
    assert(dual_verts.size() == 2);
    
      // now create the dual edge
    tmp_result = mbImpl->create_element(MBEDGE, &dual_verts[0], 2, dual_ent);
    if (MB_SUCCESS != tmp_result || 0 == dual_ent) continue;
    
      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
      // tag the primal entity with its dual entity and vica versa
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS != tmp_result) continue;

    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &dual_ent, 1, &(*rit));
    if (MB_SUCCESS != tmp_result) continue;

      // tag the edge indicating it's a dual entity
    tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
    if (MB_SUCCESS != tmp_result) continue;

      // add a graphics point to the edge; position depends on whether it's a
      // bdy face (mid-pt of dual edge) or not (mid-pt of primal face)
    if (bdy_face)
      tmp_result = add_graphics_point(dual_ent);
    else {
        // get the face's position
      tmp_result = MeshTopoUtil(mbImpl).get_average_position(*rit, avg_pos);
      if (MB_SUCCESS != tmp_result) continue;
      tmp_result = add_graphics_point(dual_ent, avg_pos);
    }
    if (MB_SUCCESS != tmp_result) continue;
  }

  return result;
}

MBErrorCode DualTool::construct_dual_faces(const MBRange &all_edges,
                                           MBRange &dual_ents) 
{
  if (all_edges.empty()) return MB_SUCCESS;
  
    // make sure they're all edges
  assert(1 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_edges.begin())) &&
         1 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_edges.rbegin())));
  
  MBRange::const_iterator rit;
  MBEntityHandle dual_ent;
  unsigned int is_dual = 0x1;
  MBErrorCode tmp_result = MB_SUCCESS;
  MBErrorCode result = MB_SUCCESS;
  
  for (rit = all_edges.begin(); rit != all_edges.end(); rit++) {
    if (tmp_result != MB_SUCCESS) 
      result = tmp_result;
    
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS == tmp_result && 0 != dual_ent) {
      dual_ents.insert(dual_ent);
      continue;
    }
    
      // no dual entity; construct one; get the dual vertices bounding the edge in radial order
    std::vector<MBEntityHandle> rad_dverts;
    bool bdy_edge;
    tmp_result = get_radial_dverts(*rit, rad_dverts, bdy_edge);
    if (MB_SUCCESS != tmp_result) 
      continue;
    
      // create the dual face
    tmp_result = mbImpl->create_element(MBPOLYGON, &rad_dverts[0], rad_dverts.size(), dual_ent);
    if (MB_SUCCESS != tmp_result || 0 == dual_ent) 
      continue;
    
      // tag it indicating it's a dual entity
    tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
    if (MB_SUCCESS != tmp_result) 
      continue;

      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
      // tag the primal entity with its dual entity and vica versa
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS != tmp_result) continue;
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &dual_ent, 1, &(*rit));
    if (MB_SUCCESS != tmp_result) continue;

      // add a graphics point to the cell; position depends on whether it's a
      // bdy cell (mid-pt of cell's vertices) or not (mid-pt of primal edge)
    double avg_pos[3];
    tmp_result = MeshTopoUtil(mbImpl).get_average_position(*rit, avg_pos);
    if (MB_SUCCESS != tmp_result) continue;
    if (bdy_edge) {
      
      tmp_result = add_graphics_point(dual_ent);
      if (MB_SUCCESS != tmp_result) continue;

        // also, if it's a bdy edge, add a new dual edge betw last 2 verts
      MBEntityHandle new_edge;
      tmp_result = mbImpl->create_element(MBEDGE, &rad_dverts[rad_dverts.size()-2], 
                                          2, new_edge);
      if (MB_SUCCESS != tmp_result) continue;

        // tag it indicating it's a dual entity
      tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &new_edge, 1, &is_dual);
      if (MB_SUCCESS != tmp_result) continue;

        // add a graphics pt, position is center of primal edge
      tmp_result = add_graphics_point(new_edge, avg_pos);
      if (MB_SUCCESS != tmp_result) continue;
    }
    
    else {
        // if inside, point goes on the 2cell, at primal edge mid-pt
      tmp_result = add_graphics_point(dual_ent, avg_pos);
    }
    if (MB_SUCCESS != tmp_result) continue;
  }

  return result;
}

MBErrorCode DualTool::construct_dual_cells(const MBRange &all_verts,
                                           MBRange &dual_ents) 
{
  if (all_verts.empty()) return MB_SUCCESS;
  
    // make sure they're all edges
  assert(0 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_verts.begin())) &&
         0 == MBCN::Dimension(TYPE_FROM_HANDLE(*all_verts.rbegin())));
  
  MBRange::const_iterator rit;
  MBEntityHandle dual_ent;
  unsigned int is_dual = 0x1;
  MBErrorCode tmp_result = MB_SUCCESS;
  MBErrorCode result = MB_SUCCESS;
  std::vector<MBEntityHandle> edges, dfaces;
  
  for (rit = all_verts.begin(); rit != all_verts.end(); rit++) {
    if (tmp_result != MB_SUCCESS) result = tmp_result;
    
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS == tmp_result && 0 != dual_ent) {
      dual_ents.insert(dual_ent);
      continue;
    }
    
      // no dual entity; construct one; get the edges bounding the vertex
    edges.clear();
    dfaces.clear();
    tmp_result = mbImpl->get_adjacencies(&(*rit), 1, 1, false, edges);
    if (MB_SUCCESS != tmp_result) continue;

      // get the dual faces corresponding to the edges
    dfaces.resize(edges.size());
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &edges[0], edges.size(), &dfaces[0]);
    if (MB_SUCCESS != tmp_result) continue;
    
      // create the dual cell from those faces
    tmp_result = mbImpl->create_element(MBPOLYHEDRON, &dfaces[0], dfaces.size(), dual_ent);
    if (MB_SUCCESS != tmp_result || 0 == dual_ent) continue;
    
      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
      // tag it indicating it's a dual entity
    tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
    if (MB_SUCCESS != tmp_result) continue;

      // tag the primal entity with its dual entity and vica versa
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS != tmp_result) continue;
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &dual_ent, 1, &(*rit));
    if (MB_SUCCESS != tmp_result) continue;
  }

  return result;
}

  //! given an edge handle, return a list of dual vertices in radial order 
  //! around the edge
MBErrorCode DualTool::get_radial_dverts(const MBEntityHandle edge,
                                        std::vector<MBEntityHandle> &rad_dverts,
                                        bool &bdy_edge) 
{
  rad_dverts.clear();
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
      break;
    }
  }

    // if no single-region faces, just pick the last one; don't take it off
    // the list, though, so that we get back to it
  if (0 == last_face)
    last_face = *all_faces.rbegin();
  
  MBRange regions;
  MBEntityHandle last_region;
  MBEntityHandle dvert;

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
    result = mbImpl->tag_get_data(dualEntity_tag(), &last_region, 1, &dvert);
    if (MB_SUCCESS != result) return result;
    assert(0 != dvert);
    rad_dverts.push_back(dvert);
  }

    // if it's a closed loop, we got all the vertices; if not, we need 2 more;
    // closed is indicated by first_face being zero
  if (0 != first_face) {
    bdy_edge = true;
    
      // not closed - get the last and first vertices
      // last is on dual of last_face not on the end of the list
      // get dual of last face
    MBEntityHandle dedge;
    result = mbImpl->tag_get_data(dualEntity_tag(), &last_face, 1, &dedge);
    if (MB_SUCCESS != result) return result;
    
      // get connectivity of that edge
    const MBEntityHandle *connect;
    int num_connect;
    result = mbImpl->get_connectivity(dedge, connect, num_connect);
    if (MB_SUCCESS != result) return result;
    
      // we want the one that's not already on the list; reuse last_face
    last_face = (connect[0] == rad_dverts.back() ? connect[1] : connect[0]);
    rad_dverts.push_back(last_face);

      // do the same for first_face
    result = mbImpl->tag_get_data(dualEntity_tag(), &first_face, 1, &dedge);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->get_connectivity(dedge, connect, num_connect);
    if (MB_SUCCESS != result) return result;
    first_face = (connect[0] == rad_dverts[0] ? connect[1] : connect[0]);
    rad_dverts.push_back(first_face);
/*
      // make a dual edge for this (reuse connect array and first_face)
    MBEntityHandle tmp_conn[3];
    tmp_conn[0] = first_face;
    tmp_conn[1] = last_face;
    result = mbImpl->create_element(MBEDGE, connect, 2, tmp_conn[2]);

      // tag it indicating it's a dual entity
    unsigned int is_dual = 0x1;
    result = mbImpl->tag_set_data(isDualCell_tag(), &tmp_conn[2], 1, &is_dual);
    if (MB_SUCCESS != result) return result;

      // add a graphics point to it too
    result = add_graphics_point(tmp_conn[2]);
    if (MB_SUCCESS != result) return result;
*/
  }
  else {
    bdy_edge = false;
  }
  
  return MB_SUCCESS;
}

  //! construct the dual entities for a hex mesh, including dual surfaces & curves
MBErrorCode DualTool::construct_hex_dual() 
{
    // really quite simple: 

    // construct the dual...
  MBErrorCode result = construct_dual();
  if (MB_SUCCESS != result) return result;
  
    // now traverse to build 1d and 2d hyperplanes
  result = construct_dual_hyperplanes(1);
  if (MB_SUCCESS != result) return result;

  result = construct_dual_hyperplanes(2);
  if (MB_SUCCESS != result) return result;
  
  result = construct_hp_parent_child();
  if (MB_SUCCESS != result) return result;

    // see?  simple, just like I said
  return MB_SUCCESS;
}

  //! get the cells of the dual
MBErrorCode DualTool::get_dual_entities(const int dim, MBRange &dual_ents) 
{
  if (0 == isDualCell_tag()) return MB_SUCCESS;
  if (0 > dim || 3 < dim) return MB_INDEX_OUT_OF_RANGE;

  unsigned int dum = 0x1;
  const void *dum_ptr = &dum;
  static MBEntityType dual_type[] = {MBVERTEX, MBEDGE, MBPOLYGON, MBPOLYHEDRON};
  
  return
    mbImpl->get_entities_by_type_and_tag(0, dual_type[dim], &isDualCellTag, &dum_ptr, 1,
                                         dual_ents);
}

  //! get the faces of the dual
MBErrorCode DualTool::get_dual_entities(const int dim, 
                                        std::vector<MBEntityHandle> &dual_ents) 
{
  MBRange tmp_range;
  MBErrorCode result = get_dual_entities(dim, tmp_range);
  if (MB_SUCCESS != result)
    return result;

    // dual_ents.insert(dual_ents.end(), tmp_range.begin(), tmp_range.end());
  dual_ents.reserve(dual_ents.size() + tmp_range.size());
  for (MBRange::const_iterator it = tmp_range.begin();
       it != tmp_range.end();
       ++it)
  {
    dual_ents.push_back(*it);
  }
  return MB_SUCCESS;
}

MBErrorCode DualTool::get_dual_hyperplanes(const MBInterface *impl, const int dim, 
                                           MBRange &dual_ents) 
{
  if (dim != 1 && dim != 2) return MB_INDEX_OUT_OF_RANGE;

  MBTag dual_tag;
  MBErrorCode result;
  
  if (dim == 1)
    result = impl->tag_get_handle(DUAL_CURVE_TAG_NAME, dual_tag );
  else 
    result = impl->tag_get_handle(DUAL_SURFACE_TAG_NAME, dual_tag );

  if (MB_SUCCESS == result)
    result = impl->get_entities_by_type_and_tag(0, MBENTITYSET, &dual_tag, NULL, 1, dual_ents,
                                                MBInterface::UNION);
      
  return result;
}

MBErrorCode DualTool::construct_dual_hyperplanes(const int dim) 
{
    // this function traverses dual faces of input dimension, constructing
    // dual hyperplanes of them in sets as it goes

    // check various inputs
  int num_quads, num_hexes;
  if (
      // this function only makes sense for dim == 1 or dim == 2
    (dim != 1 && dim != 2) ||
      // should either be quads or hexes around
    mbImpl->get_number_entities_by_type(0, MBQUAD, num_quads) != MB_SUCCESS ||
    mbImpl->get_number_entities_by_type(0, MBQUAD, num_hexes) != MB_SUCCESS ||
      // if we're asking for 1d dual ents, should be quads around
    (num_quads == 0 && dim == 1) ||
      // if we're asking for 2d dual ents, should be hexes around
    (num_hexes == 0 && dim == 2))
    return MB_FAILURE;
  
    // get tag name for this dimension hyperplane
  MBTag gid_tag, mark_tag;
  int dum = -1;
  MBErrorCode result = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, 4, 
                                                        MB_TAG_DENSE, gid_tag, &dum);

  MBTag hp_tag = (1 == dim ? dualCurve_tag() : dualSurface_tag());
  
  unsigned short mark_val = 0x0;
  result = mbImpl->tag_get_handle("__hyperplane_mark", mark_tag);
  if (MB_SUCCESS != result) {
    dum = 0x0;
    result = mbImpl->tag_create("__hyperplane_mark", 1, 
                                MB_TAG_BIT, mark_tag, &mark_val);
    if (MB_SUCCESS != result) return result;
  }
  mark_val = 0x1;
  
    // two stacks: one completely untreated entities, and the other untreated 
    // entities on the current dual hyperplane
  std::vector<MBEntityHandle> tot_untreated, hp_untreated;
  
    // put dual entities of this dimension on the untreated list
  result = get_dual_entities(dim, tot_untreated);
  if (MB_SUCCESS != result) return result;
  
    // main part of traversal loop
  MBEntityHandle this_ent;
  MBEntityHandle this_hp;
  int hp_val;
  int hp_ids = 1;
  std::vector<MBEntityHandle> parents;

  while (!tot_untreated.empty()) {
    if (debug && dim == 2 /*(tot_untreated.size()%report == 0)*/)
      std::cout << "Untreated list size " << tot_untreated.size() << "." << std::endl;
      
    this_ent = tot_untreated.back(); tot_untreated.pop_back();
    result = mbImpl->tag_get_data(hp_tag, &this_ent, 1, &this_hp);
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return result;

      // test for this entity having a hyperplane assignment already
    else if (this_hp != 0)
      continue;

    if (1 == dim && check_1d_loop_edge(this_ent)) continue;
    
      // ok, doesn't have one; make a new hyperplane
    hp_val = hp_ids++;
    result = mbImpl->create_meshset(MESHSET_SET, this_hp);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->tag_set_data(gid_tag, &this_hp, 1, &hp_val);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->tag_set_data(hp_tag, &this_hp, 1, &this_hp);
    if (MB_SUCCESS != result) return result;
    hp_untreated.clear();

      // assign a category name to these sets
    static const char dual_category_names[2][CATEGORY_TAG_NAME_LENGTH] = 
      {"Chord\0", "Sheet\0"};
    
    result = mbImpl->tag_set_data(categoryTag, &this_hp, 1, dual_category_names[dim-1]);
    if (MB_SUCCESS != result) return result;

      // inner loop: traverse the hyperplane 'till we don't have any more
    MBRange tmp_star, star, tmp_range;
    while (0 != this_ent) {
      if (debug && hp_untreated.size()%10 == 0) 
        std::cout << "Dual surface " << hp_val << ", hp_untreated list size = " 
                  << hp_untreated.size() << "." << std::endl;

        // set the hp_val for this_ent
      result = mbImpl->tag_set_data(hp_tag, &this_ent, 1, &this_hp);
      if (MB_SUCCESS != result) return result;
      result = mbImpl->add_entities(this_hp, &this_ent, 1);
      if (MB_SUCCESS != result) return result;

        // get all neighbors connected to this entity
      tmp_range.clear(); tmp_star.clear(); star.clear();
      tmp_range.insert(this_ent);
      result = mbImpl->get_adjacencies(tmp_range, dim-1, true, tmp_star);
      if (MB_SUCCESS != result) return result;
      result = mbImpl->get_adjacencies(tmp_star, dim, false, star,
                                       MBInterface::UNION);
      if (MB_SUCCESS != result) return result;
      star.erase(this_ent);
      
        // for each star entity, see if it shares a cell with this one; if so,
        // it's not opposite

        // this_ent is already in tmp_range
      for (MBRange::iterator rit = star.begin(); rit != star.end(); rit++) {
          // check for tag first, 'cuz it's probably faster than checking adjacencies
        MBEntityHandle r_val;
        result = mbImpl->tag_get_data(hp_tag, &(*rit), 1, &r_val);
        if (MB_SUCCESS == result && 0 != r_val) 
          continue;

          // assign to avoid valgrind warning
        unsigned short tmp_mark = 0x0;
        result = mbImpl->tag_get_data(mark_tag, &(*rit), 1, &tmp_mark);
        if (MB_SUCCESS == result && mark_val == tmp_mark) 
          continue;

          // if it's on the loop, it's not eligible
        if (1 == dim && check_1d_loop_edge(*rit)) continue;

          // passed tag test; check adjacencies
        tmp_star.clear();
        tmp_range.insert(*rit);
        result = mbImpl->get_adjacencies(tmp_range, dim+1, false, tmp_star);
        if (MB_SUCCESS != result) return result;
        
        if (tmp_star.empty()) {
            // have one on this hp; just put it on the hp_untreated list for now,
            // will get tagged and put in the hp set later
          hp_untreated.push_back(*rit);
          result = mbImpl->tag_set_data(mark_tag, &(*rit), 1, &mark_val);
          if (MB_SUCCESS != result) return result;
        }

          // take *rit out of tmp_range, then proceed to the next
        tmp_range.erase(*rit);
      }
      
        // end of inner loop; get the next this_ent, or set to zero
      if (hp_untreated.empty()) this_ent = 0;
      else {
        this_ent = hp_untreated.back();
        hp_untreated.pop_back();
      }
    }
  }

  if (debug)
    std::cout << "Constructed " << hp_ids << " hyperplanes of dim = " << dim << "." << std::endl;
  return MB_SUCCESS;
}

bool DualTool::check_1d_loop_edge(MBEntityHandle this_ent) 
{
    // make sure it's an edge
  if (MBEDGE != mbImpl->type_from_handle(this_ent)) return false;

    // also has to be a dual entity
  unsigned char dum;
  MBErrorCode result = mbImpl->tag_get_data(isDualCell_tag(), &this_ent, 1, &dum);
  if (MB_SUCCESS != result || dum != 0x1) return false;
  
  const MBEntityHandle *verts;
  MBEntityHandle vert_tags[2];
  int num_verts;
  result = mbImpl->get_connectivity(this_ent, verts, num_verts);
  if (MB_SUCCESS != result) return false;
  
  result = mbImpl->tag_get_data(dualEntity_tag(), verts, 2, vert_tags);
  if (MB_SUCCESS != result ||
      mbImpl->type_from_handle(vert_tags[0]) != MBQUAD ||
      mbImpl->type_from_handle(vert_tags[1]) != MBQUAD) return false;
  
  else return true;
}

MBErrorCode DualTool::construct_hp_parent_child() 
{
  MBRange dual_surfs, dual_cells, dual_edges;
  MBErrorCode result = this->get_dual_hyperplanes(mbImpl, 2, dual_surfs);
  if (MB_SUCCESS != result || dual_surfs.empty()) return result;
  std::vector<MBEntityHandle> dual_curve_sets;
  
  for (MBRange::iterator surf_it = dual_surfs.begin(); surf_it != dual_surfs.end();
       surf_it++) {
      // get all the cells, edges in those cells, and chords for those edges
    dual_cells.clear();
    result = mbImpl->get_entities_by_handle(*surf_it, dual_cells);
    if (MB_SUCCESS != result) return result;
    dual_edges.clear();
    result = mbImpl->get_adjacencies(dual_cells, 1, false, dual_edges, MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    dual_curve_sets.reserve(dual_edges.size());
    result = mbImpl->tag_get_data(dualCurve_tag(), dual_edges, &dual_curve_sets[0]);
    if (MB_SUCCESS != result) return result;

      // reuse dual_cells to get unique list of chord sets
    dual_cells.clear();
    for (unsigned int i = 0; i < dual_edges.size(); i++)
      if (dual_curve_sets[i] != 0) dual_cells.insert(dual_curve_sets[i]);
    
      // now connect up this dual surf with all the 1d ones
    for (MBRange::iterator rit = dual_cells.begin(); rit != dual_cells.end(); rit++) {
      result = mbImpl->add_parent_child(*surf_it, *rit);
      if (MB_SUCCESS != result) return result;
    }
  }

  return MB_SUCCESS;
}
    
MBErrorCode DualTool::get_graphics_points(MBEntityHandle dual_ent,
                                          std::vector<int> &npts,
                                          std::vector<GraphicsPoint> &points) 
{
    // shouldn't be a set
  assert(MBENTITYSET != mbImpl->type_from_handle(dual_ent));
  
    // get the graphics points comprising the given entity
  GraphicsPoint gp_array[DualTool::GP_SIZE];

  MBErrorCode result = MB_SUCCESS;
  
    // switch based on topological dimension
  switch (mbImpl->dimension_from_handle(dual_ent)) {
    case 0:
        // just return the vertex point
      result = mbImpl->tag_get_data(dualGraphicsPoint_tag(), &dual_ent, 1,
                                    gp_array);
      if (MB_SUCCESS == result)
        points.push_back(gp_array[0]);
        
      break;

    case 1:
        // get my graphics point then those of my vertices
      const MBEntityHandle *connect;
      int num_connect;
      result = mbImpl->get_connectivity(dual_ent, connect, num_connect);
      if (MB_SUCCESS != result) break;
      
      result = mbImpl->tag_get_data(dualGraphicsPoint_tag(), connect, 2,
                                    gp_array);
      if (MB_SUCCESS == result) {
        points.push_back(gp_array[0]);
        points.push_back(gp_array[0]);
        points.push_back(gp_array[1]);
        result = mbImpl->tag_get_data(dualGraphicsPoint_tag(), &dual_ent, 1,
                                      gp_array);
        if (MB_SUCCESS == result) points[1] = gp_array[0];
      }

      npts.push_back(3);
      
      break;
      
    case 2:
      result = get_cell_points(dual_ent, npts, points);
      break;
  }
  
  return result;
}

MBErrorCode DualTool::get_cell_points(MBEntityHandle dual_ent,
                                      std::vector<int> &npts,
                                      std::vector<GraphicsPoint> &points) 
{
  assert(MBPOLYGON == mbImpl->type_from_handle(dual_ent));
  
    // get the 1cells in this 2cell
  MBRange one_cells;
  
  MBRange tc_range; tc_range.insert(dual_ent);
  MBErrorCode result = mbImpl->get_adjacencies(tc_range, 1, false, one_cells, 
                                               MBInterface::UNION); RR;

  int num_edges = one_cells.size();
  std::vector<GraphicsPoint> dum_gps(num_edges+1);
  
    // get graphics points for 0cells and for this cell
  result = mbImpl->tag_get_data(dualGraphicsPoint_tag(), one_cells,
                                &dum_gps[0]); RR;
  result = mbImpl->tag_get_data(dualGraphicsPoint_tag(), &dual_ent, 1,
                                &(dum_gps[num_edges])); RR;

  MBRange::iterator eit;
  const MBEntityHandle *connect;
  int num_connect;
  GraphicsPoint vert_gps[2];
  int i;
  for (i = 0, eit = one_cells.begin(); i < num_edges; i++, eit++) {
      // get the vertices and the graphics points for them
    result = mbImpl->get_connectivity(*eit, connect, num_connect); RR;
    result = mbImpl->tag_get_data(dualGraphicsPoint_tag(), connect, 2, 
                                  vert_gps); RR;

      // make the 2 tris corresponding to this edge; don't worry about order
      // for now
    npts.push_back(3);
    points.push_back(dum_gps[num_edges]);
    points.push_back(vert_gps[0]);
    points.push_back(dum_gps[i]);
    
    npts.push_back(3);
    points.push_back(dum_gps[num_edges]);
    points.push_back(dum_gps[i]);
    points.push_back(vert_gps[1]);
  }

  return result;
}

MBErrorCode DualTool::get_graphics_points(const MBRange &in_range,
                                          std::vector<GraphicsPoint> &points,
                                          const bool assign_ids,
                                          const int start_id) 
{
    // return graphics points on dual entities in in_range or in entities
    // in sets in in_range
  MBErrorCode result;

    // for each dual hyperplane set:
  MBRange::const_iterator rit;
  
  MBRange two_cells, all_cells;
  for (rit = in_range.begin(); rit != in_range.end(); rit++) {
      // for each entity:
    two_cells.clear();
    MBEntityType this_type = mbImpl->type_from_handle(*rit);
    if (MBENTITYSET == this_type) {
      result = mbImpl->get_entities_by_handle(*rit, two_cells); RR;
    
      std::copy(two_cells.begin(), two_cells.end(), mb_range_inserter(all_cells));
    }
    
    else {
      two_cells.insert(*rit);
      assert(this_type == MBVERTEX || this_type == MBEDGE ||
             this_type == MBPOLYGON || this_type == MBPOLYHEDRON);
    }

    result = mbImpl->get_adjacencies(two_cells, 0, false, all_cells, 
                                     MBInterface::UNION); RR;
    result = mbImpl->get_adjacencies(two_cells, 1, false, all_cells, 
                                     MBInterface::UNION); RR;
  }
      
      // get graphics points
  points.resize(all_cells.size());
  
  result = mbImpl->tag_get_data(dualGraphicsPointTag, all_cells,
                                &points[0]); RR;

  if (assign_ids) {
    int i = start_id;
    
    for (std::vector<GraphicsPoint>::iterator vit = points.begin(); 
         vit != points.end(); vit++)
      vit->id = i++;

    result = mbImpl->tag_set_data(dualGraphicsPoint_tag(), all_cells, 
                                  &points[0]); RR;
  }

  return result;
}

MBEntityHandle DualTool::next_loop_vertex(const MBEntityHandle last_v,
                                          const MBEntityHandle this_v,
                                          const MBEntityHandle dual_surf)
{
    // given two vertices, find the next one on the loop; if one is a dual
    // surface, then just choose either one for that surface
  assert((0 == last_v || mbImpl->type_from_handle(last_v) == MBVERTEX) &&
         mbImpl->type_from_handle(this_v) == MBVERTEX &&
         mbImpl->type_from_handle(dual_surf) == MBENTITYSET);

    // get the connected vertices
  MBSkinner skinner(mbImpl);
  MBRange other_verts;
  MBErrorCode result = skinner.get_bridge_adjacencies(this_v, 1, 0, other_verts);
  if (MB_SUCCESS != result || other_verts.empty()) return 0;
  
    //if (mbImpl->type_from_handle(last_v) == MBENTITYSET) {
      // dual surface, choose either; first get a 2cell on this surface
  MBRange tcells, tcells2, verts;
  result = mbImpl->get_entities_by_type(dual_surf, MBPOLYGON, tcells);
  if (MB_SUCCESS != result || tcells.empty()) return 0;

    // ok, pay attention here: first get 2cells common to dual surface and this_v
  verts.insert(this_v);
  result = mbImpl->get_adjacencies(verts, 2, false, tcells);
  if (MB_SUCCESS != result || tcells.empty()) return 0;

    // next get 2cells common to this_v and last_v, if last_v isn't zero, and remove
    // them from tcells
  if (0 != last_v) {
    verts.insert(last_v);
    result = mbImpl->get_adjacencies(verts, 2, false, tcells2);
    if (MB_SUCCESS != result || tcells2.empty()) return 0;
    MBRange temp_range = tcells.subtract(tcells2);
    tcells.swap(temp_range);
  }
  if (tcells.empty()) return 0;
  if (tcells.size() == 2) tcells.erase(tcells.begin());
  assert(tcells.size() == 1);
  
    // next loop vertex is intersection of other_verts and vertices adjacent to first
    // 2cell on tcells
  result = mbImpl->get_adjacencies(tcells, 0, false, other_verts);
  if (MB_SUCCESS != result || other_verts.size() != 1) return 0;
  
  return *other_verts.begin();
}

