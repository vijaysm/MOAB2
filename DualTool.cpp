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
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBSkinner.hpp"
#include "MBCore.hpp"
#include "MeshTopoUtil.hpp"
#include "AEntityFactory.hpp"
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <assert.h>

#define RR if (MB_SUCCESS != result) return result

bool debug = false;
bool debug_ap = true;

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
  
  result = mbImpl->tag_create(IS_DUAL_CELL_TAG_NAME, sizeof(unsigned int), MB_TAG_SPARSE, 
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

  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                              globalIdTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
}

DualTool::~DualTool() 
{}

  //! construct the dual entities for the entire mesh
MBErrorCode DualTool::construct_dual(MBEntityHandle *entities, 
                                     const int num_entities) 
{
    // allocate a dual entity for each primal entity in the mesh, starting
    // with highest dimension and working downward; do each dimension in a separate code
    // block, since they're all handled slightly differently
  
  MBRange regions, faces, edges, vertices;
  MBErrorCode result;

  if (NULL == entities || 0 == num_entities) {
    
      // first, construct all the aentities, since they're currently needed to 
      // compute the dual
    result = mbImpl->get_entities_by_dimension(0, 0, vertices);
    if (MB_SUCCESS != result) return result;

    result = MeshTopoUtil(mbImpl).construct_aentities(vertices);
    if (MB_SUCCESS != result) return result;

      // get all edges, faces and regions now, so we don't need to filter out dual
      // entities later
  
    result = mbImpl->get_entities_by_dimension(0, 1, edges);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->get_entities_by_dimension(0, 2, faces);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->get_entities_by_dimension(0, 3, regions);
    if (MB_SUCCESS != result) return result;
  }
  else {
      // get entities of various dimensions adjacent to these
    result = mbImpl->get_adjacencies(entities, num_entities, 0, true, vertices, 
                                     MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->get_adjacencies(entities, num_entities, 1, true, edges,
                                     MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->get_adjacencies(entities, num_entities, 2, true, faces,
                                     MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    result = mbImpl->get_adjacencies(entities, num_entities, 3, true, regions,
                                     MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
  
  }
  
  MBRange dual_verts;
  result = construct_dual_vertices(regions, dual_verts);
  if (MB_SUCCESS != result || dual_verts.size() != regions.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_verts.size() << " dual vertices." << std::endl;

    // don't really need dual edges, but construct 'em anyway
  MBRange dual_edges;
  result = construct_dual_edges(faces, dual_edges);
  if (MB_SUCCESS != result || dual_edges.size() != faces.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_edges.size() << " dual edges." << std::endl;

    // construct dual faces
  MBRange dual_faces;
  result = construct_dual_faces(edges, dual_faces);
  if (MB_SUCCESS != result || dual_faces.size() != edges.size()) return result;
  if (debug)
    std::cout << "Constructed " << dual_faces.size() << " dual faces." << std::endl;

    // construct dual cells
  MBRange dual_cells;
  result = construct_dual_cells(vertices, dual_cells);
  if (MB_SUCCESS != result || dual_cells.size() != vertices.size()) return result;
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
    if (MB_SUCCESS != tmp_result || out_ents.empty()) continue;

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
  MBRange equiv_edges;
#define TRC if (MB_SUCCESS != tmp_result) {result = tmp_result; continue;}
  for (rit = all_edges.begin(); rit != all_edges.end(); rit++) {
    
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS == tmp_result && 0 != dual_ent) {
      dual_ents.insert(dual_ent);
      continue;
    }
    
      // no dual entity; construct one; get the dual vertices bounding the edge in radial order,
      // then construct the dual face
    std::vector<MBEntityHandle> rad_dverts;
    bool bdy_edge;
    tmp_result = get_radial_dverts(*rit, rad_dverts, bdy_edge);TRC
    if (rad_dverts.empty()) continue;
    
    tmp_result = mbImpl->create_element(MBPOLYGON, &rad_dverts[0], rad_dverts.size(), dual_ent);TRC

      // tag it indicating it's a dual entity, and tag primal/dual with dual/primal
    tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);TRC
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);TRC
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &dual_ent, 1, &(*rit));TRC

      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
      // add a graphics point to the cell; position depends on whether it's a
      // bdy cell (mid-pt of cell's vertices) or not (mid-pt of primal edge)
    double avg_pos[3];
    tmp_result = MeshTopoUtil(mbImpl).get_average_position(*rit, avg_pos);TRC
    if (bdy_edge) {
      
        // add a new dual edge betw last 2 verts
      MBEntityHandle new_edge;
      tmp_result = mbImpl->create_element(MBEDGE, &rad_dverts[rad_dverts.size()-2], 
                                          2, new_edge);TRC
      tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &new_edge, 1, &is_dual);TRC

        // tag the new dual edge with the primal edge as it's dual entity; primal
        // edge IS NOT likewise tagged, since it's already tagged with the 2cell
      tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &new_edge, 1, &(*rit)); TRC
      
        // add a graphics pt, position is center of primal edge
      tmp_result = add_graphics_point(dual_ent);TRC
      tmp_result = add_graphics_point(new_edge, avg_pos);TRC
    }
    
    else {
        // if inside, point goes on the 2cell, at primal edge mid-pt
      tmp_result = add_graphics_point(dual_ent, avg_pos);TRC
    }

      // check to see whether we have equiv entities; if we find any, save for later fixup
    MBRange dum_edges, dum_poly(dual_ent, dual_ent);
    tmp_result = mbImpl->get_adjacencies(dum_poly, 1, false, dum_edges);
    if (MB_MULTIPLE_ENTITIES_FOUND == tmp_result) {
        // we do - need to add adjacencies to disambiguate; use the primal
      equiv_edges.merge(dum_edges);
    }    
  }

  if (!equiv_edges.empty()) 
    result = check_dual_equiv_edges(equiv_edges);

  return result;
}

MBErrorCode DualTool::check_dual_equiv_edges(MBRange &dual_edges)
{
    // fix equivalent dual edges (i.e. edges whose vertices define multiple edges)
    // by explicitly adding adjacencies to containing polygons; adjacent polygons
    // found by going through primal
  MBErrorCode tmp_result, result = MB_SUCCESS;

  MBRange all_dedges(dual_edges);
    // first, go through all dual edges and find equivalent edges (by looking for
    // up-adjacent edges on the vertices of each edge)
  for (MBRange::iterator rit = dual_edges.begin(); rit != dual_edges.end(); rit++) {
    MBRange connect, dum_range(*rit, *rit);
    tmp_result = mbImpl->get_adjacencies(dum_range, 0, false, connect);
    if (MB_SUCCESS != tmp_result) continue;
    tmp_result = mbImpl->get_adjacencies(connect, 1, false, all_dedges, MBInterface::UNION);
    if (MB_SUCCESS != tmp_result) continue;
  }

    // save a copy for checking later
  MBRange save_all_2cells;

    // go through each edge
  while (!all_dedges.empty()) {
    MBEntityHandle this_edge = *all_dedges.begin();
    all_dedges.erase(all_dedges.begin());
    
    const MBEntityHandle *connect;
    int num_connect;
    result = mbImpl->get_connectivity(this_edge, connect, num_connect);
    if (MB_SUCCESS != result) continue;

    MBRange dum_edges, verts;
    verts.insert(connect[0]);
    verts.insert(connect[1]);
    MBErrorCode tmp_result = mbImpl->get_adjacencies(verts, 1, false, dum_edges);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      continue;
    }
    if (dum_edges.size() == 1) {
        // not an equiv edge - already removed from list, so just continue
      continue;
    }
    
      // ok, have an equiv entity - fix by looking through primal
      // pre-get the primal of these
    MBEntityHandle dedge_quad;
    tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &this_edge, 1, &dedge_quad);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      continue;
    }
  
    if (MBQUAD == mbImpl->type_from_handle(dedge_quad)) {
      
        // get the primal edges adj to quad
      MBRange dum_quad_range(dedge_quad, dedge_quad), adj_pedges;
      tmp_result = mbImpl->get_adjacencies(dum_quad_range, 1, false, adj_pedges);
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        continue;
      }
        // get the dual 2cells corresponding to those pedges
      std::vector<MBEntityHandle> dcells;
      dcells.resize(adj_pedges.size());
      tmp_result = mbImpl->tag_get_data(dualEntity_tag(), adj_pedges, &dcells[0]);
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        continue;
      }
        // now add explicit adjacencies from the dedge to those dcells
      std::vector<MBEntityHandle>::iterator vit;
      for (vit = dcells.begin(); vit != dcells.end(); vit++) {
        save_all_2cells.insert(*vit);
        
        assert(MBPOLYGON == mbImpl->type_from_handle(*vit));
        tmp_result = mbImpl->add_adjacencies(this_edge, &(*vit), 1, false);
        if (MB_SUCCESS != tmp_result) {
          result = tmp_result;
          continue;
        }
          // check that there are really adjacencies and *vit is in them
        const MBEntityHandle *adjs;
        int num_adjs;
        tmp_result = reinterpret_cast<MBCore*>(mbImpl)->a_entity_factory()->
          get_adjacencies(this_edge, adjs, num_adjs);
        if (NULL == adjs || std::find(adjs, adjs+num_adjs, *vit) == adjs+num_adjs)
          std::cout << "Add_adjacencies failed in construct_dual_faces." << std::endl;
      }
    }
    else {
        // else, have a dual edge representing a bdy edge - tie directly to
        // dual entity if its dual entity
      MBEntityHandle bdy_dcell;
      tmp_result = mbImpl->tag_get_data(dualEntity_tag(), &dedge_quad, 1, &bdy_dcell); TRC
      assert(MBPOLYGON == mbImpl->type_from_handle(bdy_dcell));
      
      tmp_result = mbImpl->add_adjacencies(this_edge, &bdy_dcell, 1, false); 
      if (MB_SUCCESS != tmp_result) {
        result = tmp_result;
        continue;
      }
    }
  }
  

    // sanity check - look for adj edges again, and check for equiv entities
  for (MBRange::iterator vit = save_all_2cells.begin(); 
       vit != save_all_2cells.end(); vit++) {
    MBRange adj_edges, dum_quad_range;
    dum_quad_range.insert(*vit);
    assert(MBPOLYGON == mbImpl->type_from_handle(*vit));
    tmp_result = mbImpl->get_adjacencies(dum_quad_range, 1, false, adj_edges);
    if (MB_MULTIPLE_ENTITIES_FOUND == tmp_result) {
      std::cout << "Multiple entities returned for polygon " << mbImpl->id_from_handle(*vit)
                << "." << std::endl;
      continue;
    }
  }
    // success!
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
  
  std::vector<MBEntityHandle> rad_faces, rad_ents;
  MBErrorCode result = MeshTopoUtil(mbImpl).star_entities(edge, rad_faces, bdy_edge, 0, 
                                                          &rad_ents);
  if (MB_SUCCESS != result) return result;

  if (bdy_edge) {
      // if we're a bdy edge, change the order back to what DualTool expects
    rad_ents.push_back(*rad_faces.rbegin());
    rad_ents.push_back(*rad_faces.begin());
  }
  
  rad_dverts.resize(rad_ents.size());
  for (unsigned int i = 0; i < rad_ents.size(); i++) {
    MBEntityHandle dual_ent;
    result = mbImpl->tag_get_data(dualEntity_tag(), &rad_ents[i], 1,
                                  &dual_ent);
    if (!bdy_edge || i < rad_ents.size()-2) rad_dverts[i] = dual_ent;
    else {
      // fix up this entry
      assert(mbImpl->type_from_handle(dual_ent) == MBEDGE);
    
        // get connectivity of that edge
      const MBEntityHandle *connect;
      int num_connect;
      result = mbImpl->get_connectivity(dual_ent, connect, num_connect);
      if (MB_SUCCESS != result) return result;
    
        // we want the one that's not already on the list; reuse last_face
      int last_hex = (i == rad_ents.size()-1 ? 0 : i-1);
      MBEntityHandle last_face = (connect[0] == rad_dverts[last_hex] ? connect[1] : connect[0]);
      rad_dverts[i] = last_face;
    }
  }

  return result;
}

  //! construct the dual entities for a hex mesh, including dual surfaces & curves
MBErrorCode DualTool::construct_hex_dual(MBRange &entities) 
{
  std::vector<MBEntityHandle> evec;
  std::copy(entities.begin(), entities.end(), std::back_inserter(evec));
  return construct_hex_dual(&evec[0], evec.size());
}
  
  //! construct the dual entities for a hex mesh, including dual surfaces & curves
MBErrorCode DualTool::construct_hex_dual(MBEntityHandle *entities,
                                         const int num_entities) 
{
    // really quite simple: 

    // construct the dual...
  MBErrorCode result = construct_dual(entities, num_entities);
  if (MB_SUCCESS != result) return result;
  
    // now traverse to build 1d and 2d hyperplanes
  result = construct_dual_hyperplanes(1, entities, num_entities);
  if (MB_SUCCESS != result) return result;

  result = construct_dual_hyperplanes(2, entities, num_entities);
  if (MB_SUCCESS != result) return result;
  
  result = construct_hp_parent_child();
  if (MB_SUCCESS != result) return result;

    // see?  simple, just like I said
  return MB_SUCCESS;
}

  //! get the cells of the dual
MBErrorCode DualTool::get_dual_entities(const int dim, 
                                        MBEntityHandle *entities, 
                                        const int num_entities, 
                                        MBRange &dual_ents) 
{
  if (0 == isDualCell_tag()) return MB_SUCCESS;
  if (0 > dim || 3 < dim) return MB_INDEX_OUT_OF_RANGE;

  unsigned int dum = 0x1;
  const void *dum_ptr = &dum;
  static MBEntityType dual_type[] = {MBVERTEX, MBEDGE, MBPOLYGON, MBPOLYHEDRON};

  MBRange dim_ents;
  
  MBErrorCode result;

  if (0 == entities || 0 == num_entities) {
      // just get all the dual entities of this dimension
    result = mbImpl->get_entities_by_type_and_tag(0, dual_type[dim], 
                                                  &isDualCellTag, &dum_ptr, 1,
                                                  dual_ents);
  }
  else {
      // else look for specific dual entities
    result = mbImpl->get_adjacencies(entities, num_entities, 3-dim, false,
                                     dim_ents, MBInterface::UNION);
    if (MB_SUCCESS != result) return result;
    std::vector<MBEntityHandle> dual_ents_vec(dim_ents.size());
    result = mbImpl->tag_get_data(dualEntity_tag(), dim_ents, &dual_ents_vec[0]);
    if (MB_SUCCESS != result) return result;
    std::copy(dual_ents_vec.begin(), dual_ents_vec.end(), 
              mb_range_inserter(dual_ents));
  }
  
  return result;
}

  //! get the faces of the dual
MBErrorCode DualTool::get_dual_entities(const int dim, 
                                        MBEntityHandle *entities, 
                                        const int num_entities, 
                                        std::vector<MBEntityHandle> &dual_ents) 
{
  MBRange tmp_range;
  MBErrorCode result = get_dual_entities(dim, entities, num_entities, tmp_range);
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

MBErrorCode DualTool::construct_dual_hyperplanes(const int dim, 
                                                 MBEntityHandle *entities, 
                                                 const int num_entities) 
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
    mbImpl->get_number_entities_by_type(0, MBHEX, num_hexes) != MB_SUCCESS ||
      // if we're asking for 1d dual ents, should be quads around
    (num_quads == 0 && dim == 1) ||
      // if we're asking for 2d dual ents, should be hexes around
    (num_hexes == 0 && dim == 2))
    return MB_FAILURE;
  
    // get tag name for this dimension hyperplane
  MBTag gid_tag;
  int dum = -1;
  MBErrorCode result = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, 4, 
                                                        MB_TAG_DENSE, gid_tag, &dum);

  MBTag hp_tag = (1 == dim ? dualCurve_tag() : dualSurface_tag());
  
    // two stacks: one completely untreated entities, and the other untreated 
    // entities on the current dual hyperplane
  std::vector<MBEntityHandle> tot_untreated;
  
    // put dual entities of this dimension on the untreated list
  result = get_dual_entities(dim, entities, num_entities, tot_untreated);
  if (MB_SUCCESS != result) 
    return result;
  
    // main part of traversal loop
  MBEntityHandle this_ent;
  MBEntityHandle this_hp;
  std::vector<MBEntityHandle> parents;

  while (!tot_untreated.empty()) {
    if (debug && dim == 2 /*(tot_untreated.size()%report == 0)*/)
      std::cout << "Untreated list size " << tot_untreated.size() << "." << std::endl;
      
    this_ent = tot_untreated.back(); tot_untreated.pop_back();
    result = mbImpl->tag_get_data(hp_tag, &this_ent, 1, &this_hp);
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) 
      return result;

      // test for this entity having a hyperplane assignment already
    else if (this_hp != 0)
      continue;

    if (1 == dim && check_1d_loop_edge(this_ent)) continue;
    
      // inner loop: traverse the hyperplane 'till we don't have any more
    result = traverse_hyperplane(hp_tag, this_hp, this_ent);
    if (MB_SUCCESS != result) {
      std::cout << "Failed to traverse hyperplane ";
      if (this_hp) std::cout << mbImpl->id_from_handle(this_hp) << "."  << std::endl;
      else std::cout << "0." << std::endl;
      return result;
    }

      // ok, now order the edges if it's a chord
    if (1 == dim) order_chord(this_hp);
  }

  return MB_SUCCESS;
}

MBErrorCode DualTool::traverse_hyperplane(const MBTag hp_tag, 
                                          MBEntityHandle &this_hp, 
                                          MBEntityHandle this_ent) 
{
  MBRange tmp_star, star, tmp_range, new_hyperplane_ents;
  std::vector<MBEntityHandle> hp_untreated;
  int dim = mbImpl->dimension_from_handle(this_ent);
  MeshTopoUtil mtu(mbImpl);
  this_hp = 0;
  MBErrorCode result;
  
  unsigned short mark_val = 0x0;
  MBTag mark_tag;
  result = mbImpl->tag_get_handle("__hyperplane_mark", mark_tag);
  if (MB_SUCCESS != result) {
    result = mbImpl->tag_create("__hyperplane_mark", 1, 
                                MB_TAG_BIT, mark_tag, &mark_val);
    if (MB_SUCCESS != result) 
      return result;
  }
  mark_val = 0x1;
  
  while (0 != this_ent) {
    MBEntityHandle tmp_hp = get_dual_hyperplane(this_ent);
    if (0 == this_hp && 0 != tmp_hp) this_hp = tmp_hp;
    
    if (0 == tmp_hp) new_hyperplane_ents.insert(this_ent);
    
    if (debug && hp_untreated.size()%10 == 0) 
      std::cout << "Dual surface " << this_hp << ", hp_untreated list size = " 
                << hp_untreated.size() << "." << std::endl;

      // get the 2nd order adjacencies through lower dimension
    tmp_range.clear(); tmp_star.clear(); star.clear();
    result = mtu.get_bridge_adjacencies(this_ent, dim-1, dim, star); RR;

      // get the bridge adjacencies through higher dimension
    result = mtu.get_bridge_adjacencies(this_ent, dim+1, dim, tmp_star); RR;
    tmp_range = star.subtract(tmp_star);
    
    for (MBRange::iterator rit = tmp_range.begin(); rit != tmp_range.end(); rit++) {
      if (new_hyperplane_ents.find(*rit) != new_hyperplane_ents.end()) continue;
      
        // check for tag first, 'cuz it's probably faster than checking adjacencies
        // assign to avoid valgrind warning
      unsigned short tmp_mark = 0x0;
      result = mbImpl->tag_get_data(mark_tag, &(*rit), 1, &tmp_mark);
      if (MB_SUCCESS == result && mark_val == tmp_mark) 
        continue;

        // if it's on the loop, it's not eligible
      if (1 == dim && check_1d_loop_edge(*rit)) continue;

        // have one on this hp; just put it on the hp_untreated list for now,
        // will get tagged and put in the hp set later
      hp_untreated.push_back(*rit);
      result = mbImpl->tag_set_data(mark_tag, &(*rit), 1, &mark_val);
      if (MB_SUCCESS != result) 
        return result;
    }

      // end of inner loop; get the next this_ent, or set to zero
    if (hp_untreated.empty()) this_ent = 0;
    else {
      this_ent = hp_untreated.back();
      hp_untreated.pop_back();
    }
  }

  if (debug_ap) {
    std::string hp_name;
    if (2 == dim) hp_name = "sheet";
    else hp_name = "chord";
    
    if (0 == this_hp) std::cout << "Constructed new " << hp_name << " with ";
    else {
      int this_id;
      result = mbImpl->tag_get_data(globalId_tag(), &this_hp, 1, &this_id); RR;
      std::cout << "Added to " << hp_name << " " << this_id << " ";
    }
    if (dim == 2) std::cout << "edges:" << std::endl;
    else std::cout << "quads:" << std::endl;
    std::vector<MBEntityHandle> pents(new_hyperplane_ents.size());
    result = mbImpl->tag_get_data(dualEntity_tag(), new_hyperplane_ents,
                                  &pents[0]); RR;
    for (std::vector<MBEntityHandle>::iterator vit = pents.begin(); 
         vit != pents.end(); vit++) {
      if (vit != pents.begin()) std::cout << ", ";
      std::cout << mbImpl->id_from_handle(*vit);
    }
    std::cout << std::endl;
  }

  if (0 == this_hp) {
      // ok, doesn't have one; make a new hyperplane
    int new_id = -1;
    result = construct_new_hyperplane(dim, this_hp, new_id);
    if (MB_SUCCESS != result) return result;

    if (debug_ap) {
      std::cout << "New ";
      if (2 == dim) std::cout << " sheet ";
      else std::cout << " chord ";
      std::cout << new_id << " constructed." << std::endl;
    }
  }

    // set the hp_val for entities which didn't have one before
  std::vector<MBEntityHandle> hp_tags(new_hyperplane_ents.size());
  std::fill(hp_tags.begin(), hp_tags.end(), this_hp);
  result = mbImpl->tag_set_data(hp_tag, new_hyperplane_ents, &hp_tags[0]);
  if (MB_SUCCESS != result) 
    return result;
  result = mbImpl->add_entities(this_hp, new_hyperplane_ents);
  if (MB_SUCCESS != result) 
    return result;

    // unmark the entities by removing the tag
  result = mbImpl->tag_delete(mark_tag);
  if (MB_SUCCESS != result) 
    return result;

  return MB_SUCCESS;
}

MBErrorCode DualTool::order_chord(MBEntityHandle chord_set) 
{
    // re-order the 1cells in the set so they are in order along the chord
    // start by finding the vertex dual to a quad
  MBRange verts, one_cells;
  MBErrorCode result = mbImpl->get_entities_by_dimension(chord_set, 1, one_cells); 
  if (MB_SUCCESS != result || one_cells.empty()) return MB_FAILURE;
  
  result = mbImpl->get_adjacencies(one_cells, 0, false, verts, MBInterface::UNION);
  if (MB_SUCCESS != result || verts.empty()) return MB_FAILURE;
  
  MBEntityHandle last_vert = 0;
  for (MBRange::iterator rit = verts.begin(); rit != verts.end(); rit++) {
    if (TYPE_FROM_HANDLE(get_dual_entity(*rit)) == MBQUAD) {
      last_vert = *rit;
      break;
    }
  }
    // if there's no vertex owned by a quad, just start with 1st one
  if (0 == last_vert) last_vert = *verts.begin();
  
    // now, skip from vertex to vertex, building a list of 1cells
  std::vector<MBEntityHandle> ordered_1cells;
  MBEntityHandle last_1cell = 0;
  MBRange dum1, dum2;
  const MBEntityHandle *connect;
  int num_connect;
  MBErrorCode tmp_result = MB_SUCCESS;
  while(ordered_1cells.size() != one_cells.size()) {
    dum1 = one_cells;
    result = mbImpl->get_adjacencies(&last_vert, 1, 1, false, dum1);
    if (0 != last_1cell) dum1.erase(last_1cell);
      // assert(1 == dum1.size());
    if (1 != dum1.size()) {
      std::cerr << "unexpected size traversing chord." << std::endl;
      tmp_result = MB_FAILURE;
    }
      
    last_1cell = *dum1.begin();
    ordered_1cells.push_back(last_1cell);
    result = mbImpl->get_connectivity(last_1cell, connect, num_connect); RR;
    if (last_vert == connect[0]) last_vert = connect[1];
    else last_vert = connect[0];
  }
  
    // now have the 1cells in order, replace them in the set
  if (MB_SUCCESS == tmp_result) {
    result = mbImpl->remove_entities(chord_set, one_cells); RR;
    result = mbImpl->add_entities(chord_set, &ordered_1cells[0], ordered_1cells.size()); RR;
  }
  
  return MB_SUCCESS;
}

MBErrorCode DualTool::construct_new_hyperplane(const int dim,
                                               MBEntityHandle &new_hyperplane,
                                               int &id) 
{
  MBErrorCode result;
  if (1 == dim)
    result = mbImpl->create_meshset((MESHSET_ORDERED | MESHSET_TRACK_OWNER), 
                                    new_hyperplane);
  else
    result = mbImpl->create_meshset((MESHSET_SET | MESHSET_TRACK_OWNER), 
                                    new_hyperplane);
  if (MB_SUCCESS != result) return result;

  if (-1 == id) {
    MBRange all_hyperplanes;
    result = get_dual_hyperplanes(mbImpl, dim, all_hyperplanes); RR;
    id = all_hyperplanes.size() + 1;
  }
    
  result = mbImpl->tag_set_data(globalId_tag(), &new_hyperplane, 1, &id); RR;
  MBTag hp_tag = (1 == dim ? dualCurve_tag() : dualSurface_tag());
  result = mbImpl->tag_set_data(hp_tag, &new_hyperplane, 1, &new_hyperplane);

    // assign a category name to these sets
  static const char dual_category_names[2][CATEGORY_TAG_NAME_LENGTH] = 
    {"Chord\0", "Sheet\0"};
    
  result = mbImpl->tag_set_data(categoryTag, &new_hyperplane, 1, dual_category_names[dim-1]);

  return result;
}

bool DualTool::check_1d_loop_edge(MBEntityHandle this_ent) 
{
    // make sure it's an edge
  if (MBEDGE != mbImpl->type_from_handle(this_ent)) return false;

    // also has to be a dual entity
  unsigned int dum;
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
  MeshTopoUtil tpu(mbImpl);
  MBRange other_verts;
  MBErrorCode result = tpu.get_bridge_adjacencies(this_v, 1, 0, other_verts);
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

    // next get vertices common to both 2cells and subtract from other_verts; also
    // remove last_v if it's non-zero
  verts.clear();
  result = mbImpl->get_adjacencies(tcells, 0, false, verts);
  if (MB_SUCCESS != result || verts.empty()) return 0;
  
  MBRange tmp_verts = other_verts.subtract(verts);
  other_verts.swap(tmp_verts);
  if (0 != last_v) other_verts.erase(last_v);

    // now get intersection of remaining vertices and 2 2cells vertices
    // look at each one successively; should find one, maybe not on both
  tmp_verts = other_verts;
  MBRange tmp_faces(*tcells.begin(), *tcells.begin());
  result = mbImpl->get_adjacencies(tmp_faces, 0, false, tmp_verts);
  if (MB_SUCCESS == result && !tmp_verts.empty()) return *tmp_verts.begin();
  tmp_faces.clear();
  tmp_faces.insert(*tcells.rbegin());
  result = mbImpl->get_adjacencies(tmp_faces, 0, false, other_verts);
  if (MB_SUCCESS == result && !other_verts.empty()) return *other_verts.begin();

    // if we got here, there isn't any
  return 0;
}

MBEntityHandle DualTool::get_dual_hyperplane(const MBEntityHandle ncell) 
{
    // get the sheet or chord it's in
  std::vector<MBEntityHandle> adj_sets;
  MBErrorCode result = mbImpl->get_adjacencies(&ncell, 1, 4, false, adj_sets);
  if (MB_SUCCESS != result) return 0;
    
  MBEntityHandle dum_set;
  for (std::vector<MBEntityHandle>::iterator vit = adj_sets.begin(); 
       vit != adj_sets.end(); vit++) {
    if (mbImpl->tag_get_data(dualCurve_tag(), &(*vit), 1, &dum_set) != MB_TAG_NOT_FOUND ||
        mbImpl->tag_get_data(dualSurface_tag(), &(*vit), 1, &dum_set) != MB_TAG_NOT_FOUND)
      return *vit;
  }
  
  return 0;
}

    //! set the dual surface or curve for an entity
MBErrorCode DualTool::set_dual_surface_or_curve(MBEntityHandle entity, 
                                                const MBEntityHandle dual_hyperplane,
                                                const int dual_entity_dimension)
{
  if (1 == dual_entity_dimension)
    mbImpl->tag_set_data(dualCurve_tag(), &entity, 1, &dual_hyperplane);
  else if (2 == dual_entity_dimension)
    mbImpl->tag_set_data(dualSurface_tag(), &entity, 1, &dual_hyperplane);
  else 
    return MB_INDEX_OUT_OF_RANGE;

  return MB_SUCCESS;
}

//! return the corresponding dual entity
MBEntityHandle DualTool::get_dual_entity(const MBEntityHandle this_ent) const
{
  MBEntityHandle dual_ent;
  MBErrorCode result = mbImpl->tag_get_data(dualEntity_tag(), &this_ent, 1, &dual_ent);
  if (MB_SUCCESS != result || MB_TAG_NOT_FOUND == result) return 0;
  else return dual_ent;
}

//! return the corresponding dual entity
MBEntityHandle DualTool::get_extra_dual_entity(const MBEntityHandle this_ent) 
{
  MBEntityHandle dual_ent;
  MBErrorCode result = mbImpl->tag_get_data(extraDualEntity_tag(), &this_ent, 1, &dual_ent);
  if (MB_SUCCESS != result || MB_TAG_NOT_FOUND == result) return 0;
  else return dual_ent;
}

MBErrorCode DualTool::atomic_pillow(MBEntityHandle odedge, MBEntityHandle &new_hp) 
{
    // perform an atomic pillow operation around dedge

    // 0. get star 2cells and 3cells around odedge (before odedge changes)
  MeshTopoUtil mtu(mbImpl);
  MBRange star_2cells, star_3cells;
  MBErrorCode result = mbImpl->get_adjacencies(&odedge, 1, 2, false, star_2cells); RR;
  result = mbImpl->get_adjacencies(&odedge, 1, 3, false, star_3cells); RR;
  
    // tear down the dual entities which will be modified by the ap first
  result = delete_dual_entities(star_3cells);RR;
  result = delete_dual_entities(star_2cells);RR;

    // grab the quad before deleting the odedge
  MBEntityHandle quad = get_dual_entity(odedge);
  assert(0 != quad);
  result = delete_dual_entities(&odedge, 1); RR;

    // now change the quad to an ap
  std::vector<MBEntityHandle> verts;
  result = mbImpl->get_connectivity(&quad, 1, verts); RR;
  
    // get average position of vertices
  double coords[12], avg[3] = {0.0, 0.0, 0.0};
  result = mbImpl->get_coords(&verts[0], verts.size(), coords); RR;
  for (int i = 0; i < 4; i++) {
    avg[0] += coords[3*i]; avg[1] += coords[3*i+1]; avg[2] += coords[3*i+2];
  }
  for (int i = 0; i < 3; i++) avg[i] *= 0.25;

    // for each position, get a corresponding position 1/2 way to avg
  double new_coords[12];
  for (int i = 0; i < 4; i++) {
    new_coords[3*i] = avg[0] + .5*(coords[3*i]-avg[0]);
    new_coords[3*i+1] = avg[1] + .5*(coords[3*i+1]-avg[1]);
    new_coords[3*i+2] = avg[2] + .5*(coords[3*i+2]-avg[2]);
  }
  
    // make the 4 new vertices; store in vector long enough for hex connectivity
  for (int i = 0; i < 4; i++) {
    verts.push_back(0);
    result = mbImpl->create_vertex(&coords[3*i], verts[4+i]); RR;
  }

    // get the hexes connected to the quad
  MBRange hexes;
  result = mbImpl->get_adjacencies(&quad, 1, 3, false, hexes); RR;
  
    // remove any explicit adjacency from the first hex, since that'll get connected
    // to the new outer quad; add adjacency between quad and other hex
  result = mbImpl->remove_adjacencies(quad, &(*hexes.begin()), 1); RR;
  if (hexes.size() == 2) {
    result = mbImpl->add_adjacencies(quad, &(*hexes.rbegin()), 1, false);
    RR;
  }
  
    // create the new, outer quad, and make it explicitly adjacent to 1st hex
  MBEntityHandle new_quad;
  result = mbImpl->create_element(MBQUAD, &verts[0], 4, new_quad); RR;
  result = mbImpl->add_adjacencies(new_quad, &(*hexes.begin()), 1, false); RR;
  
    // now make two inner hexes, connect each to one of the quads; note connectivity
    // array is flipped for the two hexes
  MBEntityHandle new_hexes[2];
  result = mbImpl->create_element(MBHEX, &verts[0], 8, new_hexes[0]); RR;
  result = mbImpl->add_adjacencies(quad, &new_hexes[0], 1, false); RR;
  
    // reverse the connectivities for the 2nd hex
  std::reverse(verts.begin(), verts.begin()+4);
  std::reverse(verts.begin()+4, verts.end());
  result = mbImpl->create_element(MBHEX, &verts[0], 8, new_hexes[1]); RR;
  result = mbImpl->add_adjacencies(new_quad, &new_hexes[1], 1, false); RR;

    // now update the dual
  result = construct_hex_dual(&new_hexes[0], 2); RR;

    // get the new dual surface, by getting one of the edges between the center
    // and outer vertex rings
  MBRange new_edge;
  verts[1] = verts[4];
  result = mbImpl->get_adjacencies(&verts[0], 2, 1, false, new_edge);
  if (MB_SUCCESS != result || new_edge.size() != 1) return result;
  new_hp = get_dual_hyperplane(get_dual_entity(*new_edge.begin()));
  
  return MB_SUCCESS;
}

  //! effect reverse atomic pillow operation
MBErrorCode DualTool::rev_atomic_pillow(MBEntityHandle pillow, MBRange &chords) 
{
    // get the dual entities associated with elements in the pillow; go through
    // the elements instead of the pillow sheet so you get all of them, not just
    // the ones on the sheet
  MBRange dverts;
  MBErrorCode result = get_dual_entities(pillow, NULL, NULL,
                                         &dverts, NULL, NULL);
  if (MB_SUCCESS != result) return result;
  assert(2 == dverts.size());
  
  MBEntityHandle hexes[2];
  result = mbImpl->tag_get_data(dualEntity_tag(), dverts, hexes); RR;
  assert(hexes[0] != 0 && hexes[1] != 0);

  std::vector<MBEntityHandle> dcells[4];
  MBRange pcells[4];
  std::copy(hexes, hexes+2, mb_range_inserter(pcells[3]));
  std::copy(dverts.begin(), dverts.end(), std::back_inserter(dcells[0]));
  for (int dim = 0; dim <= 2; dim++) {
    result = mbImpl->get_adjacencies(hexes, 2, dim, false, pcells[dim], 
                                     MBInterface::UNION); RR;
    dcells[3-dim].resize(pcells[dim].size());
    result = mbImpl->tag_get_data(dualEntity_tag(), pcells[dim], &dcells[3-dim][0]); RR;
  }
  
    // delete the dual entities which are part of the original pillow
  result = mbImpl->delete_entities(&pillow, 1);
  if (MB_SUCCESS != result) return result;
  
  result = mbImpl->delete_entities(chords);
  if (MB_SUCCESS != result) return result;

  for (int i = 3; i >= 0; i--) {
    result = delete_dual_entities(&dcells[i][0], dcells[i].size()); RR;
  }

    // delete the primal entities inserted by the ap; be careful to get the right
    // faces, edges and vertices
  MBRange del_faces, del_edges, del_verts, tmp_faces, tmp_verts;
    // faces are the shared 5 and the 1 other one with greater handle (which
    // we know will be later in the range)
  result = mbImpl->get_adjacencies(hexes, 2, 2, false, del_faces); RR;
  assert(5 == del_faces.size());
  std::copy(pcells[2].begin(), pcells[2].end(), mb_range_inserter(tmp_faces));
  tmp_faces = tmp_faces.subtract(del_faces);
  del_faces.insert(*tmp_faces.rbegin());
  result = mbImpl->get_adjacencies(tmp_faces, 0, false, tmp_verts); RR;
  std::copy(pcells[0].begin(), pcells[0].end(), mb_range_inserter(del_verts));
  del_verts = del_verts.subtract(tmp_verts);
  assert(4 == del_verts.size());
  result = mbImpl->get_adjacencies(del_verts, 1, false, del_edges, MBInterface::UNION); RR;
  assert(8 == del_edges.size());
  
  result = mbImpl->delete_entities(hexes, 2); RR;
  result = mbImpl->delete_entities(del_faces); RR;
  result = mbImpl->delete_entities(del_edges); RR;
  result = mbImpl->delete_entities(del_verts); RR;

    // recompute the dual for the hexes on either side of the quad affected
    // by the ap removal
  MBRange tmp_hexes;
  result = mbImpl->get_adjacencies(tmp_verts, 3, false, tmp_hexes, MBInterface::UNION); RR;
  result = construct_hex_dual(tmp_hexes); RR;

  return MB_SUCCESS;
}

  
  
  
  

MBErrorCode DualTool::delete_dual_entities(MBEntityHandle *entities, 
                                           const int num_entities) 
{
  MBEntityHandle null_entity = 0;
  MBErrorCode result;
  std::vector<MBEntityHandle> ents_to_delete;
  
  for (int i = 0; i < num_entities; i++) {
      // reset the primal's dual entity
    MBEntityHandle primal = get_dual_entity(entities[i]);
    if (get_dual_entity(primal) == entities[i]) {
      result = mbImpl->tag_set_data(dualEntity_tag(), &primal, 1, &null_entity); RR;
    }
    MBEntityHandle extra = get_extra_dual_entity(primal);
    if (0 != extra) {
      result = mbImpl->tag_set_data(extraDualEntity_tag(), &primal, 1, &null_entity); RR;
    }
    ents_to_delete.push_back(entities[i]);
    
      // check for extra dual entities
    if (mbImpl->type_from_handle(entities[i]) == MBPOLYGON) {
      // for 2cell, might be a loop edge
      MBRange loop_edges;
      result = mbImpl->get_adjacencies(&entities[i], 1, 1, false, loop_edges);
      for (MBRange::iterator rit = loop_edges.begin(); rit != loop_edges.end(); rit++) {
        if (check_1d_loop_edge(*rit)) {
          MBEntityHandle this_ent = *rit;
          result = delete_dual_entities(&this_ent, 1); RR;
        }
      }
    }
    else if (extra && extra != entities[i])
        // just put it on the list; primal for which we're extra has already been
        // reset to not point to extra entity
      ents_to_delete.push_back(extra);
  }

    // now delete the entities (sheets and chords will be updated automatically)
  return mbImpl->delete_entities(&ents_to_delete[0], ents_to_delete.size());
}

MBErrorCode DualTool::delete_dual_entities(MBRange &entities) 
{
  MBErrorCode result = MB_SUCCESS;
  for (MBRange::iterator rit = entities.begin(); rit != entities.end(); rit++) {
    MBEntityHandle this_ent = *rit;
    MBErrorCode tmp_result = delete_dual_entities(&this_ent, 1);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }

  return result;
}

MBErrorCode DualTool::face_open_collapse(MBEntityHandle ocl, MBEntityHandle ocr,
                                         MBEntityHandle tcm) 
{
  MBErrorCode result;

    // gather data we can get just from looking at ocl, ocr, tcm
  MBEntityHandle zclf, zclb, zcrf, zcrb;
  MBEntityHandle tclu, tclm, tcll, tcru, tcrm, tcrl;
  MBEntityHandle thclu, thcll, thcmu, thcml, thcru, thcrl;
  MBEntityHandle sl, sm, sr, cl, cr;
  
  result = foc_gather_data(ocl, ocr, tcm,
                           zclf, zclb, zcrf, zcrb,
                           tclu, tclm, tcll, tcru, tcrm, tcrl,
                           thclu, thcll, thcmu, thcml, thcru, thcrl,
                           sl, sm, sr, cl, cr); RR;

  MBEntityHandle new_ocb, new_ocf, new_cb, new_cf;
  
  result = foc_1cells(zclf, zclb, ocl, cl,
                      zcrf, zcrb, ocr, cr,
                      sm, sr,
                      new_ocb, new_ocf, new_cb, new_cf);
  
  return MB_SUCCESS;
}

MBErrorCode DualTool::foc_1cells(MBEntityHandle zclf, MBEntityHandle zclb, 
                                 MBEntityHandle ocl, MBEntityHandle cl,
                                 MBEntityHandle zcrf, MBEntityHandle zcrb, 
                                 MBEntityHandle ocr, MBEntityHandle cr,
                                 MBEntityHandle sm, MBEntityHandle sr,
                                 MBEntityHandle &new_ocb, MBEntityHandle &new_ocf,
                                 MBEntityHandle &new_cb, MBEntityHandle &new_cf) 
{
  std::vector<MBEntityHandle> cl_1cells, cl_split_1cells, cr_1cells, cr_split_1cells;
  MBErrorCode result;

    // break the chords; make them go in opposite directions, so it's easier
    // to join them up
  result = foc_break_chord(cl, ocl, zclb, cl_1cells, cl_split_1cells); RR;
  result = foc_break_chord(cr, ocr, zcrf, cr_1cells, cr_split_1cells); RR;

    // make new 1cells
  MBEntityHandle new_verts[2];
  new_verts[0] = zclb; new_verts[1] = zcrb; 
  result = mbImpl->create_element(MBEDGE, new_verts, 2, new_ocb); RR;
  new_verts[0] = zclf; new_verts[1] = zcrf; 
  result = mbImpl->create_element(MBEDGE, new_verts, 2, new_ocf); RR;
  
    // construct new chord 1cell lists
  std::vector<MBEntityHandle> cb_1cells, cf_1cells;

    // if cr was blind, insert reversed cl split and ocf onto front of cb
  if (cr_split_1cells.empty() && !cl_split_1cells.empty()) {
    std::copy(cl_split_1cells.rbegin(), cl_split_1cells.rend(), std::back_inserter(cb_1cells));
    cb_1cells.push_back(new_ocf);
  }

    // reverse cr & put on cb
  std::copy(cr_1cells.rbegin(), cr_1cells.rend(), std::back_inserter(cb_1cells));
    // add new_ocb and append cl
  cb_1cells.push_back(new_ocb);
  std::copy(cl_1cells.begin(), cl_1cells.end(), std::back_inserter(cb_1cells));
  
    // if both chords got split, add the pieces to a new list
  if (!cl_split_1cells.empty() && !cr_split_1cells.empty()) {
    std::copy(cl_split_1cells.rbegin(), cl_split_1cells.rend(), std::back_inserter(cf_1cells));
    cf_1cells.push_back(new_ocf);
    std::copy(cr_split_1cells.begin(), cr_split_1cells.end(), std::back_inserter(cf_1cells));
  }
    // if niether was split, just need to add ocf
  else if (cl_split_1cells.empty() && cr_split_1cells.empty()) {
    cb_1cells.push_back(new_ocf);
  }
    // if cl was blind, add ocf and remainder of cr split to cb
  else if (cl_split_1cells.empty()) {
    cb_1cells.push_back(new_ocf);
    std::copy(cr_split_1cells.begin(), cr_split_1cells.end(), std::back_inserter(cb_1cells));
  }
    
    // make new chords
    // cb
  result = mbImpl->create_meshset((MESHSET_ORDERED | MESHSET_TRACK_OWNER), new_cb); RR;

    // give them parents sr and sm
  result = mbImpl->add_parent_child(sr, new_cb); RR;
  result = mbImpl->add_parent_child(sm, new_cb); RR;
  
    // now add 1cells
  result = mbImpl->add_entities(new_cb, &cb_1cells[0], cb_1cells.size()); RR;

    // cf
  result = mbImpl->create_meshset((MESHSET_ORDERED | MESHSET_TRACK_OWNER), new_cf); RR;

    // give them parents sr and sm
  result = mbImpl->add_parent_child(sr, new_cf); RR;
  result = mbImpl->add_parent_child(sm, new_cf); RR;
  
    // now add 1cells
  result = mbImpl->add_entities(new_cf, &cf_1cells[0], cf_1cells.size()); RR;

    // now delete chords cl, cr
  result = mbImpl->delete_entities(&cl, 1); RR;
  result = mbImpl->delete_entities(&cr, 1); RR;
  
  return MB_SUCCESS;
}

MBErrorCode DualTool::foc_break_chord(MBEntityHandle chord,
                                      MBEntityHandle first_1cell,
                                      MBEntityHandle next_0cell,
                                      std::vector<MBEntityHandle> &chord_1cells,
                                      std::vector<MBEntityHandle> &new_chord_1cells) 
{
  std::vector<MBEntityHandle> tmp_1cells;
  MBErrorCode result;
  result = mbImpl->get_entities_by_handle(chord, tmp_1cells); RR;
  std::vector<MBEntityHandle>::iterator vit;
  
    // position list at first_1cell
  vit = std::find(tmp_1cells.begin(), tmp_1cells.end(), first_1cell);
  if (vit == tmp_1cells.end()) return MB_FAILURE;
  int index = vit - tmp_1cells.end();

    // get common vtx with next 1cell to see if we're forward or reverse
  MeshTopoUtil mtu(mbImpl);
  MBEntityHandle next_v = 
    mtu.common_entity(first_1cell, tmp_1cells[(index+1)%tmp_1cells.size()], 0);
  int direction = (next_0cell == next_v ? 1 : -1);

    // start copy one past first_1cell, so first_1cell isn't on chord
  int i = (index + direction) % tmp_1cells.size();
  for (; i != (int)tmp_1cells.size() && i != -1; i += direction)
    chord_1cells.push_back(tmp_1cells[i]);
  
  i = (i + tmp_1cells.size()) % tmp_1cells.size();
  if (is_blind(chord)) {
    // now get the others
    for (; i != index; i += direction)
      chord_1cells.push_back(tmp_1cells[i]);
  }
  else {
      // put remaining 1cells on new list, but going away from first_1cell
    i = (index - direction + tmp_1cells.size()) % tmp_1cells.size();
    for (; i != (int)tmp_1cells.size() && i != -1; i -= direction)
      new_chord_1cells.push_back(tmp_1cells[i]);
  }

    // done
  return MB_SUCCESS;
}
  
MBErrorCode DualTool::foc_gather_data(const MBEntityHandle ocl, const MBEntityHandle ocr, 
                                      const MBEntityHandle tcm, 
                                        // 0-cells, left & right
                                      MBEntityHandle &zclf, MBEntityHandle &zclb, 
                                      MBEntityHandle &zcrf, MBEntityHandle &zcrb,
                                        // 2-cells, left & right
                                      MBEntityHandle &tclu, MBEntityHandle &tclm, MBEntityHandle &tcll, 
                                      MBEntityHandle &tcru, MBEntityHandle &tcrm, MBEntityHandle &tcrl,
                                        // 3-cells, left & right
                                      MBEntityHandle &thclu, MBEntityHandle &thcll, MBEntityHandle &thcmu, 
                                      MBEntityHandle &thcml, MBEntityHandle &thcru, MBEntityHandle &thcrl,
                                        // sheets
                                      MBEntityHandle &sl, MBEntityHandle &sm, MBEntityHandle &sr,
                                        // chords
                                      MBEntityHandle &cl, MBEntityHandle &cr) 
{
    // for explanation of algorithm, see notes 7/13/05
  MBErrorCode result;
  
    // get vertices around tcm
  const MBEntityHandle *tcm_verts, *ocl_verts, *ocr_verts;
  int tcm_verts_size, oc_verts_size;
  result = mbImpl->get_connectivity(tcm, tcm_verts, tcm_verts_size); RR;
  result = mbImpl->get_connectivity(ocl, ocl_verts, oc_verts_size); RR;
  assert(2 == oc_verts_size);
  result = mbImpl->get_connectivity(ocr, ocr_verts, oc_verts_size); RR;
  assert(2 == oc_verts_size);
  
  int side_no, sense, offset;
    // check ordering such that zclb comes before zclf in tcm's vert list; if
    // ocl and tcm are same sense, then zclb,zclf are in order on ocl
  result = mbImpl->side_number(tcm, ocl, side_no, sense, offset); RR;
  zclb = ocl_verts[(1-sense)/2];
  zclf = ocl_verts[(1+sense)/2];

    // reversed for ocr
  result = mbImpl->side_number(tcm, ocr, side_no, sense, offset); RR;
  zcrb = ocl_verts[(1+sense)/2];
  zcrf = ocl_verts[(1-sense)/2];

    // thcmu, thcml
  MBRange dum_range, dum_range_2;
  result = mbImpl->get_adjacencies(&tcm, 1, 3, false, dum_range); RR;
  assert(2 == dum_range.size());
  thcmu = *dum_range.begin();
  thcml = *dum_range.rbegin();

  result = foc_get_neighbor_23cells(ocl, tcm, thcmu, thcml,
                                    tclu, tclm, tcll,
                                    thclu, thcll); RR;

  result = foc_get_neighbor_23cells(ocr, tcm, thcmu, thcml,
                                    tcru, tcrm, tcrl,
                                    thcru, thcrl); RR;

    // sm, sl, sr
  sm = get_dual_hyperplane(tcm);
  if (0 == sm) {
    std::cerr << "Couldn't get dual surface for tcm." << std::endl;
    return MB_FAILURE;
  }
  
  sl = get_dual_hyperplane(tclu);
  if (0 == sl) {
    std::cerr << "Couldn't get dual surface for tclu." << std::endl;
    return MB_FAILURE;
  }
  
  sr = get_dual_hyperplane(tcru);
  if (0 == sl) {
    std::cerr << "Couldn't get dual surface for tcru." << std::endl;
    return MB_FAILURE;
  }
  
    // cl, cr
  cl = get_dual_hyperplane(ocl);
  if (0 == sm) {
    std::cerr << "Couldn't get dual curve for ocl." << std::endl;
    return MB_FAILURE;
  }
  
    // cl, cr
  cr = get_dual_hyperplane(ocr);
  if (0 == sm) {
    std::cerr << "Couldn't get dual curve for ocr." << std::endl;
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode DualTool::foc_get_neighbor_23cells(const MBEntityHandle oc,
                                               const MBEntityHandle tcm,
                                               const MBEntityHandle thcmu,
                                               const MBEntityHandle thcml,
                                               MBEntityHandle &tcxu, 
                                               MBEntityHandle &tcxm, 
                                               MBEntityHandle &tcxl, 
                                               MBEntityHandle &thcxu, 
                                               MBEntityHandle &thcxl) 
{
    // given a 1cell, 2cell, and an upper and lower 3cell (in that order),
    // return the adjacent upper, middle and lower 2cells sharing the 1cell
    // and adjacent to the 3cells (respectively), and the adjacent 3cells
  MBErrorCode result;
  MBRange dum_range, dum_range_2;
  
    // tcu, tcl, tcm
  dum_range.insert(tcm);
  dum_range.insert(thcmu);
  dum_range.insert(oc);
  result = mbImpl->get_adjacencies(dum_range, 2, false, dum_range_2); RR;
  assert(dum_range_2.size() == 2 &&
         (*dum_range_2.begin() == tcm || *dum_range_2.rbegin() == tcm));
  if (*dum_range_2.begin() == tcm)
    tcxu = *dum_range_2.rbegin();
  else
    tcxu = *dum_range_2.begin();

  dum_range_2.clear();
  dum_range.erase(thcmu);
  dum_range.insert(thcml);
  result = mbImpl->get_adjacencies(dum_range, 2, false, dum_range_2); RR;
  assert(dum_range_2.size() == 2 &&
         (*dum_range_2.begin() == tcm || *dum_range_2.rbegin() == tcm));
  if (*dum_range_2.begin() == tcm)
    tcxl = *dum_range_2.rbegin();
  else
    tcxl = *dum_range_2.begin();
  
  dum_range.clear(); dum_range_2.clear();
  dum_range.insert(oc);
  result = mbImpl->get_adjacencies(dum_range, 2, false, dum_range_2); RR;
  assert(4 == dum_range_2.size());
  dum_range_2.erase(tcm);
  dum_range_2.erase(tcxu);
  dum_range_2.erase(tcxl);
  tcxm = *dum_range_2.begin();
  
    // thcu, thcl
  dum_range.clear(); dum_range_2.clear();
  dum_range.insert(tcxm);
  dum_range.insert(tcxu);
  result = mbImpl->get_adjacencies(dum_range, 3, false, dum_range_2); RR;
  assert(1 == dum_range_2.size());
  thcxu = *dum_range_2.begin();
  dum_range_2.clear();
  dum_range.erase(tcxu);
  dum_range.insert(tcxl);
  result = mbImpl->get_adjacencies(dum_range, 3, false, dum_range_2); RR;
  assert(1 == dum_range_2.size());
  thcxl = *dum_range_2.begin();
  
  return MB_SUCCESS;
}

//! returns true if first & last vertices are dual to hexes (not faces)
bool DualTool::is_blind(const MBEntityHandle chord) 
{
    // must be an entity set
  if (TYPE_FROM_HANDLE(chord) != MBENTITYSET) return false;
  
    // get the vertices
  MBRange verts, ents;
  MBErrorCode result = mbImpl->get_entities_by_handle(chord, ents); 
  if (MB_SUCCESS != result || ents.empty()) return false;
  
  result = mbImpl->get_adjacencies(ents, 0, false, verts, MBInterface::UNION);
  if (MB_SUCCESS != result || verts.empty()) return false;
  
  for (MBRange::iterator rit = verts.begin(); rit != verts.end(); rit++) {
      // get dual entity for this vertex
    MBEntityHandle dual_ent = get_dual_entity(*rit);
    if (0 == dual_ent) continue;
    if (TYPE_FROM_HANDLE(dual_ent) == MBQUAD) return false;
  }

    // if none of the vertices' duals were quads, chord must be blind
  return true;
}

  //! given a 1-cell and a chord, return the neighboring vertices on the
  //! chord, in the same order as the 1-cell's vertices
MBErrorCode DualTool::get_opposite_verts(const MBEntityHandle middle_edge, 
                                         const MBEntityHandle chord, 
                                         MBEntityHandle *verts) 
{
    // get the edges on the chord, in order, and move to middle_edge
  std::vector<MBEntityHandle> chord_edges;
  const MBEntityHandle *connect;
  int num_connect;
  
  MBErrorCode result = mbImpl->get_entities_by_handle(chord, chord_edges); RR;
  std::vector<MBEntityHandle>::iterator vit = std::find(chord_edges.begin(), chord_edges.end(),
                                                        middle_edge);
  result = mbImpl->get_connectivity(middle_edge, connect, num_connect); RR;

  if (
      // middle_edge isn't on this chord
    vit == chord_edges.end() || 
      // chord only has 1 edge
    chord_edges.size() == 1 ||
      // middle_edge is at beginning or end and chord isn't blind
      ((vit == chord_edges.begin() || vit == chord_edges.end()-1) && !is_blind(chord))) 
    return MB_FAILURE;

  else if (chord_edges.size() == 2) {
      // easier if it's a 2-edge blind chord, just get vertices in opposite order
    verts[0] = connect[1];
    verts[1] = connect[0];
    return MB_SUCCESS;
  }
  
    // get vertices with the prev edge & subtract vertices of 1-cell
  if (vit == chord_edges.begin())
    vit = chord_edges.end() - 1;
  else vit--;
  MBRange dum_connect, middle_connect;
  result = mbImpl->get_connectivity(&middle_edge, 1, middle_connect); RR;
  result = mbImpl->get_connectivity(&(*vit), 1, dum_connect); RR;
  dum_connect = dum_connect.subtract(middle_connect);
  assert(dum_connect.size() == 1);

    // put in verts[0]
  verts[0] = *dum_connect.begin();

    // same with prev edge
  vit++;
  if (vit == chord_edges.end()) vit = chord_edges.begin();
  else vit++;
  dum_connect.clear();
  result = mbImpl->get_connectivity(&(*vit), 1, dum_connect); RR;
  dum_connect = dum_connect.subtract(middle_connect);
  assert(dum_connect.size() == 1);

    // put in verts[1]
  verts[1] = *dum_connect.begin();

    // if verts[0] and 1st vertex of 1cell don't have common edge, switch verts
  MeshTopoUtil mtu(mbImpl);
  if (0 == mtu.common_entity(verts[0], connect[0], 1)) {
    MBEntityHandle dum_h = verts[0];
    verts[0] = verts[1];
    verts[1] = dum_h;
  }
  
  assert(0 != mtu.common_entity(verts[0], connect[0], 1));
  
  return MB_SUCCESS;
}

MBErrorCode DualTool::get_dual_entities(const MBEntityHandle dual_ent,
                                        MBRange *dcells,
                                        MBRange *dedges,
                                        MBRange *dverts,
                                        MBRange *dverts_loop,
                                        MBRange *dedges_loop)
{
  MBErrorCode result;

  if (NULL != dcells) {
    result = mbImpl->get_entities_by_type(dual_ent, MBPOLYGON, *dcells);
    if (MB_SUCCESS != result) return result;
  }
  
  if (NULL != dedges) {
    if (NULL != dcells)
      result = mbImpl->get_adjacencies(*dcells, 1, false, *dedges, MBInterface::UNION);
    else 
      result = mbImpl->get_entities_by_type(dual_ent, MBEDGE, *dedges);

    if (MB_SUCCESS != result) return result;
  }

  if (NULL != dverts) {
    if (NULL != dcells)
      result = mbImpl->get_adjacencies(*dcells, 0, false, *dverts, MBInterface::UNION);
    else if (NULL != dedges)
      result = mbImpl->get_adjacencies(*dedges, 0, false, *dverts, MBInterface::UNION);
    else {
      MBRange all_ents;
      result = mbImpl->get_entities_by_handle(dual_ent, all_ents); RR;
      result = mbImpl->get_adjacencies(all_ents, 0, false, *dverts, MBInterface::UNION);
    }

    if (MB_SUCCESS != result) return result;
  }

  if (NULL != dverts_loop && NULL != dverts) {
    static std::vector<MBEntityHandle> dual_ents;
    dual_ents.reserve(dverts->size());
    result = mbImpl->tag_get_data(dualEntity_tag(), *dverts, &dual_ents[0]);
    if (MB_SUCCESS != result) return result;
    MBRange::iterator rit;
    unsigned int i;
    for (rit = dverts->begin(), i = 0; rit != dverts->end(); rit++, i++)
      if (0 != dual_ents[i] && mbImpl->type_from_handle(dual_ents[i]) == MBQUAD)
        dverts_loop->insert(*rit);
  }
  
  if (NULL != dedges_loop && NULL != dedges) {
    static std::vector<MBEntityHandle> dual_ents;
    dual_ents.reserve(dedges->size());
    result = mbImpl->tag_get_data(dualEntity_tag(), *dedges, &dual_ents[0]);
    if (MB_SUCCESS != result) return result;
    MBRange::iterator rit;
    unsigned int i;
    for (rit = dedges->begin(), i = 0; rit != dedges->end(); rit++, i++)
      if (0 != dual_ents[i] && mbImpl->type_from_handle(dual_ents[i]) == MBEDGE)
        dedges_loop->insert(*rit);
  }

  return result;
}

MBErrorCode DualTool::list_entities(const MBEntityHandle *entities,
                                    const int num_entities) const
{
  MBRange temp_range;
  MBErrorCode result;
  if (NULL == entities && num_entities <= 0) {
      // just list the numbers of each entity type
    int num_ents;
    std::cout << std::endl;
    std::cout << "Number of entities per type: " << std::endl;
    for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++) {
      result = mbImpl->get_number_entities_by_type(0, this_type, num_ents);
      std::cout << MBCN::EntityTypeName(this_type) << ": " << num_ents << std::endl;
    }
    std::cout << std::endl;

      // if negative num_entities, list the set hierarchy too
    if (0 > num_entities) {
      MBRange sets;
      result = mbImpl->get_entities_by_type(0, MBENTITYSET, sets);
      if (MB_SUCCESS != result) return result;
      for (MBRange::iterator rit = sets.begin(); rit != sets.end(); rit++) {
        mbImpl->list_entities(&(*rit), 1);
        result = mbImpl->get_number_entities_by_handle(*rit, num_ents);
        std::cout << "(" << num_ents << " total entities)" << std::endl;
      }
    }
    
    return MB_SUCCESS;
  }
      
  else if (NULL == entities) {

      // list all entities of all types
    std::cout << std::endl;
    for (MBEntityType this_type = MBVERTEX; this_type < MBMAXTYPE; this_type++) {
      temp_range.clear();
      std::cout << "Entities of type: " << MBCN::EntityTypeName(this_type) << ": " << std::endl;
      result = mbImpl->get_entities_by_type(0, this_type, temp_range);
      result = list_entities(temp_range);
    }
    std::cout << std::endl;
    
    return MB_SUCCESS;
  }

  else {
    std::copy(entities, entities+num_entities, mb_range_inserter(temp_range));
    return list_entities(temp_range);
  }
}
  
MBErrorCode DualTool::list_entities(const MBRange &entities) const
{
  MBErrorCode result;
  std::vector<MBEntityHandle> adj_vec;

    // list entities
  for (MBRange::const_iterator iter = entities.begin(); iter != entities.end(); iter++) {
    MBEntityType this_type = TYPE_FROM_HANDLE(*iter);
    std::cout << MBCN::EntityTypeName(this_type) << " " << mbImpl->id_from_handle(*iter) << ":" << std::endl;

    if (this_type == MBVERTEX) {
      double coords[3];
      result = mbImpl->get_coords(&(*iter), 1, coords);
      if (MB_SUCCESS != result) return result;
      std::cout << "Coordinates: (" << coords[0] << ", " << coords[1] << ", " << coords[2] 
           << ")" << std::endl;
    }
    else if (this_type == MBENTITYSET)
      mbImpl->list_entities(&(*iter), 1);

    MBEntityHandle dual_ent = get_dual_entity(*iter);
    if (0 != dual_ent) {
      std::cout << "Dual to " 
                << MBCN::EntityTypeName(mbImpl->type_from_handle(dual_ent)) << " "
                << mbImpl->id_from_handle(dual_ent) << std::endl;
    }
    
    std::cout << "  Adjacencies:" << std::endl;
    bool some = false;
    int multiple = 0;
    for (int dim = 0; dim <= 3; dim++) {
      if (dim == MBCN::Dimension(this_type)) continue;
      adj_vec.clear();
      result = mbImpl->get_adjacencies(&(*iter), 1, dim, false, adj_vec);
      if (MB_FAILURE == result) continue;
      for (std::vector<MBEntityHandle>::iterator adj_it = adj_vec.begin(); 
           adj_it != adj_vec.end(); adj_it++) {
        if (adj_it != adj_vec.begin()) std::cout << ", ";
        else std::cout << "   ";
        std::cout << MBCN::EntityTypeName(mbImpl->type_from_handle(*adj_it)) << " " 
                  << mbImpl->id_from_handle(*adj_it);
      }
      if (!adj_vec.empty()) {
        std::cout << std::endl;
        some = true;
      }
      if (MB_MULTIPLE_ENTITIES_FOUND == result)
        multiple += dim;
    }
    if (!some) std::cout << "(none)" << std::endl;
    if (multiple != 0)
      std::cout << "   (MULTIPLE = " << multiple << ")" << std::endl;

    std::cout << std::endl;
  }

  return MB_SUCCESS;
}

