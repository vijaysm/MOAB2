#include "DualTool.hpp"
#include "MeshTopoUtil.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"

  //! tag name for dual surfaces
char *DualTool::DUAL_SURFACE_TAG_NAME = "DUAL_SURFACE";

  //! tag name for dual curves
char *DualTool::DUAL_CURVE_TAG_NAME = "DUAL_CURVE";

  //! tag name for dual cells
char *DualTool::IS_DUAL_CELL_TAG_NAME = "DUAL_CELL";

  //! tag name for dual entities
char *DualTool::DUAL_ENTITY_TAG_NAME = "DUAL_ENTITY";

DualTool::DualTool(MBInterface *impl) 
    : mbImpl(impl)
{
  MBTag dummy = 0;

  MBErrorCode result = mbImpl->tag_create(DUAL_SURFACE_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                                          dualSurfaceTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
  result = mbImpl->tag_create(DUAL_CURVE_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                              dualCurveTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);
  
  result = mbImpl->tag_create(IS_DUAL_CELL_TAG_NAME, 1, MB_TAG_BIT, 
                              isDualCellTag, &dummy);
  assert(MB_ALREADY_ALLOCATED == result || MB_SUCCESS == result);

  result = mbImpl->tag_create(DUAL_ENTITY_TAG_NAME, sizeof(MBEntityHandle), MB_TAG_DENSE, 
                              dualEntityTag, &dummy);
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

    // don't really need dual edges, but construct 'em anyway
  MBRange dual_edges;
  result = construct_dual_edges(all_faces, dual_edges);
  if (MB_SUCCESS != result || dual_edges.size() != all_faces.size()) return result;

    // construct dual faces
  MBRange dual_faces;
  result = construct_dual_faces(all_edges, dual_faces);
  if (MB_SUCCESS != result || dual_faces.size() != all_edges.size()) return result;

    // construct dual cells
  MBRange dual_cells;
  result = construct_dual_cells(all_vertices, dual_cells);
  if (MB_SUCCESS != result || dual_cells.size() != all_vertices.size()) return result;

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
  unsigned int is_dual = 0x1;
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
    
      // no dual entity; construct one; first need the avg coordinates
    double avg_pos[3];
    tmp_result = MeshTopoUtil(mbImpl).get_average_position(*rit, avg_pos);
    if (MB_SUCCESS != tmp_result) continue;
    
      // now construct the new dual entity
    tmp_result = mbImpl->create_vertex(avg_pos, dual_ent);
    if (MB_SUCCESS != tmp_result) continue;
    
      // tag it indicating it's a dual entity
    tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
    if (MB_SUCCESS != tmp_result) continue;
    
      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
      // tag the primal entity with its dual entity
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS != tmp_result) continue;
  }

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
    
    if (dual_verts.size() == 1) {
        // boundary face - make a dual vertex at the face center and put in list
      double avg_pos[3];      
      tmp_result = MeshTopoUtil(mbImpl).get_average_position(*rit, avg_pos);
      if (MB_SUCCESS != tmp_result) continue;
      
      tmp_result = mbImpl->create_vertex(avg_pos, dual_ent);
      if (MB_SUCCESS != tmp_result) continue;
    
        // tag it indicating it's a dual entity
      tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
      if (MB_SUCCESS != tmp_result) continue;

        // put it on vertex list
      dual_verts.push_back(dual_ent);
    }
    
    assert(dual_verts.size() == 2);
    
      // now create the dual edge
    tmp_result = mbImpl->create_element(MBEDGE, &dual_verts[0], 2, dual_ent);
    if (MB_SUCCESS != tmp_result || 0 == dual_ent) continue;
    
      // save it in the list of new dual ents
    dual_ents.insert(dual_ent);
    
      // tag the primal entity with its dual entity
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
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
    tmp_result = get_radial_dverts(*rit, rad_dverts);
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
    
      // tag the primal entity with its dual entity
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS != tmp_result) 
      continue;
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
    dfaces.reserve(edges.size());
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

      // tag the primal entity with its dual entity
    tmp_result = mbImpl->tag_set_data(dualEntity_tag(), &(*rit), 1, &dual_ent);
    if (MB_SUCCESS != tmp_result) continue;
  }

  return result;
}

  //! given an edge handle, return a list of dual vertices in radial order 
  //! around the edge
MBErrorCode DualTool::get_radial_dverts(const MBEntityHandle edge,
                                        std::vector<MBEntityHandle> &rad_dverts) 
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
  }
  
  return MB_SUCCESS;
}

/*  
  //! construct the dual entities for a hex mesh, including dual surfaces & curves
MBErrorCode DualTool::construct_hex_dual();
  
  //! get the sets representing dual surfaces (sheets) in hex dual
MBErrorCode DualTool::get_dual_surfaces(MBRange &dual_surfs);
  
  //! get the sets representing dual curves (chords) in hex dual
MBErrorCode DualTool::get_dual_curves(MBRange &dual_curves);
  
*/

  //! get the cells of the dual
MBErrorCode DualTool::get_dual_faces(MBRange &dual_faces) 
{
  if (0 == isDualCell_tag()) return MB_SUCCESS;

  unsigned int dum = 0x1;
  const void *dum_ptr = &dum;
  
  return
    mbImpl->get_entities_by_type_and_tag(0, MBPOLYGON, &isDualCellTag, &dum_ptr, 1,
                                         dual_faces);
}

  //! get the cells of the dual
MBErrorCode DualTool::get_dual_cells(MBRange &dual_cells) 
{
  if (0 == isDualCell_tag()) return MB_SUCCESS;
  
  unsigned int dum = 0x1;
  const void *dum_ptr = &dum;

  return
    mbImpl->get_entities_by_type_and_tag(0, MBPOLYHEDRON, &isDualCellTag, &dum_ptr, 1,
                                         dual_cells);
}


