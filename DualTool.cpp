#include "DualTool.hpp"
#include "MeshTopoUtil.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"
#include <string>
#include <iostream>
#include <sstream>

bool debug = false;

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

      // tag the edge indicating it's a dual entity
    tmp_result = mbImpl->tag_set_data(isDualCell_tag(), &dual_ent, 1, &is_dual);
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
  if (MB_SUCCESS != result) return result;
  
  dual_ents.insert(dual_ents.end(), tmp_range.begin(), tmp_range.end());
  return MB_SUCCESS;
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
  const char *hp_tag_names[2] = {"dual_hyperplane_1", "dual_hyperplane_2"};
  MBTag hp_tag, gid_tag, mark_tag;
  int dum = -1;
  MBErrorCode result = mbImpl->tag_get_handle(hp_tag_names[dim-1], hp_tag);
  if (MB_SUCCESS != result) result = mbImpl->tag_create(hp_tag_names[dim-1], 
                                                        4, MB_TAG_SPARSE, hp_tag, &dum);
  if (MB_SUCCESS != result) return result;
  result = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, 4, 
                                                        MB_TAG_DENSE, gid_tag, &dum);
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
  int hp_val;
  int hp_ids = 0;
    // int report = tot_untreated.size() / 10;
  while (!tot_untreated.empty()) {
    if (debug && dim == 2 /*(tot_untreated.size()%report == 0)*/)
      std::cout << "Untreated list size " << tot_untreated.size() << "." << std::endl;
      
    this_ent = tot_untreated.back(); tot_untreated.pop_back();
    result = mbImpl->tag_get_data(hp_tag, &this_ent, 1, &hp_val);
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) continue;

      // test for this entity having a hyperplane assignment already
    else if (hp_val != -1)
      continue;
    
      // ok, doesn't have one; make a new hyperplane
    hp_val = hp_ids++;
    MBEntityHandle this_hp;
    result = mbImpl->create_meshset(MESHSET_SET, this_hp);
    if (MB_SUCCESS != result) continue;
    result = mbImpl->tag_set_data(gid_tag, &this_hp, 1, &hp_val);
    if (MB_SUCCESS != result) continue;
    result = mbImpl->tag_set_data(hp_tag, &this_hp, 1, &hp_val);
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) continue;
    hp_untreated.clear();

      // inner loop: traverse the hyperplane 'till we don't have any more
    MBRange tmp_star, star, tmp_range;
    while (0 != this_ent) {
      if (debug && hp_untreated.size()%10 == 0) 
        std::cout << "Dual surface " << hp_val << ", hp_untreated list size = " 
                  << hp_untreated.size() << "." << std::endl;
      
        // set the hp_val for this_ent
      result = mbImpl->tag_set_data(hp_tag, &this_ent, 1, &hp_val);
      if (MB_SUCCESS != result) continue;
      result = mbImpl->add_entities(this_hp, &this_ent, 1);
      if (MB_SUCCESS != result) continue;

        // get all neighbors connected to this entity
      tmp_range.clear(); tmp_star.clear(); star.clear();
      tmp_range.insert(this_ent);
      result = mbImpl->get_adjacencies(tmp_range, dim-1, true, tmp_star);
      if (MB_SUCCESS != result) continue;
      result = mbImpl->get_adjacencies(tmp_star, dim, false, star,
                                       MBInterface::UNION);
      if (MB_SUCCESS != result) continue;
      star.erase(this_ent);
      
        // for each star entity, see if it shares a cell with this one; if so,
        // it's not opposite

        // this_ent is already in tmp_range
      for (MBRange::iterator rit = star.begin(); rit != star.end(); rit++) {
          // check for tag first, 'cuz it's probably faster than checking adjacencies
        int r_val;
        result = mbImpl->tag_get_data(hp_tag, &(*rit), 1, &r_val);
        if (MB_SUCCESS == result && -1 != r_val) 
          continue;

        unsigned short tmp_mark;
        result = mbImpl->tag_get_data(mark_tag, &(*rit), 1, &tmp_mark);
        if (MB_SUCCESS == result && mark_val == tmp_mark) 
          continue;

          // passed tag test; check adjacencies
        tmp_star.clear();
        tmp_range.insert(*rit);
        result = mbImpl->get_adjacencies(tmp_range, dim+1, false, tmp_star);
        if (MB_SUCCESS != result) continue;
        
        if (tmp_star.empty()) {
            // have one on this hp; just put it on the hp_untreated list for now,
            // will get tagged and put in the hp set later
          hp_untreated.push_back(*rit);
          result = mbImpl->tag_set_data(mark_tag, &(*rit), 1, &mark_val);
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

  
