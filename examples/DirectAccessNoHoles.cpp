/** @example DirectAccess.cpp \n
 * \brief Use direct access to MOAB data to avoid calling through API \n
 *
 * This example creates a 1d row of quad elements, with a user-specified number of "holes" (missing quads) in the row:
 *
 *  ----------------------      ----------------------      --------
 *  |      |      |      |      |      |      |      |      |      |       
 *  |      |      |      |(hole)|      |      |      |(hole)|      | ...
 *  |      |      |      |      |      |      |      |      |      |
 *  ----------------------      ----------------------      --------
 *
 * This makes (nholes+1) contiguous runs of quad handles in the handle space
 * This example shows how to use the xxx_iterate functions in MOAB (xxx = coords, connect, tag, adjacencies) to get 
 * direct pointer access to MOAB internal storage, which can be used without calling through the MOAB API.
 *
 *    -#  Initialize MOAB \n
 *    -#  Create a quad mesh with holes, as depicted above
 *    -#  Create 2 dense tags (tag1, tag2) for avg position to assign to quads, and # verts per quad (tag3)
 *    -#  Get connectivity, coordinate, tag1 iterators
 *    -#  Iterate through quads, computing midpoint based on vertex positions, set on quad-based tag1
 *    -#  Set up map from starting quad handle for a chunk to struct of (tag1_ptr, tag2_ptr, tag3_ptr), pointers to
 *        the dense tag storage for those tags for the chunk
 *    -#  Iterate through vertices, summing positions into tag2 on connected quads and incrementing vertex count
 *    -#  Iterate through quads, normalizing tag2 by vertex count and comparing values of tag1 and tag2
 *
 * <b>To compile</b>: \n
 *    make DirectAccessWithHoles MOAB_DIR=<installdir>  \n
 * <b>To run</b>: ./DirectAccess [-nquads <# quads>] [-holes <# holes>]\n
 *
 */

#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ReadUtilIface.hpp"
#include <map>
#include <iostream>
#include <assert.h>

// Error routines for use with MOAB API
#define CHKERR(CODE, MSG) do {if (MB_SUCCESS != (CODE)) {std::string errstr;  mbImpl->get_last_error(errstr);   \
        std::cerr << errstr << std::endl; std::cerr << MSG << std::endl; return CODE;}} while(false)

using namespace moab;

ErrorCode create_mesh_no_holes(Interface *mbImpl, int nquads);

int main(int argc, char **argv)
{
    // get MOAB instance
  Interface *mbImpl = new Core;
  if (NULL == mbImpl) return 1;

  int nquads = 1000;
  
    // parse options
  ProgOptions opts;
  opts.addOpt<int>(std::string("nquads,n"), std::string("Number of quads in the mesh (default = 1000"), &nquads);
  opts.parseCommandLine(argc, argv);

    // create simple structured mesh with hole, but using unstructured representation
  ErrorCode rval = create_mesh_no_holes(mbImpl, nquads); CHKERR(rval, "Trouble creating mesh.");
  
    // get all vertices and non-vertex entities
  Range verts, quads;
  rval = mbImpl->get_entities_by_handle(0, quads); CHKERR(rval, "Getting all entities.");
  verts = quads.subset_by_dimension(0);
  quads -= verts;

    // create tag1 (element-based avg), tag2 (vertex-based avg), tag3 (# connected verts)
  Tag tag1, tag2, tag3;
  rval = mbImpl->tag_get_handle("tag1", 3, MB_TYPE_DOUBLE, tag1, MB_TAG_DENSE | MB_TAG_CREAT); CHKERR(rval, "Creating tag1.");
  double def_val[3] = {0.0, 0.0, 0.0}; // need a default value for tag2 because we sum into it
  rval = mbImpl->tag_get_handle("tag2", 3, MB_TYPE_DOUBLE, tag2, MB_TAG_DENSE | MB_TAG_CREAT, def_val); CHKERR(rval, "Creating tag2.");
  int def_val_int = 0;  // need a default value for tag3 because we increment it
  rval = mbImpl->tag_get_handle("tag3", 1, MB_TYPE_INTEGER, tag3, MB_TAG_DENSE | MB_TAG_CREAT, &def_val_int); CHKERR(rval, "Creating tag3.");
  
    // Get pointers to connectivity, coordinate, tag, and adjacency arrays; each of these returns a count,
    // which should be compared to the # entities you expect to verify there's only one chunk (no holes)
  int count, vpere;
  EntityHandle *conn_ptr;
  rval = mbImpl->connect_iterate(quads.begin(), quads.end(), conn_ptr, vpere, count); CHKERR(rval, "connect_iterate.");
  assert(count == (int) quads.size()); // should end up with just one contiguous chunk of quads

  double *x_ptr, *y_ptr, *z_ptr;
  rval = mbImpl->coords_iterate(verts.begin(), verts.end(), x_ptr, y_ptr, z_ptr, count); CHKERR(rval, "coords_iterate.");
  assert(count == (int) verts.size()); // should end up with just one contiguous chunk of vertices

  double *tag1_ptr, *tag2_ptr;
  int *tag3_ptr;
  rval = mbImpl->tag_iterate(tag1, quads.begin(), quads.end(), count, (void*&)tag1_ptr); CHKERR(rval, "tag1_iterate.");
  assert(count == (int) quads.size()); // should end up with just one contiguous chunk of quads
  rval = mbImpl->tag_iterate(tag2, quads.begin(), quads.end(), count, (void*&)tag2_ptr); CHKERR(rval, "tag2_iterate.");
  assert(count == (int) quads.size()); // should end up with just one contiguous chunk of quads
  rval = mbImpl->tag_iterate(tag3, quads.begin(), quads.end(), count, (void*&)tag3_ptr); CHKERR(rval, "tag3_iterate.");
  assert(count == (int) quads.size()); // should end up with just one contiguous chunk of quads

  const std::vector<EntityHandle> **adjs_ptr;
  rval = mbImpl->adjacencies_iterate(verts.begin(), verts.end(), adjs_ptr, count); CHKERR(rval, "adjacencies_iterate.");
  assert(count == (int) verts.size()); // should end up with just one contiguous chunk of vertices
  
    // Start_ handles used to compute indices into vertex/quad arrays; can use direct subtraction because we know
    // there aren't any holes in the handle spaces for verts or quads
  EntityHandle start_vert = *verts.begin(), start_quad = *quads.begin();

    // iterate over elements, computing tag1 from coords positions
  for (int i = 0; i < nquads; i++) {
    tag1_ptr[3*i+0] = tag1_ptr[3*i+1] = tag1_ptr[3*i+2] = 0.0; // initialize position for this element
    for (int j = 0; j < vpere; j++) { // loop over vertices in this element
      int v_index = conn_ptr[vpere*i+j] - start_vert; // vert index is just the offset from start vertex
      tag1_ptr[3*i+0] += x_ptr[v_index];
      tag1_ptr[3*i+1] += y_ptr[v_index]; // sum vertex positions into tag1...
      tag1_ptr[3*i+2] += z_ptr[v_index];
    }
    tag1_ptr[3*i+0] /= vpere;
    tag1_ptr[3*i+1] /= vpere; // then normalize
    tag1_ptr[3*i+2] /= vpere;
  } // loop over elements in chunk
    
    // Iterate through vertices, summing positions into tag2 on connected elements and incrementing vertex count
  for (int v = 0; v < count; v++) {
    const std::vector<EntityHandle> *avec = *(adjs_ptr+v);
    for (std::vector<EntityHandle>::const_iterator ait = avec->begin(); ait != avec->end(); ait++) {
        // *ait is the quad handle, its index is computed by subtracting the start quad handle
      int a_ind = *ait - start_quad;
      tag2_ptr[3*a_ind+0] += x_ptr[v];   // tag on each element is 3 doubles, x/y/z
      tag2_ptr[3*a_ind+1] += y_ptr[v];
      tag2_ptr[3*a_ind+2] += z_ptr[v];
      tag3_ptr[a_ind]++;  // increment the vertex count
    }
  }
        
    // Normalize tag2 by vertex count (tag3); loop over elements using same approach as before
    // At the same time, compare values of tag1 and tag2
  for (Range::iterator q_it = quads.begin(); q_it != quads.end(); q_it++) {
    int i = *q_it - start_quad;
    for (int j = 0; j < 3; j++) tag2_ptr[3*i+j] /= (double)tag3_ptr[i];  // normalize by # verts
    if (tag1_ptr[3*i] != tag2_ptr[3*i] || tag1_ptr[3*i+1] != tag2_ptr[3*i+1] || tag1_ptr[3*i+2] != tag2_ptr[3*i+2]) 
      std::cout << "Tag1, tag2 disagree for element " << *q_it + i << std::endl;
  }

    // Ok, we're done, shut down MOAB
  delete mbImpl;

  return 0;
}

ErrorCode create_mesh_no_holes(Interface *mbImpl, int nquads) 
{
    // first make the mesh, a 1d array of quads with left hand side x = elem_num; vertices are numbered in layers
  ReadUtilIface *read_iface;
  ErrorCode rval = mbImpl->query_interface(read_iface); CHKERR(rval, "query_interface");
  std::vector<double *> coords;
  EntityHandle start_vert, start_elem, *connect;
    // create verts, num is 4(nquads+1) because they're in a 1d row; will initialize coords in loop over quads later
  rval = read_iface->get_node_coords (3, 2*(nquads+1), 0, start_vert, coords); CHKERR(rval, "get_node_arrays");
    // create quads
  rval = read_iface->get_element_connect(nquads, 4, MBQUAD, 0, start_elem, connect); CHKERR(rval, "get_element_connect");
  for (int i = 0; i < nquads; i++) {
    coords[0][2*i] = coords[0][2*i+1] = (double) i; // x values are all i
    coords[1][2*i] = 0.0; coords[1][2*i+1] = 1.0; // y coords
    coords[2][2*i] = coords[2][2*i+1] = (double) 0.0; // z values, all zero (2d mesh)
    EntityHandle quad_v = start_vert + 2*i;
    for (int j = 0; j < 4; j++) connect[4*i+j] = quad_v+j; // connectivity of each quad is a sequence starting from quad_v
  }
    // last two vertices
  coords[0][2*nquads] = coords[0][2*nquads+1] = (double) nquads;
  coords[1][2*nquads] = 0.0; coords[1][2*nquads+1] = 1.0; // y coords
  coords[2][2*nquads] = coords[2][2*nquads+1] = (double) 0.0; // z values, all zero (2d mesh)
  
    // call a vertex-quad adjacencies function to generate vertex-element adjacencies in MOAB
  Range dum_range;
  rval = mbImpl->get_adjacencies(&start_vert, 1, 2, false, dum_range); CHKERR(rval, "get_adjacencies");
  assert(!dum_range.empty());

  return MB_SUCCESS;
}

    
  
