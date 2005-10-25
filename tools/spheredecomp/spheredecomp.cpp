#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MeshTopoUtil.hpp"
#include "MBCN.hpp"
#include <iostream>

MBInterface *gMB;
MBTag subdiv_vertices_tag = 0;
const char *SUBDIV_VERTICES_TAG_NAME = "SUBDIV_VERTICES";
MBErrorCode compute_nodes(const int dim);
MBErrorCode build_hexes(std::vector<MBEntityHandle> &sphere_hexes,
                        std::vector<MBEntityHandle> &interstic_hexes);
MBErrorCode retrieve_subdiv_verts(MBEntityHandle tet, MBEntityHandle this_ent,
                                  const MBEntityHandle *tet_conn,
                                  const int dim, MBEntityHandle *subdiv_verts);
MBErrorCode subdivide_tet(MBEntityHandle tet, 
                          std::vector<MBEntityHandle> &sphere_hexes,
                          std::vector<MBEntityHandle> &interstic_hexes);


#define RR if (MB_SUCCESS != result) return result

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <input_mesh> <output_mesh>" << std::endl;
    return 0;
  }
  
    // create MOAB
  gMB = new MBCore();
  
    // read in mesh
  MBErrorCode result = gMB->load_mesh(argv[1]); RR;

    // need to make sure all interior edges and faces are created
  MBRange all_verts;
  result = gMB->get_entities_by_type(0, MBVERTEX, all_verts); RR;
  MeshTopoUtil mtu(gMB);
  result = mtu.construct_aentities(all_verts);
  
    // create tag to hold vertices
  result = gMB->tag_create(SUBDIV_VERTICES_TAG_NAME, 9*sizeof(MBEntityHandle), 
                           MB_TAG_DENSE, MB_TYPE_HANDLE, subdiv_vertices_tag, NULL); RR;

    // compute nodal positions for each dimension element
  result = compute_nodes(1); RR;
  result = compute_nodes(2); RR;
  result = compute_nodes(3); RR;
  
    // build hex elements
  std::vector<MBEntityHandle> sphere_hexes, interstic_hexes;
  result = build_hexes(sphere_hexes, interstic_hexes); RR;
  
    // write mesh
  result = gMB->write_mesh(argv[2]); RR;
  
  return 0;
}

MBErrorCode compute_nodes(const int dim) 
{
    // get facets of that dimension
  MBRange these_ents;
  const MBEntityType the_types[4] = {MBVERTEX, MBEDGE, MBTRI, MBTET};
  
  MBErrorCode result = gMB->get_entities_by_dimension(0, dim, these_ents); RR;
  assert(gMB->type_from_handle(*these_ents.begin()) == the_types[dim] &&
         gMB->type_from_handle(*these_ents.rbegin()) == the_types[dim]);
  
  MBEntityHandle subdiv_vertices[9];
  MeshTopoUtil mtu(gMB);
  double avg_pos[3], vert_pos[12], new_vert_pos[12], new_new_vert_pos[3];
  double radii[4];
  int num_verts = MBCN::VerticesPerEntity(the_types[dim]);
  
  for (MBRange::iterator rit = these_ents.begin(); rit != these_ents.end(); rit++) {
    
      // get vertices
    const MBEntityHandle *connect;
    int num_connect;
    result = gMB->get_connectivity(*rit, connect, num_connect); RR;

      // compute center
    result = mtu.get_average_position(connect, num_connect, avg_pos); RR;

      // create center vertex
    result = gMB->create_vertex(avg_pos, subdiv_vertices[num_verts]); RR;
    
      // get coords of other vertices
    result = gMB->get_coords(connect, num_connect, vert_pos); RR;
    
      // get radii associated with each vertex
      // for now, just set to constant
    for (int i = 0; i < num_verts; i++) radii[i] = 0.1;
    
      // compute subdiv vertex position for each vertex
    for (int i = 0; i < num_verts; i++) {
      for (int j = 0; j < 3; j++)
        new_vert_pos[3*i+j] = vert_pos[3*i+j] + radii[i] * (avg_pos[j] - vert_pos[3*i+j]);

      // create vertex at this position
      result = gMB->create_vertex(&new_vert_pos[3*i], subdiv_vertices[i]); RR;
    }
    
      // compute subdiv vertex positions for vertices inside spheres; just mid-pt between
      // previous subdiv vertex and corner vertex
    for (int i = 0; i < num_verts; i++) {
      for (int j = 0; j < 3; j++) 
        new_new_vert_pos[j] = .5 * (vert_pos[3*i+j] + new_vert_pos[3*i+j]);

      result = gMB->create_vertex(new_new_vert_pos, subdiv_vertices[num_verts+1+i]);
    }

      // set the tag
    result = gMB->tag_set_data(subdiv_vertices_tag, &(*rit), 1, subdiv_vertices); RR;
  }
  
  return result;
}

MBErrorCode build_hexes(std::vector<MBEntityHandle> &sphere_hexes,
                        std::vector<MBEntityHandle> &interstic_hexes) 
{
    // build hexes inside each tet element separately
  MBRange tets;
  MBErrorCode result = gMB->get_entities_by_type(0, MBTET, tets); RR;
  
  for (MBRange::iterator vit = tets.begin(); vit != tets.end(); vit++) {
    result = subdivide_tet(*vit, sphere_hexes, interstic_hexes); RR;
  }
  
  return MB_SUCCESS;
}

MBErrorCode subdivide_tet(MBEntityHandle tet, 
                          std::vector<MBEntityHandle> &sphere_hexes,
                          std::vector<MBEntityHandle> &interstic_hexes) 
{
    // 99: (#subdiv_verts/entity=9) * (#edges=6 + #faces=4 + 1=tet)
  MBEntityHandle subdiv_verts[99];

    // get tet connectivity
  std::vector<MBEntityHandle> tet_conn;
  MBErrorCode result = gMB->get_connectivity(&tet, 1, tet_conn); RR;
  
  for (int dim = 1; dim <= 3; dim++) {
      // get entities of this dimension
    std::vector<MBEntityHandle> ents;
    if (dim != 3) {
      result = gMB->get_adjacencies(&tet, 1, dim, false, ents); RR; 
    }
    else ents.push_back(tet);
    
      // for each, get subdiv verts & put into vector
    for (std::vector<MBEntityHandle>::iterator vit = ents.begin(); vit != ents.end(); vit++) {
      result = retrieve_subdiv_verts(tet, *vit, &tet_conn[0], dim, subdiv_verts); RR;
    }
  }

    // ok, subdiv_verts are in canonical order; now create the hexes, using pre-computed templates
    // first, interstices hexes
#define EDGE 0
#define FACE 1
#define TET 2
#define AINDEX 0
#define BINDEX 1
#define CINDEX 2
#define DINDEX 3
#define EINDEX 4
#define FINDEX 5
#define GINDEX 6
#define HINDEX 7
#define IINDEX 8
#define V0INDEX 0
#define V1INDEX 1
#define V2INDEX 2
#define V3INDEX 3
#define CV(a) tet_conn[a]
#define ESV(a,b) subdiv_verts[a*9+b]
#define FSV(a,b) subdiv_verts[54+a*9+b]
#define TSV(a,b) subdiv_verts[90+a*9+b]

  MBEntityHandle this_connect[8], this_hex;

// V0:
  int i = 0;
  this_connect[i++]=ESV(0,AINDEX); this_connect[i++]=ESV(0,CINDEX); this_connect[i++]=FSV(3,DINDEX); this_connect[i++]=FSV(3,AINDEX); 
  this_connect[i++]=FSV(0,AINDEX); this_connect[i++]=FSV(0,DINDEX); this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=TSV(0,AINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(0,AINDEX); this_connect[i++]=FSV(0,DINDEX); this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=TSV(0,AINDEX); 
  this_connect[i++]=ESV(3,AINDEX); this_connect[i++]=ESV(3,CINDEX); this_connect[i++]=FSV(2,DINDEX); this_connect[i++]=FSV(2,AINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(3,AINDEX); this_connect[i++]=FSV(3,DINDEX); this_connect[i++]=ESV(2,CINDEX); this_connect[i++]=ESV(2,BINDEX); 
  this_connect[i++]=TSV(0,AINDEX); this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=FSV(2,DINDEX); this_connect[i++]=FSV(2,AINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);


// V1:
  i = 0;
  this_connect[i++]=ESV(0,CINDEX); this_connect[i++]=ESV(0,BINDEX); this_connect[i++]=FSV(3,CINDEX); this_connect[i++]=FSV(3,DINDEX); 
  this_connect[i++]=FSV(0,DINDEX); this_connect[i++]=FSV(0,BINDEX); this_connect[i++]=TSV(0,BINDEX); this_connect[i++]=TSV(0,EINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(0,DINDEX); this_connect[i++]=FSV(0,BINDEX); this_connect[i++]=TSV(0,BINDEX); this_connect[i++]=TSV(0,EINDEX); 
  this_connect[i++]=ESV(4,CINDEX); this_connect[i++]=ESV(4,AINDEX); this_connect[i++]=FSV(1,AINDEX); this_connect[i++]=FSV(1,DINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(1,DINDEX); this_connect[i++]=FSV(1,AINDEX); this_connect[i++]=TSV(0,BINDEX); this_connect[i++]=TSV(0,EINDEX); 
  this_connect[i++]=ESV(1,CINDEX); this_connect[i++]=ESV(1,AINDEX); this_connect[i++]=FSV(3,CINDEX); this_connect[i++]=FSV(3,DINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);


// V2:
  i = 0;
  this_connect[i++]=FSV(3,DINDEX); this_connect[i++]=ESV(1,CINDEX); this_connect[i++]=ESV(1,BINDEX); this_connect[i++]=FSV(3,BINDEX); 
  this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=FSV(1,DINDEX); this_connect[i++]=FSV(1,BINDEX); this_connect[i++]=TSV(0,CINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=FSV(1,DINDEX); this_connect[i++]=FSV(1,BINDEX); this_connect[i++]=TSV(0,CINDEX); 
  this_connect[i++]=FSV(2,DINDEX); this_connect[i++]=ESV(5,CINDEX); this_connect[i++]=ESV(5,AINDEX); this_connect[i++]=FSV(2,CINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=TSV(0,CINDEX); this_connect[i++]=FSV(2,CINDEX); this_connect[i++]=ESV(2,AINDEX); this_connect[i++]=FSV(3,BINDEX); 
  this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=FSV(2,DINDEX); this_connect[i++]=ESV(2,CINDEX); this_connect[i++]=FSV(3,DINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);


// V3:
  i = 0;
  this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=FSV(1,DINDEX); this_connect[i++]=ESV(5,CINDEX); this_connect[i++]=FSV(2,DINDEX); 
  this_connect[i++]=TSV(0,DINDEX); this_connect[i++]=FSV(1,CINDEX); this_connect[i++]=ESV(5,BINDEX); this_connect[i++]=FSV(2,BINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(0,DINDEX); this_connect[i++]=ESV(4,CINDEX); this_connect[i++]=FSV(1,DINDEX); this_connect[i++]=TSV(0,EINDEX); 
  this_connect[i++]=FSV(0,CINDEX); this_connect[i++]=ESV(4,BINDEX); this_connect[i++]=FSV(1,CINDEX); this_connect[i++]=TSV(0,DINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=ESV(3,CINDEX); this_connect[i++]=FSV(0,DINDEX); this_connect[i++]=TSV(0,EINDEX); this_connect[i++]=FSV(2,DINDEX); 
  this_connect[i++]=ESV(3,BINDEX); this_connect[i++]=FSV(0,CINDEX); this_connect[i++]=TSV(0,DINDEX); this_connect[i++]=FSV(2,BINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  interstic_hexes.push_back(this_hex);


// V0:
  i = 0;
  this_connect[i++]=CV(V0INDEX); this_connect[i++]=ESV(0,DINDEX); this_connect[i++]=FSV(3,EINDEX); this_connect[i++]=ESV(2,EINDEX); 
  this_connect[i++]=ESV(3,DINDEX); this_connect[i++]=FSV(0,EINDEX); this_connect[i++]=TSV(0,FINDEX); this_connect[i++]=FSV(2,EINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=ESV(0,DINDEX); this_connect[i++]=ESV(0,AINDEX); this_connect[i++]=FSV(3,AINDEX); this_connect[i++]=FSV(3,EINDEX); 
  this_connect[i++]=FSV(0,EINDEX); this_connect[i++]=FSV(0,AINDEX); this_connect[i++]=TSV(0,AINDEX); this_connect[i++]=TSV(0,FINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(3,EINDEX); this_connect[i++]=FSV(3,AINDEX); this_connect[i++]=ESV(2,BINDEX); this_connect[i++]=ESV(2,EINDEX); 
  this_connect[i++]=TSV(0,FINDEX); this_connect[i++]=TSV(0,AINDEX); this_connect[i++]=FSV(2,AINDEX); this_connect[i++]=FSV(2,EINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=TSV(0,FINDEX); this_connect[i++]=TSV(0,AINDEX); this_connect[i++]=FSV(2,AINDEX); this_connect[i++]=FSV(2,EINDEX); 
  this_connect[i++]=FSV(0,EINDEX); this_connect[i++]=FSV(0,AINDEX); this_connect[i++]=ESV(3,AINDEX); this_connect[i++]=ESV(3,DINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);


// V1:
  i = 0;
  this_connect[i++]=CV(V1INDEX); this_connect[i++]=ESV(1,DINDEX); this_connect[i++]=FSV(3,GINDEX); this_connect[i++]=ESV(0,EINDEX); 
  this_connect[i++]=ESV(4,DINDEX); this_connect[i++]=FSV(1,EINDEX); this_connect[i++]=TSV(0,GINDEX); this_connect[i++]=FSV(0,FINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(3,GINDEX); this_connect[i++]=ESV(1,DINDEX); this_connect[i++]=ESV(1,AINDEX); this_connect[i++]=FSV(3,CINDEX); 
  this_connect[i++]=TSV(0,GINDEX); this_connect[i++]=FSV(1,EINDEX); this_connect[i++]=FSV(1,AINDEX); this_connect[i++]=TSV(0,BINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=TSV(0,GINDEX); this_connect[i++]=FSV(1,EINDEX); this_connect[i++]=FSV(1,AINDEX); this_connect[i++]=TSV(0,BINDEX); 
  this_connect[i++]=FSV(0,FINDEX); this_connect[i++]=ESV(4,DINDEX); this_connect[i++]=ESV(4,AINDEX); this_connect[i++]=FSV(0,BINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=ESV(0,BINDEX); this_connect[i++]=ESV(0,EINDEX); this_connect[i++]=FSV(3,GINDEX); this_connect[i++]=FSV(3,CINDEX); 
  this_connect[i++]=FSV(0,BINDEX); this_connect[i++]=FSV(0,FINDEX); this_connect[i++]=TSV(0,GINDEX); this_connect[i++]=TSV(0,BINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);


// V2:
  i = 0;
  this_connect[i++]=ESV(1,BINDEX); this_connect[i++]=ESV(1,EINDEX); this_connect[i++]=FSV(3,FINDEX); this_connect[i++]=FSV(3,BINDEX); 
  this_connect[i++]=FSV(1,BINDEX); this_connect[i++]=FSV(1,FINDEX); this_connect[i++]=TSV(0,HINDEX); this_connect[i++]=TSV(0,CINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(3,FINDEX); this_connect[i++]=ESV(1,EINDEX); this_connect[i++]=CV(V2INDEX); this_connect[i++]=ESV(2,DINDEX); 
  this_connect[i++]=TSV(0,HINDEX); this_connect[i++]=FSV(1,FINDEX); this_connect[i++]=ESV(5,DINDEX); this_connect[i++]=FSV(2,GINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=TSV(0,HINDEX); this_connect[i++]=FSV(1,FINDEX); this_connect[i++]=ESV(5,DINDEX); this_connect[i++]=FSV(2,GINDEX); 
  this_connect[i++]=TSV(0,CINDEX); this_connect[i++]=FSV(1,BINDEX); this_connect[i++]=ESV(5,AINDEX); this_connect[i++]=FSV(2,CINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(3,BINDEX); this_connect[i++]=FSV(3,FINDEX); this_connect[i++]=ESV(2,DINDEX); this_connect[i++]=ESV(2,AINDEX); 
  this_connect[i++]=TSV(0,CINDEX); this_connect[i++]=TSV(0,HINDEX); this_connect[i++]=FSV(2,GINDEX); this_connect[i++]=FSV(2,CINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);


// V3:
  i = 0;
  this_connect[i++]=FSV(0,CINDEX); this_connect[i++]=ESV(4,BINDEX); this_connect[i++]=FSV(1,EINDEX); this_connect[i++]=TSV(0,DINDEX); 
  this_connect[i++]=FSV(0,GINDEX); this_connect[i++]=ESV(4,EINDEX); this_connect[i++]=FSV(1,GINDEX); this_connect[i++]=TSV(0,IINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=ESV(3,BINDEX); this_connect[i++]=FSV(0,CINDEX); this_connect[i++]=TSV(0,DINDEX); this_connect[i++]=FSV(2,BINDEX); 
  this_connect[i++]=ESV(3,EINDEX); this_connect[i++]=FSV(0,GINDEX); this_connect[i++]=TSV(0,IINDEX); this_connect[i++]=FSV(2,FINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=TSV(0,DINDEX); this_connect[i++]=FSV(1,CINDEX); this_connect[i++]=ESV(5,BINDEX); this_connect[i++]=FSV(2,BINDEX); 
  this_connect[i++]=TSV(0,IINDEX); this_connect[i++]=FSV(1,GINDEX); this_connect[i++]=ESV(5,EINDEX); this_connect[i++]=FSV(2,FINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  i = 0;
  this_connect[i++]=FSV(0,GINDEX); this_connect[i++]=ESV(4,EINDEX); this_connect[i++]=FSV(1,GINDEX); this_connect[i++]=TSV(0,IINDEX); 
  this_connect[i++]=ESV(3,EINDEX); this_connect[i++]=CV(V3INDEX); this_connect[i++]=ESV(5,EINDEX); this_connect[i++]=FSV(2,FINDEX);
  result = gMB->create_element(MBHEX, this_connect, 8, this_hex); RR;
  sphere_hexes.push_back(this_hex);

  return result;
}

MBErrorCode retrieve_subdiv_verts(MBEntityHandle tet, MBEntityHandle this_ent,
                                  const MBEntityHandle *tet_conn,
                                  const int dim, MBEntityHandle *subdiv_verts) 
{
    // get the subdiv verts for this entity
  MBErrorCode result;
  
    // if it's a tet, just put them on the end & return
  if (tet == this_ent) {
    result = gMB->tag_get_data(subdiv_vertices_tag, &this_ent, 1, &subdiv_verts[90]);
    return MB_SUCCESS;
  }
  
    // if it's a sub-entity, need to find index, relative orientation, and offset
    // get connectivity of sub-entity
  std::vector<MBEntityHandle> this_conn;
  result = gMB->get_connectivity(&this_ent, 1, this_conn); RR;
  
    // get relative orientation
  int sense, side_no, offset;
  int success = MBCN::SideNumber(tet_conn, MBTET, &this_conn[0],
                                 this_conn.size(), dim, side_no, sense, offset);
  if (-1 == success) return MB_FAILURE;
  
    // start of this entity's subdiv_verts; edges go first, then preceding sides, then this one;
    // this assumes 6 edges/tet
  MBEntityHandle *subdiv_start = &subdiv_verts[((dim-1)*6 + side_no) * 9];
  
    // get subdiv_verts and put them into proper place
  result = gMB->tag_get_data(subdiv_vertices_tag, &this_ent, 1, subdiv_start);

    // could probably do this more elegantly, but isn't worth it
#define SWITCH(a,b) {MBEntityHandle tmp_handle = a; a = b; b = tmp_handle;}
  switch (dim) {
    case 1:
      if (offset != 0 || sense == -1) {
        SWITCH(subdiv_start[0], subdiv_start[1]);
        SWITCH(subdiv_start[3], subdiv_start[4]);
      }
      break;
    case 2:
        // rotate first
      if (0 != offset) {
        std::rotate(subdiv_start, subdiv_start+offset, subdiv_start+3);
        std::rotate(subdiv_start+4, subdiv_start+4+offset, subdiv_start+7);
      }
        // now flip, if necessary
      if (-1 == sense) {
        SWITCH(subdiv_start[1], subdiv_start[2]);
        SWITCH(subdiv_start[5], subdiv_start[6]);
      }
      break;
    default:
      return MB_FAILURE;
  }
  
    // ok, we're done
  return MB_SUCCESS;
}
