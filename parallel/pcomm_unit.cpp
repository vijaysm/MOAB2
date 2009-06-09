#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "MeshTopoUtil.hpp"
#include "TestUtil.hpp"
#include <algorithm>
#include <vector>
#include <set>

#ifdef USE_MPI
#  include <mpi.h>
#endif

/** Test pack/unpack of vertices */
void test_pack_vertices();
/** Test pack/unpack of elements */
void test_pack_elements();
/** Test pack/unpack of higher-order elements */
void test_pack_higher_order();
/** Test pack/unpack of polygons & polyhedra */
void test_pack_poly();
/** Test pack/unpack of entity sets */
void test_pack_sets_simple();
/** Test pack/unpack of entity sets including implicit packing of set contents */
void test_pack_set_contents();
/** Test pack/unpack of entity sets containing entity sets */
void test_pack_sets_of_sets();
/** Test pack/unpack of set parent/child relations */
void test_pack_set_parent_child();
/** Test pack/unpack tag values*/
void test_pack_tag_data_sparse();
void test_pack_tag_data_dense();
void test_pack_tag_data_default_value();
/** Test pack/unpack tag values*/
void test_pack_bit_tag_data();
/** Test pack/unpack of variable length tag values*/
void test_pack_variable_length_tag();
/** Test pack/unpack tag values*/
void test_pack_tag_handle_data();
/** Test pack/unpack of shared entities*/
void test_pack_shared_entities();
/** Test filter_pstatus function*/
void test_filter_pstatus();

int main( int argc, char* argv[] )
{
#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif

  int num_err = 0;
  num_err += RUN_TEST( test_pack_vertices );
  num_err += RUN_TEST( test_pack_elements );
  num_err += RUN_TEST( test_pack_higher_order );
  num_err += RUN_TEST( test_pack_poly );
  num_err += RUN_TEST( test_pack_sets_simple );
  num_err += RUN_TEST( test_pack_set_contents );
  num_err += RUN_TEST( test_pack_sets_of_sets );
  num_err += RUN_TEST( test_pack_set_parent_child );
  num_err += RUN_TEST( test_pack_tag_data_sparse );
  num_err += RUN_TEST( test_pack_tag_data_dense );
  num_err += RUN_TEST( test_pack_tag_data_default_value );
  //num_err += RUN_TEST( test_pack_bit_tag_data );
  num_err += RUN_TEST( test_pack_variable_length_tag );
  num_err += RUN_TEST( test_pack_tag_handle_data );
  num_err += RUN_TEST( test_pack_shared_entities );
  num_err += RUN_TEST( test_filter_pstatus );
  
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return num_err;
}

/* Utility method: pack mesh data, clear moab instance, and unpack */
void pack_unpack_noremoteh( MBCore& moab, MBRange& entities )
{
  MBErrorCode rval;
  if (entities.empty()) {
    rval = moab.get_entities_by_handle( 0, entities );
    CHECK_ERR(rval);
  }
  
  MBParallelComm *pcomm = new MBParallelComm( &moab );
  int size = 0;
  std::vector<unsigned char> buff;
  std::vector<int> addl_procs;

    // get the necessary vertices too
  MBRange tmp_range = entities.subset_by_type(MBENTITYSET);
  entities = entities.subtract(tmp_range);
  rval = moab.get_adjacencies(entities, 0, false, entities, MBInterface::UNION);
  CHECK_ERR(rval);
  entities.merge(tmp_range);
  
  rval = pcomm->pack_buffer( entities, false, true, false, 
                             -1, buff, size);
  CHECK_ERR(rval);
  
  delete pcomm;
  moab.~MBCore();

  new (&moab) MBCore();
  pcomm = new MBParallelComm( &moab);
  
  entities.clear();
  std::vector<std::vector<MBEntityHandle> > L1h;
  std::vector<std::vector<int> > L1p;
  std::vector<MBEntityHandle> L2hloc, L2hrem;
  std::vector<unsigned int> L2p;
  rval = pcomm->unpack_buffer( &buff[0], false, -1, -1, L1h, L1p, L2hloc, 
                               L2hrem, L2p, entities);
  CHECK_ERR(rval);

  delete pcomm;
}
void pack_unpack_noremoteh( MBCore& moab )
{
  MBRange empty;
  pack_unpack_noremoteh( moab, empty );
}

/* Utility method -- check expected sizes */
void check_sizes( MBInterface& moab,
                  int num_vtx,
                  int num_edge,
                  int num_tri,
                  int num_quad,
                  int num_polygon,
                  int num_tet,
                  int num_pyr,
                  int num_wedge,
                  int num_knife,
                  int num_hex,
                  int num_polyhedron )
{
  int count;
  MBErrorCode rval;
  
  rval = moab.get_number_entities_by_type( 0, MBVERTEX, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_vtx, count );

  rval = moab.get_number_entities_by_type( 0, MBEDGE, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_edge, count );

  rval = moab.get_number_entities_by_type( 0, MBTRI, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tri, count );

  rval = moab.get_number_entities_by_type( 0, MBQUAD, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_quad, count );

  rval = moab.get_number_entities_by_type( 0, MBPOLYGON, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_polygon, count );

  rval = moab.get_number_entities_by_type( 0, MBTET, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tet, count );

  rval = moab.get_number_entities_by_type( 0, MBPYRAMID, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_pyr, count );

  rval = moab.get_number_entities_by_type( 0, MBPRISM, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_wedge, count );

  rval = moab.get_number_entities_by_type( 0, MBKNIFE, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_knife, count );

  rval = moab.get_number_entities_by_type( 0, MBHEX, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_hex, count );

  rval = moab.get_number_entities_by_type( 0, MBPOLYHEDRON, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_polyhedron, count );
}

/* Create a simple mesh for use in tests */
void create_simple_grid( MBInterface& moab, unsigned x, unsigned y, unsigned z );
void create_simple_grid( MBInterface& moab, unsigned xyz = 3 )
{
  create_simple_grid( moab, xyz, xyz, xyz );
}
void create_simple_grid( MBInterface& moab, unsigned x, unsigned y, unsigned z )
{
  MBErrorCode rval;
  MBEntityHandle *verts = new MBEntityHandle[x*y*z];
  for (unsigned k = 0; k < z; ++k)
    for (unsigned j = 0; j < y; ++j)
      for (unsigned i = 0; i < x; ++i) {
        const double coords[3] = { i, j, k };
        rval = moab.create_vertex( coords, verts[x*y*k + x*j + i] );
        CHECK_ERR(rval);
      }
  
  MBEntityHandle *elems = new MBEntityHandle[(x-1)*(y-1)*(z-1)];
  for (unsigned k = 0; k < (z-1); ++k)
    for (unsigned j = 0; j < (y-1); ++j)
      for (unsigned i = 0; i < (x-1); ++i) {
        const size_t idx = (size_t)i + (size_t)j*x + (size_t)k*x*y; 
        const MBEntityHandle conn[8] = { verts[idx        ],
                                         verts[idx      +1],
                                         verts[idx    +x+1],
                                         verts[idx    +x  ],
                                         verts[idx+x*y    ],
                                         verts[idx+x*y  +1],
                                         verts[idx+x*y+x+1],
                                         verts[idx+x*y+x  ] };
        rval = moab.create_element( MBHEX, conn, 8, elems[(x-1)*(y-1)*k + (x-1)*j + i] );
        CHECK_ERR(rval);
      }
  delete [] verts;
  delete [] elems;
}

MBErrorCode set_owners(unsigned char pstatus,
                       MBParallelComm *pc0, MBEntityHandle ent0, 
                       MBParallelComm *pc1, MBEntityHandle ent1, 
                       MBParallelComm *pc2 = NULL, MBEntityHandle ent2 = 0, 
                       MBParallelComm *pc3 = NULL, MBEntityHandle ent3 = 0)
{
  int owners[MAX_SHARING_PROCS];
  MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
  
  MBErrorCode result = MB_SUCCESS;
  
  std::fill(owners, owners+MAX_SHARING_PROCS, -1);
  std::fill(tmp_handles, tmp_handles+MAX_SHARING_PROCS, 0);
  owners[0] = pc0->proc_config().proc_rank(); owners[1] = pc1->proc_config().proc_rank();
  tmp_handles[0] = ent0; tmp_handles[1] = ent1;
  int np = 2;
  if (pc2) {
    owners[2] = pc2->proc_config().proc_rank();
    tmp_handles[2] = ent2;
    np++;
  }
  if (pc3) {
    owners[3] = pc3->proc_config().proc_rank();
    tmp_handles[3] = ent3;
    np++;
  }
  if (np > 2) {
    result = pc0->get_moab()->tag_set_data(pc0->sharedps_tag(), &ent0, 1, owners);
    result = pc0->get_moab()->tag_set_data(pc0->sharedhs_tag(), &ent0, 1, tmp_handles);
    result = pc0->get_moab()->tag_set_data(pc0->pstatus_tag(), &ent0, 1, &pstatus);
      // 2nd and greater don't own them
    pstatus |= PSTATUS_NOT_OWNED;
    result = pc1->get_moab()->tag_set_data(pc1->sharedps_tag(), &ent1, 1, owners);
    result = pc1->get_moab()->tag_set_data(pc1->sharedhs_tag(), &ent1, 1, tmp_handles);
    result = pc1->get_moab()->tag_set_data(pc1->pstatus_tag(), &ent1, 1, &pstatus);
    if (pc2) {
      result = pc2->get_moab()->tag_set_data(pc2->sharedps_tag(), &ent2, 1, owners);
      result = pc2->get_moab()->tag_set_data(pc2->sharedhs_tag(), &ent2, 1, tmp_handles);
      result = pc2->get_moab()->tag_set_data(pc2->pstatus_tag(), &ent2, 1, &pstatus);
    }
    if (pc3) {
      result = pc3->get_moab()->tag_set_data(pc3->sharedps_tag(), &ent3, 1, owners);
      result = pc3->get_moab()->tag_set_data(pc3->sharedhs_tag(), &ent3, 1, tmp_handles);
      result = pc3->get_moab()->tag_set_data(pc3->pstatus_tag(), &ent3, 1, &pstatus);
    }
  }
  else {
    result = pc0->get_moab()->tag_set_data(pc0->sharedp_tag(), &ent0, 1, &owners[1]);
    result = pc0->get_moab()->tag_set_data(pc0->sharedh_tag(), &ent0, 1, &tmp_handles[1]);
    result = pc0->get_moab()->tag_set_data(pc0->pstatus_tag(), &ent0, 1, &pstatus);
      // 2nd and greater don't own them
    pstatus |= PSTATUS_NOT_OWNED;
    result = pc1->get_moab()->tag_set_data(pc1->sharedp_tag(), &ent1, 1, owners);
    result = pc1->get_moab()->tag_set_data(pc1->sharedh_tag(), &ent1, 1, tmp_handles);
    result = pc1->get_moab()->tag_set_data(pc1->pstatus_tag(), &ent1, 1, &pstatus);
  }
  
  return result;
}

MBErrorCode create_patch(MBInterface *moab, MBRange &verts, MBRange &quads,
                         unsigned int n, double *xyz, int *gids) 
{
    // create vertices/quads in square array
  MBErrorCode result = moab->create_vertices(xyz, n*n, verts);
  if (MB_SUCCESS != result) return result;
  std::vector<MBEntityHandle> connect;
  for (unsigned int j = 0; j < n-1; j++) {
    for (unsigned int i = 0; i < n-1; i++) {
      connect.push_back(verts[n*j+i]);
      connect.push_back(verts[n*j+i+1]);
      connect.push_back(verts[n*(j+1)+i+1]);
      connect.push_back(verts[n*(j+1)+i]);
    }
  }
  
  unsigned int nquads = (n-1)*(n-1);
  for (unsigned int i = 0; i < nquads; i++) {
    MBEntityHandle dum_quad;
    result = moab->create_element(MBQUAD, &connect[4*i], 4, dum_quad);
    if (MB_SUCCESS != result) return result;
    quads.insert(dum_quad);
  }
  
    // global ids
  MBTag gid_tag;
  int dum_default = -1;
  result = moab->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                           MB_TYPE_INTEGER, gid_tag, &dum_default, true);
  if (MB_SUCCESS != result) return result;
  result = moab->tag_set_data(gid_tag, verts, gids);
  if (MB_SUCCESS != result) return result;

  return result;
}

MBErrorCode set_owners(MBParallelComm **pc, MBRange *verts, MBRange *quads,
                       const bool edges_too) 
{
    // P0-P2
  unsigned char pstat = PSTATUS_SHARED | PSTATUS_INTERFACE;
  MBErrorCode rval = set_owners(pstat, pc[0], verts[0][6], pc[2], verts[2][0]);
  rval = set_owners(pstat, pc[0], verts[0][7], pc[2], verts[2][1]); CHECK_ERR(rval);

    // P1-P2
  rval = set_owners(pstat, pc[1], verts[1][7], pc[2], verts[2][5]); CHECK_ERR(rval);
  rval = set_owners(pstat, pc[1], verts[1][8], pc[2], verts[2][8]); CHECK_ERR(rval);
    // P0-P1-P3
  pstat |= PSTATUS_MULTISHARED;
  rval = set_owners(pstat, pc[0], verts[0][2], pc[1], verts[1][0], pc[3], verts[3][0]); CHECK_ERR(rval);
  rval = set_owners(pstat, pc[0], verts[0][5], pc[1], verts[1][3], pc[3], verts[3][3]); CHECK_ERR(rval);
    // P0-P1-P2-P3
  rval = set_owners(pstat, pc[0], verts[0][8], pc[1], verts[1][6], 
                    pc[2], verts[2][2], pc[3], verts[3][6]); CHECK_ERR(rval);

  if (edges_too) {
    MeshTopoUtil *mtu[4];
    MBInterface *mb[4];
    MBRange dum_range;
      // create mtu's and explicit edges
    for (unsigned int i = 0; i < 4; i++) {
      mb[i] = pc[i]->get_moab();
      assert(mb[i]);
      mtu[i] = new MeshTopoUtil(mb[i]);
      rval = mb[i]->get_adjacencies(quads[i], 1, true, dum_range, 
                                    MBInterface::UNION);
      CHECK_ERR(rval);
      dum_range.clear();
    }
    
    MBEntityHandle edge1, edge2, edge3;
    pstat = PSTATUS_SHARED | PSTATUS_INTERFACE;
      // P0-P2
    edge1 = mtu[0]->common_entity(verts[0][6], verts[0][7], 1);
    edge2 = mtu[2]->common_entity(verts[2][0], verts[2][1], 1);
    assert(edge1 && edge2);
    rval = set_owners(pstat, pc[0], edge1, pc[2], edge2); CHECK_ERR(rval);
    edge1 = mtu[0]->common_entity(verts[0][7], verts[0][8], 1);
    edge2 = mtu[2]->common_entity(verts[2][1], verts[2][2], 1);
    assert(edge1 && edge2);
    rval = set_owners(pstat, pc[0], edge1, pc[2], edge2); CHECK_ERR(rval);
      // P1-P2
    edge1 = mtu[1]->common_entity(verts[1][6], verts[1][7], 1);
    edge2 = mtu[2]->common_entity(verts[2][2], verts[2][5], 1);
    assert(edge1 && edge2);
    rval = set_owners(pstat, pc[1], edge1, pc[2], edge2); CHECK_ERR(rval);
    edge1 = mtu[1]->common_entity(verts[1][7], verts[1][8], 1);
    edge2 = mtu[2]->common_entity(verts[2][5], verts[2][8], 1);
    assert(edge1 && edge2);
    rval = set_owners(pstat, pc[1], edge1, pc[2], edge2); CHECK_ERR(rval);
      // P0-P1-P3
    pstat |= PSTATUS_MULTISHARED;
    edge1 = mtu[0]->common_entity(verts[0][2], verts[0][5], 1);
    edge2 = mtu[1]->common_entity(verts[1][0], verts[1][3], 1);
    edge3 = mtu[3]->common_entity(verts[3][0], verts[3][3], 1);
    assert(edge1 && edge2 && edge3);
    rval = set_owners(pstat, pc[0], edge1, pc[1], edge2, pc[3], edge3); CHECK_ERR(rval);
    edge1 = mtu[0]->common_entity(verts[0][5], verts[0][8], 1);
    edge2 = mtu[1]->common_entity(verts[1][3], verts[1][6], 1);
    edge3 = mtu[3]->common_entity(verts[3][3], verts[3][6], 1);
    assert(edge1 && edge2 && edge3);
    rval = set_owners(pstat, pc[0], edge1, pc[1], edge2, pc[3], edge3); CHECK_ERR(rval);

    for (unsigned int i = 0; i < 4; i++)
      delete mtu[i];
  }
    
  return MB_SUCCESS;
}

MBErrorCode create_shared_grid(MBParallelComm **pc, MBRange *verts, MBRange *quads) 
{
//          
//        P2______
//         /__/__/ /|P1
//    y:  /__/__/ /||
//     1  _____   ||/     _____  1   
//       |  |  |  |/|-1  |  |  | 
//    .5 |__|__|  ||/    |__|__| .5
//       |P0|  |  |/-.5  |P3|  | 
//     0 |__|__| z:0     |__|__| 0
//    x:-1 -.5  0       0  .5  1  
//
// nodes:   P2
//            18 16 14      P1
//          17 15 13      14               
//         6  7  8     13 12            
//                   8 11 10               
//       6  7  8     5  9      8  23 24 P3
//       3  4  5     2         5  21 22 
//   P0  0  1  2               2  19 20 

  int gids[] = {
      0, 1, 2, 3, 4, 5, 6, 7, 8, // P0
      2, 9, 10, 5, 11, 12, 8, 13, 14, // P1
      6, 7, 8, 17, 15, 13, 18, 16, 14, // P2
      2, 19, 20, 5, 21, 22, 8, 23, 24 // P3
  };
  double xyz[] =  {
      -1.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.5, 0.0,
      -0.5, 0.5, 0.0, 0.0, 0.5, 0.0, -1.0, 1.0, 0.0, -0.5, 1.0, 0.0,
      0.0, 1.0, 0.0, // n0-8
      0.0, 0.0, -0.5, 0.0, 0.0, -1.0, 0.0, 0.5, -0.5, 
      0.0, 0.5, -1.0, 0.0, 1.0, -0.5, 0.0, 1.0, -1.0, // n9-14
      -0.5, 1.0, -0.5, -0.5, 1.0, -1.0, -1.0, 1.0, -0.5, -1.0, 1.0, -1.0, // n15-18
      0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5, 0.0, 1.0, 0.5, 0.0, 
      0.5, 1.0, 0.0, 1.0, 1.0, 0.0, // n19-24
  };
  double xyztmp[27];

    // create the test mesh above
  for (unsigned int i = 0; i < 4; i++) {
    for (unsigned int j = 0; j < 9; j++) {
      xyztmp[3*j] = xyz[3*gids[9*i+j]];
      xyztmp[3*j+1] = xyz[3*gids[9*i+j]+1];
      xyztmp[3*j+2] = xyz[3*gids[9*i+j]+2];
    }
    
    create_patch(pc[i]->get_moab(), verts[i], quads[i], 3, xyztmp, &gids[9*i]);
  }

  MBErrorCode rval = set_owners(pc, verts, quads, false);
  CHECK_ERR(rval);
  return rval;
}

void test_pack_vertices()
{
  MBCore moab;
  MBErrorCode rval;
  MBRange verts;
  
  const size_t num_verts = 4;
  const double coords[3*num_verts] = { -0.5, -1./3, 0.0,
                                        0.5, -1./3, 0.0,
                                        0.0,  2./3, 0.0,
                                        0.0,  0.0,  0.745356 };
  
  rval = moab.create_vertices( coords, num_verts, verts );
  CHECK_ERR(rval);
  
  pack_unpack_noremoteh( moab, verts );
  CHECK_EQUAL( num_verts, verts.size() );
  
  double coords2[3*num_verts];
  rval = moab.get_coords( verts, coords2 );
  CHECK_ERR(rval);
  
  std::vector<bool> seen( num_verts, false );
  for (unsigned i = 0; i < num_verts; ++i) {
    unsigned j;
    for (j = 0; j < num_verts; ++j) 
      if (coords[3*j  ] == coords2[3*i  ] &&
          coords[3*j+1] == coords2[3*i+1] &&
          coords[3*j+2] == coords2[3*i+2])
        break;
    CHECK( j < num_verts );
    CHECK( !seen[j] );
    seen[j] = true;
  }
}

void test_pack_elements()
{
  MBCore moab;
  MBErrorCode rval;
  MBRange elems;
  
    // define some vertices
  const size_t num_verts = 12;
  MBEntityHandle verts[num_verts];
  const double hex_corners[3*num_verts] = { -1, -1, -1,
                                             1, -1, -1,
                                             1,  1, -1,
                                            -1,  1, -1,
                                            -1, -1,  0,
                                             1, -1,  0,
                                             1,  1,  0,
                                            -1,  1,  0,
                                            -1, -1,  1,
                                             1, -1,  1,
                                             1,  1,  1,
                                            -1,  1,  1 };
  for (size_t i = 0; i < num_verts; ++i) {
    rval = moab.create_vertex( hex_corners + 3*i, verts[i] );
    CHECK_ERR(rval);
  }
  
    // define two adjacent hexes
  const size_t num_hex = 2;
  MBEntityHandle hexes[num_hex];
  rval = moab.create_element( MBHEX, verts, 8, hexes[0] );
  CHECK_ERR(rval);
  elems.insert( hexes[0] );
  rval = moab.create_element( MBHEX, verts+4, 8, hexes[1] );
  CHECK_ERR(rval);
  elems.insert( hexes[1] );
  
    // define a single quad on the adjacent sides of the hexes
  const size_t num_quad = 1;
  MBEntityHandle quad;
  rval = moab.create_element( MBQUAD, verts+4, 4, quad );
  CHECK_ERR(rval);
  elems.insert( quad );
  
    // define a decomposition of the first hex into 5 tets
  const size_t num_tet = 5;
  MBEntityHandle tets[num_tet];
  MBEntityHandle tet_conn[num_tet][4] = 
                               { { verts[0], verts[1], verts[3], verts[4] },
                                 { verts[1], verts[2], verts[3], verts[6] },
                                 { verts[1], verts[3], verts[4], verts[6] },
                                 { verts[4], verts[6], verts[5], verts[1] },
                                 { verts[7], verts[6], verts[4], verts[3] } };
  rval = moab.create_element( MBTET, tet_conn[0], 4, tets[0] );
  CHECK_ERR(rval);
  elems.insert( tets[0] );
  rval = moab.create_element( MBTET, tet_conn[1], 4, tets[1] );
  CHECK_ERR(rval);
  elems.insert( tets[1] );
  rval = moab.create_element( MBTET, tet_conn[2], 4, tets[2] );
  CHECK_ERR(rval);
  elems.insert( tets[2] );
  rval = moab.create_element( MBTET, tet_conn[3], 4, tets[3] );
  CHECK_ERR(rval);
  elems.insert( tets[3] );
  rval = moab.create_element( MBTET, tet_conn[4], 4, tets[4] );
  CHECK_ERR(rval);
  elems.insert( tets[4] );
  
    // define the 4 shared faces of the above tets as tris
    // (the faces of the 3rd tet)
  const size_t num_tri = 4;
  MBEntityHandle tris[num_tri];
  MBEntityHandle tri_conn[num_tri][3] = { { verts[3], verts[1], verts[4] },
                                          { verts[1], verts[6], verts[4] },
                                          { verts[1], verts[3], verts[6] },
                                          { verts[3], verts[4], verts[6] } };
  rval = moab.create_element( MBTRI, tri_conn[0], 3, tris[0] );
  CHECK_ERR(rval);
  elems.insert( tris[0] );
  rval = moab.create_element( MBTRI, tri_conn[1], 3, tris[1] );
  CHECK_ERR(rval);
  elems.insert( tris[1] );
  rval = moab.create_element( MBTRI, tri_conn[2], 3, tris[2] );
  CHECK_ERR(rval);
  elems.insert( tris[2] );
  rval = moab.create_element( MBTRI, tri_conn[3], 3, tris[3] );
  CHECK_ERR(rval);
  elems.insert( tris[3] );
  
    // define a decomposition of the second hex into two wedges
  const size_t num_wedge = 2;
  MBEntityHandle wedges[num_wedge];
  MBEntityHandle wedge_conn[num_wedge][6] = {
    { verts[4], verts[5], verts[7], verts[8], verts[9], verts[11] },
    { verts[5], verts[6], verts[7], verts[9], verts[10],verts[11] } };
  rval = moab.create_element( MBPRISM, wedge_conn[0], 6, wedges[0] );
  CHECK_ERR(rval);
  elems.insert( wedges[0] );
  rval = moab.create_element( MBPRISM, wedge_conn[1], 6, wedges[1] );
  CHECK_ERR(rval);
  elems.insert( wedges[1] );
  
    // define a pyramid
  MBEntityHandle pyr;
  rval = moab.create_element( MBPYRAMID, verts, 5, pyr );
  CHECK_ERR(rval);
  elems.insert( pyr );
  
    // pack and unpack mesh
  pack_unpack_noremoteh( moab, elems );
  
    // check_counts
  check_sizes( moab, num_verts, 0, num_tri, num_quad, 0, 
                     num_tet, 1, num_wedge, 0, num_hex, 0 );

    // get connectivity for two hexes and a quad
  MBRange range;
  const MBEntityHandle *conn1, *conn2, *conn3;
  int len1, len2, len3;
  rval = moab.get_entities_by_type( 0, MBHEX, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_hex, range.size() );
  rval = moab.get_connectivity( range.front(), conn1, len1, true );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, len1 );
  rval = moab.get_connectivity( range.back(),  conn2, len2, true );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, len2 );
  range.clear();
  rval = moab.get_entities_by_type( 0, MBQUAD, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_quad, range.size() );
  rval = moab.get_connectivity( range.front(), conn3, len3, true );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len3 );
  
    // Check if hexes are reversed
  if (conn1[0] == conn2[4]) {
    std::swap( conn1, conn2 );
  }
  
    // Check consistant connectivity between hexes
  CHECK_EQUAL( conn1[4], conn2[0] );
  CHECK_EQUAL( conn1[5], conn2[1] );
  CHECK_EQUAL( conn1[6], conn2[2] );
  CHECK_EQUAL( conn1[7], conn2[3] );
    // Check connectivity of quad on shared face
  CHECK_EQUAL( conn1[4], conn3[0] );
  CHECK_EQUAL( conn1[5], conn3[1] );
  CHECK_EQUAL( conn1[6], conn3[2] );
  CHECK_EQUAL( conn1[7], conn3[3] );
    // Check coordinates
  const MBEntityHandle combined[12] = { conn1[0], conn1[1], conn1[2], conn1[3],
                                        conn3[0], conn3[1], conn3[2], conn3[3],
                                        conn2[4], conn2[5], conn2[6], conn2[7] };
  double coords[36];
  rval = moab.get_coords( combined, 12, coords );
  CHECK_ERR(rval);
  for (int i = 0; i < 36; ++i) {
    CHECK_REAL_EQUAL( hex_corners[i], coords[i], 1e-12 );
  }
}

void test_pack_higher_order()
{
  MBCore moab;
  MBErrorCode rval;
  MBRange elems;
  
    // define coordinates for a pyramid decomposed
    // into two 10-node tets
  const size_t num_vert = 14;
  const double coords[3*num_vert] = { -1, -1, 0,
                                       0, -1, 0,
                                       1, -1, 0,
                                       1,  0, 0,
                                       1,  1, 0,
                                       0,  1, 0,
                                      -1,  1, 0,
                                      -1,  0, 0,
                                       0,  0, 0, 
                                       0,  0, 1,
                                     -.5, -.5, .5,
                                      .5, -.5, .5,
                                      .5,  .5, .5,
                                     -.5,  .5, .5 };
  MBEntityHandle verts[num_vert];
  for (size_t i = 0; i < num_vert; ++i) {
    rval = moab.create_vertex( coords + 3*i, verts[i] );
    CHECK_ERR(rval);
  }
  
    // define two tets
  const size_t num_tet = 2;
  MBEntityHandle tet_conn[2][10] = {
   { verts[ 0], verts[ 4], verts[ 9], verts[2],
     verts[ 8], verts[12], verts[10],
     verts[ 1], verts[ 3], verts[11] }, 
   { verts[ 0], verts[ 9], verts[ 4], verts[6],
     verts[10], verts[12], verts[ 8],
     verts[ 7], verts[13], verts[ 5] } };
     
  MBEntityHandle tets[num_tet];
  rval = moab.create_element( MBTET, tet_conn[0], 10, tets[0] );
  CHECK_ERR(rval);
  elems.insert( tets[0] );
  rval = moab.create_element( MBTET, tet_conn[1], 10, tets[1] );
  CHECK_ERR(rval);
  elems.insert( tets[1] );
   
    // define interior tri face
  const size_t num_tri = 1;
  MBEntityHandle tri_conn[6] = 
    { verts[0], verts[4], verts[9],
      verts[8], verts[12],verts[10] };
  MBEntityHandle tri;
  rval = moab.create_element( MBTRI, tri_conn, 6, tri );
  CHECK_ERR(rval);
  elems.insert( tri );
  
    // pack and unpack mesh
  pack_unpack_noremoteh( moab, elems );
  
    // check_counts
  check_sizes( moab, num_vert, 0, num_tri, 0, 0, num_tet, 0, 0, 0, 0, 0 );

    // get connectivity for two tets and a tri
  MBRange range;
  const MBEntityHandle *conn1, *conn2, *conn3;
  int len1, len2, len3;
  rval = moab.get_entities_by_type( 0, MBTET, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tet, range.size() );
  rval = moab.get_connectivity( range.front(), conn1, len1, false );
  CHECK_ERR(rval);
  CHECK_EQUAL( 10, len1 );
  rval = moab.get_connectivity( range.back(),  conn2, len2, false );
  CHECK_ERR(rval);
  CHECK_EQUAL( 10, len2 );
  range.clear();
  rval = moab.get_entities_by_type( 0, MBTRI, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tri, range.size() );
  rval = moab.get_connectivity( range.front(), conn3, len3, false );
  CHECK_ERR(rval);
  CHECK_EQUAL( 6, len3 );
  
    // The first face of one of the tets is in the 
    // same order as the tri.
  if (conn3[1] != conn1[1])
    std::swap( conn1, conn2 );
  
    // check consistant connectivity for face shared by tets
  CHECK_EQUAL( conn1[0], conn2[0] );
  CHECK_EQUAL( conn1[1], conn2[2] );
  CHECK_EQUAL( conn1[2], conn2[1] );
  CHECK_EQUAL( conn1[4], conn2[6] );
  CHECK_EQUAL( conn1[5], conn2[5] );
  CHECK_EQUAL( conn1[6], conn2[4] );
    // check connsistant connectivity between tet face and tri
  CHECK_EQUAL( conn1[0], conn3[0] );
  CHECK_EQUAL( conn1[1], conn3[1] );
  CHECK_EQUAL( conn1[2], conn3[2] );
  CHECK_EQUAL( conn1[4], conn3[3] );
  CHECK_EQUAL( conn1[5], conn3[4] );
  CHECK_EQUAL( conn1[6], conn3[5] );
  
    // order vertex handles corresponding to original coordinate list
  const MBEntityHandle combined[num_vert] = {
    conn1[0], 
    conn1[7],
    conn1[3],
    conn1[8],
    conn1[1],
    conn2[9],
    conn2[3],
    conn2[7],
    conn1[4],
    conn1[2],
    conn1[6],
    conn1[9],
    conn1[5],
    conn2[8] };
  double coords2[3*num_vert];
  rval = moab.get_coords( combined, num_vert, coords2 );
  CHECK_ERR(rval);
  
    // check vertex coordinates
  for (int i = 0; i < 36; ++i) {
    CHECK_REAL_EQUAL( coords[i], coords2[i], 1e-12 );
  }
}


void test_pack_poly()
{
  MBCore moab;
  MBErrorCode rval;
  MBRange elems;
  
    // define a pyramid w/ a octagonal base
  const double a = 0.5;
  const double b = 1+sqrt(2.0)/2.0;
  const size_t num_vert = 9;
  const double coords[3*num_vert] = {
    -a, -b, 0,
     a, -b, 0,
     b, -a, 0,
     b,  a, 0,
     a,  b, 0,
    -a,  b, 0,
    -b,  a, 0,
    -b, -a, 0,
     0,  0, 1 };
  MBEntityHandle verts[num_vert];
  for (size_t i = 0; i < num_vert; ++i) {
    rval = moab.create_vertex( coords + 3*i, verts[i] );
    CHECK_ERR(rval);
  }
  
    // define octagonal base
  const size_t num_polygon = 1;
  MBEntityHandle octagon;
  rval = moab.create_element( MBPOLYGON, verts, 8, octagon );
  CHECK_ERR(rval);
  
    // define triangular sides
  const size_t num_tri = num_vert-1;
  MBEntityHandle tri[num_tri];
  for (size_t i = 0; i < num_tri; ++i) {
    const MBEntityHandle conn[3] = { verts[i], verts[(i+1)%num_tri], verts[num_tri] };
    rval = moab.create_element( MBTRI, conn, 3, tri[i] );
    CHECK_ERR(rval);
  }
  
    // define the octagon-based pyramid
  const size_t num_polyhedron = 1;
  MBEntityHandle polyhedron;
  MBEntityHandle all_faces[num_vert];
  all_faces[0] = octagon;
  std::copy( tri, tri+num_tri, all_faces+1 );
  rval = moab.create_element( MBPOLYHEDRON, all_faces, num_vert, polyhedron );
  CHECK_ERR(rval);
  
    // pack and unpack the mesh
  elems.clear();
  elems.insert( polyhedron );
  elems.insert( octagon );
  std::copy( tri, tri+num_tri, mb_range_inserter(elems) );
  pack_unpack_noremoteh( moab, elems );
  
    // check counts
  check_sizes( moab, num_vert, 0, num_tri, 0, num_polygon, 
                     0, 0, 0, 0, 0, num_polyhedron );
  
    // get entities
  MBRange range;
  rval = moab.get_entities_by_type( 0, MBPOLYHEDRON, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_polyhedron, range.size() );
  polyhedron = range.front();
  
  range.clear();
  rval = moab.get_entities_by_type( 0, MBPOLYGON, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_polygon, range.size() );
  octagon = range.front();
  
  range.clear();
  rval = moab.get_entities_by_type( 0, MBTRI, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tri, range.size() );
  
    // check coords of octagon vertices
  const MBEntityHandle* oct_conn;
  int eight = 0;
  rval = moab.get_connectivity( octagon, oct_conn, eight );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, eight );
  double oct_coords[3*8];
  rval = moab.get_coords( oct_conn, 8, oct_coords );
  CHECK_ERR(rval);
  for (int i = 0; i < 3*8; ++i) {
    CHECK_REAL_EQUAL( coords[i], oct_coords[i], 1e-12 );
  }
  
    // check faces of polyhedron
  std::vector<MBEntityHandle> volconn;
  rval = moab.get_connectivity( &polyhedron, 1, volconn );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tri + 1, volconn.size() );
  CHECK_EQUAL( volconn[0], octagon );
  for (MBRange::iterator i = range.begin(); i != range.end(); ++i) {
    CHECK( std::find(volconn.begin(), volconn.end(), *i) != volconn.end() );
  }
}

void test_pack_sets_simple()
{
  MBCore moab;
  MBErrorCode rval;
  create_simple_grid( moab, 3 );  

    // delete any existing sets
  MBRange sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
  
    // get all entities
  MBRange entities;
  rval = moab.get_entities_by_handle( 0, entities );
  CHECK_ERR(rval);
    // expect 8 elements and 27 vertices
  CHECK_EQUAL( (MBEntityHandle)35, entities.size() );
  CHECK_EQUAL( 27u, entities.num_of_type( MBVERTEX ) );
  CHECK_EQUAL( 8u, entities.num_of_type( MBHEX ) );
  
  
    // create five sets:
    // 1) one with all the elements and vertices,
    // 2) one with half of the elements, 
    // 3) one with the other half of the elements, 
    // 4) one with a single vertex, 
    // 5) an empty set, 
  MBEntityHandle all_set, half1_set, half2_set, vertex_set, empty_set;
  const unsigned int all_opt = MESHSET_SET | MESHSET_TRACK_OWNER,
                     half1_opt = MESHSET_SET,
                     half2_opt = MESHSET_SET,
                     vertex_opt = MESHSET_ORDERED,
                     empty_opt = MESHSET_ORDERED | MESHSET_TRACK_OWNER;
  rval = moab.create_meshset( all_opt, all_set );
  CHECK_ERR(rval);
  entities.insert( all_set );
  rval = moab.create_meshset( half1_opt, half1_set );
  CHECK_ERR(rval);
  entities.insert( half1_set );
  rval = moab.create_meshset( half2_opt, half2_set );
  CHECK_ERR(rval);
  entities.insert( half2_set );
  rval = moab.create_meshset( vertex_opt, vertex_set );
  CHECK_ERR(rval);
  entities.insert( vertex_set );
  rval = moab.create_meshset( empty_opt, empty_set );
  CHECK_ERR(rval);
  entities.insert( empty_set );

  MBRange elems, verts, half;
  rval = moab.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  rval = moab.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  
  rval = moab.add_entities( all_set, verts );
  CHECK_ERR(rval);
  rval = moab.add_entities( all_set, elems );
  CHECK_ERR(rval);
  half.merge( elems.begin(), elems.begin() += elems.size() / 2 );
  rval = moab.add_entities( half1_set, half );
  CHECK_ERR(rval);
  half.clear();
  half.merge( elems.begin() += elems.size() / 2, elems.end() );
  rval = moab.add_entities( half2_set, half );
  CHECK_ERR(rval);
  MBEntityHandle vert = verts.front();
  rval = moab.add_entities( vertex_set, &vert, 1 );
  CHECK_ERR(rval);
  
    // do pack and unpack
  pack_unpack_noremoteh( moab, entities );
  
    // get entities by type
  verts.clear();
  rval = moab.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  elems.clear();
  rval = moab.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  sets.clear();
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  
  CHECK_EQUAL( (MBEntityHandle)27, verts.size() );
  CHECK_EQUAL( (MBEntityHandle)8, elems.size() );
  CHECK_EQUAL( (MBEntityHandle)5, sets.size() );
  
    // guess which is which
  empty_set = all_set = vertex_set = half1_set = half2_set = 0;
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i) {
    int num_vtx, num_elem;
    rval = moab.get_number_entities_by_type( *i, MBVERTEX, num_vtx );
    CHECK_ERR(rval);
    rval = moab.get_number_entities_by_type( *i, MBHEX, num_elem );
    CHECK_ERR(rval);
    if (num_vtx == 0) {
      if (num_elem == 0) {
        CHECK_EQUAL( (MBEntityHandle)0, empty_set );
        empty_set = *i;
      }
      else if (!half1_set) {
        half1_set = *i;
      }
      else {
        CHECK_EQUAL( (MBEntityHandle)0, half2_set );
        half2_set = *i;
      }
    }
    else if (num_vtx == 1) {
      CHECK_EQUAL( 0, num_elem );
      CHECK_EQUAL( (MBEntityHandle)0, vertex_set );
      vertex_set = *i;
    }
    else {
      CHECK_EQUAL( 8, num_elem );
      CHECK_EQUAL( 27, num_vtx );
      CHECK_EQUAL( (MBEntityHandle)0, all_set );
      all_set = *i;
    }
  }
  
    // check set options
  unsigned opt;
  rval = moab.get_meshset_options( all_set, opt );
  CHECK_ERR(rval);
  CHECK_EQUAL( all_opt, opt );
  rval = moab.get_meshset_options( half1_set, opt );
  CHECK_ERR(rval);
  CHECK_EQUAL( half1_opt, opt );
  rval = moab.get_meshset_options( half2_set, opt );
  CHECK_ERR(rval);
  CHECK_EQUAL( half2_opt, opt );
  rval = moab.get_meshset_options( vertex_set, opt );
  CHECK_ERR(rval);
  CHECK_EQUAL( vertex_opt, opt );
  rval = moab.get_meshset_options( empty_set, opt );
  CHECK_ERR(rval);
  CHECK_EQUAL( empty_opt, opt );
}

void test_pack_set_contents()
{
  MBCore moab;
  MBErrorCode rval;
  create_simple_grid( moab, 3 );  

    // delete any existing sets
  MBRange sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
  
    // get all vertices
  MBRange vertices;
  rval = moab.get_entities_by_type( 0, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK_EQUAL( (MBEntityHandle)27, vertices.size() );
    // create meshset containing vertices
  MBEntityHandle set;
  rval = moab.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  rval = moab.add_entities( set, vertices );
  CHECK_ERR(rval);
  
    // pack and unpack range containing only set handle.
    // Will fail unless we also pass in set contents explicitly
  MBRange entities( vertices );
  entities.insert( set );
  pack_unpack_noremoteh( moab, entities );
  
    // expect single set in mesh
  entities.clear();
  rval = moab.get_entities_by_type( 0, MBENTITYSET, entities );
  CHECK_ERR(rval);
  CHECK_EQUAL( (MBEntityHandle)1, entities.size() );
  set = entities.front();
  
    // expect 27 vertices in mesh
  vertices.clear();
  rval = moab.get_entities_by_type( 0, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK_EQUAL( (MBEntityHandle)27, vertices.size() );
  
    // expect set to contain all 27 vertices
  vertices.clear();
  rval = moab.get_entities_by_type( set, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK_EQUAL( (MBEntityHandle)27, vertices.size() );
}

void test_pack_sets_of_sets()
{
  MBCore moab;
  MBErrorCode rval;
 
    // delete any existing sets
  MBRange sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
 
    // create three sets such that set2 contains set1, and set3 contains
    // both set1 and set2
  MBEntityHandle set1, set2, set3;
  sets.clear();
  rval = moab.create_meshset( MESHSET_ORDERED, set1 );
  CHECK_ERR(rval);
  sets.insert( set1 );
  rval = moab.create_meshset( MESHSET_SET, set2 );
  CHECK_ERR(rval);
  rval = moab.add_entities( set2, sets );
  CHECK_ERR(rval);
  sets.insert( set2 );
  rval = moab.create_meshset( MESHSET_SET, set3 );
  CHECK_ERR(rval);
  rval = moab.add_entities( set3, sets );
  CHECK_ERR(rval);
  sets.insert( set3 );
  
    // pack and unpack
  pack_unpack_noremoteh( moab, sets );
  
    // get sets
  sets.clear();
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( (MBEntityHandle)3, sets.size() );
  
    // figure out which is which
  set1 = set2 = set3 = 0;
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i) {
    int count;
    rval = moab.get_number_entities_by_type( *i, MBENTITYSET, count );
    CHECK_ERR(rval);
    CHECK( count >= 0 && count <= 2 );
    switch (count) {
      case 0:
        CHECK_EQUAL( (MBEntityHandle)0, set1 );
        set1 = *i;
        break;
      case 1:
        CHECK_EQUAL( (MBEntityHandle)0, set2 );
        set2 = *i;
        break;
      case 2:
        CHECK_EQUAL( (MBEntityHandle)0, set3 );
        set3 = *i;
        break;
    }
  }
  
    // check that set2 contains set1
  sets.clear();
  rval = moab.get_entities_by_type( set2, MBENTITYSET, sets );
  CHECK_EQUAL( (MBEntityHandle)1, sets.size() );
  CHECK_EQUAL( set1, sets.front() );
  
    // check that set3 contains set1 and set2
  sets.clear();
  rval = moab.get_entities_by_type( set3, MBENTITYSET, sets );
  CHECK_EQUAL( (MBEntityHandle)2, sets.size() );
  if (sets.front() == set1) {
    CHECK_EQUAL( set2, sets.back() );
  }
  else {
    CHECK_EQUAL( set2, sets.front() );
    CHECK_EQUAL( set1, sets.back() );
  }
}

void test_pack_set_parent_child()
{
  MBCore moab;
  MBErrorCode rval;
 
    // delete any existing sets
  MBRange sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
 
    // create three sets such that set3 has a child link to
    // set1, set2 has a parent link to set1, and such that set3 and
    // set2 are parent and child, respectively.
  MBEntityHandle set1, set2, set3;
  sets.clear();
  rval = moab.create_meshset( MESHSET_ORDERED, set1 );
  CHECK_ERR(rval);
  sets.insert( set1 );
  rval = moab.create_meshset( MESHSET_SET, set2 );
  CHECK_ERR(rval);
  sets.insert( set2 );
  rval = moab.create_meshset( MESHSET_SET, set3 );
  CHECK_ERR(rval);
  sets.insert( set3 );
  
  rval = moab.add_child_meshset( set3, set1 );
  CHECK_ERR(rval);
  rval = moab.add_parent_meshset( set2, set1 );
  CHECK_ERR(rval);
  rval = moab.add_parent_child( set3, set2 );
  CHECK_ERR(rval);
  
    // make sure everything is valid before doing the pack/unpack
  int count;
  rval = moab.num_child_meshsets( set1, &count );
  CHECK( MB_SUCCESS == rval && 0 == count );
  rval = moab.num_child_meshsets( set2, &count );
  CHECK( MB_SUCCESS == rval && 0 == count );
  rval = moab.num_child_meshsets( set3, &count );
  CHECK( MB_SUCCESS == rval && 2 == count );
  rval = moab.num_parent_meshsets( set1, &count );
  CHECK( MB_SUCCESS == rval && 0 == count );
  rval = moab.num_parent_meshsets( set2, &count );
  CHECK( MB_SUCCESS == rval && 2 == count );
  rval = moab.num_parent_meshsets( set3, &count );
  CHECK( MB_SUCCESS == rval && 0 == count );
  
    // pack and unpack
  pack_unpack_noremoteh( moab, sets );
  
    // get sets
  sets.clear();
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( (MBEntityHandle)3, sets.size() );
  
    // look for a set with two child links (set3)
  set1 = set2 = set3 = 0;
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i) {
    int count;
    rval = moab.num_child_meshsets( *i, &count );
    CHECK_ERR(rval);
    if (count == 2) {
      set3 = *i;
      break;
    }
  }
  CHECK( 0 != set3 );
  
    // check set relations
  std::vector<MBEntityHandle> parents, children;
  rval = moab.get_child_meshsets( set3, children );
  CHECK_ERR(rval);
  CHECK_EQUAL( (std::vector<MBEntityHandle>::size_type)2, children.size() );
  set1 = children[0];
  set2 = children[1];
  rval = moab.get_parent_meshsets( set1, parents );
  CHECK_ERR(rval);
  CHECK( parents.empty() );
  children.clear();
  rval = moab.get_parent_meshsets( set1, children );
  CHECK_ERR(rval);
  CHECK( children.empty() );
  rval = moab.get_parent_meshsets( set2, parents );
  CHECK_ERR(rval);
  CHECK_EQUAL( (std::vector<MBEntityHandle>::size_type)2, parents.size() );
  CHECK_EQUAL( set1, parents[0] );
  CHECK_EQUAL( set3, parents[1] );
}

void test_pack_tag_data_sparse()
{
  MBRange::iterator i;
  int size;
  MBDataType type;
  MBTagType storage;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
  create_simple_grid( mb, 3 );  
  MBRange elems;
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  CHECK( !elems.empty() );
 
    // Define a sparse tag containing two integers.  For every other element
    // in the mesh, set the tag value to the floor of the x and y
    // coordinates of the first vertex in the elements connectivity list.
  const char sparse_2_int_tag_name[] = "test tag 1";
  MBTag sparse_2_int_tag;
  rval = mb.tag_create( sparse_2_int_tag_name,
                        2*sizeof(int),
                        MB_TAG_SPARSE,
                        MB_TYPE_INTEGER,
                        sparse_2_int_tag,
                        0 );
  CHECK_ERR(rval);
  bool skip = false;
  for (i = elems.begin(); i != elems.end(); ++i, skip = !skip) {
    if (skip)
      continue;
    
    const MBEntityHandle* conn =0;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    double coords[3];
    rval = mb.get_coords( conn, 1, coords );
    CHECK_ERR(rval);
    const int data[2] = { (int)coords[0], (int)coords[1] };
    rval = mb.tag_set_data( sparse_2_int_tag, &*i, 1, data );
    CHECK_ERR(rval);
  }
  
    // pack and unpack
  MBRange ents;
  pack_unpack_noremoteh( moab, ents );
  elems.clear();
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  
  
    // check tag meta for sparse_2_int_tag
  rval = mb.tag_get_handle( sparse_2_int_tag_name, sparse_2_int_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_size( sparse_2_int_tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( (int)(2*sizeof(int)), size );
  rval = mb.tag_get_type( sparse_2_int_tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_SPARSE, storage );
  rval = mb.tag_get_data_type( sparse_2_int_tag, type );
  CHECK_EQUAL( MB_TYPE_INTEGER, type );
  int intdata[2];
  rval = mb.tag_get_default_value( sparse_2_int_tag, intdata );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // check tag data for sparse_2_int_tag
  MBRange tagged;
  rval = mb.get_entities_by_type_and_tag( 0, MBHEX, 
                                          &sparse_2_int_tag, 0, 1,
                                          tagged );
  CHECK_ERR(rval);
  CHECK_EQUAL( (elems.size()+1) / 2, tagged.size() );
  for (i = tagged.begin(); i != tagged.end(); ++i) {
    rval = mb.tag_get_data( sparse_2_int_tag, &*i, 1, intdata );
    CHECK_ERR(rval);
    
    const MBEntityHandle* conn =0;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    double coords[3];
    rval = mb.get_coords( conn, 1, coords );
    CHECK_ERR(rval);
    
    CHECK_EQUAL( (int)(coords[0]), intdata[0] );
    CHECK_EQUAL( (int)(coords[1]), intdata[1] );
  }
}


void test_pack_tag_data_dense()
{
  MBRange::iterator i;
  int size;
  MBDataType type;
  MBTagType storage;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
  create_simple_grid( mb, 3 );  
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );

    // Define a dense tag containing a single double-precision floating
    // point value.  For each vertex, store the distance from the origin
    // in this tag.
  const char dense_1_double_tag_name[] = "test tag 2";
  MBTag dense_1_double_tag;
  rval = mb.tag_create( dense_1_double_tag_name,
                        sizeof(double),
                        MB_TAG_DENSE,
                        MB_TYPE_DOUBLE,
                        dense_1_double_tag,
                        0 );
  CHECK_ERR(rval);
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    double val = sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
    rval = mb.tag_set_data( dense_1_double_tag, &*i, 1, &val );
    CHECK_ERR(rval);
  }
  
    // pack and unpack
  MBRange ents;
  pack_unpack_noremoteh( moab, ents );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  
    // check tag meta for dense_1_double_tag
  rval = mb.tag_get_handle( dense_1_double_tag_name, dense_1_double_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_size( dense_1_double_tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( (int)sizeof(double), size );
  rval = mb.tag_get_type( dense_1_double_tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  rval = mb.tag_get_data_type( dense_1_double_tag, type );
  CHECK_EQUAL( MB_TYPE_DOUBLE, type );
  double dval;
  rval = mb.tag_get_default_value( dense_1_double_tag, &dval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // check tag data for dense_1_double_tag
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    const double expected = sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]);
    rval = mb.tag_get_data( dense_1_double_tag, &*i, 1, &dval );
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL( expected, dval, 1e-6 );
  }
}

void test_pack_tag_data_default_value()
{
  MBRange::iterator i;
  int size;
  MBDataType type;
  MBTagType storage;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
  create_simple_grid( mb, 3 );  
  MBRange verts, elems, sets;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  CHECK( !elems.empty() );

    // Define a dense, opaque tag with a default value of "DEFLT".
    // Set the tag on one element, one vertex,and one set to "TAGGD".
  const char dense_5_opaque_tag_name[] = "This is intentionally a very long tag name in an attempt to test for an arbitrary limitations on tag name length.";
  MBTag dense_5_opaque_tag;
  rval = mb.tag_create( dense_5_opaque_tag_name,
                        5,
                        MB_TAG_DENSE,
                        dense_5_opaque_tag,
                        "DEFLT" );
  CHECK_ERR(rval);
  MBEntityHandle set;
  rval = mb.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  const MBEntityHandle handles[3] = { verts.front(), elems.front(), set };
  const char data[] = "TAGGDTAGGDTAGGD";
  rval = mb.tag_set_data( dense_5_opaque_tag, handles, 3, data );
  CHECK_ERR(rval);
  
    // pack and unpack
  MBRange ents;
  pack_unpack_noremoteh( moab, ents );
  elems.clear();
  verts.clear();
  sets.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  
    // check tag meta for dense_5_opaque_tag
  rval = mb.tag_get_handle( dense_5_opaque_tag_name, dense_5_opaque_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_size( dense_5_opaque_tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 5, size );
  rval = mb.tag_get_type( dense_5_opaque_tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  rval = mb.tag_get_data_type( dense_5_opaque_tag, type );
  CHECK_EQUAL( MB_TYPE_OPAQUE, type );
  char odata[6]; odata[5] = '\0';
  rval = mb.tag_get_default_value( dense_5_opaque_tag, odata );
  CHECK_ERR( rval );
  CHECK_EQUAL( std::string("DEFLT"), std::string(odata) );
  
    // count number of each type with tag set to non-default
  int vcount = 0, ecount = 0, scount =0;
  for (i = verts.begin(); i != verts.end(); ++i) {
    rval = mb.tag_get_data( dense_5_opaque_tag, &*i, 1, odata );
    CHECK_ERR(rval);
    if (strcmp( odata, "DEFLT" )) {
      CHECK_EQUAL( std::string("TAGGD"), std::string(odata) );
      ++vcount;
    }
  }
  CHECK_EQUAL( 1, vcount );
  for (i = elems.begin(); i != elems.end(); ++i) {
    rval = mb.tag_get_data( dense_5_opaque_tag, &*i, 1, odata );
    CHECK_ERR(rval);
    if (strcmp( odata, "DEFLT" )) {
      CHECK_EQUAL( "TAGGD", odata );
      ++ecount;
    }
  }
  CHECK_EQUAL( 1, ecount );
  for (i = sets.begin(); i != sets.end(); ++i) {
    rval = mb.tag_get_data( dense_5_opaque_tag, &*i, 1, odata );
    CHECK_ERR(rval);
    if (strcmp( odata, "DEFLT" )) {
      CHECK_EQUAL( "TAGGD", odata );
      ++scount;
    }
  }
  CHECK_EQUAL( 1, scount );
}

void test_pack_bit_tag_data()
{
  MBRange::iterator i;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
 
    // Create a bit tag
  const char tag_name[] = "test bit";
  MBTag tag;
  rval = mb.tag_create( tag_name, 3, MB_TAG_BIT, tag, 0 );
  CHECK_ERR(rval);
  
    // Set bits to 1 unless cooresponding coordinate of 
    // vertex is zero.
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    unsigned char data = 0;
    for (int j = 0; j < 3; ++j) 
      if (fabs(coords[j]) > 1e-6)
        data |= (1 << j);
    rval = mb.tag_set_data( tag, &*i, 1, &data );
    CHECK_ERR(rval);
  }
  
    // pack and unpack
  MBRange ents;
  pack_unpack_noremoteh( moab, ents );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);

    // check tag meta 
  rval = mb.tag_get_handle( tag_name, tag );
  CHECK_ERR(rval);
  
  int size;
  rval = mb.tag_get_size( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, size );
  
  MBTagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_BIT, storage );
  
  MBDataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_EQUAL( MB_TYPE_BIT, type );

    // check tag values
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    unsigned char expected = 0;
    for (int j = 0; j < 3; ++j) 
      if (fabs(coords[j]) > 1e-6)
        expected |= (1 << j);
        
    unsigned char data = (unsigned char)0xFF;
    rval = mb.tag_get_data( tag, &*i, 1, &data );
    CHECK_ERR(rval);
    
    CHECK_EQUAL( (int)expected, (int)data );
  }
}

void test_pack_variable_length_tag()
{
  MBRange::iterator i;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  
    // create a variable-length tag 
  const char* tag_name = "var_int_tag";
  const int default_val[] = { 0xBEEF, 0xFEED, 0xDEAD, 0xBAD, 0xBEAD };
  MBTag tag;
  rval = mb.tag_create_variable_length( tag_name, MB_TAG_DENSE, MB_TYPE_INTEGER,
                                        tag, default_val, sizeof(default_val) );
  CHECK_ERR(rval);
  
    // for each vertex, store in the tag an integer between 1 and 3, 
    // followed by the floor of the cooresponding number of vertex
    // coordinates, beginning with x.
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    const int num_coord = 1 + *i % 3;
    const int data_size = sizeof(int) * (num_coord + 1);
    const int data[4] = { num_coord, (int)coords[0], (int)coords[1], (int)coords[2] };
    const void* data_ptrs[1] = { data };
    rval = mb.tag_set_data( tag, &*i, 1, data_ptrs, &data_size );
    CHECK_ERR(rval);
  }
  
    // pack and unpack
  pack_unpack_noremoteh( moab );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);

    // check tag meta 
  rval = mb.tag_get_handle( tag_name, tag );
  CHECK_ERR(rval);
  
  int size;
  rval = mb.tag_get_size( tag, size );
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, rval );
  
  MBTagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  
  MBDataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_EQUAL( MB_TYPE_INTEGER, type );

  const void* defval_ptr;
  rval = mb.tag_get_default_value( tag, defval_ptr, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( sizeof(default_val), (size_t)size );
  const int* defval_arr = reinterpret_cast<const int*>(defval_ptr);
  for (size_t j = 0; j < sizeof(default_val)/sizeof(int); ++j)
    CHECK_EQUAL( default_val[j], defval_arr[j] );
  
    // check tag values
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    int size;
    const void* valptr;
    rval = mb.tag_get_data( tag, &*i, 1, &valptr, &size );
    CHECK_ERR(rval);
    CHECK( (size % sizeof(int)) == 0 );
    int count = size / sizeof(int);
    CHECK( count > 1 );
    CHECK( count <= 4 );
    
    const int* valarr = reinterpret_cast<const int*>(valptr);
    CHECK( valarr[0] >= 1 );
    CHECK( valarr[0] <= 3 );
    for (int j = 0; j < valarr[0]; ++j) {
      CHECK_EQUAL( (int)(coords[j]), valarr[j+1] );
    }
  }
}

void test_pack_tag_handle_data()
{
  MBRange::iterator i;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  MBRange verts, elems;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  CHECK( !elems.empty() );
 
    // create a tag
  const char* tag_name = "entity tag";
  MBEntityHandle default_val[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  MBTag tag;
  rval = mb.tag_create( tag_name, 8*sizeof(MBEntityHandle), MB_TAG_SPARSE, MB_TYPE_HANDLE, tag, &default_val );
  CHECK_ERR(rval);
  
    // Store on each vertex the handles of the adjacent hexes, padded
    // with NULL handles.
  MBEntityHandle tagdata[8*8];
  for (i = elems.begin(); i != elems.end(); ++i) {
    const MBEntityHandle* conn;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 8, len );
    
    rval = mb.tag_get_data( tag, conn, len, tagdata );
    CHECK_ERR(rval);
    
    for (int j = 0; j < 8; ++j) {
      MBEntityHandle* vdata = tagdata + 8*j;
      int idx = 0;
      while (vdata[idx]) {
        ++idx;
        CHECK(idx < 8);
      }
      vdata[idx] = *i;
    }
     
    rval = mb.tag_set_data( tag, conn, len, tagdata );
    CHECK_ERR(rval);
  }
  
    // pack and unpack
  pack_unpack_noremoteh( moab );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);

    // check tag meta 
  rval = mb.tag_get_handle( tag_name, tag );
  CHECK_ERR(rval);
  
  int size;
  rval = mb.tag_get_size( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8*sizeof(MBEntityHandle), (size_t)size );
  
  MBTagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_SPARSE, storage );
  
  MBDataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_EQUAL( MB_TYPE_HANDLE, type );

  rval = mb.tag_get_default_value( tag, tagdata );
  CHECK_ERR(rval);
  for (int j = 0; j < 8; ++j) {
    CHECK_EQUAL( (MBEntityHandle)0, tagdata[j] );
  }
  
    // check tag values
  for (i = elems.begin(); i != elems.end(); ++i) {
    const MBEntityHandle* conn;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 8, len );
    
    rval = mb.tag_get_data( tag, conn, len, tagdata );
    CHECK_ERR(rval);
    
    for (int j = 0; j < 8; ++j) {
      MBEntityHandle* vdata = tagdata + 8*j;
      int idx = 0;
      while (vdata[idx] != *i) {
        ++idx;
        CHECK(idx < 8);
      }
      vdata[idx] = 0;
    }
     
    rval = mb.tag_set_data( tag, conn, len, tagdata );
    CHECK_ERR(rval);
  }
  
  for (i = verts.begin(); i != verts.end(); ++i) {
    rval = mb.tag_get_data( tag, &*i, 1, tagdata );
    CHECK_ERR(rval);
    for (int j = 0; j < 8; ++j) {
      CHECK_EQUAL( (MBEntityHandle)0, tagdata[j] );
    }
  }
}
  
MBErrorCode get_entities(MBInterface *mb,
                         std::vector<MBEntityHandle> &ent_verts, 
                         int verts_per_entity, int dim, 
                         MBRange &ents) 
{
  assert(!(ent_verts.size()%verts_per_entity));
  unsigned int num_ents = ent_verts.size() / verts_per_entity;
  MBRange dum_ents;
  MBErrorCode result;
  for (unsigned int i = 0; i < num_ents; i++) {
    result = mb->get_adjacencies(&ent_verts[verts_per_entity*i], verts_per_entity,
                                dim, true, dum_ents);
    CHECK_ERR(result);
    assert(dum_ents.size() == 1);
    ents.merge(dum_ents);
    dum_ents.clear();
  }
  return MB_SUCCESS;
}
  
void test_pack_shared_entities()
{
  MBCore moab[4];
  MBInterface *mb[] = {&moab[0], &moab[1], &moab[2], &moab[3]};
  MBParallelComm *pc[4];
  for (unsigned int i = 0; i < 4; i++) {
    pc[i] = new MBParallelComm(&moab[i]);
    pc[i]->set_rank(i);
    for (unsigned int j = 0; j < 4; j++) {
      if (j == i) continue;
      else pc[i]->get_buffers(j);
    }
  }

  MBRange verts[4], quads[4];
  MBErrorCode rval = create_shared_grid(pc, verts, quads);

  std::vector<std::vector<MBEntityHandle> > L1h[4];
  std::vector<std::vector<int> > L1p[4];
  for (unsigned int i = 0; i < 4; i++) {
    L1h[i].resize(4);
    L1p[i].resize(4);
  }
  
  std::vector<MBEntityHandle> L2hloc[4];
  std::vector<MBEntityHandle> L2hrem[4];
  std::vector<unsigned int> L2p[4];
  MBRange new_ents[4];
  MBRange sent_ents;
  std::vector<MBEntityHandle> ent_verts, dum_handles;
  std::vector<int> dum_procs;
  std::vector<std::set<unsigned int> > entprocs(4);
  std::vector<unsigned char> buffer;
  buffer.reserve(1);
  unsigned char *buff_ptr = &buffer[0];

    //========================
    // interface, shared, 0->2
    // get edges that we want
  ent_verts.push_back(verts[0][6]); ent_verts.push_back(verts[0][7]);
  ent_verts.push_back(verts[0][7]); ent_verts.push_back(verts[0][8]);
  rval = get_entities(mb[0], ent_verts, 2, 1, sent_ents);
  CHECK_ERR(rval);
  assert(2 == sent_ents.size());
    // set entprocs
  entprocs[0].insert(2); entprocs[1].insert(2);
  
  rval = pc[0]->pack_entities(sent_ents, buffer, buff_ptr, true, 2,
                              true, &entprocs);
  CHECK_ERR(rval);

    // now unpack the buffer
  buff_ptr = &buffer[0];
  rval = pc[2]->unpack_entities(buff_ptr, true, 0, true,
                                L1h[2], L1p[2], L2hloc[2], L2hrem[2], L2p[2], new_ents[2]);
    // all L1h lists should be empty, since we're dealing with iface
  if (!L1h[2][0].empty() || !L1h[2][1].empty() || 
      !L1h[2][2].empty() || !L1h[2][3].empty() || 
      !L2p[2].empty() || !new_ents[2].empty()) rval = MB_FAILURE;
  CHECK_ERR(rval);

  buffer.clear(); buff_ptr = &buffer[0];
  ent_verts.clear(); sent_ents.clear();

    //========================
    // interface, multishared, 1st message, 0->1
  ent_verts.push_back(verts[0][8]); ent_verts.push_back(verts[0][5]);
  ent_verts.push_back(verts[0][5]); ent_verts.push_back(verts[0][2]);
  rval = get_entities(mb[0], ent_verts, 2, 1, sent_ents);
  CHECK_ERR(rval);
  assert(2 == sent_ents.size());
    // sending these edges to 1 and 3
  entprocs[0].insert(1); entprocs[0].insert(3);
  entprocs[1].insert(1); entprocs[1].insert(3);

  rval = pc[0]->pack_entities(sent_ents, buffer, buff_ptr, true, 1,
                              true, &entprocs);
  CHECK_ERR(rval);

    // now unpack the buffer
  buff_ptr = &buffer[0];
  rval = pc[1]->unpack_entities(buff_ptr, true, 0, true,
                                L1h[1], L1p[1], L2hloc[1], L2hrem[1], L2p[1], new_ents[1]);
    // all L1h lists should be empty, since we're dealing with iface
  if (!L1h[1][0].empty() || !L1h[1][1].empty() || 
      !L1h[1][2].empty() || !L1h[1][3].empty() || 
      !L2p[1].empty() || !new_ents[1].empty()) rval = MB_FAILURE;
  CHECK_ERR(rval);

  entprocs[0].clear(); entprocs[1].clear();
  buffer.clear(); buff_ptr = &buffer[0];
  ent_verts.clear(); sent_ents.clear();
  
    //========================
    // interface, multishared, 2nd message, 3->1
    // sending these edges to 0 and 1
  ent_verts.push_back(verts[3][6]); ent_verts.push_back(verts[3][3]);
  ent_verts.push_back(verts[3][3]); ent_verts.push_back(verts[3][0]);
  rval = get_entities(mb[3], ent_verts, 2, 1, sent_ents);
  CHECK_ERR(rval);
  assert(2 == sent_ents.size());
  entprocs[0].insert(0); entprocs[0].insert(1);
  entprocs[1].insert(0); entprocs[1].insert(1);
  rval = pc[3]->pack_entities(sent_ents, buffer, buff_ptr, true, 1,
                              true, &entprocs);
  CHECK_ERR(rval);

    // now unpack the buffer
  buff_ptr = &buffer[0];
  rval = pc[1]->unpack_entities(buff_ptr, true, 3, true,
                                L1h[1], L1p[1], L2hloc[1], L2hrem[1], L2p[1], new_ents[1]);
    // all L1h lists should be empty, since we're dealing with iface
  if (!L1h[1][0].empty() || !L1h[1][1].empty() || 
      !L1h[1][2].empty() || !L1h[1][3].empty())
    rval = MB_FAILURE;
  CHECK_ERR(rval);

  ent_verts.clear(); sent_ents.clear();
  entprocs[0].clear(); entprocs[1].clear();
  buffer.clear(); buff_ptr = &buffer[0];
  
    //========================================================================
    // prepare for ghost communication; set all iface entities ownership
    // stuff, so we don't need to do it explicitly using pack/unpack;
    // passing true makes it happen for edges too
  rval = set_owners(pc, verts, quads, true);
  CHECK_ERR(rval);

    //========================
    // ghost, unshared, 2->1
  sent_ents.insert(quads[2][3]);  
    // add vertices not already shared
  sent_ents.insert(verts[2][4]); 
  sent_ents.insert(verts[2][7]); 
    // entprocs lists are all just 1
  entprocs.resize(sent_ents.size());
  entprocs[0].insert(1); 
  entprocs[1].insert(1); 
  entprocs[2].insert(1); 

  rval = pc[2]->pack_entities(sent_ents, buffer, buff_ptr, true, 1,
                              true, &entprocs);
  CHECK_ERR(rval);

    // now unpack the buffer
  buff_ptr = &buffer[0];
  rval = pc[1]->unpack_entities(buff_ptr, true, 2, false,
                                L1h[1], L1p[1], L2hloc[1], L2hrem[1], L2p[1], new_ents[1]);
  if (
        // to P2: 2 handles per 3 entities = 6
      L1h[1][1].size() != 6 || 
        // to all others, 0
      !L1h[1][0].empty() || !L1h[1][2].empty() || !L1h[1][3].empty()) 
    rval = MB_FAILURE;
  CHECK_ERR(rval);

  buffer.clear(); buff_ptr = &buffer[0];
  rval = pc[1]->pack_remote_handles(L1h[1][1], L1p[1][1], 2, buffer, buff_ptr);
  CHECK_ERR(rval);

  buff_ptr = &buffer[0];
  rval = pc[2]->unpack_remote_handles(1, buff_ptr, false, L2hloc[2], L2hrem[2], L2p[2]);
  CHECK_ERR(rval);

  ent_verts.clear(); sent_ents.clear();
  entprocs[0].clear(); entprocs[1].clear(); entprocs[2].clear();
  buffer.clear(); buff_ptr = &buffer[0]; L1h[1][1].clear();

    //========================
    // ghost, multishared, 1st message, 2->0
    // sent: v3, v4, v5, q0, q1
  sent_ents.insert(quads[2][0]); sent_ents.insert(quads[2][1]);  
    // add vertices not already shared
  sent_ents.insert(verts[2][3]); sent_ents.insert(verts[2][4]); 
  sent_ents.insert(verts[2][5]); 
    // entprocs lists are entity-dependent
    // v3, v4, v5, q1: only to 0 (v4, v5 already on 1)
  entprocs.resize(sent_ents.size());
  entprocs[0].insert(0); 
  entprocs[1].insert(0); 
  entprocs[2].insert(0); 
  entprocs[3].insert(0); 
    // q2: P0 and P1
  entprocs[4].insert(0); entprocs[4].insert(1); 

  rval = pc[2]->pack_entities(sent_ents, buffer, buff_ptr, true, 0,
                              true, &entprocs);
  CHECK_ERR(rval);

    // now unpack the buffer
  buff_ptr = &buffer[0];
  rval = pc[0]->unpack_entities(buff_ptr, true, 2, false,
                                L1h[0], L1p[0], L2hloc[0], L2hrem[0], L2p[0], new_ents[0]);
  if (
        // 2 handles per 5 entities = 10
      L1h[0][1].size() != 10 || 
        // 2 handles per 3 entities
      L1h[0][0].size() != 6 || 
        // none received from 3, no 4th proc
      !L1h[0][2].empty() || !L1h[0][3].empty()) 
    rval = MB_FAILURE;
  CHECK_ERR(rval);

  ent_verts.clear(); sent_ents.clear();
  for (int i = 0; i < 5; i++) entprocs[i].clear();
  L1p[0][1].clear();  L1h[0][1].clear();
  L1p[0][0].clear();  L1h[0][0].clear();
  buffer.clear(); buff_ptr = &buffer[0];

    //========================
    // ghost, multishared, 2nd message, 1->0
    // sent: v1, v4, v7, q0, q2
  sent_ents.insert(quads[1][0]); sent_ents.insert(quads[1][2]);
    // add vertices not already shared
  sent_ents.insert(verts[1][1]); sent_ents.insert(verts[1][4]); 
  sent_ents.insert(verts[1][7]); 
    // entprocs lists are entity-dependent
    // v1, q0: only to 0
  entprocs.resize(sent_ents.size());
  entprocs[0].insert(0); 
  entprocs[2].insert(0);
  entprocs[3].insert(0); 
    // v4, v7, q2: P0 and P2
  entprocs[1].insert(0); entprocs[1].insert(2); 
  entprocs[4].insert(0); entprocs[4].insert(2); 

  rval = pc[1]->pack_entities(sent_ents, buffer, buff_ptr, true, 0,
                              true, &entprocs);
  CHECK_ERR(rval);

    // now unpack the buffer
  buff_ptr = &buffer[0];
  rval = pc[0]->unpack_entities(buff_ptr, true, 1, false,
                                L1h[0], L1p[0], L2hloc[0], L2hrem[0], L2p[0], new_ents[0]);
  if (
        // 2 handles per 5 entities = 10
      L1h[0][0].size() != 10 || 
        // 2 handles per 3 entities
      L1h[0][1].size() != 6 || 
        // none received from 3, no 4th proc
      !L1h[0][2].empty() || !L1h[0][3].empty()) 
    rval = MB_FAILURE;
  CHECK_ERR(rval);

    //========================
    // now pack/unpack the handles
    // P0 -> P2
    // xxx moved pack/unpack handles to after P2, P1 send
  buffer.clear(); buff_ptr = &buffer[0];
  rval = pc[0]->pack_remote_handles(L1h[0][1], L1p[0][1], 2, buffer, buff_ptr);
  CHECK_ERR(rval);

  rval = pc[2]->unpack_remote_handles(0, buff_ptr, false, L2hloc[2], L2hrem[2], L2p[2]);
  CHECK_ERR(rval);

    // P0 -> P1
  buffer.clear(); buff_ptr = &buffer[0];
  rval = pc[0]->pack_remote_handles(L1h[0][0], L1p[0][0], 1, buffer, buff_ptr);
  CHECK_ERR(rval);

  rval = pc[1]->unpack_remote_handles(0, buff_ptr, false, L2hloc[1], L2hrem[1], L2p[1]);
  CHECK_ERR(rval);

  ent_verts.clear(); sent_ents.clear();
  for (int i = 0; i < 5; i++) entprocs[i].clear();
  L1p[0][1].clear();  L1h[0][1].clear();
  buffer.clear(); buff_ptr = &buffer[0];

}

void test_filter_pstatus()
{
  MBRange::iterator i;
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  std::vector<MBEntityHandle> verts;
  MBRange dum_vertsr, vertsr;
  rval = mb.get_entities_by_type( 0, MBVERTEX, dum_vertsr );
  CHECK_ERR(rval);
  vertsr.insert(dum_vertsr[0], dum_vertsr[8]);
  for (unsigned int i = 0; i < 9; i++) verts.push_back(vertsr[i]);

  CHECK( !verts.empty() );
 
  MBParallelComm *pcomm = new MBParallelComm( &moab );

  std::vector<int> procs(70, -1);
  for (unsigned int i = 0; i < 6; i++) procs[i] = i;

  std::vector<unsigned char> pvals(verts.size(), 0);
    // interface, owned
  pvals[0] = (PSTATUS_INTERFACE | PSTATUS_SHARED); // p0
  rval = moab.tag_set_data(pcomm->sharedp_tag(), &verts[0], 1, &procs[0]); CHECK_ERR(rval);  
    // interface, not owned
  pvals[1] = (PSTATUS_NOT_OWNED | PSTATUS_INTERFACE | PSTATUS_SHARED); // p1
  rval = moab.tag_set_data(pcomm->sharedp_tag(), &verts[1], 1, &procs[1]); CHECK_ERR(rval);  
    // interface, multi-shared, owned
  pvals[2] = (PSTATUS_INTERFACE | PSTATUS_SHARED | PSTATUS_MULTISHARED); // p0, p1
  rval = moab.tag_set_data(pcomm->sharedps_tag(), &verts[2], 1, &procs[0]); CHECK_ERR(rval);  
    // interface, multi-shared, not owned
  pvals[3] = (PSTATUS_INTERFACE | PSTATUS_MULTISHARED | PSTATUS_NOT_OWNED | PSTATUS_SHARED); // p1, p2
  rval = moab.tag_set_data(pcomm->sharedps_tag(), &verts[3], 1, &procs[1]); CHECK_ERR(rval);  
    // ghost, shared
  pvals[4] = (PSTATUS_GHOST | PSTATUS_SHARED | PSTATUS_NOT_OWNED); // p2
  rval = moab.tag_set_data(pcomm->sharedp_tag(), &verts[4], 1, &procs[2]); CHECK_ERR(rval);  
    // ghost, multi-shared
  pvals[5] = (PSTATUS_GHOST | PSTATUS_MULTISHARED | PSTATUS_NOT_OWNED | PSTATUS_SHARED); // p2, p3
  rval = moab.tag_set_data(pcomm->sharedps_tag(), &verts[5], 1, &procs[2]); CHECK_ERR(rval);  
    // owned, shared
  pvals[6] = (PSTATUS_SHARED); // p4
  rval = moab.tag_set_data(pcomm->sharedp_tag(), &verts[6], 1, &procs[4]); CHECK_ERR(rval);  
    // owned, multi-shared
  pvals[7] = (PSTATUS_MULTISHARED | PSTATUS_SHARED); // p4, p5
  rval = moab.tag_set_data(pcomm->sharedps_tag(), &verts[7], 1, &procs[4]); CHECK_ERR(rval);  
    // not shared, owned
  pvals[8] = 0x0;

  rval = moab.tag_set_data(pcomm->pstatus_tag(), &verts[0], 9, &pvals[0]);
  CHECK_ERR(rval);
  

  MBRange tmp_range = vertsr;

    // interface ents
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_INTERFACE, PSTATUS_AND);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 4 && *tmp_range.begin() == verts[0] && 
        *tmp_range.rbegin() == verts[3]);
    // not interface
  tmp_range = vertsr;
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_INTERFACE, PSTATUS_NOT);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 5 && *tmp_range.begin() == verts[4] && 
        *tmp_range.rbegin() == verts[8]);
    // interface not owned
  tmp_range = vertsr;
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_INTERFACE | PSTATUS_NOT_OWNED, PSTATUS_AND);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 2 && *tmp_range.begin() == verts[1] && 
        *tmp_range.rbegin() == verts[3]);
    // ghost
  tmp_range = vertsr;
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_GHOST, PSTATUS_AND);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 2 && *tmp_range.begin() == verts[4] && 
        *tmp_range.rbegin() == verts[5]);
    // shared not multi-shared
  tmp_range = vertsr;
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_SHARED, PSTATUS_AND);
  CHECK_ERR(rval);
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_MULTISHARED, PSTATUS_NOT);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 4 && tmp_range[0] == verts[0] && 
        tmp_range[1] == verts[1] && tmp_range[2] == verts[4] && tmp_range[3] == verts[6]);
    // shared w/ p0
  tmp_range = vertsr;
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_SHARED, PSTATUS_AND, 0);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 2 && tmp_range[1] == verts[2]);
    // shared w/ p2 && not owned
  tmp_range = vertsr;
  rval = pcomm->filter_pstatus(tmp_range, PSTATUS_SHARED | PSTATUS_NOT_OWNED, PSTATUS_AND, 2);
  CHECK_ERR(rval);
  CHECK(tmp_range.size() == 3 && tmp_range[0] == verts[3] && 
        tmp_range[1] == verts[4] && tmp_range[2] == verts[5]);
  
}
