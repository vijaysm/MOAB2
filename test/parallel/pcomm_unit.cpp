#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "TestUtil.hpp"
#include <algorithm>
#include <vector>
#include <set>
#include <sstream>

using namespace moab;

#ifdef USE_MPI
#  include "moab_mpi.h"
#endif

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

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
/** Test pack/unpack of shared entities in 2d*/
void test_pack_shared_entities_2d();
/** Test pack/unpack of shared entities in 3d*/
void test_pack_shared_entities_3d();
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
  num_err += RUN_TEST( test_pack_shared_entities_2d );
  num_err += RUN_TEST( test_pack_shared_entities_3d );
  num_err += RUN_TEST( test_filter_pstatus );
  
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return num_err;
}

/* Utility method: pack mesh data, clear moab instance, and unpack */
void pack_unpack_noremoteh( Core& moab, Range& entities )
{
  ErrorCode rval;
  if (entities.empty()) {
    rval = moab.get_entities_by_handle( 0, entities );
    CHECK_ERR(rval);
  }
  
  ParallelComm *pcomm = new ParallelComm( &moab );
  std::vector<int> addl_procs;

    // get the necessary vertices too
  Range tmp_range = entities.subset_by_type(MBENTITYSET);
  entities = subtract( entities, tmp_range);
  rval = moab.get_adjacencies(entities, 0, false, entities, Interface::UNION);
  CHECK_ERR(rval);
  entities.merge(tmp_range);
  
  ParallelComm::Buffer buff(ParallelComm::INITIAL_BUFF_SIZE);
  buff.reset_ptr(sizeof(int));
  rval = pcomm->pack_buffer( entities, false, true, false, 
                             -1, &buff);
  CHECK_ERR(rval);
  buff.set_stored_size();
  
  delete pcomm;
  moab.~Core();

  new (&moab) Core();
  pcomm = new ParallelComm( &moab);
  
  entities.clear();
  std::vector<std::vector<EntityHandle> > L1hloc, L1hrem;
  std::vector<std::vector<int> > L1p;
  std::vector<EntityHandle> L2hloc, L2hrem;
  std::vector<unsigned int> L2p;
  buff.reset_ptr(sizeof(int));
  std::vector<EntityHandle> entities_vec(entities.size());
  std::copy(entities.begin(), entities.end(), entities_vec.begin());
  rval = pcomm->unpack_buffer(buff.buff_ptr, false, -1, -1, L1hloc, L1hrem, L1p, L2hloc, 
                              L2hrem, L2p, entities_vec);
  CHECK_ERR(rval);
  std::copy(entities_vec.begin(), entities_vec.end(), range_inserter(entities));

  delete pcomm;
}
void pack_unpack_noremoteh( Core& moab )
{
  Range empty;
  pack_unpack_noremoteh( moab, empty );
}

/* Utility method -- check expected sizes */
void check_sizes( Interface& moab,
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
  ErrorCode rval;
  
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
void create_simple_grid( Interface& moab, unsigned x, unsigned y, unsigned z );
void create_simple_grid( Interface& moab, unsigned xyz = 3 )
{
  create_simple_grid( moab, xyz, xyz, xyz );
}
void create_simple_grid( Interface& moab, unsigned x, unsigned y, unsigned z )
{
  ErrorCode rval;
  EntityHandle *verts = new EntityHandle[x*y*z];
  for (unsigned k = 0; k < z; ++k)
    for (unsigned j = 0; j < y; ++j)
      for (unsigned i = 0; i < x; ++i) {
        const double coords[3] = { i, j, k };
        rval = moab.create_vertex( coords, verts[x*y*k + x*j + i] );
        CHECK_ERR(rval);
      }
  
  EntityHandle *elems = new EntityHandle[(x-1)*(y-1)*(z-1)];
  for (unsigned k = 0; k < (z-1); ++k)
    for (unsigned j = 0; j < (y-1); ++j)
      for (unsigned i = 0; i < (x-1); ++i) {
        const size_t idx = (size_t)i + (size_t)j*x + (size_t)k*x*y; 
        const EntityHandle conn[8] = { verts[idx        ],
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

ErrorCode create_patch(Interface *moab, Range &verts, Range &quads,
                         unsigned int n, double *xyz, int *gids) 
{
    // create vertices/quads in square array
  ErrorCode result = moab->create_vertices(xyz, n*n, verts);
  if (MB_SUCCESS != result) return result;
  std::vector<EntityHandle> connect;
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
    EntityHandle dum_quad;
    result = moab->create_element(MBQUAD, &connect[4*i], 4, dum_quad);
    if (MB_SUCCESS != result) return result;
    quads.insert(dum_quad);
  }
  
    // global ids
  Tag gid_tag;
  int dum_default = 0;
  result = moab->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag,
                                MB_TAG_DENSE|MB_TAG_CREAT, &dum_default);
  if (MB_SUCCESS != result) return result;
  result = moab->tag_set_data(gid_tag, verts, gids);
  if (MB_SUCCESS != result) return result;

  return result;
}

ErrorCode create_shared_grid_2d(ParallelComm **pc, Range *verts, Range *quads) 
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

  ErrorCode rval = ParallelComm::resolve_shared_ents(pc, 4, 2);
  CHECK_ERR(rval);

  return rval;
}

ErrorCode create_shared_grid_3d(ParallelComm **pc, Range *verts, Range *hexes) 
{
//    
//   4 _____   _____ 
//    |  |  | |  |  |
//   3|__|__| |__|__|          GIDS               HANDLES
//    |P1|  | |P2|  |
//    |__|__| |__|__|   20 21 22   22 23 24      7  8  9    7  8  9
//    2 ___________     15 16 17   17 18 19      4  5  6    4  5  6
//     |  |  |  |  |    10 11 12   12 13 14      1  2  3    1  2  3
// /  1|__|__|__|__|         
// |   |  |P0|  |  |     10 11 12 13 14           10 11 12 13 14     
// J  0|__|__|__|__|      5  6  7  8  9            5  6  7  8  9
// I-> 0  1  2  3  4      0  1  2  3  4            0  1  2  3  4
//
//  P3 - k = 2..4
  
    // create structured meshes
    // ijkmin[p][ijk], ijkmax[p][ijk]
#define P 4
  int ijkmin[P][3] = { {0, 0, 0}, {0, 2, 0}, {2, 2, 0}, {0, 0, 2}};
  int ijkmax[P][3] = { {4, 2, 2}, {2, 4, 2}, {4, 4, 2}, {4, 4, 4}};

  int nijk[P][3];
  int NIJK[3] = {0, 0, 0};
#define INDEXG(i, j, k) (k * NIJK[1] * NIJK[0] + j * NIJK[0] + i)
#define INDEXL(i, j, k) ((k-ijkmin[p][2])*nijk[p][1]*nijk[p][0] + \
                         (j-ijkmin[p][1])*nijk[p][0] + (i - ijkmin[p][0]))

  int p, i, j, k;
  for (int p = 0; p < P; p++) {
    for (int i = 0; i < 3; i++) {
      nijk[p][i] =  ijkmax[p][i] - ijkmin[p][i] + 1;
      NIJK[i] = std::max(NIJK[i], nijk[p][i]);
    }
  }
    
  std::vector<int> gids;
  std::vector<double> xyz;
  ErrorCode rval;
  Tag gid_tag;
  int dum_default = -1;

  for (p = 0; p < P; p++) {
    rval = pc[p]->get_moab()->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER,
                                             gid_tag, MB_TAG_DENSE|MB_TAG_CREAT,
                                             &dum_default);
    if (MB_SUCCESS != rval) return rval;

      // make vertices
    int nverts = nijk[p][0] * nijk[p][1] * nijk[p][2];
    xyz.resize(3*nverts);
    gids.resize(nverts);

      // set vertex gids
    int nv = 0;
    for (k = ijkmin[p][2]; k <= ijkmax[p][2]; k++) 
      for (j = ijkmin[p][1]; j <= ijkmax[p][1]; j++) 
        for (i = ijkmin[p][0]; i <= ijkmax[p][0]; i++) {
            // xyz
          xyz[3*nv] = i;
          xyz[3*nv+1] = j;
          xyz[3*nv+2] = k;
          
            // gid
          gids[nv++] = INDEXG(i, j, k);
        }
    

    rval = pc[p]->get_moab()->create_vertices(&xyz[0], nverts, verts[p]);
    CHECK_ERR(rval);

    rval = pc[p]->get_moab()->tag_set_data(gid_tag, verts[p], &gids[0]);
    if (MB_SUCCESS != rval) return rval;

      // make elements
    nv = 0;
    EntityHandle connect[8], dum_hex;
    for (k = ijkmin[p][2]; k < ijkmax[p][2]; k++) 
      for (j = ijkmin[p][1]; j < ijkmax[p][1]; j++) 
        for (i = ijkmin[p][0]; i < ijkmax[p][0]; i++) {
            // gid
          connect[0] = verts[p][INDEXL(i, j, k)];
          connect[1] = verts[p][INDEXL(i+1, j, k)];
          connect[2] = verts[p][INDEXL(i+1, j+1, k)];
          connect[3] = verts[p][INDEXL(i, j+1, k)];
          connect[4] = verts[p][INDEXL(i, j, k+1)];
          connect[5] = verts[p][INDEXL(i+1, j, k+1)];
          connect[6] = verts[p][INDEXL(i+1, j+1, k+1)];
          connect[7] = verts[p][INDEXL(i, j+1, k+1)];
          rval = pc[p]->get_moab()->create_element(MBHEX, connect, 8, dum_hex);
          hexes[p].insert(dum_hex);
          gids[nv++] = INDEXG(i, j, k);
        }
    rval = pc[p]->get_moab()->tag_set_data(gid_tag, hexes[p], &gids[0]);
    if (MB_SUCCESS != rval) return rval;

    std::ostringstream fname;
    fname << "tmp" << p << ".h5m";
    rval = pc[p]->get_moab()->write_file(fname.str().c_str());
    if (MB_SUCCESS != rval) return rval;
  }
  rval = ParallelComm::resolve_shared_ents(pc, 4, 3);
  CHECK_ERR(rval);
  return rval;
}

void test_pack_vertices()
{
  Core moab;
  ErrorCode rval;
  Range verts;
  
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
  Core moab;
  ErrorCode rval;
  Range elems;
  
    // define some vertices
  const size_t num_verts = 12;
  EntityHandle verts[num_verts];
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
  EntityHandle hexes[num_hex];
  rval = moab.create_element( MBHEX, verts, 8, hexes[0] );
  CHECK_ERR(rval);
  elems.insert( hexes[0] );
  rval = moab.create_element( MBHEX, verts+4, 8, hexes[1] );
  CHECK_ERR(rval);
  elems.insert( hexes[1] );
  
    // define a single quad on the adjacent sides of the hexes
  const size_t num_quad = 1;
  EntityHandle quad;
  rval = moab.create_element( MBQUAD, verts+4, 4, quad );
  CHECK_ERR(rval);
  elems.insert( quad );
  
    // define a decomposition of the first hex into 5 tets
  const size_t num_tet = 5;
  EntityHandle tets[num_tet];
  EntityHandle tet_conn[num_tet][4] = 
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
  EntityHandle tris[num_tri];
  EntityHandle tri_conn[num_tri][3] = { { verts[3], verts[1], verts[4] },
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
  EntityHandle wedges[num_wedge];
  EntityHandle wedge_conn[num_wedge][6] = {
    { verts[4], verts[5], verts[7], verts[8], verts[9], verts[11] },
    { verts[5], verts[6], verts[7], verts[9], verts[10],verts[11] } };
  rval = moab.create_element( MBPRISM, wedge_conn[0], 6, wedges[0] );
  CHECK_ERR(rval);
  elems.insert( wedges[0] );
  rval = moab.create_element( MBPRISM, wedge_conn[1], 6, wedges[1] );
  CHECK_ERR(rval);
  elems.insert( wedges[1] );
  
    // define a pyramid
  EntityHandle pyr;
  rval = moab.create_element( MBPYRAMID, verts, 5, pyr );
  CHECK_ERR(rval);
  elems.insert( pyr );
  
    // pack and unpack mesh
  pack_unpack_noremoteh( moab, elems );
  
    // check_counts
  check_sizes( moab, num_verts, 0, num_tri, num_quad, 0, 
                     num_tet, 1, num_wedge, 0, num_hex, 0 );

    // get connectivity for two hexes and a quad
  Range range;
  const EntityHandle *conn1, *conn2, *conn3;
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
  const EntityHandle combined[12] = { conn1[0], conn1[1], conn1[2], conn1[3],
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
  Core moab;
  ErrorCode rval;
  Range elems;
  
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
  EntityHandle verts[num_vert];
  for (size_t i = 0; i < num_vert; ++i) {
    rval = moab.create_vertex( coords + 3*i, verts[i] );
    CHECK_ERR(rval);
  }
  
    // define two tets
  const size_t num_tet = 2;
  EntityHandle tet_conn[2][10] = {
   { verts[ 0], verts[ 4], verts[ 9], verts[2],
     verts[ 8], verts[12], verts[10],
     verts[ 1], verts[ 3], verts[11] }, 
   { verts[ 0], verts[ 9], verts[ 4], verts[6],
     verts[10], verts[12], verts[ 8],
     verts[ 7], verts[13], verts[ 5] } };
     
  EntityHandle tets[num_tet];
  rval = moab.create_element( MBTET, tet_conn[0], 10, tets[0] );
  CHECK_ERR(rval);
  elems.insert( tets[0] );
  rval = moab.create_element( MBTET, tet_conn[1], 10, tets[1] );
  CHECK_ERR(rval);
  elems.insert( tets[1] );
   
    // define interior tri face
  const size_t num_tri = 1;
  EntityHandle tri_conn[6] = 
    { verts[0], verts[4], verts[9],
      verts[8], verts[12],verts[10] };
  EntityHandle tri;
  rval = moab.create_element( MBTRI, tri_conn, 6, tri );
  CHECK_ERR(rval);
  elems.insert( tri );
  
    // pack and unpack mesh
  pack_unpack_noremoteh( moab, elems );
  
    // check_counts
  check_sizes( moab, num_vert, 0, num_tri, 0, 0, num_tet, 0, 0, 0, 0, 0 );

    // get connectivity for two tets and a tri
  Range range;
  const EntityHandle *conn1, *conn2, *conn3;
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
  const EntityHandle combined[num_vert] = {
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
  Core moab;
  ErrorCode rval;
  Range elems;
  
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
  EntityHandle verts[num_vert];
  for (size_t i = 0; i < num_vert; ++i) {
    rval = moab.create_vertex( coords + 3*i, verts[i] );
    CHECK_ERR(rval);
  }
  
    // define octagonal base
  const size_t num_polygon = 1;
  EntityHandle octagon;
  rval = moab.create_element( MBPOLYGON, verts, 8, octagon );
  CHECK_ERR(rval);
  
    // define triangular sides
  const size_t num_tri = num_vert-1;
  EntityHandle tri[num_tri];
  for (size_t i = 0; i < num_tri; ++i) {
    const EntityHandle conn[3] = { verts[i], verts[(i+1)%num_tri], verts[num_tri] };
    rval = moab.create_element( MBTRI, conn, 3, tri[i] );
    CHECK_ERR(rval);
  }
  
    // define the octagon-based pyramid
  const size_t num_polyhedron = 1;
  EntityHandle polyhedron;
  EntityHandle all_faces[num_vert];
  all_faces[0] = octagon;
  std::copy( tri, tri+num_tri, all_faces+1 );
  rval = moab.create_element( MBPOLYHEDRON, all_faces, num_vert, polyhedron );
  CHECK_ERR(rval);
  
    // pack and unpack the mesh
  elems.clear();
  elems.insert( polyhedron );
  elems.insert( octagon );
  std::copy( tri, tri+num_tri, range_inserter(elems) );
  pack_unpack_noremoteh( moab, elems );
  
    // check counts
  check_sizes( moab, num_vert, 0, num_tri, 0, num_polygon, 
                     0, 0, 0, 0, 0, num_polyhedron );
  
    // get entities
  Range range;
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
  const EntityHandle* oct_conn;
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
  std::vector<EntityHandle> volconn;
  rval = moab.get_connectivity( &polyhedron, 1, volconn );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_tri + 1, volconn.size() );
  CHECK_EQUAL( volconn[0], octagon );
  for (Range::iterator i = range.begin(); i != range.end(); ++i) {
    CHECK( std::find(volconn.begin(), volconn.end(), *i) != volconn.end() );
  }
}

void test_pack_sets_simple()
{
  Core moab;
  ErrorCode rval;
  create_simple_grid( moab, 3 );  

    // delete any existing sets
  Range sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
  
    // get all entities
  Range entities;
  rval = moab.get_entities_by_handle( 0, entities );
  CHECK_ERR(rval);
    // expect 8 elements and 27 vertices
  CHECK_EQUAL( (EntityHandle)35, entities.size() );
  CHECK_EQUAL( 27u, entities.num_of_type( MBVERTEX ) );
  CHECK_EQUAL( 8u, entities.num_of_type( MBHEX ) );
  
  
    // create five sets:
    // 1) one with all the elements and vertices,
    // 2) one with half of the elements, 
    // 3) one with the other half of the elements, 
    // 4) one with a single vertex, 
    // 5) an empty set, 
  EntityHandle all_set, half1_set, half2_set, vertex_set, empty_set;
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

  Range elems, verts, half;
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
  EntityHandle vert = verts.front();
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
  
  CHECK_EQUAL( (EntityHandle)27, verts.size() );
  CHECK_EQUAL( (EntityHandle)8, elems.size() );
  CHECK_EQUAL( (EntityHandle)5, sets.size() );
  
    // guess which is which
  empty_set = all_set = vertex_set = half1_set = half2_set = 0;
  for (Range::iterator i = sets.begin(); i != sets.end(); ++i) {
    int num_vtx, num_elem;
    rval = moab.get_number_entities_by_type( *i, MBVERTEX, num_vtx );
    CHECK_ERR(rval);
    rval = moab.get_number_entities_by_type( *i, MBHEX, num_elem );
    CHECK_ERR(rval);
    if (num_vtx == 0) {
      if (num_elem == 0) {
        CHECK_EQUAL( (EntityHandle)0, empty_set );
        empty_set = *i;
      }
      else if (!half1_set) {
        half1_set = *i;
      }
      else {
        CHECK_EQUAL( (EntityHandle)0, half2_set );
        half2_set = *i;
      }
    }
    else if (num_vtx == 1) {
      CHECK_EQUAL( 0, num_elem );
      CHECK_EQUAL( (EntityHandle)0, vertex_set );
      vertex_set = *i;
    }
    else {
      CHECK_EQUAL( 8, num_elem );
      CHECK_EQUAL( 27, num_vtx );
      CHECK_EQUAL( (EntityHandle)0, all_set );
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
  Core moab;
  ErrorCode rval;
  create_simple_grid( moab, 3 );  

    // delete any existing sets
  Range sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
  
    // get all vertices
  Range vertices;
  rval = moab.get_entities_by_type( 0, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK_EQUAL( (EntityHandle)27, vertices.size() );
    // create meshset containing vertices
  EntityHandle set;
  rval = moab.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  rval = moab.add_entities( set, vertices );
  CHECK_ERR(rval);
  
    // pack and unpack range containing only set handle.
    // Will fail unless we also pass in set contents explicitly
  Range entities( vertices );
  entities.insert( set );
  pack_unpack_noremoteh( moab, entities );
  
    // expect single set in mesh
  entities.clear();
  rval = moab.get_entities_by_type( 0, MBENTITYSET, entities );
  CHECK_ERR(rval);
  CHECK_EQUAL( (EntityHandle)1, entities.size() );
  set = entities.front();
  
    // expect 27 vertices in mesh
  vertices.clear();
  rval = moab.get_entities_by_type( 0, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK_EQUAL( (EntityHandle)27, vertices.size() );
  
    // expect set to contain all 27 vertices
  vertices.clear();
  rval = moab.get_entities_by_type( set, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK_EQUAL( (EntityHandle)27, vertices.size() );
}

void test_pack_sets_of_sets()
{
  Core moab;
  ErrorCode rval;
 
    // delete any existing sets
  Range sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
 
    // create three sets such that set2 contains set1, and set3 contains
    // both set1 and set2
  EntityHandle set1, set2, set3;
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
  CHECK_EQUAL( (EntityHandle)3, sets.size() );
  
    // figure out which is which
  set1 = set2 = set3 = 0;
  for (Range::iterator i = sets.begin(); i != sets.end(); ++i) {
    int count;
    rval = moab.get_number_entities_by_type( *i, MBENTITYSET, count );
    CHECK_ERR(rval);
    CHECK( count >= 0 && count <= 2 );
    switch (count) {
      case 0:
        CHECK_EQUAL( (EntityHandle)0, set1 );
        set1 = *i;
        break;
      case 1:
        CHECK_EQUAL( (EntityHandle)0, set2 );
        set2 = *i;
        break;
      case 2:
        CHECK_EQUAL( (EntityHandle)0, set3 );
        set3 = *i;
        break;
    }
  }
  
    // check that set2 contains set1
  sets.clear();
  rval = moab.get_entities_by_type( set2, MBENTITYSET, sets );
  CHECK_EQUAL( (EntityHandle)1, sets.size() );
  CHECK_EQUAL( set1, sets.front() );
  
    // check that set3 contains set1 and set2
  sets.clear();
  rval = moab.get_entities_by_type( set3, MBENTITYSET, sets );
  CHECK_EQUAL( (EntityHandle)2, sets.size() );
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
  Core moab;
  ErrorCode rval;
 
    // delete any existing sets
  Range sets;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  if (!sets.empty()) {
    rval = moab.delete_entities( sets );
    CHECK_ERR(rval);
  }
 
    // create three sets such that set3 has a child link to
    // set1, set2 has a parent link to set1, and such that set3 and
    // set2 are parent and child, respectively.
  EntityHandle set1, set2, set3;
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
  CHECK_EQUAL( (EntityHandle)3, sets.size() );
  
    // look for a set with two child links (set3)
  set1 = set2 = set3 = 0;
  for (Range::iterator i = sets.begin(); i != sets.end(); ++i) {
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
  std::vector<EntityHandle> parents, children;
  rval = moab.get_child_meshsets( set3, children );
  CHECK_ERR(rval);
  CHECK_EQUAL( (std::vector<EntityHandle>::size_type)2, children.size() );
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
  CHECK_EQUAL( (std::vector<EntityHandle>::size_type)2, parents.size() );
  CHECK_EQUAL( set1, parents[0] );
  CHECK_EQUAL( set3, parents[1] );
}

void test_pack_tag_data_sparse()
{
  Range::iterator i;
  int size;
  DataType type;
  TagType storage;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
  create_simple_grid( mb, 3 );  
  Range elems;
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  CHECK( !elems.empty() );
 
    // Define a sparse tag containing two integers.  For every other element
    // in the mesh, set the tag value to the floor of the x and y
    // coordinates of the first vertex in the elements connectivity list.
  const char sparse_2_int_tag_name[] = "test tag 1";
  Tag sparse_2_int_tag;
  rval = mb.tag_get_handle( sparse_2_int_tag_name,
                        2, MB_TYPE_INTEGER,
                        sparse_2_int_tag,
                        MB_TAG_SPARSE|MB_TAG_CREAT);
  CHECK_ERR(rval);
  bool skip = false;
  for (i = elems.begin(); i != elems.end(); ++i, skip = !skip) {
    if (skip)
      continue;
    
    const EntityHandle* conn =0;
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
  Range ents;
  pack_unpack_noremoteh( moab, ents );
  elems.clear();
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  
  
    // check tag meta for sparse_2_int_tag
  rval = mb.tag_get_handle( sparse_2_int_tag_name, 2, MB_TYPE_INTEGER, sparse_2_int_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_length( sparse_2_int_tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, size );
  rval = mb.tag_get_type( sparse_2_int_tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_SPARSE, storage );
  rval = mb.tag_get_data_type( sparse_2_int_tag, type );
  CHECK_EQUAL( MB_TYPE_INTEGER, type );
  int intdata[2];
  rval = mb.tag_get_default_value( sparse_2_int_tag, intdata );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // check tag data for sparse_2_int_tag
  Range tagged;
  rval = mb.get_entities_by_type_and_tag( 0, MBHEX, 
                                          &sparse_2_int_tag, 0, 1,
                                          tagged );
  CHECK_ERR(rval);
  CHECK_EQUAL( (elems.size()+1) / 2, tagged.size() );
  for (i = tagged.begin(); i != tagged.end(); ++i) {
    rval = mb.tag_get_data( sparse_2_int_tag, &*i, 1, intdata );
    CHECK_ERR(rval);
    
    const EntityHandle* conn =0;
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
  Range::iterator i;
  int size;
  DataType type;
  TagType storage;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
  create_simple_grid( mb, 3 );  
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );

    // Define a dense tag containing a single double-precision floating
    // point value.  For each vertex, store the distance from the origin
    // in this tag.
  const char dense_1_double_tag_name[] = "test tag 2";
  Tag dense_1_double_tag;
  rval = mb.tag_get_handle( dense_1_double_tag_name,
                        1, MB_TYPE_DOUBLE,
                        dense_1_double_tag,
                        MB_TAG_DENSE|MB_TAG_EXCL );
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
  Range ents;
  pack_unpack_noremoteh( moab, ents );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  
    // check tag meta for dense_1_double_tag
  rval = mb.tag_get_handle( dense_1_double_tag_name, 1, MB_TYPE_DOUBLE, dense_1_double_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_length( dense_1_double_tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, size );
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
  Range::iterator i;
  int size;
  DataType type;
  TagType storage;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
  create_simple_grid( mb, 3 );  
  Range verts, elems, sets;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  CHECK( !elems.empty() );

    // Define a dense, opaque tag with a default value of "DEFLT".
    // Set the tag on one element, one vertex,and one set to "TAGGD".
  const char dense_5_opaque_tag_name[] = "This is intentionally a very long tag name in an attempt to test for an arbitrary limitations on tag name length.";
  Tag dense_5_opaque_tag;
  rval = mb.tag_get_handle( dense_5_opaque_tag_name,
                            5, MB_TYPE_OPAQUE,
                            dense_5_opaque_tag,
                            MB_TAG_DENSE|MB_TAG_EXCL,
                            "DEFLT" );
  CHECK_ERR(rval);
  EntityHandle set;
  rval = mb.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  const EntityHandle handles[3] = { verts.front(), elems.front(), set };
  const char data[] = "TAGGDTAGGDTAGGD";
  rval = mb.tag_set_data( dense_5_opaque_tag, handles, 3, data );
  CHECK_ERR(rval);
  
    // pack and unpack
  Range ents;
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
  rval = mb.tag_get_handle( dense_5_opaque_tag_name, 5, MB_TYPE_OPAQUE, dense_5_opaque_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_length( dense_5_opaque_tag, size );
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
  Range::iterator i;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
 
    // Create a bit tag
  const char tag_name[] = "test bit";
  Tag tag;
  rval = mb.tag_get_handle( tag_name, 3, MB_TYPE_BIT, tag, MB_TAG_EXCL );
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
  Range ents;
  pack_unpack_noremoteh( moab, ents );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);

    // check tag meta 
  rval = mb.tag_get_handle( tag_name, 3, MB_TYPE_BIT, tag );
  CHECK_ERR(rval);
  
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, size );
  
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_BIT, storage );
  
  DataType type;
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
  Range::iterator i;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  
    // create a variable-length tag 
  const char* tag_name = "var_int_tag";
  const int defval_size = 5;
  const int default_val[defval_size] = { 0xBEEF, 0xFEED, 0xDEAD, 0xBAD, 0xBEAD };
  Tag tag;
  rval = mb.tag_get_handle( tag_name, defval_size, MB_TYPE_INTEGER, tag, 
                             MB_TAG_DENSE|MB_TAG_VARLEN|MB_TAG_EXCL,
                            default_val );
  CHECK_ERR(rval);
  
    // for each vertex, store in the tag an integer between 1 and 3, 
    // followed by the floor of the cooresponding number of vertex
    // coordinates, beginning with x.
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    const int num_coord = 1 + *i % 3;
    const int data_size = num_coord + 1;
    const int data[4] = { num_coord, (int)coords[0], (int)coords[1], (int)coords[2] };
    const void* data_ptrs[1] = { data };
    rval = mb.tag_set_by_ptr( tag, &*i, 1, data_ptrs, &data_size );
    CHECK_ERR(rval);
  }
  
    // pack and unpack
  pack_unpack_noremoteh( moab );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);

    // check tag meta 
  rval = mb.tag_get_handle( tag_name, 0, MB_TYPE_INTEGER, tag );
  CHECK_ERR(rval);
  
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, rval );
  
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  
  DataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_EQUAL( MB_TYPE_INTEGER, type );

  const void* defval_ptr;
  rval = mb.tag_get_default_value( tag, defval_ptr, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( defval_size, size );
  const int* defval_arr = reinterpret_cast<const int*>(defval_ptr);
  for (int j = 0; j < size; ++j)
    CHECK_EQUAL( default_val[j], defval_arr[j] );
  
    // check tag values
  for (i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    int size;
    const void* valptr;
    rval = mb.tag_get_by_ptr( tag, &*i, 1, &valptr, &size );
    CHECK_ERR(rval);
    CHECK( size > 1 );
    CHECK( size <= 4 );
    
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
  Range::iterator i;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  Range verts, elems;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  rval = mb.get_entities_by_type( 0, MBHEX, elems );
  CHECK_ERR(rval);
  CHECK( !elems.empty() );
 
    // create a tag
  const char* tag_name = "entity tag";
  EntityHandle default_val[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
  Tag tag;
  rval = mb.tag_get_handle( tag_name, 8, MB_TYPE_HANDLE, tag, MB_TAG_SPARSE|MB_TAG_EXCL, &default_val );
  CHECK_ERR(rval);
  
    // Store on each vertex the handles of the adjacent hexes, padded
    // with NULL handles.
  EntityHandle tagdata[8*8];
  for (i = elems.begin(); i != elems.end(); ++i) {
    const EntityHandle* conn;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 8, len );
    
    rval = mb.tag_get_data( tag, conn, len, tagdata );
    CHECK_ERR(rval);
    
    for (int j = 0; j < 8; ++j) {
      EntityHandle* vdata = tagdata + 8*j;
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
  rval = mb.tag_get_handle( tag_name, 8, MB_TYPE_HANDLE, tag );
  CHECK_ERR(rval);
  
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, size );
  
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_SPARSE, storage );
  
  DataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_EQUAL( MB_TYPE_HANDLE, type );

  rval = mb.tag_get_default_value( tag, tagdata );
  CHECK_ERR(rval);
  for (int j = 0; j < 8; ++j) {
    CHECK_EQUAL( (EntityHandle)0, tagdata[j] );
  }
  
    // check tag values
  for (i = elems.begin(); i != elems.end(); ++i) {
    const EntityHandle* conn;
    int len;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 8, len );
    
    rval = mb.tag_get_data( tag, conn, len, tagdata );
    CHECK_ERR(rval);
    
    for (int j = 0; j < 8; ++j) {
      EntityHandle* vdata = tagdata + 8*j;
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
      CHECK_EQUAL( (EntityHandle)0, tagdata[j] );
    }
  }
}
  
ErrorCode get_entities(Interface *mb,
                         std::vector<EntityHandle> &ent_verts, 
                         int verts_per_entity, int dim, 
                         Range &ents) 
{
  assert(!(ent_verts.size()%verts_per_entity));
  unsigned int num_ents = ent_verts.size() / verts_per_entity;
  Range dum_ents;
  ErrorCode result;
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
  
void test_pack_shared_entities_2d()
{
  Core moab[4];
  ParallelComm *pc[4];
  for (unsigned int i = 0; i < 4; i++) {
    pc[i] = new ParallelComm(&moab[i]);
    pc[i]->set_rank(i);
  }

  Range verts[4], quads[4];
  ErrorCode rval = create_shared_grid_2d(pc, verts, quads);

    //moab[0].list_entities(0,1);
  
    // exchange interface cells
  rval = ParallelComm::exchange_ghost_cells(pc, 4, -1, -1, 0, 0, true);
  CHECK_ERR(rval);
  
    // now 1 layer of hex ghosts
  rval = ParallelComm::exchange_ghost_cells(pc, 4, 2, 0, 1, 0, true);
  CHECK_ERR(rval);

    // now 1 layer of hex ghosts w/ edges
  rval = ParallelComm::exchange_ghost_cells(pc, 4, 2, 0, 1, 1, true);
  CHECK_ERR(rval);

  for (unsigned int i = 0; i < 4; i++)
    delete pc[i];
}

void test_pack_shared_entities_3d()
{
  Core moab[4];
  ParallelComm *pc[4];
  for (unsigned int i = 0; i < 4; i++) {
    pc[i] = new ParallelComm(&moab[i]);
    pc[i]->set_rank(i);
    for (unsigned int j = 0; j < 4; j++) {
      if (j == i) continue;
      else pc[i]->get_buffers(j);
    }
  }

  Range verts[4], hexes[4];
  ErrorCode rval = create_shared_grid_3d(pc, verts, hexes);

    // exchange interface cells
  rval = ParallelComm::exchange_ghost_cells(pc, 4, -1, -1, 0, 0, true);
  CHECK_ERR(rval);
  
    // now 1 layer of hex ghosts
  rval = ParallelComm::exchange_ghost_cells(pc, 4, 3, 0, 1, 0, true);
  CHECK_ERR(rval);

    // now 1 layer of hex ghosts w/ faces, edges
  rval = ParallelComm::exchange_ghost_cells(pc, 4, 3, 0, 1, 3, true);
  CHECK_ERR(rval);

  for (unsigned int i = 0; i < 4; i++)
    delete pc[i];
}

void test_filter_pstatus()
{
  Range::iterator i;
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create some mesh
  create_simple_grid( mb, 3 );  
  std::vector<EntityHandle> verts;
  Range dum_vertsr, vertsr;
  rval = mb.get_entities_by_type( 0, MBVERTEX, dum_vertsr );
  CHECK_ERR(rval);
  vertsr.insert(dum_vertsr[0], dum_vertsr[8]);
  for (unsigned int i = 0; i < 9; i++) verts.push_back(vertsr[i]);

  CHECK( !verts.empty() );
 
  ParallelComm *pcomm = new ParallelComm( &moab );

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
  

  Range tmp_range = vertsr;

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
  
  delete pcomm;
}
