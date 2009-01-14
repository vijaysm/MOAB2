#include "TestUtil.hpp"
#include "MBCore.hpp"
#include "Tqdcfr.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"
#include "MBCN.hpp"
#include <math.h>
#include <algorithm>

/**\brief Input test file.
 * Cubit 10.2 file.
 * File contains:
 * Two merged 10x10x10 bricks sharing a single surface (surface 6).
 * - Each brick is meshed with 8 5x5x5 hexes.
 * - Each surface is meshed with 4 5x5 quads.
 * - Each curve is meshed with 2 5-unit edges.
 * A single block containing both bricks.
 * Two side sets
 * - sideset 1: surfaces 1 and 7
 * - sideset 2: surfaces 5 and 11
 * Two node sets:
 * - nodeset 1: surfaces 2 and 8
 * - nodeset 2: surfaces 3 and 9
 *
 * Surfaces:
 *                    2          8
 *                   /          /
 *          o----------o----------o
 *         /.      /  /.      /  /|
 *        / .  (5)/  / . (11)/  / |
 *       /  .    L  /  .    L  /  |
 *      o----------o----------o   |
 *      |   .      |   .      |(12)
 *   4--|-> o . . .|. .o. . . | . o
 *      |  . (1)   |  . (9)   |  /
 *      | .   ^    | .   ^    | /
 *      |.    |    |.    |    |/
 *      o----------o----------o
 *            |          |
 *            3          9
 *
 * Curves:
 *
 *     o----8-----o----20----o
 *    /.         /.         /|
 *  12 .       11 .       24 |
 *  /  7       /  5       /  17
 * o----2-----o----14----o   |
 * |   .      |   .      |   |
 * |   o . .6.|. .o. . 18| . o
 * 3  .       1  .      13  /
 * | 9        | 10       | 22
 * |.         |.         |/
 * o----4-----o----16----o
 *
 * Vertices:
 *
 *     8----------5----------13
 *    /.         /.         /|
 *   / .        / .        / |
 *  /  .       /  .       /  |
 * 3----------2----------10  |
 * |   .      |   .      |   |
 * |   7 . . .|. .6. . . | . 14
 * |  .       |  .       |  /
 * | .        | .        | /
 * |.         |.         |/
 * 4----------1----------9
*/
#ifdef SRCDIR
static const char input_file_1[] = STRINGIFY(SRCDIR) "/test.cub";
static const char quad9_file[] = STRINGIFY(SRCDIR) "/quad9.cub";
static const char hex27_file[] = STRINGIFY(SRCDIR) "/hex27.cub";
#else
static const char input_file_1[] = "test.cub";
static const char quad9_file[] = "quad9.cub";
static const char hex27_file[] = "hex27.cub";
#endif

void read_file( MBInterface& moab, 
                const char* input_file,
                MBEntityHandle* file_set = 0 );

void test_vertices();

void test_edges();

void test_quads();

void test_hexes();

void test_geometric_topology();

void test_geometric_sets();

void test_blocks();

void test_side_sets();

void test_node_sets();

void test_file_set();

void test_quad9();

void test_hex27();

// Check that adjacent lower-order entities have
// higher-order nodes consitent with input entities.
void check_adj_ho_nodes( MBInterface& moab,
                         const MBEntityHandle* entities,
                         int num_entities );

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_vertices);
  result += RUN_TEST(test_edges);
  result += RUN_TEST(test_quads);
  result += RUN_TEST(test_hexes);
  result += RUN_TEST(test_geometric_topology);
  result += RUN_TEST(test_geometric_sets);
  result += RUN_TEST(test_blocks);
  result += RUN_TEST(test_side_sets);
  result += RUN_TEST(test_node_sets);
  result += RUN_TEST(test_file_set);
  result += RUN_TEST(test_hex27);
  result += RUN_TEST(test_quad9);
  
  return result;
}

void read_file( MBInterface& moab, 
                const char* input_file,
                MBEntityHandle* file_set )
{
  MBErrorCode rval;
  MBEntityHandle set;
  Tqdcfr reader( &moab );
  FileOptions opts("");
  rval = reader.load_file( input_file, set, opts, 0, 0 );
  if (file_set)
    *file_set = set;
  CHECK_ERR(rval);
}

void test_vertices()
{
  /* Node coordinates, in order of node ID, beginning with 1. */
  const double node_coords[] = { 5, -5,  5, // 1
                                 5,  5,  5, 
                                 5,  0,  5,
                                -5,  5,  5,
                                 0,  5,  5, // 5
                                -5, -5,  5,
                                -5,  0,  5,
                                 0, -5,  5,
                                 0,  0,  5,
                                 5,  5, -5, // 10
                                 5, -5, -5,
                                 5,  0, -5,
                                -5, -5, -5,
                                 0, -5, -5,
                                -5,  5, -5, // 15
                                -5,  0, -5,
                                 0,  5, -5,
                                 0,  0, -5,
                                -5, -5,  0,
                                 5, -5,  0, // 20
                                 0, -5,  0,
                                -5,  5,  0,
                                -5,  0,  0,
                                 5,  5,  0,
                                 0,  5,  0, // 25
                                 5,  0,  0,
                                 0,  0,  0,
                                15, -5,  5,
                                15,  5,  5,
                                15,  0,  5, // 30
                                10,  5,  5, 
                                10, -5,  5,
                                10,  0,  5,
                                15,  5, -5,
                                15, -5, -5, // 35
                                15,  0, -5,
                                10, -5, -5,
                                10,  5, -5,
                                10,  0, -5,
                                15, -5,  0, // 40
                                10, -5,  0,
                                15,  5,  0,
                                10,  5,  0,
                                15,  0,  0,
                                10,  0,  0  // 45
                                };

  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, input_file_1 );
  
    // get vertex handles and check correct number of vertices
  const size_t num_nodes = sizeof(node_coords)/(3*sizeof(double));
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts ); CHECK_ERR(rval);
  CHECK_EQUAL( num_nodes, (size_t)verts.size() );
  
    // check global ids (should be 1 to 45 for vertices.)
  MBTag gid_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", gid_tag ); CHECK_ERR(rval);
  std::vector<int> ids(num_nodes);
  rval = mb.tag_get_data( gid_tag, verts, &ids[0] ); CHECK_ERR(rval);
  std::vector<int> sorted(ids);
  std::sort( sorted.begin(), sorted.end() );
  for (size_t i = 0; i < num_nodes; ++i)
    CHECK_EQUAL( (int)(i+1), sorted[i] );
    
    // check coordinates of each vertex
  std::vector<double> coords(3*num_nodes);
  rval = mb.get_coords( verts, &coords[0] ); CHECK_ERR(rval);
  for (size_t i = 0; i < num_nodes; ++i) {
    const double* exp = node_coords + 3*(ids[i]-1);
    const double* act = &coords[3*i];
    CHECK_REAL_EQUAL( exp[0], act[0], 1e-8 );
    CHECK_REAL_EQUAL( exp[1], act[1], 1e-8 );
    CHECK_REAL_EQUAL( exp[2], act[2], 1e-8 );
  }
}


void test_element( const char* filename,
                   MBEntityType type, 
                   int num_elem,
                   int node_per_elem,
                   const int* conn_list )
{
  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, filename );
  
  MBRange elems;
  rval = mb.get_entities_by_type( 0, type, elems ); CHECK_ERR(rval);
  CHECK_EQUAL( num_elem, (int)elems.size() );
  
    // get global ids
  MBTag gid_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", gid_tag ); CHECK_ERR(rval);
  std::vector<int> ids(num_elem);
  rval = mb.tag_get_data( gid_tag, elems, &ids[0] ); CHECK_ERR(rval);
  
    // check that global ids are consecutive, beginning with 1
  std::vector<int> sorted(ids);
  std::sort( sorted.begin(), sorted.end() );
  for (int i = 0; i < num_elem; ++i)
    CHECK_EQUAL( i+1, sorted[i] );
  
    // check connectivity of each element
  std::vector<int> conn_ids(node_per_elem);
  std::vector<MBEntityHandle> conn_h;
  MBRange::iterator j = elems.begin();
  for (int i = 0; i < num_elem; ++i, ++j) {
    conn_h.clear();
    rval = mb.get_connectivity( &*j, 1, conn_h ); CHECK_ERR(rval);
    CHECK_EQUAL( node_per_elem, (int)conn_h.size() );
    rval = mb.tag_get_data( gid_tag, &conn_h[0], node_per_elem, &conn_ids[0] );
    CHECK_ERR(rval);
    const int* exp = conn_list + node_per_elem * (ids[i]-1);
    for (int k = 0; k < node_per_elem; ++k)
      CHECK_EQUAL( exp[k], conn_ids[k] );
  }
}
    

void test_edges()
{
  const int edge_conn[] = { 1, 3, // 1
                            3, 2,
                            2, 5, 
                            5, 4,
                            4, 7, 
                            7, 6,
                            6, 8,
                            8, 1,
                            3, 9, 
                            9, 8, // 10
                            5, 9,
                            9, 7, 
                            10, 12,
                            12, 11,
                            11, 14,
                            14, 13,
                            13, 16,
                            16, 15,
                            15, 17,
                            17, 10, // 20
                            12, 18,
                            18, 17,
                            14, 18,
                            18, 16,
                            6, 19,
                            19, 13,
                            1, 20,
                            20, 11,
                            19, 21,
                            21, 8, // 30
                            14, 21,
                            21, 20,
                            4, 22,
                            22, 15,
                            22, 23,
                            23, 7,
                            16, 23,
                            23, 19,
                            2, 24,
                            24, 10, // 40
                            24, 25,
                            25, 5,
                            17, 25,
                            25, 22,
                            20, 26,
                            26, 3,
                            12, 26,
                            26, 24,
                            28, 30,
                            30, 29, // 50
                            29, 31,
                            31, 2,
                            1, 32,
                            32, 28,
                            30, 33,
                            33, 32,
                            31, 33,
                            33, 3,
                            34, 36,
                            36, 35, // 60
                            35, 37,
                            37, 11,
                            10, 38,
                            38, 34,
                            36, 39,
                            39, 38,
                            37, 39,
                            39, 12,
                            28, 40,
                            40, 35, // 70
                            20, 41,
                            41, 32,
                            37, 41,
                            41, 40,
                            29, 42,
                            42, 34,
                            42, 43,
                            43, 31,
                            38, 43,
                            43, 24, // 80
                            40, 44,
                            44, 30,
                            36, 44,
                            44, 42 };
  test_element( input_file_1, MBEDGE, 84, 2, edge_conn );
}

void test_quads()
{
  const int quad_conn[] = { 1, 3, 9, 8, // 1
                            3, 2, 5, 9,
                            8, 9, 7, 6,
                            9, 5, 4, 7,
                           10,12,18,17,
                           12,11,14,18,
                           17,18,16,15,
                           18,14,13,16,
                            6,19,21, 8,
                           19,13,14,21, // 10
                            8,21,20, 1,
                           21,14,11,20,
                            4,22,23, 7,
                           22,15,16,23,
                            7,23,19, 6,
                           23,16,13,19,
                            2,24,25, 5,
                           24,10,17,25,
                            5,25,22, 4,
                           25,17,15,22, // 20
                            1,20,26, 3,
                           20,11,12,26,
                            3,26,24, 2,
                           26,12,10,24,
                           28,30,33,32,
                           30,29,31,33,
                           32,33, 3, 1,
                           33,31, 2, 3,
                           34,36,39,38,
                           36,35,37,39, // 30
                           38,39,12,10,
                           39,37,11,12,
                            1,20,41,32,
                           20,11,37,41,
                           32,41,40,28,
                           41,37,35,40,
                           29,42,43,31,
                           42,34,38,43,
                           31,43,24, 2,
                           43,38,10,24, // 40
                           28,40,44,30,
                           40,35,36,44,
                           30,44,42,29,
                           44,36,34,42,
                          };
  test_element( input_file_1, MBQUAD, 44, 4, quad_conn );
}

void test_hexes()
{
  const int hex_conn[] = {  6, 19, 23,  7,  8, 21, 27,  9,
                           19, 13, 16, 23, 21, 14, 18, 27,
                            7, 23, 22,  4,  9, 27, 25,  5,
                           23, 16, 15, 22, 27, 18, 17, 25,
                            8, 21, 27,  9,  1, 20, 26,  3,
                           21, 14, 18, 27, 20, 11, 12, 26,
                            9, 27, 25,  5,  3, 26, 24,  2,
                           27, 18, 17, 25, 26, 12, 10, 24,
                            1, 20, 26,  3, 32, 41, 45, 33,
                           20, 11, 12, 26, 41, 37, 39, 45,
                            3, 26, 24,  2, 33, 45, 43, 31,
                           26, 12, 10, 24, 45, 39, 38, 43,
                           32, 41, 45, 33, 28, 40, 44, 30,
                           41, 37, 39, 45, 40, 35, 36, 44,
                           33, 45, 43, 31, 30, 44, 42, 29,
                           45, 39, 38, 43, 44, 36, 34, 42 };
  test_element( input_file_1, MBHEX, 16, 8, hex_conn );
}

template <int L>
std::vector<int> find_parents( const int parent_conn[][L], int num_parent, int id )
{
  std::vector<int> results;
  for (int i = 0; i < num_parent; ++i) {
    for (int j = 0; j < L; ++j) {
      if (parent_conn[i][j] == id)
        results.push_back( i+1 );
    }
  }
  return results;
}

int check_geometric_set( MBInterface& moab,
                         int dim, int id, 
                         const int* children, int num_children, 
                         std::vector<int> parents )
{
  MBErrorCode rval;
  MBTag gid_tag, dim_tag;
  
  rval = moab.tag_get_handle( "GLOBAL_ID", gid_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( "GEOM_DIMENSION", dim_tag ); CHECK_ERR(rval);
  void* tag_vals[] = { &dim, &id };
  MBTag tags[] = { dim_tag, gid_tag };
  MBRange ents;
  rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, tag_vals, 2, ents );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1u, (unsigned)ents.size() );
  
  const MBEntityHandle geom = ents.front();
  std::vector<int> exp_rel, act_rel;
  std::vector<MBEntityHandle> rel;
  
  if (num_children) {
    exp_rel.resize( num_children );
    std::copy( children, children+num_children, exp_rel.begin() );
    std::sort( exp_rel.begin(), exp_rel.end() );
    rel.clear();
    rval = moab.get_child_meshsets( geom, rel ); CHECK_ERR(rval);
    CHECK_EQUAL( num_children, (int)rel.size() );
    act_rel.resize( rel.size() );
    rval = moab.tag_get_data( gid_tag, &rel[0], rel.size(), &act_rel[0] ); 
    CHECK_ERR(rval);
    std::sort( act_rel.begin(), act_rel.end() );
    CHECK( exp_rel == act_rel );
  }
  
  if (!parents.empty()) {
    exp_rel = parents;
    std::sort( exp_rel.begin(), exp_rel.end() );
    rel.clear();
    rval = moab.get_parent_meshsets( geom, rel ); CHECK_ERR(rval);
    CHECK_EQUAL( parents.size(), rel.size() );
    act_rel.resize( rel.size() );
    rval = moab.tag_get_data( gid_tag, &rel[0], rel.size(), &act_rel[0] ); 
    CHECK_ERR(rval);
    std::sort( act_rel.begin(), act_rel.end() );
    CHECK( exp_rel == act_rel );
  }
  
  return 0;
}

void test_geometric_topology()
{
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, input_file_1 );
  // expected geometric vertices, specified by global ID
  const int vertex_ids[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14 };
  // List of global IDs of surfacs in geometric volumes, indexed by ID-1  
  const int volume_surfs[2][6] = { { 1, 2, 3, 4, 5, 6 },
                                   { 7, 8, 9, 6,11,12 } };
  // List of global IDs of curves in geometric surfaces, indexed by ID-1  
  // Curve IDs of zero indicates that corresponding surface doesn't exist.
  const int surf_curves[12][4] = { { 1, 2, 3, 4 },
                                   { 5, 6, 7, 8 },
                                   { 9, 6,10, 4 },
                                   {11, 7, 9, 3 },
                                   {12, 8,11, 2 },
                                   {10, 5,12, 1 },
                                   {13,14, 1,16 },
                                   {17,18, 5,20 },
                                   {10,18,22,16 },
                                   { 0, 0, 0, 0 }, // no surf 10
                                   {24,20,12,14 },
                                   {22,17,24,13 } };
  // List of global IDs of vertices in geometric curves, indexed by ID-1  
  // Vertex IDs of zero indicates that corresponding curve doesn't exist.
  const int curve_verts[24][2] =  { { 1, 2 },
                                    { 2, 3 },
                                    { 3, 4 },
                                    { 4, 1 },
                                    { 5, 6 }, // 5
                                    { 6, 7 },
                                    { 7, 8 },
                                    { 8, 5 },
                                    { 4, 7 },
                                    { 1, 6 }, // 10
                                    { 3, 8 },
                                    { 2, 5 },
                                    { 9,10 }, // 13
                                    {10, 2 },
                                    { 0, 0 }, // no curve 15
                                    { 1, 9 },
                                    {13,14 },
                                    {14, 6 }, 
                                    { 0, 0 }, // no curve 19
                                    { 5,13 }, 
                                    { 0, 0 }, // no curve 21
                                    { 9,14 }, 
                                    { 0, 0 }, // no curve 23
                                    {10,13 } };
  
    // check all vertices
  for (unsigned i = 0; i < (sizeof(vertex_ids)/sizeof(vertex_ids[0])); ++i)
    check_geometric_set( mb, 0, vertex_ids[i], 0, 0, 
                         find_parents<2>(curve_verts,24,vertex_ids[i]) );
  
    // check all curves
  for (int i = 1; i <= 24; ++i)
    if (curve_verts[i-1][0])
      check_geometric_set( mb, 1, i, curve_verts[i-1], 2, 
                           find_parents<4>(surf_curves,12,i) );
   
    // check all surfs
  for (int i = 1; i <= 12; ++i)
    if (surf_curves[i-1][0])
      check_geometric_set( mb, 2, i, surf_curves[i-1], 4, 
                          find_parents<6>(volume_surfs,2,i) );
    
    // check all volumes
  std::vector<int> empty;
  for (int i = 1; i <= 2; ++i)
    check_geometric_set( mb, 3, i, volume_surfs[i-1], 6, empty );
}

void test_geometric_sets()
{
  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, input_file_1 );
  MBTag gid_tag, dim_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", gid_tag ); CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GEOM_DIMENSION", dim_tag ); CHECK_ERR(rval);

    // verify mesh entity counts
  MBRange verts, curves, surfs, vols;
  int dim = 0;
  const void* vals[] = {&dim};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, vals, 1, verts );
  CHECK_ERR(rval);
  dim = 1;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, vals, 1, curves );
  CHECK_ERR(rval);
  dim = 2;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, vals, 1, surfs );
  CHECK_ERR(rval);
  dim = 3;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, vals, 1, vols );
  CHECK_ERR(rval);
  
  CHECK_EQUAL( 12u, (unsigned)verts.size() );
  CHECK_EQUAL( 20u, (unsigned)curves.size() );
  CHECK_EQUAL( 11u, (unsigned)surfs.size() );
  CHECK_EQUAL( 2u,  (unsigned)vols.size() );
  
    // check that each vertex has a single node, and that the 
    // node is also contained in any parent curve
  MBRange ents;
  MBRange::iterator i;
  for (i = verts.begin(); i != verts.end(); ++i) {
    ents.clear();
    rval = mb.get_entities_by_handle( *i, ents ); CHECK_ERR(rval);
    CHECK_EQUAL( 1u, (unsigned)ents.size() );
    CHECK( ents.all_of_type(MBVERTEX) );
  }

    // check that each curve has one node and two edges
  for (i = curves.begin(); i != curves.end(); ++i) {
    ents.clear();
    rval = mb.get_entities_by_handle( *i, ents ); CHECK_ERR(rval);
    CHECK_EQUAL( 1u, (unsigned)ents.num_of_type(MBVERTEX) );
    CHECK_EQUAL( 2u, (unsigned)ents.num_of_type(MBEDGE) );
    CHECK_EQUAL( 3u, (unsigned)ents.size() );
  }

    // check that each surface has 1 node, 4 edges, 4 quads
  for (i = surfs.begin(); i != surfs.end(); ++i) {
    ents.clear();
    rval = mb.get_entities_by_handle( *i, ents ); CHECK_ERR(rval);
    CHECK_EQUAL( 1u, (unsigned)ents.num_of_type(MBVERTEX) );
    CHECK_EQUAL( 4u, (unsigned)ents.num_of_type(MBEDGE) );
    CHECK_EQUAL( 4u, (unsigned)ents.num_of_type(MBQUAD) );
    CHECK_EQUAL( 9u, (unsigned)ents.size() );
  }

    // check that each volume has 1 node and 8 hexes.
  for (i = vols.begin(); i != vols.end(); ++i) {
    ents.clear();
    rval = mb.get_entities_by_handle( *i, ents ); CHECK_ERR(rval);
    CHECK_EQUAL( 1u, (unsigned)ents.num_of_type(MBVERTEX) );
    CHECK_EQUAL( 8u, (unsigned)ents.num_of_type(MBHEX) );
    CHECK_EQUAL( 9u, (unsigned)ents.size() );
  }
  
    // Check that for each geometric entity, any contained vertices
    // are adjacent to some entity in one of its parents.
  MBRange parents, geom, nodes, tmp;
  for (int d = 0; d < 3; ++d) {
    const void* vals[] = {&d};
    rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, vals, 1, geom );
    CHECK_ERR(rval);
    
    for (i = geom.begin(); i != geom.end(); ++i) {
      nodes.clear();
      ents.clear();
      parents.clear();
      rval = mb.get_entities_by_type( *i, MBVERTEX, nodes ); CHECK_ERR(rval);
      rval = mb.get_parent_meshsets( *i, parents ); CHECK_ERR(rval);
      for (MBRange::iterator j = parents.begin(); j != parents.end(); ++j) {
        tmp.clear();
        rval = mb.get_entities_by_dimension( *j, d+1, tmp ); CHECK_ERR(rval);
        ents.merge( tmp );
      }
      tmp.clear();
      rval = mb.get_adjacencies( ents, 0, false, tmp, MBInterface::UNION );
      CHECK_ERR(rval);
      nodes = nodes.subtract( tmp );
      CHECK( nodes.empty() );
    }
  }
}

// expect one block containing entire mesh, with id == 1
void test_blocks()
{
  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, input_file_1 );
  MBTag mat_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, mat_tag ); CHECK_ERR(rval);

  MBRange blocks;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &mat_tag, 0, 1, blocks );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1u, (unsigned)blocks.size() );
  MBEntityHandle block = blocks.front();
  int id;
  rval = mb.tag_get_data( mat_tag, &block, 1, &id );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, id );
  
  MBRange block_hexes, mesh_hexes;
  rval = mb.get_entities_by_dimension( 0, 3, mesh_hexes ); CHECK_ERR(rval);
  rval = mb.get_entities_by_dimension( block, 3, block_hexes, true ); CHECK_ERR(rval);
  CHECK( mesh_hexes == block_hexes );  
}

// Common code for test_side_sets and test_node sets
//\param tag_name  NEUMANN_SET_TAG_NAME or DIRICHLET_SET_TAG_NAME
//\param count     Number of expected sets
//\param ids       Expected IDs of sets
//\param set_surfs One list for each id in "ids" containing the
//                 ids of the geometric surfaces expected to be
//                 contained in the boundary condition set.
void test_bc_sets( const char* tag_name, unsigned count, 
                   const int* ids,
                   const std::vector<int> set_surfs[] )
{
  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, input_file_1 );
  MBTag ss_tag, gid_tag, dim_tag;
  rval = mb.tag_get_handle( tag_name, ss_tag ); CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GLOBAL_ID", gid_tag ); CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GEOM_DIMENSION", dim_tag ); CHECK_ERR(rval);

    // check number of sidesets and IDs
  MBRange sidesets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &ss_tag, 0, 1, sidesets );
  CHECK_ERR(rval);
  CHECK_EQUAL( count, (unsigned)sidesets.size() );
  std::vector<MBEntityHandle> handles(count,0);
  for (MBRange::iterator i = sidesets.begin(); i != sidesets.end(); ++i) {
    int id;
    rval = mb.tag_get_data( ss_tag, &*i, 1, &id ); CHECK_ERR(rval);
    unsigned idx;
    for (idx = 0; idx < count; ++idx) 
      if (ids[idx] == id)
        break;
    CHECK( idx != count );
    CHECK( handles[idx] == 0 );
    handles[idx] = *i;
  }

    
    // get surface faces
  std::vector<MBRange> exp(count);
  MBRange surfs, tmp;
  MBTag tags[] = { dim_tag, gid_tag };
  for (unsigned i = 0; i < count; ++i) {
    exp[i].clear();
    surfs.clear();
    const int two = 2;
    for (unsigned j = 0; j < set_surfs[i].size(); ++j) {
      const void* vals[] = { &two, &set_surfs[i][j] };
      surfs.clear();
      rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, surfs );
      CHECK_ERR(rval);
      CHECK_EQUAL( 1u, (unsigned)surfs.size() );
      tmp.clear();
      rval = mb.get_entities_by_dimension( surfs.front(), 2, tmp, true ); CHECK_ERR(rval);
      exp[i].merge( tmp );
    }
  }
  
    // check each bc set
  MBRange act;
  for (unsigned i = 0; i < count; ++i) {
    act.clear();
    rval = mb.get_entities_by_dimension( handles[i], 2, act, true ); CHECK_ERR(rval);
    CHECK( exp[i] == act );
  }
}


// expect two sidesets containing two geometric surfaces each:
// sideset 1 : surfaces 1 and 7
// sideset 2 : surfaces 5 and 11
void test_side_sets()
{
  int ids[] = { 1, 2 };
  std::vector<int> surfs[2];
  surfs[0].push_back( 1 );
  surfs[0].push_back( 7 );
  surfs[1].push_back( 5 );
  surfs[1].push_back( 11 );
  test_bc_sets( NEUMANN_SET_TAG_NAME, 2, ids, surfs );
}

void test_node_sets()
{
  int ids[] = { 1, 2 };
  std::vector<int> surfs[2];
  surfs[0].push_back( 2 );
  surfs[0].push_back( 8 );
  surfs[1].push_back( 3 );
  surfs[1].push_back( 9 );
  test_bc_sets( DIRICHLET_SET_TAG_NAME, 2, ids, surfs );
}

void test_file_set()
{
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  MBEntityHandle set;
  read_file( mb, input_file_1, &set );
  
  MBErrorCode rval;
  MBRange exp, act;
  rval = mb.get_entities_by_handle( 0, exp ); CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( set, act ); CHECK_ERR(rval);
  exp.erase( set );
  CHECK( exp == act );
}

static int q9_idx_from_coord( double coord )
{
    // 10x10 element grid has 21 vertices in each direction,
    // equally spaced in [-5,5].
  return (int)round( 2*coord + 10 );
}

void test_quad9()
{
  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, quad9_file );
  
    // we're expecting a 10x10 element grid in the xy plane.
  MBEntityHandle verts[21][21], quads[10][10];
  memset( verts, 0, sizeof(verts) );
  memset( quads, 0, sizeof(quads) );
  
    // get all quads
  MBRange range;
  rval = mb.get_entities_by_type( 0, MBQUAD, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)(10*10), range.size() );
  
    // for each quad
  for (MBRange::iterator i = range.begin(); i != range.end(); ++i) {
      // get connectivity
    const MBEntityHandle* conn = 0;
    int conn_len = 0;
    rval = mb.get_connectivity( *i, conn, conn_len, false );
    CHECK_ERR(rval);
    CHECK_EQUAL( 9, conn_len );
    
      // translate connectivity to indices into vertex grid
    double coords[3*9];
    rval = mb.get_coords( conn, 9, coords ); CHECK_ERR(rval);
    int x[9], y[9];
    for (int j = 0; j < 9; ++j) {
      x[j] = q9_idx_from_coord( coords[3*j] );
      y[j] = q9_idx_from_coord( coords[3*j+1] );
      CHECK_REAL_EQUAL( 0.0, coords[3*j+2], 1e-12 );
      CHECK( 0 <= x[j] && x[j] <= 20 );
      CHECK( 0 <= y[j] && y[j] <= 20 );
    }
    
      // calculate element location from mid-node location
    int ex = (x[8]-1)/2;
    int ey = (y[8]-1)/2;
    CHECK_EQUAL( (MBEntityHandle)0, quads[ex][ey] );
    quads[ex][ey] = *i;
    
      // check each vertex
    bool seen[3][3] = { {false, false, false},
                        {false, false, false},
                        {false, false, false} };
    for (int j = 0; j < 9; ++j) {
        // check for duplicate vertices
      if (verts[x[j]][y[j]] == 0)
        verts[x[j]][y[j]] = conn[j];
      else {
        CHECK_EQUAL( verts[x[j]][y[j]], conn[j] );
      }
        // check that vertices have correct coordinates
      int sx = x[j] - x[8] + 1;
      int sy = y[j] - y[8] + 1;
      CHECK( 0 <= sx && sx <= 2 );
      CHECK( 0 <= sy && sy <= 2 );
      seen[sx][sy] = true;
    }
    
    CHECK( seen[0][0] );
    CHECK( seen[0][1] );
    CHECK( seen[0][2] );
    CHECK( seen[1][0] );
    CHECK( seen[1][1] );
    CHECK( seen[1][2] );
    CHECK( seen[2][0] );
    CHECK( seen[2][1] );
    CHECK( seen[2][2] );
  }
  
    // check that we saw all the elements and vertices we expected to
  for (int i = 0; i < 21; ++i)
    for (int j = 0; j < 21; ++j)
      CHECK( verts[i][j] != 0 );
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      CHECK( quads[i][j] != 0 );
        
  check_adj_ho_nodes( mb, &quads[0][0], 10*10 );
}


static int h27_idx_from_coord( double coord )
{
    // 2x2 element grid has 7 vertices in each direction,
    // equally spaced in [-1,1].
  return (int)round( 2*coord + 2 );
}

void test_hex27()
{
  MBErrorCode rval;
  MBCore mb_impl;
  MBInterface& mb = mb_impl;
  read_file( mb, hex27_file );
  
    // we're expecting a 2x2x2 element grid
  MBEntityHandle verts[5][5][5], hexes[2][2][2];
  memset( verts, 0, sizeof(verts) );
  memset( hexes, 0, sizeof(hexes) );
  
    // get all hexes
  MBRange range;
  rval = mb.get_entities_by_type( 0, MBHEX, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)(2*2*2), range.size() );
  
    // for each quad
  for (MBRange::iterator i = range.begin(); i != range.end(); ++i) {
      // get connectivity
    const MBEntityHandle* conn = 0;
    int conn_len = 0;
    rval = mb.get_connectivity( *i, conn, conn_len, false );
    CHECK_ERR(rval);
    CHECK_EQUAL( 27, conn_len );
    
      // translate connectivity to indices into vertex grid
    double coords[3*27];
    rval = mb.get_coords( conn, 27, coords ); CHECK_ERR(rval);
    int x[27], y[27], z[27];
    for (int j = 0; j < 27; ++j) {
      x[j] = h27_idx_from_coord( coords[3*j] );
      y[j] = h27_idx_from_coord( coords[3*j+1] );
      z[j] = h27_idx_from_coord( coords[3*j+2] );
      CHECK( 0 <= x[j] && x[j] <= 4 );
      CHECK( 0 <= y[j] && y[j] <= 4 );
      CHECK( 0 <= z[j] && z[j] <= 4 );
    }
    
      // calculate element location from mid-node location
    int ex = (x[26]-1)/2;
    int ey = (y[26]-1)/2;
    int ez = (z[26]-1)/2;
    CHECK_EQUAL( (MBEntityHandle)0, hexes[ex][ey][ez] );
    hexes[ex][ey][ez] = *i;
    
      // check each vertex
    bool seen[3][3][3];
    std::fill( &seen[0][0][0], &seen[0][0][0] + 3*3*3, false );
    for (int j = 0; j < 27; ++j) {
        // check for duplicate vertices
      if (verts[x[j]][y[j]][z[j]] == 0)
        verts[x[j]][y[j]][z[j]] = conn[j];
      else {
        CHECK_EQUAL( verts[x[j]][y[j]][z[j]], conn[j] );
      }
        // check that vertices have correct coordinates
      int sx = x[j] - x[26] + 1;
      int sy = y[j] - y[26] + 1;
      int sz = z[j] - z[26] + 1;
      CHECK( 0 <= sx && sx <= 2 );
      CHECK( 0 <= sy && sy <= 2 );
      CHECK( 0 <= sz && sz <= 2 );
      seen[sx][sy][sz] = true;
    }
    
    for (int k = 0; k < 3; ++k) 
      for (int l = 0; l < 3; ++l)
        for (int m = 0; m < 3; ++m)
          CHECK( seen[k][l][m] );
  }
  
    // check that we saw all the elements and vertices we expected to
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k)
        CHECK( verts[i][j][k] != 0 );
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        CHECK( hexes[i][j][k] != 0 );
        
  check_adj_ho_nodes( mb, &hexes[0][0][0], 2*2*2 );
}

static MBEntityHandle find_side( MBInterface& moab, 
                                 MBEntityHandle entity,
                                 int side_dim,
                                 int side_num )
{
  MBErrorCode rval;
  
  std::vector<MBEntityHandle> adj;
  rval = moab.get_adjacencies( &entity, 1, side_dim, false, adj );
  CHECK_ERR(rval);
  
  int sub_ent_indices[4];
  MBCN::SubEntityVertexIndices( TYPE_FROM_HANDLE(entity), side_dim, side_num, 
                                sub_ent_indices );
  MBEntityType subtype = MBCN::SubEntityType( TYPE_FROM_HANDLE(entity),
                                              side_dim, side_num );
  int sub_ent_corners = MBCN::VerticesPerEntity(subtype);
  
  const MBEntityHandle* conn;
  int conn_len;
  rval = moab.get_connectivity( entity, conn, conn_len );
  CHECK_ERR(rval);
  
  for (size_t i = 0; i < adj.size(); ++i) {
    if (TYPE_FROM_HANDLE(adj[i]) != subtype)
      continue;
  
    const MBEntityHandle* sub_conn;
    int sub_len;
    rval = moab.get_connectivity( adj[i], sub_conn, sub_len );
    CHECK_ERR(rval);
    
    int n = std::find( sub_conn, sub_conn+sub_len, conn[sub_ent_indices[0]] ) 
                - sub_conn;
    if (n == sub_len) // no vertex in common
      continue;
  
      // check forward direction
    int j;
    for (j = 1; j < sub_ent_corners; ++j)
      if (conn[sub_ent_indices[j]] != sub_conn[(j+n)%sub_ent_corners])
        break;
    if (j == sub_ent_corners)
      return adj[i];
  
      // check reverse direction
    for (j = 1; j < sub_ent_corners; ++j)
      if (conn[sub_ent_indices[j]] != sub_conn[(n+sub_ent_corners-j)%sub_ent_corners])
        break;
    if (j == sub_ent_corners)
      return adj[i];
  }
  
  // no match
  return 0;
}
      
  
void check_adj_ho_nodes( MBInterface& moab,
                         const MBEntityHandle* entities,
                         int num_entities )
{
  for (int i = 0; i < num_entities; ++i) {
    MBEntityType type = TYPE_FROM_HANDLE(entities[i]);
    const MBEntityHandle* conn;
    int conn_len;
    MBErrorCode rval = moab.get_connectivity( entities[i], conn, conn_len );
    CHECK_ERR(rval);
    
    int ho[4];
    MBCN::HasMidNodes( type, conn_len, ho );
    for (int dim = MBCN::Dimension(type)-1; dim > 0; --dim) {
      if (!ho[dim])
        continue;
        
      for (int j = 0; j < MBCN::NumSubEntities( type, dim ); ++j) {
        MBEntityHandle side = find_side( moab, entities[i], dim, j );
        if (!side)
          continue;
        
        const MBEntityHandle* side_conn;
        int side_len;
        rval = moab.get_connectivity( side, side_conn, side_len );
        CHECK_ERR(rval);
        
        int this_idx = MBCN::HONodeIndex( type, conn_len, dim, j );
        int side_idx = MBCN::HONodeIndex( TYPE_FROM_HANDLE(side), side_len, dim, 0 );
        CHECK_EQUAL( side_conn[side_idx], conn[this_idx] );
      }
    }
  }
}

