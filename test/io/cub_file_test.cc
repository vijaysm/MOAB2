#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include "moab/Range.hpp"
#include "moab/GeomTopoTool.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

/**\brief Input test file: test.cub
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

/* Input test file: ho_test.cub
 * 
 * File is expected to contain at least one block for every
 * supported higher-order element type.  The coordinates of
 * every higher-order node are expected to be the mean of the
 * adjacent corner vertices of the element.
 */
static const char* input_file_1 = std::string( TestDir + "/io/test.cub" ).c_str();
static const char* ho_file = std::string( TestDir + "/io/ho_test.cub" ).c_str();
static const char* cubit12_file = std::string( TestDir + "/io/cubtest12.cub" ).c_str();
static const char* cubit14_file = std::string( TestDir + "/io/cubtest14.cub" ).c_str();

void read_file( Interface& moab, const char* input_file );


// Check that adjacent lower-order entities have
// higher-order nodes consitent with input entity.
void check_adj_ho_nodes( Interface& moab,
                         EntityHandle entity );

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( Interface& moab, 
                       EntityHandle entity,
                       int mid_nodes[4] );

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( EntityType type, int num_nodes );

void test_vertices();

void test_edges();

void test_quads();

void test_hexes();

void test_geometric_topology();

void test_geometric_sets();

void test_blocks();

void test_side_sets();

void test_node_sets();

void test_tri6 () { test_ho_elements(MBTRI, 6); }
void test_tri7 () { test_ho_elements(MBTRI, 7); }

void test_quad5() { test_ho_elements(MBQUAD, 5); }
void test_quad8() { test_ho_elements(MBQUAD, 8); }
void test_quad9() { test_ho_elements(MBQUAD, 9); }

void test_tet8 () { test_ho_elements(MBTET,  8); }
void test_tet10() { test_ho_elements(MBTET, 10); }
void test_tet14() { test_ho_elements(MBTET, 14); }

void test_hex9 () { test_ho_elements(MBHEX,  9); }
void test_hex20() { test_ho_elements(MBHEX, 20); }
void test_hex27() { test_ho_elements(MBHEX, 27); }

void test_multiple_files();                    

void test_cubit12();
void test_cubit14();

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
  result += RUN_TEST(test_tri6 );
  result += RUN_TEST(test_tri7 );
  result += RUN_TEST(test_quad5);
  result += RUN_TEST(test_quad8);
  result += RUN_TEST(test_quad9);
  result += RUN_TEST(test_tet8 );
  result += RUN_TEST(test_tet10);
  result += RUN_TEST(test_tet14);
  result += RUN_TEST(test_hex9 );
  result += RUN_TEST(test_hex20);
  result += RUN_TEST(test_hex27);
  result += RUN_TEST(test_multiple_files);
  result += RUN_TEST(test_cubit12);
  result += RUN_TEST(test_cubit14);
  return result;
}

void read_file( Interface& moab, const char* input_file )
{
  ErrorCode rval = moab.load_file( input_file );
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

  ErrorCode rval;
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, input_file_1 );
  
    // get vertex handles and check correct number of vertices
  const size_t num_nodes = sizeof(node_coords)/(3*sizeof(double));
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts ); CHECK_ERR(rval);
  CHECK_EQUAL( num_nodes, (size_t)verts.size() );
  
    // check global ids (should be 1 to 45 for vertices.)
  Tag gid_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid_tag ); CHECK_ERR(rval);
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
                   EntityType type, 
                   int num_elem,
                   int node_per_elem,
                   const int* conn_list )
{
  ErrorCode rval;
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, filename );
  
  Range elems;
  rval = mb.get_entities_by_type( 0, type, elems ); CHECK_ERR(rval);
  CHECK_EQUAL( num_elem, (int)elems.size() );
  
    // get global ids
  Tag gid_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid_tag ); CHECK_ERR(rval);
  std::vector<int> ids(num_elem);
  rval = mb.tag_get_data( gid_tag, elems, &ids[0] ); CHECK_ERR(rval);
  
    // check that global ids are consecutive, beginning with 1
  std::vector<int> sorted(ids);
  std::sort( sorted.begin(), sorted.end() );
  for (int i = 0; i < num_elem; ++i)
    CHECK_EQUAL( i+1, sorted[i] );
  
    // check connectivity of each element
  std::vector<int> conn_ids(node_per_elem);
  std::vector<EntityHandle> conn_h;
  Range::iterator j = elems.begin();
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

int check_geometric_set( Interface& moab,
                         int dim, int id, 
                         const int* children, int num_children, 
                         std::vector<int> parents )
{
  ErrorCode rval;
  Tag gid_tag, dim_tag;
  
  rval = moab.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( "GEOM_DIMENSION", 1, MB_TYPE_INTEGER, dim_tag ); CHECK_ERR(rval);
  void* tag_vals[] = { &dim, &id };
  Tag tags[] = { dim_tag, gid_tag };
  Range ents;
  rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, tag_vals, 2, ents );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1u, (unsigned)ents.size() );
  
  const EntityHandle geom = ents.front();
  std::vector<int> exp_rel, act_rel;
  std::vector<EntityHandle> rel;
  
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
  Core mb_impl;
  Interface& mb = mb_impl;
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
  ErrorCode rval;
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, input_file_1 );
  Tag gid_tag, dim_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid_tag ); CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GEOM_DIMENSION", 1, MB_TYPE_INTEGER, dim_tag ); CHECK_ERR(rval);

    // verify mesh entity counts
  Range verts, curves, surfs, vols;
  int dim = 0;
  // Cppcheck warning (false positive): variable dim is assigned a value that is never used
  // Cppcheck warning (false positive): variable dim is reassigned a value before the old one has been used
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
  Range ents;
  Range::iterator i;
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
  Range parents, geom, nodes, tmp;
  for (int d = 0; d < 3; ++d) {
    const void* vals1[] = {&d};
    rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, vals1, 1, geom );
    CHECK_ERR(rval);
    
    for (i = geom.begin(); i != geom.end(); ++i) {
      nodes.clear();
      ents.clear();
      parents.clear();
      rval = mb.get_entities_by_type( *i, MBVERTEX, nodes ); CHECK_ERR(rval);
      rval = mb.get_parent_meshsets( *i, parents ); CHECK_ERR(rval);
      for (Range::iterator j = parents.begin(); j != parents.end(); ++j) {
        tmp.clear();
        rval = mb.get_entities_by_dimension( *j, d+1, tmp ); CHECK_ERR(rval);
        ents.merge( tmp );
      }
      tmp.clear();
      rval = mb.get_adjacencies( ents, 0, false, tmp, Interface::UNION );
      CHECK_ERR(rval);
      nodes = subtract( nodes, tmp );
      CHECK( nodes.empty() );
    }
  }
}

// expect one block containing entire mesh, with id == 1
void test_blocks()
{
  ErrorCode rval;
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, input_file_1 );
  Tag mat_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, mat_tag ); CHECK_ERR(rval);

  Range blocks;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &mat_tag, 0, 1, blocks );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1u, (unsigned)blocks.size() );
  EntityHandle block = blocks.front();
  int id;
  rval = mb.tag_get_data( mat_tag, &block, 1, &id );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, id );
  
  Range block_hexes, mesh_hexes;
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
  ErrorCode rval;
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, input_file_1 );
  Tag ss_tag, gid_tag, dim_tag;
  rval = mb.tag_get_handle( tag_name, 1, MB_TYPE_INTEGER, ss_tag ); CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid_tag ); CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GEOM_DIMENSION", 1, MB_TYPE_INTEGER, dim_tag ); CHECK_ERR(rval);

    // check number of sidesets and IDs
  Range sidesets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &ss_tag, 0, 1, sidesets );
  CHECK_ERR(rval);
  CHECK_EQUAL( count, (unsigned)sidesets.size() );
  std::vector<EntityHandle> handles(count,0);
  for (Range::iterator i = sidesets.begin(); i != sidesets.end(); ++i) {
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
  std::vector<Range> exp(count);
  Range surfs, tmp;
  Tag tags[] = { dim_tag, gid_tag };
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
  Range act;
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

static EntityHandle find_side( Interface& moab, 
                                 EntityHandle entity,
                                 int side_dim,
                                 int side_num )
{
  ErrorCode rval;
  
  std::vector<EntityHandle> adj;
  rval = moab.get_adjacencies( &entity, 1, side_dim, false, adj );
  CHECK_ERR(rval);
  
  int sub_ent_indices[4];
  CN::SubEntityVertexIndices( TYPE_FROM_HANDLE(entity), side_dim, side_num, 
                                sub_ent_indices );
  EntityType subtype = CN::SubEntityType( TYPE_FROM_HANDLE(entity),
                                              side_dim, side_num );
  int sub_ent_corners = CN::VerticesPerEntity(subtype);
  
  const EntityHandle* conn;
  int conn_len;
  rval = moab.get_connectivity( entity, conn, conn_len );
  CHECK_ERR(rval);
  
  for (size_t i = 0; i < adj.size(); ++i) {
    if (TYPE_FROM_HANDLE(adj[i]) != subtype)
      continue;
  
    const EntityHandle* sub_conn;
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
      
  
void check_adj_ho_nodes( Interface& moab,
                         EntityHandle entity )
{
    EntityType type = TYPE_FROM_HANDLE(entity);
    const EntityHandle* conn;
    int conn_len;
    ErrorCode rval = moab.get_connectivity( entity, conn, conn_len );
    CHECK_ERR(rval);
    
    int ho[4];
    CN::HasMidNodes( type, conn_len, ho );
    for (int dim = CN::Dimension(type)-1; dim > 0; --dim) {
      if (!ho[dim])
        continue;
        
      for (int j = 0; j < CN::NumSubEntities( type, dim ); ++j) {
        EntityHandle side = find_side( moab, entity, dim, j );
        if (!side)
          continue;
        
        const EntityHandle* side_conn;
        int side_len;
        rval = moab.get_connectivity( side, side_conn, side_len );
        CHECK_ERR(rval);
        
        int this_idx = CN::HONodeIndex( type, conn_len, dim, j );
        int side_idx = CN::HONodeIndex( TYPE_FROM_HANDLE(side), side_len, dim, 0 );
        CHECK_EQUAL( side_conn[side_idx], conn[this_idx] );
      }
    }
}

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( Interface& moab, 
                       EntityHandle entity,
                       int mid_nodes[4] )
{
    // get element info
  const EntityType type = TYPE_FROM_HANDLE(entity);
  const EntityHandle* conn;
  int conn_len;
  ErrorCode rval = moab.get_connectivity( entity, conn, conn_len );
  CHECK_ERR(rval);
  std::vector<double> coords(3*conn_len);
  rval = moab.get_coords( conn, conn_len, &coords[0] );
  CHECK_ERR(rval);
  
    // calculate and verify expected number of mid nodes
  int num_nodes = CN::VerticesPerEntity(type);
  for (int d = 1; d <= CN::Dimension(type); ++d)
    if (mid_nodes[d])
      num_nodes += CN::NumSubEntities(type, d);
  CHECK_EQUAL( num_nodes, conn_len );
  
    // verify that each higher-order node is at the center
    // of its respective sub-entity.
  for (int i = CN::VerticesPerEntity(type); i < num_nodes; ++i) {
      // get sub-entity owning ho-node  
    int sub_dim, sub_num;
    CN::HONodeParent( type, num_nodes, i, sub_dim, sub_num );
      // get corner vertex indices
    int sub_conn[8], num_sub;
    if (sub_dim < CN::Dimension(type)) {
      CN::SubEntityVertexIndices( type, sub_dim, sub_num, sub_conn );
      EntityType sub_type = CN::SubEntityType( type, sub_dim, sub_num );
      num_sub = CN::VerticesPerEntity( sub_type );
    }
    else {
      num_sub = CN::VerticesPerEntity(type);
      for (int j = 0; j < num_sub; ++j)
        sub_conn[j] = j;
    }
      // calculate mean of corner vertices
    double mean[3] = {0,0,0};
    for (int j = 0; j < num_sub; ++j) {
      int co = 3*sub_conn[j];
      mean[0] += coords[co  ];
      mean[1] += coords[co+1];
      mean[2] += coords[co+2];
    }
    mean[0] /= num_sub;
    mean[1] /= num_sub;
    mean[2] /= num_sub;
      // verify that higher-order node is at expected location
    CHECK_REAL_EQUAL( mean[0], coords[3*i  ], 1e-6 );
    CHECK_REAL_EQUAL( mean[1], coords[3*i+1], 1e-6 );
    CHECK_REAL_EQUAL( mean[2], coords[3*i+2], 1e-6 );
  }
}

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( EntityType type, int num_nodes )
{
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, ho_file );

  ErrorCode rval;
  Tag ho_tag, block_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, block_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( HAS_MID_NODES_TAG_NAME, 4, MB_TYPE_INTEGER, ho_tag );
  CHECK_ERR(rval);
  
  // get material sets with expected higher-order nodes
  Range blocks;
  int ho_flags[4];
  CN::HasMidNodes( type, num_nodes, ho_flags );
  Tag tags[2] = {ho_tag, block_tag};
  void* vals[2] = {ho_flags, NULL};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, blocks );
  CHECK_ERR(rval);
  
  Range::iterator i;
  Range entities;
  for (i = blocks.begin(); i != blocks.end(); ++i) {
    rval = mb.get_entities_by_type( *i, type, entities, true );
    CHECK_ERR(rval);
  }
    // test file should contain all types of HO entities
  CHECK(!entities.empty());
    // test ho node positions -- expect to be at center of adj corners
  for (i = entities.begin(); i != entities.end(); ++i)
    check_ho_element( mb, *i, ho_flags );
    // test ho node handles consistent with adjacent entities
  for (i = entities.begin(); i != entities.end(); ++i)
    check_adj_ho_nodes( mb, *i );
}

void test_multiple_files()
{
    // load two surface meshes, one at z=+5 at z=-5.
  ErrorCode rval;
  Core mb_impl;
  Interface& mb = mb_impl;
  Range file1_ents, file2_ents;
  read_file( mb, input_file_1 );
  mb.get_entities_by_handle( 0, file1_ents );
  read_file( mb, input_file_1 );
  mb.get_entities_by_handle( 0, file2_ents );
  file2_ents = subtract( file2_ents, file1_ents );
  EntityHandle file1, file2;
  mb.create_meshset( MESHSET_SET, file1 );
  mb.create_meshset( MESHSET_SET, file2 );
  mb.add_entities( file1, file1_ents );
  mb.add_entities( file2, file2_ents );
  
    // first check that we get the same number of verts from 
    // each file and that they are distinct vertices
  Range file1_verts, file2_verts;
  rval = mb.get_entities_by_type( file1, MBVERTEX, file1_verts );
  CHECK_ERR(rval);
  CHECK( !file1_verts.empty() );
  rval = mb.get_entities_by_type( file2, MBVERTEX, file2_verts );
  CHECK_ERR(rval);
  CHECK( !file2_verts.empty() );
  CHECK_EQUAL( file1_verts.size(), file2_verts.size() );
  CHECK( intersect( file1_verts,  file2_verts ).empty() );
  
    // now check that we get the same number of elements from 
    // each file and that they are distinct
  Range file1_elems, file2_elems;
  rval = mb.get_entities_by_dimension( file1, 3, file1_elems );
  CHECK_ERR(rval);
  CHECK( !file1_elems.empty() );
  rval = mb.get_entities_by_dimension( file2, 3, file2_elems );
  CHECK_ERR(rval);
  CHECK( !file2_elems.empty() );
  CHECK_EQUAL( file1_elems.size(), file2_elems.size() );
  CHECK( intersect( file1_elems,  file2_elems ).empty() );

    // now check that the connectivity for each element is
    // defined using the appropriate vertex instances
  Range file1_elem_verts, file2_elem_verts;
  rval = mb.get_adjacencies( file1_elems, 0, false, file1_elem_verts, Interface::UNION );
  CHECK_ERR(rval);
  rval = mb.get_adjacencies( file2_elems, 0, false, file2_elem_verts, Interface::UNION );
  CHECK_ERR(rval);
  CHECK_EQUAL( file1_elem_verts.size(), file2_elem_verts.size() );
  CHECK( intersect( file1_elem_verts,  file1_verts ) == file1_elem_verts );
  CHECK( intersect( file2_elem_verts,  file2_verts ) == file2_elem_verts );
}

void test_cubit12() 
{
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, cubit12_file);
}

void test_cubit14()
{
  Core mb_impl;
  Interface& mb = mb_impl;
  read_file( mb, cubit14_file);
  // check the global id for some geometry sets
  GeomTopoTool gtt(&mb_impl);
  Range ranges[5];
  ErrorCode rval = gtt.find_geomsets(ranges);
  CHECK_ERR(rval);
  EntityHandle set0=ranges[0][0]; // does it have a global id > 0?
  Tag gid_tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, gid_tag ); CHECK_ERR(rval);

  int val;
  rval = mb.tag_get_data(gid_tag, &set0, 1, &val );
  CHECK ( val!=0 );
}
