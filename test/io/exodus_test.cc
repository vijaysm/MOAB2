#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/MBTagConventions.hpp"
#include "moab/CN.hpp"
#define IS_BUILDING_MB
#include "ReadNCDF.hpp"
#include "WriteNCDF.hpp"
#include "FileOptions.hpp"
#include "ExoIIUtil.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

/* Input test file: ho_test.g
 * 
 * File is expected to contain at least one block for every
 * supported higher-order element type.  The coordinates of
 * every higher-order node are expected to be the mean of the
 * adjacent corner vertices of the element.
 */
#ifdef SRCDIR
static const char ho_file[] = STRINGIFY(SRCDIR) "/ho_test.g";
#else
static const char ho_file[] = "ho_test.g";
#endif

void read_file( Interface& moab, 
                const char* input_file );

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( Interface& moab, 
                       EntityHandle entity,
                       int mid_nodes[4] );

void test_read_side( int sideset_id,
                     EntityType sideset_type,
                     int sideset_nodes_per_elem,
                     bool shell_side = false );

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( EntityType type, int num_nodes );

void test_types();

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

void test_read_tri6_side()    { test_read_side( 1, MBEDGE, 3 ); }  // sideset 1
void test_read_shell_side()   { test_read_side( 3, MBQUAD, 9, true ); } // sideset 3
void test_read_shell_edge()   { test_read_side( 4, MBEDGE, 3 ); } // sideset 4
void test_read_hex20_side()   { test_read_side( 2, MBQUAD, 8 ); }  // sideset 2

void test_read_block_ids();
void test_read_sideset_ids();
void test_read_nodeset_ids();

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_types);
  
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
  
  result += RUN_TEST(test_read_tri6_side );
  result += RUN_TEST(test_read_shell_side);
  result += RUN_TEST(test_read_shell_edge);
  result += RUN_TEST(test_read_hex20_side);
  
  result += RUN_TEST(test_read_block_ids );
  result += RUN_TEST(test_read_sideset_ids);
  result += RUN_TEST(test_read_nodeset_ids);
  
  return result;
}

struct TestType {
  EntityType moab_type;
  ExoIIElementType exo_type;
  int num_nodes;
  std::string name;
};

void check_type( const TestType& type )
{
  int has_mid_nodes[4];
  CN::HasMidNodes( type.moab_type, type.num_nodes, has_mid_nodes );
  
  CHECK_EQUAL( type.moab_type, ExoIIUtil::ExoIIElementMBEntity[type.exo_type] );
  CHECK_EQUAL( type.name, std::string( ExoIIUtil::ElementTypeNames[type.exo_type] ) );
  CHECK_EQUAL( type.num_nodes, ExoIIUtil::VerticesPerElement[type.exo_type] );
  switch (CN::Dimension(type.moab_type)) {
    case 3: CHECK_EQUAL( has_mid_nodes[3], ExoIIUtil::HasMidNodes[type.exo_type][3] );
    case 2: CHECK_EQUAL( has_mid_nodes[2], ExoIIUtil::HasMidNodes[type.exo_type][2] );
    case 1: CHECK_EQUAL( has_mid_nodes[1], ExoIIUtil::HasMidNodes[type.exo_type][1] );
  }
  
  Core moab;
  ExoIIUtil tool(&moab);
  CHECK_EQUAL( type.exo_type, tool.element_name_to_type( type.name.c_str() ) );
  CHECK_EQUAL( type.name, std::string(tool.element_type_name( type.exo_type ) ) );
}
  
void test_types()
{
  const TestType types[] = {
    { MBVERTEX,  EXOII_SPHERE,     1, "SPHERE" },
    { MBEDGE,    EXOII_SPRING,     1, "SPRING" },
    { MBEDGE,    EXOII_BAR,        2, "BAR" },
    { MBEDGE,    EXOII_BAR2,       2, "BAR2" },
    { MBEDGE,    EXOII_BAR3,       3, "BAR3" },
    { MBEDGE,    EXOII_BEAM,       2, "BEAM" },
    { MBEDGE,    EXOII_BEAM2,      2, "BEAM2" },
    { MBEDGE,    EXOII_BEAM3,      3, "BEAM3" },
    { MBEDGE,    EXOII_TRUSS,      2, "TRUSS" },
    { MBEDGE,    EXOII_TRUSS2,     2, "TRUSS2" },
    { MBEDGE,    EXOII_TRUSS3,     3, "TRUSS3" },
    { MBTRI,     EXOII_TRI,        3, "TRI" },
    { MBTRI,     EXOII_TRI3,       3, "TRI3" },
    { MBTRI,     EXOII_TRI6,       6, "TRI6" },
    { MBTRI,     EXOII_TRI7,       7, "TRI7" },
    { MBQUAD,    EXOII_QUAD,       4, "QUAD" },
    { MBQUAD,    EXOII_QUAD4,      4, "QUAD4" },
    { MBQUAD,    EXOII_QUAD5,      5, "QUAD5" },
    { MBQUAD,    EXOII_QUAD8,      8, "QUAD8" },
    { MBQUAD,    EXOII_QUAD9,      9, "QUAD9" },
    { MBQUAD,    EXOII_SHELL,      4, "SHELL" },
    { MBQUAD,    EXOII_SHELL4,     4, "SHELL4" },
    { MBQUAD,    EXOII_SHELL5,     5, "SHELL5" },
    { MBQUAD,    EXOII_SHELL8,     8, "SHELL8" },
    { MBQUAD,    EXOII_SHELL9,     9, "SHELL9" },
    { MBTET,     EXOII_TETRA,      4, "TETRA" },
    { MBTET,     EXOII_TETRA4,     4, "TETRA4" },
    { MBTET,     EXOII_TETRA8,     8, "TETRA8" },
    { MBTET,     EXOII_TETRA10,   10, "TETRA10" },
    { MBTET,     EXOII_TETRA14,   14, "TETRA14" },
    { MBPYRAMID, EXOII_PYRAMID,    5, "PYRAMID" },
    { MBPYRAMID, EXOII_PYRAMID5,   5, "PYRAMID5" },
    { MBPYRAMID, EXOII_PYRAMID10, 10, "PYRAMID10" },
    { MBPYRAMID, EXOII_PYRAMID13, 13, "PYRAMID13" },
    { MBPYRAMID, EXOII_PYRAMID18, 18, "PYRAMID18" },
    { MBPRISM,   EXOII_WEDGE,      6, "WEDGE" },
    { MBKNIFE,   EXOII_KNIFE,      7, "KNIFE" },
    { MBHEX,     EXOII_HEX,        8, "HEX" },
    { MBHEX,     EXOII_HEX8,       8, "HEX8" },
    { MBHEX,     EXOII_HEX9,       9, "HEX9" },
    { MBHEX,     EXOII_HEX20,     20, "HEX20" },
    { MBHEX,     EXOII_HEX27,     27, "HEX27" },
    { MBHEX,     EXOII_HEXSHELL,  12, "HEXSHELL" } };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (int i = 0; i < num_types; ++i) 
    check_type( types[i] );
}

void read_file( Interface& moab, 
                const char* input_file )
{
  ErrorCode rval;
  ReadNCDF reader( &moab );
  FileOptions opts("");
  rval = reader.load_file( input_file, 0, opts, 0, 0, 0 );
  CHECK_ERR(rval);
}

void write_and_read( Interface& write_mb,
                     Interface& read_mb,
                     EntityHandle block = 0 )
{
  const char* tmp_file = "exodus_test_tmp.g";
  ErrorCode rval;
  ReadNCDF reader( &read_mb );
  WriteNCDF writer( &write_mb );
  FileOptions opts("");
  
  EntityHandle* write_set_list = &block;
  int write_set_list_len = 0;//(block != 0);
  std::vector<std::string> qa_records;
  rval = writer.write_file( tmp_file, true, opts, 
                            write_set_list, write_set_list_len,
                            qa_records, NULL, 0, 3 );
  if (MB_SUCCESS != rval) 
    remove(tmp_file);
  CHECK_ERR(rval);
  
  rval = reader.load_file( tmp_file, 0, opts, 0, 0, 0 );
  remove( tmp_file );
  CHECK_ERR(rval);
}

void check_ho_elements( Interface& moab, 
                        EntityHandle block,
                        EntityType type,
                        int mid_nodes[4] )
{
  ErrorCode rval;
  Range elems;
  rval = moab.get_entities_by_handle( block, elems );
  CHECK_ERR(rval);
  CHECK(!elems.empty());
  CHECK(elems.all_of_type(type));
  for (Range::const_iterator i = elems.begin(); i != elems.end(); ++i)
    check_ho_element( moab, *i, mid_nodes );
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


EntityHandle find_block( Interface& mb, EntityType type, const int has_mid_nodes[4] )
{

  ErrorCode rval;
  Tag ho_tag, block_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, block_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( HAS_MID_NODES_TAG_NAME, ho_tag );
  CHECK_ERR(rval);
  
  // get material sets with expected higher-order nodes
  Range blocks;
  Tag tags[2] = {ho_tag, block_tag};
  const void* vals[2] = {has_mid_nodes, NULL};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, blocks );
  CHECK_ERR(rval);
  
  for (Range::iterator i = blocks.begin(); i != blocks.end(); ++i) {
    int n;
    rval = mb.get_number_entities_by_type( *i, type, n );
    CHECK_ERR(rval);
    if (n > 0)
      return *i;
  }
  
  CHECK(false); // no block matching element type description
  return 0;
}

EntityHandle find_sideset( Interface& mb, 
                             int sideset_id,
                             EntityType side_type )
{
  ErrorCode rval;
  Tag ss_tag;
  rval = mb.tag_get_handle( NEUMANN_SET_TAG_NAME, ss_tag );
  CHECK_ERR(rval);
  
  const void* tag_vals[] = { &sideset_id };
  Range side_sets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &ss_tag, tag_vals, 1, side_sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)side_sets.size() );
  return side_sets.front();
}

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( EntityType type, int num_nodes )
{
  Core mb_impl1, mb_impl2;
  Interface &mb1 = mb_impl1, &mb2 = mb_impl2;
  int ho_flags[4];
  CN::HasMidNodes( type, num_nodes, ho_flags );

    // read file 
  read_file( mb1, ho_file );
    // test element connectivity order
  EntityHandle block = find_block( mb1, type, ho_flags );
  CHECK(block != 0);
  check_ho_elements( mb1, block, type, ho_flags );
  
    // write block and read it back in
  write_and_read( mb1, mb2, block );
    // test element connectivity order on re-read data
  block = find_block( mb2, type, ho_flags );
  CHECK(block != 0);
  check_ho_elements( mb2, block, type, ho_flags );
}

void test_read_side( int id,
                     EntityType sideset_type,
                     int sideset_nodes_per_elem,
                     bool shell_side )
{
  // read test file
  Core mb_impl;
  Interface& moab = mb_impl;
  read_file( moab, ho_file );
  
  // get side set 
  EntityHandle set = find_sideset( moab, id, sideset_type );
  CHECK(set != 0);
  
  // check expected element connectivity
  int ho_flags[4];
  CN::HasMidNodes( sideset_type, sideset_nodes_per_elem, ho_flags );
  check_ho_elements( moab, set, sideset_type, ho_flags );
  
  if (shell_side)
    return;
  
  // check that each element is on the boundary of the mesh
  Range elems;
  ErrorCode rval = mb_impl.get_entities_by_handle( set, elems );
  CHECK_ERR(rval);
  
  int dim = CN::Dimension( sideset_type );
  for (Range::iterator i= elems.begin(); i != elems.end(); ++i) {
    Range adj;
    rval = mb_impl.get_adjacencies( &*i, 1, dim+1, false, adj, Interface::UNION );
    CHECK_ERR(rval);
    CHECK_EQUAL( 1, (int)adj.size() );
  }
    
}

void test_read_ids_common( const char* file_name,
                           const char* tag_name,
                           const int* expected_vals,
                           int num_expected )
{
  Core mb;
  ReadNCDF reader(&mb);
  
  FileOptions opts("");
  std::vector<int> values;
  ErrorCode rval = reader.read_tag_values( file_name, tag_name, opts, values );
  CHECK_ERR(rval);
  
  std::vector<int> expected( expected_vals, expected_vals+num_expected );
  std::sort( values.begin(), values.end() );
  std::sort( expected.begin(), expected.end() );
  CHECK_EQUAL( expected, values );
}

void test_read_block_ids() {
  const int expected[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
  test_read_ids_common( ho_file, MATERIAL_SET_TAG_NAME, expected, sizeof(expected)/sizeof(expected[0]) );
}

void test_read_sideset_ids() {
  const int expected[] = { 1, 2, 3, 4 };
  test_read_ids_common( ho_file, NEUMANN_SET_TAG_NAME, expected, sizeof(expected)/sizeof(expected[0]) );
}

void test_read_nodeset_ids() {
  test_read_ids_common( ho_file, DIRICHLET_SET_TAG_NAME, 0, 0 );
}

