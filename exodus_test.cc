#include "TestUtil.hpp"
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBCN.hpp"
#include "ReadNCDF.hpp"
#include "WriteNCDF.hpp"
#include "FileOptions.hpp"
#include <math.h>
#include <algorithm>

/* Input test file: ho_test.g
 * 
 * File is expected to contain at least one block for every
 * supported higher-order element type.  The coordinates of
 * every higher-order node are expected to be the mean of the
 * adjacent corner vertices of the element.
 */
#ifdef SRCDIR
static const char ho_file[] = STRINGIFY(SRCDIR) "/test/ho_test.g";
#else
static const char ho_file[] = "test/ho_test.g";
#endif

void read_file( MBInterface& moab, 
                const char* input_file );

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( MBInterface& moab, 
                       MBEntityHandle entity,
                       int mid_nodes[4] );

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( MBEntityType type, int num_nodes );

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

int main()
{
  int result = 0;
  
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
  
  return result;
}

void read_file( MBInterface& moab, 
                const char* input_file )
{
  MBErrorCode rval;
  MBEntityHandle set;
  ReadNCDF reader( &moab );
  FileOptions opts("");
  rval = reader.load_file( input_file, set, opts, 0, 0 );
  CHECK_ERR(rval);
}

void write_and_read( MBInterface& write_mb,
                     MBInterface& read_mb,
                     MBEntityHandle block = 0 )
{
  const char* tmp_file = "exodus_test_tmp.g";
  MBErrorCode rval;
  MBEntityHandle set;
  ReadNCDF reader( &read_mb );
  WriteNCDF writer( &write_mb );
  FileOptions opts("");
  
  MBEntityHandle* write_set_list = &block;
  int write_set_list_len = 0;//(block != 0);
  std::vector<std::string> qa_records;
  rval = writer.write_file( tmp_file, true, opts, 
                            write_set_list, write_set_list_len,
                            qa_records, NULL, 0, 3 );
  if (MB_SUCCESS != rval) 
    remove(tmp_file);
  CHECK_ERR(rval);
  
  rval = reader.load_file( tmp_file, set, opts, 0, 0 );
  remove( tmp_file );
  CHECK_ERR(rval);
}

void check_ho_elements( MBInterface& moab, 
                        MBEntityHandle block,
                        MBEntityType type,
                        int mid_nodes[4] )
{
  MBErrorCode rval;
  MBRange elems;
  rval = moab.get_entities_by_handle( block, elems );
  CHECK_ERR(rval);
  CHECK(!elems.empty());
  CHECK(elems.all_of_type(type));
  for (MBRange::const_iterator i = elems.begin(); i != elems.end(); ++i)
    check_ho_element( moab, *i, mid_nodes );
}

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( MBInterface& moab, 
                       MBEntityHandle entity,
                       int mid_nodes[4] )
{
    // get element info
  const MBEntityType type = TYPE_FROM_HANDLE(entity);
  const MBEntityHandle* conn;
  int conn_len;
  MBErrorCode rval = moab.get_connectivity( entity, conn, conn_len );
  CHECK_ERR(rval);
  std::vector<double> coords(3*conn_len);
  rval = moab.get_coords( conn, conn_len, &coords[0] );
  CHECK_ERR(rval);
  
    // calculate and verify expected number of mid nodes
  int num_nodes = MBCN::VerticesPerEntity(type);
  for (int d = 1; d <= MBCN::Dimension(type); ++d)
    if (mid_nodes[d])
      num_nodes += MBCN::NumSubEntities(type, d);
  CHECK_EQUAL( num_nodes, conn_len );
  
    // verify that each higher-order node is at the center
    // of its respective sub-entity.
  for (int i = MBCN::VerticesPerEntity(type); i < num_nodes; ++i) {
      // get sub-entity owning ho-node  
    int sub_dim, sub_num;
    MBCN::HONodeParent( type, num_nodes, i, sub_dim, sub_num );
      // get corner vertex indices
    int sub_conn[8], num_sub;
    if (sub_dim < MBCN::Dimension(type)) {
      MBCN::SubEntityVertexIndices( type, sub_dim, sub_num, sub_conn );
      MBEntityType sub_type = MBCN::SubEntityType( type, sub_dim, sub_num );
      num_sub = MBCN::VerticesPerEntity( sub_type );
    }
    else {
      num_sub = MBCN::VerticesPerEntity(type);
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


MBEntityHandle find_block( MBInterface& mb, MBEntityType type, const int has_mid_nodes[4] )
{

  MBErrorCode rval;
  MBTag ho_tag, block_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, block_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( HAS_MID_NODES_TAG_NAME, ho_tag );
  CHECK_ERR(rval);
  
  // get material sets with expected higher-order nodes
  MBRange blocks;
  MBTag tags[2] = {ho_tag, block_tag};
  const void* vals[2] = {has_mid_nodes, NULL};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, blocks );
  CHECK_ERR(rval);
  
  for (MBRange::iterator i = blocks.begin(); i != blocks.end(); ++i) {
    int n;
    rval = mb.get_number_entities_by_type( *i, type, n );
    CHECK_ERR(rval);
    if (n > 0)
      return *i;
  }
  
  CHECK(false); // no block matching element type description
  return 0;
}

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( MBEntityType type, int num_nodes )
{
  MBCore mb_impl1, mb_impl2;
  MBInterface &mb1 = mb_impl1, &mb2 = mb_impl2;
  int ho_flags[4];
  MBCN::HasMidNodes( type, num_nodes, ho_flags );

    // read file 
  read_file( mb1, ho_file );
    // test element connectivity order
  MBEntityHandle block = find_block( mb1, type, ho_flags );
  CHECK(block != 0);
  check_ho_elements( mb1, block, type, ho_flags );
  
    // write block and read it back in
  write_and_read( mb1, mb2, block );
    // test element connectivity order on re-read data
  block = find_block( mb2, type, ho_flags );
  CHECK(block != 0);
  check_ho_elements( mb2, block, type, ho_flags );
}
