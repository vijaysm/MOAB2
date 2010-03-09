#include "TestUtil.hpp"
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBCN.hpp"
#define IS_BUILDING_MB
#include "ReadNASTRAN.hpp"
#include "MBRange.hpp"
#include "FileOptions.hpp"
#include <math.h>
#include <algorithm>

#ifdef SRCDIR
static const char example[] = STRINGIFY(SRCDIR) "/test.nas";
#else
static const char example[] = "test.nas";
#endif

void read_file( MBInterface& moab, const char* input_file );
void test_read_nodes();
void test_read_tets();
void test_read_prisms();
void test_read_hexes();
void test_read_material_set1();
void test_read_material_set2();

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_read_nodes);
  result += RUN_TEST(test_read_tets);
  result += RUN_TEST(test_read_prisms);
  result += RUN_TEST(test_read_hexes);
  result += RUN_TEST(test_read_material_set1);
  result += RUN_TEST(test_read_material_set2);
  
  return result;
}


void read_file( MBInterface& moab, const char* input_file )
{
  MBErrorCode rval;
  ReadNASTRAN reader( &moab );
  FileOptions opts("");
  rval = reader.load_file( input_file, 0, opts, 0, 0, 0 );
  CHECK_ERR(rval);
}

void test_read_nodes()
{
  const double eps = 1e-100;
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  read_file( moab, example );
  
  std::vector<MBEntityHandle> nodes;
  rval = mb.get_entities_by_type( 0, MBVERTEX, nodes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)19, nodes.size() );
  
  MBTag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(nodes.size());
  rval = mb.tag_get_data( id_tag, &nodes[0], nodes.size(), &ids[0] );
  CHECK_ERR(rval);
  
  std::vector<int> sorted_ids( ids );
  std::sort( sorted_ids.begin(), sorted_ids.end() );
  
  std::vector<double> coords(3*nodes.size());
  rval = mb.get_coords( &nodes[0], nodes.size(), &coords[0] );
  CHECK_ERR(rval);
  
  int idx, pos = 0;
  // shared between 2 tets and 2 prisms
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], -2.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], -1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], -1.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], -1.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], -1.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], -1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], -1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );
  // hex element
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 5.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 5.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 5.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 5.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 10.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 10.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 10.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 5.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 10.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 10.0, eps );

}  


void test_read_tets()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  read_file( moab, example );
  
  std::vector<MBEntityHandle> tets;
  rval = mb.get_entities_by_type( 0, MBTET, tets );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)2, tets.size() );
  
  MBTag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(tets.size());
  rval = mb.tag_get_data( id_tag, &tets[0], tets.size(), &ids[0] );
  CHECK_ERR(rval);
  
  if (ids[0] != 1) {
    std::swap( ids[0], ids[1] );
    std::swap( tets[0], tets[1] );
  }
  
  int vtx_ids[4];
  const MBEntityHandle* conn;
  int len;
  
  const int conn1[] = { 8, 9, 10, 11 };
  int pos = 0;
  CHECK_EQUAL( pos+1, ids[pos] );
  rval = mb.get_connectivity( tets[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn1, 4, vtx_ids, len );
  
  const int conn2[] = { 4, 3, 2, 1 };
  ++pos;
  CHECK_EQUAL( pos+1, ids[pos] );
  rval = mb.get_connectivity( tets[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn2, 4, vtx_ids, len );
}  

void test_read_prisms()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  read_file( moab, example );
  
  std::vector<MBEntityHandle> prisms;
  rval = mb.get_entities_by_type( 0, MBPRISM, prisms );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)2, prisms.size() );
  
  MBTag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(prisms.size());
  rval = mb.tag_get_data( id_tag, &prisms[0], prisms.size(), &ids[0] );
  CHECK_ERR(rval);
  
  if (ids[0] != 3) {
    std::swap( ids[0], ids[1] );
    std::swap( prisms[0], prisms[1] );
  }
  
  int vtx_ids[6];
  const MBEntityHandle* conn;
  int len;
  
  const int conn1[] = { 2, 3, 4, 5, 6, 7 };
  int pos = 0;
  // Element ids 1 and 2 are the two tet elements.
  // Element ids 3 and 4 are the two prism elements.
  CHECK_EQUAL( pos+3, ids[pos] );
  rval = mb.get_connectivity( prisms[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 6, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn1, 6, vtx_ids, len );
  
  const int conn2[] = { 5, 6, 7, 8, 9, 10 };
  ++pos;
  CHECK_EQUAL( pos+3, ids[pos] );
  rval = mb.get_connectivity( prisms[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 6, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn2, 6, vtx_ids, len );
}  


void test_read_hexes()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  read_file( moab, example );
  
  std::vector<MBEntityHandle> hexes;
  rval = mb.get_entities_by_type( 0, MBHEX, hexes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)1, hexes.size() );
  
  MBTag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(hexes.size());
  rval = mb.tag_get_data( id_tag, &hexes[0], hexes.size(), &ids[0] );
  CHECK_ERR(rval);
  
  int vtx_ids[8];
  const MBEntityHandle* conn;
  int len;
  
  const int conn1[] = { 12, 13, 14, 15, 16, 17, 18, 19 };
  int pos = 0;
  // Element id 5 is the hex
  CHECK_EQUAL( pos+5, ids[pos] );
  rval = mb.get_connectivity( hexes[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn1, 8, vtx_ids, len );
} 
 
// The tets are in material set 1.
void test_read_material_set1()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  read_file( moab, example );
  
  MBTag mat_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, mat_tag );
  CHECK_ERR(rval);
  
  MBRange mat_set_one;  
  const int one = 1;
  const void* const one_val[] = {&one};  
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &mat_tag, one_val, 1, mat_set_one );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)mat_set_one.size() );
  
  std::vector<MBEntityHandle> tets, contents;
  rval = mb.get_entities_by_type( 0, MBTET, tets );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( mat_set_one.front(), contents );
  CHECK_ERR(rval);
  std::sort( tets.begin(), tets.end() );
  std::sort( contents.begin(), contents.end() );
  CHECK_EQUAL( tets, contents );
}

// The prisms are in material set 2.
void test_read_material_set2()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  read_file( moab, example );
  
  MBTag mat_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, mat_tag );
  CHECK_ERR(rval);
  
  MBRange mat_set_two;  
  const int two = 2;
  const void* const two_val[] = {&two};  
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &mat_tag, two_val, 1, mat_set_two );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)mat_set_two.size() );
  
  std::vector<MBEntityHandle> prisms, contents;
  rval = mb.get_entities_by_type( 0, MBPRISM, prisms );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( mat_set_two.front(), contents );
  CHECK_ERR(rval);
  std::sort( prisms.begin(), prisms.end() );
  std::sort( contents.begin(), contents.end() );
  CHECK_EQUAL( prisms, contents );
}
