#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#define IS_BUILDING_MB
#include "ReadIDEAS.hpp"
#include "moab/Range.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/test.unv";
#else
static const char example[] = "test.unv";
#endif

void read_file( Interface& moab, const char* input_file );
void test_read_nodes();
void test_read_tets();
void test_read_hexes();
void test_read_material_set();
void test_read_physical_set();

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_read_nodes);
  result += RUN_TEST(test_read_tets);
  result += RUN_TEST(test_read_hexes);
  result += RUN_TEST(test_read_material_set);
  result += RUN_TEST(test_read_physical_set);
  
  return result;
}


void read_file( Interface& moab, const char* input_file )
{
  ErrorCode rval = moab.load_file( input_file );
  CHECK_ERR(rval);
}

void test_read_nodes()
{
  const double eps = 1e-6;
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> nodes;
  rval = mb.get_entities_by_type( 0, MBVERTEX, nodes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)17, nodes.size() );
  
  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
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
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
  
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
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );
  
  // id=4
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );

  // id=8
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
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  // id=12
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0000022760448197, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], -2.2760448196157412e-06, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], -2.2760448196157412e-6, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0000022760448197, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], -2.2760448197267635e-6, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], -2.2760448196157412e-6, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );

  // id=16
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 5.0e-1, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 5.0e-1, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 3.0, eps );

  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0000022760448197, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0000022760448197, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 2.0, eps );
}  


void test_read_tets()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> tets;
  rval = mb.get_entities_by_type( 0, MBTET, tets );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)2, tets.size() );
  
  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(tets.size());
  rval = mb.tag_get_data( id_tag, &tets[0], tets.size(), &ids[0] );
  CHECK_ERR(rval);
  
  if (ids[0] != 3) {
    std::swap( ids[0], ids[1] );
    std::swap( tets[0], tets[1] );
  }
  
  int vtx_ids[4];
  const EntityHandle* conn;
  int len;
  
  // The first tet has id=3
  const int conn1[] = { 13, 14, 15, 16 };
  int pos = 0, offset = 3;
  CHECK_EQUAL( pos+offset, ids[pos] );
  rval = mb.get_connectivity( tets[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn1, 4, vtx_ids, len );
  
  // The second tet has id=4
  const int conn2[] = { 13, 17, 14, 16 };
  ++pos;
  CHECK_EQUAL( pos+offset, ids[pos] );
  rval = mb.get_connectivity( tets[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn2, 4, vtx_ids, len );
}  

void test_read_hexes()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> hexes;
  rval = mb.get_entities_by_type( 0, MBHEX, hexes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)2, hexes.size() );
  
  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(hexes.size());
  rval = mb.tag_get_data( id_tag, &hexes[0], hexes.size(), &ids[0] );
  CHECK_ERR(rval);
  
  int vtx_ids[8];
  const EntityHandle* conn;
  int len;
  
  const int conn1[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  int pos = 0, offset = 1;
  // Element id 1 is a hex
  CHECK_EQUAL( pos+offset, ids[pos] );
  rval = mb.get_connectivity( hexes[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn1, 8, vtx_ids, len );
  
  const int conn2[] = { 2, 9, 10, 3, 6, 11, 12, 7 };
  ++pos;
  // Element id 2 is a hex
  CHECK_EQUAL( pos+offset, ids[pos] );
  rval = mb.get_connectivity( hexes[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn2, 8, vtx_ids, len );
} 

// Two tets and two hexes are in material set 100.
void test_read_material_set()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  Tag mat_tag;
  rval = mb.tag_get_handle( MAT_PROP_TABLE_TAG, 1, MB_TYPE_INTEGER, mat_tag );
  CHECK_ERR(rval);
  
  Range mat_set;  
  const int mat_set_id = 100;
  const void* const mat_set_id_val[] = {&mat_set_id};  
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &mat_tag, mat_set_id_val, 
                                          1, mat_set );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)mat_set.size() );

  std::vector<EntityHandle> elements, contents;
  rval = mb.get_entities_by_type( 0, MBTET, elements );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type( 0, MBHEX, elements );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( mat_set.front(), contents );
  CHECK_ERR(rval);
  std::sort( elements.begin(), elements.end() );
  std::sort( contents.begin(), contents.end() );
  CHECK_EQUAL( elements, contents );
}

// The tets are in physical set 4, which corresponds to volume 4 in Cubit.
void test_read_physical_set()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  Tag phys_tag;
  rval = mb.tag_get_handle( PHYS_PROP_TABLE_TAG, 1, MB_TYPE_INTEGER, phys_tag );
  CHECK_ERR(rval);
  
  Range phys_set;  
  const int phys_set_id = 4;
  const void* const phys_set_id_val[] = {&phys_set_id};  
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &phys_tag, phys_set_id_val,
                                          1, phys_set );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)phys_set.size() );
  
  std::vector<EntityHandle> tets, contents;
  rval = mb.get_entities_by_type( 0, MBTET, tets );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( phys_set.front(), contents );
  CHECK_ERR(rval);
  std::sort( tets.begin(), tets.end() );
  std::sort( contents.begin(), contents.end() );
  CHECK_EQUAL( tets, contents );
}
