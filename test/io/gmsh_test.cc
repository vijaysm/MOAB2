#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

/* Input test file: gmsh2.msh
 * 
 * Example version 2.0 ASCII input file from Gmsh 2.4 manual.
 */
#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/gmsh2.msh";
#else
static const char example[] = "gmsh2.msh";
#endif

void test_read_nodes();
void test_read_quads();
void test_read_material_set();
void test_read_geom_set();

void read_file( Interface& moab, const char* input_file );

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_read_nodes);
  result += RUN_TEST(test_read_quads);
  result += RUN_TEST(test_read_material_set);
  result += RUN_TEST(test_read_geom_set);
  
  return result;
}


void read_file( Interface& moab, const char* input_file )
{
  ErrorCode rval;
  rval = moab.load_file( input_file );
  CHECK_ERR(rval);
}

void test_read_nodes()
{
  const double eps = 1e-100;
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> nodes;
  rval = mb.get_entities_by_type( 0, MBVERTEX, nodes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)6, nodes.size() );
  
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
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
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
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
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
  CHECK_REAL_EQUAL( coords[3*idx+0], 2.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
  
  ++pos;
  CHECK_EQUAL( pos+1, sorted_ids[pos] );
  idx = std::find( ids.begin(), ids.end(),pos+1 ) - ids.begin();
  CHECK_REAL_EQUAL( coords[3*idx+0], 2.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
}  


void test_read_quads()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> quads;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)2, quads.size() );
  
  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
  CHECK_ERR(rval);
  
  std::vector<int> ids(quads.size());
  rval = mb.tag_get_data( id_tag, &quads[0], quads.size(), &ids[0] );
  CHECK_ERR(rval);
  
  if (ids[0] != 1) {
    std::swap( ids[0], ids[1] );
    std::swap( quads[0], quads[1] );
  }
  
  int vtx_ids[4];
  const EntityHandle* conn;
  int len;
  
  const int conn1[] = { 1, 2, 3, 4 };
  int pos = 0;
  CHECK_EQUAL( pos+1, ids[pos] );
  rval = mb.get_connectivity( quads[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn1, 4, vtx_ids, len );
  
  const int conn2[] = { 2, 5, 6, 3 };
  ++pos;
  CHECK_EQUAL( pos+1, ids[pos] );
  rval = mb.get_connectivity( quads[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4, len );
  rval = mb.tag_get_data( id_tag, conn, len, vtx_ids );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( conn2, 4, vtx_ids, len );
}  

void test_read_material_set()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  Tag mat_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, mat_tag );
  CHECK_ERR(rval);
  
  Range sets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &mat_tag, 0, 1, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)sets.size() );
  EntityHandle set = sets.front();
  
  int id;
  rval = mb.tag_get_data( mat_tag, &set, 1, &id );
  CHECK_ERR(rval);
  CHECK_EQUAL( 99, id );
  
  std::vector<EntityHandle> quads, contents;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( set, contents );
  CHECK_ERR(rval);
  std::sort( quads.begin(), quads.end() );
  std::sort( contents.begin(), contents.end() );
  CHECK_EQUAL( quads, contents );
}


void test_read_geom_set()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  Tag dim_tag, id_tag;
  rval = mb.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, dim_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
  CHECK_ERR(rval);
  
  Range sets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &dim_tag, 0, 1, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)sets.size() );
  EntityHandle set = sets.front();
  
  int dim;
  rval = mb.tag_get_data( dim_tag, &set, 1, &dim );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, dim );
  
  int id;
  rval = mb.tag_get_data( id_tag, &set, 1, &id );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, id );
  
  std::vector<EntityHandle> quads, contents;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_handle( set, contents );
  CHECK_ERR(rval);
  std::sort( quads.begin(), quads.end() );
  std::sort( contents.begin(), contents.end() );
  CHECK_EQUAL( quads, contents );
}


