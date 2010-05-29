#include "TestUtil.hpp"
#include "moab/Core.hpp"
#define IS_BUILDING_MB
#include "moab/Range.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/io/three.smf";
#else
static const char example[] = "three.smf";
#endif

void read_file( Interface& moab, const char* input_file );
void test_read_nodes();
void test_read_triangles();

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_read_nodes);
  result += RUN_TEST(test_read_triangles);
  
  return result;
}


void read_file( Interface& moab, const char* input_file )
{
  ErrorCode rval = moab.load_file( input_file );
  CHECK_ERR(rval);
}

void test_read_nodes()
{
  const double eps = 1e-10;
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> nodes;
  rval = mb.get_entities_by_type( 0, MBVERTEX, nodes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)24, nodes.size() );
  
  std::vector<double> coords(3*nodes.size());
  rval = mb.get_coords( &nodes[0], nodes.size(), &coords[0] );
  CHECK_ERR(rval);
  
  int idx = 0;
  CHECK_REAL_EQUAL( coords[3*idx+0], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
  
  ++idx;
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );
 
  idx = 8; // second cube 
  CHECK_REAL_EQUAL( coords[3*idx+0], 1.2, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 0.0, eps );

  idx = 15; // last node of second cube
  CHECK_REAL_EQUAL( coords[3*idx+0], 2.2, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 1.0, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], 1.0, eps );
 
  idx = 16; // first node of third cube 
  CHECK_REAL_EQUAL( coords[3*idx+0], -1.2, eps );
  CHECK_REAL_EQUAL( coords[3*idx+1], 0.5, eps );
  CHECK_REAL_EQUAL( coords[3*idx+2], -0.2071067812, eps );
  
}  


void test_read_triangles()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  read_file( moab, example );
  
  std::vector<EntityHandle> triangles;
  rval = mb.get_entities_by_type( 0, MBTRI, triangles );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)36, triangles.size() );
  
  int vtx_ids[3];
  const EntityHandle* conn;
  int len;
 
  const int conn1[] = { 1, 4, 2}; 
  int pos =0;  
  rval = mb.get_connectivity( triangles[pos], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, len );
  int i=0;
  for ( i=0; i<3; i++)
       vtx_ids[i] = mb.id_from_handle(conn[i]);
  CHECK_ARRAYS_EQUAL( conn1, 3, vtx_ids, len );
  

  // last triangle
  const int conn2[] = { 19, 21, 23 };
  rval = mb.get_connectivity( triangles[35], conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, len );
  for ( i=0; i<3; i++)
       vtx_ids[i] = mb.id_from_handle(conn[i]);
  CHECK_ARRAYS_EQUAL( conn2, 3, vtx_ids, len );

}  

