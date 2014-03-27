#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

/* Input test file: test/sample.stl
 * 
 * Example ASCII STL file from: 
 *  http://people.sc.fsu.edu/~burkardt/data/stla/stla.html
 */
#ifdef MESHDIR
static const char sample[] = STRINGIFY(MESHDIR) "/io/sample.stl";
#else
static const char sample[] = "sample.stl";
#endif

// Use static keyword so that tmp_file can only be accessed within this file
static const char* tmp_file = "test.stl";

void test_read_ascii();
void test_write_ascii();
void test_type_option();
void test_detect_type();
void test_endian_option();
void test_big_endian();
void test_little_endian();
void test_detect_byte_order();

void read_file( Interface& moab, 
                const char* input_file,
                const char* options = "" );
void convert_file( const char* source_file,
                   const char* dest_file,
                   const char* options = "" );
// check that the mesh constains the simple tetrahedron defined
// in test/sample.stl
void check_mesh_is_tet( Interface& moab );

int main()
{
  int result = 0;
  
  result += RUN_TEST(test_read_ascii);
  result += RUN_TEST(test_write_ascii);
  result += RUN_TEST(test_type_option);
  result += RUN_TEST(test_detect_type);
  result += RUN_TEST(test_endian_option);
  result += RUN_TEST(test_big_endian);
  result += RUN_TEST(test_little_endian);
  result += RUN_TEST(test_detect_byte_order);
  
  remove( tmp_file );
  return result;
}

ErrorCode read_file_( Interface& moab, const char* input_file, const char* options = "" )
{
  ErrorCode rval = moab.load_file( input_file, 0, options );
  return rval;
}

void read_file( Interface& moab, const char* input_file, const char* options )
{
  ErrorCode rval = read_file_( moab, input_file, options );
  CHECK_ERR(rval);
}

void convert_file( const char* input_file, const char* output_file, const char* options )
{
  ErrorCode rval;
  Core moab;

  rval = moab.load_file( input_file );
  CHECK_ERR(rval);
  
  rval = moab.write_file( output_file, "STL", options );
  CHECK_ERR(rval);
}

void test_read_ascii()
{
  Core moab;
  read_file( moab, sample, "ASCII" );
  check_mesh_is_tet( moab );
}

void test_write_ascii()
{
  convert_file( sample, tmp_file, "ASCII" );
  Core moab;
  read_file( moab, tmp_file, "ASCII" );
  remove( tmp_file );
  check_mesh_is_tet( moab );
}

void test_type_option()
{
  ErrorCode rval;
  Core moab;

  rval = read_file_( moab, sample, "BINARY" );
  CHECK( MB_SUCCESS != rval );
  
  convert_file( sample, tmp_file, "BINARY" );
  rval = read_file_( moab, tmp_file, "ASCII" );
  CHECK( MB_SUCCESS != rval );
  
  remove( tmp_file );
}

void test_detect_type()
{
  Core moab;

  read_file( moab, sample );
  
  convert_file( sample, tmp_file, "BINARY" );
  read_file_( moab, tmp_file );
  
  remove( tmp_file );
}

void test_endian_option()
{
  ErrorCode rval;
  Core moab;

  convert_file( sample, tmp_file, "BINARY;BIG_ENDIAN" );
  rval = read_file_( moab, tmp_file, "BINARY;LITTLE_ENDIAN" );
  CHECK( MB_SUCCESS != rval );
  
  convert_file( sample, tmp_file, "BINARY;LITTLE_ENDIAN" );
  rval = read_file_( moab, tmp_file, "BINARY;BIG_ENDIAN" );
  CHECK( MB_SUCCESS != rval );
  
  remove( tmp_file );
}

void test_big_endian()
{
  Core moab;
  convert_file( sample, tmp_file, "BINARY;BIG_ENDIAN" );
  read_file( moab, tmp_file, "BINARY;BIG_ENDIAN" );
  check_mesh_is_tet( moab );
  remove( tmp_file );
}

void test_little_endian()
{
  Core moab;
  convert_file( sample, tmp_file, "BINARY;LITTLE_ENDIAN" );
  read_file( moab, tmp_file, "BINARY;LITTLE_ENDIAN" );
  check_mesh_is_tet( moab );
  remove( tmp_file );
}

void test_detect_byte_order()
{
  Core moab;

  convert_file( sample, tmp_file, "BINARY;LITTLE_ENDIAN" );
  read_file( moab, tmp_file, "BINARY" );

  convert_file( sample, tmp_file, "BINARY;BIG_ENDIAN" );
  read_file( moab, tmp_file, "BINARY" );
  
  remove( tmp_file );
}


void check_mesh_is_tet( Interface& moab )
{
  ErrorCode rval;
  Range verts, tris, other;
  rval = moab.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  rval = moab.get_entities_by_type( 0, MBTRI, tris );
  CHECK_ERR(rval);
  rval = moab.get_entities_by_handle( 0, other );
  CHECK_ERR(rval);
  
  CHECK_EQUAL( 4, (int)verts.size() );
  CHECK_EQUAL( 4, (int)tris.size() );
  other = subtract( other, verts );
  other = subtract( other, tris );
  CHECK( other.all_of_type(MBENTITYSET) );
  
  const double expt_coords[4][3] = { { 0, 0, 0 },
                                     { 1, 0, 0 },
                                     { 0, 1, 0 },
                                     { 0, 0, 1 } };
  EntityHandle vert_handles[4] = { 0, 0, 0, 0 };
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = moab.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    
    bool found = false;
    for (int j = 0; j < 4; ++j) {
      double ds = 0;
      for (int d = 0; d < 3; ++d) {
        double dl = expt_coords[j][d] - coords[d];
        ds += dl*dl;
      }
      
      if (ds < 1e-6) {
        CHECK_EQUAL( (EntityHandle)0, vert_handles[j] );
        vert_handles[j] = *i;
        found = true;
        break;
      }
    }
    CHECK(found);
  }
  
  const int expt_conn[4][3] = { { 0, 1, 3 },
                                { 0, 2, 1 },
                                { 0, 3, 2 },
                                { 1, 2, 3 } };
  EntityHandle tri_handles[4] = { 0, 0, 0, 0 };
  for (Range::iterator i = tris.begin(); i != tris.end(); ++i) {
    const EntityHandle* conn = 0;
    int len = 0;
    rval = moab.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 3, len );
    
    int conn_idx[3] = { 
        static_cast<int>(std::find( vert_handles, vert_handles + 4, conn[0] ) - vert_handles),
        static_cast<int>(std::find( vert_handles, vert_handles + 4, conn[1] ) - vert_handles),
        static_cast<int>(std::find( vert_handles, vert_handles + 4, conn[2] ) - vert_handles) };
    CHECK( conn_idx[0] != 4 );
    CHECK( conn_idx[1] != 4 );
    CHECK( conn_idx[2] != 4 );
    
    bool found = false;
    for (int j = 0; j < 4; ++j) {
      int k = std::find( expt_conn[j], expt_conn[j]+3, conn_idx[0] ) - expt_conn[j];
      if (k == 3)
        continue;
      
      if (expt_conn[j][(k+1)%3] == conn_idx[1] && 
          expt_conn[j][(k+2)%3] == conn_idx[2]) {
        CHECK_EQUAL( (EntityHandle)0, tri_handles[j] );
        tri_handles[j] = *i;
        found = true;
        break;
      }
    }
    CHECK(found);
  }
}
