/** \file   h5portable.cpp
 *  \author Jason Kraftcheck 
 *  \date   2010-09-23
 * 
 * Tests for HDF5 portability because we call H5Tconvert ourselves
 * to work around parallel performance issues in the HDF5 library.
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "TestRunner.hpp"
#include "ReadHDF5.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace moab;

#define FILE_NAME_BASE "portable"
#ifdef SRCDIR
#  define BASE_NAME STRINGIFY(SRCDIR) "/" FILE_NAME_BASE
#else 
#  define BASE_NAME FILE_NAME_BASE
#endif
#define NATIVE_NAME FILE_NAME_BASE ".h5m"
#define READ_OPTS "BUFFER_SIZE=256"

const int default_int_tag[] = { -100, -99 };
const int mesh_int_tag[] = {1134, -1134};
const double default_real_tag[] = { -1, -2, -3, -4, -5, -6, -7, -8 };
const double mesh_real_tag[] = { 8, 7, 6, 5, 4, 3, 2, 1 };

const size_t INTERVALS = 8;
const size_t NUM_VERT = (INTERVALS+1)*(INTERVALS+1);
const size_t NUM_QUAD = INTERVALS * INTERVALS;
const double Z = 0.0;
const double EPS = std::numeric_limits<double>::epsilon();

// This function is not actually called from any test.
// It is the utility code used to generate the initial
// input file.
//
// Create a mesh with the following properties:
//   A 8x8 grid of quadrilateral elements in the Z plane
//   Lower left vertex is at origin
//   Vertex coordinates are incremented by 1 in increasing X and Y directions
//   Each element's connectivity is ordered starting with the lower-left vertex
//      counter-clockwise about the positive Z axis.
//   Each element is tagged with the integer value of it's lower-left vertex 
//      coordinates ("int_tag", default = [-100,-99], mesh = [1134,-1134])
//   Each element is tagged with the array of its vertex coordinates (interleaved)
//      ("real_tag", default = [-1, -2, -3, -4, -5, -6, -7, -8], mesh = [8, 7, 5, 5, 4, 3, 2, 1])
//   Each vertex is tagged with the handle of the element for which it is the
//      lower left vertex ("handle_tag")
//   A set containing each row of quads that is a child of the row below it
//   Each set tagged with a bit tag indicating whether it is an even or odd row,
//     with the first row even. ("bit_tag");
//   Each set will also be tagged with an integer tag named 'rowset' containing
//     the same value as the "bit_tag" for use in testing partial reads.
//
//  |      |      |      |      |      |      |      |      |
//  |  17  |  18  |  19  |  20  |  21  |  22  |  23  |  24  |
//  |      |      |      |      |      |      |      |      |
// 19-----20-----21-----22-----23-----24-----25-----26-----27
//  |      |      |      |      |      |      |      |      |
//  |   9  |  10  |  11  |  12  |  13  |  14  |  15  |  16  |
//  |      |      |      |      |      |      |      |      |
// 10-----11-----12-----13-----14-----15-----16-----17-----18
//  |      |      |      |      |      |      |      |      |
//  |   1  |   2  |   3  |   4  |   5  |   6  |   7  |   8  |
//  |      |      |      |      |      |      |      |      |
//  1------2------3------4------5------6------7------8------9
void create_mesh( const char* filename )
{
  EntityHandle verts[NUM_VERT];
  EntityHandle quads[NUM_QUAD];
  
  ErrorCode rval;
  Core core;
  Interface& mb = core;
  const EntityHandle root = 0;
  
  Tag int_tag, real_tag, handle_tag, bit_tag, row_tag;

  rval = mb.tag_get_handle( "int_tag", 2, MB_TYPE_INTEGER, int_tag, MB_TAG_DENSE|MB_TAG_EXCL, default_int_tag );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( int_tag, &root, 1, mesh_int_tag );
  CHECK_ERR(rval);

  rval = mb.tag_get_handle( "real_tag", 8, MB_TYPE_DOUBLE, real_tag, MB_TAG_DENSE|MB_TAG_EXCL, default_real_tag );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( real_tag, &root, 1, mesh_real_tag );
  CHECK_ERR(rval);
  
  rval = mb.tag_get_handle( "handle_tag", 1, MB_TYPE_HANDLE, handle_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
  rval = mb.tag_get_handle( "bit_tag", 1, MB_TYPE_BIT, bit_tag, MB_TAG_EXCL );
  CHECK_ERR(rval);

  rval = mb.tag_get_handle( "rowset", 1, MB_TYPE_INTEGER, row_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
  
  for (size_t i = 0; i < NUM_VERT; ++i) {
    double coords[] = { i % 9, i / 9, Z };
    rval = mb.create_vertex( coords, verts[i] );
    CHECK_ERR(rval);
  }
  
  for (size_t i = 0; i < NUM_QUAD; ++i) {
    size_t j = (i / 8) * 9 + i % 8;
    EntityHandle conn[] = { verts[j], verts[j+1], verts[j+10], verts[j+9] };
    rval = mb.create_element( MBQUAD, conn, 4, quads[i] );
    CHECK_ERR(rval);
    
    double coords[4*3];
    rval = mb.get_coords( conn, 4, coords );
    CHECK_ERR(rval);
    
    int int_val[] = { (int)coords[0], (int)coords[1] };
    rval = mb.tag_set_data( int_tag, quads+i, 1, int_val );
    CHECK_ERR(rval);
    
    double real_val[] = { coords[0], coords[1],
                          coords[3], coords[4],
                          coords[6], coords[7],
                          coords[9], coords[10] };
    rval = mb.tag_set_data( real_tag, quads+i, 1, real_val );
    CHECK_ERR(rval);
    
    rval = mb.tag_set_data( handle_tag, conn, 1, quads+i );
    CHECK_ERR(rval);
  }
  
  EntityHandle prev = 0;
  for (size_t i = 0; i < INTERVALS; ++i) {
    EntityHandle set;
    int flag = i < 4 ? MESHSET_SET : MESHSET_ORDERED;
    rval = mb.create_meshset( flag, set );
    CHECK_ERR(rval);
    
    rval = mb.add_entities( set, quads + i*INTERVALS, INTERVALS );
    CHECK_ERR(rval);
    
    char bit = i % 2;
    rval = mb.tag_set_data( bit_tag, &set, 1, &bit );
    CHECK_ERR(rval);
    
    int ival = bit;
    rval = mb.tag_set_data( row_tag, &set, 1, &ival );
    CHECK_ERR(rval);
    
    if (prev != 0) {
      rval = mb.add_parent_child( prev, set );
      CHECK_ERR(rval);
    }    
    prev = set;
  }
  
  rval = mb.write_file( filename, "MOAB" );
  CHECK_ERR(rval);
}

void test_load_file( const char* filename )
{
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
}

void test_read_vertices( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
    // get the vertices
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK_EQUAL( NUM_VERT, verts.size() );
    // check vertex coordinates
  bool seen[INTERVALS+1][INTERVALS+1];
  for (size_t i = 0; i <= INTERVALS; ++i)
    for (size_t j = 0; j <= INTERVALS; ++j)
      seen[i][j] = false;
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL( Z, coords[2], EPS );
    size_t x = (size_t)coords[0];
    size_t y = (size_t)coords[1];
    CHECK_REAL_EQUAL( (double)x, coords[0], EPS );
    CHECK_REAL_EQUAL( (double)y, coords[1], EPS );
    CHECK( x >= 0 && x <= INTERVALS );
    CHECK( y >= 0 && y <= INTERVALS );
    CHECK( !seen[x][y] );
    seen[x][y] = true;
  }
}

void test_elements( Interface& mb, bool odd_only )
{
    // get the elements
  Range quads;
  ErrorCode rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  CHECK_EQUAL( NUM_QUAD/(1+odd_only), quads.size() );
    // check quad connectivity
  bool seen[INTERVALS][INTERVALS];
  for (size_t i = 0; i < INTERVALS; ++i)
    for (size_t j = 0; j < INTERVALS; ++j)
      seen[i][j] = false;
  for (Range::iterator i = quads.begin(); i != quads.end(); ++i) {
    const EntityHandle* conn = 0;
    int len = 0;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 4, len );
    double coords[4*3];
    rval = mb.get_coords( conn, len, coords );
    CHECK_ERR(rval);

    size_t x = (size_t)coords[0];
    size_t y = (size_t)coords[1];
    CHECK_REAL_EQUAL( (double)x, coords[0], EPS );
    CHECK_REAL_EQUAL( (double)y, coords[1], EPS );
    CHECK( x >= 0 && x < INTERVALS );
    CHECK( y >= 0 && y < INTERVALS );
    
    CHECK_REAL_EQUAL( Z,     coords[2], EPS );
    CHECK_REAL_EQUAL( 1.0+x, coords[3], EPS );
    CHECK_REAL_EQUAL( 0.0+y, coords[4], EPS );
    CHECK_REAL_EQUAL( Z,     coords[5], EPS );
    CHECK_REAL_EQUAL( 1.0+x, coords[6], EPS );
    CHECK_REAL_EQUAL( 1.0+y, coords[7], EPS );
    CHECK_REAL_EQUAL( Z,     coords[8], EPS );
    CHECK_REAL_EQUAL( 0.0+x, coords[9], EPS );
    CHECK_REAL_EQUAL( 1.0+y, coords[10], EPS );
    CHECK_REAL_EQUAL( Z,     coords[11], EPS );
  }
  
  if (odd_only) {
    for (size_t i = 0; i < INTERVALS; i += 2)
      for (size_t j = 0; j < INTERVALS; ++j)
        CHECK(!seen[i][j]);
  }
}

void test_read_elements( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
    // check the elements
  test_elements( mb, false );
}

void read_sets( Interface& mb, EntityHandle rows[INTERVALS] )
{
    // get the sets
  Range sets;
  ErrorCode rval = mb.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( INTERVALS, sets.size() );
    // check set contents
  for (Range::iterator i = sets.begin(); i != sets.end(); ++i) {
    Range contents;
    rval = mb.get_entities_by_handle( *i, contents );
    CHECK_ERR(rval);
    CHECK_EQUAL( INTERVALS, contents.size() );
    CHECK( contents.all_of_type(MBQUAD) );
    
    const EntityHandle* conn = 0;
    int len = 0;
    rval = mb.get_connectivity( contents.front(), conn, len );
    CHECK_ERR(rval);
    double coords[3];
    rval = mb.get_coords( conn, 1, coords );
    CHECK_ERR(rval);
    
    size_t y = (size_t)coords[1];
    CHECK( y >= 0 && y < INTERVALS );
    rows[y] = *i;
  }
}

void test_read_set_contents( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
    // get the sets
  EntityHandle rows[INTERVALS];
  read_sets( mb, rows );
    // check set contents
  for (size_t i = 0; i < INTERVALS; ++i) {
    Range contents;
    rval = mb.get_entities_by_handle( rows[i], contents );
    CHECK_ERR(rval);
    CHECK_EQUAL( INTERVALS, contents.size() );
    CHECK( contents.all_of_type(MBQUAD) );
    
    for (Range::iterator j = contents.begin(); j != contents.end(); ++j) {
      const EntityHandle* conn = 0;
      int len = 0;
      rval = mb.get_connectivity( *j, conn, len );
      CHECK_ERR(rval);
      double coords[3];
      rval = mb.get_coords( conn, 1, coords );
      CHECK_ERR(rval);
    
      size_t y = (size_t)coords[1];
      CHECK_EQUAL( i, y );
    }
  }
}

void test_read_set_parent_child( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
    // get the sets
  EntityHandle rows[INTERVALS];
  read_sets( mb, rows );
    // check set children
  std::vector<EntityHandle> list;
  for (size_t i = 0; i < INTERVALS-1; ++i) {
    list.clear();
    rval = mb.get_child_meshsets( rows[i], list );
    CHECK_ERR(rval);
    CHECK_EQUAL( (size_t)1, list.size() );
    CHECK_EQUAL( rows[i+1], list.front() );
  }
  list.clear();
  rval = mb.get_child_meshsets( rows[INTERVALS-1], list );
  CHECK_ERR(rval);
  CHECK( list.empty() );
    // check set parents
  list.clear();
  rval = mb.get_parent_meshsets( rows[0], list );
  CHECK_ERR(rval);
  CHECK(list.empty());
  for (size_t i = 1; i < INTERVALS; ++i) {
    list.clear();
    rval = mb.get_parent_meshsets( rows[i], list );
    CHECK_ERR(rval);
    CHECK_EQUAL( (size_t)1, list.size() );
    CHECK_EQUAL( rows[i-1], list.front() );
  }
}

void test_read_int_tag( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
  
    // get the tag
  Tag tag;
  rval = mb.tag_get_handle( "int_tag", 2, MB_TYPE_INTEGER, tag );
  CHECK_ERR(rval);
  
    // check tag properties
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, size );
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  DataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TYPE_INTEGER, type );
  
    // check default value
  int value[2];
  rval = mb.tag_get_default_value( tag, value );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( default_int_tag, 2, value, 2 );
  
    // check mesh value
  const EntityHandle root = 0;
  rval = mb.tag_get_data( tag, &root, 1, value );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( mesh_int_tag, 2, value, 2 );
  
    // check entity values
  Range quads;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  for (Range::iterator i = quads.begin(); i != quads.end(); ++i) {
    const EntityHandle* conn = 0;
    int len = 0;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    double coords[3];
    rval = mb.get_coords( conn, 1, coords );
    CHECK_ERR(rval);
    int exp[] = { (int)coords[0], (int)coords[1] };

    rval = mb.tag_get_data( tag, &*i, 1, value );
    CHECK_ERR(rval);
    CHECK_ARRAYS_EQUAL( exp, 2, value, 2 );
  }
}

void test_read_real_tag( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
  
    // get the tag
  Tag tag;
  rval = mb.tag_get_handle( "real_tag", 8, MB_TYPE_DOUBLE, tag );
  CHECK_ERR(rval);
  
    // check tag properties
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 8, size );
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  DataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TYPE_DOUBLE, type );
  
    // check default value
  double value[8];
  rval = mb.tag_get_default_value( tag, value );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( default_real_tag, 8, value, 8 );
  
    // check mesh value
  const EntityHandle root = 0;
  rval = mb.tag_get_data( tag, &root, 1, value );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( mesh_real_tag, 8, value, 8 );
  
    // check entity values
  Range quads;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  for (Range::iterator i = quads.begin(); i != quads.end(); ++i) {
    const EntityHandle* conn = 0;
    int len = 0;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( 4, len );
    double coords[3*4];
    rval = mb.get_coords( conn, len, coords );
    CHECK_ERR(rval);
    double exp[] = { coords[0], coords[1],
                     coords[3], coords[4],
                     coords[6], coords[7],
                     coords[9], coords[10] };

    rval = mb.tag_get_data( tag, &*i, 1, value );
    CHECK_ERR(rval);
    CHECK_ARRAYS_EQUAL( exp, 8, value, 8 );
  }
}

void test_read_handle_tag( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
  
    // get the tag
  Tag tag;
  rval = mb.tag_get_handle( "handle_tag", 1, MB_TYPE_HANDLE, tag );
  CHECK_ERR(rval);
  
    // check tag properties
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, size );
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_SPARSE, storage );
  DataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TYPE_HANDLE, type );
  
    // check entity values
  Range quads;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  for (Range::iterator i = quads.begin(); i != quads.end(); ++i) {
    const EntityHandle* conn = 0;
    int len = 0;
    rval = mb.get_connectivity( *i, conn, len );
    CHECK_ERR(rval);
    
    EntityHandle value;
    rval = mb.tag_get_data( tag, conn, 1, &value );
    CHECK_ERR(rval);
    CHECK_EQUAL( *i, value );    
  }
}

void test_read_bit_tag( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  ErrorCode rval = mb.load_file( filename, 0, READ_OPTS );
  CHECK_ERR(rval);
  
    // get the tag
  Tag tag;
  rval = mb.tag_get_handle( "bit_tag", 1, MB_TYPE_BIT, tag );
  CHECK_ERR(rval);
  
    // check tag properties
  int size;
  rval = mb.tag_get_length( tag, size );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, size );
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_BIT, storage );
  DataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TYPE_BIT, type );
  
    // check entity values
  EntityHandle sets[INTERVALS];
  read_sets( mb, sets );
  char values[INTERVALS], expected[INTERVALS];
  rval = mb.tag_get_data( tag, sets, INTERVALS, values );
  CHECK_ERR(rval);
  for (size_t i = 0; i < INTERVALS; ++i) 
    expected[i] = (char)(i%2);
  CHECK_ARRAYS_EQUAL( expected, INTERVALS, values, INTERVALS );
}

void test_read_partial( const char* filename )
{
    // load the file
  Core core;
  Interface& mb = core;
  const int odd = 1;
  ErrorCode rval = mb.load_file( filename, 0, "CHILDREN=NONE;SETS=NONE;" READ_OPTS, "rowset", &odd, 1 );
  CHECK_ERR(rval);
    // check the elements
  test_elements( mb, true );
}

// Make sure everything works without data conversion just to verify
// that the tests themselves are valid
void test_native_read()
{
  const char* filename = NATIVE_NAME;
  create_mesh( filename );
  test_read_vertices( filename );
  test_read_elements( filename );
  test_read_set_contents( filename );
  test_read_set_parent_child( filename );
  test_read_int_tag( filename );
  test_read_real_tag( filename );
  test_read_handle_tag( filename );
  test_read_bit_tag( filename );
  test_read_partial( filename );
  if (!getenv("KEEP_FILE"))
    remove(filename);
}

#define DEFINE_TEST_SET( NAME ) \
  void test_load_file_ ## NAME () { test_load_file( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_vertices_ ## NAME () { test_read_vertices( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_elements_ ## NAME () { test_read_elements( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_set_contents_ ## NAME () { test_read_set_contents( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_set_parent_child_ ## NAME () { test_read_set_parent_child( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_int_tag_ ## NAME () { test_read_int_tag( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_real_tag_ ## NAME () { test_read_real_tag( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_handle_tag_ ## NAME () { test_read_handle_tag( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_bit_tag_ ## NAME () { test_read_bit_tag( BASE_NAME "_" #NAME ".h5m" ); } \
  void test_read_partial_ ## NAME () { test_read_partial( BASE_NAME "_" #NAME ".h5m" ); } 
  
#define REGISTER_TEST_SET( NAME ) \
  REGISTER_DEP_TEST( test_load_file_ ## NAME , test_native_read ); \
  REGISTER_DEP_TEST( test_read_vertices_ ## NAME , test_load_file_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_elements_ ## NAME , test_read_vertices_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_set_contents_ ## NAME , test_read_elements_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_set_parent_child_ ## NAME , test_read_set_contents_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_int_tag_ ## NAME , test_read_elements_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_real_tag_ ## NAME , test_read_elements_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_handle_tag_ ## NAME ,test_read_elements_ ## NAME  ); \
  REGISTER_DEP_TEST( test_read_bit_tag_ ## NAME , test_read_set_contents_ ## NAME ); \
  REGISTER_DEP_TEST( test_read_partial_ ## NAME , test_read_elements_ ## NAME ); 

DEFINE_TEST_SET( x86_64 )
DEFINE_TEST_SET( x86_32 )
DEFINE_TEST_SET( power_32 )

int main( int argc, char* argv[] )
{ 
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif
  REGISTER_TEST( test_native_read );
  REGISTER_TEST_SET( x86_64 );
  REGISTER_TEST_SET( x86_32 );
  REGISTER_TEST_SET( power_32 );
  int result = RUN_TESTS( argc, argv );
#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif
  return result;
}
