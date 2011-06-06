#include "moab/Core.hpp"
#include "moab/Range.hpp"

using namespace moab;

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>

#define DECLARE_TEST(A) \
  bool test_ ## A(); \
  int A ## _reg_var = register_test( &test_ ## A, #A );

typedef bool (*test_ptr)();
struct test_data { test_ptr test; const char* name; bool result; };
size_t num_tests = 0;
test_data *test_array = 0;
int register_test( test_ptr test, const char* name )
{
  test_array = (test_data*)realloc( test_array, sizeof(test_data)*(num_tests+1) );
  test_array[num_tests].test = test;
  test_array[num_tests].name = name;
  test_array[num_tests].result = true;
  ++num_tests;
  return 0;
}   

DECLARE_TEST(edge2)
DECLARE_TEST(edge3)
DECLARE_TEST(tri3)
DECLARE_TEST(tri6)
DECLARE_TEST(quad4)
DECLARE_TEST(quad8)
DECLARE_TEST(quad9)
DECLARE_TEST(polygon)
DECLARE_TEST(tet4)
DECLARE_TEST(tet10)
DECLARE_TEST(hex8)
DECLARE_TEST(hex20)
DECLARE_TEST(hex27)
DECLARE_TEST(wedge)
DECLARE_TEST(wedge15)
DECLARE_TEST(pyramid)
DECLARE_TEST(pyramid13)

DECLARE_TEST(structured_points_2d)
DECLARE_TEST(structured_grid_2d)
DECLARE_TEST(rectilinear_grid_2d)
DECLARE_TEST(structured_points_3d)
DECLARE_TEST(structured_grid_3d)
DECLARE_TEST(rectilinear_grid_3d)

DECLARE_TEST(scalar_attrib_1_bit)
DECLARE_TEST(scalar_attrib_1_uchar)
DECLARE_TEST(scalar_attrib_1_char)
DECLARE_TEST(scalar_attrib_1_ushort)
DECLARE_TEST(scalar_attrib_1_short)
DECLARE_TEST(scalar_attrib_1_uint)
DECLARE_TEST(scalar_attrib_1_int)
DECLARE_TEST(scalar_attrib_1_ulong)
DECLARE_TEST(scalar_attrib_1_long)
DECLARE_TEST(scalar_attrib_1_float)
DECLARE_TEST(scalar_attrib_1_double)

DECLARE_TEST(scalar_attrib_4_bit)
DECLARE_TEST(scalar_attrib_4_uchar)
DECLARE_TEST(scalar_attrib_4_char)
DECLARE_TEST(scalar_attrib_4_ushort)
DECLARE_TEST(scalar_attrib_4_short)
DECLARE_TEST(scalar_attrib_4_uint)
DECLARE_TEST(scalar_attrib_4_int)
DECLARE_TEST(scalar_attrib_4_ulong)
DECLARE_TEST(scalar_attrib_4_long)
DECLARE_TEST(scalar_attrib_4_float)
DECLARE_TEST(scalar_attrib_4_double)

DECLARE_TEST(vector_attrib_bit)
DECLARE_TEST(vector_attrib_uchar)
DECLARE_TEST(vector_attrib_char)
DECLARE_TEST(vector_attrib_ushort)
DECLARE_TEST(vector_attrib_short)
DECLARE_TEST(vector_attrib_uint)
DECLARE_TEST(vector_attrib_int)
DECLARE_TEST(vector_attrib_ulong)
DECLARE_TEST(vector_attrib_long)
DECLARE_TEST(vector_attrib_float)
DECLARE_TEST(vector_attrib_double)

DECLARE_TEST(tensor_attrib_uchar)
DECLARE_TEST(tensor_attrib_char)
DECLARE_TEST(tensor_attrib_ushort)
DECLARE_TEST(tensor_attrib_short)
DECLARE_TEST(tensor_attrib_uint)
DECLARE_TEST(tensor_attrib_int)
DECLARE_TEST(tensor_attrib_ulong)
DECLARE_TEST(tensor_attrib_long)
DECLARE_TEST(tensor_attrib_float)
DECLARE_TEST(tensor_attrib_double)

DECLARE_TEST(subset)

DECLARE_TEST(unstructured_field)

int main( int argc, char* argv[] )
{
  int *test_indices = (int*)malloc( sizeof(int) * num_tests );
  int test_count;
    // if no arguments, do all tests
  if (argc == 1) {
    for (unsigned i = 0; i < num_tests; ++i)
      test_indices[i] = i;
    test_count = num_tests;
  }
    // otherwise run only specified tests
  else {
    test_count = 0;
    for (int i = 1; i < argc; ++i) 
      for (unsigned j = 0; j < num_tests; ++j)
        if (!strcmp( test_array[j].name, argv[i]))
          test_indices[test_count++] = j;
  }
  
  int fail_count = 0;
  for (int i = 0; i < test_count; ++i)
  {
    test_data& test = test_array[test_indices[i]];
    printf("Testing %s...\n", test.name );
    if (!(test.result = test.test()))
      ++fail_count;
  }
  
  printf("\n\n");
  if (fail_count) {
    printf("FAILED TESTS:\n");
    for (int i = 0; i < test_count; ++i)
    {
      test_data& test = test_array[test_indices[i]];
      if (!test.result)
        printf("\t%s\n", test.name);
    }
  }
    
  
  if (test_count == 0)
    printf("0 VTK tests run\n");
  else if (fail_count == 0)
    printf("%d tests passed\n", test_count );
  else
    printf("%d of %d tests failed\n", fail_count, test_count );
  printf("\n");
  
  return fail_count;
}

#define CHECK(A) if (is_error((A))) return do_error( #A, __LINE__ )
static bool do_error( const char* string, int line )
{
  fprintf(stderr, "Check failed at line %d: %s\n", line, string );
  return false;
}
static inline bool is_error( bool b )
  { return !b; }

//static bool do_error( ErrorCode err, int line )
//{
//  Core tmp_core;
//  fprintf(stderr, "API failed at line %d: %s (%d)\n", 
//    line, tmp_core.get_error_string(err).c_str(), (int)err );
//  return false;
//}
static inline bool is_error( ErrorCode b )
  { return MB_SUCCESS != b; }

bool read_file( Interface* iface, const char* file );
bool write_and_read( Interface* iface1, Interface* iface2 );


bool test_read_write_element( const double* coords, unsigned num_coords,
                              const int* vtk_conn, const int* moab_conn,
                              unsigned num_conn,
                              unsigned num_elem, unsigned vtk_type,
                              EntityType moab_type );

bool test_edge2()
{
  const double coords[] =
   { 0, 0, 0,
     1, 0, 0,
     1, 1, 0,
     0, 1, 0,
     0, 0, 1 };
  const int conn[] =
    { 0, 1,
      1, 2,
      2, 3, 
      3, 4, 
      0, 4 };
      
  return test_read_write_element( coords, 5, conn, conn, 10, 5, 3, MBEDGE );
}
  
  
bool test_edge3()
{
  const double coords[] =
    { -1, -1, 2,
       1, -1, 2,
       1,  1, 2,
      -1,  1, 2,
       0.000, -0.707, 2,
       0.707,  0.000, 2,
       0.000,  0.707, 2,
      -0.707,  0.000, 2 };
  const int conn[] = 
    { 0, 1, 4,
      1, 2, 5,
      2, 3, 6,
      3, 0, 7 };
  
  return test_read_write_element( coords, 8, conn, conn, 12, 4, 21, MBEDGE );
}


bool test_tri3()
{
  const double coords[] = 
    {  0,  0,  0,
       5,  0,  0,
       0,  5,  0,
      -5,  0,  0,
       0, -5,  0 };
  const int conn[] = 
    { 0, 1, 2,
      0, 2, 3,
      0, 3, 4,
      0, 4, 1 };
  
  return test_read_write_element( coords, 5, conn, conn, 12, 4, 5, MBTRI );
}


bool test_tri6()
{
  const double coords[] = 
    {  0,  2,  0,
       0,  0,  2,
       0, -2,  0,
       0,  0, -2,
       0,  1,  1,
       0, -1,  1,
       0, -1, -1,
       0,  1, -1,
       0,  0,  0 };
  const int conn[] = 
    { 0, 1, 3, 4, 5, 8, 
      1, 2, 3, 6, 7, 8 };
  
  return test_read_write_element( coords, 9, conn, conn, 12, 2, 22, MBTRI );
}

const double grid_3x3[] =
  { 0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    3, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    3, 1, 0,
    0, 2, 0,
    1, 2, 0,
    2, 2, 0,
    3, 2, 0,
    0, 3, 0,
    1, 3, 0,
    2, 3, 0,
    3, 3, 0 };
    
const int quad_structured_conn[] =
  { 0, 1, 5, 4,
    1, 2, 6, 5,
    2, 3, 7, 6,
    4, 5, 9, 8,
    5, 6, 10, 9,
    6, 7, 11, 10,
    8, 9, 13, 12,
    9, 10, 14, 13,
    10, 11, 15, 14 };

bool test_quad4()
{
  // test read as quads
  bool rval1 = test_read_write_element( grid_3x3, 16, quad_structured_conn, quad_structured_conn, 36, 9, 9, MBQUAD );
  
  // test read as pixels
  const int conn2[] = 
    { 0, 1, 4, 5,
      1, 2, 5, 6, 
      2, 3, 6, 7,
      4, 5, 8, 9,
      5, 6, 9, 10,
      6, 7, 10, 11,
      8, 9, 12, 13,
      9, 10, 13, 14,
      10, 11, 14, 15 };
  bool rval2 = test_read_write_element( grid_3x3, 16, conn2, quad_structured_conn, 36, 9, 8, MBQUAD );
  
  return rval1 && rval2;
}

bool test_quad8()
{
  const double coords[] = 
    { 0, 0, 0,
      0, 2, 0,
      0, 4, 0,
      0, 0, 4,
      0, 2, 4, 
      0, 4, 4,
      4, 0, 0,
      4, 2, 0,
      4, 4, 0,
      2, 0, 0,
      2, 4, 0,
      0, 0, 2,
      0, 4, 2
    };
  const int conn[] = 
    { 0, 2, 5, 3, 1, 12, 4, 11,
      2, 0, 6, 8, 1, 9, 7, 10 };
  
  return test_read_write_element( coords, 13, conn, conn, 16, 2, 23, MBQUAD );
}

bool test_quad9()
{
  const double coords[] = 
    { 0, 0, 0,
      0, 2, 0,
      0, 4, 0,
      0, 0, 4,
      0, 2, 4, 
      0, 4, 4,
      4, 0, 0,
      4, 2, 0,
      4, 4, 0,
      2, 0, 0,
      2, 2, 0,
      2, 4, 0,
      0, 0, 2,
      0, 2, 2,
      0, 4, 2
    };
  const int conn[] = 
    { 0, 2, 5, 3, 1, 14, 4, 12, 12,
      2, 0, 6, 8, 1, 9, 7, 11, 10 };
  
  return test_read_write_element( coords, 15, conn, conn, 18, 2, 28, MBQUAD );
}

bool test_polygon()
{
  const double coords[] = 
    { 0, 0, 0,
      0, 2, 0,
      0, 4, 0,
      0, 0, 4,
      0, 2, 4, 
      0, 4, 4,
      4, 0, 0,
      4, 2, 0,
      4, 4, 0,
      2, 0, 0,
      2, 4, 0,
      0, 0, 2,
      0, 4, 2
    };
  const int conn[] = 
    { 0, 1, 2, 12, 5, 4, 3, 11,
      2, 1, 0, 9, 6, 7, 8, 10 };
  
  return test_read_write_element( coords, 13, conn, conn, 16, 2, 7, MBPOLYGON );
}

bool test_tet4()
{
  const double coords[] = 
    {  1, -1,  0,
       1,  1,  0,
      -1,  1,  0,
      -1, -1,  0,
       0,  0, -1,
       0,  0,  1 };
  const int conn[] = 
    { 0, 1, 3, 5,
      1, 2, 3, 5,
      0, 1, 4, 3, 
      1, 2, 4, 3 };
  
  return test_read_write_element( coords, 6, conn, conn, 16, 4, 10, MBTET );
}

bool test_tet10()
{
  const double coords[] =
    { 4, 0, 0,
      0, 2, 0,
      0,-2, 0,
      0, 0,-2,
      0, 0, 2,
      0, 1, 1,
      2, 0, 1,
      0,-1, 1,
      0, 0, 0,
      2, 1, 0,
      2,-1, 0,
      2, 0,-1,
      0,-1,-1,
      0, 1,-1 };
  const int conn[] = 
    { 0, 1, 2, 4,  9, 8, 10, 6, 5, 7,
      2, 1, 0, 3,  8, 9, 10, 12, 13, 11 };
  
  return test_read_write_element( coords, 14, conn, conn, 20, 2, 24, MBTET );
}


const double grid_2x2x2[] =
  { 0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    0, 2, 0,
    1, 2, 0,
    2, 2, 0,
    0, 0, 1,
    1, 0, 1,
    2, 0, 1,
    0, 1, 1,
    1, 1, 1,
    2, 1, 1,
    0, 2, 1,
    1, 2, 1,
    2, 2, 1,
    0, 0, 2,
    1, 0, 2,
    2, 0, 2,
    0, 1, 2,
    1, 1, 2,
    2, 1, 2,
    0, 2, 2,
    1, 2, 2,
    2, 2, 2 };

const int hex_structured_conn[] =
  { 0, 1, 4, 3, 9,10,13,12,
    1, 2, 5, 4,10,11,14,13,
    3, 4, 7, 6,12,13,16,15,
    4, 5, 8, 7,13,14,17,16,
    9,10,13,12,18,19,22,21,
   10,11,14,13,19,20,23,22,
   12,13,16,15,21,22,25,24,
   13,14,17,16,22,23,26,25 };

bool test_hex8()
{
    // check vtk hexes
  bool rval1 = test_read_write_element( grid_2x2x2, 27, hex_structured_conn, hex_structured_conn, 64, 8, 12, MBHEX );
  CHECK(rval1);
  
  const int conn2[] = 
    { 0, 1, 3, 4, 9,10,12,13,
      1, 2, 4, 5,10,11,13,14,
      3, 4, 6, 7,12,13,15,16,
      4, 5, 7, 8,13,14,16,17,
      9,10,12,13,18,19,21,22,
     10,11,13,14,19,20,22,23,
     12,13,15,16,21,22,24,25,
     13,14,16,17,22,23,25,26 };
  
    // check with vtk voxels
  bool rval2 = test_read_write_element( grid_2x2x2, 27, conn2, hex_structured_conn, 64, 8, 11, MBHEX );
  CHECK(rval2);
  
  return true;
}

bool test_hex20()
{
  const int vtk_conn[] = 
    {  0,  2,  8,  6, 
      18, 20, 26, 24, 
       1,  5,  7,  3, 
      19, 23, 25, 21,
       9, 11, 17, 15 };
  const int exo_conn[] =
    {  0,  2,  8,  6, 
      18, 20, 26, 24, 
       1,  5,  7,  3, 
       9, 11, 17, 15,
      19, 23, 25, 21 };
  
  return test_read_write_element( grid_2x2x2, 27, vtk_conn, exo_conn, 20, 1, 25, MBHEX );       
}

bool test_hex27()
{
  const int vtk_conn[] = 
    {  0,  2,  8,  6, 
      18, 20, 26, 24, 
       1,  5,  7,  3, 
      19, 23, 25, 21,
       9, 11, 17, 15,
      10, 16, 14, 12,
       4, 22, 13 };
  const int moab_conn[] =
    {  0,  2,  8,  6, 
      18, 20, 26, 24, 
       1,  5,  7,  3, 
       9, 11, 17, 15,
      19, 23, 25, 21,
      14, 16, 12, 10,
       4, 22, 13 };
  
  return test_read_write_element( grid_2x2x2, 27, vtk_conn, moab_conn, 27, 1, 29, MBHEX );       
}

bool test_wedge()
{
  const double coords[] =
    { 1, 0, 0,
      1, 1, 0,
      0, 1, 0,
      0, 0, 0,
      1, 0, 1,
      1, 1, 1,
      0, 1, 1,
      0, 0, 1 };
  const int exo_conn[] = 
    { 0, 1, 3, 4, 5, 7,
      1, 2, 3, 5, 6, 7 };
  const int vtk_conn[] =
    { 0, 3, 1, 4, 7, 5,
      1, 3, 2, 5, 7, 6 };
  return test_read_write_element( coords, 8, vtk_conn, exo_conn, 12, 2, 13, MBPRISM ); 
}

bool test_wedge15()
{
  const double coords[] =
    { 2, 0, 0,
      2, 2, 0,
      0, 2, 0,
      0, 0, 0,
      2, 0, 2,
      2, 2, 2,
      0, 2, 2,
      0, 0, 2,
      2, 1, 0,
      1, 2, 0,
      0, 1, 0,
      1, 0, 0,
      1, 1, 0,
      2, 1, 2,
      1, 2, 2,
      0, 1, 2,
      1, 0, 2,
      1, 1, 2,
      2, 0, 1,
      2, 2, 1,
      0, 2, 1,
      0, 0, 1 };
  const int exo_conn[] = 
    { 0, 1, 3, 4, 5, 7, 8, 12, 11, 18, 19, 21, 13, 17, 16,
      1, 2, 3, 5, 6, 7, 9, 10, 12, 19, 20, 21, 14, 15, 17 };
  const int vtk_conn[] =
    { 0, 3, 1, 4, 7, 5, 11, 12, 8, 16, 17, 13, 18, 21, 19,
      1, 3, 2, 5, 7, 6, 12, 10, 9, 17, 15, 14, 19, 21, 20 };
  return test_read_write_element( coords, 22, vtk_conn, exo_conn, 30, 2, 26, MBPRISM ); 
}

bool test_pyramid()
{
  const double coords[] = 
    {  1, -1,  0,
       1,  1,  0,
      -1,  1,  0,
      -1, -1,  0,
       0,  0, -1,
       0,  0,  1 };
  const int conn[] = 
    { 0, 1, 2, 3, 5,
      3, 2, 1, 0, 4 };
  
  return test_read_write_element( coords, 6, conn, conn, 10, 2, 14, MBPYRAMID );
}

bool test_pyramid13()
{
  const double coords[] = 
    {  2, -2,  0,
       2,  2,  0,
      -2,  2,  0,
      -2, -2,  0,
       0,  0, -2,
       0,  0,  2,
       2,  0,  0,
       0,  2,  0,
      -2,  0,  0,
       0, -2,  0,
       1, -1, -1,
       1,  1, -1,
      -1,  1, -1,
      -1, -1, -1,
       1, -1,  1,
       1,  1,  1,
      -1,  1,  1,
      -1, -1,  1 };
  const int conn[] = 
    { 0, 1, 2, 3, 5, 6, 7, 8, 9, 14, 15, 16, 17,
      3, 2, 1, 0, 4, 8, 7, 6, 9, 13, 12, 11, 10 };
  
  return test_read_write_element( coords, 18, conn, conn, 26, 2, 27, MBPYRAMID );
}

bool test_structured_2d( const char* file );
bool test_structured_3d( const char* file );

bool test_structured_points_2d()
{
  const char file[] = 
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET STRUCTURED_POINTS\n"
   "DIMENSIONS 4 4 1\n"
   "ORIGIN 0 0 0\n"
   "SPACING 1 1 1\n";
  bool rval1 = test_structured_2d( file );
  
    // test again w/ old 1.0 ASPECT_RATIO keyword
  const char file2[] = 
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET STRUCTURED_POINTS\n"
   "DIMENSIONS 4 4 1\n"
   "ORIGIN 0 0 0\n"
   "ASPECT_RATIO 1 1 1\n";
  bool rval2 = test_structured_2d( file2 );
  
  return rval1 && rval2;
}

bool test_structured_grid_2d()
{
  char file[4096] =
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET STRUCTURED_GRID\n"
   "DIMENSIONS 4 4 1\n"
   "POINTS 16 double\n";
  int len = strlen(file);
  for (unsigned i = 0; i < 16; ++i)
    len += sprintf(file+len, "%f %f %f\n", grid_3x3[3*i], grid_3x3[3*i+1], grid_3x3[3*i+2]);
   
  return test_structured_2d( file );
}

bool test_rectilinear_grid_2d()
{
  const char file[] =
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET RECTILINEAR_GRID\n"
   "DIMENSIONS 4 4 1\n"
   "X_COORDINATES 4 float 0 1 2 3\n"
   "Y_COORDINATES 4 float 0 1 2 3\n"
   "Z_COORDINATES 1 float 0\n";
   
  return test_structured_2d( file );
}

bool test_structured_points_3d()
{
  const char file[] = 
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET STRUCTURED_POINTS\n"
   "DIMENSIONS 3 3 3\n"
   "ORIGIN 0 0 0\n"
   "SPACING 1 1 1\n";
  return test_structured_3d( file );
}

bool test_structured_grid_3d()
{
  char file[4096] =
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET STRUCTURED_GRID\n"
   "DIMENSIONS 3 3 3\n"
   "POINTS 27 double\n";
   
  int len = strlen(file);
  for (unsigned i = 0; i < 27; ++i)
    len += sprintf(file+len, "%f %f %f\n", grid_2x2x2[3*i], grid_2x2x2[3*i+1], grid_2x2x2[3*i+2]);
   
  return test_structured_3d( file );
}

bool test_rectilinear_grid_3d()
{
  const char file[] =
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET RECTILINEAR_GRID\n"
   "DIMENSIONS 3 3 3\n"
   "X_COORDINATES 3 float 0 1 2\n"
   "Y_COORDINATES 3 float 0 1 2\n"
   "Z_COORDINATES 3 float 0 1 2\n";
   
  return test_structured_3d( file );
}

bool test_scalar_attrib(const char* vtk_type, DataType mb_type, int count);
bool test_vector_attrib(const char* vtk_type, DataType mb_type);
bool test_tensor_attrib(const char* vtk_type, DataType mb_type);

bool test_scalar_attrib_1_bit()
{
  return test_scalar_attrib("bit", MB_TYPE_BIT, 1);
}

bool test_scalar_attrib_1_uchar()
{
  return test_scalar_attrib("unsigned_char", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_char()
{
  return test_scalar_attrib("char", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_ushort()
{
  return test_scalar_attrib("unsigned_short", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_short()
{
  return test_scalar_attrib("short", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_uint()
{
  return test_scalar_attrib("unsigned_int", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_int()
{
  return test_scalar_attrib("int", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_ulong()
{
  return test_scalar_attrib("unsigned_long", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_long()
{
  return test_scalar_attrib("long", MB_TYPE_INTEGER, 1);
}

bool test_scalar_attrib_1_float()
{
  return test_scalar_attrib("float", MB_TYPE_DOUBLE, 1);
}

bool test_scalar_attrib_1_double()
{
  return test_scalar_attrib("double", MB_TYPE_DOUBLE, 1);
}

bool test_scalar_attrib_4_bit()
{
  return test_scalar_attrib("bit", MB_TYPE_BIT, 4);
}

bool test_scalar_attrib_4_uchar()
{
  return test_scalar_attrib("unsigned_char", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_char()
{
  return test_scalar_attrib("char", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_ushort()
{
  return test_scalar_attrib("unsigned_short", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_short()
{
  return test_scalar_attrib("short", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_uint()
{
  return test_scalar_attrib("unsigned_int", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_int()
{
  return test_scalar_attrib("int", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_ulong()
{
  return test_scalar_attrib("unsigned_long", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_long()
{
  return test_scalar_attrib("long", MB_TYPE_INTEGER, 4);
}

bool test_scalar_attrib_4_float()
{
  return test_scalar_attrib("float", MB_TYPE_DOUBLE, 4);
}

bool test_scalar_attrib_4_double()
{
  return test_scalar_attrib("double", MB_TYPE_DOUBLE, 4);
}

bool test_vector_attrib_bit()
{
  return test_vector_attrib("bit", MB_TYPE_BIT);
}

bool test_vector_attrib_uchar()
{
  return test_vector_attrib("unsigned_char", MB_TYPE_INTEGER);
}

bool test_vector_attrib_char()
{
  return test_vector_attrib("char", MB_TYPE_INTEGER);
}

bool test_vector_attrib_ushort()
{
  return test_vector_attrib("unsigned_short", MB_TYPE_INTEGER);
}

bool test_vector_attrib_short()
{
  return test_vector_attrib("short", MB_TYPE_INTEGER);
}

bool test_vector_attrib_uint()
{
  return test_vector_attrib("unsigned_int", MB_TYPE_INTEGER);
}

bool test_vector_attrib_int()
{
  return test_vector_attrib("int", MB_TYPE_INTEGER);
}

bool test_vector_attrib_ulong()
{
  return test_vector_attrib("unsigned_long", MB_TYPE_INTEGER);
}

bool test_vector_attrib_long()
{
  return test_vector_attrib("long", MB_TYPE_INTEGER);
}

bool test_vector_attrib_float()
{
  return test_vector_attrib("float", MB_TYPE_DOUBLE);
}

bool test_vector_attrib_double()
{
  return test_vector_attrib("double", MB_TYPE_DOUBLE);
}


bool test_tensor_attrib_uchar()
{
  return test_tensor_attrib("unsigned_char", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_char()
{
  return test_tensor_attrib("char", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_ushort()
{
  return test_tensor_attrib("unsigned_short", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_short()
{
  return test_tensor_attrib("short", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_uint()
{
  return test_tensor_attrib("unsigned_int", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_int()
{
  return test_tensor_attrib("int", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_ulong()
{
  return test_tensor_attrib("unsigned_long", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_long()
{
  return test_tensor_attrib("long", MB_TYPE_INTEGER);
}

bool test_tensor_attrib_float()
{
  return test_tensor_attrib("float", MB_TYPE_DOUBLE);
}

bool test_tensor_attrib_double()
{
  return test_tensor_attrib("double", MB_TYPE_DOUBLE);
}




bool read_file( Interface* iface, const char* file )
{
  char fname[] = "tmp_file.vtk";
  FILE* fptr = fopen( fname, "w" );
  fputs( file, fptr ); 
  fclose( fptr );
  
  ErrorCode rval = iface->load_mesh( fname );
  remove( fname );
  CHECK(rval);
  return true;
}

bool write_and_read( Interface* iface1, Interface* iface2 )
{
  const char fname[] = "tmp_file.vtk";
  ErrorCode rval1 = iface1->write_mesh( fname );
  ErrorCode rval2 = iface2->load_mesh( fname );
  remove( fname );
  CHECK(rval1);
  CHECK(rval2);
  return true;
}

bool compare_connectivity( EntityType ,
                           const int* conn1,
                           const int* conn2,
                           unsigned len )
{
  for (unsigned i = 0; i < len; ++i)
    if (conn1[i] != conn2[i])
      return false;
  return true;
}

bool match_vertices_and_elements( Interface* iface,
                                  EntityType moab_type,
                                  unsigned num_vert,
                                  unsigned num_elem,
                                  unsigned vert_per_elem,
                                  const double* coords,
                                  const int* connectivity,
                                  EntityHandle* vert_handles,
                                  EntityHandle* elem_handles )
{
  ErrorCode rval;
  
    // get vertices and check count
  Range verts;
  rval = iface->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK(rval);
  CHECK(verts.size() == num_vert);
  
    // get elements and check count
  Range elems;
  rval = iface->get_entities_by_type( 0, moab_type, elems );
  CHECK(rval);
  CHECK(elems.size() == num_elem);
  
    // get vertex coordinates
  std::vector<EntityHandle> vert_array(num_vert);
  std::copy(verts.begin(), verts.end(), vert_array.begin());
  std::vector<double> mb_coords(3*num_vert);
  rval = iface->get_coords( &vert_array[0], num_vert, &mb_coords[0] );
  CHECK(rval);
  
    // compare vertex coordinates to construct map from 
    // EntityHandle to index in input coordinate list
  std::map<EntityHandle,int> vert_map;
  std::vector<bool> seen(num_vert, false);
  for (unsigned i = 0; i < num_vert; ++i) {
    double* vert_coords = &mb_coords[3*i];
    bool found = false;
    for (unsigned j = 0; j < num_vert; ++j) {
      const double* file_coords = &coords[3*j];
      double dsqr = 0;
      for (unsigned k = 0; k < 3; ++k) {
        double diff = vert_coords[k] - file_coords[k];
        dsqr += diff*diff;
      }
      if (dsqr < 1e-6) {
        CHECK(!seen[j]); // duplicate vertex
        seen[j] = found = true;
        vert_map[vert_array[i]] = j;
        vert_handles[j] = vert_array[i];
        break;
      }
    }
    CHECK(found); // not found?
  }
  
    // check element connectivity
  seen.clear(); seen.resize( num_elem, false );
  Range::iterator iter = elems.begin();
  for (unsigned i = 0; i < num_elem; ++i) {
      // get element connectivity
    EntityHandle elem = *iter; ++iter;
    std::vector<EntityHandle> elem_conn;
    rval = iface->get_connectivity( &elem, 1, elem_conn );
    CHECK(rval);
    CHECK( elem_conn.size() == vert_per_elem );
    
      // convert to input vertex ordering
    std::vector<int> elem_conn2(vert_per_elem);
    for (unsigned j = 0; j < vert_per_elem; ++j) {
      std::map<EntityHandle,int>::iterator k = vert_map.find(elem_conn[j]);
      CHECK( k != vert_map.end() );
      elem_conn2[j] = k->second;
    }
    
      // search input list for matching element
    bool found = false;
    for (unsigned j = 0; j < num_elem; ++j) {
      const int* conn_arr = connectivity + j*vert_per_elem;
      if (!seen[j] &&
          compare_connectivity( moab_type, conn_arr, &elem_conn2[0], vert_per_elem))
      {
        seen[j] = found = true;
        elem_handles[j] = elem;
        break;
      }
    }
    CHECK(found);
  }
  
  return true;
}


bool check_elements( Interface* iface,
                     EntityType moab_type, unsigned num_elem, unsigned vert_per_elem,
                     const double* coords, unsigned num_vert, 
                     const int* connectivity )
{
  std::vector<EntityHandle> junk1(num_vert), junk2(num_elem);
  bool rval = match_vertices_and_elements( iface, moab_type, num_vert, num_elem, 
                                         vert_per_elem, coords, connectivity,
                                         &junk1[0], &junk2[0] );
  CHECK(rval);
  return true;
}
      
                     

bool test_read_write_element( const double* coords, unsigned num_verts,
                              const int* vtk_conn, const int* moab_conn,
                              unsigned num_conn,
                              unsigned num_elem, unsigned vtk_type,
                              EntityType moab_type )
      
{
    // construct VTK file
  char file[4096] = 
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET UNSTRUCTURED_GRID\n";
  size_t len = strlen(file);

  len += sprintf(file+len, "POINTS %u double\n", num_verts);
  for (unsigned i = 0; i < num_verts; ++i)
    len += sprintf(file+len, "%f %f %f\n", coords[3*i], coords[3*i+1], coords[3*i+2] );

  len += sprintf(file+len, "CELLS %u %u\n", num_elem, num_conn+num_elem);
  assert( num_conn % num_elem == 0 );
  unsigned conn_len = num_conn / num_elem;
  for (unsigned i = 0; i < num_elem; ++i) {
    len += sprintf(file+len, "%u", conn_len );
    for (unsigned j = 0; j < conn_len; ++j) 
      len += sprintf(file+len, " %u", vtk_conn[conn_len*i+j]);
    len += sprintf(file+len,"\n");
  }
  
  len += sprintf(file+len,"CELL_TYPES %u\n", num_elem);
  for (unsigned i =0; i < num_elem; ++i)
    len += sprintf(file+len, "%u\n", vtk_type);
  
    // read VTK file and check results
  Core instance1, instance2;
  bool bval = read_file( &instance1, file ); CHECK(bval);
  bval = check_elements( &instance1, moab_type, num_elem, conn_len, coords, num_verts, moab_conn );
  CHECK(bval);
  
    // write, re-read, and check results
  bval = write_and_read( &instance1, &instance2 ); CHECK(bval);
  bval = check_elements( &instance2, moab_type, num_elem, conn_len, coords, num_verts, moab_conn );
  CHECK(bval);
  
  return true;
}
  
bool test_structured_2d( const char* file )
{
    // read VTK file and check results
  Core instance;
  bool bval = read_file( &instance, file ); CHECK(bval);
  bval = check_elements( &instance, MBQUAD, 9, 4, grid_3x3, 16, quad_structured_conn );
  CHECK(bval);
  
  return true;
}

bool test_structured_3d( const char* file )
{
    // read VTK file and check results
  Core instance;
  bool bval = read_file( &instance, file ); CHECK(bval);
  bval = check_elements( &instance, MBHEX, 8, 8, grid_2x2x2, 27, hex_structured_conn );
  CHECK(bval);
  
  return true;
}

const char two_quad_mesh[] =
   "# vtk DataFile Version 3.0\n"
   "MOAB Version 1.00\n"
   "ASCII\n"
   "DATASET UNSTRUCTURED_GRID\n"
   "POINTS 6 float\n"
   "-1 0 0\n"
   " 0 0 0\n"
   " 1 0 0\n"
   "-1 1 0\n"
   " 0 1 0\n"
   " 1 1 0\n"
   "CELLS 2 10\n"
   "4 0 1 4 3\n"
   "4 1 2 5 4\n"
   "CELL_TYPES 2\n"
   "9 9\n";
   
const double two_quad_mesh_coords[] = {
  -1, 0, 0,
   0, 0, 0,
   1, 0, 0,
  -1, 1, 0,
   0, 1, 0,
   1, 1, 0 };
const int two_quad_mesh_conn[] = {
  0, 1, 4, 3,
  1, 2, 5, 4 };

const int vertex_values[] = { 9, 3, 8, 2, 0, 6, 
                              4, 1, 4, 1, 0, 3,
                              8, 6, 6, 4, 0, 2,
                              1, 2, 3, 4, 5, 6,
                              6, 5, 4, 3, 2, 1,
                              0, 6, 1, 5, 2, 4,
                              3, 6, 9, 2, 5, 8,
                              1, 3, 5, 7, 1, 3,
                              5, 8, 1, 9, 7, 4 };
const int element_values[] = { 1001, 1002, 1004, 1003, 50, 60, 51, 61,
                               1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 
                               0, 9, 8, 7, 6, 5, 4, 3, 2, 1 };

void write_data( char* file, size_t& len, DataType type, unsigned count, const int* vals )
{
  switch(type) {
    case MB_TYPE_BIT:
      for (unsigned i = 0; i < count; ++i)
        len += sprintf(file+len, "%d\n", abs(vals[i])%2);
      break;
    case MB_TYPE_INTEGER: 
      for (unsigned i = 0; i < count; ++i)
        len += sprintf(file+len, "%d\n", vals[i]);
      break;
    case MB_TYPE_DOUBLE: 
      for (unsigned i = 0; i < count; ++i)
        len += sprintf(file+len, "%f\n", (double)vals[i]);
      break;
    case MB_TYPE_OPAQUE:
      for (unsigned i = 0; i < count; ++i)
        len += sprintf(file+len, "%d\n", abs(vals[i]%256));
      break;
    default:
      assert( false /* VTK files cannot handle this type */ );
  }
}

bool check_tag_values( Interface* iface,
                       DataType tag_type, int tag_length,
                       int num_entities, const EntityHandle* entities,
                       const int* values )
{
  Tag tag;
  ErrorCode rval = iface->tag_get_handle( "data", tag_length, tag_type, tag ); CHECK(rval);
  
  int size, *intptr;
  double* dblptr;
  rval = iface->tag_get_bytes( tag, size ); CHECK(rval);
  std::vector<unsigned char> data( size * num_entities );
  
  switch (tag_type) {
    case MB_TYPE_BIT:
      rval = iface->tag_get_length( tag, size ); CHECK(rval);
      CHECK( tag_length == size );
      for (int i = 0; i < num_entities; ++i) {
        unsigned char val;
        rval = iface->tag_get_data( tag, entities + i, 1, &val ); CHECK(rval);
        for (int j = 0; j < tag_length; ++j) {
          int bitval = !!(val & (1 << j));
          int expval = abs(*values) % 2;
          CHECK( bitval == expval );
          ++values;
        }
      }
      break;
    case MB_TYPE_OPAQUE:
      rval = iface->tag_get_data( tag, entities, num_entities, &data[0] ); CHECK(rval);
      CHECK( tag_length == size );
      for (int i = 0; i < num_entities; ++i) 
        for (int j = 0; j < tag_length; ++j, ++values) 
          CHECK( (unsigned)(*values % 256) == data[i*tag_length+j] );
      break;
    case MB_TYPE_INTEGER:
      rval = iface->tag_get_data( tag, entities, num_entities, &data[0] ); CHECK(rval);
      CHECK( tag_length*sizeof(int) == (unsigned)size );
      intptr = reinterpret_cast<int*>(&data[0]);
      for (int i = 0; i < num_entities; ++i) 
        for (int j = 0; j < tag_length; ++j, ++values) 
          CHECK( *values == intptr[i*tag_length+j] );
      break;
    case MB_TYPE_DOUBLE:
      rval = iface->tag_get_data( tag, entities, num_entities, &data[0] ); CHECK(rval);
      CHECK( tag_length*sizeof(double) == (unsigned)size );
      dblptr = reinterpret_cast<double*>(&data[0]);
      for (int i = 0; i < num_entities; ++i) 
        for (int j = 0; j < tag_length; ++j, ++values) 
          CHECK( *values == dblptr[i*tag_length+j] );
      break;
    default:
      assert(false);
      return false;
  }
  return true;
}


bool check_tag_values( Interface* iface, DataType type, int vals_per_ent )
{
  EntityHandle vert_handles[6], elem_handles[2];
  bool rval = match_vertices_and_elements( iface, MBQUAD, 6, 2, 4, 
                           two_quad_mesh_coords, two_quad_mesh_conn,
                           vert_handles, elem_handles ); CHECK(rval);
  
  rval = check_tag_values( iface, type, vals_per_ent, 6, vert_handles, vertex_values );
  CHECK(rval);
  rval = check_tag_values( iface, type, vals_per_ent, 2, elem_handles, element_values );
  CHECK(rval);
  return rval;
}

bool check_tag_data( const char* file, DataType type, int vals_per_ent )
{
  bool bval;
  Core instance1, instance2;
  
  bval = read_file( &instance1, file ); CHECK(bval);
  bval = check_tag_values( &instance1, type, vals_per_ent ); CHECK(bval);
  bval = write_and_read( &instance1, &instance2 ); CHECK(bval);
  bval = check_tag_values( &instance2, type, vals_per_ent ); CHECK(bval);
  return true;
}

bool test_scalar_attrib(const char* vtk_type, DataType mb_type, int count)
{
  char file[4096];
  strcpy( file, two_quad_mesh );
  size_t len = strlen(file);
  len += sprintf(file+len, "POINT_DATA 6\n");
  len += sprintf(file+len, "SCALARS data %s %d\n", vtk_type, count );
  len += sprintf(file+len, "LOOKUP_TABLE default\n");
  write_data( file, len, mb_type, 6*count, vertex_values );
  len += sprintf(file+len, "CELL_DATA 2\n");
  len += sprintf(file+len, "SCALARS data %s %d\n", vtk_type, count );
  len += sprintf(file+len, "LOOKUP_TABLE default\n");
  write_data( file, len, mb_type, 2*count, element_values );
  
  return check_tag_data( file, mb_type, count );
} 

bool test_vector_attrib( const char* vtk_type, DataType mb_type )
{
  char file[4096];
  strcpy( file, two_quad_mesh );
  size_t len = strlen(file);
  len += sprintf(file+len, "POINT_DATA 6\n");
  len += sprintf(file+len, "VECTORS data %s\n", vtk_type );
  write_data( file, len, mb_type, 6*3, vertex_values );
  len += sprintf(file+len, "CELL_DATA 2\n");
  len += sprintf(file+len, "VECTORS data %s\n", vtk_type );
  write_data( file, len, mb_type, 2*3, element_values );
  
  return check_tag_data( file, mb_type, 3 );
} 

bool test_tensor_attrib( const char* vtk_type, DataType mb_type )
{
  char file[4096];
  strcpy( file, two_quad_mesh );
  size_t len = strlen(file);
  len += sprintf(file+len, "POINT_DATA 6\n");
  len += sprintf(file+len, "TENSORS data %s\n", vtk_type );
  write_data( file, len, mb_type, 6*9, vertex_values );
  len += sprintf(file+len, "CELL_DATA 2\n");
  len += sprintf(file+len, "TENSORS data %s\n", vtk_type );
  write_data( file, len, mb_type, 2*9, element_values );
  
  return check_tag_data( file, mb_type, 9 );
} 

bool test_subset( )
{
  Core moab_inst;
  Interface& moab = moab_inst;
  ErrorCode rval;
  
    // create 9 nodes in grid pattern
  EntityHandle verts[9];
  const double coords[][3] = { { 0, 0, 0 },
                               { 1, 0, 0 },
                               { 2, 0, 0 },
                               { 0, 1, 0 },
                               { 1, 1, 0 },
                               { 2, 1, 0 },
                               { 0, 2, 0 },
                               { 1, 2, 0 },
                               { 2, 2, 0 } };
  for (unsigned i = 0; i < 9; ++i) {
    rval = moab.create_vertex(coords[i], verts[i]);
    assert(MB_SUCCESS == rval);
  }
  
    // create 4 quad elements in grid pattern
  const int conn[][4] = { { 0, 1, 4, 3 },
                          { 1, 2, 5, 4 },
                          { 3, 4, 7, 6 },
                          { 4, 5, 8, 7 } };
  EntityHandle econn[4], elems[4];
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 4; ++j)
      econn[j] = verts[conn[i][j]];
    rval = moab.create_element( MBQUAD, econn, 4, elems[i] );
    assert(MB_SUCCESS == rval);
  }
  
    // create 3 meshsets
  EntityHandle sets[3];
  for (unsigned i = 0;i < 3; ++i) {
    rval = moab.create_meshset( 0, sets[i] );
    assert(MB_SUCCESS == rval);
  }
  
    // add element 3 to set 0
  rval = moab.add_entities( sets[0], elems+3, 1 );
  assert(MB_SUCCESS == rval);
    // add node 2 to set 1
  rval = moab.add_entities( sets[1], verts+2, 1 );
  assert(MB_SUCCESS == rval);
    // add element 2 and 3 to set 2
  rval = moab.add_entities( sets[2], elems+2, 2 );
  assert(MB_SUCCESS == rval);
  
    // make set 2 a child of set 1
  rval = moab.add_child_meshset( sets[1], sets[2] );
  assert(MB_SUCCESS == rval);
    // put set 1 in set 0
  rval = moab.add_entities( sets[0], sets+1, 1 );
  assert(MB_SUCCESS == rval);
  
    // write sets[0] to vtk file
  rval = moab.write_mesh(  "tmp_file.vtk", sets, 1 );
  CHECK(rval);
   
    // read data back in
  moab.delete_mesh();
  rval = moab.load_mesh( "tmp_file.vtk" );
  remove( "tmp_file.vtk" );
  CHECK(rval);
  
    // writer should have written all three sets,
    // so the resulting mesh should be elems[2], elems[3], 
    // and verts[2]
  Range new_elems, new_verts;
  rval = moab.get_entities_by_type( 0, MBQUAD, new_elems );
  CHECK(rval);
  CHECK( new_elems.size() == 2 );
  rval = moab.get_entities_by_type( 0, MBVERTEX, new_verts );
  CHECK(rval);
  CHECK( new_verts.size() == 7 );
  
    // vertex not in element closure should have coords of 2,0,0
  Range elem_verts;
  rval = moab.get_adjacencies( new_elems, 0, false, elem_verts, Interface::UNION );
  CHECK(rval);
  CHECK(elem_verts.size() == 6);
  Range free_verts( subtract( new_verts, elem_verts ) );
  CHECK(free_verts.size() == 1 );
  double vcoords[3];
  rval = moab.get_coords( free_verts, vcoords );
  CHECK( vcoords[0] == 2 );
  CHECK( vcoords[1] == 0 );
  CHECK( vcoords[2] == 0 );
  
  return true;
}

// Test technically invalid but somewhat common insertion of
// FIELD blocks within an UNSTRUCTURED_GRID dataset
bool test_unstructured_field()
{
    // Use existing file defined in 'two_quad_mesh', but
    // insert a few field data blocks
  std::istringstream base_data( two_quad_mesh );
  std::ostringstream file_data;
  std::string line;
  while (getline( base_data, line )) {
    if (0 == line.find("POINTS")) {
      file_data << "FIELD FieldData 2" << std::endl
                << "avtOriginalBounds 1 6 float" << std::endl
                << "-10 10 -10 10 -10 10 " << std::endl
                << "TIME 1 1 double" << std::endl
                << "10.543" << std::endl;
    }
    else if (0 == line.find("CELLS")) {
      file_data << "FIELD more_data 2" << std::endl
                << "first_array 3 2 int" << std::endl
                << "0 1 2" << std::endl
                << "3 4 5" << std::endl
                << "second_array 4 3 bit" << std::endl
                << "0 0 0 0" << std::endl
                << "1 1 1 1" << std::endl
                << "1 0 1 0" << std::endl;
    }
    file_data << line << std::endl;
  }
  
  Core core;
  Interface& mb = core;
  bool rval = read_file(&mb, file_data.str().c_str());
  CHECK(rval);

  EntityHandle vert_handles[6], elem_handles[2];
  rval = match_vertices_and_elements( &mb, MBQUAD, 6, 2, 4, 
                       two_quad_mesh_coords, two_quad_mesh_conn,
                       vert_handles, elem_handles ); 
  CHECK(rval);
  
  return true;
}
