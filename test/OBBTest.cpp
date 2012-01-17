#include "OrientedBox.hpp"
#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"

using namespace moab;

#include <assert.h>

const double TOL = 1e-6;
int error_count = 0;

static void test_basic();         // test basic properties (volume(), etc.)
static void test_contained();     
static void test_ray_intersect(); 
static void test_closest_point();
static void test_build_from_tri();
static void test_build_from_pts();
static void test_save();

#include <iostream>
#define ASSERT_VECTOR_ELEMENT( A, B ) assert_vector_element( (A), (B), #A, #B, __LINE__ )
#define ASSERT_VECTORS_EQUAL(A, B) assert_vectors_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT_DOUBLES_EQUAL(A, B) assert_doubles_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT(B) assert_bool( (B), #B, __LINE__ )

static void assert_vector_element( const CartVect& a, const CartVect* b, const char* sa, const char* sb, int lineno );
static void assert_vectors_equal( const CartVect& a, const CartVect& b, const char* sa, const char* sb, int lineno );
static void assert_doubles_equal( double a, double b, const char* sa, const char* sb, int lineno );
static void assert_bool( bool b, const char* sb, int lineno );

int main()
{
  test_basic();
  test_contained();
  test_ray_intersect();
  test_closest_point();
  test_build_from_tri();
  test_build_from_pts();
  test_save();
  
  return error_count;
}

/********************* Declare some boxes to test ***************************/

  // define unit box centered at origin
const CartVect origin( 0.0, 0.0, 0.0 );
const CartVect unitaxes[3] = { CartVect(0.5, 0.0, 0.0),
                                 CartVect(0.0, 0.5, 0.0),
                                 CartVect(0.0, 0.0, 0.5) };
const OrientedBox unitbox( unitaxes, origin );

  // define axis-aligned unit box outside origin
const CartVect unitcenter( 10, 20, 30 );
const OrientedBox offsetbox( unitaxes, unitcenter );

  // define non-unit centered at origin
const CartVect origaxes[3] = { 5*unitaxes[0],
                                10*unitaxes[1],
                                .1*unitaxes[2] };
const OrientedBox oblongbox( origaxes, origin );

  // define non-axis-aligned unit box at origin
const CartVect rotaxes[3] = { unit( CartVect( 1.0, 1.0, 0.0 ) ),
                                unit( CartVect( 1.0,-1.0, 1.0 ) ),
                                unit( CartVect( 1.0, 1.0, 0.0 )*CartVect( 1.0,-1.0, 1.0 ) ) };
const OrientedBox rotbox( rotaxes, origin );

/********************* Utility methods for tests ***************************/

// return point at specified fraction between box center and specified box corner
static CartVect scaled_corner( const OrientedBox& box, int corner, double factor )
{
  static const int signs[][3] = { { 1, 1,-1},
                                  {-1, 1,-1},
                                  {-1,-1,-1},
                                  { 1,-1,-1},
                                  { 1, 1, 1},
                                  {-1, 1, 1},
                                  {-1,-1, 1},
                                  { 1,-1, 1} };
  return box.center 
       + signs[corner][0]*factor*box.scaled_axis(0)
       + signs[corner][1]*factor*box.scaled_axis(1)
       + signs[corner][2]*factor*box.scaled_axis(2);
}

// return point at specified fraction between box center and specified box face
static CartVect scaled_face( const OrientedBox& box, int face, double factor )
{
  assert(face >= 0 && face <= 6);
  int sign = face % 2 ? -1 : 1;
  return box.center + factor * sign * box.scaled_axis(face/2);
}

// get vector containing axis lengths, ordered from smallest to largest
static void axis_dims( const CartVect axis[3], CartVect& dims )
{
  dims = CartVect(axis[0].length(), axis[1].length(), axis[2].length());
  if (dims[0] > dims[1]) 
    std::swap(dims[0], dims[1]);
  if (dims[1] > dims[2])
    std::swap(dims[1], dims[2]);
  if (dims[0] > dims[1]) 
    std::swap(dims[0], dims[1]);
} 
  

/********************* The Actual Tests ***************************/

static void test_basic()
{
  CartVect dims;
  
  axis_dims( unitaxes, dims );
  ASSERT_VECTORS_EQUAL( unitbox.center, origin );
  ASSERT_VECTOR_ELEMENT( unitbox.scaled_axis(0), unitaxes );
  ASSERT_VECTOR_ELEMENT( unitbox.scaled_axis(1), unitaxes );
  ASSERT_VECTOR_ELEMENT( unitbox.scaled_axis(2), unitaxes );
  ASSERT_DOUBLES_EQUAL( unitbox.inner_radius(), dims[0] );
  ASSERT_DOUBLES_EQUAL( unitbox.outer_radius(), dims.length() );
  ASSERT_DOUBLES_EQUAL( unitbox.volume(), 8.0*dims[0]*dims[1]*dims[2] );
  ASSERT_VECTORS_EQUAL( unitbox.dimensions(), 2*dims );
  
  axis_dims( unitaxes, dims );
  ASSERT_VECTORS_EQUAL( offsetbox.center, unitcenter );
  ASSERT_VECTOR_ELEMENT( offsetbox.scaled_axis(0), unitaxes );
  ASSERT_VECTOR_ELEMENT( offsetbox.scaled_axis(1), unitaxes );
  ASSERT_VECTOR_ELEMENT( offsetbox.scaled_axis(2), unitaxes );
  ASSERT_DOUBLES_EQUAL( offsetbox.inner_radius(), dims[0] );
  ASSERT_DOUBLES_EQUAL( offsetbox.outer_radius(), dims.length() );
  ASSERT_DOUBLES_EQUAL( offsetbox.volume(), 8.0*dims[0]*dims[1]*dims[2] );
  ASSERT_VECTORS_EQUAL( offsetbox.dimensions(), 2*dims );
  
  axis_dims( origaxes, dims );
  ASSERT_VECTORS_EQUAL( oblongbox.center, origin );
  ASSERT_VECTOR_ELEMENT( oblongbox.scaled_axis(0), origaxes );
  ASSERT_VECTOR_ELEMENT( oblongbox.scaled_axis(1), origaxes );
  ASSERT_VECTOR_ELEMENT( oblongbox.scaled_axis(2), origaxes );
  ASSERT_DOUBLES_EQUAL( oblongbox.inner_radius(), dims[0] );
  ASSERT_DOUBLES_EQUAL( oblongbox.outer_radius(), dims.length() );
  ASSERT_DOUBLES_EQUAL( oblongbox.volume(), 8.0*dims[0]*dims[1]*dims[2] );
  ASSERT_VECTORS_EQUAL( oblongbox.dimensions(), 2*dims );
  
  axis_dims( rotaxes, dims );
  ASSERT_VECTORS_EQUAL( rotbox.center, origin );
  ASSERT_VECTOR_ELEMENT( rotbox.scaled_axis(0), rotaxes );
  ASSERT_VECTOR_ELEMENT( rotbox.scaled_axis(1), rotaxes );
  ASSERT_VECTOR_ELEMENT( rotbox.scaled_axis(2), rotaxes );
  ASSERT_DOUBLES_EQUAL( rotbox.inner_radius(), dims[0] );
  ASSERT_DOUBLES_EQUAL( rotbox.outer_radius(), dims.length() );
  ASSERT_DOUBLES_EQUAL( rotbox.volume(), 8.0*dims[0]*dims[1]*dims[2] );
  ASSERT_VECTORS_EQUAL( rotbox.dimensions(), 2*dims );
}
  

static void test_contained() 
{
    // first do tests of unit box
    
    // test points inside box
  ASSERT( unitbox.contained( unitbox.center, TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 0, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 1, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 2, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 3, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 4, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 5, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 6, 0.6 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 7, 0.6 ), TOL ) );
    
    // test points at box corners
  ASSERT( unitbox.contained( scaled_corner( unitbox, 0, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 1, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 2, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 3, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 4, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 5, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 6, 1.0 ), TOL ) );
  ASSERT( unitbox.contained( scaled_corner( unitbox, 7, 1.0 ), TOL ) );
    
    // test points at center of each face
  ASSERT( unitbox.contained( scaled_face( unitbox, 0, 1.0), TOL ) );
  ASSERT( unitbox.contained( scaled_face( unitbox, 1, 1.0), TOL ) );
  ASSERT( unitbox.contained( scaled_face( unitbox, 2, 1.0), TOL ) );
  ASSERT( unitbox.contained( scaled_face( unitbox, 3, 1.0), TOL ) );
  ASSERT( unitbox.contained( scaled_face( unitbox, 4, 1.0), TOL ) );
  ASSERT( unitbox.contained( scaled_face( unitbox, 5, 1.0), TOL ) );

    // test points beyond each corner
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 0, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 1, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 2, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 3, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 4, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 5, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 6, 1.2 ), TOL ) );
  ASSERT( ! unitbox.contained( scaled_corner( unitbox, 7, 1.2 ), TOL ) );
   
    // test points beyond the center of each face
  ASSERT( ! unitbox.contained( scaled_face( unitbox, 0, 1.3), TOL ) );
  ASSERT( ! unitbox.contained( scaled_face( unitbox, 1, 1.3), TOL ) );
  ASSERT( ! unitbox.contained( scaled_face( unitbox, 2, 1.3), TOL ) );
  ASSERT( ! unitbox.contained( scaled_face( unitbox, 3, 1.3), TOL ) );
  ASSERT( ! unitbox.contained( scaled_face( unitbox, 4, 1.3), TOL ) );
  ASSERT( ! unitbox.contained( scaled_face( unitbox, 5, 1.3), TOL ) );

    // now do offset box 
    
    // test points inside box
  ASSERT( offsetbox.contained( offsetbox.center, TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 0, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 1, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 2, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 3, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 4, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 5, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 6, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 7, 0.6 ), TOL ) );
    
    // test points at box corners
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 0, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 1, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 2, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 3, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 4, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 5, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 6, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 7, 1.0 ), TOL ) );
    
    // test points at center of each face
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 0, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 1, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 2, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 3, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 4, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 5, 1.0), TOL ) );

    // test points beyond each corner
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 0, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 1, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 2, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 3, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 4, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 5, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 6, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 7, 1.2 ), TOL ) );
   
    // test points beyond the center of each face
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 0, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 1, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 2, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 3, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 4, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 5, 1.3), TOL ) );

    // now do oblong box 
    
    // test points inside box
  ASSERT( oblongbox.contained( oblongbox.center, TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 0, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 1, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 2, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 3, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 4, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 5, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 6, 0.6 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 7, 0.6 ), TOL ) );
    
    // test points at box corners
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 0, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 1, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 2, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 3, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 4, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 5, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 6, 1.0 ), TOL ) );
  ASSERT( oblongbox.contained( scaled_corner( oblongbox, 7, 1.0 ), TOL ) );
    
    // test points at center of each face
  ASSERT( oblongbox.contained( scaled_face( oblongbox, 0, 1.0), TOL ) );
  ASSERT( oblongbox.contained( scaled_face( oblongbox, 1, 1.0), TOL ) );
  ASSERT( oblongbox.contained( scaled_face( oblongbox, 2, 1.0), TOL ) );
  ASSERT( oblongbox.contained( scaled_face( oblongbox, 3, 1.0), TOL ) );
  ASSERT( oblongbox.contained( scaled_face( oblongbox, 4, 1.0), TOL ) );
  ASSERT( oblongbox.contained( scaled_face( oblongbox, 5, 1.0), TOL ) );

    // test points beyond each corner
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 0, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 1, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 2, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 3, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 4, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 5, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 6, 1.2 ), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_corner( oblongbox, 7, 1.2 ), TOL ) );
   
    // test points beyond the center of each face
  ASSERT( ! oblongbox.contained( scaled_face( oblongbox, 0, 1.3), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_face( oblongbox, 1, 1.3), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_face( oblongbox, 2, 1.3), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_face( oblongbox, 3, 1.3), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_face( oblongbox, 4, 1.3), TOL ) );
  ASSERT( ! oblongbox.contained( scaled_face( oblongbox, 5, 1.3), TOL ) );

    // now do offset box 
    
    // test points inside box
  ASSERT( offsetbox.contained( offsetbox.center, TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 0, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 1, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 2, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 3, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 4, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 5, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 6, 0.6 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 7, 0.6 ), TOL ) );
    
    // test points at box corners
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 0, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 1, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 2, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 3, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 4, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 5, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 6, 1.0 ), TOL ) );
  ASSERT( offsetbox.contained( scaled_corner( offsetbox, 7, 1.0 ), TOL ) );
    
    // test points at center of each face
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 0, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 1, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 2, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 3, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 4, 1.0), TOL ) );
  ASSERT( offsetbox.contained( scaled_face( offsetbox, 5, 1.0), TOL ) );

    // test points beyond each corner
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 0, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 1, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 2, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 3, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 4, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 5, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 6, 1.2 ), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_corner( offsetbox, 7, 1.2 ), TOL ) );
   
    // test points beyond the center of each face
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 0, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 1, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 2, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 3, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 4, 1.3), TOL ) );
  ASSERT( ! offsetbox.contained( scaled_face( offsetbox, 5, 1.3), TOL ) );

    // now do rotated box 
    
    // test points inside box
  ASSERT( rotbox.contained( rotbox.center, TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 0, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 1, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 2, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 3, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 4, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 5, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 6, 0.6 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 7, 0.6 ), TOL ) );
    
    // test points at box corners
  ASSERT( rotbox.contained( scaled_corner( rotbox, 0, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 1, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 2, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 3, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 4, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 5, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 6, 1.0 ), TOL ) );
  ASSERT( rotbox.contained( scaled_corner( rotbox, 7, 1.0 ), TOL ) );
    
    // test points at center of each face
  ASSERT( rotbox.contained( scaled_face( rotbox, 0, 1.0), TOL ) );
  ASSERT( rotbox.contained( scaled_face( rotbox, 1, 1.0), TOL ) );
  ASSERT( rotbox.contained( scaled_face( rotbox, 2, 1.0), TOL ) );
  ASSERT( rotbox.contained( scaled_face( rotbox, 3, 1.0), TOL ) );
  ASSERT( rotbox.contained( scaled_face( rotbox, 4, 1.0), TOL ) );
  ASSERT( rotbox.contained( scaled_face( rotbox, 5, 1.0), TOL ) );

    // test points beyond each corner
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 0, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 1, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 2, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 3, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 4, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 5, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 6, 1.2 ), TOL ) );
  ASSERT( ! rotbox.contained( scaled_corner( rotbox, 7, 1.2 ), TOL ) );
   
    // test points beyond the center of each face
  ASSERT( ! rotbox.contained( scaled_face( rotbox, 0, 1.3), TOL ) );
  ASSERT( ! rotbox.contained( scaled_face( rotbox, 1, 1.3), TOL ) );
  ASSERT( ! rotbox.contained( scaled_face( rotbox, 2, 1.3), TOL ) );
  ASSERT( ! rotbox.contained( scaled_face( rotbox, 3, 1.3), TOL ) );
  ASSERT( ! rotbox.contained( scaled_face( rotbox, 4, 1.3), TOL ) );
  ASSERT( ! rotbox.contained( scaled_face( rotbox, 5, 1.3), TOL ) );
}

static void test_closest_point()
{
  CartVect result;
  
    // start with unit box
    
    // test locations inside box
  unitbox.closest_location_in_box( unitbox.center, result );
  ASSERT_VECTORS_EQUAL( result, unitbox.center );
  double f = 0.5;
  unitbox.closest_location_in_box( scaled_corner(unitbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,0,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,1,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,2,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,3,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,4,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,5,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,6,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,7,f) );
    
    // test each corner
  f = 1.0;
  unitbox.closest_location_in_box( scaled_corner(unitbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,0,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,1,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,2,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,3,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,4,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,5,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,6,f) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,7,f) );

    // test outside each corner
  f = 1.5;
  unitbox.closest_location_in_box( scaled_corner(unitbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,0,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,1,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,2,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,3,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,4,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,5,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,6,1) );
  unitbox.closest_location_in_box( scaled_corner(unitbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(unitbox,7,1) );
  
    // test on each face
  f = 1.0;
  unitbox.closest_location_in_box( scaled_face(unitbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,0,f) );
  unitbox.closest_location_in_box( scaled_face(unitbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,1,f) );
  unitbox.closest_location_in_box( scaled_face(unitbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,2,f) );
  unitbox.closest_location_in_box( scaled_face(unitbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,3,f) );
  unitbox.closest_location_in_box( scaled_face(unitbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,4,f) );
  unitbox.closest_location_in_box( scaled_face(unitbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,5,f) );
  
    // test outside each face
  f = 1.5;
  unitbox.closest_location_in_box( scaled_face(unitbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,0,1) );
  unitbox.closest_location_in_box( scaled_face(unitbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,1,1) );
  unitbox.closest_location_in_box( scaled_face(unitbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,2,1) );
  unitbox.closest_location_in_box( scaled_face(unitbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,3,1) );
  unitbox.closest_location_in_box( scaled_face(unitbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,4,1) );
  unitbox.closest_location_in_box( scaled_face(unitbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(unitbox,5,1) );
  
    // next offset box
    
    // test locations inside box
  offsetbox.closest_location_in_box( offsetbox.center, result );
  ASSERT_VECTORS_EQUAL( result, offsetbox.center );
  f = 0.5;
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,0,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,1,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,2,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,3,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,4,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,5,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,6,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,7,f) );
    
    // test each corner
  f = 1.0;
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,0,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,1,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,2,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,3,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,4,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,5,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,6,f) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,7,f) );

    // test outside each corner
  f = 1.5;
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,0,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,1,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,2,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,3,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,4,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,5,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,6,1) );
  offsetbox.closest_location_in_box( scaled_corner(offsetbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(offsetbox,7,1) );
  
    // test on each face
  f = 1.0;
  offsetbox.closest_location_in_box( scaled_face(offsetbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,0,f) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,1,f) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,2,f) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,3,f) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,4,f) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,5,f) );
  
    // test outside each face
  f = 1.5;
  offsetbox.closest_location_in_box( scaled_face(offsetbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,0,1) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,1,1) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,2,1) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,3,1) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,4,1) );
  offsetbox.closest_location_in_box( scaled_face(offsetbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(offsetbox,5,1) );
  
    // next oblong box
    
    // test locations inside box
  oblongbox.closest_location_in_box( oblongbox.center, result );
  ASSERT_VECTORS_EQUAL( result, oblongbox.center );
  f = 0.5;
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,0,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,1,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,2,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,3,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,4,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,5,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,6,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,7,f) );
    
    // test each corner
  f = 1.0;
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,0,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,1,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,2,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,3,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,4,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,5,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,6,f) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,7,f) );

    // test outside each corner
  f = 1.5;
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,0,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,1,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,2,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,3,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,4,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,5,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,6,1) );
  oblongbox.closest_location_in_box( scaled_corner(oblongbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(oblongbox,7,1) );
  
    // test on each face
  f = 1.0;
  oblongbox.closest_location_in_box( scaled_face(oblongbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,0,f) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,1,f) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,2,f) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,3,f) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,4,f) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,5,f) );
  
    // test outside each face
  f = 1.5;
  oblongbox.closest_location_in_box( scaled_face(oblongbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,0,1) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,1,1) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,2,1) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,3,1) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,4,1) );
  oblongbox.closest_location_in_box( scaled_face(oblongbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(oblongbox,5,1) );
  
    // next rotated box
    
    // test locations inside box
  rotbox.closest_location_in_box( rotbox.center, result );
  ASSERT_VECTORS_EQUAL( result, rotbox.center );
  f = 0.5;
  rotbox.closest_location_in_box( scaled_corner(rotbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,0,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,1,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,2,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,3,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,4,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,5,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,6,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,7,f) );
    
    // test each corner
  f = 1.0;
  rotbox.closest_location_in_box( scaled_corner(rotbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,0,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,1,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,2,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,3,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,4,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,5,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,6,f) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,7,f) );

    // test outside each corner
  f = 1.5;
  rotbox.closest_location_in_box( scaled_corner(rotbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,0,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,1,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,2,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,3,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,4,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,5,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,6,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,6,1) );
  rotbox.closest_location_in_box( scaled_corner(rotbox,7,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_corner(rotbox,7,1) );
  
    // test on each face
  f = 1.0;
  rotbox.closest_location_in_box( scaled_face(rotbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,0,f) );
  rotbox.closest_location_in_box( scaled_face(rotbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,1,f) );
  rotbox.closest_location_in_box( scaled_face(rotbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,2,f) );
  rotbox.closest_location_in_box( scaled_face(rotbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,3,f) );
  rotbox.closest_location_in_box( scaled_face(rotbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,4,f) );
  rotbox.closest_location_in_box( scaled_face(rotbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,5,f) );
  
    // test outside each face
  f = 1.5;
  rotbox.closest_location_in_box( scaled_face(rotbox,0,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,0,1) );
  rotbox.closest_location_in_box( scaled_face(rotbox,1,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,1,1) );
  rotbox.closest_location_in_box( scaled_face(rotbox,2,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,2,1) );
  rotbox.closest_location_in_box( scaled_face(rotbox,3,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,3,1) );
  rotbox.closest_location_in_box( scaled_face(rotbox,4,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,4,1) );
  rotbox.closest_location_in_box( scaled_face(rotbox,5,f), result );
  ASSERT_VECTORS_EQUAL( result, scaled_face(rotbox,5,1) );
}

void test_ray_intersect()
{
  CartVect dir, pt;
  
  // start with unit box
  
    // test ray from box center towards each face
  dir = scaled_face( unitbox, 0, 1.0 ) - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( unitbox.center, dir, TOL ) );
  dir = scaled_face( unitbox, 1, 1.0 ) - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( unitbox.center, dir, TOL ) );
  dir = scaled_face( unitbox, 2, 1.0 ) - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( unitbox.center, dir, TOL ) );
  dir = scaled_face( unitbox, 3, 1.0 ) - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( unitbox.center, dir, TOL ) );
  dir = scaled_face( unitbox, 4, 1.0 ) - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( unitbox.center, dir, TOL ) );
  dir = scaled_face( unitbox, 5, 1.0 ) - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( unitbox.center, dir, TOL ) );
  
    // test ray starting outside each face and pointing towards box center
  pt = scaled_face( unitbox, 0, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 1, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 2, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 3, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 4, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 5, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray from outside first face toward box center, with nonnegative ray length
  pt = scaled_face( unitbox, 0, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  const double short_pos =  0.99;
  const double short_neg = -0.99;
  const double long_pos  =  1.01;
  const double long_neg  = -1.01;
  ASSERT(  unitbox.intersect_ray( pt, dir, TOL, &long_pos ) );
  ASSERT( !unitbox.intersect_ray( pt, dir, TOL, &short_pos ) );

    // test ray from outside first face away from box center, with negative ray length
  ASSERT(  unitbox.intersect_ray( pt, -dir, TOL, NULL, &long_neg ) );
  ASSERT( !unitbox.intersect_ray( pt, -dir, TOL, NULL, &short_neg ) );

    // test ray from outside first face toward box center, with both ray lengths
    // Note that box.intersect_ray requires neg_ray_len<0 and -neg_ray_len<=nonneg_ray_len
  ASSERT(  unitbox.intersect_ray( pt, dir, TOL, &long_pos,  &long_neg ) );
  ASSERT( !unitbox.intersect_ray( pt, dir, TOL, &short_pos, &long_neg ) );
  ASSERT(  unitbox.intersect_ray( pt, dir, TOL, &long_pos,  &short_neg ) );
  ASSERT( !unitbox.intersect_ray( pt, dir, TOL, &short_pos, &short_neg ) );

    // test ray starting inside box and passing through a corner
  pt = scaled_corner( unitbox, 0, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 1, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 2, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 3, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 4, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 5, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 6, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 7, 0.3 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside box and passing through opposite corners
  pt = scaled_corner( unitbox, 0, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 1, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 2, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 3, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 4, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 5, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 6, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( unitbox, 7, 3.0 ); dir = unitbox.center - pt; dir.normalize();
  ASSERT( unitbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and pointing away from box
  pt = scaled_face( unitbox, 0, 3.0 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 1, 3.0 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 2, 3.0 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 3, 3.0 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 4, 3.0 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 5, 3.0 ); dir = pt - unitbox.center; dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and parallel to face
  pt = scaled_face( unitbox, 0, 3.0 ); dir = unitbox.scaled_axis(1); dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 1, 3.0 ); dir = unitbox.scaled_axis(2); dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 2, 3.0 ); dir = unitbox.scaled_axis(0); dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 3, 3.0 ); dir = unitbox.scaled_axis(2); dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 4, 3.0 ); dir = unitbox.scaled_axis(0); dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( unitbox, 5, 3.0 ); dir = unitbox.scaled_axis(1); dir.normalize();
  ASSERT( ! unitbox.intersect_ray( pt, dir, TOL ) );
  
  // next do offset box
  
    // test ray from box center towards each face
  dir = scaled_face( offsetbox, 0, 1.0 ) - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( offsetbox.center, dir, TOL ) );
  dir = scaled_face( offsetbox, 1, 1.0 ) - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( offsetbox.center, dir, TOL ) );
  dir = scaled_face( offsetbox, 2, 1.0 ) - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( offsetbox.center, dir, TOL ) );
  dir = scaled_face( offsetbox, 3, 1.0 ) - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( offsetbox.center, dir, TOL ) );
  dir = scaled_face( offsetbox, 4, 1.0 ) - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( offsetbox.center, dir, TOL ) );
  dir = scaled_face( offsetbox, 5, 1.0 ) - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( offsetbox.center, dir, TOL ) );
  
    // test ray starting outside each face and pointing towards box center
  pt = scaled_face( offsetbox, 0, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 1, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 2, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 3, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 4, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 5, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting inside box and passing through a corner
  pt = scaled_corner( offsetbox, 0, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 1, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 2, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 3, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 4, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 5, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 6, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 7, 0.3 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside box and passing through opposite corners
  pt = scaled_corner( offsetbox, 0, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 1, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 2, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 3, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 4, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 5, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 6, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( offsetbox, 7, 3.0 ); dir = offsetbox.center - pt; dir.normalize();
  ASSERT( offsetbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and pointing away from box
  pt = scaled_face( offsetbox, 0, 3.0 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 1, 3.0 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 2, 3.0 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 3, 3.0 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 4, 3.0 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 5, 3.0 ); dir = pt - offsetbox.center; dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and parallel to face
  pt = scaled_face( offsetbox, 0, 3.0 ); dir = offsetbox.scaled_axis(1); dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 1, 3.0 ); dir = offsetbox.scaled_axis(2); dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 2, 3.0 ); dir = offsetbox.scaled_axis(0); dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 3, 3.0 ); dir = offsetbox.scaled_axis(2); dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 4, 3.0 ); dir = offsetbox.scaled_axis(0); dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( offsetbox, 5, 3.0 ); dir = offsetbox.scaled_axis(1); dir.normalize();
  ASSERT( ! offsetbox.intersect_ray( pt, dir, TOL ) );
  
  // next do oblong box
  
    // test ray from box center towards each face
  dir = scaled_face( oblongbox, 0, 1.0 ) - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( oblongbox.center, dir, TOL ) );
  dir = scaled_face( oblongbox, 1, 1.0 ) - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( oblongbox.center, dir, TOL ) );
  dir = scaled_face( oblongbox, 2, 1.0 ) - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( oblongbox.center, dir, TOL ) );
  dir = scaled_face( oblongbox, 3, 1.0 ) - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( oblongbox.center, dir, TOL ) );
  dir = scaled_face( oblongbox, 4, 1.0 ) - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( oblongbox.center, dir, TOL ) );
  dir = scaled_face( oblongbox, 5, 1.0 ) - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( oblongbox.center, dir, TOL ) );
  
    // test ray starting outside each face and pointing towards box center
  pt = scaled_face( oblongbox, 0, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 1, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 2, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 3, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 4, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 5, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting inside box and passing through a corner
  pt = scaled_corner( oblongbox, 0, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 1, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 2, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 3, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 4, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 5, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 6, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 7, 0.3 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside box and passing through opposite corners
  pt = scaled_corner( oblongbox, 0, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 1, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 2, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 3, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 4, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 5, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 6, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( oblongbox, 7, 3.0 ); dir = oblongbox.center - pt; dir.normalize();
  ASSERT( oblongbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and pointing away from box
  pt = scaled_face( oblongbox, 0, 3.0 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 1, 3.0 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 2, 3.0 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 3, 3.0 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 4, 3.0 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 5, 3.0 ); dir = pt - oblongbox.center; dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and parallel to face
  pt = scaled_face( oblongbox, 0, 3.0 ); dir = oblongbox.scaled_axis(1); dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 1, 3.0 ); dir = oblongbox.scaled_axis(2); dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 2, 3.0 ); dir = oblongbox.scaled_axis(0); dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 3, 3.0 ); dir = oblongbox.scaled_axis(2); dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 4, 3.0 ); dir = oblongbox.scaled_axis(0); dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( oblongbox, 5, 3.0 ); dir = oblongbox.scaled_axis(1); dir.normalize();
  ASSERT( ! oblongbox.intersect_ray( pt, dir, TOL ) );
  
  // next do rotated box
  
    // test ray from box center towards each face
  dir = scaled_face( rotbox, 0, 1.0 ) - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( rotbox.center, dir, TOL ) );
  dir = scaled_face( rotbox, 1, 1.0 ) - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( rotbox.center, dir, TOL ) );
  dir = scaled_face( rotbox, 2, 1.0 ) - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( rotbox.center, dir, TOL ) );
  dir = scaled_face( rotbox, 3, 1.0 ) - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( rotbox.center, dir, TOL ) );
  dir = scaled_face( rotbox, 4, 1.0 ) - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( rotbox.center, dir, TOL ) );
  dir = scaled_face( rotbox, 5, 1.0 ) - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( rotbox.center, dir, TOL ) );
  
    // test ray starting outside each face and pointing towards box center
  pt = scaled_face( rotbox, 0, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 1, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 2, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 3, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 4, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 5, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting inside box and passing through a corner
  pt = scaled_corner( rotbox, 0, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 1, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 2, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 3, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 4, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 5, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 6, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 7, 0.3 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside box and passing through opposite corners
  pt = scaled_corner( rotbox, 0, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 1, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 2, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 3, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 4, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 5, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 6, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_corner( rotbox, 7, 3.0 ); dir = rotbox.center - pt; dir.normalize();
  ASSERT( rotbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and pointing away from box
  pt = scaled_face( rotbox, 0, 3.0 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 1, 3.0 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 2, 3.0 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 3, 3.0 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 4, 3.0 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 5, 3.0 ); dir = pt - rotbox.center; dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  
    // test ray starting outside face and parallel to face
  pt = scaled_face( rotbox, 0, 3.0 ); dir = rotbox.scaled_axis(1); dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 1, 3.0 ); dir = rotbox.scaled_axis(2); dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 2, 3.0 ); dir = rotbox.scaled_axis(0); dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 3, 3.0 ); dir = rotbox.scaled_axis(2); dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 4, 3.0 ); dir = rotbox.scaled_axis(0); dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
  pt = scaled_face( rotbox, 5, 3.0 ); dir = rotbox.scaled_axis(1); dir.normalize();
  ASSERT( ! rotbox.intersect_ray( pt, dir, TOL ) );
}

void test_build_from_tri()
{
  ErrorCode rval;
  Core moab;
  int i;
  
  // define a planar patch of triangles
  const double coords[] = { 0, 0, 0,
                            5, 0, 0,
                            5, 5, 0,
                            0, 5, 0,
                           -5, 5, 0,
                           -5, 0, 0,
                           -5,-5, 0,
                            0,-5, 0,
                            5,-5, 0,
                           10, 0, 0,
                            8, 5, 0,
                            5, 8, 0,
                            0,10, 0,
                           -5, 8, 0,
                           -8, 5, 0,
                          -10, 0, 0,
                           -8,-5, 0,
                           -5,-8, 0,
                           0,-10, 0,
                            5,-8, 0,
                            8,-5, 0 };
 const int conn[] = { 3,12,13,
                      3,11,12,
                      4,13,14,
                      4, 3,13,
                      3, 2,11,
                      2,10,11,
                      
                      5,14,15,
                      5, 4,14,
                      3, 4, 5,
                      0, 3, 5,
                      0, 1, 3,
                      1, 2, 3, 
                      2, 1,10,
                      1, 9,10,
                      
                      5,15,16,
                      6, 5,16,
                      5, 6, 7,
                      0, 5, 7,
                      1, 0, 7,
                      1, 7, 8,
                      1, 8,20,
                      9, 1,20,
                      
                      6,16,17,
                      7, 6,17,
                      7,17,18,
                      7,18,19,
                      8, 7,19,
                      8,19,20 };
  
    // build triangle mesh
  std::vector<EntityHandle> vertices(21);
  for (i = 0; i < 21; ++i) {
    rval = moab.create_vertex( coords + 3*i, vertices[i] );
    ASSERT(MB_SUCCESS == rval);
  }
  Range tris;
  for (i = 0; i < 28; ++i) {
    EntityHandle tri;
    EntityHandle c[3] = { vertices[conn[3*i  ]],
                            vertices[conn[3*i+1]],
                            vertices[conn[3*i+2]] };
    rval = moab.create_element( MBTRI, c, 3, tri );
    ASSERT(MB_SUCCESS == rval);
    tris.insert( tri );
  }
  
    // create box from triangles
  OrientedBox box;
  rval = OrientedBox::compute_from_2d_cells( box, &moab, tris );
  ASSERT( MB_SUCCESS == rval );
  
    // compute range along each box axis for input vertices
  const CartVect axis[3] = { box.scaled_axis(0),
                               box.scaled_axis(1),
                               box.scaled_axis(2) };
  double min[3], max[3];
  CartVect v = CartVect(coords) - box.center;
  min[0] = max[0] = box.scaled_axis(0) % v;
  min[1] = max[1] = box.scaled_axis(1) % v;
  min[2] = max[2] = box.scaled_axis(2) % v;
  for (i = 1; i < 21; ++i) {
    CartVect vi( coords + 3*i );
    CartVect v = vi - box.center;
    for (int j = 0; j < 3; ++j) {
      double d = (axis[j] % v) / (axis[j] % axis[j]);
      if (d < min[j])
        min[j] = d;
      if (d > max[j])
        max[j] = d;
    }
  }
  
    // Vrify that all points are contained in box 
    // and that box fits points tightly.
    // Triangles line in xy plane, so assuming axes are
    // sorted by length, first axis should be zero.
  ASSERT_DOUBLES_EQUAL( min[1], -1 );
  ASSERT_DOUBLES_EQUAL( min[2], -1 );
  ASSERT_DOUBLES_EQUAL( max[1],  1 );
  ASSERT_DOUBLES_EQUAL( max[2],  1 );
    
    // verify that the box is flat along first axis
  ASSERT( box.dimensions()[0] <= TOL );
    // verify that other two axes are in XY plane
  const CartVect z_axis(0.0,0.0,1.0);
  ASSERT( fabs(box.axis[1] % z_axis) <= TOL );
  ASSERT( fabs(box.axis[2] % z_axis) <= TOL );
}
                         
    

void test_build_from_pts()
{
  ErrorCode rval;
  Core moab;
  Interface* const gMB = &moab;
  
  const double vertex_coords[] =
    { 0, 0, 0,
      1, 0, 0,
      1, 1, 0,
      0, 1, 0,
      0, 0, 1,
      1, 0, 1,
      1, 1, 1,
      0, 1, 1 };
  const int num_double = sizeof(vertex_coords) / (sizeof(double));
  const int num_vertex = num_double / 3;
  assert (0 == num_double%3);
  
    // create some vertices
  Range vertices;
  double min[3] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double max[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < num_vertex; ++i) {
    EntityHandle h;
    rval = gMB->create_vertex( vertex_coords + 3*i, h );
    ASSERT( MB_SUCCESS == rval );
    vertices.insert(h);
    for (int j = 0; j < 3; ++j) {
      const double c = vertex_coords[3*i+j];
      if (c < min[j]) min[j] = c;
      if (c > max[j]) max[j] = c;
    }
  }
  
    // create box containing vertices
  OrientedBox box;
  rval = OrientedBox::compute_from_vertices( box, gMB, vertices );
  ASSERT( MB_SUCCESS == rval );
  
  for (int i = 0; i < num_vertex; ++i) {
    ASSERT( box.contained( CartVect(vertex_coords[3*i]), 1e-6 ) );
  }
  
    // define a set of points otside of the box to test
  double center[3] = { 0.5 * (max[0] + min[0]),
                       0.5 * (max[1] + min[1]),
                       0.5 * (max[2] + min[2]) };
  double diag[3] = { max[0] - min[0],
                     max[1] - min[1],
                     max[2] - min[2] };
  double outside1[3] = { center[0] + diag[0],
                         center[1] + diag[1],
                         center[2] + diag[2] };
  double outside2[3] = { center[0] + diag[0],
                         center[1] + diag[1],
                         center[2] };
  double outside3[3] = { center[0] + diag[0],
                         center[1],
                         center[2] };
  double *outside[3] = {outside1, outside2, outside3};
    // test 'contained' method for all points
  for (int i = 0; i < 3; ++i) {
    ASSERT( ! box.contained( CartVect(outside[i]), 1e-6 ) );
  } 
}

static void test_save()
{
  ErrorCode rval;
  Core moab;
  OrientedBox box;
  
    // create a tag to store the data in
  Tag tag;
  rval = OrientedBox::tag_handle( tag, &moab, "FOO" );
  ASSERT( MB_SUCCESS == rval );
  
    // check tag size
  int size;
  rval= moab.tag_get_bytes( tag, size );
  ASSERT( MB_SUCCESS == rval );
  ASSERT( size == sizeof(OrientedBox) );
}


/********************* Error Checking Code ***************************/

static void assert_vector_element( const CartVect& a, 
                                   const CartVect b[3], 
                                   const char* sa, 
                                   const char* sb, 
                                   int lineno )
{
  int i;
  for (i = 0; i < 3; ++i) 
    if (fabs(a[0] - b[i][0]) <= TOL
     && fabs(a[1] - b[i][1]) <= TOL
     && fabs(a[2] - b[i][2]) <= TOL)
      return;
      
  ++error_count;
  std::cerr << "Assertion failed at line " << lineno << std::endl
            << "\t" << sa << " in " << sb << std::endl
            << "\t" << sa << " = " << a << std::endl
            << "\t" << sb << " = " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;
}


void assert_vectors_equal( const CartVect& a, const CartVect& b, 
                           const char* sa, const char* sb,
                           int lineno )
{
  if (fabs(a[0] - b[0]) > TOL ||
      fabs(a[1] - b[1]) > TOL ||
      fabs(a[2] - b[2]) > TOL) {
    std::cerr << "Assertion failed at line " << lineno << std::endl
              << "\t" << sa << " == " << sb << std::endl
              << "\t[" << a[0] << ", " << a[1] << ", " << a[2] << "] == ["
              << b[0] << ", " << b[1] << ", " << b[2] << "]" << std::endl;
    ++error_count;
  }
}

void assert_doubles_equal( double a, double b, const char* sa, const char* sb, int lineno )
{
  if (fabs(a - b) > TOL) {
    std::cerr << "Assertion failed at line " << lineno << std::endl
              << "\t" << sa << " == " << sb << std::endl
              << "\t" << a << " == " << b << std::endl;
    ++error_count;
  }
}

void assert_bool( bool b, const char* sb, int lineno )
{
  if (!b) {
    std::cerr << "Assertion failed at line " << lineno << std::endl
              << "\t" << sb << std::endl;
    ++error_count;
  }
}
