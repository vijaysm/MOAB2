#include "MBCore.hpp"
#include "MBGeomUtil.hpp"

using namespace MBGeomUtil;

#include <iostream>
#define ASSERT_VECTORS_EQUAL(A, B) assert_vectors_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT_DOUBLES_EQUAL(A, B) assert_doubles_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT(B) assert_bool( (B), #B, __LINE__ )

const double TOL = 1e-6;

int error_count = 0;

void assert_vectors_equal( const MBCartVect& a, const MBCartVect& b, 
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



void test_ray_tri_intersect()
{
  bool xsect;
  double t;

    // define a triangle
  const MBCartVect tri[3] = { MBCartVect(1.0, 0.0, 0.0), 
                              MBCartVect(0.0, 1.0, 0.0),
                              MBCartVect(0.0, 0.0, 1.0) };
  
    // try a ray through the center of the triangle
  xsect = ray_tri_intersect( tri, 
                             MBCartVect( 0.0, 0.0, 0.0 ),
                             MBCartVect( 1.0, 1.0, 1.0 ),
                             TOL, t );
  ASSERT(xsect);
  ASSERT_DOUBLES_EQUAL( 1.0/3.0, t );
  
    // try a same ray, but move base point above triangle
  xsect = ray_tri_intersect( tri, 
                             MBCartVect( 1.0, 1.0, 1.0 ),
                             MBCartVect( 1.0, 1.0, 1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
    // try a same ray the other direction with base point below triangle
  xsect = ray_tri_intersect( tri, 
                             MBCartVect( 0.0, 0.0, 0.0 ),
                             MBCartVect(-1.0,-1.0,-1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
  
    // try a ray that passes above the triangle
  xsect = ray_tri_intersect( tri, 
                             MBCartVect( 1.0, 1.0, 1.0 ),
                             MBCartVect(-1.0,-1.0, 1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
    // try a skew ray
  xsect = ray_tri_intersect( tri, 
                             MBCartVect( 0.0, 0.0, 0.0 ),
                             MBCartVect( 1.0, 1.0,-0.1 ),
                             TOL, t );
  ASSERT(!xsect);
}

void test_closest_location_on_tri()
{
  MBCartVect result, input;
  
    // define a triangle
  const MBCartVect tri[3] = { MBCartVect(1.0, 0.0, 0.0), 
                              MBCartVect(0.0, 1.0, 0.0),
                              MBCartVect(0.0, 0.0, 1.0) };
  
    // try point at triangle centroid
  input = MBCartVect( 1.0/3.0, 1.0/3.0, 1.0/3.0 );
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, input );
  
    // try point at each vertex
  closest_location_on_tri( tri[0], tri, result );
  ASSERT_VECTORS_EQUAL( result, tri[0] );
  closest_location_on_tri( tri[1], tri, result );
  ASSERT_VECTORS_EQUAL( result, tri[1] );
  closest_location_on_tri( tri[2], tri, result );
  ASSERT_VECTORS_EQUAL( result, tri[2] );
  
    // try point at center of each edge
  input = 0.5 * (tri[0] + tri[1]);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.5 * (tri[0] + tri[2]);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.5 * (tri[2] + tri[1]);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, input );
  
    // try a point above the center of the triangle
  input = MBCartVect(1.0,1.0,1.0);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, MBCartVect( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) );
  
    // try a point below the center of the triangle
  input = MBCartVect(0.0,0.0,0.0);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, MBCartVect( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) );
  
    // try a point closest to each vertex and 'outside' of both adjacent edges.
  input = 2*tri[0];
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, tri[0] );
  input = 2*tri[1];
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, tri[1] );
  input = 2*tri[2];
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, tri[2] );
  
    // try a point outside and closest to each edge
  input = tri[0] + tri[1];
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, 0.5 * input );
  input = tri[2] + tri[1];
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, 0.5 * input );
  input = tri[0] + tri[2];
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, 0.5 * input );
  
    // define an equilateral triangle in the xy-plane
  const MBCartVect tri_xy[3] = { MBCartVect( 0.0, sqrt(3.0)/2.0, 0.0), 
                                 MBCartVect( 0.5, 0.0, 0.0),
                                 MBCartVect(-0.5, 0.0, 0.0) };
  
    // for each vertex, test point that is
    // - outside triangle
    // - closest to vertex
    // - 'inside' one of the adjacent edges
    // - 'outside' the other adjacent edge
  closest_location_on_tri( MBCartVect(-0.3, 1.2, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[0] );
  closest_location_on_tri( MBCartVect( 0.3, 1.2, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[0] );
  closest_location_on_tri( MBCartVect( 1.0, 0.1, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[1] );
  closest_location_on_tri( MBCartVect( 0.6,-0.5, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[1] );
  closest_location_on_tri( MBCartVect(-0.6,-0.5, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[2] );
  closest_location_on_tri( MBCartVect(-1.0, 0.1, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[2] );
}

void test_closest_location_on_polygon()
{
  MBCartVect result, input;
  
    // define a unit square in xy plane
  const MBCartVect quad[4] = { MBCartVect( 0.0, 0.0, 0.0), 
                               MBCartVect( 1.0, 0.0, 0.0),
                               MBCartVect( 1.0, 1.0, 0.0),
                               MBCartVect( 0.0, 1.0, 0.0) };
  
    // test input in center of square
  closest_location_on_polygon(  MBCartVect( 0.5, 0.5, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, MBCartVect( 0.5, 0.5, 0.0 ) );
    // test above center of square  
  closest_location_on_polygon(  MBCartVect( 0.5, 0.5, 1.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, MBCartVect( 0.5, 0.5, 0.0 ) );
    // test below center of square  
  closest_location_on_polygon(  MBCartVect( 0.5, 0.5,-1.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, MBCartVect( 0.5, 0.5, 0.0 ) );

    // test points within square, but not at center
  input = MBCartVect( 0.25, 0.25, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = MBCartVect( 0.75, 0.25, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = MBCartVect( 0.75, 0.75, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = MBCartVect( 0.25, 0.75, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );

    // test at each corner
  closest_location_on_polygon(  quad[0], quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[0] );
  closest_location_on_polygon(  quad[1], quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[1] );
  closest_location_on_polygon(  quad[2], quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[2] );
  closest_location_on_polygon(  quad[3], quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[3] );
  
    // test at point on each edge
  input = 0.5 * quad[0] + 0.5 * quad[1];
  closest_location_on_polygon(  input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.2 * quad[1] + 0.8 * quad[2];
  closest_location_on_polygon(  input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.7 * quad[2] + 0.3 * quad[3];
  closest_location_on_polygon(  input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.6 * quad[3] + 0.4 * quad[0];
  closest_location_on_polygon(  input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  
    // test at point outside and closest to each corner
  closest_location_on_polygon( MBCartVect(-1.0,-1.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[0] );
  closest_location_on_polygon( MBCartVect( 2.0,-1.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[1] );
  closest_location_on_polygon( MBCartVect( 2.0, 2.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[2] );
  closest_location_on_polygon( MBCartVect(-1.0, 2.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[3] );
  
    // test at point outside and closest to an edge
  MBCartVect x(1.0,0.0,0.0), y(0.0,1.0,0.0);
  input = 0.5 * quad[0] + 0.5 * quad[1];
  closest_location_on_polygon(  input-y, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.2 * quad[1] + 0.8 * quad[2];
  closest_location_on_polygon(  input+x, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.7 * quad[2] + 0.3 * quad[3];
  closest_location_on_polygon(  input+y, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = 0.6 * quad[3] + 0.4 * quad[0];
  closest_location_on_polygon(  input-x, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
}

int main()
{
  test_ray_tri_intersect();
  test_closest_location_on_tri();
  test_closest_location_on_polygon();
  return error_count;
}
