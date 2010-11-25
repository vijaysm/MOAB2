#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"

using namespace moab;
using namespace moab::GeomUtil;

#include <iostream>

#include "TestUtil.hpp"
const double TOL = 1e-6;
#define ASSERT_VECTORS_EQUAL(A, B) assert_vectors_equal( (A), (B), #A, #B, __LINE__ )
#define ASSERT_DOUBLES_EQUAL(A, B) CHECK_REAL_EQUAL( A, B, TOL )
#define ASSERT(B) CHECK(B)


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
    FLAG_ERROR;
  }
}

void test_box_plane_norm( CartVect norm, 
                          CartVect min,
                          CartVect max )
{
  CartVect c_lower = min;
  CartVect c_upper = max;
  for (int i = 0; i < 3; ++i)
    if (norm[i] < 0.0)
      std::swap(c_lower[i],c_upper[i]);
  
  CartVect p_below = c_lower - norm;
  CartVect p_lower = c_lower + norm;
  CartVect p_upper = c_upper - norm;
  CartVect p_above = c_upper + norm;
  
  double below = -(p_below % norm);
  double lower = -(p_lower % norm);
  double upper = -(p_upper % norm);
  double above = -(p_above % norm);
  
  ASSERT( !box_plane_overlap( norm, below, min, max ) );
  ASSERT(  box_plane_overlap( norm, lower, min, max ) );
  ASSERT(  box_plane_overlap( norm, upper, min, max ) );
  ASSERT( !box_plane_overlap( norm, above, min, max ) );
}

void test_box_plane_axis( int axis, double ns, 
                          const CartVect& min, 
                          const CartVect& max )
{
  CartVect norm(0.0);
  norm[axis] = ns;
  test_box_plane_norm( norm, min, max );
}

void test_box_plane_edge( int axis1, int axis2, bool flip_axis2,
                          CartVect min, CartVect max )
{
  CartVect norm(0.0);
  norm[axis1] = max[axis1] - min[axis1];
  if (flip_axis2)
    norm[axis2] = min[axis2] - max[axis2];
  else
    norm[axis2] = max[axis2] - min[axis2];
  norm.normalize();
  
  test_box_plane_norm( norm, min, max );
}

void test_box_plane_corner( int xdir, int ydir, int zdir, 
                            CartVect min, CartVect max )
{
  CartVect norm(max - min);
  norm[0] *= xdir;
  norm[1] *= ydir;
  norm[2] *= zdir;
  test_box_plane_norm( norm, min, max );
}

void test_box_plane_overlap()
{
  const CartVect min( -1, -2, -3 );
  const CartVect max(  6,  4,  2 );
  
    // test with planes orthogonal to Z axis
  test_box_plane_axis( 2, 2.0, min, max );
    // test with planes orthogonal to X axis
  test_box_plane_axis( 1,-2.0, min, max );
    // test with planes orthogonal to Y axis
  test_box_plane_axis( 1, 1.0, min, max );

    // test with plane orthogonal to face diagonals
  test_box_plane_edge( 0, 1, true,  min, max );
  test_box_plane_edge( 0, 1, false, min, max );
  test_box_plane_edge( 0, 2, true,  min, max );
  test_box_plane_edge( 0, 2, false, min, max );
  test_box_plane_edge( 2, 1, true,  min, max );
  test_box_plane_edge( 2, 1, false, min, max );
  
    // test with plane orthogonal to box diagonals
  test_box_plane_corner( 1, 1, 1, min, max );
  test_box_plane_corner( 1, 1,-1, min, max );
  test_box_plane_corner( 1,-1,-1, min, max );
  test_box_plane_corner( 1,-1, 1, min, max );
}     


class ElemOverlapTest {
  public:
  
    virtual bool operator()( const CartVect* coords, 
                             const CartVect& box_center, 
                             const CartVect& box_dims ) const = 0;
};
class LinearElemOverlapTest : public ElemOverlapTest {
  public:
    const EntityType type;
    LinearElemOverlapTest(EntityType t) : type(t) {}
    bool operator()( const CartVect* coords, 
                     const CartVect& box_center, 
                     const CartVect& box_dims ) const
      { return box_linear_elem_overlap( coords, type, box_center, box_dims ); }
};
class TypeElemOverlapTest : public ElemOverlapTest {
  public:
    bool (*func)( const CartVect*, const CartVect&, const CartVect& );
    TypeElemOverlapTest( bool (*f)( const CartVect*, const CartVect&, const CartVect& ) )
      : func(f){}
    bool operator()( const CartVect* coords, 
                     const CartVect& box_center, 
                     const CartVect& box_dims ) const
      { return (*func)( coords, box_center, box_dims ); }
};    

void general_box_tri_overlap_test( const ElemOverlapTest& overlap )
{
  CartVect coords[3];
  CartVect center, dims;
  
    // test box projection within triangle, z-plane
  coords[0] = CartVect( 0, 0, 0 );
  coords[1] = CartVect( 0, 4, 0 );
  coords[2] = CartVect(-4, 0, 0 );
  center = CartVect( -2, 1, 0 );
  dims = CartVect( 1, 0.5, 3 );
  ASSERT(  overlap( coords, center, dims ) );
    // move box below plane of triangle
  center[2] = -4;
  ASSERT( !overlap( coords, center, dims ) );
    // move box above plane of triangle
  center[2] =  4;
  ASSERT( !overlap( coords, center, dims ) );
  
    // test box projection within triangle, x-plane
  coords[0] = CartVect( 3, 3, 0 );
  coords[1] = CartVect( 3, 3, 1 );
  coords[2] = CartVect( 3, 0, 0 );
  center = CartVect( 3, 2.5, .25 );
  dims = CartVect( 0.001, 0.4, .2 );
  ASSERT(  overlap( coords, center, dims ) );
    // move box below plane of triangle
  center[0] = 2;
  ASSERT( !overlap( coords, center, dims ) );
    // move box above plane of triangle
  center[0] = 4;
  ASSERT( !overlap( coords, center, dims ) );
  
    // test tri slices corner at +x,+y,+z
  coords[0] = CartVect(3,1,1);
  coords[1] = CartVect(1,3,1);
  coords[2] = CartVect(1,1,3);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri above the corner
  ASSERT( !overlap( coords, CartVect(0,0,0), CartVect(1,1,1) ) );
    // test tri slices corner at -x,-y,-z
  coords[0] = CartVect(-1,1,1);
  coords[1] = CartVect(1,-1,1);
  coords[2] = CartVect(1,1,-1);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri below the corner
  ASSERT( !overlap( coords, CartVect(2,2,2),CartVect(1,1,1) ) );
  
    // test tri slices corner at -x,+y,+z
  coords[0] = CartVect( 0.5, 0.0, 2.5);
  coords[1] = CartVect( 0.5, 2.5, 0.0);
  coords[2] = CartVect(-0.5, 0.0, 0.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri above the corner
  ASSERT( !overlap( coords, CartVect(2,1,1), CartVect(1,1,1) ) );
  
    // test tri slices corner at +x,-y,-z
  coords[0] = CartVect( 0.5, 0.0,-1.5);
  coords[1] = CartVect( 0.5,-1.5, 0.0);
  coords[2] = CartVect( 1.5, 0.0, 0.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri above the corner
  ASSERT( !overlap( coords, CartVect(0,1,1), CartVect(1,1,1) ) );

    // test tri slices corner at +x,-y,+z
  coords[0] = CartVect( 1.0, 1.0, 2.5 );
  coords[1] = CartVect( 2.5, 1.0, 1.0 );
  coords[2] = CartVect( 1.0,-0.5, 1.0 );
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri above the corner
  ASSERT( !overlap( coords, CartVect(-1,1,1), CartVect(1,1,1) ) );

    // test tri slices corner at -x,+y,-z  
  coords[0] = CartVect( 1.0,  1.0,-0.5 );
  coords[1] = CartVect(-0.5,  1.0, 1.0 );
  coords[2] = CartVect( 1.0,  2.5, 1.0 );
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri above the corner
  ASSERT( !overlap( coords, CartVect(3,1,1), CartVect(1,1,1) ) );

    // test tri slices corner at +x,+y,-z
  coords[0] = CartVect(-0.1, 1.0, 1.0);
  coords[1] = CartVect( 1.0,-0.1, 1.0);
  coords[2] = CartVect( 1.0, 1.0,-0.1);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,3), CartVect(1,1,1) ) );
  
    // test tri slices corner at -x,-y,+z
  coords[0] = CartVect( 2.1, 1.0, 1.0);
  coords[1] = CartVect( 1.0, 2.1, 1.0);
  coords[2] = CartVect( 1.0, 1.0, 2.1);
  ASSERT(  box_tri_overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,-1), CartVect(1,1,1) ) );
  
    // box edge parallel to x at +y,+z passes through triangle
  coords[0] = CartVect( 1.0, 1.0, 3.0);
  coords[1] = CartVect( 1.0, 3.0, 3.0);
  coords[2] = CartVect( 1.0, 3.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,0.3), CartVect(1,1,1) ) );
  
    // box edge parallel to x at +y,-z passes through triangle
  coords[0] = CartVect( 1.0, 3.0, 1.0);
  coords[1] = CartVect( 1.0, 3.0,-1.0);
  coords[2] = CartVect( 1.0, 1.0,-1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,1.7), CartVect(1,1,1) ) );
  
    // box edge parallel to x at -y,-z passes through triangle
  coords[0] = CartVect( 1.0,-1.0, 1.0);
  coords[1] = CartVect( 1.0,-1.0,-1.0);
  coords[2] = CartVect( 1.0, 1.0,-1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,1.7), CartVect(1,1,1) ) );
  
    // box edge parallel to x at -y,+z passes through triangle
  coords[0] = CartVect( 1.0,-1.0, 1.0);
  coords[1] = CartVect( 1.0,-1.0, 3.0);
  coords[2] = CartVect( 1.0, 1.0, 3.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,0.3), CartVect(1,1,1) ) );
  
    // box edge parallel to y at +x,+z passes through triangle
  coords[0] = CartVect( 1.0, 1.0, 3.0);
  coords[1] = CartVect( 3.0, 1.0, 3.0);
  coords[2] = CartVect( 3.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,0.3), CartVect(1,1,1) ) );
  
    // box edge parallel to y at -x,+z passes through triangle
  coords[0] = CartVect( 1.0, 1.0, 3.0);
  coords[1] = CartVect(-1.0, 1.0, 3.0);
  coords[2] = CartVect(-1.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,0.3), CartVect(1,1,1) ) );
  
    // box edge parallel to y at +x,-z passes through triangle
  coords[0] = CartVect( 1.0, 1.0,-1.0);
  coords[1] = CartVect( 3.0, 1.0,-1.0);
  coords[2] = CartVect( 3.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,1.7), CartVect(1,1,1) ) );
  
    // box edge parallel to y at -x,-z passes through triangle
  coords[0] = CartVect( 1.0, 1.0,-1.0);
  coords[1] = CartVect(-1.0, 1.0,-1.0);
  coords[2] = CartVect(-1.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1,1,1.7), CartVect(1,1,1) ) );
  
    // box edge parallel to z at +x,+y passes through triangle
  coords[0] = CartVect( 1.0, 3.0, 1.0);
  coords[1] = CartVect( 3.0, 3.0, 1.0);
  coords[2] = CartVect( 3.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(0.3,1,1), CartVect(1,1,1) ) );
  
    // box edge parallel to z at +x,-y passes through triangle
  coords[0] = CartVect( 1.0,-1.0, 1.0);
  coords[1] = CartVect( 3.0,-1.0, 1.0);
  coords[2] = CartVect( 3.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(0.3,1,1), CartVect(1,1,1) ) );
  
    // box edge parallel to z at -x,+y passes through triangle
  coords[0] = CartVect( 1.0, 3.0, 1.0);
  coords[1] = CartVect(-1.0, 3.0, 1.0);
  coords[2] = CartVect(-1.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1.7,1,1), CartVect(1,1,1) ) );
  
    // box edge parallel to z at -x,-y passes through triangle
  coords[0] = CartVect( 1.0,-1.0, 1.0);
  coords[1] = CartVect(-1.0,-1.0, 1.0);
  coords[2] = CartVect(-1.0, 1.0, 1.0);
  ASSERT(  overlap( coords, CartVect(1,1,1), CartVect(1,1,1) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(1.7,1,1), CartVect(1,1,1) ) );
  
    // triangle penetrates +x face
  coords[0] = CartVect( 2.0, 2.0, 2.0 );
  coords[1] = CartVect( 5.0, 3.0, 2.0 );
  coords[2] = CartVect( 5.0, 1.0, 2.0 );
  ASSERT(  overlap( coords, CartVect(2,2,2), CartVect(2,2,2) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(-1,2,2), CartVect(2,2,2) ) );
  
    // triangle penetrates -x face
  coords[0] = CartVect( 2.0, 2.0, 2.0 );
  coords[1] = CartVect(-1.0, 3.0, 2.0 );
  coords[2] = CartVect(-1.0, 1.0, 2.0 );
  ASSERT(  overlap( coords, CartVect(2,2,2), CartVect(2,2,2) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(5,2,2), CartVect(2,2,2) ) );
  
    // triangle penetrates +y face
  coords[0] = CartVect( 2.0, 2.0, 2.0 );
  coords[1] = CartVect( 3.0, 5.0, 2.0 );
  coords[2] = CartVect( 1.0, 5.0, 2.0 );
  ASSERT(  overlap( coords, CartVect(2,2,2), CartVect(2,2,2) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(2,-1,2), CartVect(2,2,2) ) );
  
    // triangle penetrates -y face
  coords[0] = CartVect( 2.0, 2.0, 2.0 );
  coords[1] = CartVect( 3.0,-1.0, 2.0 );
  coords[2] = CartVect( 1.0,-1.0, 2.0 );
  ASSERT(  overlap( coords, CartVect(2,2,2), CartVect(2,2,2) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(2,5,2), CartVect(2,2,2) ) );
  
    // triangle penetrates +z face
  coords[0] = CartVect( 2.0, 2.0, 2.0 );
  coords[1] = CartVect( 2.0, 3.0, 5.0 );
  coords[2] = CartVect( 2.0, 1.0, 5.0 );
  ASSERT(  overlap( coords, CartVect(2,2,2), CartVect(2,2,2) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(2,2,-1), CartVect(2,2,2) ) );
  
    // triangle penetrates -z face
  coords[0] = CartVect( 2.0, 2.0, 2.0 );
  coords[1] = CartVect( 2.0, 3.0,-1.0 );
  coords[2] = CartVect( 2.0, 1.0,-1.0 );
  ASSERT(  overlap( coords, CartVect(2,2,2), CartVect(2,2,2) ) );
    // test with tri outside box
  ASSERT( !overlap( coords, CartVect(2,2,5), CartVect(2,2,2) ) );
}

void general_box_hex_overlap_test( const ElemOverlapTest& overlap )
{
  CartVect coords[8];

    // test against axis-aligned rectilinear hex
  coords[0] = CartVect(-0.5,-0.5,-0.5);
  coords[1] = CartVect( 0.5,-0.5,-0.5);
  coords[2] = CartVect( 0.5, 0.5,-0.5);
  coords[3] = CartVect(-0.5, 0.5,-0.5);
  coords[4] = CartVect(-0.5,-0.5, 0.5);
  coords[5] = CartVect( 0.5,-0.5, 0.5);
  coords[6] = CartVect( 0.5, 0.5, 0.5);
  coords[7] = CartVect(-0.5, 0.5, 0.5);

  ASSERT( overlap( coords, CartVect( 0, 0, 0), CartVect(1,1,1) ) );

  ASSERT( overlap( coords, CartVect( 1, 0, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0, 1, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0, 0, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1, 0, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0,-1, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0, 0,-1), CartVect(1,1,1) ) );

  ASSERT( overlap( coords, CartVect( 1, 1, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1, 1, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1,-1, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 1,-1, 0), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 1, 0, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1, 0, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1, 0,-1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 1, 0,-1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0, 1, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0,-1, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0,-1,-1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 0, 1,-1), CartVect(1,1,1) ) );

  ASSERT( overlap( coords, CartVect( 1, 1, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1, 1, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1,-1, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 1,-1, 1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 1, 1,-1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1, 1,-1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect(-1,-1,-1), CartVect(1,1,1) ) );
  ASSERT( overlap( coords, CartVect( 1,-1,-1), CartVect(1,1,1) ) );

  ASSERT(!overlap( coords, CartVect( 3, 0, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0, 3, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0, 0, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3, 0, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0,-3, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0, 0,-3), CartVect(1,1,1) ) );

  ASSERT(!overlap( coords, CartVect( 3, 3, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3, 3, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3,-3, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 3,-3, 0), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 3, 0, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3, 0, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3, 0,-3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 3, 0,-3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0, 3, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0,-3, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0,-3,-3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 0, 3,-3), CartVect(1,1,1) ) );

  ASSERT(!overlap( coords, CartVect( 3, 3, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3, 3, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3,-3, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 3,-3, 3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 3, 3,-3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3, 3,-3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect(-3,-3,-3), CartVect(1,1,1) ) );
  ASSERT(!overlap( coords, CartVect( 3,-3,-3), CartVect(1,1,1) ) );

    // test against rectilinear hex rotated 45 degrees about z axis
  const double r = sqrt(2.0)/2.0;
  coords[0] = CartVect( r, 0,-0.5);
  coords[1] = CartVect( 0, r,-0.5);
  coords[2] = CartVect(-r, 0,-0.5);
  coords[3] = CartVect( 0,-r,-0.5);
  coords[4] = CartVect( r, 0, 0.5);
  coords[5] = CartVect( 0, r, 0.5);
  coords[6] = CartVect(-r, 0, 0.5);
  coords[7] = CartVect( 0,-r, 0.5);

  ASSERT( overlap( coords, CartVect( 1, 0, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT( overlap( coords, CartVect(-1, 0, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT( overlap( coords, CartVect( 0, 1, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT( overlap( coords, CartVect( 0,-1, 0 ), CartVect(0.5,0.5,0.5) ) );

  ASSERT(!overlap( coords, CartVect( 1, 0, 2 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect(-1, 0, 2 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect( 0, 1, 2 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect( 0,-1, 2 ), CartVect(0.5,0.5,0.5) ) );

  ASSERT(!overlap( coords, CartVect( 2, 0, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect(-2, 0, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect( 0, 2, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect( 0,-2, 0 ), CartVect(0.5,0.5,0.5) ) );

  ASSERT(!overlap( coords, CartVect( 1, 1, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect(-1, 1, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect(-1,-1, 0 ), CartVect(0.5,0.5,0.5) ) );
  ASSERT(!overlap( coords, CartVect( 1,-1, 0 ), CartVect(0.5,0.5,0.5) ) );

  ASSERT( overlap( coords, CartVect( 1, 1, 0 ), CartVect(0.75,0.75,0.5) ) );
  ASSERT( overlap( coords, CartVect(-1, 1, 0 ), CartVect(0.75,0.75,0.5) ) );
  ASSERT( overlap( coords, CartVect(-1,-1, 0 ), CartVect(0.75,0.75,0.5) ) );
  ASSERT( overlap( coords, CartVect( 1,-1, 0 ), CartVect(0.75,0.75,0.5) ) );
}

void general_box_tet_overlap_test( const ElemOverlapTest& overlap )
{
  CartVect coords[4];
  
    // Octant I
  coords[0] = CartVect(0,0,0);
  coords[1] = CartVect(1,0,0);
  coords[2] = CartVect(0,1,0);
  coords[3] = CartVect(0,0,1);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect(-1,-1,-1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(0.2,0.2,0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(0.5,0.5,0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(0.5,0.5,0.5), CartVect(0.15,0.15,0.15) ) );
  
    // Octant II
  coords[0] = CartVect(0,1,0);
  coords[1] = CartVect(-1,0,0);
  coords[2] = CartVect(0,0,0);
  coords[3] = CartVect(0,0,1);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect( 1,-1,-1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(-0.2,0.2,0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(-0.5,0.5,0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(-0.5,0.5,0.5), CartVect(0.15,0.15,0.15) ) );
  
    // Octant III
  coords[0] = CartVect(0,-1,0);
  coords[1] = CartVect(0,0,0);
  coords[2] = CartVect(-1,0,0);
  coords[3] = CartVect(0,0,1);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect( 1, 1,-1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(-0.2,-0.2,0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(-0.5,-0.5,0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(-0.5,-0.5,0.5), CartVect(0.15,0.15,0.15) ) );
  
    // Octant IV
  coords[0] = CartVect(1,0,0);
  coords[1] = CartVect(0,-1,0);
  coords[2] = CartVect(0,0,1);
  coords[3] = CartVect(0,0,0);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect(-1, 1,-1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(0.2,-0.2,0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(0.5,-0.5,0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(0.5,-0.5,0.5), CartVect(0.15,0.15,0.15) ) );
   
    // Octant V
  coords[0] = CartVect(0,0,0);
  coords[1] = CartVect(0,1,0);
  coords[2] = CartVect(1,0,0);
  coords[3] = CartVect(0,0,-1);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect(-1,-1, 1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(0.2,0.2,-0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(0.5,0.5,-0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(0.5,0.5,-0.5), CartVect(0.15,0.15,0.15) ) );
  
    // Octant VI
  coords[0] = CartVect(-1,0,0);
  coords[1] = CartVect(0,1,0);
  coords[2] = CartVect(0,0,0);
  coords[3] = CartVect(0,0,-1);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect( 1,-1, 1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(-0.2,0.2,-0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(-0.5,0.5,-0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(-0.5,0.5,-0.5), CartVect(0.15,0.15,0.15) ) );
  
    // Octant VII
  coords[0] = CartVect(0,0,0);
  coords[1] = CartVect(0,-1,0);
  coords[2] = CartVect(-1,0,0);
  coords[3] = CartVect(0,0,-1);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect( 1, 1, 1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(-0.2,-0.2,-0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(-0.5,-0.5,-0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(-0.5,-0.5,-0.5), CartVect(0.15,0.15,0.15) ) );
  
    // Octant VIII
  coords[0] = CartVect(0,-1,0);
  coords[1] = CartVect(1,0,0);
  coords[2] = CartVect(0,0,-1);
  coords[3] = CartVect(0,0,0);
    // tet entirely within box
  ASSERT( overlap( coords, CartVect(-1, 1, 1), CartVect(3,3,3) ) );
    // box entirely within tet
  ASSERT( overlap( coords, CartVect(0.2,-0.2,-0.2), CartVect(0.1,0.1,0.1) ) );
    // box corner penetrates tet face
  ASSERT( overlap( coords, CartVect(0.5,-0.5,-0.5), CartVect(0.2,0.2,0.2) ) );
    // box corner does not penetrate face
  ASSERT( !overlap( coords, CartVect(0.5,-0.5,-0.5), CartVect(0.15,0.15,0.15) ) );
 
  
    // Box edge -x,-z
  coords[0] = CartVect( 0, 0, 0);
  coords[1] = CartVect( 2,-1, 0);
  coords[2] = CartVect( 2, 1, 0);
  coords[3] = CartVect( 0, 0, 2);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(1.5,0.0,1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(2.5,0.0,2.5), CartVect(1,1,1) ) );
  
    // Box edge -y,-z
  coords[0] = CartVect( 1, 2, 0);
  coords[1] = CartVect(-1, 2, 0);
  coords[2] = CartVect( 0, 0, 0);
  coords[3] = CartVect( 0, 0, 2);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(0.0,1.5,1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(0.0,2.5,2.5), CartVect(1,1,1) ) );
  
    // Box edge +x,-z
  coords[0] = CartVect(-2,-1, 0);
  coords[1] = CartVect(-2, 1, 0);
  coords[2] = CartVect( 0, 0, 2);
  coords[3] = CartVect( 0, 0, 0);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(-1.5,0.0,1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(-2.5,0.0,2.5), CartVect(1,1,1) ) );
  
    // Box edge +y,-z
  coords[0] = CartVect( 2,-1, 0);
  coords[1] = CartVect( 0, 0, 0);
  coords[2] = CartVect(-2,-1, 0);
  coords[3] = CartVect( 0, 0, 2);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(0.0,-1.5,1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(0.0,-2.5,2.5), CartVect(1,1,1) ) );
  
    // Box edge -x,+z
  coords[0] = CartVect( 2,-1, 0);
  coords[1] = CartVect( 0, 0, 0);
  coords[2] = CartVect( 2, 1, 0);
  coords[3] = CartVect( 0, 0,-2);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(1.5,0.0,-1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(2.5,0.0,-2.5), CartVect(1,1,1) ) );
  
    // Box edge -y,+z
  coords[0] = CartVect(-1, 2, 0);
  coords[1] = CartVect( 1, 2, 0);
  coords[2] = CartVect( 0, 0, 0);
  coords[3] = CartVect( 0, 0,-2);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(0.0,1.5,-1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(0.0,2.5,-2.5), CartVect(1,1,1) ) );
  
    // Box edge +x,+z
  coords[0] = CartVect(-2, 1, 0);
  coords[1] = CartVect(-2,-1, 0);
  coords[2] = CartVect( 0, 0,-2);
  coords[3] = CartVect( 0, 0, 0);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(-1.5,0.0,-1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(-2.5,0.0,-2.5), CartVect(1,1,1) ) );
  
    // Box edge +y,+z
  coords[0] = CartVect( 0, 0, 0);
  coords[1] = CartVect( 2,-1, 0);
  coords[2] = CartVect(-2,-1, 0);
  coords[3] = CartVect( 0, 0,-2);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(0.0,-1.5,-1.5), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(0.0,-2.5,-2.5), CartVect(1,1,1) ) );
    
    // Box edge -x,-y
  coords[0] = CartVect( 0, 0, 0);
  coords[1] = CartVect( 0, 2,-1);
  coords[2] = CartVect( 0, 2, 1);
  coords[3] = CartVect( 2, 0, 0);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(1.5,1.5,0.0), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(2.5,2.5,0.0), CartVect(1,1,1) ) );
    
    // Box edge +x,-y
  coords[0] = CartVect( 0, 2,-1);
  coords[1] = CartVect( 0, 0, 0);
  coords[2] = CartVect( 0, 2, 1);
  coords[3] = CartVect(-2, 0, 0);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(-1.5,1.5,0.0), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(-2.5,2.5,0.0), CartVect(1,1,1) ) );
    
    // Box edge -x,+y
  coords[0] = CartVect( 0,-2, 1);
  coords[1] = CartVect( 0,-2,-1);
  coords[2] = CartVect( 0, 0, 0);
  coords[3] = CartVect( 2, 0, 0);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(1.5,-1.5,0.0), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(2.5,-2.5,0.0), CartVect(1,1,1) ) );
    
    // Box edge +x,+y
  coords[0] = CartVect( 0,-2,-1);
  coords[1] = CartVect(-2, 0, 0);
  coords[2] = CartVect( 0,-2, 1);
  coords[3] = CartVect( 0, 0, 0);
    // box edge passes through tet
  ASSERT( overlap( coords, CartVect(-1.5,-1.5,0.0), CartVect(1,1,1) ) );
    // box edge does not pass through tet
  ASSERT( !overlap( coords, CartVect(-2.5,-2.5,0.0), CartVect(1,1,1) ) );
  
  
    // Test tet edge through box 
  coords[0] = CartVect( -0.13369421660900116, -2.9871494770050049,  0.0526076555252075 );
  coords[1] = CartVect( -0.00350524857640266, -3.3236153125762939,  0.2924639880657196 );
  coords[2] = CartVect(  0.16473215818405151, -2.9966945648193359, -0.1936169415712357 );
  coords[3] = CartVect(  0.26740345358848572, -2.8492588996887207,  0.1519143134355545 );
  ASSERT( overlap( coords, CartVect( -2.5, -2.8, -2.5 ), CartVect( 2.5, 0.31, 2.5 ) ) );
}

void test_box_tri_overlap()
{
  general_box_tri_overlap_test( TypeElemOverlapTest(&box_tri_overlap) );
}

void test_box_linear_elem_overlap_tri()
{
  general_box_tri_overlap_test( LinearElemOverlapTest(MBTRI) );
}

void test_box_hex_overlap()
{
  general_box_hex_overlap_test( TypeElemOverlapTest(&box_hex_overlap) );
}

void test_box_linear_elem_overlap_hex()
{
  general_box_hex_overlap_test( LinearElemOverlapTest(MBHEX) );
}

void test_box_tet_overlap()
{
  general_box_tet_overlap_test( TypeElemOverlapTest(&box_tet_overlap) );
}

void test_box_linear_elem_overlap_tet()
{
  general_box_tet_overlap_test( LinearElemOverlapTest(MBTET) );
}


void test_ray_tri_intersect()
{
  bool xsect;
  double t;

    // define a triangle
  const CartVect tri[3] = { CartVect(1.0, 0.0, 0.0), 
                              CartVect(0.0, 1.0, 0.0),
                              CartVect(0.0, 0.0, 1.0) };
  
    // try a ray through the center of the triangle
  xsect = ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0, 1.0 ),
                             TOL, t );
  ASSERT(xsect);
  ASSERT_DOUBLES_EQUAL( 1.0/3.0, t );
  
    // try a same ray, but move base point above triangle
  xsect = ray_tri_intersect( tri, 
                             CartVect( 1.0, 1.0, 1.0 ),
                             CartVect( 1.0, 1.0, 1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
    // try a same ray the other direction with base point below triangle
  xsect = ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect(-1.0,-1.0,-1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
  
    // try a ray that passes above the triangle
  xsect = ray_tri_intersect( tri, 
                             CartVect( 1.0, 1.0, 1.0 ),
                             CartVect(-1.0,-1.0, 1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
    // try a skew ray
  xsect = ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0,-0.1 ),
                             TOL, t );
  ASSERT(!xsect);
}

void test_plucker_ray_tri_intersect()
{
  bool xsect;
  double t;

    // define a triangle
  const CartVect tri[3] = { CartVect(1.0, 0.0, 0.0), 
                              CartVect(0.0, 1.0, 0.0),
                              CartVect(0.0, 0.0, 1.0) };
  
    // try a ray through the center of the triangle
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0, 1.0 ),
                             TOL, t );
  ASSERT(xsect);
  ASSERT_DOUBLES_EQUAL( 1.0/3.0, t );
  
    // try a same ray, but move base point above triangle
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 1.0, 1.0, 1.0 ),
                             CartVect( 1.0, 1.0, 1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
    // try a same ray the other direction with base point below triangle
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect(-1.0,-1.0,-1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
  
    // try a ray that passes above the triangle
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 1.0, 1.0, 1.0 ),
                             CartVect(-1.0,-1.0, 1.0 ),
                             TOL, t );
  ASSERT(!xsect);
  
    // try a skew ray
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0,-0.1 ),
                             TOL, t );
  ASSERT(!xsect);

    // try a ray that intersects with wrong orientation
  const int orientation = -1.0;
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0, 1.0 ),
				     TOL, t, NULL, NULL, &orientation );
  ASSERT(!xsect);

    // try a ray that intersects beyond the nonneg_ray_len
  const double nonneg_ray_len = 0.25;
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0, 1.0 ),
				     TOL, t, &nonneg_ray_len );
  ASSERT(!xsect);

    // try a ray that intersects behind the origin
  const double neg_ray_len = -2.0;
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( -1.0, -1.0, -1.0 ),
				     TOL, t, NULL, &neg_ray_len );
  ASSERT(xsect);

    // try a ray that intersects a node
  intersection_type int_type;
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 0.0, 0.0 ),
				     TOL, t, NULL, NULL, NULL, &int_type);
  ASSERT(xsect);
  ASSERT(NODE0 == int_type);

    // try a ray that intersects an edge
  xsect = plucker_ray_tri_intersect( tri, 
                             CartVect( 0.0, 0.0, 0.0 ),
                             CartVect( 1.0, 1.0, 0.0 ),
				     TOL, t, NULL, NULL, NULL, &int_type);
  ASSERT(xsect);
  ASSERT(EDGE0 == int_type);
}

void test_closest_location_on_tri()
{
  CartVect result, input;
  
    // define a triangle
  const CartVect tri[3] = { CartVect(1.0, 0.0, 0.0), 
                              CartVect(0.0, 1.0, 0.0),
                              CartVect(0.0, 0.0, 1.0) };
  
    // try point at triangle centroid
  input = CartVect( 1.0/3.0, 1.0/3.0, 1.0/3.0 );
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
  input = CartVect(1.0,1.0,1.0);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, CartVect( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) );
  
    // try a point below the center of the triangle
  input = CartVect(0.0,0.0,0.0);
  closest_location_on_tri( input, tri, result );
  ASSERT_VECTORS_EQUAL( result, CartVect( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) );
  
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
  const CartVect tri_xy[3] = { CartVect( 0.0, sqrt(3.0)/2.0, 0.0), 
                                 CartVect( 0.5, 0.0, 0.0),
                                 CartVect(-0.5, 0.0, 0.0) };
  
    // for each vertex, test point that is
    // - outside triangle
    // - closest to vertex
    // - 'inside' one of the adjacent edges
    // - 'outside' the other adjacent edge
  closest_location_on_tri( CartVect(-0.3, 1.2, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[0] );
  closest_location_on_tri( CartVect( 0.3, 1.2, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[0] );
  closest_location_on_tri( CartVect( 1.0, 0.1, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[1] );
  closest_location_on_tri( CartVect( 0.6,-0.5, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[1] );
  closest_location_on_tri( CartVect(-0.6,-0.5, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[2] );
  closest_location_on_tri( CartVect(-1.0, 0.1, 0.0), tri_xy, result );
  ASSERT_VECTORS_EQUAL( result, tri_xy[2] );
}

void test_closest_location_on_polygon()
{
  CartVect result, input;
  
    // define a unit square in xy plane
  const CartVect quad[4] = { CartVect( 0.0, 0.0, 0.0), 
                               CartVect( 1.0, 0.0, 0.0),
                               CartVect( 1.0, 1.0, 0.0),
                               CartVect( 0.0, 1.0, 0.0) };
  
    // test input in center of square
  closest_location_on_polygon(  CartVect( 0.5, 0.5, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, CartVect( 0.5, 0.5, 0.0 ) );
    // test above center of square  
  closest_location_on_polygon(  CartVect( 0.5, 0.5, 1.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, CartVect( 0.5, 0.5, 0.0 ) );
    // test below center of square  
  closest_location_on_polygon(  CartVect( 0.5, 0.5,-1.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, CartVect( 0.5, 0.5, 0.0 ) );

    // test points within square, but not at center
  input = CartVect( 0.25, 0.25, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = CartVect( 0.75, 0.25, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = CartVect( 0.75, 0.75, 0 );
  closest_location_on_polygon( input, quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, input );
  input = CartVect( 0.25, 0.75, 0 );
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
  closest_location_on_polygon( CartVect(-1.0,-1.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[0] );
  closest_location_on_polygon( CartVect( 2.0,-1.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[1] );
  closest_location_on_polygon( CartVect( 2.0, 2.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[2] );
  closest_location_on_polygon( CartVect(-1.0, 2.0, 0.0 ), quad, 4, result );
  ASSERT_VECTORS_EQUAL( result, quad[3] );
  
    // test at point outside and closest to an edge
  CartVect x(1.0,0.0,0.0), y(0.0,1.0,0.0);
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

void test_segment_box_intersect()
{
  const double box_min = 0.0;
  const double box_max = 2.0;
  const double box_wid = box_max - box_min;
  const double box_mid = 0.5 * (box_min + box_max);
  const CartVect min( box_min );
  const CartVect max( box_max );
  const CartVect X(1,0,0), Y(0,1,0), Z(0,0,1);
  CartVect pt;
  double start, end;
  bool r;
  
    // test line through box in +x direction
  double offset = 1;
  pt = CartVect( box_min - offset, box_mid, box_mid );
  start = -HUGE_VAL; end = HUGE_VAL;
  r = segment_box_intersect( min, max, pt, X, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, box_min + offset );
  ASSERT_DOUBLES_EQUAL( end - start, box_wid);
  
    // test with ray ending left of the box
  start = -HUGE_VAL; end = 0;
  r = segment_box_intersect( min, max, pt, X, start, end );
  ASSERT( !r );
  
    // test with ray ending within box
  start = -HUGE_VAL; end = box_mid + offset;
  r = segment_box_intersect( min, max, pt, X, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, box_min + offset );
  ASSERT_DOUBLES_EQUAL( end, box_mid + offset );
  
    // test with ray beginning within box
  start = box_mid + offset; end = HUGE_VAL;
  r = segment_box_intersect( min, max, pt, X, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, box_mid + offset );
  ASSERT_DOUBLES_EQUAL( end, box_max + offset );
    
    // test with ray right of box
  start = offset + offset + box_max; end = HUGE_VAL;
  r = segment_box_intersect( min, max, pt, X, start, end );
  ASSERT( !r );
  
   
    // test line through box in -y direction
  offset = 1;
  pt = CartVect( box_mid, box_min - offset, box_mid );
  start = -HUGE_VAL; end = HUGE_VAL;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( end - start, box_wid );
  ASSERT_DOUBLES_EQUAL( end, box_min - offset);
  
    // test with ray ending left of the box
  start = box_min; end = HUGE_VAL;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( !r );
  
    // test with ray beginning within box
  start = -box_mid - offset; end = HUGE_VAL;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, -box_mid - offset );
  ASSERT_DOUBLES_EQUAL( end, box_min - offset );
  
    // test with ray ending within box
  start = -HUGE_VAL; end = -box_mid - offset;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, -box_max - offset );
  ASSERT_DOUBLES_EQUAL( end, -box_mid - offset );
    
    // test with ray right of box
  start = -HUGE_VAL; end = -box_max - offset - offset;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( !r );
 
    // test ray outside in Z direction, parallel to Z plane, and
    // intersecting in projections into other planes
  pt = CartVect( box_mid, box_mid, box_max + 1 );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt,  X, start, end );
  ASSERT( !r );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt, -X, start, end );
  ASSERT( !r );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt,  Y, start, end );
  ASSERT( !r );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( !r );
  
    // try the other side (less than the min Z);
  pt = CartVect( box_mid, box_mid, box_min - 1 );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt,  X, start, end );
  ASSERT( !r );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt, -X, start, end );
  ASSERT( !r );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt,  Y, start, end );
  ASSERT( !r );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( !r );
  
    // now move the ray such that it lies exactly on the side of the box
  pt = CartVect( box_mid, box_mid, box_min );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt,  X, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, 0 );
  ASSERT_DOUBLES_EQUAL( end, 0.5 * box_wid );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt, -X, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, 0 );
  ASSERT_DOUBLES_EQUAL( end, 0.5 * box_wid );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt,  Y, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, 0 );
  ASSERT_DOUBLES_EQUAL( end, 0.5 * box_wid );
  start = 0; end = box_wid;
  r = segment_box_intersect( min, max, pt, -Y, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, 0 );
  ASSERT_DOUBLES_EQUAL( end, 0.5 * box_wid );
  
    // try a skew line segment
  pt = CartVect( box_min - 0.25 * box_wid, box_mid, box_mid );
  CartVect dir( 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0 );
  start = 0; end = 1.5 / sqrt(2.0) * box_wid;
  r = segment_box_intersect( min, max, pt, dir, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, 0.5 / sqrt(2.0) * box_wid );
  ASSERT_DOUBLES_EQUAL( end, box_wid / sqrt(2.0) );
  
    // try with skew line segment that just touches edge of box
  pt = CartVect( box_min - 0.5 * box_wid, box_mid, box_mid );
  start = 0; end = 3.0 / sqrt(2.0) * box_wid;
  r = segment_box_intersect( min, max, pt, dir, start, end );
  ASSERT( r );
  ASSERT_DOUBLES_EQUAL( start, box_wid / sqrt(2.0) );
  ASSERT_DOUBLES_EQUAL( end, box_wid / sqrt(2.0) );

    // try with skew line segment outside of box
  pt = CartVect( box_min - 0.75 * box_wid, box_mid, box_mid );
  start = 0; end = 3.0 / sqrt(2.0) * box_wid;
  r = segment_box_intersect( min, max, pt, dir, start, end );
  ASSERT( !r );
}

void test_closest_location_on_box()
{
  const CartVect min(0,0,0), max(1,2,3);
  CartVect pt;
  
    // inside
  closest_location_on_box( min, max, CartVect(0.5,0.5,0.5), pt );
  ASSERT_VECTORS_EQUAL( CartVect(0.5,0.5,0.5), pt );
  
    // closest to min x side
  closest_location_on_box( min, max, CartVect(-1.0,0.5,0.5), pt );
  ASSERT_VECTORS_EQUAL( CartVect(0.0,0.5,0.5), pt );
  
    // closest to max x side
  closest_location_on_box( min, max, CartVect(2.0,0.5,0.5), pt );
  ASSERT_VECTORS_EQUAL( CartVect(1.0,0.5,0.5), pt );
  
    // closest to min y side
  closest_location_on_box( min, max, CartVect(0.5,-1.0,0.5), pt );
  ASSERT_VECTORS_EQUAL( CartVect(0.5,0.0,0.5), pt );
  
    // closest to max y side
  closest_location_on_box( min, max, CartVect(0.5,2.5,0.5), pt );
  ASSERT_VECTORS_EQUAL( CartVect(0.5,2.0,0.5), pt );
  
    // closest to min z side
  closest_location_on_box( min, max, CartVect(0.5,0.5,-0.1), pt );
  ASSERT_VECTORS_EQUAL( CartVect(0.5,0.5,0.0), pt );
  
    // closest to max z side
  closest_location_on_box( min, max, CartVect(0.5,0.5,100.0), pt );
  ASSERT_VECTORS_EQUAL( CartVect(0.5,0.5,3.0), pt );
  
    // closest to min corner
  closest_location_on_box( min, max, CartVect(-1,-1,-1), pt );
  ASSERT_VECTORS_EQUAL( min, pt );
  
    // closest to max corner
  closest_location_on_box( min, max, CartVect(2,3,4), pt );
  ASSERT_VECTORS_EQUAL( max, pt );
}

int main()
{
  int error_count = 0;
  error_count += RUN_TEST(test_box_plane_overlap);
  error_count += RUN_TEST(test_box_linear_elem_overlap_tri);
  error_count += RUN_TEST(test_box_linear_elem_overlap_tet);
  error_count += RUN_TEST(test_box_linear_elem_overlap_hex);
  error_count += RUN_TEST(test_box_tri_overlap);
  error_count += RUN_TEST(test_box_tet_overlap);
  error_count += RUN_TEST(test_box_hex_overlap);
  error_count += RUN_TEST(test_ray_tri_intersect);
  error_count += RUN_TEST(test_plucker_ray_tri_intersect);
  error_count += RUN_TEST(test_closest_location_on_tri);
  error_count += RUN_TEST(test_closest_location_on_polygon);
  error_count += RUN_TEST(test_segment_box_intersect);
  error_count += RUN_TEST(test_closest_location_on_box);
  return error_count;
}
