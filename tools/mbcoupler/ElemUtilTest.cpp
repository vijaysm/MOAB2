#include "TestUtil.hpp"
#include "ElemUtil.hpp"
#include <iostream>

using namespace moab;

void test_hex_nat_coords();

int main()
{
  int rval = 0;
  rval += RUN_TEST(test_hex_nat_coords);
  return rval;
}

const CartVect cube_corners[8] = { CartVect( 0, 0, 0 ),
                                     CartVect( 1, 0, 0 ),
                                     CartVect( 1, 1, 0 ),
                                     CartVect( 0, 1, 0 ),
                                     CartVect( 0, 0, 1 ),
                                     CartVect( 1, 0, 1 ),
                                     CartVect( 1, 1, 1 ),
                                     CartVect( 0, 1, 1 ) };
                                    

const CartVect hex_corners[8] = { CartVect( 1.0, 0.0, 0.0 ),
                                    CartVect( 1.0, 1.0, 0.3 ),
                                    CartVect( 0.0, 2.0, 0.6 ),
                                    CartVect( 0.2, 1.1, 0.4 ),
                                    CartVect( 1.5, 0.3, 1.0 ),
                                    CartVect( 1.5, 1.3, 1.0 ),
                                    CartVect( 0.5, 2.3, 1.0 ),
                                    CartVect( 0.7, 1.4, 1.0 ) };

/** shape function for trilinear hex */
CartVect hex_map( const CartVect& xi, const CartVect* corners )
{
  CartVect x(0.0);
  x += (1 - xi[0]) * (1 - xi[1]) * (1 - xi[2]) * corners[0];
  x += (1 + xi[0]) * (1 - xi[1]) * (1 - xi[2]) * corners[1];
  x += (1 + xi[0]) * (1 + xi[1]) * (1 - xi[2]) * corners[2];
  x += (1 - xi[0]) * (1 + xi[1]) * (1 - xi[2]) * corners[3];
  x += (1 - xi[0]) * (1 - xi[1]) * (1 + xi[2]) * corners[4];
  x += (1 + xi[0]) * (1 - xi[1]) * (1 + xi[2]) * corners[5];
  x += (1 + xi[0]) * (1 + xi[1]) * (1 + xi[2]) * corners[6];
  x += (1 - xi[0]) * (1 + xi[1]) * (1 + xi[2]) * corners[7];
  return x *= 0.125;
}

static void hex_bounding_box( const CartVect* corners, CartVect& min, CartVect& max  )
{
  min = max = corners[0];
  for (unsigned i = 1; i < 8; ++i)
    for (unsigned d = 0; d < 3; ++d)
      if (corners[i][d] < min[d])
        min[d] = corners[i][d];
      else if (corners[i][d] > max[d])
        max[d] = corners[i][d];
}

static bool in_range( const CartVect& xi )
  { return fabs(xi[0]) <= 1 
        && fabs(xi[1]) <= 1 
        && fabs(xi[2]) <= 1; 
  }        

void test_hex_nat_coords()
{
  CartVect xi, result_xi;
  bool valid;
  const double EPS = 1e-6;
  
    // first test with cube because it's easier to debug failures
  for (xi[0] = -1; xi[0] <= 1; xi[0] += 0.2) {
    for (xi[1] = -1; xi[1] <= 1; xi[1] += 0.2) {
      for (xi[2] = -1; xi[2] <= 1; xi[2] += 0.2) {
        const CartVect pt = hex_map(xi, cube_corners);
        valid = ElemUtil::nat_coords_trilinear_hex( cube_corners, pt, result_xi, EPS/10 );
        CHECK(valid);
        CHECK_REAL_EQUAL( xi[0], result_xi[0], EPS );
        CHECK_REAL_EQUAL( xi[1], result_xi[1], EPS );
        CHECK_REAL_EQUAL( xi[2], result_xi[2], EPS );
      }
    }
  }
  
    // now test with distorted hex
  for (xi[0] = -1; xi[0] <= 1; xi[0] += 0.2) {
    for (xi[1] = -1; xi[1] <= 1; xi[1] += 0.2) {
      for (xi[2] = -1; xi[2] <= 1; xi[2] += 0.2) {
        const CartVect pt = hex_map(xi, hex_corners);
        valid = ElemUtil::nat_coords_trilinear_hex( hex_corners, pt, result_xi, EPS/10 );
        CHECK(valid);
        CHECK_REAL_EQUAL( xi[0], result_xi[0], EPS );
        CHECK_REAL_EQUAL( xi[1], result_xi[1], EPS );
        CHECK_REAL_EQUAL( xi[2], result_xi[2], EPS );
      }
    }
  }
  
    // test points outside of element
  CartVect x, min, max;
  hex_bounding_box( cube_corners, min, max );
  for (x[0] = -1; x[0] <= 2; x[0] += 0.4) {
    for (x[1] = -1; x[1] <= 2; x[1] += 0.4) {
      for (x[2] = -1; x[2] <= 2; x[2] += 0.4) {
        bool in_box = x[0] >= min[0] && x[0] <= max[0] 
                   && x[1] >= min[1] && x[1] <= max[1]
                   && x[2] >= min[2] && x[2] <= max[2];
        if (in_box)
          continue;
        valid = ElemUtil::nat_coords_trilinear_hex( cube_corners, x, result_xi, EPS/10 );
//std::cout << (valid ? 'y' : 'n');
        CHECK(!valid || !in_range(result_xi));
      }
    }
  }
//std::cout << std::endl;

  hex_bounding_box( hex_corners, min, max );
  for (x[0] = -1; x[0] <= 3; x[0] += 0.5) {
    for (x[1] = -2; x[1] <= 4; x[1] += 0.5) {
      for (x[2] = -1; x[2] <= 2; x[2] += 0.4) {
        bool in_box = x[0] >= min[0] && x[0] <= max[0] 
                   && x[1] >= min[1] && x[1] <= max[1]
                   && x[2] >= min[2] && x[2] <= max[2];
        if (in_box)
          continue;
        valid = ElemUtil::nat_coords_trilinear_hex( hex_corners, x, result_xi, EPS/10 );
//std::cout << (valid ? 'y' : 'n');
        CHECK(!valid || !in_range(result_xi));
      }
    }
  }
//std::cout << std::endl;
}
  
