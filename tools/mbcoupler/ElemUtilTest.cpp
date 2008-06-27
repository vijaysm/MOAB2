#include "TestUtil.hpp"
#include "MBElemUtil.hpp"
#include <iostream>

void test_hex_nat_coords();

int main()
{
  int rval = 0;
  rval += RUN_TEST(test_hex_nat_coords);
  return rval;
}

const MBCartVect cube_corners[8] = { MBCartVect( 0, 0, 0 ),
                                     MBCartVect( 1, 0, 0 ),
                                     MBCartVect( 1, 1, 0 ),
                                     MBCartVect( 0, 1, 0 ),
                                     MBCartVect( 0, 0, 1 ),
                                     MBCartVect( 1, 0, 1 ),
                                     MBCartVect( 1, 1, 1 ),
                                     MBCartVect( 0, 1, 1 ) };
                                    

const MBCartVect hex_corners[8] = { MBCartVect( 1.0, 0.0, 0.0 ),
                                    MBCartVect( 1.0, 1.0, 0.3 ),
                                    MBCartVect( 0.0, 2.0, 0.6 ),
                                    MBCartVect( 0.2, 1.1, 0.4 ),
                                    MBCartVect( 1.5, 0.3, 1.0 ),
                                    MBCartVect( 1.5, 1.3, 1.0 ),
                                    MBCartVect( 0.5, 2.3, 1.0 ),
                                    MBCartVect( 0.7, 1.4, 1.0 ) };

/** shape function for trilinear hex */
MBCartVect hex_map( const MBCartVect& xi, const MBCartVect* corners )
{
  MBCartVect x(0.0);
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


void test_hex_nat_coords()
{
  MBCartVect xi, result_xi;
  bool valid;
  const double EPS = 1e-6;
  
    // first test with cube because it's easier to debug failures
  for (xi[0] = -1; xi[0] <= 1; xi[0] += 0.2) {
    for (xi[1] = -1; xi[1] <= 1; xi[1] += 0.2) {
      for (xi[2] = -1; xi[2] <= 1; xi[2] += 0.2) {
        const MBCartVect pt = hex_map(xi, cube_corners);
        valid = MBElemUtil::nat_coords_trilinear_hex( cube_corners, pt, result_xi, EPS/10 );
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
        const MBCartVect pt = hex_map(xi, hex_corners);
        valid = MBElemUtil::nat_coords_trilinear_hex( hex_corners, pt, result_xi, EPS/10 );
        CHECK(valid);
        CHECK_REAL_EQUAL( xi[0], result_xi[0], EPS );
        CHECK_REAL_EQUAL( xi[1], result_xi[1], EPS );
        CHECK_REAL_EQUAL( xi[2], result_xi[2], EPS );
      }
    }
  }
}
  
