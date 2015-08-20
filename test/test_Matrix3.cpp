
#include <vector>
#include <assert.h>
#include "moab/Matrix3.hpp"
#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "TestUtil.hpp"

using namespace moab;


void test_EigenDecomp();

int main ()
{

  int result = 0;

  result += RUN_TEST(test_EigenDecomp);
  
  return result; 

}

// test to ensure the Eigenvalues/vectors are calculated correctly and returned properly
// from the Matrix3 class for a simple case
void test_EigenDecomp()
{
  //Create a matrix
  moab::Matrix3 mat;

  mat(0) = 3;
  mat(1) = 2;
  mat(2) = 4;
  mat(3) = 2;
  mat(4) = 0;
  mat(5) = 2;
  mat(6) = 4;
  mat(7) = 2;
  mat(8) = 3;

  //now do the Eigen Decomposition of this Matrix

  double lamda[3]; 
  moab::CartVect vectors[3];
  moab::ErrorCode rval = moab::Matrix::EigenDecomp( mat, lamda, vectors);
  CHECK_ERR(rval);

  //Hardcoded check values for the results
  double lamda_check[3];
  lamda_check[0] = 8; lamda_check[1] = -1; lamda_check[2] = -1;

  moab::CartVect vec0_check( 0.666666, 0.3333333, 0.66666666);
  moab::CartVect vec1_check( 0.596285, 0.298142, -0.745356);
  moab::CartVect vec2_check( -0.447214, 0.894427, 0); 
  
  //now verfy that the returns Eigenvalues and Eigenvectors are correct (within some tolerance)
  double tol = 1e-04; 

  //check that the correct Eigenvalues are returned correctly (in order)
  CHECK_REAL_EQUAL( lamda[0], lamda_check[0], tol);
  CHECK_REAL_EQUAL( lamda[1], lamda_check[1], tol);
  CHECK_REAL_EQUAL( lamda[2], lamda_check[2], tol);

  //check the Eigenvector values (order should correspond to the Eigenvalues)
  //first vector
  CHECK_REAL_EQUAL( vectors[0][0], vec0_check[0], tol );
  CHECK_REAL_EQUAL( vectors[0][1], vec0_check[1], tol );
  CHECK_REAL_EQUAL( vectors[0][2], vec0_check[2], tol );
  
  //sceond vector
  CHECK_REAL_EQUAL( vectors[1][0], vec1_check[0], tol );
  CHECK_REAL_EQUAL( vectors[1][1], vec1_check[1], tol );
  CHECK_REAL_EQUAL( vectors[1][2], vec1_check[2], tol );

  //third vector
  CHECK_REAL_EQUAL( vectors[2][0], vec2_check[0], tol );
  CHECK_REAL_EQUAL( vectors[2][1], vec2_check[1], tol );
  CHECK_REAL_EQUAL( vectors[2][2], vec2_check[2], tol );

  //another check to ensure the result is valid (AM-kM = 0)
  unsigned int i;
  for(i=0; i<3; i++){
    moab::CartVect v = moab::Matrix::matrix_vector(mat, vectors[i])-lamda[i]*vectors[i];
    CHECK_REAL_EQUAL( v.length(), 0, tol );
  }
  
  return;
}
