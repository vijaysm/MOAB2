
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

  mat(0) = 2;
  mat(1) = -1;
  mat(2) = 0;
  mat(3) = -1;
  mat(4) = 2;
  mat(5) = -1;
  mat(6) = 0;
  mat(7) = -1;
  mat(8) = 2;

  //now do the Eigen Decomposition of this Matrix

  double lamda[3]; 
  moab::CartVect vectors[3];
  moab::ErrorCode rval = moab::Matrix::EigenDecomp( mat, lamda, vectors);
  CHECK_ERR(rval);

  //Hardcoded check values for the results
  double lamda_check[3];
  lamda_check[0] = 3.41421; lamda_check[1] = 2; lamda_check[2] = 0.585786;

  moab::CartVect vec0_check(0.5, -0.707107, 0.5);
  moab::CartVect vec1_check(0.707107, 3.37748e-17, -0.707107);
  moab::CartVect vec2_check(0.5, 0.707107, 0.5); 
  
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

  //for a real, symmetric matrix the Eigenvectors should be orthogonal
  CHECK_REAL_EQUAL( vectors[0]%vectors[1], 0 , tol );
  CHECK_REAL_EQUAL( vectors[0]%vectors[2], 0 , tol );
  
  return;
}


