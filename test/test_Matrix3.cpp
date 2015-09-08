
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
  for (int i=0; i < 3; ++i)
    vectors[i].normalize();

  // check to ensure that the eigenvalues are ordered from highest to lowest
  CHECK( lamda[0] >= lamda[1]);
  CHACK( lamda[1] >= lamda[2]);
  
  // check to ensure the result is valid (AM-kM = 0)
  for(unsigned i=0; i<3; ++i) {
    moab::CartVect v = moab::Matrix::matrix_vector(mat, vectors[i])-lamda[i]*vectors[i];
    CHECK_REAL_EQUAL( v.length(), 0, tol );
  }

  //for a real, symmetric matrix the Eigenvectors should be orthogonal
  CHECK_REAL_EQUAL( vectors[0]%vectors[1], 0 , tol );
  CHECK_REAL_EQUAL( vectors[0]%vectors[2], 0 , tol );
  CHECK_REAL_EQUAL( vectors[1]%vectors[2], 0 , tol );
  
  return;
}


