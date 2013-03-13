/**
 * \file elem_eval_test.cpp
 *
 * \brief test ElemEvaluator and the various element types in MOAB
 *
 */
#include "moab/Core.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/LinearHex.hpp"
#include "moab/LinearTet.hpp"
#include "moab/QuadraticHex.hpp"
#include "moab/CartVect.hpp"
#include "TestUtil.hpp"

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

using namespace moab;

void test_linear_hex();
void test_linear_tet();
void test_quadratic_hex();

CartVect hex_verts[] = { 
      // corners
    CartVect( -1, -1, -1 ), CartVect( 1, -1, -1 ), CartVect( 1, 1, -1 ), CartVect( -1, 1, -1 ),
    CartVect( -1, -1, 1 ), CartVect( 1, -1, 1 ), CartVect( 1, 1, 1 ), CartVect( -1, 1, 1 ),
      // mid-edge (bottom, middle, top)
    CartVect( 0, -1, -1 ), CartVect( 1, 0, -1 ), CartVect( 0, 1, -1 ), CartVect( -1, 0, -1 ),
    CartVect( -1, -1, 0 ), CartVect( 1, -1, 0 ), CartVect( 1, 1, 0 ), CartVect( -1, 1, 0 ),
    CartVect( 0, -1, 1 ), CartVect( 1, 0, 1 ), CartVect( 0, 1, 1 ), CartVect( -1, 0, 1 ),
      // mid-face (middle, bottom, top)
    CartVect( 0, -1, 0 ), CartVect( 1, 0, 0 ), CartVect( 0, 1, 0 ), CartVect( -1, 0, 0 ),
    CartVect( 0, 0, -1 ), CartVect( 0, 0, 1 ), 
      // mid-element
    CartVect( 0, 0, 0 )
};

void test_eval(ElemEvaluator &ee, bool test_integrate) 
{
  
  CartVect params, posn, params2;
  bool is_inside;
  const double EPS1 = 1e-6;
  Matrix3 jacob;
  ErrorCode rval;
  
  for (params[0] = -1; params[0] <= 1; params[0] += 0.2) {
    for (params[1] = -1; params[1] <= 1; params[1] += 0.2) {
      for (params[2] = -1; params[2] <= 1; params[2] += 0.2) {

          // forward/reverse evaluation should get back to the same point, within tol
        rval = ee.eval(params, posn.array()); CHECK_ERR(rval);
        rval = ee.reverse_eval(posn, EPS1, params2, &is_inside); CHECK_ERR(rval);
        CHECK_REAL_EQUAL(0.0, (params - params2).length(), EPS1);

          // jacobian should be >= 0
        rval = ee.jacobian(params, jacob); CHECK_ERR(rval);
        CHECK(jacob.determinant() >= 0.0);
        
      }
    }
  }

    // tag equal to coordinates should integrate to avg position, since volume is 1
  Tag tag;
  rval = ee.get_moab()->tag_get_handle(NULL, 3, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE | MB_TAG_CREAT); CHECK_ERR(rval);
  rval = ee.get_moab()->tag_set_data(tag, ee.get_vert_handles(), ee.get_num_verts(), hex_verts[0].array()); CHECK_ERR(rval);
  
  rval = ee.set_tag_handle(tag, 0); CHECK_ERR(rval);

  if (test_integrate) {
    CartVect integral, avg(0.0);
    rval = ee.integrate(integral.array()); CHECK_ERR(rval);
    CHECK_REAL_EQUAL(0.0, (avg - integral).length(), EPS1);
  }
  
}

int main()
{
  int failures = 0;
  
  failures += RUN_TEST(test_linear_hex);
  failures += RUN_TEST(test_quadratic_hex);
  failures += RUN_TEST(test_linear_tet);

  return failures;
}

void test_linear_hex() 
{
  Core mb;
  Range verts;
  ErrorCode rval = mb.create_vertices((double*)hex_verts, 8, verts); CHECK_ERR(rval);
  EntityHandle hex;
  std::vector<EntityHandle> connect;
  std::copy(verts.begin(), verts.end(), std::back_inserter(connect));
  rval = mb.create_element(MBHEX, connect.data(), 8, hex); CHECK_ERR(rval);
  
  ElemEvaluator ee(&mb, hex, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBHEX, LinearHex::eval_set());

  test_eval(ee, true);
}

void test_quadratic_hex() 
{
  Core mb;
  Range verts;
  ErrorCode rval = mb.create_vertices((double*)hex_verts, 27, verts); CHECK_ERR(rval);
  EntityHandle hex;
  std::vector<EntityHandle> connect;
  std::copy(verts.begin(), verts.end(), std::back_inserter(connect));
  rval = mb.create_element(MBHEX, connect.data(), 27, hex); CHECK_ERR(rval);
  
  ElemEvaluator ee(&mb, hex, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBHEX, QuadraticHex::eval_set());

  test_eval(ee, false);
}

void test_linear_tet() 
{
  Core mb;
  Range verts;
  ErrorCode rval = mb.create_vertices((double*)hex_verts[1].array(), 4, verts); CHECK_ERR(rval);
  EntityHandle tet;
  std::vector<EntityHandle> connect;
  std::copy(verts.begin(), verts.end(), std::back_inserter(connect));
  rval = mb.create_element(MBTET, connect.data(), 4, tet); CHECK_ERR(rval);
  
  ElemEvaluator ee(&mb, tet, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBTET, LinearTet::eval_set());

  test_eval(ee, false);
}
