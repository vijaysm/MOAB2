/**
 * \file elem_eval_test.cpp
 *
 * \brief test ElemEvaluator and the various element types in MOAB
 *
 */
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
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

const double EPS1 = 1e-6;

void test_eval(ElemEvaluator &ee, bool test_integrate) 
{
  
  CartVect params, posn, params2;
  bool is_inside;
  Matrix3 jacob;
  ErrorCode rval;
  
  for (params[0] = -1; params[0] <= 1; params[0] += 0.2) {
    for (params[1] = -1; params[1] <= 1; params[1] += 0.2) {
      for (params[2] = -1; params[2] <= 1; params[2] += 0.2) {

          // forward/reverse evaluation should get back to the same point, within tol
        rval = ee.eval(params.array(), posn.array()); CHECK_ERR(rval);
        rval = ee.reverse_eval(posn.array(), EPS1, params2.array(), &is_inside); CHECK_ERR(rval);
        CHECK_REAL_EQUAL(0.0, (params - params2).length(), EPS1);

          // jacobian should be >= 0
        rval = ee.jacobian(params.array(), jacob.array()); CHECK_ERR(rval);
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

void test_evals(ElemEvaluator &ee, bool test_integrate, EntityHandle *starth, int num_ents, double total_vol) 
{
  for (int i = 0; i < num_ents; i++)
    test_eval(ee, false);
    
  if (!test_integrate) return;
  if (!num_ents || !starth) CHECK_ERR(MB_FAILURE);

  double tot_vol = 0.0;
  EntityHandle eh = *starth;
  for (int i = 0; i < num_ents; i++) {
    double tmp_vol;
    ee.set_ent_handle(eh++);
    ErrorCode rval = ee.integrate(&tmp_vol); CHECK_ERR(rval);
    tot_vol += tmp_vol;
  }
  CHECK_REAL_EQUAL(total_vol, tot_vol, EPS1);
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
  Range verts, tets;
  ErrorCode rval = mb.create_vertices((double*)hex_verts[0].array(), 8, verts); CHECK_ERR(rval);
  EntityHandle starth = 1, *conn;
  int conn_inds[] = {1, 6, 4, 5,    1, 4, 6, 3,    0, 1, 3, 4,    1, 2, 3, 6,    3, 4, 6, 7};
  ReadUtilIface *ru;
  rval = mb.query_interface(ru); CHECK_ERR(rval);
  rval = ru->get_element_connect(5, 4, MBTET, 1, starth, conn); CHECK_ERR(rval);
  for (unsigned int i = 0; i < 20; i++) conn[i] = verts[conn_inds[i]];
  
  ElemEvaluator ee(&mb, starth, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBTET, LinearTet::eval_set());

  test_evals(ee, true, &starth, 5, 8.0);
}
