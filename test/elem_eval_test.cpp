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

const double EPS1 = 1.0e-6;

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
        if (!ee.inside(params.array(), EPS1)) continue;
        
        rval = ee.eval(params.array(), posn.array()); CHECK_ERR(rval);
        rval = ee.reverse_eval(posn.array(), EPS1, params2.array(), &is_inside); CHECK_ERR(rval);
        if ((params - params2).length() > 3*EPS1) 
          std::cerr << params << std::endl;
        CHECK_REAL_EQUAL(0.0, (params - params2).length(), 3*EPS1);

          // jacobian should be >= 0
        rval = ee.jacobian(params.array(), jacob.array()); CHECK_ERR(rval);
        CHECK(jacob.determinant() >= 0.0);
        
      }
    }
  }

  if (!test_integrate) return;
  
    // tag equal to coordinates should integrate to avg position, since volume is 1
  Tag tag;
    // make a temporary tag and set it on vertices to the vertex positions
  rval = ee.get_moab()->tag_get_handle(NULL, 3, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE | MB_TAG_CREAT); CHECK_ERR(rval);
  rval = ee.get_moab()->tag_set_data(tag, ee.get_vert_handles(), ee.get_num_verts(), hex_verts[0].array()); CHECK_ERR(rval);
    // set that temporary tag on the evaluator so that's what gets integrated
  rval = ee.set_tag_handle(tag, 0); CHECK_ERR(rval);

  CartVect integral, avg(0.0);
  rval = ee.integrate(integral.array()); CHECK_ERR(rval);
  CHECK_REAL_EQUAL(0.0, (avg - integral).length(), EPS1);

    // delete the temporary tag
  rval = ee.get_moab()->tag_delete(tag);
  CHECK_ERR(rval);
    // set the ee's tag back to coords
  rval = ee.set_tag_handle(0, 0); CHECK_ERR(rval);
}

void test_evals(ElemEvaluator &ee, bool test_integrate, EntityHandle *ents, int num_ents, Tag onetag, double total_vol) 
{
  for (int i = 0; i < num_ents; i++) {
    ee.set_ent_handle(ents[i]);
    test_eval(ee, false);
  }
    
  if (!test_integrate) return;
  if (!num_ents || !ents) CHECK_ERR(MB_FAILURE);

  ErrorCode rval = ee.set_tag_handle(onetag, 0); CHECK_ERR(rval);
  
  double tot_vol = 0.0;
  for (int i = 0; i < num_ents; i++) {
    double tmp_vol;
    ee.set_ent_handle(ents[i]);
    rval = ee.integrate(&tmp_vol); CHECK_ERR(rval);
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
  rval = mb.create_element(MBHEX, &connect[0], 8, hex); CHECK_ERR(rval);
  
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
  rval = mb.create_element(MBHEX, &connect[0], 27, hex); CHECK_ERR(rval);
  
  ElemEvaluator ee(&mb, hex, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBHEX, QuadraticHex::eval_set());

  test_eval(ee, false);
}

void test_linear_tet() 
{
  Core mb;
  Range verts;
  ErrorCode rval = mb.create_vertices((double*)hex_verts[0].array(), 8, verts); CHECK_ERR(rval);
  EntityHandle starth = 1, *conn;
  int conn_inds[] = {1, 6, 4, 5,    1, 4, 6, 3,    0, 1, 3, 4,    1, 2, 3, 6,    3, 4, 6, 7};
  ReadUtilIface *ru;
  rval = mb.query_interface(ru); CHECK_ERR(rval);
  rval = ru->get_element_connect(5, 4, MBTET, 1, starth, conn); CHECK_ERR(rval);
  for (unsigned int i = 0; i < 20; i++) conn[i] = verts[conn_inds[i]];
  EntityHandle tets[5];
  for (unsigned int i = 0; i < 5; i++) tets[i] = starth + i;
  ElemEvaluator ee(&mb, 0, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBTET, LinearTet::eval_set());

    // make a tag whose value is one on all vertices
  Tag tag;
  rval = mb.tag_get_handle(NULL, 1, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE | MB_TAG_CREAT); CHECK_ERR(rval);
  std::vector<double> vals(verts.size(), 1.0);
  rval = mb.tag_set_data(tag, verts, &vals[0]); CHECK_ERR(rval);

  test_evals(ee, true, tets, 5, tag, 8.0);
}
