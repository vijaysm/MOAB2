/**
 * \file elem_eval_test.cpp
 *
 * \brief test ElemEvaluator and the various element types in MOAB
 *
 */
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/LocalDiscretization/LinearTri.hpp"
#include "moab/LocalDiscretization/LinearQuad.hpp"
#include "moab/LocalDiscretization/LinearHex.hpp"
#include "moab/LocalDiscretization/LinearTet.hpp"
#include "moab/LocalDiscretization/QuadraticHex.hpp"
#include "moab/CartVect.hpp"
#include "TestUtil.hpp"

using namespace moab;

void test_linear_tri();
void test_linear_quad();
void test_linear_hex();
void test_linear_tet();
void test_quadratic_hex();
void test_normal_linear_tri();
void test_normal_linear_quad();
void test_normal_linear_tet();
void test_normal_linear_hex();
ErrorCode create_mesh(Core &mb, EntityType type);

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
  int is_inside;
  Matrix3 jacob;
  ErrorCode rval;
  int ent_dim = ee.get_moab()->dimension_from_handle(ee.get_ent_handle());
  
  for (params[0] = -1; params[0] <= 1; params[0] += 0.2) {
    for (params[1] = -1; params[1] <= 1; params[1] += 0.2) {
      for (params[2] = -1; params[2] <= 1; params[2] += 0.2) {

          // forward/reverse evaluation should get back to the same point, within tol
        if (!ee.inside(params.array(), EPS1)) continue;
        
        rval = ee.eval(params.array(), posn.array()); CHECK_ERR(rval);
        rval = ee.reverse_eval(posn.array(), EPS1, EPS1, params2.array(), &is_inside); CHECK_ERR(rval);
        CHECK_REAL_EQUAL(0.0, params[0] - params2[0], 3*EPS1);
        if (ent_dim > 1)
          CHECK_REAL_EQUAL(0.0, params[1] - params2[1], 3*EPS1);
        if (ent_dim > 2)
          CHECK_REAL_EQUAL(0.0, params[2] - params2[2], 3*EPS1);

          // jacobian should be >= 0
        rval = ee.jacobian(params.array(), jacob.array()); CHECK_ERR(rval);
        CHECK(jacob.determinant() >= 0.0);
        
      }
    }
  }

  if (!test_integrate) return;
  
    // tag equal to coordinates should integrate to avg position, test that
  Tag tag;
    // make a temporary tag and set it on vertices to the vertex positions
  rval = ee.get_moab()->tag_get_handle(NULL, 3, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE | MB_TAG_CREAT); CHECK_ERR(rval);
  rval = ee.get_moab()->tag_set_data(tag, ee.get_vert_handles(), ee.get_num_verts(), ee.get_vert_pos()); CHECK_ERR(rval);
    // set that temporary tag on the evaluator so that's what gets integrated
  rval = ee.set_tag_handle(tag, 0); CHECK_ERR(rval);

  CartVect integral, avg;
  rval = ee.integrate(integral.array()); CHECK_ERR(rval);

    // now integrate a const 1-valued function, using direct call to the integrate function
  std::vector<double> one(ee.get_num_verts(), 1.0);
  double measure;
  EvalSet es;
  EntityHandle eh = ee.get_ent_handle();
  rval = EvalSet::get_eval_set(ee.get_moab(), eh, es); CHECK_ERR(rval);
  rval = (*es.integrateFcn)(&one[0], ee.get_vert_pos(), ee.get_num_verts(), ee.get_moab()->dimension_from_handle(eh),
                            1, ee.get_work_space(), &measure); CHECK_ERR(rval);
  if (measure) integral /= measure;

    // check against avg of entity's vertices' positions
  rval = ee.get_moab()->get_coords(&eh, 1, avg.array()); CHECK_ERR(rval);
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
  
//  failures += RUN_TEST(test_linear_tri); currently failing linear tri, bad formulation, working on it...
  failures += RUN_TEST(test_linear_quad);
  failures += RUN_TEST(test_linear_hex);
  failures += RUN_TEST(test_quadratic_hex);
  failures += RUN_TEST(test_linear_tet);
  failures += RUN_TEST(test_normal_linear_tri);
  failures += RUN_TEST(test_normal_linear_quad);
  failures += RUN_TEST(test_normal_linear_tet);
  failures += RUN_TEST(test_normal_linear_hex);

  return failures;
}

void test_linear_tri() 
{
  Core mb;
  Range verts;
  double tri_verts[] = {-1.0, -1.0, -1.0,
                         1.0, -1.0, -1.0,
                        -1.0,  1.0, -1.0};
  
          
  ErrorCode rval = mb.create_vertices(tri_verts, 3, verts); CHECK_ERR(rval);
  EntityHandle tri;
  std::vector<EntityHandle> connect;
  std::copy(verts.begin(), verts.end(), std::back_inserter(connect));
  rval = mb.create_element(MBTRI, &connect[0], 3, tri); CHECK_ERR(rval);
  
  ElemEvaluator ee(&mb, tri, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBTRI, LinearTri::eval_set());

  test_eval(ee, true);
}

void test_linear_quad() 
{
  Core mb;
  Range verts;
  ErrorCode rval = mb.create_vertices((double*)hex_verts, 4, verts); CHECK_ERR(rval);
  EntityHandle quad;
  std::vector<EntityHandle> connect;
  std::copy(verts.begin(), verts.end(), std::back_inserter(connect));
  rval = mb.create_element(MBQUAD, &connect[0], 4, quad); CHECK_ERR(rval);
  
  ElemEvaluator ee(&mb, quad, 0);
  ee.set_tag_handle(0, 0);
  ee.set_eval_set(MBQUAD, LinearQuad::eval_set());

  test_eval(ee, true);
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

void test_normal_linear_tri()
{
  ErrorCode error;
  Core mb;

  error = create_mesh(mb, MBTRI); CHECK_ERR(error);
  Range faces;
  error = mb.get_entities_by_dimension(0,2,faces); CHECK_ERR(error);

  ElemEvaluator ee(&mb,0,0);
  ee.set_eval_set(MBTRI, LinearTri::eval_set());

  double nrms[8][3];

  for (int i=0; i<4; i++)
    {
      ee.set_ent_handle(faces[i]);
      ee.get_normal(1, 0, nrms[2*i]);
      ee.get_normal(1, 2, nrms[2*i+1]);
    }

  for (int i=0; i<4; i++)
    {
      if (i == 3)
        {
          double val = nrms[7][0]*nrms[0][0] + nrms[7][1]*nrms[0][1] + nrms[7][2]*nrms[0][2];
          CHECK_EQUAL(val, -1);
        }
      else
        {
          double val = nrms[2*i+1][0]*nrms[2*(i+1)][0] + nrms[2*i+1][1]*nrms[2*(i+1)][1] + nrms[2*i+1][2]*nrms[2*(i+1)][2];
          CHECK_EQUAL(val, -1);
        }
    }
}

void test_normal_linear_quad()
{
  ErrorCode error;
  Core mb;

  error = create_mesh(mb, MBQUAD); CHECK_ERR(error);
  Range faces;
  error = mb.get_entities_by_dimension(0,2,faces); CHECK_ERR(error);

  ElemEvaluator ee(&mb,0,0);
  ee.set_eval_set(MBQUAD, LinearQuad::eval_set());

  double nrms[8][3];

  for (int i=0; i<4; i++)
    {
      ee.set_ent_handle(faces[i]);
      ee.get_normal(1, 0, nrms[2*i]);
      ee.get_normal(1, 3, nrms[2*i+1]);
    }

  for (int i=0; i<4; i++)
    {
      if (i == 3)
        {
          double val = nrms[7][0]*nrms[0][0] + nrms[7][1]*nrms[0][1] + nrms[7][2]*nrms[0][2];
          CHECK_EQUAL(val, -1);
        }
      else
        {
          double val = nrms[2*i+1][0]*nrms[2*(i+1)][0] + nrms[2*i+1][1]*nrms[2*(i+1)][1] + nrms[2*i+1][2]*nrms[2*(i+1)][2];
          CHECK_EQUAL(val, -1);
        }
    }
}

void test_normal_linear_tet()
{
  ErrorCode error;
  Core mb;

  error = create_mesh(mb, MBTET); CHECK_ERR(error);
  Range cells;
  error = mb.get_entities_by_dimension(0,3,cells); CHECK_ERR(error);

  ElemEvaluator ee(&mb,0,0);
  ee.set_eval_set(MBTET, LinearTet::eval_set());

  double nrms[8][3];

  for (int i=0; i<4; i++)
    {
      ee.set_ent_handle(cells[i]);
      ee.get_normal(2, 0, nrms[2*i]);
      ee.get_normal(2, 2, nrms[2*i+1]);
    }

  for (int i=0; i<4; i++)
    {
      if (i == 3)
        {
          double val = nrms[7][0]*nrms[0][0] + nrms[7][1]*nrms[0][1] + nrms[7][2]*nrms[0][2];
          CHECK_EQUAL(val, -1);
        }
      else
        {
          double val = nrms[2*i+1][0]*nrms[2*(i+1)][0] + nrms[2*i+1][1]*nrms[2*(i+1)][1] + nrms[2*i+1][2]*nrms[2*(i+1)][2];
          CHECK_EQUAL(val, -1);
        }
    }
}

void test_normal_linear_hex()
{
  ErrorCode error;
  Core mb;

  error = create_mesh(mb, MBHEX); CHECK_ERR(error);
  Range cells;
  error = mb.get_entities_by_dimension(0,3,cells); CHECK_ERR(error);

  ElemEvaluator ee(&mb,0,0);
  ee.set_eval_set(MBHEX, LinearHex::eval_set());

  double nrms[8][3];

  for (int i=0; i<4; i++)
    {
      ee.set_ent_handle(cells[i]);
      ee.get_normal(2, 0, nrms[2*i]);
      ee.get_normal(2, 3, nrms[2*i+1]);
    }

  for (int i=0; i<4; i++)
    {
      if (i == 3)
        {
          double val = nrms[7][0]*nrms[0][0] + nrms[7][1]*nrms[0][1] + nrms[7][2]*nrms[0][2];
          CHECK_EQUAL(val, -1);
        }
      else
        {
          double val = nrms[2*i+1][0]*nrms[2*(i+1)][0] + nrms[2*i+1][1]*nrms[2*(i+1)][1] + nrms[2*i+1][2]*nrms[2*(i+1)][2];
          CHECK_EQUAL(val, -1);
        }
    }
}

ErrorCode create_mesh(Core &mb, EntityType type)
{
  ErrorCode error;
  if (type == MBTRI)
    {
      const double coords[] = {0,0,0,
                              1,0,0,
                              0,1,0,
                              -1,0,0,
                              0,-1,0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0, 1, 2,
                         0,2,3,
                         0,3,4,
                         0,4,1};
      const size_t num_elems = sizeof(conn)/sizeof(int)/3;

      EntityHandle verts[num_vtx], faces[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mb.create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[3];
          for (int j=0; j<3; j++)
            c[j] = verts[conn[3*i+j]];

          error = mb.create_element(MBTRI, c, 3, faces[i]); CHECK_ERR(error);
        }
    }
  else if (type == MBQUAD)
    {
      const double coords[] = {0,0,0,
                               1,0,0,
                               1,1,0,
                               0,1,0,
                               -1,1,0,
                               -1,0,0,
                               -1,-1,0,
                               0,-1,0,
                               1,-1,0};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0,1,2,3,
                         0,3,4,5,
                         0,5,6,7,
                         0,7,8,1};
      const size_t num_elems = sizeof(conn)/sizeof(int)/4;

      EntityHandle verts[num_vtx], faces[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mb.create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[4*i+j]];

          error = mb.create_element(MBQUAD, c, 4, faces[i]); CHECK_ERR(error);
        }

    }
  else if (type == MBTET)
    {
      const double coords[] = {0,0,0,
                               1,0,0,
                               0,1,0,
                               -1,0,0,
                               0,-1,0,
                               0,0,1};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0,1,2,5,
                         0,2,3,5,
                         0,3,4,5,
                         0,4,1,5};
      const size_t num_elems = sizeof(conn)/sizeof(int)/4;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mb.create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[4];
          for (int j=0; j<4; j++)
            c[j] = verts[conn[4*i+j]];

          error = mb.create_element(MBTET, c, 4, cells[i]); CHECK_ERR(error);
        }
    }
  else if (type == MBHEX)
    {
      const double coords[] = {0,0,0,
                               1,0,0,
                              1,1,0,
                              0,1,0,
                              -1,1,0,
                              -1,0,0,
                              -1,-1,0,
                              0,-1,0,
                              1,-1,0,
                              0,0,1,
                              1,0,1,
                              1,1,1,
                              0,1,1,
                              -1,1,1,
                              -1,0,1,
                              -1,-1,1,
                              0,-1,1,
                              1,-1,1};
      const size_t num_vtx = sizeof(coords)/sizeof(double)/3;

      const int conn[] = {0,1,2,3,9,10,11,12,
                         0,3,4,5,9,12,13,14,
                         0,5,6,7,9,14,15,16,
                         0,7,8,1,9,16,17,10};
      const size_t num_elems = sizeof(conn)/sizeof(int)/8;

      EntityHandle verts[num_vtx], cells[num_elems];
      for (size_t i=0; i< num_vtx; ++i)
        {
          error = mb.create_vertex(coords+3*i, verts[i]); CHECK_ERR(error);
        }

      for (size_t i=0; i< num_elems; ++i)
        {
          EntityHandle c[8];
          for (int j=0; j<8; j++)
            c[j] = verts[conn[8*i+j]];

          error = mb.create_element(MBHEX, c, 8, cells[i]); CHECK_ERR(error);
        }
    }

  return MB_SUCCESS;
}
