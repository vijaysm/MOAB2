#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "moab/BoundBox.hpp"
#include "moab/ScdInterface.hpp"

using namespace moab;

#include <iostream>

//! test add_entities, get_entities, num_entities
void test_bound_box();

int main()
{
  int err = 0;
  
  err += RUN_TEST(test_bound_box);
  
  if (!err) 
    printf("ALL TESTS PASSED\n");
  else
    printf("%d TESTS FAILED\n",err);
  
  return err;
}

void test_bound_box() 
{
  BoundBox box;
  CartVect vec(0.0);
  double tol = 0.0;
  std::vector<double> vals(6, 0.0);
  BoundBox other_box(&vals[0]);
  
    // test for contains point failure
  bool result = box.contains_point(vec.array(), tol);
  CHECK(!result);
  result = box.contains_box(other_box, tol);
  CHECK(!result);
  
  box = other_box;
  tol = 1.0e-10;

    // test for success
  result = box.contains_point(&vals[0], tol);
  CHECK(result);
  result = box.contains_box(other_box, tol);
  CHECK(result);
  
    // check update functions
  CartVect three(3.0);
  box.update_max(three.array());
  result = box.contains_point(three.array(), tol);
  CHECK(result);
  result = box.contains_point((three*1.1).array(), tol);
  CHECK(!result);
  result = box.contains_box(BoundBox(three, three), tol);
  CHECK(result);
  result = box.contains_box(BoundBox(three, 1.1*three), tol);
  CHECK(!result);
  
  CartVect negthree(-3.0);
  box.update_min(negthree.array());
  result = box.contains_point(negthree.array(), tol);
  CHECK(result);
  result = box.contains_point((negthree*1.1).array(), tol);
  CHECK(!result);
  result = box.contains_box(BoundBox(negthree, negthree), tol);
  CHECK(result);
  result = box.contains_box(BoundBox(1.1*negthree, negthree), tol);
  CHECK(!result);
  
  for (int i = 0; i < 3; i++) {
    vals[i] = -4.0; vals[3+i] = 4.0;
  }
  box.update(&vals[0]);
  result = box.contains_point(&vals[0], tol);
  CHECK(result);
  result = box.contains_point((CartVect(&vals[0])*1.1).array(), tol);
  CHECK(!result);
  result = box.contains_box(BoundBox(&vals[0]), tol);
  CHECK(result);
  result = box.contains_box(BoundBox(1.1*CartVect(&vals[0]), CartVect(&vals[3])), tol);
  CHECK(!result);
  
    // check length functions
  tol = 1.0e-6;
  
    // box should be 8 on a side, or 3*(2^3)^2 or 192 for diag length squared
  double diagsq = box.diagonal_squared();
  CHECK_REAL_EQUAL(diagsq, 192.0, tol);
  double diag = box.diagonal_length();
  CHECK_REAL_EQUAL(diag, sqrt(3)*8.0, tol);
  
    // check distance function

    // face-centered
  vals[0] = vals[1] = 0.0;
  vals[2] = 6.0;
  double dist = box.distance_squared(CartVect(&vals[0]).array());
  CHECK_REAL_EQUAL(dist, 4.0, tol);
  dist = box.distance(CartVect(&vals[0]).array());
  CHECK_REAL_EQUAL(dist, 2.0, tol);

    // edge-centered
  vals[0] = 0.0;
  vals[1] = vals[2] = 6.0;
  dist = box.distance_squared(CartVect(&vals[0]).array());
  CHECK_REAL_EQUAL(dist, 8.0, tol);
  dist = box.distance(CartVect(&vals[0]).array());
  CHECK_REAL_EQUAL(dist, sqrt(8.0), tol);

    // vertex-centered
  vals[0] = vals[1] = vals[2] = 6.0;
  dist = box.distance_squared(CartVect(&vals[0]).array());
  CHECK_REAL_EQUAL(dist, 12.0, tol);
  dist = box.distance(CartVect(&vals[0]).array());
  CHECK_REAL_EQUAL(dist, sqrt(12.0), tol);
  
    // check entity-based functions
  Core mb;
  ScdInterface *scdi;
  ErrorCode rval = mb.query_interface(scdi);
  CHECK_ERR(rval);
  ScdBox *scd_box;
    // create a 10x10x10 box
  rval = scdi->construct_box(HomCoord(0, 0, 0), HomCoord(10, 10, 10), NULL, 0, scd_box);
  CHECK_ERR(rval);
  
  EntityHandle vert = scd_box->start_vertex(), elem = scd_box->start_element();
  rval = box.update(mb, vert);
  CHECK_ERR(rval);
  rval = mb.get_coords(&vert, 1, &vals[0]);
  CHECK_ERR(rval);
  tol = 1.0e-10;
  result = box.contains_point(&vals[0], tol);
  CHECK(result);
  rval = box.update(mb, elem);
  CHECK_ERR(rval);
  rval = mb.get_coords(&elem, 1, &vals[0]);
  CHECK_ERR(rval);
  result = box.contains_point(&vals[0], tol);
  CHECK(result);
  
  Range all_elems(elem, elem + scd_box->num_elements() - 1);
  rval = box.update(mb, all_elems);
  CHECK_ERR(rval);
  CartVect center;
  box.compute_center(center);
  tol = 1.0e-6;
  CHECK_REAL_EQUAL((center - CartVect(5.0, 5.0, 5.0)).length(), 0.0, tol);
}

  
