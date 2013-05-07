#include "moab/Core.hpp"
#include "moab/SpatialLocator.hpp"
#include "moab/Tree.hpp"
#include "moab/HomXform.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/CartVect.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include "TestUtil.hpp"

#include <cstdlib>

using namespace moab;

const int npoints = 1000;

void test_just_tree();
ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n = 10, int dim = 3);

int main(int argc, char **argv)
{
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#else
  // silence the warning of parameters not used, in serial; there should be a smarter way :(
  argv[0]=argv[argc-argc];
#endif

  RUN_TEST(test_just_tree);
  
#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return 0;
}

void test_just_tree() 
{
  ErrorCode rval;
  Core mb;
  
    // create a simple mesh to test
  Range elems;
  rval = create_hex_mesh(mb, elems); CHECK_ERR(rval);

    // initialize spatial locator with the elements and the default tree type
  SpatialLocator *sl = new SpatialLocator(&mb, elems);
  CartVect box_min, box_max, box_del, test_pt, test_res;
  rval = sl->get_bounding_box(box_min, box_max); CHECK_ERR(rval);
  box_del = box_max - box_min;

  double denom = 1.0 / (double)RAND_MAX;
  bool is_in;
  EntityHandle ent;
  for (int i = 0; i < npoints; i++) {    
      // generate a small number of random point to test
    double rx = (double)rand() * denom, ry = (double)rand() * denom, rz = (double)rand() * denom;
    test_pt = box_min + CartVect(rx*box_del[0], ry*box_del[1], rz*box_del[2]);

    // call spatial locator to locate points
    rval = sl->locate_points(test_pt.array(), 1, &ent, test_res.array(), 0.0, 0.0, &is_in); CHECK_ERR(rval);

    // verify that the point was found
    CHECK_EQUAL(is_in, true);
  }

    // destroy spatial locator, and tree along with it
  delete sl;
}

ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n, int dim) 
{
  ScdInterface *scdi;
  ErrorCode rval = mb.query_interface(scdi); CHECK_ERR(rval);
  HomCoord high(n-1, -1, -1);
  if (dim > 1) high[1] = n-1;
  if (dim > 2) high[2] = n-1;
  ScdBox *new_box;
  rval = scdi->construct_box(HomCoord(0, 0, 0), high, NULL, 0, new_box); CHECK_ERR(rval);
  rval = mb.release_interface(scdi); CHECK_ERR(rval);

  rval = mb.get_entities_by_dimension(0, dim, elems); CHECK_ERR(rval);

  return rval;
}
