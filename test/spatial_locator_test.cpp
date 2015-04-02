#include "moab/Core.hpp"
#include "moab/SpatialLocator.hpp"
#include "moab/Tree.hpp"
#include "moab/HomXform.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/CartVect.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/BVHTree.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include "TestUtil.hpp"

#include <cstdlib>
#include <sstream>

using namespace moab;

void test_kd_tree();
void test_bvh_tree();
void test_locator(SpatialLocator *sl);

ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n, int dim);

int max_depth = 30;
int npoints = 1000;
int leaf = 6;
bool print_tree = false;
int ints = 10;

int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#else
  // silence the warning of parameters not used, in serial; there should be a smarter way :(
  argv[0]=argv[argc-argc];
#endif

  ProgOptions po("spatial_locator_test options" );
  po.addOpt<int>( "ints,i", "Number of intervals on each side of scd mesh", &ints);
  po.addOpt<int>( "leaf,l", "Maximum number of elements per leaf", &leaf);
  po.addOpt<int>( "max_depth,m", "Maximum depth of tree", &max_depth);
  po.addOpt<int>( "npoints,n", "Number of query points", &npoints);
  po.addOpt<void>( "print,p", "Print tree details", &print_tree);
  po.parseCommandLine(argc, argv);

  RUN_TEST(test_kd_tree);
  RUN_TEST(test_bvh_tree);
  
#ifdef MOAB_HAVE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return 0;
}

void test_kd_tree() 
{
  ErrorCode rval;
  Core mb;
  
    // create a simple mesh to test
  Range elems;
  rval = create_hex_mesh(mb, elems, ints, 3); CHECK_ERR(rval);

    // initialize spatial locator with the elements and KDtree
  AdaptiveKDTree kd(&mb);
  std::ostringstream opts;
  opts << "MAX_DEPTH=" << max_depth << ";MAX_PER_LEAF=" << leaf;
  FileOptions fo(opts.str().c_str());
  rval = kd.parse_options(fo);
  SpatialLocator *sl = new SpatialLocator(&mb, elems, &kd);

  test_locator(sl);

    // test with an evaluator
  ElemEvaluator eval(&mb);
  kd.set_eval(&eval);
  test_locator(sl);

    // destroy spatial locator, and tree along with it
  delete sl;
}

void test_bvh_tree() 
{
  ErrorCode rval;
  Core mb;
  
    // create a simple mesh to test
  Range elems;
  rval = create_hex_mesh(mb, elems, ints, 3); CHECK_ERR(rval);

    // initialize spatial locator with the elements and a BVH tree
  BVHTree bvh(&mb);
  std::ostringstream opts;
  opts << "MAX_DEPTH=" << max_depth << ";MAX_PER_LEAF=" << leaf;
  FileOptions fo(opts.str().c_str());
  rval = bvh.parse_options(fo);
  SpatialLocator *sl = new SpatialLocator(&mb, elems, &bvh);
  test_locator(sl);
  
    // test with an evaluator
  ElemEvaluator eval(&mb);
  bvh.set_eval(&eval);
  test_locator(sl);

    // destroy spatial locator, and tree along with it
  delete sl;
}

void test_locator(SpatialLocator *sl) 
{
  CartVect box_del, test_pt, test_res;
  BoundBox box = sl->local_box();
  box_del = box.bMax - box.bMin;

  double denom = 1.0 / (double)RAND_MAX;
  int is_in;
  EntityHandle ent;
  ErrorCode rval;
  for (int i = 0; i < npoints; i++) {    
      // generate a small number of random point to test
    double rx = (double)rand() * denom, ry = (double)rand() * denom, rz = (double)rand() * denom;
    test_pt = box.bMin + CartVect(rx*box_del[0], ry*box_del[1], rz*box_del[2]);

    // call spatial locator to locate points
    rval = sl->locate_points(test_pt.array(), 1, &ent, test_res.array(), &is_in); CHECK_ERR(rval);

    // verify that the point was found
    CHECK_EQUAL(is_in, true);
  }

  std::cout << "Traversal stats:" << std::endl;
  sl->get_tree()->tree_stats().output_trav_stats();

  if (print_tree) {
    std::cout << "Tree information: " << std::endl;
    rval = sl->get_tree()->print();
    CHECK_ERR(rval);
  }
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
