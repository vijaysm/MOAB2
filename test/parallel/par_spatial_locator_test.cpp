#include "moab/Core.hpp"
#include "moab/SpatialLocator.hpp"
#include "moab/Tree.hpp"
#include "moab/HomXform.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/CartVect.hpp"
#include "moab/BVHTree.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ParallelComm.hpp"
#include "moab_mpi.h"

#include "TestUtil.hpp"

#include <cstdlib>
#include <sstream>
#include <string>

using namespace moab;

void test_kd_tree();
void test_bvh_tree();
void test_locator(SpatialLocator *sl);

ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n, int dim);
ErrorCode load_file(Interface &mb, std::string &fn, Range &elems);

int max_depth = 30;
int npoints = 1000;
int leaf = 6;
int tree = -1;
bool print_tree = false;
int ints = 10;
std::string fname;

int main(int argc, char **argv)
{
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;

  ProgOptions po("spatial_locator_test options" );
  po.addOpt<std::string>( "file,f", "Input filename", &fname);
  po.addOpt<int>( "ints,i", "Number of intervals on each side of scd mesh", &ints);
  po.addOpt<int>( "leaf,l", "Maximum number of elements per leaf", &leaf);
  po.addOpt<int>( "max_depth,m", "Maximum depth of tree", &max_depth);
  po.addOpt<int>( "npoints,n", "Number of query points", &npoints);
  po.addOpt<void>( "print,p", "Print tree details", &print_tree);
  po.addOpt<int>( "tree,t", "Tree type (-1=all (default), 0=AdaptiveKD, 1=BVH", &tree);
  po.parseCommandLine(argc, argv);

  if (-1 == tree || 0 == tree)
    RUN_TEST(test_kd_tree);
  if (-1 == tree || 1 == tree)
  RUN_TEST(test_bvh_tree);
  
  fail = MPI_Finalize();
  if (fail) return fail;

  return 0;
}

void test_kd_tree() 
{
  ErrorCode rval;
  Core mb;
  
    // create a simple mesh to test
  Range elems;
  if (fname.empty()) {
    rval = create_hex_mesh(mb, elems, ints, 3); CHECK_ERR(rval);
  }
  else {
    rval = load_file(mb, fname, elems); CHECK_ERR(rval);
  }

    // initialize spatial locator with the elements and the default tree type
  SpatialLocator *sl = new SpatialLocator(&mb, elems);

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
  if (fname.empty()) {
    rval = create_hex_mesh(mb, elems, ints, 3); CHECK_ERR(rval);
  }
  else {
    rval = load_file(mb, fname, elems); CHECK_ERR(rval);
  }

    // initialize spatial locator with the elements and a BVH tree
  BVHTree bvh(&mb);
  std::ostringstream opts;
  opts << "MAX_DEPTH=" << max_depth << ";MAX_PER_LEAF=" << leaf;
  FileOptions fo(opts.str().c_str());
  rval = bvh.parse_options(fo);
  SpatialLocator *sl = new SpatialLocator(&mb, elems, &bvh);
  test_locator(sl);
  
    // destroy spatial locator, and tree along with it
  delete sl;
}

bool is_neg(int is_neg) 
{
  return (is_neg < 0);
}

void test_locator(SpatialLocator *sl) 
{
  CartVect box_del;
  BoundBox box = sl->local_box();
  box_del = box.bMax - box.bMin;

  double denom = 1.0 / (double)RAND_MAX;
  std::vector<CartVect> test_pts(npoints);
  
  for (int i = 0; i < npoints; i++) {    
      // generate a small number of random point to test
    double rx = (double)rand() * denom, ry = (double)rand() * denom, rz = (double)rand() * denom;
    test_pts[i] = box.bMin + CartVect(rx*box_del[0], ry*box_del[1], rz*box_del[2]);
  }
  
    // call spatial locator to locate points
  ParallelComm *pc = ParallelComm::get_pcomm(sl->moab(), 0);
  CHECK(pc != NULL);
  
  ErrorCode rval = sl->par_locate_points(pc, test_pts[0].array(), npoints); CHECK_ERR(rval);
  if (pc->rank() == 0) {
    int num_out = std::count_if(sl->par_loc_table().vi_rd, sl->par_loc_table().vi_rd+2*npoints, is_neg);
    num_out /= 2;
  
    std::cout << "Number of points inside an element = " << npoints-num_out << "/" << npoints 
              << " (" << 100.0*((double)npoints-num_out)/npoints << "%)" << std::endl;
    std::cout << "Traversal stats:" << std::endl;
    sl->get_tree()->tree_stats().output();

    if (print_tree) {
      std::cout << "Tree information: " << std::endl;
      rval = sl->get_tree()->print();
      CHECK_ERR(rval);
    }
  }
}

ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n, int dim) 
{
  ScdInterface *scdi;
  ErrorCode rval = mb.query_interface(scdi); CHECK_ERR(rval);
  ScdParData spd;
  spd.gDims[0] = spd.gDims[1] = spd.gDims[2] = 0;
  spd.gDims[3] = n;
  spd.partMethod = ScdParData::SQIJK;
  if (dim > 1) spd.gDims[4] = n;
  if (dim > 2) spd.gDims[5] = n;
  ScdBox *new_box;
  rval = scdi->construct_box(HomCoord(0, 0, 0), HomCoord(0, 0, 0), 
                             NULL, 0, new_box, NULL, &spd, false, 0); CHECK_ERR(rval);

  rval = mb.get_entities_by_dimension(0, dim, elems); CHECK_ERR(rval);

  return rval;
}

ErrorCode load_file(Interface &mb, std::string &fn, Range &elems) 
{
  std::string options;
  ErrorCode rval;

  options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;";
  rval = mb.load_file(fn.c_str(), NULL, options.c_str());
  if (MB_SUCCESS != rval) return rval;
  
  rval = mb.get_entities_by_dimension(0, 3, elems);
  if (MB_SUCCESS != rval) return rval;

  if (elems.empty()) {
    rval = mb.get_entities_by_dimension(0, 2, elems);
    if (MB_SUCCESS != rval) return rval;
  }
  
  return MB_SUCCESS;
}

