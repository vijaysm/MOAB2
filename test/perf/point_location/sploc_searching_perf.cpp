#include "moab/Core.hpp"
#include "moab/SpatialLocator.hpp"
#include "moab/Tree.hpp"
#include "moab/HomXform.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/CartVect.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/BVHTree.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "moab/ElemEvaluator.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include <cstdlib>
#include <sstream>

using namespace moab;

ErrorCode test_locator(SpatialLocator &sl, int npoints, double rtol, double &cpu_time, double &percent_outside);
ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n, int dim);

int main(int argc, char **argv)
{
#ifdef MOAB_HAVE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#else
  // silence the warning of parameters not used, in serial; there should be a smarter way :(
  argv[0]=argv[argc-argc];
#endif

  int npoints = 100, dim = 3;
  int dints = 1, dleafs = 1, ddeps = 1;
  bool eval = false;
  double rtol = 1.0e-10;
  
  ProgOptions po("tree_searching_perf options" );
  po.addOpt<void>( ",e", "Use ElemEvaluator in tree search", &eval);
  po.addOpt<int>( "ints,i", "Number of doublings of intervals on each side of scd mesh", &dints);
  po.addOpt<int>( "leaf,l", "Number of doublings of maximum number of elements per leaf", &dleafs);
  po.addOpt<int>( "max_depth,m", "Number of 5-intervals on maximum depth of tree", &ddeps);
  po.addOpt<int>( "npoints,n", "Number of query points", &npoints);
  po.addOpt<double>( "tol,t", "Relative tolerance of point search", &rtol);
//  po.addOpt<void>( "print,p", "Print tree details", &print_tree);
  po.parseCommandLine(argc, argv);

  std::vector<int> ints, deps, leafs;
  ints.push_back(10);
  for (int i = 1; i < dints; i++) ints.push_back(2*ints[i-1]);
  deps.push_back(30);
  for (int i = 1; i < ddeps; i++) deps.push_back(deps[i-1]-5);
  leafs.push_back(6);
  for (int i = 1; i < dleafs; i++) leafs.push_back(2*leafs[i-1]);

  ErrorCode rval = MB_SUCCESS;
  std::cout << "Tree_type" << " "
            << "Elems_per_leaf" << " "
            << "Tree_depth" << " "
            << "Ints_per_side" << " "
            << "N_elements" << " "
            << "search_time" << " "
            << "perc_outside" << " "
            << "initTime" << " "
            << "nodesVisited" << " "
            << "leavesVisited" << " "
            << "numTraversals" << " "
            << "leafObjectTests" << std::endl;

// outermost iteration: # elements
  for (std::vector<int>::iterator int_it = ints.begin(); int_it != ints.end(); int_it++) {
    Core mb;
    Range elems;
    rval = create_hex_mesh(mb, elems, *int_it, 3);
    if (MB_SUCCESS != rval) return rval;
    
      // iteration: tree depth
    for (std::vector<int>::iterator dep_it = deps.begin(); dep_it != deps.end(); dep_it++) {
  
        // iteration: tree max elems/leaf
      for (std::vector<int>::iterator leafs_it = leafs.begin(); leafs_it != leafs.end(); leafs_it++) {
  
          // iteration: tree type
        for (int tree_tp = 0; tree_tp < 2; tree_tp++) {
            // create tree
          Tree *tree;
          if (0 == tree_tp)
            tree = new BVHTree(&mb);
          else
            tree = new AdaptiveKDTree(&mb);

          ElemEvaluator *eeval = NULL;
          if (eval) {
              // create an element evaluator
            eeval = new ElemEvaluator(&mb);
            rval = eeval->set_eval_set(*elems.begin());
            if (MB_SUCCESS != rval) return rval;
          }
          
          std::ostringstream opts;
          opts << "MAX_DEPTH=" << *dep_it << ";MAX_PER_LEAF=" << *leafs_it;
          FileOptions fo(opts.str().c_str());
          rval = tree->parse_options(fo);
          SpatialLocator sl(&mb, elems, tree, eeval);

            // call evaluation
          double cpu_time, perc_outside;
          rval = test_locator(sl, npoints, rtol, cpu_time, perc_outside);
          if (MB_SUCCESS != rval) return rval;

          std::cout << (tree_tp == 0 ? "BVH" : "KD") << " "
                    << *leafs_it << " "
                    << *dep_it << " "
                    << *int_it << " "
                    << (*int_it)*(*int_it)*(dim == 3 ? *int_it : 1) << " "
                    << cpu_time << " "
                    << perc_outside << " ";

          tree->tree_stats().output_all_stats();

          if (eeval) delete eeval;

        } // tree_tp

      } // max elems/leaf

    } // max depth

  } // # elements

  
#ifdef MOAB_HAVE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return 0;
}

ErrorCode test_locator(SpatialLocator &sl, int npoints, double rtol, double &cpu_time, double &percent_outside) 
{
  BoundBox box = sl.local_box();
  CartVect box_del = box.bMax - box.bMin;

  std::vector<CartVect> test_pts(npoints), test_res(npoints);
  std::vector<EntityHandle> ents(npoints);
  int *is_in = new int[npoints];

  double denom = 1.0 / (double)RAND_MAX;
  for (int i = 0; i < npoints; i++) {    
      // generate a small number of random point to test
    double rx = (double)rand() * denom, ry = (double)rand() * denom, rz = (double)rand() * denom;
    test_pts[i] = box.bMin + CartVect(rx*box_del[0], ry*box_del[1], rz*box_del[2]);
  }
  
  CpuTimer ct;
  
    // call spatial locator to locate points
  ErrorCode rval = sl.locate_points(test_pts[0].array(), npoints, &ents[0], test_res[0].array(), &is_in[0], rtol, 0.0);
  if (MB_SUCCESS != rval) return rval;

  cpu_time = ct.time_elapsed();

  int num_out = std::count(is_in, is_in+npoints, false);
  percent_outside = ((double)num_out)/npoints;
  delete [] is_in;
  
  return rval;
}

ErrorCode create_hex_mesh(Interface &mb, Range &elems, int n, int dim) 
{
  ScdInterface *scdi;
  ErrorCode rval = mb.query_interface(scdi);
  if (MB_SUCCESS != rval) return rval;
  
  HomCoord high(n-1, -1, -1);
  if (dim > 1) high[1] = n-1;
  if (dim > 2) high[2] = n-1;
  ScdBox *new_box;
  rval = scdi->construct_box(HomCoord(0, 0, 0), high, NULL, 0, new_box);
  if (MB_SUCCESS != rval) return rval;
  rval = mb.release_interface(scdi);
  if (MB_SUCCESS != rval) return rval;

  rval = mb.get_entities_by_dimension(0, dim, elems);
  if (MB_SUCCESS != rval) return rval;

  return rval;
}
