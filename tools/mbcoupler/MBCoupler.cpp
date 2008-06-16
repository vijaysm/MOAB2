#include "MBCoupler.hpp"
#include "MBParallelComm.hpp"

#include "assert.h"

MBCoupler::MBCoupler(MBInterface *impl,
                     MBParallelComm *pc,
                     MBRange &local_elems,
                     int coupler_id,
                     bool init_tree) 
    : mbImpl(impl), myPc(pc), myId(coupler_id)
{
  assert(NULL != impl && NULL != myPc);
  
    // keep track of the local points, at least for now
  myRange = local_elems;
  
    // now initialize the tree
  if (init_tree) initialize_tree();
}

  /* destructor
   */
MBCoupler::~MBCoupler() 
{}


MBErrorCode MBCoupler::locate_points(double *xyz, int num_points,
                                     tuple_list *tl,
                                     bool store_local) 
{
    // for each point, find box(es) containing the point, 
    // appending results to vector

    // allocate tuple_list to hold these points, procs

    // perform scatter/gather, to gather points to source mesh procs

    // find leaf node(s) in local kdtree containing point

    // find natural coordinates of point in element(s) in that leaf

    // for any point/element with nat coords within bounds, store
    // handle/nat coords in vector, and proc/index in outgoing tuple
    // (-1 for index if no elements containing that point)

    // perform scatter/gather to send proc/index tuples back to procs

    // store proc/index tuples locally, and/or pass back to application

    // done
  return MB_SUCCESS;
}

MBErrorCode MBCoupler::interpolate(MBCoupler::Method method,
                                   MBTag tag,
                                   double *interp_vals,
                                   tuple_list *tl,
                                   bool normalize) 
{
  if (LINEAR_FE != method) return MB_FAILURE;
  
    // if no tl, gather procs/indices of locations being interpolated

    // scatter/gather interpolation points

    // perform interpolation on local source mesh

    // normalize interpolation

    // scatter/gather interpolation data

    // done
  return MB_SUCCESS;
}

MBErrorCode MBCoupler::interpolate(MBCoupler::Method method,
                                   const char *tag_name,
                                   double *interp_vals,
                                   tuple_list *tl,
                                   bool normalize) 
{
  MBTag this_tag;
  MBErrorCode result = mbImpl->tag_get_handle(tag_name, this_tag);
  if (MB_SUCCESS != result) return result;
  
  return interpolate(method, this_tag, interp_vals, tl, normalize);
}

MBErrorCode MBCoupler::locate_points(MBRange &ents,
                                     tuple_list *tl,
                                     bool store_local) 
{
    // only do vertices for now
  MBRange verts = ents.subset_by_type(MBVERTEX);
  if (verts.size() != ents.size()) return MB_FAILURE;
  
    // get coordinates
  std::vector<double> coords(3*verts.size());
  MBErrorCode result = mbImpl->get_coords(verts, &coords[0]);
  if (MB_SUCCESS != result) return result;
  
  return locate_points(&coords[0], verts.size(), tl, store_local);
}
  
MBErrorCode MBCoupler::initialize_tree() 
{
    // initialize the tree for the local mesh

    // allocate box extents storage array
  allBoxes = new double[6*myPc->proc_config().proc_size()];
  
    // read local tree minimum/maximum vectors into array starting at
    // allBoxes[6*myPc->proc_rank()], allBoxes[6*myPc->proc_rank()+3]

    // now allgather so each proc has copy of box minima/maxima
  int mpi_err = MPI_Allgather(allBoxes+6*myPc->proc_config().proc_rank(), 6, MPI_DOUBLE,
                              allBoxes, 6, MPI_DOUBLE, myPc->proc_config().proc_comm());
  if (MPI_SUCCESS == mpi_err) return MB_SUCCESS;
  else return MB_FAILURE;
}

  
