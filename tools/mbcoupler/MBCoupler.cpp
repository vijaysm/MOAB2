#include "MBCoupler.hpp"
#include "MBParallelComm.hpp"
#include "types.h"
#include "errmem.h"
#include "minmax.h"
#include "sort.h"
#include "tuple_list.h"
#include "transfer.h"

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

    // initialize tuple lists to indicate not initialized
  mappedPts = NULL;
  targetPts = NULL;
}

  /* destructor
   */
MBCoupler::~MBCoupler() 
{}


MBErrorCode MBCoupler::locate_points(double *xyz, int num_points,
                                     tuple_list *tl,
                                     bool store_local) 
{
  assert(tl || store_local);
  
    // allocate tuple_list to hold point data: (p, i, , xyz), i = point index
  tuple_list target_pts;
  tuple_list_init_max(&target_pts, 2, 0, 0, 3, 3*num_points);

    // for each point, find box(es) containing the point, 
    // appending results to tuple_list; 
    // keep local points separately, in local_pts, which has pairs
    // of <local_index, mapped_index>, where mapped_index is the index
    // into the mappedPts tuple list
  int threen = 0;
  for (int i = 0; i < 3*num_points; i+=3) {
    for (;;/*  <marco - loop over boxes> */ ) {
      if (/* marco - test point i/3 in box j */ false) {
          // check size, grow if we're at max
        if (target_pts.n == target_pts.max) 
          tuple_list_grow(&target_pts);

        target_pts.vi[2*target_pts.n] = -1 /* <marco - proc rank j> */;
        target_pts.vi[2*target_pts.n+1] = i/3;
        target_pts.vr[threen] = xyz[i];
        target_pts.vr[threen] = xyz[i+1];
        target_pts.vr[threen] = xyz[i+2];
        target_pts.n++;
        threen += 3;
      }
    }
  }

    // perform scatter/gather, to gather points to source mesh procs
  gs_transfer(1, &target_pts, 0, myPc->proc_config().crystal_router());

    // after scatter/gather: 
    // target_pts.n = # points local proc has to map
    // target_pts.vi[2*i] = proc sending point i
    // target_pts.vi[2*i+1] = index of point i on sending proc
    // target_pts.vr[3*i..3*i+2] = xyz of point i
    // 
    // Mapping builds the tuple list:
    // source_pts.n = target_pts.n
    // source_pts.vi[3*i] = target_pts.vi[2*i] = sending proc
    // source_pts.vi[3*i+1] = index of point i on sending proc
    // source_pts.vi[3*i+2] = index of mapped point (-1 if not mapped)
    //
    // Also, mapping builds local tuple_list mappedPts:
    // mappedPts->n = # mapped points
    // mappedPts->vul[i] = local handle of mapped entity
    // mappedPts->vr[3*i..3*i+2] = natural coordinates in mapped entity

    // initialize source_pts and local_pts
  tuple_list source_pts;
  mappedPts = new tuple_list;
  tuple_list_init_max(&source_pts, 3, 0, 0, 0, target_pts.n);
  tuple_list_init_max(mappedPts, 0, 0, 1, 3, target_pts.n);

  for (unsigned i = 0; i < target_pts.n; i++) {
      // find leaf node(s) in local kdtree containing point i
      // <marco...>

      // find natural coordinates of point in element(s) in that leaf
      // <marco...>

      // for any point/element with nat coords within bounds, store
      // handle/nat coords in mappedPts, and proc/index in outgoing tuple
      // (-1 for index if no elements containing that point)
      // <marco...>
    
  }

    // no longer need target_pts
  tuple_list_free(&target_pts);

    // perform scatter/gather to send proc/index tuples back to procs
  gs_transfer(1, &source_pts, 0, myPc->proc_config().crystal_router());

    // store proc/index tuples in targetPts, and/or pass back to application;
    // the tuple this gets stored to looks like:
    // tl.n = # mapped points
    // tl.vi[3*i] = remote proc mapping point
    // tl.vi[3*i+1] = local index of mapped point
    // tl.vi[3*i+2] = remote index of mapped point
    //
    // Local index is mapped into either myRange, holding the handles of
    // local mapped entities, or myXyz, holding locations of mapped pts

    // go through and count non-negatives
  int num_pts = 0;
  for (unsigned int i = 0; i < source_pts.n; i++) 
    if (-1 != source_pts.vi[3*i+2]) num_pts++;

  targetPts = new tuple_list;
  tuple_list *tl_tmp = targetPts;
  if (!store_local) tl_tmp = tl;
  tuple_list_init_max(tl_tmp, 3, 0, 0, 1, num_pts);
  tl_tmp->n = num_pts;
  for (unsigned int i = 0; i < source_pts.n; i++) {
    if (-1 != source_pts.vi[3*i+2]) {
      tl_tmp->vi[3*i] = source_pts.vi[3*i];
      tl_tmp->vi[3*i+1] = source_pts.vi[3*i+1];
      tl_tmp->vi[3*i+2] = source_pts.vi[3*i+2];
    }
  }
  
    // no longer need source_pts
  tuple_list_free(&source_pts);

    // copy into tl if passed in and storing locally
  if (tl && store_local) {
    tuple_list_init_max(tl, 3, 0, 0, 1, num_pts);
    memcpy(tl->vi, tl_tmp->vi, 3*num_pts*sizeof(int));
    tl->n = tl_tmp->n;
  }
  
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

  tuple_list *tl_tmp = (tl ? tl : targetPts);
  
    // scatter/gather interpolation points
  gs_transfer(1, tl_tmp, 0, myPc->proc_config().crystal_router());

    // perform interpolation on local source mesh; put results into
    // tl_tmp->vr[i]

    // interpolate locally mapped points too, putting results directly 
    // into interp_vals

    // normalize interpolation

    // scatter/gather interpolation data
  gs_transfer(1, tl_tmp, 0, myPc->proc_config().crystal_router());
  for (unsigned int i = 0; i < tl_tmp->n; i++)
    interp_vals[tl_tmp->vi[3*i+1]] = tl_tmp->vr[i];

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

  
