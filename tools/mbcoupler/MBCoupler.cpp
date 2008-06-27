#include "MBCoupler.hpp"
#include "MBParallelComm.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBGeomUtil.hpp"
#include "MBElemUtil.hpp"
#include "MBCN.hpp"
#include "iostream"

extern "C" 
{
#include "types.h"
#include "errmem.h"
#include "minmax.h"
#include "sort.h"
#include "tuple_list.h"
#include "transfer.h"
}

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


MBErrorCode MBCoupler::initialize_tree()
{
  
  MBRange local_ents;
  MBAdaptiveKDTree::Settings settings;
  settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION;
  allBoxes.resize(6*myPc->proc_config().proc_size());

    //get entities on the local part
  MBErrorCode result = myPc->get_part_entities(local_ents, 3);
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting entities by dimension" << std::endl;
    return result;
  }

    // build the tree for local processor
  myTree = new MBAdaptiveKDTree(mbImpl);
  result = myTree->build_tree(local_ents, localRoot, &settings);
  if (MB_SUCCESS != result) {
    std::cout << "Problems building tree" << std::endl;
    return result;
  }

    // get the bounding box for local tree
  allBoxes.resize(6*myPc->proc_config().proc_size());
  unsigned int my_rank = myPc->proc_config().proc_rank();
  result = myTree->get_tree_box(localRoot, &allBoxes[6*my_rank], &allBoxes[6*my_rank+3]);

    // now communicate to get all boxes
  int mpi_err = MPI_Allgather(&allBoxes[6*my_rank], 6, MPI_DOUBLE,
                              &allBoxes[0], 6, MPI_DOUBLE, 
                              myPc->proc_config().proc_comm());

#ifndef NDEBUG
  double min[3] = {0,0,0}, max[3] = {0,0,0};
  unsigned int dep;
  myTree->get_info(localRoot, min, max, dep);
  std::cout << "Proc " << my_rank << ": box min/max, tree depth = ("
            << min[0] << "," << min[1] << "," << min[2] << "), ("
            << max[0] << "," << max[1] << "," << max[2] << "), "
            << dep << std::endl;
#endif  

  if (MPI_SUCCESS == mpi_err) return MB_SUCCESS;
  else return MB_FAILURE;
}

MBErrorCode MBCoupler::locate_points(double *xyz, int num_points,
                                     tuple_list *tl,
                                     bool store_local)
{
  assert(tl || store_local);

    // allocate tuple_list to hold point data: (p, i, , xyz), i = point index
  tuple_list target_pts;
  tuple_list_init_max(&target_pts, 2, 0, 0, 3, num_points);

    // initialize source_pts and local_pts
  tuple_list source_pts;
  mappedPts = new tuple_list;
  tuple_list_init_max(&source_pts, 3, 0, 0, 0, target_pts.max); 
  tuple_list_init_max(mappedPts, 0, 0, 1, 3, target_pts.max); 

  mappedPts->n = 0;
  source_pts.n = 0;
  MBErrorCode result;

    // keep track of which points have been located
  std::vector<unsigned char> located_pts(num_points, 0);

    // for each point, find box(es) containing the point,
    // appending results to tuple_list;
    // keep local points separately, in local_pts, which has pairs
    // of <local_index, mapped_index>, where mapped_index is the index
    // of <local_index, mapped_index>, where mapped_index is the index
    // into the mappedPts tuple list

  unsigned int my_rank = myPc->proc_config().proc_rank();
  bool point_located;
  
  for (int i = 0; i < 3*num_points; i+=3) 
  {
      // test point locally first
    result = test_local_box(xyz+i, my_rank, i/3, i/3, point_located);
    if (MB_SUCCESS != result) return result;
    if (point_located) {
      located_pts[i/3] = 0x1;
      continue;
    }

      // if not located locally, test other procs' boxes
    for (unsigned int j = 0; j < myPc->proc_config().proc_size(); j++)
    {
      if (j == my_rank) continue;
      
        // test if point is in proc's box
      if (allBoxes[6*j] <= xyz[i] && xyz[i] <= allBoxes[6*j+3] && 
          allBoxes[6*j+1] <= xyz[i+1] && xyz[i+1] <= allBoxes[6*j+4] && 
          allBoxes[6*j+2] <= xyz[i+2] && xyz[i+2] <= allBoxes[6*j+5])
      {
          // if in this proc's box, will send to proc to test further
          // check size, grow if we're at max
        if (target_pts.n == target_pts.max)
          tuple_list_grow(&target_pts);
  
        target_pts.vi[2*target_pts.n] = j;
        target_pts.vi[2*target_pts.n+1] = i/3;
        target_pts.vr[3*target_pts.n] = xyz[i];
        target_pts.vr[3*target_pts.n+1] = xyz[i+1];
        target_pts.vr[3*target_pts.n+2] = xyz[i+2];
        target_pts.n++;
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

    // test target points against my elements
  for (unsigned i = 0; i < target_pts.n; i++) 
  {
    result = test_local_box(target_pts.vr+3*i, 
                            target_pts.vi[2*i], target_pts.vi[2*i+1], i, 
                            point_located, &source_pts);
    if (MB_SUCCESS != result) return result;
  }

  // no longer need target_pts
  tuple_list_free(&target_pts);

    // send target points back to target procs
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

  // count non-negatives
  int num_pts = 0;
  for (unsigned int i = 0; i < source_pts.n; i++)
    if (-1 != source_pts.vi[3*i+2]) num_pts++;

    // store information about located points
  targetPts = new tuple_list;
  tuple_list *tl_tmp = targetPts;
  if (!store_local) tl_tmp = tl;
  tuple_list_init_max(tl_tmp, 3, 0, 0, 1, num_pts);
  for (unsigned int i = 0; i < source_pts.n; i++) {
    if (-1 != source_pts.vi[3*i+2] && !located_pts[3*i+1]) {
      tl_tmp->vi[3*i] = source_pts.vi[3*i];
      tl_tmp->vi[3*i+1] = source_pts.vi[3*i+1];
      tl_tmp->vi[3*i+2] = source_pts.vi[3*i+2];
      tl_tmp->n++;
    }
  }

  assert(tl_tmp->n + localMappedPts.size()/2 == (unsigned int) num_points);
  
    // no longer need source_pts
  tuple_list_free(&source_pts);

    // copy into tl if passed in and storing locally
  if (tl && store_local) {
    tuple_list_init_max(tl, 3, 0, 0, 1, num_pts);
    memcpy(tl->vi, tl_tmp->vi, 3*tl_tmp->n*sizeof(int));
    tl->n = tl_tmp->n;
  }

    // done
  return MB_SUCCESS;
}

MBErrorCode MBCoupler::test_local_box(double *xyz, 
                                      int from_proc, int remote_index, int index, 
                                      bool &point_located,
                                      tuple_list *tl)
{
  
  std::vector<MBEntityHandle> entities;
  std::vector<MBCartVect> nat_coords;
          
  MBErrorCode result = nat_param(xyz, entities, nat_coords);
  if (MB_SUCCESS != result) return result;

    // grow if we know we'll exceed size
  if (mappedPts->n+entities.size() >= mappedPts->max)
    tuple_list_grow(mappedPts);

  if (entities.empty() && tl) {
    tl->vi[3*index] = from_proc;
    tl->vi[3*index+1] = remote_index;
    tl->vi[3*index+2] = -1;
    point_located = false;
    return MB_SUCCESS;
  }
  
  std::vector<MBEntityHandle>::iterator eit = entities.begin();
  std::vector<MBCartVect>::iterator ncit = nat_coords.begin();
  for (; eit != entities.end(); eit++, ncit++) {
      // store in tuple mappedPts
    mappedPts->vr[3*mappedPts->n] = (*ncit)[0];
    mappedPts->vr[3*mappedPts->n+1] = (*ncit)[1];
    mappedPts->vr[3*mappedPts->n+2] = (*ncit)[2];
    mappedPts->vul[mappedPts->n] = *eit;
    mappedPts->n++;

      // also store local point, mapped point indices
    if (tl) 
    {
      if (tl->n == tl->max) tuple_list_grow(tl);

        // store in tuple source_pts
      tl->vi[3*tl->n] = from_proc;
      tl->vi[3*tl->n+1] = remote_index;
      tl->vi[3*tl->n+2] = mappedPts->n-1;
      tl->n++;
    }
    else {
      localMappedPts.push_back(index);
      localMappedPts.push_back(mappedPts->n-1);
    }
  }

  point_located = true;
  
  return MB_SUCCESS;
}

MBErrorCode MBCoupler::interpolate(MBCoupler::Method method,
                                   std::string &interp_tag,
                                   double *interp_vals,
                                   tuple_list *tl,
                                   bool normalize)
{
  MBTag tag;
  MBErrorCode result = mbImpl->tag_get_handle(interp_tag.c_str(), tag);
  if (MB_SUCCESS != result) return result;
  return interpolate(method, tag, interp_vals, tl, normalize);
}
  
MBErrorCode MBCoupler::interpolate(MBCoupler::Method method,
                                   MBTag tag,
                                   double *interp_vals,
                                   tuple_list *tl,
                                   bool normalize)
{
  if (LINEAR_FE != method) return MB_FAILURE;

  tuple_list *tl_tmp = (tl ? tl : targetPts);

    // remote pts first
  
    // scatter/gather interpolation points
  gs_transfer(1, tl_tmp, 0, myPc->proc_config().crystal_router());

    // perform interpolation on local source mesh; put results into
    // tl_tmp->vr[i]
  MBErrorCode result;
  for (unsigned int i = 0; i < tl_tmp->n; i++) {
    int mindex = tl_tmp->vi[3*i+2];
    result = interp_field_for_hex(mappedPts->vul[mindex],
                                  MBCartVect(mappedPts->vr+3*mindex), 
                                              tag, tl_tmp->vr[i]);
    if (MB_SUCCESS != result) return result;
  }
  
    // scatter/gather interpolation data
  gs_transfer(1, tl_tmp, 0, myPc->proc_config().crystal_router());

  if (!tl) {
      // mapped whole targetPts tuple; put into proper place in interp_vals
    for (unsigned int i = 0; i < tl_tmp->n; i++)
      interp_vals[tl_tmp->vi[3*i+1]] = tl_tmp->vr[i];

      // now do locally-contained pts, since we're mapping everything
    for (std::vector<unsigned int>::iterator vit = localMappedPts.begin();
         vit != localMappedPts.end(); vit += 2) {
      int mindex = *(vit+1);
      result = interp_field_for_hex(mappedPts->vul[mindex],
                                    MBCartVect(mappedPts->vr+3*mindex), 
                                    tag, interp_vals[*vit]);
      if (MB_SUCCESS != result) return result;
    }
  }
  
    // done
  return MB_SUCCESS;
}

MBErrorCode MBCoupler::nat_param(double xyz[3], 
                                 std::vector<MBEntityHandle> &entities, 
                                 std::vector<MBCartVect> &nat_coords)
{
  MBAdaptiveKDTreeIter treeiter;
  MBErrorCode result = myTree->get_tree_iterator(localRoot, treeiter); 
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting iterator" << std::endl;
    return result;
  }

  result = myTree->leaf_containing_point(localRoot, xyz, treeiter);
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting leaf " << std::endl;
    return result;
  }

    // find natural coordinates of point in element(s) in that leaf
  MBCartVect tmp_nat_coords; 
  MBRange range_leaf;
  result = mbImpl->get_entities_by_dimension(treeiter.handle(), 3, range_leaf, false);
  if(result != MB_SUCCESS) std::cout << "Problem getting leaf in a range" << std::endl;

    // loop over the range_leaf 
  for(MBRange::iterator iter = range_leaf.begin(); iter != range_leaf.end(); iter++)
  {
    const MBEntityHandle *connect;
    int num_connect;

      //get connectivity
    result = mbImpl->get_connectivity(*iter, connect, num_connect, true);

      //get coordinates of the vertices
    std::vector<MBCartVect> coords_vert(num_connect);
    std::vector<double> coords(3*num_connect);
    result = mbImpl->get_coords(connect, num_connect, &coords[0]);
	
    for(int j = 0; j < num_connect; j++)
    {
      coords_vert[j][0] = coords[3*j];
      coords_vert[j][1] = coords[3*j+1];
      coords_vert[j][2] = coords[3*j+2];
    }

      //test to find out in which hex the point is
		
      // get natural coordinates
    MBElemUtil::nat_coords_trilinear_hex2(&coords_vert[0], MBCartVect(xyz), 
                                         tmp_nat_coords, 1e-10);
    if (-1.0 <= tmp_nat_coords[0] && tmp_nat_coords[0] <= 1.0 &&
        -1.0 <= tmp_nat_coords[1] && tmp_nat_coords[1] <= 1.0 &&
        -1.0 <= tmp_nat_coords[2] && tmp_nat_coords[2] <= 1.0) {
      entities.push_back(*iter);
      nat_coords.push_back(tmp_nat_coords);
      return MB_SUCCESS;
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBCoupler::interp_field_for_hex(MBEntityHandle elem,
                                            MBCartVect nat_coord, 
                                            MBTag tag,
                                            double &field)
{
    //set the vertices coordinates in the natural system

  const double xi[8] = {-1,1,1,-1,-1,1,1,-1};
  const double etha[8] = {-1,-1,1,1,-1,-1,1,1};
  const double mu[8] = {-1,-1,-1,-1,1,1,1,1};
  double vfields[MB_MAX_SUB_ENTITIES*MB_MAX_SUB_ENTITY_VERTICES];

    // get the tag values at the vertices
  const MBEntityHandle *connect;
  int num_connect;
  MBErrorCode result = mbImpl->get_connectivity(elem, connect, num_connect);
  if (MB_SUCCESS != result) return result;
  result = mbImpl->tag_get_data(tag, connect, num_connect, vfields);
  if (MB_SUCCESS != result) return result;
  
    //function for the interpolation
  field = 0;

    //calculate the field

    // just hexes for now
  assert(num_connect <= 8);
  
  for(int i = 0; i < num_connect; i++)
  {
    field += 0.125 * vfields[i] *
      (1+xi[i]*nat_coord[0]) * (1+etha[i]*nat_coord[1]) * (1+mu[i]*nat_coord[2]);
  }

  return MB_SUCCESS;
}
