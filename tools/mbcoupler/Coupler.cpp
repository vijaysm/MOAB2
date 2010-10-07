#include "Coupler.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/GeomUtil.hpp"
#include "ElemUtil.hpp"
#include "moab/CN.hpp"
#include "iMesh_extensions.h"
#include "iostream"

extern "C" 
{
#include "types.h"
#include "errmem.h"
#include "minmax.h"
#include "sort.h"
#include "tuple_list.h"
#include "crystal.h"
#include "transfer.h"
}

#include "assert.h"

#define ERROR(a) {if (iBase_SUCCESS != err) std::cerr << a << std::endl;}
#define ERRORR(a,b) {if (iBase_SUCCESS != err) {std::cerr << a << std::endl; return b;}}
#define ERRORMPI(a,b) {if (MPI_SUCCESS != err) {std::cerr << a << std::endl; return b;}}

#define MASTER_PROC 0

namespace moab {

bool debug = false;

Coupler::Coupler(Interface *impl,
                     ParallelComm *pc,
                     Range &local_elems,
                     int coupler_id,
                     bool init_tree)
    : mbImpl(impl), myPc(pc), myId(coupler_id), numIts(3)
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
Coupler::~Coupler()
{}


ErrorCode Coupler::initialize_tree()
{
  
  Range local_ents;
  AdaptiveKDTree::Settings settings;
  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;
  allBoxes.resize(6*myPc->proc_config().proc_size());

    //get entities on the local part
  ErrorCode result = myPc->get_part_entities(local_ents, 3);
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting entities by dimension" << std::endl;
    return result;
  }

    // build the tree for local processor
  for (int i = 0; i < numIts; i++) {
    myTree = new AdaptiveKDTree(mbImpl);
    result = myTree->build_tree(local_ents, localRoot, &settings);
    if (MB_SUCCESS != result) {
      std::cout << "Problems building tree";
      if (numIts != i) {
        delete myTree;
        settings.maxEntPerLeaf *= 2;
        std::cout << "; increasing elements/leaf to " 
                  << settings.maxEntPerLeaf << std::endl;;
      }
      else {
        std::cout << "; exiting" << std::endl;
        return result;
      }
    }
  }

    // get the bounding box for local tree
  allBoxes.resize(6*myPc->proc_config().proc_size());
  unsigned int my_rank = myPc->proc_config().proc_rank();
  result = myTree->get_tree_box(localRoot, &allBoxes[6*my_rank], &allBoxes[6*my_rank+3]);

    // now communicate to get all boxes
//   int mpi_err = MPI_Allgather(&allBoxes[6*my_rank], 6, MPI_DOUBLE,
//                               &allBoxes[0], 6, MPI_DOUBLE, 
//                               myPc->proc_config().proc_comm());
    // Change to use "in place" option
  int mpi_err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
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

ErrorCode Coupler::locate_points(double *xyz, int num_points,
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
  ErrorCode result;

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

ErrorCode Coupler::test_local_box(double *xyz, 
                                      int from_proc, int remote_index, int index, 
                                      bool &point_located,
                                      tuple_list *tl)
{
  
  std::vector<EntityHandle> entities;
  std::vector<CartVect> nat_coords;
          
  ErrorCode result = nat_param(xyz, entities, nat_coords);
  if (MB_SUCCESS != result) return result;

    // if we didn't find any ents and we're looking locally, nothing more to do
  if (entities.empty() && !tl) 
    return result;
  
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
  
  std::vector<EntityHandle>::iterator eit = entities.begin();
  std::vector<CartVect>::iterator ncit = nat_coords.begin();
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

ErrorCode Coupler::interpolate(Coupler::Method method,
                                   std::string &interp_tag,
                                   double *interp_vals,
                                   tuple_list *tl,
                                   bool normalize)
{
  Tag tag;
  ErrorCode result = mbImpl->tag_get_handle(interp_tag.c_str(), tag);
  if (MB_SUCCESS != result) return result;
  return interpolate(method, tag, interp_vals, tl, normalize);
}
  
ErrorCode Coupler::interpolate(Coupler::Method method,
                                   Tag tag,
                                   double *interp_vals,
                                   tuple_list *tl,
                                   bool normalize)
{
  if (!((LINEAR_FE == method) || (PLAIN_FE == method)))
    return MB_FAILURE;

  tuple_list *tl_tmp = (tl ? tl : targetPts);

    // remote pts first
  
    // scatter/gather interpolation points
  gs_transfer(1, tl_tmp, 0, myPc->proc_config().crystal_router());

    // perform interpolation on local source mesh; put results into
    // tl_tmp->vr[i]
  ErrorCode result;

  for (unsigned int i = 0; i < tl_tmp->n; i++) {
    int mindex = tl_tmp->vi[3*i+2];

    result = MB_FAILURE;
    if(LINEAR_FE == method){
      result = interp_field_for_hex(mappedPts->vul[mindex],
				    CartVect(mappedPts->vr+3*mindex), 
				    tag, tl_tmp->vr[i]);
    }else if (PLAIN_FE == method){
      result = plain_field_map(mappedPts->vul[mindex],
			       tag, tl_tmp->vr[i]);
    }

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

      result = MB_FAILURE;
      if(LINEAR_FE == method){
	result = interp_field_for_hex(mappedPts->vul[mindex],
				      CartVect(mappedPts->vr+3*mindex), 
				      tag, interp_vals[*vit]);
      }else if (PLAIN_FE == method){
	result = plain_field_map(mappedPts->vul[mindex],
				 tag, interp_vals[*vit]);
      }

      if (MB_SUCCESS != result) return result;

    }
  }
  
    // done
  return MB_SUCCESS;
}

ErrorCode Coupler::nat_param(double xyz[3], 
                                 std::vector<EntityHandle> &entities, 
                                 std::vector<CartVect> &nat_coords)
{
  AdaptiveKDTreeIter treeiter;
  ErrorCode result = myTree->get_tree_iterator(localRoot, treeiter); 
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
  CartVect tmp_nat_coords; 
  Range range_leaf;
  result = mbImpl->get_entities_by_dimension(treeiter.handle(), 3, range_leaf, false);
  if(result != MB_SUCCESS) std::cout << "Problem getting leaf in a range" << std::endl;

    // loop over the range_leaf 
  for(Range::iterator iter = range_leaf.begin(); iter != range_leaf.end(); iter++)
  {
    const EntityHandle *connect;
    int num_connect;

      //get connectivity
    result = mbImpl->get_connectivity(*iter, connect, num_connect, true);

      //get coordinates of the vertices
    std::vector<CartVect> coords_vert(num_connect);
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
    if (ElemUtil::nat_coords_trilinear_hex(&coords_vert[0], CartVect(xyz), 
                                             tmp_nat_coords, 1e-10)) {
      entities.push_back(*iter);
      nat_coords.push_back(tmp_nat_coords);
      return MB_SUCCESS;
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode Coupler::interp_field_for_hex(EntityHandle elem,
                                            CartVect nat_coord, 
                                            Tag tag,
                                            double &field)
{
    //set the vertices coordinates in the natural system

  const double xi[8] = {-1,1,1,-1,-1,1,1,-1};
  const double etha[8] = {-1,-1,1,1,-1,-1,1,1};
  const double mu[8] = {-1,-1,-1,-1,1,1,1,1};
  double vfields[MAX_SUB_ENTITIES*MAX_SUB_ENTITY_VERTICES];

    // get the tag values at the vertices
  const EntityHandle *connect;
  int num_connect;
  ErrorCode result = mbImpl->get_connectivity(elem, connect, num_connect);
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


//Simplest "interpolation" for element-based source fields. Set the value of the field
//at the target point to that of the field in the source element it lies in.
ErrorCode Coupler::plain_field_map(EntityHandle elem,
				       Tag tag,
				       double &field)
{
  double tempField;

  // get the tag values at the vertices
  ErrorCode result = mbImpl->tag_get_data(tag, &elem, 1, &tempField);
  if (MB_SUCCESS != result) return result;

  field = tempField;

  return MB_SUCCESS;
}

// Normalize a field over the subset of entities identified by the tags and values passed
int Coupler::normalize_subset(iBase_EntitySetHandle &root_set,
                              const char            *norm_tag,
                              const char            **tag_names,
                              int                   num_tags,
                              const char            **tag_values,
                              Coupler::IntegType    integ_type,
                              int                   num_integ_pts)
{
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  // Lookup tag handles from tag names
  for (int t = 0; t < num_tags; t++) {
    // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(iMeshInst, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return normalize_subset(root_set, 
                          norm_tag, 
                          &tag_handles[0], 
                          num_tags, 
                          tag_values, 
                          integ_type, 
                          num_integ_pts);
}

int Coupler::normalize_subset(iBase_EntitySetHandle &root_set,
                              const char            *norm_tag,
                              iBase_TagHandle       *tag_handles,
                              int                   num_tags,
                              const char            **tag_values,
                              Coupler::IntegType    integ_type,
                              int                   num_integ_pts)
{
  //  ErrorCode result = MB_SUCCESS;
  int err = iBase_SUCCESS;

//  // Get an iMesh_Instance from MBCoupler::mbImpl.
//  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  // MASTER/SLAVE START #########################################################
  // Broadcast norm_tag, tag_handles, tag_values to procs
  // This is not needed as we are assuming that each proc invokes this function
  // after the data has been distributed.
  // MASTER/SLAVE END   #########################################################

  // SLAVE START ****************************************************************
  // Search for entities based on tag_handles and tag_values

  std::vector< std::vector<iBase_EntitySetHandle> > entitySets;
  std::vector< std::vector<iBase_EntityHandle> >    entityGroups;
  err = get_matching_entities(root_set, tag_handles, tag_values, num_tags, 
                              &entitySets, &entityGroups);
  ERRORR("Failed to get matching entities.", err);

  // Get the integrated field value for each group(vector) of entities.
  // If no entities are in a group then a zero will be put in the list
  // of return values.
  unsigned int numEntGrps = entityGroups.size();
  std::vector<double> integVals(numEntGrps);

  err = get_group_integ_vals(entityGroups, integVals, norm_tag, num_integ_pts, integ_type);
  ERRORR("Failed to get integrated field values for groups in mesh.", err);
  // SLAVE END   ****************************************************************

  // SLAVE/MASTER START #########################################################
  // Send list of integrated values back to master proc.  The ordering of the 
  // values will match the ordering of the entity groups (i.e. vector of vectors)
  // sent from master to slaves earlier.  The values for each entity group will
  // be summed during the transfer.
  std::vector<double> sumIntegVals(numEntGrps);
  err = MPI_Reduce(&integVals[0], &sumIntegVals[0], numEntGrps, MPI_DOUBLE, 
                   MPI_SUM, MASTER_PROC, myPc->proc_config().proc_comm());
  ERRORMPI("Transfer and reduction of integrated values failed.", err);
  // SLAVE/MASTER END   #########################################################

  // MASTER START ***************************************************************
  // Calculate the normalization factor for each group by taking the
  // inverse of each integrated field value.  Put the normalization factor 
  // for each group back into the list in the same order.

  for (unsigned int i = 0; i < numEntGrps; i++) {
    double val = sumIntegVals[i];
    if (val != 0) sumIntegVals[i] = 1.0/val;
  }
  // MASTER END   ***************************************************************

  // MASTER/SLAVE START #########################################################
  // Broadcast the normalization factors to the procs.

  err = MPI_Bcast(&sumIntegVals[0], numEntGrps, MPI_DOUBLE, MASTER_PROC, 
                  myPc->proc_config().proc_comm());
  ERRORMPI("", err);
  // MASTER/SLAVE END   #########################################################

  // SLAVE START ****************************************************************
  // Save the normalization factors to a new tag with name of norm_tag's value
  // and the string "_normF" appended.  This new tag will be created on the entity
  // set that contains all of the entities from a group.

  err = apply_group_norm_factor(entitySets, sumIntegVals, norm_tag, integ_type);
  ERRORR("Failed to set the normalization factor for groups in mesh.", err);
  // SLAVE END   ****************************************************************

  return err;
}

// Functions supporting the subset normalization function

// Retrieve groups of entities matching tags and values if present
int Coupler::get_matching_entities(iBase_EntitySetHandle                             root_set,
                                   const char                                        **tag_names,
                                   const char                                        **tag_values,
                                   int                                               num_tags,
                                   std::vector< std::vector<iBase_EntitySetHandle> > *entity_sets,
                                   std::vector< std::vector<iBase_EntityHandle> >    *entity_groups)
{
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  for (int t = 0; t < num_tags; t++) {
    // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(iMeshInst, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return get_matching_entities(root_set, &tag_handles[0], tag_values, num_tags,
                               entity_sets, entity_groups);
}

// Retrieve groups of entities matching tags and values if present
int Coupler::get_matching_entities(iBase_EntitySetHandle                          root_set,
                                   iBase_TagHandle                                *tag_handles,
                                   const char                                     **tag_values,
                                   int                                            num_tags,
                                   std::vector< std::vector<iBase_EntitySetHandle> > *entity_sets,
                                   std::vector< std::vector<iBase_EntityHandle> > *entity_groups)
{                                        

  // SLAVE START ****************************************************************
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;
  int entSetsSize;
  int entSetsAlloc=0;
  iBase_EntitySetHandle *entSets = NULL;  // free at end

  // Get Entity Sets that match the tags and values.
  iMesh_getEntSetsByTagsRec(iMeshInst, root_set, tag_handles, 
                            tag_values, num_tags, 0,
                            &entSets, &entSetsAlloc, &entSetsSize, &err);
  ERRORR("iMesh_getEntSetsByTagsRec failed.", err);

  tuple_list *tagList = NULL;
  err = create_tuples(entSets, entSetsSize, tag_handles, num_tags, &tagList);
  ERRORR("Failed to create tuples from entity sets.", err);

  // Free up array memory from iMesh call
  free(entSets);
  entSets = NULL;
  entSetsAlloc = 0;
  entSetsSize = 0;
  // SLAVE END   ****************************************************************

  // SLAVE/MASTER START #########################################################
  // Send tuple list back to master proc for consolidation
  
  // pack the tuple_list in a buffer.
  

  // Send all buffers to the master proc for consolidation
//  MPI_Gatherv();

  // unpack the tuple_list from the buffer.

  // SLAVE/MASTER END   #########################################################

  // MASTER START ***************************************************************
  // Master proc receives all lists and consolidates to one set with no duplicates.
//   tuple_list *tuples = tagList;
//   tuple_list *consTuples = NULL;
  // TEMPORARY - assignment for serial version
  tuple_list *consTuples = tagList;
//   int numTuples = 1;
//   err = consolidate_tuples(&tuples, numTuples, &consTuples);
//   ERRORR("Failed to consolidate tuples.", err);
  // MASTER END   ***************************************************************

  // MASTER/SLAVE START #########################################################
  // Broadcast condensed tuple list back to all procs.
  // ***TODO***
  // MASTER/SLAVE END   #########################################################

  // SLAVE START ****************************************************************
  // Loop over the tuple list getting the entities with the tags in the tuple_list entry
  for (unsigned int i = 0; i < consTuples->n; i++) {
    // Get Entity Sets that match the tags and values.

    // Convert the data in the tuple_list to an array of pointers to the data
    // in the tuple_list as that is what the iMesh API call is expecting.
    int **vals = new int*[consTuples->mi];
    for (unsigned int j = 0; j < consTuples->mi; j++)
      vals[j] = &(consTuples->vi[(i*consTuples->mi) + j]);

    iMesh_getEntSetsByTagsRec(iMeshInst, root_set, tag_handles, 
                              (const char * const *) vals,
                              consTuples->mi, 0,
                              &entSets, &entSetsAlloc, &entSetsSize, &err);
    ERRORR("iMesh_getEntSetsByTagsRec failed.", err);
    if (debug) std::cout << "entSetsSize=" << entSetsSize << std::endl;

    // Free up the array of pointers
    free(vals);

    // Loop over the entity sets and then free the memory for entSets.
    std::vector<iBase_EntitySetHandle> entSetHandles;
    std::vector<iBase_EntityHandle> entHandles;
    for (int j = 0; j < entSetsSize; j++) {
      // Save the entity set
      entSetHandles.push_back(entSets[j]);

      // Get all entities for the entity set
      iBase_EntityHandle *ents = NULL;
      int entsAlloc = 0;
      int entsSize = 0;

      iMesh_getEntities(iMeshInst, entSets[j], iBase_ALL_TYPES, iMesh_ALL_TOPOLOGIES,
                        &ents, &entsAlloc, &entsSize, &err);
      ERRORR("iMesh_getEntities failed.", err);
      if (debug) std::cout << "entsSize=" << entsSize << std::endl;

      // Save all of the entities from the entity set and free the memory for ents.
      for (int k = 0; k < entsSize; k++) {
        entHandles.push_back(ents[k]);
      }
      free(ents);
      if (debug) std::cout << "entHandles.size=" << entHandles.size() << std::endl;
    }

    // Free the entity set list for next tuple iteration.
    free(entSets);
    entSets = NULL;
    entSetsAlloc = 0;
    entSetsSize = 0;

    // Push entSetHandles onto entity_sets, entHandles onto entity_groups
    // and clear both entSetHandles and entHandles.
    entity_sets->push_back(entSetHandles);
    entSetHandles.clear();
    entity_groups->push_back(entHandles);
    entHandles.clear();
    if (debug) std::cout << "entity_sets->size=" << entity_sets->size() 
                         << ", entity_groups->size=" << entity_groups->size() << std::endl;
  }
  // SLAVE END   ****************************************************************

  return err;
}


// Return a tuple_list containing  tag values for each Entity Set
// The tuple_list will have a column for each tag and a row for each
// Entity Set.  It is assumed all of the tags are integer tags.
int Coupler::create_tuples(iBase_EntitySetHandle *ent_sets,
                           int                   num_sets, 
                           const char            **tag_names,
                           int                   num_tags,
                           tuple_list            **tuple_list)
{
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);
  int err;
  std::vector<iBase_TagHandle> tag_handles;
  
  for (int t = 0; t < num_tags; t++) {
    // get tag handle & size
    iBase_TagHandle th;
    iMesh_getTagHandle(iMeshInst, tag_names[t], &th, &err, strlen(tag_names[t]));
    ERRORR("Failed to get tag handle.", err);
    tag_handles.push_back(th);
  }

  return create_tuples(ent_sets, num_sets, &tag_handles[0], num_tags, tuple_list);
}

// Return a tuple_list containing  tag values for each Entity Set
// The tuple_list will have a column for each tag and a row for each
// Entity Set.  It is assumed all of the tags are integer tags.
int Coupler::create_tuples(iBase_EntitySetHandle *ent_sets, 
                           int                   num_sets, 
                           iBase_TagHandle       *tag_handles,
                           int                   num_tags,
                           tuple_list            **tuples)
{
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;

  // ASSUMPTION: All tags are of type integer.  This may need to be expanded in future.

  // Allocate a tuple_list for the number of entity sets passed in
  tuple_list *tag_tuples = new tuple_list;
  tuple_list_init_max(tag_tuples, num_tags, 0, 0, 0, num_sets);
  if (tag_tuples->mi == 0)
    ERRORR("Failed to initialize tuple_list.", iBase_FAILURE);

  // Loop over the filtered entity sets retrieving each matching tag value one by one.
  int val;
  for (int i = 0; i < num_sets; i++) {
    for (int j = 0; j < num_tags; j++) {
      iMesh_getEntSetIntData(iMeshInst, ent_sets[i], tag_handles[j], &val, &err);
      ERRORR("Failed to get integer tag data.", err);
      tag_tuples->vi[i*tag_tuples->mi + j] = val;
    }

    // If we get here there was no error so increment n in the tuple_list
    tag_tuples->n++;
  }

  *tuples = tag_tuples;

  return err;
}

// Consolidate tuple_lists into one list with no duplicates
int Coupler::consolidate_tuples(tuple_list **all_tuples, 
                                int        num_tuples,
                                tuple_list **unique_tuples)
{
  int err = iBase_SUCCESS;

  const unsigned intSize = sizeof(sint);
  const unsigned numTags = all_tuples[0]->mi;
  const unsigned intWidth = numTags * intSize;

  int totalRcvTuples = 0;
  int offset = 0, copysz = 0;

  // Get the total size of all of the tuple_lists in all_tuples.
  for (int i = 0; i < num_tuples; i++) {
    totalRcvTuples += all_tuples[i]->n;
  }

  // Copy the tuple_lists into a single tuple_list.
  tuple_list *newTupleList = new tuple_list;
  tuple_list_init_max(newTupleList, numTags, 0, 0, 0, totalRcvTuples);
  for (int i = 0; i < num_tuples; i++) {
    copysz = all_tuples[i]->n * intWidth;
    memcpy(newTupleList->vi+offset, all_tuples[i]->vi, copysz);
    offset = offset + (all_tuples[i]->n * all_tuples[i]->mi);
    newTupleList->n = newTupleList->n + all_tuples[i]->n;
  }

  // Sort the new tuple_list.  Use a radix type sort, starting with the last (or least significant)
  // tag column in the vi array and working towards the first (or most significant) tag column.
  buffer sort_buffer;
  buffer_init(&sort_buffer, 2 * totalRcvTuples * intWidth);
  for (int i = numTags - 1; i >= 0; i--) {
    tuple_list_sort(newTupleList, i, &sort_buffer);
  }

  // Cycle through the sorted list eliminating duplicates.
  // Keep counters to the current end of the tuple_list (w/out dups) and the last tuple examined.
  unsigned int endIdx = 0, lastIdx = 1;
  while (lastIdx < newTupleList->n) {
    if (memcmp(newTupleList->vi+(endIdx*numTags), newTupleList->vi+(lastIdx*numTags), intWidth) == 0) {
      // Values equal - skip
      lastIdx += 1;
    }
    else {
      // Values different - copy
      // Move up the end index
      endIdx += 1;
      memcpy(newTupleList->vi+(endIdx*numTags), newTupleList->vi+(lastIdx*numTags), intWidth);
      lastIdx += 1;
    }
  }
  // Update the count in newTupleList
  newTupleList->n = endIdx + 1;

  // Resize the tuple_list
  tuple_list_resize(newTupleList, newTupleList->n);

  // Set the output parameter
  *unique_tuples = newTupleList;

  return err;
}

// Calculate integrated field values for groups of entities
int Coupler::get_group_integ_vals(std::vector< std::vector<iBase_EntityHandle> > &groups,
                                  std::vector<double> &integ_vals,
                                  const char *norm_tag,
                                  int num_integ_vals,
                                  Coupler::IntegType integ_type)
{
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;

  std::vector< std::vector<iBase_EntityHandle> >::iterator iter_i;
  std::vector<iBase_EntityHandle>::iterator iter_j;
  double grpIntgrVal, intgrVal;

  // Get the tag handle for norm_tag
  iBase_TagHandle norm_hdl;
  iMesh_getTagHandle(iMeshInst, norm_tag, &norm_hdl, &err, strlen(norm_tag));
  ERRORR("Failed to get norm_tag handle.", err);

  int norm_type;
  iMesh_getTagType(iMeshInst, norm_hdl, &norm_type, &err);
  ERRORR("Failed to get norm_tag type.", err);

  // Check size of integ_vals vector
  if (integ_vals.size() != groups.size())
    integ_vals.resize(groups.size());

  // Loop over the groups(vectors) of entities
  unsigned int i;
  for (i = 0, iter_i = groups.begin(); iter_i != groups.end(); i++, iter_i++) {
    grpIntgrVal = 0;

    // Loop over the all the entities in the group, integrating 
    // the field_fn over the entity in iter_j
    for (iter_j = (*iter_i).begin(); iter_j != (*iter_i).end(); iter_j++) {
      // Check that the entity in iter_j is of the same dimension as the 
      // integ_type we are performing
      int j_type;
      iMesh_getEntType(iMeshInst, (*iter_j), &j_type, &err);
      ERRORR("Failed to get entity type.", err);
      if (((integ_type == VOLUME) && (j_type != iBase_REGION)) ||
          ((integ_type == AREA)   && (j_type != iBase_FACE)))
        continue;
      
      intgrVal = 0;

      // Currently assumes we have a hexahedral element.  Will need code
      // to choose between a hex and a tet element.

      // Get coordinates of all corner vertices (in normal order) and
      // put in array of CartVec.
      CartVect corners[8];

      // Retrieve the vertices from the element
      iBase_EntityHandle *verts = NULL;
      int vertsAlloc = 0;
      int vertsSize = 0;

      iMesh_getEntAdj(iMeshInst, (*iter_j), iBase_VERTEX, &verts, &vertsAlloc, &vertsSize, &err);
      ERRORR("Failed to get vertices from entity.", err);

      if (vertsSize < 8)
        ERRORR("Failed to get 8 vertices.", iBase_FAILURE);

      // Put the vertices into a CartVect array
      double x[3];
      double vfield[8];
      for (int i = 0; i < vertsSize; i++) {
        iMesh_getVtxCoord(iMeshInst, verts[i], &x[0], &x[1], &x[2], &err);
        corners[i] = x;
        if (norm_type == iBase_INTEGER) {
          int data;
          iMesh_getIntData(iMeshInst, verts[i], norm_hdl, &data, &err);
          ERRORR("Failed to get vertex int data for norm_tag.", err);
          vfield[i] = (double) data;
        } else if (norm_type == iBase_DOUBLE) {
          double data;
          iMesh_getDblData(iMeshInst, verts[i], norm_hdl, &data, &err);
          ERRORR("Failed to get vertex double data for norm_tag.", err);
          vfield[i] = data;
        } else {
          err = iBase_FAILURE;
          ERRORR("Bad type for norm_tag.", err);
        }
      }

      // Perform the integration
      if (!ElemUtil::integrate_trilinear_hex(corners, vfield, intgrVal, num_integ_vals))
          ERRORR("Bad parameter to ElemUtil::integrate_trilinear_hex.", iBase_FAILURE);

      // Combine the result with those of the group
      grpIntgrVal += intgrVal;
    }

    // Set the group integrated value in the vector
    integ_vals[i] = grpIntgrVal;
  }

  return err;
}

// Apply a normalization factor to group of entities
int Coupler::apply_group_norm_factor(std::vector< std::vector<iBase_EntitySetHandle> > &entity_sets,
                                     std::vector<double> &norm_factors, 
                                     const char *norm_tag,
                                     Coupler::IntegType integ_type)
{
  // Get an iMesh_Instance from MBCoupler::mbImpl.
  iMesh_Instance iMeshInst = reinterpret_cast<iMesh_Instance>(mbImpl);

  int err = iBase_SUCCESS;

  // Construct the new tag for the normalization factor from the norm_tag name
  // and "_normf".
  int norm_tag_len = strlen(norm_tag);
  const char* normf_appd = "_normf";
  int normf_appd_len = strlen(normf_appd);

  char* normf_tag = (char*) malloc(norm_tag_len + normf_appd_len + 1);
  char* tmp_ptr = normf_tag;

  memcpy(tmp_ptr, norm_tag, norm_tag_len);
  tmp_ptr += norm_tag_len;
  memcpy(tmp_ptr, normf_appd, normf_appd_len);
  tmp_ptr += normf_appd_len;
  *tmp_ptr = '\0';

  iBase_TagHandle normf_hdl;
  iMesh_createTag(iMeshInst, normf_tag, 1, iBase_DOUBLE, &normf_hdl, &err, strlen(normf_tag));
  ERRORR("Failed to create normalization factor tag.", err);

  std::vector< std::vector<iBase_EntitySetHandle> >::iterator iter_i;
  std::vector<iBase_EntitySetHandle>::iterator iter_j;
  std::vector<double>::iterator iter_f;
  double grpNormFactor = 0.0;

  // Loop over the entity sets
  for (iter_i = entity_sets.begin(), iter_f = norm_factors.begin(); 
       (iter_i != entity_sets.end()) && (iter_f != norm_factors.end()); 
       iter_i++, iter_f++) {
    grpNormFactor = *iter_f;

    // Loop over the all the entity sets in iter_i and set the 
    // new normf_tag with the norm factor value on each
    for (iter_j = (*iter_i).begin(); iter_j != (*iter_i).end(); iter_j++) {

      iMesh_setEntSetDblData(iMeshInst, (*iter_j), normf_hdl, grpNormFactor, &err);
      ERRORR("Failed to set normalization factor on entity set.", err);
    }
  }

  return err;
}

} // namespace_moab
