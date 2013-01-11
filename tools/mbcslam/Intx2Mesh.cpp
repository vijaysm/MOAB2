/*
 * Intx2Mesh.cpp
 *
 *  Created on: Oct 2, 2012
 */

#include "Intx2Mesh.hpp"
//#include "Coupler.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "MBParallelConventions.h"
#include <queue>

namespace moab {

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {std::cout << str << "\n"; return rval;}

#define ERRORV(rval, str) \
    if (MB_SUCCESS != rval) {std::cout << str << "\n"; return ;}

Intx2Mesh::Intx2Mesh(Interface * mbimpl):mb(mbimpl), parcomm(NULL), myTree(NULL)
{
  dbg_1=0;
  box_error=0;
}

Intx2Mesh::~Intx2Mesh()
{
  // TODO Auto-generated destructor stub
}
void Intx2Mesh::createTags()
{
  unsigned char def_data_bit = 0; // unused by default
  ErrorCode rval = mb->tag_get_handle("blueFlag", 1, MB_TYPE_BIT, BlueFlagTag,
      MB_TAG_CREAT, &def_data_bit);
  if (MB_SUCCESS != rval)
    return;
  // maybe the red tag is better to be deleted every time, and recreated;
  // or is it easy to set all values to something again? like 0?
  rval = mb->tag_get_handle("redFlag", 1, MB_TYPE_BIT, RedFlagTag, MB_TAG_CREAT,
      &def_data_bit);
  ERRORV(rval, "can't get red flag tag");

  // assume that the edges are on the red triangles
  Range redQuads;
  //Range redEdges;
  rval = mb->get_entities_by_type(mbs2, type, redQuads, false);
  ERRORV(rval, "can't get ents by type");

  // create red edges if they do not exist yet; so when they are looked upon, they are found
  // this is the only call that is potentially NlogN, in the whole method
  rval = mb->get_adjacencies(redQuads, 1, true, RedEdges, Interface::UNION);
  ERRORV(rval, "can't get adjacent red edges");

  // now, create a map from each edge to a list of potential new nodes on a red edge
  // this memory has to be cleaned up
  // change it to a vector, and use the index in range of red edges
  int indx = 0;
  extraNodesVec.reserve(RedEdges.size());
  for (Range::iterator eit = RedEdges.begin(); eit != RedEdges.end();
      eit++, indx++)
  {
    //EntityHandle edge = *eit;
    //extraNodesMap[edge] = new std::vector<EntityHandle>;
    std::vector<EntityHandle> * nv = new std::vector<EntityHandle>;
    extraNodesVec.push_back(nv);
  }

  int defaultInt = 0;

  rval = mb->tag_get_handle("Positive", 1, MB_TYPE_INTEGER, redParentTag,
      MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);
  ERRORV(rval, "can't create positive tag");

  rval = mb->tag_get_handle("Negative", 1, MB_TYPE_INTEGER, blueParentTag,
      MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);
  ERRORV(rval, "can't create negative tag");

  rval = mb->tag_get_handle("Counting", 1, MB_TYPE_INTEGER, countTag,
        MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);
  ERRORV(rval, "can't create Counting tag");

  return;
}


// specify also desired set; we are interested only in neighbors in the set!
// we should always get manifold mesh, each edge is adjacent to 2 quads
ErrorCode Intx2Mesh::GetOrderedNeighbors(EntityHandle set, EntityHandle quad,
    EntityHandle neighbors[4])
{
  int nsides = 3;
  if (type == MBQUAD)
    nsides = 4;
  // will get the 4 ordered neighbors;
  // first quad is for nodes 0, 1, second to 1, 2, third to 2, 3,
  int nnodes;
  const EntityHandle * conn4;
  ErrorCode rval = mb->get_connectivity(quad, conn4, nnodes);
  ERRORR(rval, "can't get connectivity on an element");
  for (int i = 0; i < nsides; i++)
  {
    EntityHandle v[2];
    v[0] = conn4[i];
    v[1] = conn4[(i + 1) % nsides];
    // get quads adjacent to vertices
    std::vector<EntityHandle> quads;
    std::vector<EntityHandle> quadsInSet;
    rval = mb->get_adjacencies(v, 2, 2, false, quads, Interface::INTERSECT);
    ERRORR(rval, "can't get adjacencies on 2 nodes");
    size_t siz = quads.size();
    for (size_t j = 0; j < siz; j++)
      if (mb->contains_entities(set, &(quads[j]), 1))
        quadsInSet.push_back(quads[j]);
    siz = quadsInSet.size();

    if (siz > 2)
    {
      std::cout << "non manifold mesh, error"
          << mb->list_entities(&(quadsInSet[0]), quadsInSet.size()) << "\n";
      return MB_FAILURE; // non-manifold
    }
    if (siz == 1)
    {
      // it must be the border,
      neighbors[i] = 0; // we are guaranteed that ids are !=0; this is marking a border
      continue;
    }
    // here siz ==2, it is either the first or second
    if (quad == quads[0])
      neighbors[i] = quads[1];
    else
      neighbors[i] = quads[0];
  }
  return MB_SUCCESS;
}
// main interface; this will do the advancing front trick
// some are triangles, some are quads
ErrorCode Intx2Mesh::intersect_meshes(EntityHandle mbset1, EntityHandle mbset2,
     EntityHandle & outputSet)
{

  int nsides=3;
  if (type == MBQUAD)
    nsides=4;
  ErrorCode rval;
  mbs1 = mbset1; // is this the arrival or departure? I don't know yet
  mbs2 = mbset2;
  outSet = outputSet;

  // really, should be something from t1 and t2; blue is 1 (lagrange), red is 2 (euler)
  createTags(); //
  EntityHandle startBlue, startRed;

  Range rs1;
  Range rs2;
  mb->get_entities_by_type(mbs1, type, rs1);
  mb->get_entities_by_type(mbs2, type, rs2);
  for (Range::iterator it = rs1.begin(); it != rs1.end(); it++)
  {
    startBlue = *it;
    int found = 0;
    for (Range::iterator it2 = rs2.begin(); it2 != rs2.end() && !found; it2++)
    {
      startRed = *it2;
      double area = 0;
      // if area is > 0 , we have intersections
      double P[48]; // max 8 intx points + 8 more in the polygon
      // the red quad is convex, always, while the blue can be concave
      int nP = 0;
      int nb[4], nr[4]; // sides 3 or 4?
      computeIntersectionBetweenRedAndBlue(startRed, startBlue, P, nP, area, nb, nr);
      if (area > 0)
      {
        found = 1;
        break; // found 2 elements that intersect; these will be the seeds
      }
    }
    if (found)
      break;
  }

  std::queue<EntityHandle> blueQueue; // these are corresponding to Ta,
  blueQueue.push(startBlue);
  std::queue<EntityHandle> redQueue;
  redQueue.push(startRed);

  Range toResetReds; // will be used to reset red flags for every blue quad
  // processed

  unsigned char used = 1;
  unsigned char unused = 0; // for red flags
  // mark the start blue quad as used, so it will not come back again
  mb->tag_set_data(BlueFlagTag, &startBlue, 1, &used);
  while (!blueQueue.empty())
  {
    // flags for the side : 0 means a red quad not found on side
    // a paired red not found yet for the neighbors of blue
    EntityHandle n[4] = { EntityHandle(0) };

    EntityHandle currentBlue = blueQueue.front();
    blueQueue.pop();
    //        for (k=0; k<m_numPos; k++)
    //          redFlag[k] = 0;
    //        redFlag[m_numPos] = 1; // to guard for the boundary
    // all reds that were tagged, are now cleared
    for (Range::iterator itr = toResetReds.begin(); itr != toResetReds.end();
        itr++)
    {
      EntityHandle ttt = *itr;
      rval = mb->tag_set_data(RedFlagTag, &ttt, 1, &unused);
      ERRORR(rval, "can't set red tag");
    }
    //rval = mb2->tag_set_data(RedFlagTag, toResetReds, &unused);
    if (dbg_1)
    {
      std::cout << "reset reds: ";
      for (Range::iterator itr = toResetReds.begin(); itr != toResetReds.end();
          itr++)
        std::cout << mb->id_from_handle(*itr) << " ";
      std::cout << std::endl;
    }
    EntityHandle currentRed = redQueue.front(); // where do we check for redQueue????
    // red and blue queues are parallel
    redQueue.pop(); // mark the current red
    //redFlag[currentRed] = 1; //
    toResetReds.clear(); // empty the range of used reds, will have to be set unused again,
    // at the end of blue triangle processing
    toResetReds.insert(currentRed);
    rval = mb->tag_set_data(RedFlagTag, &currentRed, 1, &used);
    ERRORR(rval, "can't set red tag");
    //mb2->set_tag_data
    std::queue<EntityHandle> localRed;
    localRed.push(currentRed);
    while (!localRed.empty())
    {
      //
      EntityHandle redT = localRed.front();
      localRed.pop();
      double P[48], area; // area is in 2d, points are in 3d (on a sphere), back-projected
      int nP = 0; // intersection points (could include the vertices of initial quads)
      int nb[4] = { 0, 0, 0, 0 }; // means no intersection on the side (markers)
      int nr[4] = { 0, 0, 0, 0 }; // means no intersection on the side (markers)
      // nc [j] = 1 means that the side j (from j to j+1) of blue quad intersects the
      // red quad.  A potential next quad is the red quad that is adjacent to this side
      computeIntersectionBetweenRedAndBlue(/* red */redT, currentBlue, P, nP,
          area, nb, nr);
      if (nP > 0)
      {
        // intersection found: output P and original triangles if nP > 2

        EntityHandle neighbors[4];
        rval = GetOrderedNeighbors(mbs2, redT, neighbors);
        if (rval != MB_SUCCESS)
        {
          std::cout << " can't get the neighbors for red element "
              << mb->id_from_handle(redT);
          return MB_FAILURE;
        }

        // add neighbors to the localRed queue, if they are not marked
        for (int nn = 0; nn < nsides; nn++)
        {
          EntityHandle neighbor = neighbors[nn];
          if (neighbor > 0 && nr[nn]>0) // advance across red boundary n
          {
            //n[nn] = redT; // start from 0!!
            unsigned char status = 0;
            mb->tag_get_data(RedFlagTag, &neighbor, 1, &status);
            if (status == 0)
            {
              localRed.push(neighbor);
              if (dbg_1)
              {
                std::cout << " local red elem " << mb->id_from_handle(neighbor)
                    << " for blue:" << mb->id_from_handle(currentBlue) << "\n"
                    << mb->list_entities(&neighbor, 1) << "\n";
              }
              rval = mb->tag_set_data(RedFlagTag, &neighbor, 1, &used);
              //redFlag[neighbor] = 1; // flag it to not be added anymore
              toResetReds.insert(neighbor); // this is used to reset the red flag
            }
          }
          // n(find(nc>0))=ac;        % ac is starting candidate for neighbor
          if (nb[nn] > 0)
            n[nn] = redT;

        }
        if (nP > 1) // this will also construct triangles/polygons in the new mesh, if needed
          findNodes(redT, currentBlue, P, nP);
      }
      else if (dbg_1)
      {
        std::cout << " red, blue, do not intersect: "
            << mb->id_from_handle(redT) << " "
            << mb->id_from_handle(currentBlue) << "\n";
      }

    }

    EntityHandle blueNeighbors[4];
    rval = GetOrderedNeighbors(mbs1, currentBlue, blueNeighbors);
    ERRORR(rval, "can't get neighbors");
    if (dbg_1)
    {
      std::cout << "Next: neighbors for blue T ";
      for (int kk = 0; kk < nsides; kk++)
      {
        if (blueNeighbors[kk] > 0)
          std::cout << mb->id_from_handle(blueNeighbors[kk]) << " ";
        else
          std::cout << 0 << " ";
      }
      std::cout << std::endl;
    }
    for (int j = 0; j < nsides; j++)
    {
      EntityHandle blueNeigh = blueNeighbors[j];
      unsigned char status = 1;
      if (blueNeigh == 0)
        continue;
      mb->tag_get_data(BlueFlagTag, &blueNeigh, 1, &status); // status 0 is unused
      if (status == 0 && n[j] > 0) // not treated yet and marked as a neighbor
      {
        // we identified red quad n[j] as intersecting with neighbor j of the blue quad
        blueQueue.push(blueNeigh);
        redQueue.push(n[j]);
        if (dbg_1)
          std::cout << "new quads pushed: blue, red:"
              << mb->id_from_handle(blueNeigh) << " "
              << mb->id_from_handle(n[j]) << std::endl;
        mb->tag_set_data(BlueFlagTag, &blueNeigh, 1, &used);
      }
    }

  }

  if (dbg_1)
  {
    for (int k = 0; k < 6; k++)
      mout_1[k].close();
  }
  //
  clean();
  return MB_SUCCESS;
}

// clean some memory allocated
void Intx2Mesh::clean()
{
  //
  int indx = 0;
  for (Range::iterator eit = RedEdges.begin(); eit != RedEdges.end();
      eit++, indx++)
  {
    //EntityHandle edge = *eit;
    //delete extraNodesMap[edge];
    delete extraNodesVec[indx];
  }
  //extraNodesMap.clear();
  extraNodesVec.clear();
}
ErrorCode Intx2Mesh::initialize_local_kdtree(EntityHandle euler_set)
{
  if (myTree)
    return MB_SUCCESS;
  // get the DP tag, and values for each node in the euler set
    // get forst all entities in the set
  parcomm = ParallelComm::get_pcomm(mb, 0);
  if (NULL==parcomm)
    return MB_FAILURE;

  Range local_ents;
  ErrorCode rval = mb->get_entities_by_dimension(euler_set, 2, local_ents);
  ERRORR(rval, "can't get ents by dimension");

  AdaptiveKDTree::Settings settings;
  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;

    //get entities on the local part
  ErrorCode result = MB_SUCCESS;


      // build the tree for local processor
  int numIts = 3;
  for (int i = 0; i < numIts; i++) {
    myTree = new AdaptiveKDTree(mb);
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
    else
      break; // get out of tree building
  }

    // get the bounding box for local tree
  if (parcomm)
    allBoxes.resize(6*parcomm->proc_config().proc_size());
  else allBoxes.resize(6);
  unsigned int my_rank = (parcomm ? parcomm->proc_config().proc_rank() : 0);
  result = myTree->get_tree_box(localRoot, &allBoxes[6*my_rank], &allBoxes[6*my_rank+3]);
  if (MB_SUCCESS != result) return result;

    // now communicate to get all boxes
    // use "in place" option
  if (parcomm) {
    int mpi_err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                                &allBoxes[0], 6, MPI_DOUBLE,
                                parcomm->proc_config().proc_comm());
    if (MPI_SUCCESS != mpi_err) return MB_FAILURE;
  }


  #ifndef NDEBUG
    double min[3] = {0,0,0}, max[3] = {0,0,0};
    unsigned int dep;
    myTree->get_info(localRoot, min, max, dep);
    std::cout << "Proc " << my_rank << ": box min/max, tree depth = ("
              << min[0] << "," << min[1] << "," << min[2] << "), ("
              << max[0] << "," << max[1] << "," << max[2] << "), "
              << dep << std::endl;
  #endif

    return result;
  return MB_SUCCESS;
}
// this will work in parallel only
ErrorCode Intx2Mesh::locate_departure_points(EntityHandle euler_set)
{
  // get the par comm ...
  // create local search tree if not created yet
  ErrorCode  rval = initialize_local_kdtree(euler_set);
  if (MB_SUCCESS!=rval)
    return rval;
  // get the DP tag, and values for each node in the euler set
  // get first all entities in the set
  Range local_ents;
  rval = mb->get_entities_by_dimension(euler_set, 2, local_ents);
  ERRORR(rval, "can't get ents by dimension");

  // get all the owned nodes on this euler set
  Range local_verts;
  rval = mb->get_connectivity(local_ents, local_verts);
  if (parcomm)
  {
    // get only those that are owned by this proc
    rval = parcomm->filter_pstatus(local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    ERRORR(rval, "can't filter verts");
    rval = parcomm->filter_pstatus(local_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    ERRORR(rval, "can't filter elements");
  }

  rval = locate_departure_points(local_ents, local_verts);
  ERRORR(rval, "can't locate departure points");

  return MB_SUCCESS;
}
// store
ErrorCode Intx2Mesh::locate_departure_points(Range & local_elements, Range & local_verts)
{
  // get the DP tag , and get the departure points
  Tag tagh = 0;
  std::string tag_name("DP");
  ErrorCode rval = mb->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE);
  ERRORR(rval, "can't get DP tag");
  int num_owned_verts = (int)local_verts.size();

  std::vector<double> dep_points(3*num_owned_verts);
  rval = mb->tag_get_data(tagh, local_verts, (void*)&dep_points[0]);
  ERRORR(rval, "can't get DP tag values");

  // create LOC tag
  std::string loc_tag("LOC");
  Tag tag_loc = 0;
  rval = mb->tag_get_handle(loc_tag.c_str(), 2, MB_TYPE_INTEGER, tag_loc, MB_TAG_DENSE|MB_TAG_CREAT);
  ERRORR(rval, "can't create LOC tag");

  // replicate some of coupler algorithm locate_points
  // differences: work in 2d, so we will be interested in interior or boundary of a quad/triangle
  // we will not store yet the parametric position; sphere geometry involves gnomonic projection
    // target_pts: TL(to_proc, tgt_index, x, y, z): tuples sent to source mesh procs representing pts to be located
    // source_pts: TL(from_proc, tgt_index, src_index): results of source mesh proc point location, ready to send
    //             back to tgt procs; src_index of -1 indicates point not located (arguably not useful...)
  TupleList target_pts;
  target_pts.initialize(2, 0, 0, 3, num_owned_verts);
  target_pts.enableWriteAccess();

  TupleList source_pts;

  // store the remote handle (element) where the point is located in
  // we do not need another array like mappedPts or locs
  source_pts.initialize(2, 0, 1, 0, target_pts.get_max());
  source_pts.enableWriteAccess();

  source_pts.set_n(0);
  ErrorCode result;

    // for each point, find box(es) containing the point,
    // appending results to tuple_list;
    // keep local points separately, in local_pts, which has pairs
    // of <local_index, mapped_index>, where mapped_index is the index
    // of <local_index, mapped_index>, where mapped_index is the index
    // into the locs list of elements where the points were located

  unsigned int my_rank = (parcomm ? parcomm->proc_config().proc_rank() : 0);


  for (int i = 0; i < 3*num_owned_verts; i+=3)
  {
    for (unsigned int j = 0; j < (parcomm ? parcomm->proc_config().proc_size() : 0); j++)
    {
        // test if point is in proc's box
      if ( (allBoxes[6*j] <= dep_points[i] + box_error) && ( dep_points[i] <= allBoxes[6*j+3]+ box_error) &&
           (allBoxes[6*j+1] <= dep_points[i+1]+ box_error) && (dep_points[i+1] <= allBoxes[6*j+4]+ box_error) &&
           (allBoxes[6*j+2] <= dep_points[i+2]+ box_error) && (dep_points[i+2] <= allBoxes[6*j+5]+ box_error))
      {
          // if in this proc's box, will send to proc to test further
          // check size, grow if we're at max
        if (target_pts.get_n() == target_pts.get_max())
          target_pts.resize(std::max(10.0, 1.5*target_pts.get_max()));

        target_pts.vi_wr[2*target_pts.get_n()] = j; // send to processor j, because its box has the point
        target_pts.vi_wr[2*target_pts.get_n()+1] = i/3; // index in the local_verts range
        target_pts.vr_wr[3*target_pts.get_n()] = dep_points[i];  // departure position, of the node local_verts[i]
        target_pts.vr_wr[3*target_pts.get_n()+1] = dep_points[i+1];
        target_pts.vr_wr[3*target_pts.get_n()+2] = dep_points[i+2];
        target_pts.inc_n();
      }
      else
      {
        continue;
      }
    }
  }

  int num_to_me = 0;
  for (unsigned int i = 0; i < target_pts.get_n(); i++)
    if (target_pts.vi_rd[2*i] == (int)my_rank) num_to_me++;

  printf("rank: %d local points: %d, nb sent target pts: %d  num to me: %d \n",
         my_rank, num_owned_verts, target_pts.get_n(), num_to_me);
    // perform scatter/gather, to gather points to source mesh procs
  if (parcomm) {
    (parcomm->proc_config().crystal_router())->gs_transfer(1, target_pts, 0);

    int num_from_me = 0;
    for (unsigned int i = 0; i < target_pts.get_n(); i++)
      if (target_pts.vi_rd[2*i] == (int)my_rank) num_from_me++;

    printf("rank: %d after first gs nb received_pts: %d; num_from_me = %d\n",
           my_rank, target_pts.get_n(), num_from_me);
      // after scatter/gather:
      // target_pts.set_n( # points local proc has to map );
      // target_pts.vi_wr[2*i] = proc sending point i
      // target_pts.vi_wr[2*i+1] = index of point i on sending proc
      // target_pts.vr_wr[3*i..3*i+2] = xyz of point i
      //
      // Mapping builds the tuple list:
      // source_pts.set_n (target_pts.get_n() )
      // source_pts.vi_wr[3*i] = target_pts.vi_wr[2*i] = sending proc
      // source_pts.vi_wr[3*i+1] = index of point i on sending proc
      // source_pts.vi_wr[3*i+2] = index of mapped point (-1 if not mapped)
      //

      // test target points against my elements
    for (unsigned i = 0; i < target_pts.get_n(); i++)
    {
      // now, for each point, see if we have a local element that has it
      result = test_local_box(target_pts.vr_wr+3*i, // position of point to be located on local proc
                              target_pts.vi_rd[2*i], // proc sending point i; we will send back the info about the entity it was located in
                              target_pts.vi_rd[2*i+1], // index of point i on sending proc
                              &source_pts);
      if (MB_SUCCESS != result) return result;
    }

      // no longer need target_pts
    target_pts.reset();

    printf("rank: %d nb sent source pts: %d\n",
           my_rank, source_pts.get_n() );
      // send target points back to target procs
    (parcomm->proc_config().crystal_router())->gs_transfer(1, source_pts, 0);

    printf("rank: %d nb received source pts: %d\n",
           my_rank, source_pts.get_n());
  }
  // so now, after second router call, source_pts contain:
  // proc were it was located, index of the point in the local_list (remote index is now again local)
  // and the index of the element in the locs array it was located in
  // if this last index is -1, it means the point was not located on the proc it was sent to
  // we should order by the local index first, and see if we have one point that was not located
  int max_size = 3*source_pts.get_n();
  TupleList::buffer sort_buffer;
  sort_buffer.buffer_init(max_size);
  //source_pts.print("print1");
  source_pts.sort(1, &sort_buffer);// this is the local index; for each point we need at least one element and one proc
  //source_pts.print("print2");

  // now, on each of the local points, find the first location, and set it
  // now, we have one each target point, the proc and corresponding element local handle
  // source_pts.get_n( #  );/ we will also get 0 in el handle, which is not good
  // after sorting by index in the local_verts, we loop the source_pts
  int N=(int)source_pts.get_n();
  Range found_points;
  for (int i=0; i<N; i++)
  {
    // we will look at the first point
    EntityHandle local_element=source_pts.vul_rd[i];
    if (0==local_element)
      continue;
    Range::iterator it = local_elements.find(local_element);
    if (it==local_elements.end())
    {
      //big problem
      ERRORR(MB_FAILURE, "element is not locally owned, not found in the local range");
    }
    int index_in_el_range = (int)(it-local_elements.begin());


    int location[2]={source_pts.vi_rd[2*i], index_in_el_range};// processor, index in local range
    EntityHandle local_vert=local_verts[source_pts.vi_rd[2*i+1]];
    rval = mb->tag_set_data(tag_loc, &local_vert, 1, (void*)location);
    ERRORR(rval, "can't set the location tag");
    found_points.insert(local_vert);
  }
  printf("found points on this proc: %d\n", (int)found_points.size());
  return MB_SUCCESS;
      // tl.vi_wr[2*i] = remote proc mapping point
      // tl.vi_wr[2*i+1] = local index of mapped point
      // tl.vul_wr[i] = element handle in local ents; find the index
#if 0
    // source_pts.set_n( # mapped points );
    // tl.vi_wr[2*i] = remote proc mapping point
    // tl.vi_wr[2*i+1] = local index of mapped point
    // tl.vul_wr[i] = element handle in local ents; find the index
    //
    // Local index is mapped into either myRange, holding the handles of
    // local mapped entities, or myXyz, holding locations of mapped pts

    // store information about located points
  TupleList *tl_tmp;
  if (!store_local)
    tl_tmp = tl;
  else {
    targetPts = new TupleList();
    tl_tmp = targetPts;
  }

  tl_tmp->initialize(3, 0, 0, 0, num_points);
  tl_tmp->set_n(num_points); // automatically sets tl to write_enabled
    // initialize so we know afterwards how many pts weren't located
  std::fill(tl_tmp->vi_wr, tl_tmp->vi_wr+3*num_points, -1);

  unsigned int local_pts = 0;
  for (unsigned int i = 0; i < source_pts.get_n(); i++) {
    if (-1 != source_pts.vi_rd[3*i+2]) { //why bother sending message saying "i don't have the point" if it gets discarded?
      if (source_pts.vi_rd[3*i] == (int)my_rank) local_pts++;
      int tgt_index = 3*source_pts.vi_rd[3*i+1];
      tl_tmp->vi_wr[tgt_index]   = source_pts.vi_rd[3*i];
      tl_tmp->vi_wr[tgt_index+1] = source_pts.vi_rd[3*i+1];
      tl_tmp->vi_wr[tgt_index+2] = source_pts.vi_rd[3*i+2];
    }
  }

    // count missing points
  unsigned int missing_pts = 0;
  for (int i = 0; i < num_points; i++)
    if (tl_tmp->vi_rd[3*i+1] == -1) missing_pts++;

  printf("rank: %d point location: wanted %d got %u locally, %d remote, missing %d\n",
         my_rank, num_points, local_pts,  num_points-missing_pts-local_pts, missing_pts);
  assert(0==missing_pts); //will litely break on curved geometries

    // no longer need source_pts
  source_pts.reset();

    // copy into tl if passed in and storing locally
  if (tl && store_local) {
    tl = new TupleList(3, 0, 0, 0, num_points);
    tl->enableWriteAccess();
    memcpy(tl->vi_wr, tl_tmp->vi_rd, 3*tl_tmp->get_n()*sizeof(int));
    tl->set_n( tl_tmp->get_n() );
    tl->disableWriteAccess();
  }

  tl_tmp->disableWriteAccess();
#endif
    // done
  return MB_SUCCESS;
}

ErrorCode Intx2Mesh::test_local_box(double *xyz,
                                      int from_proc, int remote_index,
                                  TupleList *tl)
{
  std::vector<EntityHandle> entities;

  ErrorCode result = inside_entities(xyz, entities);
  if (MB_SUCCESS != result) return result;

    // if we didn't find any ents and we're looking locally, nothing more to do
  if (entities.empty())
  {
    if (tl->get_n() == tl->get_max())
      tl->resize(std::max(10.0, 1.5*tl->get_max()));

    tl->vi_wr[2 * tl->get_n()] = from_proc;
    tl->vi_wr[2 * tl->get_n() + 1] = remote_index;
    tl->vul_wr[ tl->get_n() ] = 0 ;
    tl->inc_n();

    return MB_SUCCESS;
  }

  std::vector<EntityHandle>::iterator eit = entities.begin();

  for (; eit != entities.end(); eit++) {

    if (tl->get_n() == tl->get_max())
      tl->resize(std::max(10.0, 1.5*tl->get_max()));

      // store in tuple source_pts ; send back to source processor:
    //  remote_index, which is the index of the original point in the local list
    // and the index in the locs vector, which is actually locations (elements)
    tl->vi_wr[2*tl->get_n()] = from_proc;// send it back to the processor that send this point
    tl->vi_wr[2*tl->get_n()+1] = remote_index;
    tl->vul_wr[tl->get_n()] = *eit;// so this will be the element in that has the point
    tl->inc_n();
  }

  //point_located = true;

  return MB_SUCCESS;
}

ErrorCode Intx2Mesh::inside_entities(double xyz[3],
                                 std::vector<EntityHandle> &entities)
{
  AdaptiveKDTreeIter treeiter;
  ErrorCode result = myTree->get_tree_iterator(localRoot, treeiter);
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting iterator" << std::endl;
    return result;
  }

  double epsilon = this->box_error;//1.e-10; try to get some close boxes, because of the
  // the curvature that can appear on a sphere
  std::vector<EntityHandle> leaves;
  if (epsilon) {
    std::vector<double> dists;

    result = myTree->leaves_within_distance(localRoot, xyz, epsilon, leaves, &dists);
    if (leaves.empty())
      // not found returns success here, with empty list, just like case with no epsilon
      return MB_SUCCESS;

  }
  else {
    result = myTree->leaf_containing_point(localRoot, xyz, treeiter);
    if(MB_ENTITY_NOT_FOUND==result) //point is outside of myTree's bounding box
      return MB_SUCCESS;
    else if (MB_SUCCESS != result) {
      std::cout << "Problems getting leaf \n";
      return result;
    }
    leaves.push_back(treeiter.handle());
  }

  Range range_leaf;
  for (unsigned int i=0; i< leaves.size(); i++)
  {
    // add to the range_leaf all candidates
    result = mb->get_entities_by_dimension(leaves[i], 2, range_leaf, false);
    if(result != MB_SUCCESS) std::cout << "Problem getting leaf in a range" << std::endl;
  }

  // loop over the range_leaf
  for(Range::iterator iter = range_leaf.begin(); iter != range_leaf.end(); iter++)
  {
    //test to find out in which entity the point is
    //get the EntityType and create the appropriate Element::Map subtype
    // if spectral, do not need coordinates, just the GL points
    EntityHandle eh=*iter;
    if (is_inside_element(xyz, eh))
    {
      entities.push_back(eh);
    }

  }

  //didn't find any elements containing the point
  return MB_SUCCESS;
}

} /* namespace moab */
