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
#include "MBTagConventions.hpp"
// this is for DBL_MAX
#include <float.h>
#include <queue>
#include <sstream>
#include "moab/GeomUtil.hpp"

namespace moab {

Intx2Mesh::Intx2Mesh(Interface * mbimpl):mb(mbimpl), parcomm(NULL), myTree(NULL)
{
  dbg_1=0;
  box_error=0;
  my_rank=0;
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
  Range redElements;
  //Range redEdges;
  rval = mb->get_entities_by_dimension(mbs2, 2, redElements, false); // so all tri, quad, poly
  ERRORV(rval, "can't get ents by dimension");

  // create red edges if they do not exist yet; so when they are looked upon, they are found
  // this is the only call that is potentially NlogN, in the whole method
  rval = mb->get_adjacencies(redElements, 1, true, RedEdges, Interface::UNION);
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

  rval = mb->tag_get_handle("RedParent", 1, MB_TYPE_INTEGER, redParentTag,
      MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);
  ERRORV(rval, "can't create positive tag");

  rval = mb->tag_get_handle("BlueParent", 1, MB_TYPE_INTEGER, blueParentTag,
      MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);
  ERRORV(rval, "can't create negative tag");

  rval = mb->tag_get_handle("Counting", 1, MB_TYPE_INTEGER, countTag,
        MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);
  ERRORV(rval, "can't create Counting tag");

  return;
}


// specify also desired set; we are interested only in neighbors in the set!
// we should always get manifold mesh, each edge is adjacent to 2 cell
// maybe we should check that first, just in case
ErrorCode Intx2Mesh::GetOrderedNeighbors(EntityHandle set, EntityHandle cell,
    EntityHandle neighbors[MAXEDGES])
{
  int nnodes = 3;
  // will get the nnodes ordered neighbors;
  // first cell is for nodes 0, 1, second to 1, 2, third to 2, 3, last to nnodes-1,
  const EntityHandle * conn4;
  ErrorCode rval = mb->get_connectivity(cell, conn4, nnodes);
  int nsides = nnodes; // just keep it for historical purposes; it is indeed nnodes
  ERRORR(rval, "can't get connectivity on an element");
  for (int i = 0; i < nsides; i++)
  {
    EntityHandle v[2];
    v[0] = conn4[i];
    v[1] = conn4[(i + 1) % nsides];
    // get quads adjacent to vertices
    std::vector<EntityHandle> cells;
    std::vector<EntityHandle> cellsInSet;
    rval = mb->get_adjacencies(v, 2, 2, false, cells, Interface::INTERSECT);
    ERRORR(rval, "can't get adjacencies on 2 nodes");
    size_t siz = cells.size();
    for (size_t j = 0; j < siz; j++)
      if (mb->contains_entities(set, &(cells[j]), 1))
        cellsInSet.push_back(cells[j]);
    siz = cellsInSet.size();

    if (siz > 2)
    {
      std::cout << "non manifold mesh, error"
          << mb->list_entities(&(cellsInSet[0]), cellsInSet.size()) << "\n";
      return MB_FAILURE; // non-manifold
    }
    if (siz == 1)
    {
      // it must be the border,
      neighbors[i] = 0; // we are guaranteed that ids are !=0; this is marking a border
      // borders do not appear for a sphere in serial, but they do appear for
      // parallel processing anyway
      continue;
    }
    // here siz ==2, it is either the first or second
    if (cell == cellsInSet[0])
      neighbors[i] = cellsInSet[1];
    else
      neighbors[i] = cellsInSet[0];
  }
  return MB_SUCCESS;
}
// main interface; this will do the advancing front trick
// some are triangles, some are quads, some are polygons ...
ErrorCode Intx2Mesh::intersect_meshes(EntityHandle mbset1, EntityHandle mbset2,
     EntityHandle & outputSet)
{

  ErrorCode rval;
  mbs1 = mbset1; // set 1 is departure, and it is completely covering the euler set on proc
  mbs2 = mbset2;
  outSet = outputSet;

  // really, should be something from t1 and t2; blue is 1 (lagrange), red is 2 (euler)
  createTags(); //
  EntityHandle startBlue, startRed;

  mb->get_entities_by_dimension(mbs1, 2, rs1);
  mb->get_entities_by_dimension(mbs2, 2, rs2);
  Range rs22=rs2; // a copy of the initial range; we will remove from it elements as we
                 // advance ; rs2 is needed for marking the polygon to the red parent
  while (!rs22.empty())
  {
    for (Range::iterator it = rs1.begin(); it != rs1.end(); it++)
    {
      startBlue = *it;
      int found = 0;
      for (Range::iterator it2 = rs22.begin(); it2 != rs22.end() && !found; it2++)
      {
        startRed = *it2;
        double area = 0;
        // if area is > 0 , we have intersections
        double P[10*MAXEDGES]; // max 8 intx points + 8 more in the polygon
        //
        int nP = 0;
        int nb[MAXEDGES], nr[MAXEDGES]; // sides 3 or 4? also, check boxes first
        int nsRed, nsBlue;
        computeIntersectionBetweenRedAndBlue(startRed, startBlue, P, nP, area, nb, nr,
            nsBlue, nsRed, true);
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

    Range toResetBlues; // will be used to reset blue flags for every red quad
    // processed

    /*if (my_rank==0)
      dbg_1 = 1;*/
    unsigned char used = 1;
    unsigned char unused = 0; // for red flags
    // mark the start blue quad as used, so it will not come back again
    mb->tag_set_data(RedFlagTag, &startRed, 1, &used);
    while (!redQueue.empty())
    {
      // flags for the side : 0 means a blue quad not found on side
      // a paired blue not found yet for the neighbors of red
      EntityHandle n[MAXEDGES] = { EntityHandle(0) };

      int nsidesRed; // will be initialized later, when we compute intersection
      EntityHandle currentRed = redQueue.front();

      redQueue.pop();
      //        for (k=0; k<m_numPos; k++)
      //          redFlag[k] = 0;
      //        redFlag[m_numPos] = 1; // to guard for the boundary
      // all reds that were tagged, are now cleared
      for (Range::iterator itr = toResetBlues.begin(); itr != toResetBlues.end();
          itr++)
      {
        EntityHandle ttt = *itr;
        rval = mb->tag_set_data(BlueFlagTag, &ttt, 1, &unused);
        ERRORR(rval, "can't set blue unused tag");
      }
      //rval = mb2->tag_set_data(RedFlagTag, toResetReds, &unused);
      /*if (dbg_1)
      {
        std::cout << "reset blues: ";
        mb->list_entities(toResetBlues);
      }*/
      //rval = mb2->tag_set_data(RedFlagTag, toResetReds, &unused);
      if (dbg_1)
      {
        std::cout << "reset blues: ";
        for (Range::iterator itr = toResetBlues.begin(); itr != toResetBlues.end();
            itr++)
          std::cout << mb->id_from_handle(*itr) << " ";
        std::cout << std::endl;
      }
      EntityHandle currentBlue = blueQueue.front(); // where do we check for redQueue????
      // red and blue queues are parallel
      blueQueue.pop(); // mark the current red
      //redFlag[currentRed] = 1; //
      toResetBlues.clear(); // empty the range of used blues, will have to be set unused again,
      // at the end of red element processing
      toResetBlues.insert(currentBlue);
      rval = mb->tag_set_data(BlueFlagTag, &currentBlue, 1, &used);
      ERRORR(rval, "can't set blue tag");
      //mb2->set_tag_data
      std::queue<EntityHandle> localBlue;
      localBlue.push(currentBlue);
      while (!localBlue.empty())
      {
        //
        EntityHandle blueT = localBlue.front();
        localBlue.pop();
        double P[10*MAXEDGES], area; //
        int nP = 0;
        int nb[MAXEDGES] = {0};
        int nr[MAXEDGES] = {0};

        int nsidesBlue; ///
        // area is in 2d, points are in 3d (on a sphere), back-projected, or in a plane
        // intersection points could include the vertices of initial elements
        // nb [j] = 0 means no intersection on the side k for element blue (markers)
        // nb [j] = 1 means that the side j (from j to j+1) of blue poly intersects the
        // red poly.  A potential next poly is the red poly that is adjacent to this side
        computeIntersectionBetweenRedAndBlue(/* red */currentRed, blueT, P, nP,
            area, nb, nr, nsidesBlue, nsidesRed);
        if (nP > 0)
        {
          if (dbg_1)
          {
            for (int k=0; k<3; k++)
            {
              std::cout << " nb, nr: " << k << " " << nb[k] << " " << nr[k] << "\n";
            }
          }
          // intersection found: output P and original triangles if nP > 2

          EntityHandle neighbors[MAXEDGES];
          rval = GetOrderedNeighbors(mbs1, blueT, neighbors);
          if (rval != MB_SUCCESS)
          {
            std::cout << " can't get the neighbors for blue element "
                << mb->id_from_handle(blueT);
            return MB_FAILURE;
          }

          // add neighbors to the localBlue queue, if they are not marked
          for (int nn = 0; nn < nsidesBlue; nn++)
          {
            EntityHandle neighbor = neighbors[nn];
            if (neighbor > 0 && nb[nn]>0) // advance across blue boundary n
            {
              //n[nn] = redT; // start from 0!!
              unsigned char status = 0;
              mb->tag_get_data(BlueFlagTag, &neighbor, 1, &status);
              if (status == 0)
              {
                localBlue.push(neighbor);
                /*if (dbg_1)
                {
                  std::cout << " local blue elem "
                      << mb->list_entities(&neighbor, 1) << "\n  for red:"
                      << mb->list_entities(&currentRed, 1) << "\n";
                }*/
                if (dbg_1)
                {
                  std::cout << " local blue elem " << mb->id_from_handle(neighbor)
                      << " for red:" << mb->id_from_handle(currentRed) << "\n";
                  mb->list_entities(&neighbor, 1);
                }
                rval = mb->tag_set_data(BlueFlagTag, &neighbor, 1, &used);
                //redFlag[neighbor] = 1; // flag it to not be added anymore
                toResetBlues.insert(neighbor); // this is used to reset the red flag
              }
            }
            // n(find(nc>0))=ac;        % ac is starting candidate for neighbor
            if (nr[nn] > 0)
              n[nn] = blueT;

          }
          if (nP > 1) // this will also construct triangles/polygons in the new mesh, if needed
            findNodes(currentRed, nsidesRed, blueT, nsidesBlue, P, nP);
        }
        else if (dbg_1)
        /*{
          std::cout << " red, blue, do not intersect: "
              << mb->list_entities(&currentRed, 1) << " "
              << mb->list_entities(&blueT, 1) << "\n";
        }*/
        {
          std::cout << " red, blue, do not intersect: "
              << mb->id_from_handle(currentRed) << " "
              << mb->id_from_handle(blueT) << "\n";
        }

      } // end while (!localBlue.empty())
      // here, we are finished with redCurrent, take it out of the rs22 range (red, arrival mesh)
      rs22.erase(currentRed);
      // also, look at its neighbors, and add to the seeds a next one

      EntityHandle redNeighbors[MAXEDGES];
      rval = GetOrderedNeighbors(mbs2, currentRed, redNeighbors);
      ERRORR(rval, "can't get neighbors");
      if (dbg_1)
      {
        std::cout << "Next: neighbors for current red ";
        for (int kk = 0; kk < nsidesRed; kk++)
        {
          if (redNeighbors[kk] > 0)
            std::cout << mb->id_from_handle(redNeighbors[kk]) << " ";
          else
            std::cout << 0 << " ";
        }
        std::cout << std::endl;
      }
      for (int j = 0; j < nsidesRed; j++)
      {
        EntityHandle redNeigh = redNeighbors[j];
        unsigned char status = 1;
        if (redNeigh == 0)
          continue;
        mb->tag_get_data(RedFlagTag, &redNeigh, 1, &status); // status 0 is unused
        if (status == 0 && n[j] > 0) // not treated yet and marked as a neighbor
        {
          // we identified red quad n[j] as intersecting with neighbor j of the blue quad
          redQueue.push(redNeigh);
          blueQueue.push(n[j]);
          if (dbg_1)
            std::cout << "new polys pushed: blue, red:"
                << mb->id_from_handle(redNeigh) << " "
                << mb->id_from_handle(n[j]) << std::endl;
          mb->tag_set_data(RedFlagTag, &redNeigh, 1, &used);
        }
      }

    } // end while (!redQueue.empty())
  }
  if (dbg_1)
  {
    for (int k = 0; k < 6; k++)
      mout_1[k].close();
  }
  // before cleaning up , we need to settle the position of the intersection points
  // on the boundary edges
  // this needs to be collective, so we should maybe wait something
  rval = correct_intersection_points_positions();
  if (rval!=MB_SUCCESS)
  {
    std::cout << "can't correct position, Intx2Mesh.cpp \n";
  }
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
  // also, delete some bit tags, used to mark processed reds and blues
  mb->tag_delete(RedFlagTag);// to mark blue quads already considered
  mb->tag_delete(BlueFlagTag);

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

    //get entities on the local part
  ErrorCode result = MB_SUCCESS;

      // build the tree for local processor
  int numIts = 3;
  int max_per_leaf = 6;
  for (int i = 0; i < numIts; i++) {
    std::ostringstream ostr;
    ostr << "CANDIDATE_PLANE_SET=SUBDIVISION;MAX_PER_LEAF=" << max_per_leaf;
    FileOptions opts(ostr.str().c_str());
    myTree = new AdaptiveKDTree(mb);
    result = myTree->build_tree(local_ents, &localRoot, &opts);
    if (MB_SUCCESS != result) {
      std::cout << "Problems building tree";
      if (numIts != i) {
        delete myTree;
        max_per_leaf *= 2;
        std::cout << "; increasing elements/leaf to "
                  << max_per_leaf << std::endl;;
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
  my_rank = (parcomm ? parcomm->proc_config().proc_rank() : 0);
  BoundBox box;
  result = myTree->get_bounding_box(box);
  if (MB_SUCCESS != result) return result;
  box.bMin.get(&allBoxes[6*my_rank]);
  box.bMax.get(&allBoxes[6*my_rank+3]);

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
  localEnts.clear(); // just to make sure :)
  rval = mb->get_entities_by_dimension(euler_set, 2, localEnts);
  ERRORR(rval, "can't get ents by dimension");

  // get all the owned nodes on this euler set
  Range local_verts;
  rval = mb->get_connectivity(localEnts, local_verts);
  Range all_local_verts=local_verts;
  if (parcomm)
  {
    // get only those that are owned by this proc
    rval = parcomm->filter_pstatus(local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    ERRORR(rval, "can't filter verts");
    rval = parcomm->filter_pstatus(localEnts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    ERRORR(rval, "can't filter elements");
    // just in case we do not have global ids for vertices and elements :)
    rval = parcomm->check_global_ids(euler_set, 2);
    ERRORR(rval, "can't get global ids for elements and vertices");
  }

  rval = locate_departure_points(local_verts);
  ERRORR(rval, "can't locate departure points");
  // experience the exchange tags;
  // the LOC tag should be set on the owned nodes only, but try to
  if(parcomm)
  {
    Tag loctag;
    rval = mb->tag_get_handle("LOC", 2, MB_TYPE_INTEGER, loctag, MB_TAG_DENSE);
    ERRORR(rval, "can't get LOC tag");
    std::vector<Tag> tags;
    tags.push_back(loctag);

    rval = parcomm->exchange_tags(tags, tags, all_local_verts);
    ERRORR(rval, "can't exchange tags");
    /*for (Range::iterator it=all_local_verts.begin(); it!=all_local_verts.end(); it++)
    {
      int loc[2];
      EntityHandle vert=*it;
      rval = mb->tag_get_data(loctag, &vert, 1, loc);
      if (rval==MB_SUCCESS)
      {
        std::cout<< " after: vertex: "<< vert << " loc: "<< loc[0] << " " << loc[1]<< " \n";
      }
    }*/
  }
  return MB_SUCCESS;
}
// store
ErrorCode Intx2Mesh::locate_departure_points(Range & local_verts)
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
  int def_val[2]={-1,-1};
  rval = mb->tag_get_handle(loc_tag.c_str(), 2, MB_TYPE_INTEGER, tag_loc, MB_TAG_DENSE|MB_TAG_CREAT, (void*)def_val);
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
  source_pts.initialize(3, 0, 0, 0, target_pts.get_max());
  source_pts.enableWriteAccess();

  source_pts.set_n(0);
  ErrorCode result;

    // for each point, find box(es) containing the point,
    // appending results to tuple_list;
    // keep local points separately, in local_pts, which has pairs
    // of <local_index, mapped_index>, where mapped_index is the index
    // of <local_index, mapped_index>, where mapped_index is the index
    // into the locs list of elements where the points were located

  //my_rank = (parcomm ? parcomm->proc_config().proc_rank() : 0);


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
      // target_pts.vi_wr[3*i] = proc sending point i
      // target_pts.vi_wr[3*i+1] = index of point i on sending proc
      // target_pts.vr_wr[3*i..3*i+2] = xyz of point i
      //
      // Mapping builds the tuple list:
      // source_pts.set_n (target_pts.get_n() )
      // source_pts.vi_wr[3*i] = target_pts.vi_wr[2*i] = sending proc
      // source_pts.vi_wr[3*i+1] = index of point i on sending proc
      // source_pts.vi_wr[3*i+2] = index of element in localEnts (-1 if not found)
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
    int index_in_el_range=source_pts.vi_rd[3*i+2];
    if (-1==index_in_el_range)
      continue;

    int location[2]={source_pts.vi_rd[3*i], index_in_el_range};// processor, index in local range
    EntityHandle local_vert=local_verts[source_pts.vi_rd[3*i+1]];
    rval = mb->tag_set_data(tag_loc, &local_vert, 1, (void*)location);
    ERRORR(rval, "can't set the location tag");
    found_points.insert(local_vert);
  }
  printf("found points on this proc: %d\n", (int)found_points.size());
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

    tl->vi_wr[3 * tl->get_n()] = from_proc;
    tl->vi_wr[3 * tl->get_n() + 1] = remote_index;
    tl->vi_wr[3 * tl->get_n() + 2] = -1 ;
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
    tl->vi_wr[3*tl->get_n()] = from_proc;// send it back to the processor that send this point
    tl->vi_wr[3*tl->get_n()+1] = remote_index;
    EntityHandle eh=*eit;
    Range::iterator it = localEnts.find(eh);
    if (it == localEnts.end())
    {
      ERRORR(MB_FAILURE, "can't find the element in local entities");
    }
    tl->vi_wr[3* tl->get_n() + 2] = (int) (it-localEnts.begin());// so this will be the element in that has the point
    tl->inc_n();
  }

  //point_located = true;

  return MB_SUCCESS;
}

ErrorCode Intx2Mesh::inside_entities(double xyz[3],
                                 std::vector<EntityHandle> &entities)
{
  double epsilon = this->box_error;//1.e-10; try to get some close boxes, because of the
  // the curvature that can appear on a sphere
  std::vector<EntityHandle> leaves;
  ErrorCode result;
  if (epsilon) {
    std::vector<double> dists;

    result = myTree->distance_search(xyz, epsilon, leaves, 0.0, &dists);
    if (leaves.empty())
      // not found returns success here, with empty list, just like case with no epsilon
      return MB_SUCCESS;

  }
  else {
    EntityHandle leaf;
    result = myTree->point_search(xyz, leaf);
    if(MB_ENTITY_NOT_FOUND==result) //point is outside of myTree's bounding box
      return MB_SUCCESS;
    else if (MB_SUCCESS != result) {
      std::cout << "Problems getting leaf \n";
      return result;
    }
    leaves.push_back(leaf);
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
ErrorCode Intx2Mesh::create_departure_mesh(EntityHandle & covering_lagr_set)
{
  // the elements are already in localEnts; get all vertices, and DP and LOC tag
  // get the DP tag , and get the departure points
  std::map<int, EntityHandle> globalID_to_handle;
  std::map<int, EntityHandle> globalID_to_eh;
  std::map<int, Range> rs;
  std::set<int> ps; // processors working with my_rank
  Tag dpTag = 0;
  std::string tag_name("DP");
  ErrorCode rval = mb->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, dpTag, MB_TAG_DENSE);
  ERRORR(rval, "can't get DP tag");
  // get all local verts
  Range local_verts;
  rval = mb->get_connectivity(localEnts, local_verts);
  int num_local_verts = (int) local_verts.size();
  ERRORR(rval, "can't get local vertices");
  int num_owned_verts = (int)local_verts.size();
  std::vector<double> dep_points(3*num_owned_verts);
  rval = mb->tag_get_data(dpTag, local_verts, (void*)&dep_points[0]);
  ERRORR(rval, "can't get DP tag values");

  // get LOC tag, that was created before
  std::string loc_tag("LOC");
  Tag locTag = 0;
  rval = mb->tag_get_handle(loc_tag.c_str(), 2, MB_TYPE_INTEGER, locTag, MB_TAG_DENSE);
  ERRORR(rval, "can't get LOC tag");

  Tag gid;
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE);
  ERRORR(rval,"can't get global ID tag" );
  // create vertices for DP points in this rank
  // count how many to send; if LOC(0) != rank, we have to send to other LOC(0)
  unsigned int v_to_send =0, q_to_send =0;
  my_rank = 0;
  std::vector<int> locs(num_local_verts*2);// fill it with loc info for vertex
  if (this->parcomm)
  {
    //my_rank = parcomm->proc_config().proc_rank();
    rval = mb->tag_get_data(locTag, local_verts, &locs[0]);
    ERRORR(rval, "can't get value for LOC tag");

    // look at every owned quad
    for (Range::iterator itq=localEnts.begin(); itq!=localEnts.end(); itq++ )
    {
      const EntityHandle * conn4;
      EntityHandle q=*itq;
      int num_nodes;
      std::set<int> ptmp;
      rval = mb->get_connectivity(q, conn4, num_nodes);
      ERRORR(rval, "can't get connectivity for quad");
      for (int i=0; i<num_nodes; i++)
      {
        int v=conn4[i];
        int index = (int)(local_verts.find(v)-local_verts.begin());
        int proc = locs[2*index]; // where it was located
        ptmp.insert(proc);
      }
      for (std::set<int>::iterator sit=ptmp.begin(); sit!=ptmp.end(); sit++)
      {
        int proc = *sit;
        rs[proc].insert(q);
        std::copy(conn4, conn4+num_nodes, range_inserter(rs[proc]) );

      }
      ps.insert(ptmp.begin(), ptmp.end());
    }
  }
  std::vector<int> gids(num_local_verts);
  rval = mb->tag_get_data(gid, local_verts, &gids[0]);
  ERRORR(rval, "can't get gid tag");
  // count now how many vertices and quads are to be sent to other processors
  for (std::set<int>::iterator st=ps.begin(); st!=ps.end(); st++)
  {
    int to_proc=*st;
    if ((int)my_rank==to_proc)
      continue;
    Range & procRange = rs[to_proc];
    v_to_send += procRange.num_of_type( MBVERTEX);
    q_to_send += procRange.num_of_type(MBQUAD);
  }
  TupleList TLv;
  TupleList TLq;
  TLv.initialize(2, 0, 0, 3, v_to_send); // to proc, GLOBAL ID, DP points
  TLv.enableWriteAccess();

  TLq.initialize(6, 0, 0, 0, q_to_send); // to proc, elem GLOBAL ID, connectivity[4] (global ID v)
  TLq.enableWriteAccess();
  // is the global ID of the element needed? maybe not...
  // at least one of the nodes is in the interior of the proc area, so the whole element needs
  // to be copied on that processor ...
  if (parcomm)
  {
    Range local=rs[my_rank];
    Range local_dp_verts = local.subset_by_type(MBVERTEX);
    for (Range::iterator it=local_dp_verts.begin(); it!=local_dp_verts.end(); it++)
    {
      int lv = *it;
      unsigned int i = local_verts.find(lv)-local_verts.begin();
      double posDP[3]={dep_points[3*i], dep_points[3*i+1], dep_points[3*i+2] };

      // create a vertex at DP point
      EntityHandle new_vert;
      rval = mb->create_vertex(posDP, new_vert);
      ERRORR(rval, "can't create new vertex");
      globalID_to_handle[gids[i]]=new_vert;

    }
  }
  else
  {
    // need to send everything to local proc; why worry about nothing?
  }
  if (parcomm)
  {
    for (std::set<int>::iterator st=ps.begin(); st!=ps.end(); st++)
    {
      int to_proc=*st;
      if (to_proc==(int)my_rank)
        continue;
      Range & procRange = rs[to_proc];
      Range V = procRange.subset_by_type(MBVERTEX);

      for (Range::iterator it=V.begin(); it!=V.end(); it++)
      {
        EntityHandle v=*it;
        unsigned int index = local_verts.find(v)-local_verts.begin();
        int n=TLv.get_n();
        TLv.vi_wr[2*n] = to_proc; // send to processor
        TLv.vi_wr[2*n+1] = gids[index]; // global id needs index in the local_verts range
        TLv.vr_wr[3*n] = dep_points[3*index];  // departure position, of the node local_verts[i]
        TLv.vr_wr[3*n+1] = dep_points[3*index+1];
        TLv.vr_wr[3*n+2] = dep_points[3*index+2];
        TLv.inc_n();
      }
      // also, prep the quad for sending ...
      Range Q = procRange.subset_by_type(MBQUAD);
      for (Range::iterator it=Q.begin(); it!=Q.end(); it++)
      {
        EntityHandle q=*it;
        int global_id;
        rval = mb->tag_get_data(gid, &q, 1, &global_id);
        ERRORR(rval, "can't get gid for quad");
        int n=TLq.get_n();
        TLq.vi_wr[6*n] = to_proc; //
        TLq.vi_wr[6*n+1] = global_id; // global id of element, used to identify it ...
        const EntityHandle * conn4;
        int num_nodes;
        rval = mb->get_connectivity(q, conn4, num_nodes);
        ERRORR(rval, "can't get connectivity for quad");
        for (int i=0; i<num_nodes; i++)
        {
          EntityHandle v = conn4[i];
          unsigned int index = local_verts.find(v)-local_verts.begin();
          TLq.vi_wr[6*n+2+i] = gids[index];
        }
        TLq.inc_n();

      }

    }
    // now we are done populating the tuples; route them to the appropriate processors
    (parcomm->proc_config().crystal_router())->gs_transfer(1, TLv, 0);
    (parcomm->proc_config().crystal_router())->gs_transfer(1, TLq, 0);
    // now, look at every TLv, and see if we have to create a vertex there or not
    int n=TLv.get_n();// the size of the points received
    for (int i=0; i<n; i++)
    {
      int globalId = TLv.vi_rd[2*i+1];
      if (globalID_to_handle.find(globalId)==globalID_to_handle.end())
      {
        EntityHandle new_vert;
        double coords[3]= {TLv.vr_wr[3*i], TLv.vr_wr[3*i+1],  TLv.vr_wr[3*i+2]};
        rval = mb->create_vertex(coords, new_vert);
        ERRORR(rval, "can't create new vertex ");
        globalID_to_handle[globalId]= new_vert;
      }
    }
    // now, all dep points should be at their place (procs)
    // look in the local list of q for this proc, and create all those quads
    Range local=rs[my_rank];
    Range local_q = local.subset_by_type(MBQUAD);
    // the local should have all the vertices by now
    for (Range::iterator it=local_q.begin(); it!=local_q.end(); it++)
    {
      EntityHandle q=*it;
      int nnodes;
      const EntityHandle * conn4;
      rval = mb->get_connectivity(q, conn4, nnodes);
      ERRORR(rval, "can't get connectivity of local q ");
      EntityHandle new_conn[4];
      for (int i=0; i<nnodes; i++)
      {
        EntityHandle v1=conn4[i];
        unsigned int index = local_verts.find(v1)-local_verts.begin();
        assert(globalID_to_handle.find(gids[index])!=globalID_to_handle.end());
        new_conn[i] = globalID_to_handle[gids[index]];
      }
      EntityHandle new_quad;
      rval = mb->create_element(MBQUAD, new_conn, 4, new_quad);
      ERRORR(rval, "can't create new quad ");
      rval = mb->add_entities(covering_lagr_set, &new_quad, 1);
      ERRORR(rval, "can't add new quad to dep set");
      int gid_el;
      // get the global ID of the initial quad
      rval=mb->tag_get_data(gid, &q, 1, &gid_el);
      ERRORR(rval, "can't get element global ID ");
      globalID_to_eh[gid_el]=new_quad;
    }
    // now look at all quads received through; we do not want to duplicate them
    n=TLq.get_n();// number of quads received by this processor
    for (int i=0; i<n; i++)
    {
      int globalIdEl = TLq.vi_rd[6*i+1];
      // do we already have a quad with this global ID, represented?
      if (globalID_to_eh.find(globalIdEl)==globalID_to_eh.end())
      {
        // construct the conn quad
        EntityHandle new_conn[4];
        for (int j=0; j<4; j++)
        {
          int vgid = TLq.vi_rd[6*i+2+j];// vertex global ID
          assert(globalID_to_handle.find(vgid)!=globalID_to_handle.end());
          new_conn[j]=globalID_to_handle[vgid];
        }
        EntityHandle new_quad;
        rval = mb->create_element(MBQUAD, new_conn, 4, new_quad);
        ERRORR(rval, "can't create new quad ");
        globalID_to_eh[globalIdEl]=new_quad;
        rval = mb->add_entities(covering_lagr_set, &new_quad, 1);
        ERRORR(rval, "can't add new quad to dep set");
      }
    }
  }
  // send global ID of the euler vertex and DP position


  return MB_SUCCESS;
}

ErrorCode Intx2Mesh::create_departure_mesh_2nd_alg(EntityHandle & euler_set, EntityHandle & covering_lagr_set)
{
  // compute the bounding box on each proc
  parcomm = ParallelComm::get_pcomm(mb, 0);
  if (NULL==parcomm)
    return MB_FAILURE;

  localEnts.clear();
  ErrorCode rval = mb->get_entities_by_dimension(euler_set, 2, localEnts);
  ERRORR(rval, "can't get ents by dimension");

  Tag dpTag = 0;
  std::string tag_name("DP");
  rval = mb->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, dpTag, MB_TAG_DENSE);
  ERRORR(rval, "can't get DP tag");
  // get all local verts
  Range local_verts;
  rval = mb->get_connectivity(localEnts, local_verts);
  int num_local_verts = (int) local_verts.size();
  ERRORR(rval, "can't get local vertices");

  Tag gid;
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE);
  ERRORR(rval,"can't get global ID tag" );
  std::vector<int> gids(num_local_verts);
  rval = mb->tag_get_data(gid, local_verts, &gids[0]);
  ERRORR(rval, "can't get local vertices gids");

  // get the position of local vertices, and decide local boxes (allBoxes...)
  double bmin[3]={DBL_MAX, DBL_MAX, DBL_MAX};
  double bmax[3] ={-DBL_MAX, -DBL_MAX, -DBL_MAX};

  std::vector<double> coords(3*num_local_verts);
  rval = mb->get_coords(local_verts, &coords[0]);

  for (int i=0; i< num_local_verts; i++)
  {
    for (int k=0; k<3; k++)
    {
      double val=coords[3*i+k];
      if (val < bmin[k])
        bmin[k]=val;
      if (val > bmax[k])
        bmax[k] = val;
    }
  }
  int numprocs=parcomm->proc_config().proc_size();
  allBoxes.resize(6*numprocs);

  my_rank = parcomm->proc_config().proc_rank() ;
  for (int k=0; k<3; k++)
  {
    allBoxes[6*my_rank+k]=bmin[k];
    allBoxes[6*my_rank+3+k] = bmax[k];
  }

   // now communicate to get all boxes
   // use "in place" option
  int mpi_err = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                               &allBoxes[0], 6, MPI_DOUBLE,
                               parcomm->proc_config().proc_comm());
  if (MPI_SUCCESS != mpi_err) return MB_FAILURE;

  if (my_rank==0)
  {
    for (int i=0; i<numprocs; i++)
    {
      std::cout<<"proc: " << i << " box min: " << allBoxes[6*i  ] << " " <<allBoxes[6*i+1] << " " << allBoxes[6*i+2]  << " \n";
      std::cout<<          "        box max: " << allBoxes[6*i+3] << " " <<allBoxes[6*i+4] << " " << allBoxes[6*i+5]  << " \n";
    }
  }


  // now see the departure points; to what boxes should we send them?
  std::vector<double> dep_points(3*num_local_verts);
  rval = mb->tag_get_data(dpTag, local_verts, (void*)&dep_points[0]);
  ERRORR(rval, "can't get DP tag values");
  // ranges to send to each processor; will hold vertices and elements (quads?)
  // will look if the box of the dep quad covers box of of euler mesh on proc (with tolerances)
  std::map<int, Range> Rto;

  for (Range::iterator eit = localEnts.begin(); eit!=localEnts.end(); eit++)
  {
    EntityHandle q=*eit;
    const EntityHandle * conn4;
    int num_nodes;
    rval= mb->get_connectivity(q, conn4, num_nodes);
    ERRORR(rval, "can't get DP tag values");
    CartVect qbmin(DBL_MAX);
    CartVect qbmax(-DBL_MAX);
    for (int i=0; i<num_nodes; i++)
    {
      EntityHandle v=conn4[i];
      size_t index=local_verts.find(v)-local_verts.begin();
      CartVect dp( &dep_points[3*index] ); // will use constructor
      for (int j=0; j<3; j++)
      {
        if (qbmin[j]>dp[j])
          qbmin[j]=dp[j];
        if (qbmax[j]<dp[j])
          qbmax[j]=dp[j];
      }
    }
    for (int p=0; p<numprocs; p++)
    {
      CartVect bbmin(&allBoxes[6*p]);
      CartVect bbmax(&allBoxes[6*p+3]);
      if ( GeomUtil::boxes_overlap( bbmin, bbmax, qbmin, qbmax, box_error) )
      {
        Rto[p].insert(q);
      }
    }
  }

  // now, build TLv and TLq, for each p
  size_t numq=0;
  size_t numv=0;
  for (int p=0; p<numprocs; p++)
  {
    if (p==(int)my_rank)
      continue; // do not "send" it, because it is already here
    Range & range_to_P = Rto[p];
    // add the vertices to it
    if (range_to_P.empty())
      continue;// nothing to send to proc p
    Range vertsToP;
    rval = mb->get_connectivity(range_to_P, vertsToP);
    ERRORR(rval, "can't get connectivity");
    numq=numq+range_to_P.size();
    numv=numv+vertsToP.size();
    range_to_P.merge(vertsToP);
  }
  TupleList TLv;
  TupleList TLq;
  TLv.initialize(2, 0, 0, 3, numv); // to proc, GLOBAL ID, DP points
  TLv.enableWriteAccess();

  int sizeTuple = 2+MAXEDGES; // max edges is now 10 :) for polygons
  TLq.initialize(12, 0, 0, 0, numq); // to proc, elem GLOBAL ID, connectivity[10] (global ID v)
  TLq.enableWriteAccess();
  std::cout << "from proc " << my_rank << " send " << numv << " vertices and " << numq << " elements\n";

  for (int to_proc=0; to_proc<numprocs; to_proc++)
  {
    if (to_proc==(int)my_rank)
      continue;
    Range & range_to_P = Rto[to_proc];
    Range V = range_to_P.subset_by_type(MBVERTEX);

    for (Range::iterator it=V.begin(); it!=V.end(); it++)
    {
      EntityHandle v=*it;
      unsigned int index = local_verts.find(v)-local_verts.begin();
      int n=TLv.get_n();
      TLv.vi_wr[2*n] = to_proc; // send to processor
      TLv.vi_wr[2*n+1] = gids[index]; // global id needs index in the local_verts range
      TLv.vr_wr[3*n] = dep_points[3*index];  // departure position, of the node local_verts[i]
      TLv.vr_wr[3*n+1] = dep_points[3*index+1];
      TLv.vr_wr[3*n+2] = dep_points[3*index+2];
      TLv.inc_n();
    }
    // also, prep the quad for sending ...
    Range Q = range_to_P.subset_by_dimension(2);
    for (Range::iterator it=Q.begin(); it!=Q.end(); it++)
    {
      EntityHandle q=*it;
      int global_id;
      rval = mb->tag_get_data(gid, &q, 1, &global_id);
      ERRORR(rval, "can't get gid for polygon");
      int n=TLq.get_n();
      TLq.vi_wr[sizeTuple*n] = to_proc; //
      TLq.vi_wr[sizeTuple*n+1] = global_id; // global id of element, used to identify it ...
      const EntityHandle * conn4;
      int num_nodes;
      rval = mb->get_connectivity(q, conn4, num_nodes);// could be up to 10;
      ERRORR(rval, "can't get connectivity for quad");
      if (num_nodes > MAXEDGES)
        ERRORR(MB_FAILURE, "too many nodes in a polygon");
      for (int i=0; i<num_nodes; i++)
      {
        EntityHandle v = conn4[i];
        unsigned int index = local_verts.find(v)-local_verts.begin();
        TLq.vi_wr[sizeTuple*n+2+i] = gids[index];
      }
      for (int k=num_nodes; k<MAXEDGES; k++)
      {
        TLq.vi_wr[sizeTuple*n+2+k] = 0; // fill the rest of node ids with 0; we know that the node ids start from 1!
      }
      TLq.inc_n();

    }

  }
  // now we can route them to each processor
  // now we are done populating the tuples; route them to the appropriate processors
  (parcomm->proc_config().crystal_router())->gs_transfer(1, TLv, 0);
  (parcomm->proc_config().crystal_router())->gs_transfer(1, TLq, 0);
  // the elements are already in localEnts;

  // maps from global ids to new vertex and quad handles, that are added
  std::map<int, EntityHandle> globalID_to_handle;
  std::map<int, EntityHandle> globalID_to_eh;
  // now, look at every TLv, and see if we have to create a vertex there or not
  int n=TLv.get_n();// the size of the points received
  for (int i=0; i<n; i++)
  {
    int globalId = TLv.vi_rd[2*i+1];
    if (globalID_to_handle.find(globalId)==globalID_to_handle.end())
    {
      EntityHandle new_vert;
      double dp_pos[3]= {TLv.vr_wr[3*i], TLv.vr_wr[3*i+1],  TLv.vr_wr[3*i+2]};
      rval = mb->create_vertex(dp_pos, new_vert);
      ERRORR(rval, "can't create new vertex ");
      globalID_to_handle[globalId]= new_vert;
    }
  }

  // now, all dep points should be at their place
  // look in the local list of q for this proc, and create all those quads and vertices if needed
  // it may be an overkill, but because it does not involve communication, we do it anyway
  Range & local=Rto[my_rank];
  Range local_q = local.subset_by_dimension(2);
  // the local should have all the vertices in local_verts
  for (Range::iterator it=local_q.begin(); it!=local_q.end(); it++)
  {
    EntityHandle q=*it;
    int nnodes;
    const EntityHandle * conn4;
    rval = mb->get_connectivity(q, conn4, nnodes);
    ERRORR(rval, "can't get connectivity of local q ");
    EntityHandle new_conn[MAXEDGES];
    for (int i=0; i<nnodes; i++)
    {
      EntityHandle v1=conn4[i];
      unsigned int index = local_verts.find(v1)-local_verts.begin();
      int globalId=gids[index];
      if(globalID_to_handle.find(globalId)==globalID_to_handle.end())
      {
        // we need to create that vertex, at this position dep_points
        double dp_pos[3]={dep_points[3*index], dep_points[3*index+1], dep_points[3*index+2]};
        EntityHandle new_vert;
        rval = mb->create_vertex(dp_pos, new_vert);
        ERRORR(rval, "can't create new vertex ");
        globalID_to_handle[globalId]= new_vert;
      }
      new_conn[i] = globalID_to_handle[gids[index]];
    }
    EntityHandle new_element;
    //
    EntityType entType = MBQUAD;
    if (nnodes >4)
      entType = MBPOLYGON;
    if (nnodes <4)
      entType = MBTRI;

    rval = mb->create_element(entType, new_conn, nnodes, new_element);
    ERRORR(rval, "can't create new quad ");
    rval = mb->add_entities(covering_lagr_set, &new_element, 1);
    ERRORR(rval, "can't add new element to dep set");
    int gid_el;
    // get the global ID of the initial quad
    rval=mb->tag_get_data(gid, &q, 1, &gid_el);
    ERRORR(rval, "can't get element global ID ");
    globalID_to_eh[gid_el]=new_element;
  }
  // now look at all elements received through; we do not want to duplicate them
  n=TLq.get_n();// number of elements received by this processor
  for (int i=0; i<n; i++)
  {
    int globalIdEl = TLq.vi_rd[sizeTuple*i+1];
    // do we already have a quad with this global ID, represented?
    if (globalID_to_eh.find(globalIdEl)==globalID_to_eh.end())
    {
      // construct the conn quad
      EntityHandle new_conn[MAXEDGES];
      int nnodes;
      for (int j=0; j<MAXEDGES; j++)
      {
        int vgid = TLq.vi_rd[sizeTuple*i+2+j];// vertex global ID
        if (vgid==0)
          new_conn[j] = 0;
        else
        {
          assert(globalID_to_handle.find(vgid)!=globalID_to_handle.end());
          new_conn[j]=globalID_to_handle[vgid];
          nnodes = j+1;// nodes are at the beginning, and are variable number
        }
      }
      EntityHandle new_element;
      //
      EntityType entType = MBQUAD;
      if (nnodes >4)
        entType = MBPOLYGON;
      if (nnodes <4)
        entType = MBTRI;
      rval = mb->create_element(entType, new_conn, nnodes, new_element);
      ERRORR(rval, "can't create new element ");
      globalID_to_eh[globalIdEl]=new_element;
      rval = mb->add_entities(covering_lagr_set, &new_element, 1);
      ERRORR(rval, "can't add new element to dep set");
    }
  }
  return MB_SUCCESS;
}
// this method will reduce number of nodes, collapse edges that are of length 0
  // so a polygon like 428 431 431 will become a line 428 431
  // or something like 428 431 431 531 -> 428 431 531
void Intx2Mesh::correct_polygon(EntityHandle * nodes, int & nP)
{
  int i = 0;
  while(i<nP)
  {
    int nextIndex = (i+1)%nP;
    if (nodes[i]==nodes[nextIndex])
    {
      // we need to reduce nP, and collapse nodes
      if (dbg_1)
      {
        std::cout<<" nodes duplicated in list: " ;
        for (int j=0; j<nP; j++)
          std::cout<<nodes[j] << " " ;
        std::cout<<"\n";
        std::cout<<" node " << nodes[i] << " at index " << i << " is duplicated" << "\n";
      }
      // this will work even if we start from 1 2 3 1; when i is 3, we find nextIndex is 0, then next thing does nothing
      //  (nP-1 is 3, so k is already >= nP-1); it will result in nodes -> 1, 2, 3
      for (int k=i; k<nP-1; k++)
        nodes[k] = nodes[k+1];
      nP--; // decrease the number of nodes; also, decrease i, just if we may need to check again
      i--;
    }
    i++;
  }
  return;
}

ErrorCode Intx2Mesh::correct_intersection_points_positions()
{
  if (parcomm)
  {
    // first, find out the edges that are shared between processors, and owned by the current processor
    Range shared_edges_owned;
    ErrorCode rval = parcomm->get_shared_entities(-1, // all other proc
        shared_edges_owned,
        1,
        true, // only on the interface
        true); // only the edges owned by the current processor
    ERRORR(rval, "can't get shared edges owned");

    rval = parcomm->settle_intersection_points(RedEdges, shared_edges_owned, extraNodesVec, epsilon_1);
    ERRORR(rval, "can't settle intx points");
  }
  return MB_SUCCESS;
}
} /* namespace moab */
