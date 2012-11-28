/*
 * Intx2Mesh.cpp
 *
 *  Created on: Oct 2, 2012
 */

#include "Intx2Mesh.hpp"
#include "Coupler.hpp"
#include "MBParallelConventions.h"
#include <queue>

namespace moab {

#define ERRORR(rval, str) \
    if (MB_SUCCESS != rval) {std::cout << str << "\n"; return rval;}

#define ERRORV(rval, str) \
    if (MB_SUCCESS != rval) {std::cout << str << "\n"; return ;}

Intx2Mesh::Intx2Mesh(Interface * mbimpl):mb(mbimpl), parcomm(NULL)
{
  dbg_1=0;
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

// this will work in parallel only
ErrorCode Intx2Mesh::locate_departure_points(EntityHandle euler_set)
{
  // get the par comm ...
  parcomm = ParallelComm::get_pcomm(mb, 0);
  if (NULL==parcomm)
    return MB_FAILURE;
  // get the DP tag, and values for each node in the euler set
  // get forst all entities in the set
  Range local_ents;
  ErrorCode rval = mb->get_entities_by_dimension(euler_set, 2, local_ents);
  ERRORR(rval, "can't get ents by dimension");

  // instantiate a coupler, which also initializes the tree
  Coupler mbc(mb, parcomm, local_ents, 0, true, 2);// init tree = true, and max_dim is 2

  // get all the owned nodes on this euler set
  Range local_verts;
  rval = mb->get_connectivity(local_ents, local_verts);
  if (parcomm)
  {
    // get only those that are owned by this proc
    rval = parcomm->filter_pstatus(local_verts, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    ERRORR(rval, "can't filter verts");
  }

  // get the DP tag , and get the departure points
  Tag tagh = 0;
  std::string tag_name("DP");
  rval = mb->tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, tagh, MB_TAG_DENSE);
  ERRORR(rval, "can't get DP tag");
  int num_owned_verts = (int)local_verts.size();

  std::vector<double> dep_points(3*num_owned_verts);
  rval = mb->tag_get_data(tagh, local_verts, (void*)&dep_points[0]);
  ERRORR(rval, "can't get DP tag values");

  mbc.locate_points(&dep_points[0], num_owned_verts, 0, 0.2);


  return MB_SUCCESS;
}
} /* namespace moab */
