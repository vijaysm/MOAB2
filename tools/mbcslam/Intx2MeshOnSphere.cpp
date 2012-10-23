/*
 * Intx2MeshOnSphere.cpp
 *
 *  Created on: Oct 3, 2012
 */

#include "Intx2MeshOnSphere.hpp"
#include "IntxUtils.hpp"

#include <queue>

namespace moab {


Intx2MeshOnSphere::Intx2MeshOnSphere(Interface * mbimpl):Intx2Mesh(mbimpl)
{
  // TODO Auto-generated constructor stub

}

Intx2MeshOnSphere::~Intx2MeshOnSphere()
{
  // TODO Auto-generated destructor stub
}

/* the red quad is convex for sure, intersect it with the blue quad
 * if the blue is convex, intersection is convex, it is pretty easy to order them and eliminate doubles
 *  first decide the projection plane, then do a gnomonic projection of both,
 *  compute intersection in the plane, then go back to the sphere for the points
 *  For the time being, blue is convex too, but later on, we should separate the case of concave blue
 */
int Intx2MeshOnSphere::computeIntersectionBetweenRedAndBlue(EntityHandle red, EntityHandle blue,
    double * P, int & nP, double & area, int markb[4], int markr[4])
{
  // the points will be at most 16; they will describe a convex patch, after the points will be ordered and
  // collapsed (eliminate doubles)
  // the area is not really required, except to see if it is greater than 0

  // gnomonic projection
  // int plane = 0;
  // get coordinates of the red quad, to decide the gnomonic plane

  int num_nodes;
  ErrorCode rval = mb->get_connectivity(red, redConn, num_nodes);
  if (MB_SUCCESS != rval || num_nodes != 4)
    return 1;
  //CartVect coords[4];
  rval = mb->get_coords(redConn, 4, &(redCoords[0][0]));
  if (MB_SUCCESS != rval)
    return 1;
  CartVect middle = 0.25
      * (redCoords[0] + redCoords[1] + redCoords[2] + redCoords[3]);

  decide_gnomonic_plane(middle, plane);// output the plane
  //CartVect bluecoords[4];
  rval = mb->get_connectivity(blue, blueConn, num_nodes);
  if (MB_SUCCESS != rval || num_nodes != 4)
    return 1;
  rval = mb->get_coords(blueConn, 4, &(blueCoords[0][0]));
  if (MB_SUCCESS != rval)
    return 1;

  if (dbg_1)
  {
    std::cout << "red " << mb->id_from_handle(red) << "\n";
    for (int j = 0; j < 4; j++)
    {
      std::cout << redCoords[j] << "\n";
    }
    std::cout << "blue " << mb->id_from_handle(blue) << "\n";
    for (int j = 0; j < 4; j++)
    {
      std::cout << blueCoords[j] << "\n";
    }
    mb->list_entities(&red, 1);
    mb->list_entities(&blue, 1);
    std::cout << "middle " << middle << "  plane:" << plane << "\n";
  }
  for (int j = 0; j < 4; j++)
  {
    // populate coords in the plane for intersection
    // they should be oriented correctly, positively
    int rc = gnomonic_projection(redCoords[j],  R, plane, redQuad[2 * j],
        redQuad[2 * j + 1]);
    if (rc != 0)
      return 1;
    rc = gnomonic_projection(blueCoords[j], R, plane, blueQuad[2 * j],
        blueQuad[2 * j + 1]);
    if (rc != 0)
      return 1;
  }
  if (dbg_1)
  {
    std::cout << "gnomonic plane: " << plane << "\n";
    std::cout << " red                                blue\n";
    for (int j = 0; j < 4; j++)
    {
      std::cout << redQuad[2 * j] << " " << redQuad[2 * j + 1] << " "
          << blueQuad[2 * j] << " " << blueQuad[2 * j + 1] << "\n";
    }
  }
  nP = 0; // number of intersection points we are marking the boundary of blue!
  int ret = EdgeIntersections2(blueQuad, redQuad, markb, markr, P, nP);
  if (ret != 0)
    return 1; // some unforeseen error

  int side[4] = { 0 };
  int extraPoints = borderPointsOfXinY2(blueQuad, redQuad, &(P[2 * nP]), side);
  if (extraPoints >= 1)
  {
    for (int k = 0; k < 4; k++)
    {
      if (side[k])
      {
        markb[k] = 1;
        markb[(k + 3) % 4] = 1; // it is the previous edge, actually, but instead of doing -1, it is
        // better to do modulo +3 (modulo 4)
        // null side b for next call
        side[k]=0;
      }
    }
  }
  nP += extraPoints;

  extraPoints = borderPointsOfXinY2(redQuad, blueQuad, &(P[2 * nP]), side);
  if (extraPoints >= 1)
  {
    for (int k = 0; k < 4; k++)
    {
      if (side[k])
      {
        markr[k] = 1;
        markr[(k + 3) % 4] = 1; // it is the previous edge, actually, but instead of doing -1, it is
        // better to do modulo +3 (modulo 4)
        // null side b for next call
      }
    }
  }
  nP += extraPoints;

  // now sort and orient the points in P, such that they are forming a convex polygon
  // this will be the foundation of our new mesh
  // this works if the polygons are convex
  SortAndRemoveDoubles2(P, nP, epsilon_1); // nP should be at most 8 in the end ?
  // if there are more than 3 points, some area will be positive
  area = 0.;
  if (nP >= 3)
  {
    for (int k = 1; k < nP - 1; k++)
      area += area2D(P, &P[2 * k], &P[2 * k + 2]);
  }

  return 0; // no error
}


// this method will also construct the triangles/polygons in the new mesh
// if we accept planar polygons, we just save them
// also, we could just create new vertices every time, and merge only in the end;
// could be too expensive, and the tolerance for merging could be an
// interesting topic
int Intx2MeshOnSphere::findNodes(EntityHandle red, EntityHandle blue, double * iP, int nP)
{
  // first of all, check against red and blue vertices
  //
  if (dbg_1)
  {
    std::cout << "red, blue, nP, P " << mb->id_from_handle(red) << " "
        << mb->id_from_handle(blue) << " " << nP << "\n";
    for (int n = 0; n < nP; n++)
      std::cout << " \t" << iP[2 * n] << "\t" << iP[2 * n + 1] << "\n";

  }

  // get the edges for the red triangle; the extra points will be on those edges, saved as
  // lists (unordered)
  EntityHandle redEdges[4];
  int i = 0;
  for (i = 0; i < 4; i++)
  {
    EntityHandle v[2] = { redConn[i], redConn[(i + 1) % 4] };
    std::vector<EntityHandle> adj_entities;
    ErrorCode rval = mb->get_adjacencies(v, 2, 1, false, adj_entities,
        Interface::INTERSECT);
    if (rval != MB_SUCCESS || adj_entities.size() < 1)
      return 0; // get out , big error
    redEdges[i] = adj_entities[0]; // should be only one edge between 2 nodes
  }
  // these will be in the new mesh, mbOut
  // some of them will be handles to the initial vertices from blue or red meshes (lagr or euler)

  EntityHandle * foundIds = new EntityHandle[nP];
  for (i = 0; i < nP; i++)
  {
    double * pp = &iP[2 * i]; // iP+2*i
    // project the point back on the sphere
    CartVect pos;
    reverse_gnomonic_projection(pp[0], pp[1], R, plane, pos);
    int found = 0;
    // first, are they on vertices from red or blue?
    // priority is the red mesh (mb2?)
    int j = 0;
    EntityHandle outNode = (EntityHandle) 0;
    for (j = 0; j < 4 && !found; j++)
    {
      //int node = redTri.v[j];
      double d2 = dist2(pp, &redQuad[2 * j]);
      if (d2 < epsilon_1)
      {

        foundIds[i] = redConn[j]; // no new node
        found = 1;
        if (dbg_1)
          std::cout << "  red node j:" << j << " id:"
              << mb->id_from_handle(redConn[j]) << " 2d coords:" << redCoords[2 * j] << "  "
              << redCoords[2 * j + 1] << " d2: " << d2 << " \n";
      }
    }

    for (j = 0; j < 4 && !found; j++)
    {
      //int node = blueTri.v[j];
      double d2 = dist2(pp, &blueQuad[2 * j]);
      if (d2 < epsilon_1)
      {
        // suspect is blueConn[j] corresponding in mbOut

        foundIds[i] = blueConn[j]; // no new node
        found = 1;
        if (dbg_1)
          std::cout << "  blue node " << j << " "
              << mb->id_from_handle(blueConn[j]) << " d2:" << d2 << " \n";
      }

    }
    if (!found)
    {
      // find the edge it belongs, first, on the red quad
      //
      for (j = 0; j < 4; j++)
      {
        int j1 = (j + 1) % 4;
        double area = area2D(&redQuad[2 * j], &redQuad[2 * j1], pp);
        if (dbg_1)
          std::cout << "   edge " << j << ": "
              << mb->id_from_handle(redEdges[j]) << " " << redConn[j] << " "
              << redConn[j1] << "  area : " << area << "\n";
        if (fabs(area) < epsilon_1)
        {
          // found the edge; now find if there is a point in the list here
          //std::vector<EntityHandle> * expts = extraNodesMap[redEdges[j]];
          int indx = -1;
          indx = RedEdges.index(redEdges[j]);
          std::vector<EntityHandle> * expts = extraNodesVec[indx];
          // if the points pp is between extra points, then just give that id
          // if not, create a new point, (check the id)
          // get the coordinates of the extra points so far
          int nbExtraNodesSoFar = expts->size();
          CartVect * coords1 = new CartVect[nbExtraNodesSoFar];
          mb->get_coords(&(*expts)[0], nbExtraNodesSoFar, &(coords1[0][0]));
          //std::list<int>::iterator it;
          for (int k = 0; k < nbExtraNodesSoFar && !found; k++)
          {
            //int pnt = *it;
            double d2 = (pos - coords1[k]).length_squared();
            if (d2 < epsilon_1)
            {
              found = 1;
              foundIds[i] = (*expts)[k];
              if (dbg_1)
                std::cout << " found node:" << foundIds[i] << std::endl;
            }
          }
          if (!found)
          {
            // create a new point in 2d (at the intersection)
            //foundIds[i] = m_num2dPoints;
            //expts.push_back(m_num2dPoints);
            // need to create a new node in mbOut
            // this will be on the edge, and it will be added to the local list
            mb->create_vertex(pos.array(), outNode);
            (*expts).push_back(outNode);
            foundIds[i] = outNode;
            found = 1;
            if (dbg_1)
              std::cout << " new node: " << outNode << std::endl;
          }
          delete[] coords1;
        }
      }
    }
    if (!found)
    {
      std::cout << " red quad: ";
      for (int j = 0; j < 4; j++)
      {
        std::cout << redQuad[2 * j] << " " << redQuad[2 * j + 1] << "\n";
      }
      std::cout << " a point pp is not on a red quad " << *pp << " " << pp[1]
          << " red quad " << mb->id_from_handle(red) << " \n";
      return 1;
    }
  }
  // now we can build the triangles, from P array, with foundIds
  // we will put them in the out set
  if (nP >= 3)
  {

    EntityHandle polyNew;
    mb->create_element(MBPOLYGON, foundIds, nP, polyNew);
    mb->add_entities(outSet, &polyNew, 1);

    // tag it with the ids from red and blue
    int id = mb->id_from_handle(blue);
    mb->tag_set_data(blueParentTag, &polyNew, 1, &id);
    id = mb->id_from_handle(red);
    mb->tag_set_data(redParentTag, &polyNew, 1, &id);

    static int count=0;
    count++;
    mb->tag_set_data(countTag, &polyNew, 1, &count);

    if (dbg_1)
    {

      std::cout << "Count: " << count+1 << "\n";
      std::cout << " polygon " << mb->id_from_handle(polyNew) << "  nodes: " << nP << " :";
      for (int i = 0; i < nP; i++)
        std::cout << " " << mb->id_from_handle(foundIds[i]);
      std::cout << " plane: " << plane << "\n";
      std::vector<CartVect> posi(nP);
      mb->get_coords(foundIds, nP, &(posi[0][0]));
      for (int i = 0; i < nP; i++)
        std::cout << iP[2 * i] << " " << iP[2 * i + 1] << " " << posi[i] << "\n";

      std::stringstream fff;
      fff << "file0" <<  count<< ".vtk";
          mb->write_mesh(fff.str().c_str(), &outSet, 1);
    }




  }
  delete[] foundIds;
  foundIds = NULL;
  return 0;
}

// main interface
ErrorCode Intx2MeshOnSphere::intersect_meshes(EntityHandle mbset1, EntityHandle mbset2,
    EntityHandle & outputSet)
{

  ErrorCode rval;
  mbs1 = mbset1; // is this the arrival or departure? I don't know yet
  mbs2 = mbset2;
  outSet = outputSet;

  // really, should be something from t1 and t2; blue is 1 (lagrange), red is 2 (euler)
  createTags(); //
  EntityHandle startBlue, startRed; // first triangles from mb1 and mb2
  //ErrorCode rval = mb1->handle_from_id(MBTRI, 1, startBlue);
  // we need to start somewhere; we will do an expensive search for one intersection
  //mb2->handle_from_id(MBTRI, 1, startRed);
  // this could be an expensive search
  // maybe we should do some KDtrees, for the worst case
  Range rs1;
  Range rs2;
  mb->get_entities_by_type(mbs1, MBQUAD, rs1);
  mb->get_entities_by_type(mbs2, MBQUAD, rs2);
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
      int nb[4], nr[4]; // sides
      computeIntersectionBetweenRedAndBlue(startRed, startBlue, P, nP, area, nb, nr);
      if (area > 0)
      {
        found = 1;
        break; // found 2 quads that intersect; these will be the seeds
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
  int k;

  if (dbg_1)
  {
    mout_1[0].open("patches1.m");
    mout_1[1].open("patches2.m");
    mout_1[2].open("patches3.m");
    mout_1[3].open("patches4.m");
    mout_1[4].open("patches5.m");
    mout_1[5].open("patches6.m");
  }
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
        if (dbg_1)
        {
          std::cout << "gnomonic plane: " << plane << "\n";
          std::cout << "area: " << area << " nP:" << nP << std::endl;
          std::cout << "nb: " << nb[0] << nb[1] << nb[2] << nb[3] << "\n";
          std::cout << "nr: " << nr[0] << nr[1] << nr[2] << nr[3] << "\n";

          mout_1[plane - 1] << "pa=[\n";

          for (k = 0; k < nP; k++)
          {

            mout_1[plane - 1] << P[2 * k] << "\t ";
          }

          mout_1[plane - 1] << "\n";
          for (k = 0; k < nP; k++)
          {

            mout_1[plane - 1] << P[2 * k + 1] << "\t ";
          }

          mout_1[plane - 1] << " ]; \n";
          mout_1[plane - 1] << " patch(pa(1,:),pa(2,:),'m');       \n";
          mout_1[plane - 1] << " pause(1);\n";
        }
        EntityHandle neighbors[4];
        rval = GetOrderedNeighbors(mbs2, redT, neighbors);
        if (rval != MB_SUCCESS)
        {
          std::cout << " can't get the neighbors for red quad "
              << mb->id_from_handle(redT);
          return MB_FAILURE;
        }

        if (dbg_1)
        {
          std::cout << " neighbors for redT " << mb->id_from_handle(redT)
              << " \n";
          for (int kk = 0; kk < 4; kk++)
          {
            if (neighbors[kk] > 0)
              std::cout << mb->id_from_handle(neighbors[kk]) << " ";
            else
              std::cout << 0 << " ";
          }
          std::cout << std::endl;
          //mb->list_entities(neighbors, 4);
        }
        // add neighbors to the localRed queue, if they are not marked
        for (int nn = 0; nn < 4; nn++)
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
                std::cout << " local red quad " << mb->id_from_handle(neighbor)
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
    if (dbg_1)
    {
      std::cout << "Next: neighbors for blue T ";
      for (int kk = 0; kk < 4; kk++)
      {
        if (blueNeighbors[kk] > 0)
          std::cout << mb->id_from_handle(blueNeighbors[kk]) << " ";
        else
          std::cout << 0 << " ";
      }
      std::cout << std::endl;
    }
    for (int j = 0; j < 4; j++)
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

} /* namespace moab */
