/*
 * Intx2MeshInPlane.cpp
 *
 *  Created on: Oct 24, 2012
 *      Author: iulian
 */

#include "Intx2MeshInPlane.hpp"
#include "moab/GeomUtil.hpp"
namespace moab {
Intx2MeshInPlane::Intx2MeshInPlane(Interface * mbimpl):Intx2Mesh(mbimpl){

}

Intx2MeshInPlane::~Intx2MeshInPlane() {
  // TODO Auto-generated destructor stub
}

int Intx2MeshInPlane::computeIntersectionBetweenRedAndBlue(EntityHandle red, EntityHandle blue,
    double * P, int & nP, double & area, int markb[4], int markr[4], bool check_boxes_first)
{

   // the points will be at most 9; they will describe a convex patch, after the points will be ordered and
   // collapsed (eliminate doubles)
   // the area is not really required
   int nsides = 3;
   if (type == MBQUAD)
     nsides = 4;
   int num_nodes;
   ErrorCode rval = mb->get_connectivity(red, redConn, num_nodes);

   if (MB_SUCCESS != rval || num_nodes != nsides)
     return 1;

   //CartVect coords[4];
   rval = mb->get_coords(redConn, num_nodes, &(redCoords[0][0]));
   if (MB_SUCCESS != rval)
     return 1;


   rval = mb->get_connectivity(blue, blueConn, num_nodes);
   if (MB_SUCCESS != rval || num_nodes != nsides)
     return 1;
   rval = mb->get_coords(blueConn, num_nodes, &(blueCoords[0][0]));
   if (MB_SUCCESS != rval)
     return 1;

   if (dbg_1)
   {
     std::cout << "red " << mb->id_from_handle(red) << "\n";
     for (int j = 0; j < num_nodes; j++)
     {
       std::cout << redCoords[j] << "\n";
     }
     std::cout << "blue " << mb->id_from_handle(blue) << "\n";
     for (int j = 0; j < num_nodes; j++)
     {
       std::cout << blueCoords[j] << "\n";
     }
     mb->list_entities(&red, 1);
     mb->list_entities(&blue, 1);
   }
   area = 0.;
   nP = 0; // number of intersection points we are marking the boundary of blue!
   if (check_boxes_first)
   {
     // look at the boxes formed with vertices; if they are far away, return false early
     if (!GeomUtil::bounding_boxes_overlap(redCoords, num_nodes, blueCoords, num_nodes, box_error))
       return 0; // no error, but no intersection, decide early to get out
   }
   for (int j = 0; j < nsides; j++)
   {
     // populate coords in the plane for intersection
     // they should be oriented correctly, positively
     redQuad[2 * j]=redCoords[j][0]; // x coordinate,
     redQuad[2 * j + 1] = redCoords[j][1]; // y coordinate
     blueQuad[2 * j]=blueCoords[j][0]; // x coordinate,
     blueQuad[2 * j + 1] = blueCoords[j][1]; // y coordinate
   }
   if (dbg_1)
   {
     //std::cout << "gnomonic plane: " << plane << "\n";
     std::cout << " red                                blue\n";
     for (int j = 0; j < 4; j++)
     {
       std::cout << redQuad[2 * j] << " " << redQuad[2 * j + 1] << " "
           << blueQuad[2 * j] << " " << blueQuad[2 * j + 1] << "\n";
     }
   }

   nP = 0; // number of intersection points
   int ret = EdgeIntersections2(blueQuad, redQuad, nsides, markb,
       markr, P, nP);
   if (ret != 0)
      return 1;// some unforeseen error

   // start copy
   int side[4] = { 0 };
   int extraPoints = borderPointsOfXinY2(blueQuad, redQuad, nsides,
       &(P[2 * nP]), side);
   if (extraPoints >= 1)
   {
     for (int k = 0; k < nsides; k++)
     {
       if (side[k])
       {
         markb[k] = 1;
         markb[(k + nsides-1) % nsides] = 1; // it is the previous edge, actually, but instead of doing -1, it is
         // better to do modulo +3 (modulo 4)
         // null side b for next call
         side[k]=0;
       }
     }
   }
   nP += extraPoints;

   extraPoints = borderPointsOfXinY2(redQuad, blueQuad, nsides, &(P[2 * nP]), side);
   if (extraPoints >= 1)
   {
     for (int k = 0; k < 4; k++)
     {
       if (side[k])
       {
         markr[k] = 1;
         markr[(k + nsides-1) % nsides] = 1; // it is the previous edge, actually, but instead of doing -1, it is
         // better to do modulo +3 (modulo 4)
         // null side b for next call
       }
     }
   }
   nP += extraPoints;

   //
   SortAndRemoveDoubles2(P, nP, epsilon_1);  // nP should be at most 6 in the end
   // if there are more than 3 points, some area will be positive
   area = 0.;
   if (nP >= 3) {
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
int Intx2MeshInPlane::findNodes(EntityHandle red, EntityHandle blue, double * iP, int nP)
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

  int nsides = 3;
  if (type == MBQUAD)
    nsides = 4;
  // get the edges for the red triangle; the extra points will be on those edges, saved as
  // lists (unordered)
  EntityHandle redEdges[4];
  int i = 0;
  for (i = 0; i < nsides; i++)
  {
    EntityHandle v[2] = { redConn[i], redConn[(i + 1) % nsides] };
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
    CartVect pos(pp[0], pp[1], 0.); // constructor with 3 doubles
    //reverse_gnomonic_projection(pp[0], pp[1], R, plane, pos);
    int found = 0;
    // first, are they on vertices from red or blue?
    // priority is the red mesh (mb2?)
    int j = 0;
    EntityHandle outNode = (EntityHandle) 0;
    for (j = 0; j < nsides && !found; j++)
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

    for (j = 0; j < nsides && !found; j++)
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
      for (j = 0; j < nsides; j++)
      {
        int j1 = (j + 1) % nsides;
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
          delete[] coords1;

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

        }
      }
    }
    if (!found)
    {
      std::cout << " red quad: ";
      for (int j1 = 0; j1 < 4; j1++)
      {
        std::cout << redQuad[2 * j1] << " " << redQuad[2 * j1 + 1] << "\n";
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
      for (int i1 = 0; i1 < nP; i1++)
        std::cout << " " << mb->id_from_handle(foundIds[i1]);
      //std::cout << " plane: " << plane << "\n";
      std::vector<CartVect> posi(nP);
      mb->get_coords(foundIds, nP, &(posi[0][0]));
      for (int i1 = 0; i1 < nP; i1++)
        std::cout << iP[2 * i1] << " " << iP[2 * i1 + 1] << " " << posi[i1] << "\n";

      std::stringstream fff;
      fff << "file0" <<  count<< ".vtk";
          mb->write_mesh(fff.str().c_str(), &outSet, 1);
    }




  }
  delete[] foundIds;
  foundIds = NULL;
  return 0;
}
bool Intx2MeshInPlane::is_inside_element(double xyz[3], EntityHandle eh)
{
  int num_nodes;
  ErrorCode rval = mb->get_connectivity(eh, redConn, num_nodes);

  if (MB_SUCCESS != rval || num_nodes != 4 || num_nodes!=3)
    return false;
  int nsides = num_nodes;

  //CartVect coords[4];
  rval = mb->get_coords(redConn, num_nodes, &(redCoords[0][0]));
  if (MB_SUCCESS != rval)
    return 1;

  for (int j = 0; j < nsides; j++)
  {
    // populate coords in the plane for decision making
    // they should be oriented correctly, positively
    redQuad[2 * j]     = redCoords[j][0];
    redQuad[2 * j + 1] = redCoords[j][1];
  }

  double pt[2]={xyz[0], xyz[1]};// xy plane only
  // now, is the projected point inside the red quad?
  // cslam utils
  if (point_in_interior_of_convex_polygon (redQuad, num_nodes, pt))
    return true;
  return false;
}
} // end namespace moab
