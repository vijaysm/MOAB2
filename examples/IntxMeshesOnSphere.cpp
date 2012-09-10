/*
 * IntxMeshesOnSphere.cpp
 *
 *  This test is for intersection of 2d meshes on the same sphere; the meshes are formed by quads
 *
 *  inputs are 2 meshes, vtk format, output will be another mesh, m3, with the common part
 *    intersected
 *
 *  Created on: Sep 2, 2012
 *      Author: iulian
 */

#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"

using namespace std;
using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#include <queue>
//#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

// we should maybe check that the meshes are correct, completely
// covering spheres and not leaving gaps
EntityHandle mbs1; // set 1, blue triangles
EntityHandle mbs2;
EntityHandle outSet;
Interface * mb; // global
double R; // radius; should be input to the problem;

Tag BlueFlagTag; // to mark blue quads already considered

Tag RedFlagTag; // to mark red quads already considered

Range RedEdges; //
// for each red edge, we keep a vector of extra nodes, coming from intersections
// use the index in RedEdges range, instead of a map, as before
// std::map<EntityHandle, std::vector<EntityHandle> *> extraNodesMap;
std::vector<std::vector<EntityHandle> *> extraNodesVec;

// red parent and blue parent tags
// these will be on the out mesh
Tag redParentTag;
Tag blueParentTag;
Tag countTag;

// some vars for current process
// after intersection, these var are still needed by findNodes, to establish the
// connections with the mesh
const EntityHandle * redConn;
const EntityHandle * blueConn;
CartVect redCoords[4];
CartVect blueCoords[4];
double redQuad[8]; // these are in plane
double blueQuad[8]; // these are in plane
int plane; // store the current gnomonic plane

int dbg_1 = 0;
ofstream mout_1[6]; // some debug files
#define epsilon_1 R*1.e-8 // cm, for coincident points in P, the intersection area
// will create the red edges too, if not existing yet, and some tags
// the nodes are unique now,
void createTags()
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
  if (MB_SUCCESS != rval)
    return;

  // assume that the edges are on the red triangles
  Range redQuads;
  //Range redEdges;
  rval = mb->get_entities_by_type(mbs2, MBQUAD, redQuads, false);
  if (MB_SUCCESS != rval)
    return;
  // create red edges if they do not exist yet; so when they are looked upon, they are found
  // this is the only call that is potentially NlogN, in the whole method
  rval = mb->get_adjacencies(redQuads, 1, true, RedEdges, Interface::UNION);

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
  rval = mb->tag_get_handle("Negative", 1, MB_TYPE_INTEGER, blueParentTag,
      MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);

  rval = mb->tag_get_handle("Counting", 1, MB_TYPE_INTEGER, countTag,
        MB_TAG_SPARSE | MB_TAG_CREAT, &defaultInt);

  return;
}
// vec utils

/*
 *
 * position on a sphere of radius R
 * if plane specified, use it; if not, return the plane, and the point in the plane
 * there are 6 planes, numbered 1 to 6
 * plane 1: x=R, plane 2: y=R, 3: x=-R, 4: y=-R, 5: z=-R, 6: z=R
 *
 * projection on the plane will preserve the orientation, such that a triangle, quad pointing
 * outside the sphere will have a positive orientation in the projection plane
 */
void decide_gnomonic_plane(const CartVect & pos)
{
  // decide plane, based on max x, y, z
  if (fabs(pos[0]) < fabs(pos[1]))
  {
    if (fabs(pos[2]) < fabs(pos[1]))
    {
      // pos[1] is biggest
      if (pos[1] > 0)
        plane = 2;
      else
        plane = 4;
    }
    else
    {
      // pos[2] is biggest
      if (pos[2] < 0)
        plane = 5;
      else
        plane = 6;
    }
  }
  else
  {
    if (fabs(pos[2]) < fabs(pos[0]))
    {
      // pos[0] is the greatest
      if (pos[0] > 0)
        plane = 1;
      else
        plane = 3;
    }
    else
    {
      // pos[2] is biggest
      if (pos[2] < 0)
        plane = 5;
      else
        plane = 6;
    }
  }
  return;
}
// point on a sphere is projected on one of six planes, decided earlier
int gnomonic_projection(const CartVect & pos, double & c1, double & c2)
{
  double alfa = 1.; // the new point will be on line alfa*pos

  switch (plane)
  {
  case 1:
  {
    // the plane with x = R; x>y, x>z
    // c1->y, c2->z
    alfa = R / pos[0];
    c1 = alfa * pos[1];
    c2 = alfa * pos[2];
    break;
  }
  case 2:
  {
    // y = R -> zx
    alfa = R / pos[1];
    c1 = alfa * pos[2];
    c2 = alfa * pos[0];
    break;
  }
  case 3:
  {
    // x=-R, -> yz
    alfa = -R / pos[0];
    c1 = -alfa * pos[1]; // the sign is to preserve orientation
    c2 = alfa * pos[2];
    break;
  }
  case 4:
  {
    // y = -R
    alfa = -R / pos[1];
    c1 = -alfa * pos[2]; // the sign is to preserve orientation
    c2 = alfa * pos[0];
    break;
  }
  case 5:
  {
    // z = -R
    alfa = -R / pos[2];
    c1 = -alfa * pos[0]; // the sign is to preserve orientation
    c2 = alfa * pos[1];
    break;
  }
  case 6:
  {
    alfa = R / pos[2];
    c1 = alfa * pos[0];
    c2 = alfa * pos[1];
    break;
  }
  default:
    return 1; // error
  }

  return 0; // no error
}
// given the position on plane (one out of 6), find out the position on sphere
int reverse_gnomonic_projection(const double & c1, const double & c2,
    CartVect & pos)
{

  // the new point will be on line beta*pos
  double len = sqrt(c1 * c1 + c2 * c2 + R * R);
  double beta = R / len; // it is less than 1, in general

  switch (plane)
  {
  case 1:
  {
    // the plane with x = R; x>y, x>z
    // c1->y, c2->z
    pos[0] = beta * R;
    pos[1] = c1 * beta;
    pos[2] = c2 * beta;
    break;
  }
  case 2:
  {
    // y = R -> zx
    pos[1] = R * beta;
    pos[2] = c1 * beta;
    pos[0] = c2 * beta;
    break;
  }
  case 3:
  {
    // x=-R, -> yz
    pos[0] = -R * beta;
    pos[1] = -c1 * beta; // the sign is to preserve orientation
    pos[2] = c2 * beta;
    break;
  }
  case 4:
  {
    // y = -R
    pos[1] = -R * beta;
    pos[2] = -c1 * beta; // the sign is to preserve orientation
    pos[0] = c2 * beta;
    break;
  }
  case 5:
  {
    // z = -R
    pos[2] = -R * beta;
    pos[0] = -c1 * beta; // the sign is to preserve orientation
    pos[1] = c2 * beta;
    break;
  }
  case 6:
  {
    pos[2] = R * beta;
    pos[0] = c1 * beta;
    pos[1] = c2 * beta;
    break;
  }
  default:
    return 1; // error
  }

  return 0; // no error
}
double dist2(double * a, double * b)
{
  double abx = b[0] - a[0], aby = b[1] - a[1];
  return sqrt(abx * abx + aby * aby);
}
double area2D(double *a, double *b, double *c)
{
  // (b-a)x(c-a) / 2
  return ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2;
}
int borderPointsOfXinY2(double * X, double * Y, double * P, int side[4])
{
  // 2 triangles, 3 corners, is the corner of X in Y?
  // Y must have a positive area
  /*
   */
  int extraPoint = 0;
  for (int i = 0; i < 4; i++)
  {
    // compute twice the area of all 4 triangles formed by a side of Y and a corner of X; if one is negative, stop
    double * A = X + 2 * i;

    int inside = 1;
    for (int j = 0; j < 4; j++)
    {
      double * B = Y + 2 * j;

      int j1 = (j + 1) % 4;
      double * C = Y + 2 * j1; // no copy of data

      double area2 = (B[0] - A[0]) * (C[1] - A[1])
          - (C[0] - A[0]) * (B[1] - A[1]);
      if (area2 < 0.)
      {
        inside = 0;
        break;
      }
    }
    if (inside)
    {
      side[i] = 1;
      P[extraPoint * 2] = A[0];
      P[extraPoint * 2 + 1] = A[1];
      extraPoint++;
    }
  }
  return extraPoint;
}

int swap2(double * p, double * q)
{
  double tmp = *p;
  *p = *q;
  *q = tmp;
  return 0;
}
int SortAndRemoveDoubles2(double * P, int & nP)
{
  if (nP < 2)
    return 0; // nothing to do
  // center of gravity for the points
  double c[2] = { 0., 0. };
  int k = 0;
  for (k = 0; k < nP; k++)
  {
    c[0] += P[2 * k];
    c[1] += P[2 * k + 1];
  }
  c[0] /= nP;
  c[1] /= nP;
  double angle[24]; // could be at most 24 points; much less usually
  for (k = 0; k < nP; k++)
  {
    double x = P[2 * k] - c[0], y = P[2 * k + 1] - c[1];
    if (x != 0. || y != 0.)
      angle[k] = atan2(y, x);
    else
    {
      angle[k] = 0;
      // this means that the points are on a line, or all coincident // degenerate case
    }
  }
  // sort according to angle; also eliminate close points
  // this is a bubble sort, for np = 24 it could be pretty big
  // maybe a better sort is needed here (qsort?)
  int sorted = 1;
  do
  {
    sorted = 1;
    for (k = 0; k < nP - 1; k++)
    {
      if (angle[k] > angle[k + 1])
      {
        sorted = 0;
        swap2(angle + k, angle + k + 1);
        swap2(P + (2 * k), P + (2 * k + 2));
        swap2(P + (2 * k + 1), P + (2 * k + 3));
      }
    }
  } while (!sorted);
  // eliminate doubles

  int i = 0, j = 1; // the next one; j may advance faster than i
  // check the unit
  // double epsilon_1 = 1.e-5; // these are cm; 2 points are the same if the distance is less than 1.e-5 cm
  while (j < nP)
  {
    double d2 = dist2(&P[2 * i], &P[2 * j]);
    if (d2 > epsilon_1)
    {
      i++;
      P[2 * i] = P[2 * j];
      P[2 * i + 1] = P[2 * j + 1];
    }
    j++;
  }
  // test also the last point with the first one (index 0)

  double d2 = dist2(P, &P[2 * i]); // check the first and last points (ordered from -pi to +pi)
  if (d2 > epsilon_1)
  {
    nP = i + 1;
  }
  else
    nP = i; // effectively delete the last point (that would have been the same with first)
  if (nP == 0)
    nP = 1; // we should be left with at least one point we already tested if nP is 0 originally
  return 0;
}

// the marks will show what edges of blue intersect the red

int EdgeIntersections2(double * blue, double * red, int markb[4], int markr[4],
    double * points, int & nPoints)
{
  /* EDGEINTERSECTIONS computes edge intersections of two triangles
   [P,n]=EdgeIntersections(X,Y) computes for the two given triangles  * red
   and blue ( stored column wise )
   (point coordinates are stored column-wise, in counter clock
   order) the points P where their edges intersect. In addition,
   in n the indices of which neighbors of red  are also intersecting
   with blue are given.
   */

  // points is an array with 48 slots   (24 * 2 doubles)
  nPoints = 0;
  markb[0] = markb[1] = markb[2] = markb[3] = 0; // no neighbors of red involved yet
  markr[0] = markr[1] = markr[2] = markr[3] = 0;
  /*for i=1:3                            % find all intersections of edges
   for j=1:3
   b=Y(:,j)-X(:,i);
   A=[X(:,mod(i,3)+1)-X(:,i) -Y(:,mod(j,3)+1)+Y(:,j)];
   if rank(A)==2                   % edges not parallel
   r=A\b;
   if r(1)>=0 & r(1)<=1 & r(2)>=0 & r(2)<=1,  % intersection found
   k=k+1; P(:,k)=X(:,i)+r(1)*(X(:,mod(i,3)+1)-X(:,i)); n(i)=1;
   end;
   end;
   end;
   end;*/
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      double b[2];
      double a[2][2]; // 2*2
      int iPlus1 = (i + 1) % 4;
      int jPlus1 = (j + 1) % 4;
      for (int k = 0; k < 2; k++)
      {
        b[k] = red[2 * j + k] - blue[2 * i + k];
        // row k of a: a(k, 0), a(k, 1)
        a[k][0] = blue[2 * iPlus1 + k] - blue[2 * i + k];
        a[k][1] = red[2 * j + k] - red[2 * jPlus1 + k];

      }
      double delta = a[0][0] * a[1][1] - a[0][1] * a[1][0];
      if (delta != 0.)
      {
        // not parallel
        double alfa = (b[0] * a[1][1] - a[0][1] * b[1]) / delta;
        double beta = (-b[0] * a[1][0] + b[1] * a[0][0]) / delta;
        if (0 <= alfa && alfa <= 1. && 0 <= beta && beta <= 1.)
        {
          // the intersection is good
          for (int k = 0; k < 2; k++)
          {
            points[2 * nPoints + k] = blue[2 * i + k]
                + alfa * (blue[2 * iPlus1 + k] - blue[2 * i + k]);
          }
          markb[i] = 1; // so neighbor number i of blue will be considered too.
          markr[j] = 1; // this will be used in advancing red around blue quad
          nPoints++;
        }
      }
      // the case delta == 0. will be considered by the interior points logic

    }
  }
  return 0;
}

/* the red quad is convex for sure, intersect it with the blue quad
 * if the blue is convex, intersection is convex, it is pretty easy to order them and eliminate doubles
 *  first decide the projection plane, then do a gnomonic projection of both,
 *  compute intersection in the plane, then go back to the sphere for the points
 *  For the time being, blue is convex too, but later on, we should separate the case of concave blue
 */
int computeIntersectionBetweenRedAndBlue(EntityHandle red, EntityHandle blue,
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

  decide_gnomonic_plane(middle);
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
    int rc = gnomonic_projection(redCoords[j], redQuad[2 * j],
        redQuad[2 * j + 1]);
    if (rc != 0)
      return 1;
    rc = gnomonic_projection(blueCoords[j], blueQuad[2 * j],
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
  SortAndRemoveDoubles2(P, nP); // nP should be at most 8 in the end ?
  // if there are more than 3 points, some area will be positive
  area = 0.;
  if (nP >= 3)
  {
    for (int k = 1; k < nP - 1; k++)
      area += area2D(P, &P[2 * k], &P[2 * k + 2]);
  }

  return 0; // no error
}

// specify also desired set; we are interested only in neighbors in the set!
// we should always get manifold mesh, each edge is adjacent to 2 quads
ErrorCode GetOrderedNeighbors(EntityHandle set, EntityHandle quad,
    EntityHandle neighbors[4])
{
  // will get the 4 ordered neighbors;
  // first quad is for nodes 0, 1, second to 1, 2, third to 2, 3,
  int nnodes;
  const EntityHandle * conn4;
  ErrorCode rval = mb->get_connectivity(quad, conn4, nnodes);
  if (MB_SUCCESS != rval || nnodes != 4)
    return MB_FAILURE;
  for (int i = 0; i < 4; i++)
  {
    EntityHandle v[2];
    v[0] = conn4[i];
    v[1] = conn4[(i + 1) % 4];
    // get quads adjacent to vertices
    std::vector<EntityHandle> quads;
    std::vector<EntityHandle> quadsInSet;
    rval = mb->get_adjacencies(v, 2, 2, false, quads, Interface::INTERSECT);
    if (MB_SUCCESS != rval)
      return rval;
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

// this method will also construct the triangles in the new mesh
// if we accept planar polygons, we just save them
// also, we could just create new vertices every time, and merge only in the end;
// could be too expensive, and the tolerance for merging could be an
// interesting topic
int findNodes(EntityHandle red, EntityHandle blue, double * iP, int nP)
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
    reverse_gnomonic_projection(pp[0], pp[1], pos);
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
        cout << iP[2 * i] << " " << iP[2 * i + 1] << " " << posi[i] << "\n";

      std::stringstream fff;
      fff << "file0" <<  count<< ".vtk";
          mb->write_mesh(fff.str().c_str(), &outSet, 1);
    }




  }
  delete[] foundIds;
  foundIds = NULL;
  return 0;
}

// clean some memory allocated
void clean()
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
// mbs1 is the set for lagrangian, mbs2 is for euler (arrival)
ErrorCode intersect_meshes(EntityHandle mbs1, EntityHandle mbs2,
    EntityHandle & outputSet)
{

  ErrorCode rval = MB_SUCCESS;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  if (MB_SUCCESS != rval)
    return rval;

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

int main(int argc, char* argv[])
{
  // check command line arg// Euler grid is red, arrival, Lagrangian is blue, departure
  // will will keep the
  const char *filename_mesh1 = STRINGIFY(SRCDIR) "/lagrangeHomme.vtk";
  const char *filename_mesh2 = STRINGIFY(SRCDIR) "/eulerHomme.vtk";
  R = 6. * sqrt(3.) / 2; // input
  const char *newFile = "intx.vtk";
  if (argc == 5)
  {
    filename_mesh1 = argv[1];
    filename_mesh2 = argv[2];
    R = atof(argv[3]);
    newFile = argv[4];
  }
  else
  {
    printf("Usage: %s <mesh_filename1> <mesh_filename2> <radius>  <newFile>\n",
        argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s  %s  %f %s\n",
        filename_mesh1, filename_mesh2, R, newFile);
  }

  // read meshes in 2 file sets
  ErrorCode rval = MB_SUCCESS;
  Core moab;
  /*Interface**/mb = &moab; // global
  EntityHandle sf1, sf2;
  rval = mb->create_meshset(MESHSET_SET, sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->create_meshset(MESHSET_SET, sf2);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->load_file(filename_mesh1, &sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->load_file(filename_mesh2, &sf2);
  if (MB_SUCCESS != rval)
    return 1;

  EntityHandle outputSet;
  rval = intersect_meshes(sf1, sf2, outputSet);
  if (MB_SUCCESS != rval)
    std::cout << " failed to intersect meshes\n";
  rval = mb->write_mesh(newFile, &outputSet, 1);
  if (MB_SUCCESS != rval)
    return 1;
  return 0;

}
