/*
 * IntersectMesh_test.cpp
 *
 *  This test is for intersection of 2d meshes, on top of each other, for
 *   a domain decomposition algorithm; it leverages the method from ProjectShell algorithm.
 *
 *  inputs are 2 meshes, vtk format, output will be another mesh, m3, with the common part
 *    intersected
 *
 *  Created on: Aug 19, 2010
 *      Author: iulian
 */

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

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

// we should really work only with triangles in their sets;
// can we get adjacencies limited to a set?
// how do we check if a triangle is really part of a set?
EntityHandle mbs1;// set 1, blue triangles
EntityHandle mbs2;
EntityHandle outSet;
Interface * mb;// global

Tag BlueFlagTag; // to mark blue triangles already considered

Tag RedFlagTag; // to mark blue triangles already considered

Tag RedNodeTag; // value will be the node in new mbOut mb2; it will
// be used to mark the corresponding node in mbOut

Tag BlueNodeTag; // value will be the node handle in mbOut mb; it will
// be used to mark the corresponding node in mbOut

Range RedEdges;//
// for each red edge, we keep a vector of extra nodes, coming from intersections
// use the index in RedEdges range, instead of a map, as before
// std::map<EntityHandle, std::vector<EntityHandle> *> extraNodesMap;
std::vector< std::vector<EntityHandle> *> extraNodesVec;

// red parent and blue parent tags
// these will be on the out mesh
Tag redParentTag;
Tag blueParentTag;

int dbg_1 = 1;
ofstream mout_1;// some debug file
double epsilon_1 = 1.e-5; // cm, for coincident points in P, the intersection area

void createTags(Interface * mb) {
   unsigned char def_data_bit = 0;// unused by default
   ErrorCode rval = mb->tag_get_handle("blueFlag", 1, MB_TYPE_BIT, BlueFlagTag,
        MB_TAG_CREAT, &def_data_bit);
   if (MB_SUCCESS != rval)
      return;
   // maybe the red tag is better to be deleted every time, and recreated;
   // or is it easy to set all values to something again? like 0?
   rval = mb->tag_get_handle("redFlag", 1, MB_TYPE_BIT, RedFlagTag, MB_TAG_CREAT, &def_data_bit);
   if (MB_SUCCESS != rval)
      return;
   // create a tag in mb for eh nodes in mbOut
   const EntityHandle def_val = (EntityHandle) 0;
   //rval = m2->tag_create("redNodeTag", sizeof(EntityHandle), MB_TAG_DENSE,
   //      MB_TYPE_HANDLE, RedNodeTag, &def_val);
   rval = mb->tag_get_handle("redNodeTag", 1, MB_TYPE_HANDLE, RedNodeTag,
                             MB_TAG_DENSE|MB_TAG_CREAT, &def_val);
   if (MB_SUCCESS != rval)
      return;
   //rval = m1->tag_create("blueNodeTag", 1, MB_TYPE_HANDLE, BlueNodeTag,
   //                      MB_TAG_DENSE|MB_TAG_CREAT, &def_val);
   rval = mb->tag_get_handle("blueNodeTag", 1, MB_TYPE_HANDLE, BlueNodeTag,
                             MB_TAG_DENSE|MB_TAG_CREAT, &def_val);
   if (MB_SUCCESS != rval)
      return;
   // assume that the edges are on the red triangles
   Range redTriangles;
   //Range redEdges;
   rval = mb->get_entities_by_type(mbs2, MBTRI, redTriangles, false);
   if (MB_SUCCESS != rval)
      return;
   // create red edges if they do not exist yet; so when they are looked upon, they are found
   // this is the only call that is potentially NlogN, in the whole method
   rval = mb ->get_adjacencies(redTriangles, 1, true, RedEdges,
         Interface::UNION);

   // now, create a map from each edge to a list of potential new nodes on a red edge
   // this memory has to be cleaned up
   // change it to a vector, and use the index in range of red edges
   int indx = 0;
   extraNodesVec.reserve(RedEdges.size());
   for (Range::iterator eit = RedEdges.begin(); eit != RedEdges.end(); eit++, indx++) {
      //EntityHandle edge = *eit;
      //extraNodesMap[edge] = new std::vector<EntityHandle>;
      std::vector<EntityHandle> * nv = new std::vector<EntityHandle>;
      extraNodesVec.push_back(nv);
   }

   int defaultInt = 0;

   rval = mb->tag_get_handle( "Positive", 1, MB_TYPE_INTEGER, redParentTag,
                               MB_TAG_SPARSE|MB_TAG_CREAT, &defaultInt );
   rval = mb->tag_get_handle( "Negative", 1, MB_TYPE_INTEGER, blueParentTag,
                               MB_TAG_SPARSE|MB_TAG_CREAT, &defaultInt );

   return;
}
// vec utils

double dist2(double * a, double * b)
{
  double abx=b[0]-a[0], aby=b[1]-a[1];
  return sqrt( abx*abx+aby*aby);
}
double area2D(double *a, double *b, double *c)
{
  // (b-a)x(c-a) / 2
  return  ( (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) )/2;
}
int borderPointsOfXinY2(double * X, double * Y, double * P) {
   // 2 triangles, 3 corners, is the corner of X in Y?
   // Y must have a positive area
   /*
    function P=PointsOfXInY(X,Y);
    % POINTSOFXINY finds corners of one triangle within another one
    %   P=PointsOfXInY(X,Y); computes for the two given triangles X
    %   and Y (point coordinates are stored column-wise, in counter clock
    %   order) the corners P of X which lie in the interior of Y.

    k=0;P=[];
    v0=Y(:,2)-Y(:,1); v1=Y(:,3)-Y(:,1);  % find interior points of X in Y
    d00=v0'*v0; d01=v0'*v1; d11=v1'*v1;  % using baricentric coordinates
    id=1/(d00*d11-d01*d01);
    for i=1:3
    v2=X(:,i)-Y(:,1); d02=v0'*v2; d12=v1'*v2;
    u=(d11*d02-d01*d12)*id; v=(d00*d12-d01*d02)*id;
    if u>=0 & v>=0 & u+v<=1            % also include nodes on the boundary
    k=k+1; P(:,k)=X(:,i);
    end;
    end;
    */
   int extraPoint = 0;
   for (int i = 0; i < 3; i++) {
      // compute twice the area of all 3 triangles formed by a side of Y and a corner of X; if one is negative, stop
      double A[2];
      for (int k = 0; k < 2; k++)
         A[k] = X[2 * i + k];
      int inside = 1;
      for (int j = 0; j < 3; j++) {
         double B[2], C[2];
         for (int k = 0; k < 2; k++) {
            B[k] = Y[2 * j + k];
            int j1 = (j + 1) % 3;
            C[k] = Y[2 * j1 + k];
         }

         double area2 = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1]
               - A[1]);
         if (area2 < 0.) {
            inside = 0;
            break;
         }
      }
      if (inside) {
         P[extraPoint * 2] = A[0];
         P[extraPoint * 2 + 1] = A[1];
         extraPoint++;
      }
   }
   return extraPoint;
}

int swap2(double * p, double * q) {
   double tmp = *p;
   *p = *q;
   *q = tmp;
   return 0;
}
int SortAndRemoveDoubles2(double * P, int & nP) {
   if (nP < 2)
      return 0; // nothing to do
   // center of gravity for the points
   double c[2] = { 0., 0. };
   int k = 0;
   for (k = 0; k < nP; k++) {
      c[0] += P[2 * k];
      c[1] += P[2 * k + 1];
   }
   c[0] /= nP;
   c[1] /= nP;
   double angle[12]; // could be at most 12 points; much less usually
   for (k = 0; k < nP; k++) {
      double x = P[2 * k] - c[0], y = P[2 * k + 1] - c[1];
      if (x != 0. || y != 0.)
         angle[k] = atan2(y, x);
      else {
         angle[k] = 0;
         // this means that the points are on a line, or all coincident // degenerate case
      }
   }
   // sort according to angle; also eliminate close points
   int sorted = 1;
   do {
      sorted = 1;
      for (k = 0; k < nP - 1; k++) {
         if (angle[k] > angle[k + 1]) {
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
   while (j < nP) {
      double d2 = dist2(&P[2 * i], &P[2 * j]);
      if (d2 > epsilon_1) {
         i++;
         P[2 * i] = P[2 * j];
         P[2 * i + 1] = P[2 * j + 1];
      }
      j++;
   }
   // test also the last point with the first one (index 0)

   double d2 = dist2(P, &P[2 * i]); // check the first and last points (ordered from -pi to +pi)
   if (d2 > epsilon_1) {
      nP = i + 1;
   } else
      nP = i; // effectively delete the last point (that would have been the same with first)
   if (nP == 0)
      nP = 1; // we should be left with at least one point we already tested if nP is 0 originally
   return 0;
}

// some local methods, utilities
// get the xy coordinates of a triangle
int getXYcoords(Interface * MB, EntityHandle tr, double * coords) {
   const EntityHandle * conn;
   int num_nodes;
   ErrorCode rval = MB->get_connectivity(tr, conn, num_nodes);
   if (MB_SUCCESS != rval || num_nodes != 3)
      return 1;
   double lCoords[9];
   rval = MB->get_coords(conn, 3, lCoords);
   if (MB_SUCCESS != rval)
      return 1;
   coords[0] = lCoords[0];
   coords[1] = lCoords[1];
   coords[2] = lCoords[3];
   coords[3] = lCoords[4];
   coords[4] = lCoords[6];
   coords[5] = lCoords[7];
   return 0;
}

int EdgeIntersections2(double * red, double * blue, int mark[3],
      double * points, int & nPoints) {
   /* EDGEINTERSECTIONS computes edge intersections of two triangles
    [P,n]=EdgeIntersections(X,Y) computes for the two given triangles  * red
    and blue ( stored column wise )
    (point coordinates are stored column-wise, in counter clock
    order) the points P where their edges intersect. In addition,
    in n the indices of which neighbors of red  are also intersecting
    with blue are given.
    */

   // points is an array with 12 slots   (12 * 2 doubles)
   nPoints = 0;
   mark[0] = mark[1] = mark[2] = 0; // no neighbors of red involved yet
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
   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
         double b[2];
         double a[2][2]; // 2*2
         int iPlus1 = (i + 1) % 3;
         int jPlus1 = (j + 1) % 3;
         for (int k = 0; k < 2; k++) {
            b[k] = blue[2 * j + k] - red[2 * i + k];
            // row k of a: a(k, 0), a(k, 1)
            a[k][0] = red[2 * iPlus1 + k] - red[2 * i + k];
            a[k][1] = blue[2 * j + k] - blue[2 * jPlus1 + k];

         }
         double delta = a[0][0] * a[1][1] - a[0][1] * a[1][0];
         if (delta != 0.) {
            // not parallel
            double alfa = (b[0] * a[1][1] - a[0][1] * b[1]) / delta;
            double beta = (-b[0] * a[1][0] + b[1] * a[0][0]) / delta;
            if (0 <= alfa && alfa <= 1. && 0 <= beta && beta <= 1.) {
               // the intersection is good
               for (int k = 0; k < 2; k++) {
                  points[2 * nPoints + k] = red[2 * i + k] + alfa * (red[2
                        * iPlus1 + k] - red[2 * i + k]);
               }
               mark[i] = 1; // so neighbor number i will be considered too.
               nPoints++;
            }
         }

      }
   }
   return 0;
}

int computeIntersectionBetweenRedAndBlue(EntityHandle red,
      EntityHandle blue, double * P, int & nP, double & area, int mark[3]) {
   // the points will be at most 9; they will describe a convex patch, after the points will be ordered and
   // collapsed (eliminate doubles)
   // the area is not really required

   double redTriangle[6];// column wise
   double blueTriangle[6];
   getXYcoords(mb, blue, blueTriangle);
   getXYcoords(mb, red, redTriangle);

   //we do not really need the mortar matrix

   //int n[3]={0, 0, 0};// no intersection of side red with blue
   //double area= 0.;
   // X corresponds to blue, Y to red
   nP = 0; // number of intersection points
   int ret = EdgeIntersections2(blueTriangle, redTriangle, mark, P, nP);
   if (ret != 0)
      return 1;// some unforeseen error

   int extraPoints = borderPointsOfXinY2(blueTriangle, redTriangle,
         &(P[2 * nP]));
   if (extraPoints > 1) {
      mark[0] = mark[1] = mark[2] = 1;
   }
   nP += extraPoints;
   extraPoints = borderPointsOfXinY2(redTriangle, blueTriangle, &(P[2 * nP]));
   nP += extraPoints;

   // now sort and orient the points in P, such that they are forming a convex polygon
   // this will be the foundation of our new mesh
   //
   SortAndRemoveDoubles2(P, nP); // nP should be at most 6 in the end
   // if there are more than 3 points, some area will be positive
   area = 0.;
   if (nP >= 3) {
      for (int k = 1; k < nP - 1; k++)
         area += area2D(P, &P[2 * k], &P[2 * k + 2]);
   }
   return 0; // no error
}

// specify also desired set; we are interested only in neighbors in the set!
ErrorCode GetOrderedNeighbors(Interface * MB, EntityHandle set,
      EntityHandle tri, EntityHandle neighbors[3]) {
   // will get the 3 ordered neighbors;
   // first triangle is for nodes 0, 1, second to 1, 2, third to 2, 0
   int nnodes;
   const EntityHandle * conn3;
   ErrorCode rval = MB->get_connectivity(tri, conn3, nnodes);
   if (MB_SUCCESS != rval || nnodes != 3)
      return MB_FAILURE;
   for (int i = 0; i < 3; i++) {
      EntityHandle v[2];
      v[0] = conn3[i];
      v[1] = conn3[(i + 1) % 3];
      // get triangles adjacent to vertices
      std::vector<EntityHandle> tris;
      rval = MB->get_adjacencies(v, 2, 2, false, tris, Interface::INTERSECT);
      if (MB_SUCCESS != rval)
         return rval;
      int siz = tris.size();
      if (siz > 2)
         return MB_FAILURE; // non-manifold
      if (siz == 1) {
         // it must be the border,
         neighbors[i] = 0; // we are guaranteed that ids are !=0
         continue;
      }
      // here siz ==2, it is either the first or second
      if (tri == tris[0])
         neighbors[i] = tris[1];
      else
         neighbors[i] = tris[0];
      // check if it is in the set; otherwise make it 0, like border
      if (!MB->contains_entities(set, &(neighbors[i]), 1)) {
         neighbors[i] = (EntityHandle) 0;
      }

   }
   return MB_SUCCESS;
}

int getXYcoordsAndNodeHandles(Interface * MB, EntityHandle tr,
      double * coords, const EntityHandle *& conn) {
   //const EntityHandle * conn;
   int num_nodes;
   ErrorCode rval = MB->get_connectivity(tr, conn, num_nodes);
   if (MB_SUCCESS != rval || num_nodes != 3)
      return 1;
   double lCoords[9];
   rval = MB->get_coords(conn, 3, lCoords);
   if (MB_SUCCESS != rval)
      return 1;
   coords[0] = lCoords[0];
   coords[1] = lCoords[1];
   coords[2] = lCoords[3];
   coords[3] = lCoords[4];
   coords[4] = lCoords[6];
   coords[5] = lCoords[7];
   return 0;
}

void CreateANewNode(double * coord, EntityHandle & outNode) {
   // this new node will have coordinate z as 0;
   // maybe we should offer more options;
   double coords[3] = { coord[0], coord[1], 0. };
   mb->create_vertex(coords, outNode);
}

// this method will also construct the triangles in the new mesh
// if we accept planar polygons, we just save them
// also, we could just create new vertices every time, and merge only in the end;
// could be too expensive, and the tolerance for merging could be an
// interesting topic
int findNodes(EntityHandle red, EntityHandle blue, double * iP, int nP) {
   // first of all, check against red and blue vertices
   //
   if (dbg_1) {
      std::cout << "red, blue, nP, P " << mb->id_from_handle(red) << " "
            << mb->id_from_handle(blue) << " " << nP << "\n";
      for (int n = 0; n < nP; n++)
         std::cout << " \t" << iP[2 * n] << "\t" << iP[2 * n + 1] << "\n";

   }
   // get coords and nodes for each of those 2 triangles
   double blueCoords[6];
   const EntityHandle * blueConn;
   getXYcoordsAndNodeHandles(mb, blue, blueCoords, blueConn);
   double redCoords[6];
   const EntityHandle * redConn;
   getXYcoordsAndNodeHandles(mb, red, redCoords, redConn);

   // get the edges for the red triangle; the extra points will be on those edges, saved as
   // lists (unordered)
   EntityHandle redEdges[3];
   int i = 0;
   for (i = 0; i < 3; i++) {
      EntityHandle v[2] = { redConn[i], redConn[(i + 1) % 3] };
      std::vector<EntityHandle> adj_entities;
      ErrorCode rval = mb->get_adjacencies(v, 2, 1, false, adj_entities,
            Interface::INTERSECT);
      if (rval != MB_SUCCESS || adj_entities.size() < 1)
         return 0; // get out , big error
      redEdges[i] = adj_entities[0]; // should be only one edge between 2 nodes
   }
   // these will be in the new mesh, mbOut

   EntityHandle * foundIds = new EntityHandle[nP];
   for (i = 0; i < nP; i++) {
      double * pp = &iP[2 * i];// iP+2*i
      int found = 0;
      // first, are they on vertices from red or blue?
      // priority is the red mesh (mb2?)
      int j = 0;
      EntityHandle outNode = (EntityHandle) 0;
      for (j = 0; j < 3 && !found; j++) {
         //int node = redTri.v[j];
         double d2 = dist2(pp, &redCoords[2 * j]);
         if (d2 < epsilon_1) {
            // suspect is redConn[j] corresponding in mbOut

            mb->tag_get_data(RedNodeTag, &redConn[j], 1, &outNode);
            if (outNode == (EntityHandle) 0) {
               // need to create a new node in mbOut
               CreateANewNode(&redCoords[2 * j], outNode);
               mb->tag_set_data(RedNodeTag, &redConn[j], 1, &outNode);
            }
            foundIds[i] = outNode; // no new node
            found = 1;
            if (dbg_1)
               std::cout << "  red node " << j << " " << mb->id_from_handle(
                     redConn[j]) << " new node " << mb->id_from_handle(
                     outNode) << " " << redCoords[2 * j] << "  " << redCoords[2
                     * j + 1] << " d2: " << d2 << " \n";
         }
      }
      // PSTriangle2D & blueTri = m_blueMesh[blue];
      for (j = 0; j < 3 && !found; j++) {
         //int node = blueTri.v[j];
         double d2 = dist2(pp, &blueCoords[2 * j]);
         if (d2 < epsilon_1) {
            // suspect is blueConn[j] corresponding in mbOut

            mb->tag_get_data(BlueNodeTag, &blueConn[j], 1, &outNode);
            if (outNode == (EntityHandle) 0) {
               // need to create a new node in mbOut
               CreateANewNode(&blueCoords[2 * j], outNode);
               mb->tag_set_data(BlueNodeTag, &blueConn[j], 1, &outNode);
            }
            foundIds[i] = outNode; // no new node
            found = 1;
            if (dbg_1)
               std::cout << "  blue node " << j << " " << mb->id_from_handle(
                     blueConn[j]) << " new node " << mb->id_from_handle(
                     outNode) << " " << blueCoords[2 * j] << blueCoords[2 * j
                     + 1] << " d2:" << d2 << " \n";
         }

      }
      if (!found) {
         // find the edge it belongs, first, on the red triangle
         //
         for (j = 0; j < 3; j++) {
            int j1 = (j + 1) % 3;
            double area = area2D(&redCoords[2 * j], &redCoords[2 * j1], pp);
            if (dbg_1)
               std::cout << "   edge " << j << ": " << mb->id_from_handle(
                     redEdges[j]) << " " << redConn[j] << " " << redConn[j1]
                     << "  area : " << area << "\n";
            if (fabs(area) < epsilon_1 * epsilon_1) {
               // found the edge; now find if there is a point in the list here
               //std::vector<EntityHandle> * expts = extraNodesMap[redEdges[j]];
               int indx =-1;
               indx = RedEdges.index(redEdges[j]);
               std::vector<EntityHandle> * expts = extraNodesVec[indx];
               // if the points pp is between extra points, then just give that id
               // if not, create a new point, (check the id)
               // get the coordinates of the extra points so far
               int nbExtraNodesSoFar = expts->size();
               double * coords1 = new double[3 * nbExtraNodesSoFar];
               mb->get_coords(&(*expts)[0], nbExtraNodesSoFar, coords1);
               //std::list<int>::iterator it;
               for (int k = 0; k < nbExtraNodesSoFar && !found; k++) {
                  //int pnt = *it;
                  double d2 = dist2(pp, &coords1[3 * k]);
                  if (d2 < epsilon_1) {
                     found = 1;
                     foundIds[i] = (*expts)[k];
                     if (dbg_1)
                        std::cout << " found node:" << foundIds[i] << std::endl;
                  }
               }
               if (!found) {
                  // create a new point in 2d (at the intersection)
                  //foundIds[i] = m_num2dPoints;
                  //expts.push_back(m_num2dPoints);
                  // need to create a new node in mbOut
                  // this will be on the edge, and it will be added to the local list
                  CreateANewNode(pp, outNode);
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
      if (!found) {
         std::cout << " a point pp is not on a red triangle " << *pp << " "
               << pp[1] << " red triangle " << mb->id_from_handle(red)
               << " \n";
         return 1;
      }
   }
   // now we can build the triangles, from P array, with foundIds
   // we will put them in the out set
   if (nP >= 3) {

       EntityHandle polyNew;
       mb->create_element(MBPOLYGON, foundIds, nP, polyNew);
       mb->add_entities(outSet, &polyNew, 1);
       if (dbg_1) {
          std::cout << " polygon " << mb->id_from_handle(polyNew)
                << "  nodes: " ;
          for (int i=0; i<nP; i++)
            std::cout << " " << mb->id_from_handle(foundIds[i]);
            std::cout    << "\n";
       }
       // tag it with the ids from red and blue
       int id = mb->id_from_handle(blue);
       mb->tag_set_data(blueParentTag, &polyNew, 1, &id);
       id = mb->id_from_handle(red);
       mb->tag_set_data(redParentTag, &polyNew, 1, &id);

   }
   delete[] foundIds;
   foundIds = NULL;
   return 0;
}

// clean some memory allocated
void clean() {
   //
   int indx=0;
   for (Range::iterator eit = RedEdges.begin(); eit != RedEdges.end(); eit++, indx++) {
      //EntityHandle edge = *eit;
      //delete extraNodesMap[edge];
      delete extraNodesVec[indx];
   }
   //extraNodesMap.clear();
   extraNodesVec.clear();
}

ErrorCode intersect_meshes(EntityHandle sf1, EntityHandle sf2, EntityHandle & outputSet)
{

  ErrorCode rval = MB_SUCCESS;
  rval = mb->create_meshset(MESHSET_SET, outputSet);
  if (MB_SUCCESS!=rval)
    return rval;

  mbs1 = (EntityHandle) sf1;
  mbs2 = (EntityHandle) sf2;
  outSet = outputSet;

  // really, should be something from t1 and t2; blue is 1, red is 2
  createTags(mb); //
  EntityHandle startBlue, startRed;// first triangles from mb1 and mb2
   //ErrorCode rval = mb1->handle_from_id(MBTRI, 1, startBlue);
   // we need to start somewhere; we will do an expensive search for one intersection
   //mb2->handle_from_id(MBTRI, 1, startRed);
   // this could be an expensive search
   // maybe we should do some KDtrees, for the worst case
   Range rs1;
   Range rs2;
   mb->get_entities_by_type(mbs1, MBTRI, rs1);
   mb->get_entities_by_type(mbs2, MBTRI, rs2);
   for (Range::iterator it = rs1.begin(); it != rs1.end(); it++) {
      startBlue = *it;
      int found = 0;
      for (Range::iterator it2 = rs2.begin(); it2 != rs2.end() && !found; it2++) {
         startRed = *it2;
         double area = 0;
         // if area is > 0 , we have intersections
         double P[24]; // max 6 points, but it may grow bigger; why worry
         int nP = 0;
         int n[3];// sides
         computeIntersectionBetweenRedAndBlue(startRed, startBlue, P, nP, area,
               n);
         if (area > 0) {
            found = 1;
            break; // found 2 triangles that intersect; these will be the seeds
         }
      }
      if (found)
         break;
   }

   std::queue<EntityHandle> blueQueue; // these are corresponding to Ta,
   blueQueue.push(startBlue);
   std::queue<EntityHandle> redQueue;
   redQueue.push(startRed);

   Range toResetReds;// will be used to reset red flags for every blue triangle
   // processed
   int k;

   if (dbg_1)
      mout_1.open("patches.m");
   unsigned char used = 1;
   unsigned char unused = 0; // for red flags
   while (!blueQueue.empty()) {
      // flags for the side : 0 means a red triangle not found on side
      // a paired red not found yet for the neighbors of blue
      EntityHandle n[3] = { EntityHandle(0), EntityHandle(0),
            EntityHandle(0) };

      EntityHandle currentBlue = blueQueue.front();
      blueQueue.pop();
      //        for (k=0; k<m_numPos; k++)
      //          redFlag[k] = 0;
      //        redFlag[m_numPos] = 1; // to guard for the boundary
      // all reds that were tagged, are now cleared
      for (Range::iterator itr = toResetReds.begin(); itr
            != toResetReds.end(); itr++) {
         EntityHandle ttt = *itr;
         rval = mb->tag_set_data(RedFlagTag, &ttt, 1, &unused);
      }
      //rval = mb2->tag_set_data(RedFlagTag, toResetReds, &unused);
      if (dbg_1) {
         std::cout << "reset reds: ";
         for (Range::iterator itr = toResetReds.begin(); itr
               != toResetReds.end(); itr++)
            std::cout << mb->id_from_handle(*itr) << " ";
         std::cout << std::endl;
      }
      EntityHandle currentRed = redQueue.front(); // where do we check for redQueue????
      // red and blue queues are parallel
      redQueue.pop();// mark the current red
      //redFlag[currentRed] = 1; //
      toResetReds.clear();// empty the range of used reds, will have to be set unused again,
      // at the end of blue triangle processing
      toResetReds.insert(currentRed);
      rval = mb->tag_set_data(RedFlagTag, &currentRed, 1, &used);
      //mb2->set_tag_data
      std::queue<EntityHandle> localRed;
      localRed.push(currentRed);
      while (!localRed.empty()) {
         //
         EntityHandle redT = localRed.front();
         localRed.pop();
         double P[24], area;
         int nP = 0; // intersection points
         int nc[3] = { 0, 0, 0 }; // means no intersection on the side (markers)
         computeIntersectionBetweenRedAndBlue(/* red */redT, currentBlue, P,
               nP, area, nc);
         if (nP > 0) {
            // intersection found: output P and original triangles if nP > 2
            if (dbg_1) {
               std::cout << "area: " << area << " nP:" << nP << std::endl;
               mout_1 << "pa=[\n";

               for (k = 0; k < nP; k++) {

                  mout_1 << P[2 * k] << "\t ";
               }

               mout_1 << "\n";
               for (k = 0; k < nP; k++) {

                  mout_1 << P[2 * k + 1] << "\t ";
               }

               mout_1 << " ]; \n";
               mout_1 << " patch(pa(1,:),pa(2,:),'m');       \n";
            }
            EntityHandle neighbors[3];
            rval = GetOrderedNeighbors(mb, mbs2, redT, neighbors);

            if (dbg_1) {
               std::cout << " neighbors for redT ";
               for (int kk = 0; kk < 3; kk++) {
                  if (neighbors[kk] > 0)
                     std::cout << mb->id_from_handle(neighbors[kk]) << " ";
                  else
                     std::cout << 0 << " ";
               }
               std::cout << std::endl;
            }
            // add neighbors to the localRed queue, if they are not marked
            for (int nn = 0; nn < 3; nn++) {
               EntityHandle neighbor = neighbors[nn];
               if (neighbor > 0) {
                  unsigned char status = 0;
                  mb->tag_get_data(RedFlagTag, &neighbor, 1, &status);
                  if (status == 0) {
                     localRed.push(neighbor);
                     rval = mb->tag_set_data(RedFlagTag, &neighbor, 1, &used);
                     //redFlag[neighbor] = 1; // flag it to not be added anymore
                     toResetReds.insert(neighbor); // this is used to reset the red flag
                  }
               }
               // n(find(nc>0))=ac;        % ac is starting candidate for neighbor
               if (nc[nn] > 0) // intersected
                  n[nn] = redT;// start from 0!!
            }
            if (nP > 1) // this will also construct triangles in the new mesh, if needed
               findNodes(redT, currentBlue, P, nP);
         }
         else
         {
           if (dbg_1)
             std::cout << " redT and currentBlue do not intersect " << mb->id_from_handle(redT)
                 <<" "<< mb->id_from_handle(currentBlue) << "\n";

         }

      }

      EntityHandle blueNeighbors[3];
      rval = GetOrderedNeighbors(mb, mbs1, currentBlue, blueNeighbors);
      if (dbg_1) {
         std::cout << " neighbors for blue T ";
         for (int kk = 0; kk < 3; kk++) {
            if (blueNeighbors[kk] > 0)
               std::cout << mb->id_from_handle(blueNeighbors[kk]) << " ";
            else
               std::cout << 0 << " ";
         }
         std::cout << std::endl;
      }
      for (int j = 0; j < 3; j++) {
         EntityHandle blueNeigh = blueNeighbors[j];
         unsigned char status = 1;
         if (blueNeigh == 0)
            continue;
         mb->tag_get_data(BlueFlagTag, &blueNeigh, 1, &status);// status 0 is unused
         if (status == 0 && n[j] > 0) // not treated yet and marked as a neighbor
         {
            // we identified triangle n[j] as intersecting with neighbor j of the blue triangle
            blueQueue.push(blueNeigh);
            redQueue.push(n[j]);
            if (dbg_1)
               std::cout << "new triangles pushed: blue, red:"
                     << mb->id_from_handle(blueNeigh) << " "
                     << mb->id_from_handle(n[j]) << std::endl;
            mb->tag_set_data(BlueFlagTag, &blueNeigh, 1, &used);
         }
      }

   }

   if (dbg_1)
      mout_1.close();
   //
   clean();
   return MB_SUCCESS;
}


int main(int argc, char* argv[])
{
  // check command line arg
  const char *filename_mesh1 = STRINGIFY(SRCDIR) "/m1.vtk";
  const char *filename_mesh2 = STRINGIFY(SRCDIR) "/m2.vtk";
  const char *newFile = "intx.vtk";
  if (argc == 4)
  {
    filename_mesh1 = argv[1];
    filename_mesh2 = argv[2];
    newFile = argv[3];
  }
  else
  {
    printf("Usage: %s <mesh_filename1> <mesh_filename2>  <newFile>\n", argv[0]);
    if (argc != 1)
      return 1;
    printf("No files specified.  Defaulting to: %s  %s  %s\n",
        filename_mesh1, filename_mesh2, newFile);
  }

  // read meshes in 2 file sets
  ErrorCode rval = MB_SUCCESS;
  Core moab;
  /*Interface**/ mb = &moab;// global
  EntityHandle sf1, sf2;
  rval = mb->create_meshset(MESHSET_SET, sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval = mb->create_meshset(MESHSET_SET, sf2);
  if (MB_SUCCESS != rval)
    return 1;
  rval=mb->load_file(filename_mesh1, &sf1);
  if (MB_SUCCESS != rval)
    return 1;
  rval=mb->load_file(filename_mesh2, &sf2);
  if (MB_SUCCESS != rval)
   return 1;

  EntityHandle outputSet;
  rval = intersect_meshes(sf1, sf2, outputSet);

  rval = mb->write_mesh(newFile, &outputSet, 1);
   return 0;


}
