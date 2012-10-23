/*
 * IntxUtils.cpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */

#include "IntxUtils.hpp"
#include <math.h>

namespace moab {
// vec utilities that could be common between quads on a plane or sphere
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
int SortAndRemoveDoubles2(double * P, int & nP, double epsilon_1)
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


// vec utils related to gnomonic projection on a sphere

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
void decide_gnomonic_plane(const CartVect & pos, int & plane)
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
int gnomonic_projection(const CartVect & pos, double R, int plane, double & c1, double & c2)
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
int reverse_gnomonic_projection(const double & c1, const double & c2, double R, int plane,
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

} //namespace moab
