/*
 * CslamUtils.cpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */

#include "CslamUtils.hpp"
#include <math.h>
// this is from mbcoupler; maybe it should be moved somewhere in moab src
// right now, add a dependency to mbcoupler
#include "ElemUtil.hpp"
#include "moab/MergeMesh.hpp"

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
int borderPointsOfXinY2(double * X, double * Y, int nsides, double * P, int side[4])
{
  // 2 triangles, 3 corners, is the corner of X in Y?
  // Y must have a positive area
  /*
   */
  int extraPoint = 0;
  for (int i = 0; i < nsides; i++)
  {
    // compute twice the area of all 4 triangles formed by a side of Y and a corner of X; if one is negative, stop
    double * A = X + 2 * i;

    int inside = 1;
    for (int j = 0; j < nsides; j++)
    {
      double * B = Y + 2 * j;

      int j1 = (j + 1) % nsides;
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

int EdgeIntersections2(double * blue, double * red, int nsides, int markb[4], int markr[4],
    double * points, int & nPoints)
{
  /* EDGEINTERSECTIONS computes edge intersections of two elements
   [P,n]=EdgeIntersections(X,Y) computes for the two given elements  * red
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
  for (int i = 0; i < nsides; i++)
  {
    for (int j = 0; j < nsides; j++)
    {
      double b[2];
      double a[2][2]; // 2*2
      int iPlus1 = (i + 1) % nsides;
      int jPlus1 = (j + 1) % nsides;
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

/*
 *
    use physical_constants, only : dd_pi
    type(cartesian3D_t), intent(in) :: cart
    type(spherical_polar_t)         :: sphere

    sphere%r=distance(cart)
    sphere%lat=ASIN(cart%z/sphere%r)
    sphere%lon=0

    ! ==========================================================
    ! enforce three facts:
    !
    ! 1) lon at poles is defined to be zero
    !
    ! 2) Grid points must be separated by about .01 Meter (on earth)
    !    from pole to be considered "not the pole".
    !
    ! 3) range of lon is { 0<= lon < 2*pi }
    !
    ! ==========================================================

    if (distance(cart) >= DIST_THRESHOLD) then
       sphere%lon=ATAN2(cart%y,cart%x)
       if (sphere%lon<0) then
          sphere%lon=sphere%lon+2*DD_PI
       end if
    end if

  end function cart_to_spherical
 */
SphereCoords cart_to_spherical(CartVect & cart3d)
{
  SphereCoords res;
  res.R = cart3d.length();
  if (res.R < 0)
  {
    res.lon=res.lat = 0.;
    return res;
  }
  res.lat = asin(cart3d[2]/res.R);
  res.lon=atan2(cart3d[1], cart3d[0]);
  if (res.lon<0)
    res.lon+=2*M_PI;// M_PI is defined in math.h? it seems to be true, although
  // there are some defines it depends on :(
  // #if defined __USE_BSD || defined __USE_XOPEN ???

  return res;
}
/*
 * ! ===================================================================
  ! spherical_to_cart:
  ! converts spherical polar {lon,lat}  to 3D cartesian {x,y,z}
  ! on unit sphere
  ! ===================================================================

  function spherical_to_cart(sphere) result (cart)

    type(spherical_polar_t), intent(in) :: sphere
    type(cartesian3D_t)                 :: cart

    cart%x=sphere%r*COS(sphere%lat)*COS(sphere%lon)
    cart%y=sphere%r*COS(sphere%lat)*SIN(sphere%lon)
    cart%z=sphere%r*SIN(sphere%lat)

  end function spherical_to_cart
 */
CartVect spherical_to_cart (SphereCoords & sc)
{
  CartVect res;
  res[0] = sc.R * cos(sc.lat)*cos(sc.lon) ; // x coordinate
  res[1] = sc.R * cos(sc.lat)*sin(sc.lon); // y
  res[2] = sc.R * sin(sc.lat);             // z
  return res;
}

ErrorCode SpectralVisuMesh(Interface * mb, Range & input, int NP, EntityHandle & outputSet)
{
  ErrorCode rval = MB_SUCCESS;

  // this will compute the gl points on parametric element
  moab::Element::SpectralQuad specQuad(NP);
  Range fineElem;

  // get all edges? or not? Form all gl points + quads, then merge, then output
  for (Range::iterator it=input.begin(); it!=input.end(); it++)
  {
    const moab::EntityHandle * conn4 = NULL;
    int num_nodes = 0;
    rval = mb->get_connectivity(*it, conn4, num_nodes);
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get conectivity for quad \n";
      return rval;
    }
    assert(num_nodes==4);

    std::vector<CartVect> verts(4);
    rval = mb->get_coords(conn4, num_nodes, &(verts[0][0]));
    if (moab::MB_SUCCESS != rval)
    {
      std::cout << "can't get coords for quad \n";
      return rval;
    }

    specQuad.set_vertices(verts);
    specQuad.compute_gl_positions();
    double * xyz[3];
    int size;
    specQuad.get_gl_points(xyz[0], xyz[1], xyz[2], size);
    assert(NP*NP==size);
    // these points are in lexicographic order;
    std::vector<EntityHandle> nodes(size);
    for (int i=0; i<size; i++)
    {
      double coords[3]={ xyz[0][i], xyz[1][i], xyz[2][i] };
      rval = mb->create_vertex(coords, nodes[i] );
      if (rval!=moab::MB_SUCCESS)
        return rval;
    }
    // create (NP-1)^2 quads, one by one (sic)
    for (int i=0; i<NP-1; i++)
    {
      for (int j=0; j<NP-1; j++)
      {
        EntityHandle connec[4]={ nodes[ i*NP+j],      nodes[i*NP+j+1],
                                  nodes[(i+1)*NP+j+1], nodes[(i+1)*NP+j] };
        EntityHandle element_handle;
        rval = mb->create_element(MBQUAD, connec, 4, element_handle);
        if (rval!=moab::MB_SUCCESS)
          return rval;
        fineElem.insert(element_handle);
      }
    }

  }

  mb->add_entities(outputSet, fineElem);
  // merge all nodes
  MergeMesh mgtool(mb);
  rval = mgtool.merge_entities(fineElem, 0.0001);

  return rval;
}

ErrorCode ProjectOnSphere(Interface * mb, EntityHandle set, double R)
{
  Range ents;
  ErrorCode rval = mb->get_entities_by_handle(set, ents);
  if (rval!=moab::MB_SUCCESS)
    return rval;

  Range nodes;
  rval = mb->get_connectivity( ents, nodes);
  if (rval!=moab::MB_SUCCESS)
    return rval;

  // one by one, get the node and project it on the sphere, with a radius given
  // the center of the sphere is at 0,0,0
  for (Range::iterator nit=nodes.begin(); nit!=nodes.end(); nit++)
  {
    EntityHandle nd=*nit;
    CartVect pos;
    rval = mb->get_coords(&nd, 1, (double*)&(pos[0]));
    if (rval!=moab::MB_SUCCESS)
      return rval;
    double len=pos.length();
    if (len==0.)
      return MB_FAILURE;
    pos = R/len*pos;
    rval = mb->set_coords(&nd, 1, (double*)&(pos[0]));
    if (rval!=moab::MB_SUCCESS)
       return rval;
  }
  return MB_SUCCESS;
}
//
bool point_in_interior_of_convex_polygon (double * points, int np, double pt[2])
{
  bool inside = true;
  // assume points are counterclockwise
  for (int i = 0; i < np; i++)
  {
    // compute the area of triangles formed by a side of polygon and the pt; if one is negative, stop
    double * A = points + 2 * i;

    int i1=(i+1)%np;
    double * B = points + 2 * i1; // no copy of data

    double area2 = (B[0] - A[0]) * (pt[1] - A[1])
        - (pt[0] - A[0]) * (B[1] - A[1]);
    if (area2 < 0.)
    {
      inside = false;
      break;
    }
  }
  return inside;
}
// assume they are one the same sphere
double spherical_angle(double * A, double * B, double * C, double Radius)
{
  // the angle by definition is between the planes OAB and OBC
  CartVect a(A);
  CartVect b(B);
  CartVect c(C);
  double err1=a.length_squared()-Radius*Radius;
  if (fabs(err1)>0.0001)
  {
    std::cout << " error in input " << a << " radius: " << Radius << " error:" << err1<< "\n";
  }
  CartVect normalOAB = a * b;
  CartVect normalOCB = c * b;
  return angle(normalOAB, normalOCB);
}
// could be bigger than M_PI;
// angle at B could be bigger than M_PI, if the orientation is such that ABC points toward the interior
double oriented_spherical_angle(double * A, double * B, double * C)
{
  // assume the same radius, sphere at origin
  CartVect a(A), b(B), c(C);
  CartVect normalOAB = a * b;
  CartVect normalOCB = c * b;
  CartVect orient = (b-a)*(c-a);
  double ang = angle(normalOAB, normalOCB); // this is between 0 and M_PI
  if (ang!=ang)
  {
    // signal of a nan
    std::cout << a << " " << b << " " << c <<"\n";
    std::cout << ang << "\n";
  }
  if (orient%a < 0)
    return (2*M_PI-ang);// the other angle

  return ang;

}
double area_spherical_triangle(double *A, double *B, double *C, double Radius)
{
  double correction = spherical_angle(A, B, C, Radius)+spherical_angle(B, C, A, Radius)+
      spherical_angle(C, A, B, Radius)-M_PI;
  double area = Radius*Radius*correction;
  // now, is it negative or positive? is it pointing toward the center or outward?
  CartVect a(A), b(B), c(C);
  CartVect abc = (b-a)*(c-a);
  if (abc%a > 0) // dot product positive, means ABC points out
    return area;
  else
    return -area;

}

double area_spherical_polygon (double * A, int N, double Radius)
{
  // this should work for non-convex polygons too
  // assume that the A, A+3, ..., A+3*(N-1) are the coordinates
  //
  if (N<=2)
    return 0.;
  double sum_angles = 0.;
  for (int i=0; i<N; i++)
  {
    int i1 = (i+1)%N;
    int i2 = (i+2)%N;
    sum_angles += oriented_spherical_angle( A+3*i, A+3*i1, A+3*i2);
  }
  double correction = sum_angles-(N-2)*M_PI;
  return Radius*Radius *correction;

}
/*
 *  l'Huiller's formula for spherical triangle
 *  http://williams.best.vwh.net/avform.htm
 *  a, b, c are arc measures in radians, too
 *  A, B, C are angles on the sphere, for which we already have formula
 *               c
 *         A -------B
 *          \       |
 *           \      |
 *            \b    |a
 *             \    |
 *              \   |
 *               \  |
 *                \C|
 *                 \|
 *
 *  (The angle at B is not necessarily a right angle)
 *
 *    sin(a)  sin(b)   sin(c)
 *    ----- = ------ = ------
 *    sin(A)  sin(B)   sin(C)
 *
 * In terms of the sides (this is excess, as before, but numerically stable)
 *
 *  E = 4*atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)))
 */

//! Interior angle between two vectors
inline double angleSafe( const CartVect& u, const CartVect& v )
  {
    double a = ( (u % v) / std::sqrt((u % u) * (v % v)) );
    if (a> 1.) a= 1.;
    if (a<-1.) a =-1.;
    return std::acos(a);
  }

double area_spherical_triangle_lHuiller(double * ptA, double * ptB, double * ptC, double Radius)
{
  // now, a is angle BOC, O is origin
  CartVect vA(ptA), vB(ptB), vC(ptC);
  double a = angleSafe(vB, vC);
  double b = angleSafe(vC, vA);
  double c = angleSafe(vA, vB);
  int sign =1;
  if ((vA*vB)%vC < 0)
    sign = -1;
  double s = (a+b+c)/2;
  double tmp = tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2);
  if (tmp<0.)
    tmp = 0.;
  double E = 4*atan(sqrt(tmp));
  if (E!=E)
    std::cout << "bad nan here \n";
  return sign * E * Radius*Radius;
}


double area_spherical_polygon_lHuiller (double * A, int N, double Radius)
{
  // this should work for non-convex polygons too
  // assume that the A, A+3, ..., A+3*(N-1) are the coordinates
  //
  // assume that the orientation is positive;
  // if negative orientation, the area will be negative
  if (N<=2)
    return 0.;
  double area = 0.;
  for (int i=1; i<N-1; i++)
  {
    int i1 = i+1;
    area += area_spherical_triangle_lHuiller( A, A+3*i, A+3*i1, Radius);
  }
  return area;
}

double area_on_sphere(Interface * mb, EntityHandle set, double R)
{
  // get all entities of dimension 2
  // then get the connectivity, etc
  Range inputRange;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, inputRange);
  if (MB_SUCCESS != rval)
    return -1;

  // compare total area with 4*M_PI * R^2
  double total_area = 0.;
  for (Range::iterator eit=inputRange.begin(); eit!= inputRange.end(); eit++ )
  {
    EntityHandle eh=*eit;
    // get the nodes, then the coordinates
    const EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(eh, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return -1;
    std::vector<double> coords(3*num_nodes);
    // get coordinates
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return -1;
    total_area+=area_spherical_polygon (&coords[0], num_nodes, R);
  }
  return total_area;
}
double area_on_sphere_lHuiller(Interface * mb, EntityHandle set, double R)
{
  // get all entities of dimension 2
  // then get the connectivity, etc
  Range inputRange;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, inputRange);
  if (MB_SUCCESS != rval)
    return -1;

  double total_area = 0.;
  for (Range::iterator eit = inputRange.begin(); eit != inputRange.end(); eit++)
  {
    EntityHandle eh = *eit;
    // get the nodes, then the coordinates
    const EntityHandle * verts;
    int num_nodes;
    rval = mb->get_connectivity(eh, verts, num_nodes);
    if (MB_SUCCESS != rval)
      return -1;
    std::vector<double> coords(3 * num_nodes);
    // get coordinates
    rval = mb->get_coords(verts, num_nodes, &coords[0]);
    if (MB_SUCCESS != rval)
      return -1;
    total_area += area_spherical_polygon_lHuiller(&coords[0], num_nodes, R);
  }
  return total_area;
}
/*
 *
 */
double distance_on_great_circle(CartVect & p1, CartVect & p2)
{
  SphereCoords sph1 = cart_to_spherical(p1);
  SphereCoords sph2 = cart_to_spherical(p2);
  // radius should be the same
  return sph1.R * acos(sin (sph1.lon)* sin(sph2.lon) + cos(sph1.lat)*cos(sph2.lat)*cos(sph2.lon-sph2.lon) );
}
/*
 * based on paper A class of deformational flow test cases for linear transport problems on the sphere
 *  longitude lambda [0.. 2*pi) and latitude theta (-pi/2, pi/2)
 *  lambda: -> lon (0, 2*pi)
 *  theta : -> lat (-pi/2, po/2)
 *  Radius of the sphere is 1 (if not, everything gets multiplied by 1)
 *
 *  cosine bell: center lambda0, theta0:
 */
void departure_point_case1(CartVect & arrival_point, double t, double delta_t, CartVect & departure_point)
{
  // always assume radius is 1 here?
  SphereCoords sph = cart_to_spherical(arrival_point);
  double T=5; // duration of integration (5 units)
  double k = 2.4; //flow parameter
  /*     radius needs to be within some range   */
  double  sl2 = sin(sph.lon/2);
  double pit = M_PI * t / T;
  double omega = M_PI/T;
  double costetha = cos(sph.lat);
  //double u = k * sl2*sl2 * sin(2*sph.lat) * cos(pit);
  double v = k * sin(sph.lon) * costetha * cos(pit);
  //double psi = k * sl2 * sl2 *costetha * costetha * cos(pit);
  double u_tilda = 2*k*sl2*sl2*sin(sph.lat)*cos(pit);

  // formula 35, page 8
  // this will approximate dep point using a Taylor series with up to second derivative
  // this will be O(delta_t^3) exact.
  double lon_dep = sph.lon - delta_t * u_tilda -delta_t*delta_t * k * sl2 *
      ( sl2 * sin (sph.lat) * sin(pit) * omega
          - u_tilda * sin(sph.lat) * cos(pit) * cos (sph.lon/2)
          - v * sl2 * costetha * cos(pit)   );
  // formula 36, page 8 again
  double lat_dep = sph.lat - delta_t*v - delta_t * delta_t/4* k *
      ( sin(sph.lon)* cos(sph.lat) * sin(pit) * omega
          - u_tilda * cos(sph.lon) * cos(sph.lat) * cos(pit)
          + v * sin(sph.lon) * sin(sph.lat) * cos(pit)  );
  SphereCoords sph_dep;
  sph_dep.R = 1.; // radius
  sph_dep.lat = lat_dep;
  sph_dep.lon = lon_dep;

  departure_point = spherical_to_cart(sph_dep);
  return;
}
} //namespace moab
