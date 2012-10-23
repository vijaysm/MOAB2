/*
 * IntxUtils.hpp
 *
 *  Created on: Oct 3, 2012
 */

#ifndef INTXUTILS_HPP_
#define INTXUTILS_HPP_

#include "moab/CartVect.hpp"
namespace moab
{
double dist2(double * a, double * b);
double area2D(double *a, double *b, double *c);
int borderPointsOfXinY2(double * X, double * Y, double * P, int side[4]);
int SortAndRemoveDoubles2(double * P, int & nP, double epsilon);
// the marks will show what edges of blue intersect the red

int EdgeIntersections2(double * blue, double * red, int markb[4], int markr[4],
    double * points, int & nPoints);

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
void decide_gnomonic_plane(const CartVect & pos, int & oPlane);
// point on a sphere is projected on one of six planes, decided earlier
int gnomonic_projection(const CartVect & pos, double R, int plane, double & c1, double & c2);
// given the position on plane (one out of 6), find out the position on sphere
int reverse_gnomonic_projection(const double & c1, const double & c2, double R, int plane,
    CartVect & pos);
}
#endif /* INTXUTILS_HPP_ */
