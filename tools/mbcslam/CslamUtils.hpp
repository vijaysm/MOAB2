/*
 * CslamUtils.hpp
 *
 *  Created on: Oct 3, 2012
 */

#ifndef CSLAMUTILS_HPP_
#define CSLAMUTILS_HPP_

#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"

namespace moab
{
double dist2(double * a, double * b);
double area2D(double *a, double *b, double *c);
int borderPointsOfXinY2(double * X, double * Y, int nsides, double * P, int side[4]);
int SortAndRemoveDoubles2(double * P, int & nP, double epsilon);
// the marks will show what edges of blue intersect the red

int EdgeIntersections2(double * blue, double * red, int nsides, int markb[4], int markr[4],
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
 * this is similar logic to /cesm1_0_4/models/atm/cam/src/dynamics/homme/share/coordinate_systems_mod.F90
 *    method: function cart2face(cart3D) result(face_no)
 */
void decide_gnomonic_plane(const CartVect & pos, int & oPlane);
// point on a sphere is projected on one of six planes, decided earlier
int gnomonic_projection(const CartVect & pos, double R, int plane, double & c1, double & c2);
// given the position on plane (one out of 6), find out the position on sphere
int reverse_gnomonic_projection(const double & c1, const double & c2, double R, int plane,
    CartVect & pos);

/*
 *   other methods to convert from spherical coord to cartesian, and back
 *   A spherical coordinate triple is (R, lon, lat)
 *   should we store it as a CartVect? probably not ...
 *   /cesm1_0_4/models/atm/cam/src/dynamics/homme/share/coordinate_systems_mod.F90
 *
     enforce three facts:
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
 */

struct SphereCoords{
  double R, lon, lat;
};

SphereCoords cart_to_spherical(CartVect &) ;

CartVect spherical_to_cart (SphereCoords &) ;

/*
 *  create a mesh used mainly for visualization for now, with nodes corresponding to
 *   GL points, a so-called refined mesh, with NP nodes in each direction.
 *   input: a range of quads (coarse), and  a desired order (NP is the number of points), so it
 *   is order + 1
 *
 *   output: a set with refined elements; with proper input, it should be pretty
 *   similar to a Homme mesh read with ReadNC
 */
ErrorCode SpectralVisuMesh(Interface * mb, Range & input, int NP, EntityHandle & outputSet);

/*
 * given an entity set, get all nodes and project them on a sphere with given radius
 */
ErrorCode ProjectOnSphere(Interface * mb, EntityHandle set, double R);

bool point_in_interior_of_convex_polygon (double * points, int np, double pt[2]);

/*
 * utilities to compute area of a polygon on which all edges are arcs of great circles on a sphere
 */
/*
 * this will compute the spherical angle ABC, when A, B, C are on a sphere of radius R
 *  the radius will not be needed, usually, just for verification the points are indeed on that sphere
 *  the center of the sphere is at origin (0,0,0)
 *  this angle can be used in Girard's theorem to compute the area of a spherical polygon
 */
double spherical_angle(double * A, double * B, double * C, double Radius);

// this could be larger than PI, because of orientation; useful for non-convex polygons
double oriented_spherical_angle(double * A, double * B, double * C);

double area_spherical_triangle(double *A, double *B, double *C, double Radius);

double area_spherical_polygon (double * A, int N, double Radius);

double area_on_sphere(Interface * mb, EntityHandle set, double R);
}
#endif /* CSLAMUTILS_HPP_ */
