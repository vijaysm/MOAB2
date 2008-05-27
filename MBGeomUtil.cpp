/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**\file MBGeometry.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-27
 */

#include "MBCartVect.hpp"
#include "MBCN.hpp"
#include "MBGeomUtil.hpp"
#include "MBElemUtil.hpp"
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <iostream>

#ifdef _MSC_VER
#  include <float.h>
#  define finite(A) _finite(A)
#endif

namespace MBGeomUtil {

bool segment_box_intersect( MBCartVect box_min,
                            MBCartVect box_max,
                            const MBCartVect& seg_pt,
                            const MBCartVect& seg_unit_dir,
                            double& seg_start, double& seg_end )
{
    // translate so that seg_pt is at origin
  box_min -= seg_pt;
  box_max -= seg_pt;
  
  for (unsigned i = 0; i < 3; ++i) {  // X, Y, and Z slabs

      // intersect line with slab planes
    const double t_min = box_min[i] / seg_unit_dir[i];
    const double t_max = box_max[i] / seg_unit_dir[i];
    
      // check if line is parallel to planes
    if (!finite(t_min)) {
      if (box_min[i] > 0.0 || box_max[i] < 0.0)
        return false; 
      continue;
    }

    if (seg_unit_dir[i] < 0) {
      if (t_min < seg_end) 
        seg_end = t_min;
      if (t_max > seg_start)
        seg_start = t_max;
    }
    else { // seg_unit_dir[i] > 0
      if (t_min > seg_start)
        seg_start = t_min; 
      if (t_max < seg_end)
        seg_end = t_max;
    }
  }

  return seg_start <= seg_end;
}


/* Impelementation copied from cgmMC ray_tri_contact (overlap.C) */
bool ray_tri_intersect( const MBCartVect vertices[3],
                        const MBCartVect& b,
                        const MBCartVect& v,
                        double /*tolerance*/,
                        double& t_out,
                        const double* ray_length)
{
  const MBCartVect p0 = vertices[0] - vertices[1]; // abc
  const MBCartVect p1 = vertices[0] - vertices[2]; // def
                                                   // ghi<-v
  const MBCartVect p = vertices[0] - b;            // jkl
  const MBCartVect c = p1 * v;                     // eiMinushf,gfMinusdi,dhMinuseg
  const double mP = p0 % c;
  const double betaP = p % c;
  if (mP > 0) {
    if (betaP < 0)
      return false;
  }
  else if (mP < 0) {
    if (betaP > 0)
      return false;
  }
  else {
    return false;
  }
  
  const MBCartVect d = p0 * p; // jcMinusal,blMinuskc,akMinusjb
  double gammaP = v % d;
  if (mP > 0) {
    if (gammaP < 0 || betaP + gammaP > mP)
      return false;
  }
  else if (betaP + gammaP < mP || gammaP > 0)
    return false;
  
  const double tP = p1 % d;
  const double m = 1.0 / mP;
  const double beta = betaP * m;
  const double gamma = gammaP * m;
  const double t = -tP * m;
  if (ray_length && t > *ray_length)
    return false;
  
  if (beta < 0 || gamma < 0 ||
      beta + gamma > 1 ||
      t < 0.0)
    return false;
  
  t_out = t;
  return true;
}

bool box_plane_overlap( const MBCartVect& normal,
                        double d,
                        MBCartVect min,
                        MBCartVect max )
{
  if (normal[0] < 0.0)
    std::swap( min[0], max[0] );
  if (normal[1] < 0.0)
    std::swap( min[1], max[1] );
  if (normal[2] < 0.0)
    std::swap( min[2], max[2] );
  
  return (normal % min <= -d) && (normal % max >= -d);
}


#define CHECK_RANGE( A, B, R ) \
  if ((A) < (B)) { \
    if ((A) > (R) || (B) < -(R)) \
      return false; \
  } \
  else if ((B) > (R) || (A) < -(R)) \
    return false

/* Adapted from: http://jgt.akpeters.com/papers/AkenineMoller01/tribox.html
 * Use separating axis theorem to test for overlap between triangle
 * and axis-aligned box.
 *
 * Test for overlap in these directions:
 * 1) {x,y,z}-directions 
 * 2) normal of triangle
 * 3) crossprod of triangle edge with {x,y,z}-direction
 */
bool box_tri_overlap( const MBCartVect vertices[3],
                      const MBCartVect& box_center,
                      const MBCartVect& box_dims )
{
    // translate everthing such that box is centered at origin
  const MBCartVect v0( vertices[0] - box_center );
  const MBCartVect v1( vertices[1] - box_center );
  const MBCartVect v2( vertices[2] - box_center );

  // do case 1) tests
  if (v0[0] > box_dims[0] && v1[0] > box_dims[0] && v2[0] > box_dims[0])
    return false;
  if (v0[1] > box_dims[1] && v1[1] > box_dims[1] && v2[1] > box_dims[1])
    return false;
  if (v0[2] > box_dims[2] && v1[2] > box_dims[2] && v2[2] > box_dims[2])
    return false;
  if (v0[0] < -box_dims[0] && v1[0] < -box_dims[0] && v2[0] < -box_dims[0])
    return false;
  if (v0[1] < -box_dims[1] && v1[1] < -box_dims[1] && v2[1] < -box_dims[1])
    return false;
  if (v0[2] < -box_dims[2] && v1[2] < -box_dims[2] && v2[2] < -box_dims[2])
    return false;
  
    // compute triangle edge vectors
  const MBCartVect e0( vertices[1] - vertices[0] );
  const MBCartVect e1( vertices[2] - vertices[1] );
  const MBCartVect e2( vertices[0] - vertices[2] );
  
    // do case 3) tests 
  double fex, fey, fez, p0, p1, p2, rad;
  fex = fabs(e0[0]);
  fey = fabs(e0[1]);
  fez = fabs(e0[2]);
  
  p0 = e0[2]*v0[1] - e0[1]*v0[2];
  p2 = e0[2]*v2[1] - e0[1]*v2[2];
  rad = fez * box_dims[1] + fey * box_dims[2];
  CHECK_RANGE( p0, p2, rad );
  
  p0 = -e0[2]*v0[0] + e0[0]*v0[2];
  p2 = -e0[2]*v2[0] + e0[0]*v2[2];
  rad = fez * box_dims[0] + fex * box_dims[2];
  CHECK_RANGE( p0, p2, rad );
    
  p1 = e0[1]*v1[0] - e0[0]*v1[1];
  p2 = e0[1]*v2[0] - e0[0]*v2[1];
  rad = fey * box_dims[0] + fex * box_dims[1];
  CHECK_RANGE( p1, p2, rad );
  
  fex = fabs(e1[0]);
  fey = fabs(e1[1]);
  fez = fabs(e1[2]);
  
  p0 = e1[2]*v0[1] - e1[1]*v0[2];
  p2 = e1[2]*v2[1] - e1[1]*v2[2];
  rad = fez * box_dims[1] + fey * box_dims[2];
  CHECK_RANGE( p0, p2, rad );
  
  p0 = -e1[2]*v0[0] + e1[0]*v0[2];
  p2 = -e1[2]*v2[0] + e1[0]*v2[2];
  rad = fez * box_dims[0] + fex * box_dims[2];
  CHECK_RANGE( p0, p2, rad );
  
  p0 = e1[1]*v0[0] - e1[0]*v0[1];
  p1 = e1[1]*v1[0] - e1[0]*v1[1];
  rad = fey * box_dims[0] + fex * box_dims[1];
  CHECK_RANGE( p0, p1, rad );
  
  fex = fabs(e2[0]);
  fey = fabs(e2[1]);
  fez = fabs(e2[2]);
  
  p0 = e2[2]*v0[1] - e2[1]*v0[2];
  p1 = e2[2]*v1[1] - e2[1]*v1[2];
  rad = fez * box_dims[1] + fey * box_dims[2];
  CHECK_RANGE( p0, p1, rad );
  
  p0 = -e2[2]*v0[0] + e2[0]*v0[2];
  p1 = -e2[2]*v1[0] + e2[0]*v1[2];
  rad = fez * box_dims[0] + fex * box_dims[2];
  CHECK_RANGE( p0, p1, rad );
  
  p1 = e2[1]*v1[0] - e2[0]*v1[1];
  p2 = e2[1]*v2[0] - e2[0]*v2[1];
  rad = fey * box_dims[0] + fex * box_dims[1];
  CHECK_RANGE( p1, p2, rad );
  
  // do case 2) test
  MBCartVect n = e0 * e1;
  return box_plane_overlap( n, -(n % v0), -box_dims, box_dims );
}
  

bool box_tri_overlap( const MBCartVect  triangle_corners[3],
                      const MBCartVect& box_min_corner,
                      const MBCartVect& box_max_corner,
                      double            tolerance )
{
  const MBCartVect box_center = 0.5 * (box_max_corner + box_min_corner);
  const MBCartVect box_hf_dim = 0.5 * (box_max_corner - box_min_corner);
  return box_tri_overlap( triangle_corners,
                          box_center,
                          box_hf_dim + MBCartVect(tolerance) );
} 

bool box_elem_overlap( const MBCartVect *elem_corners,
                       MBEntityType elem_type,
                       const MBCartVect& center,
                       const MBCartVect& dims )
{

  switch (elem_type) {
    case MBTRI:
      return box_tri_overlap( elem_corners, center, dims );
    case MBHEX:
      return box_hex_overlap( elem_corners, center, dims );
    case MBPOLYGON:
    case MBPOLYHEDRON:
      assert(false);
      return false;
    default:
      return box_linear_elem_overlap( elem_corners, elem_type, center, dims );
  }
}

static inline MBCartVect quad_norm( const MBCartVect& v1,
                                    const MBCartVect& v2,
                                    const MBCartVect& v3,
                                    const MBCartVect& v4 )
{ return (-v1+v2+v3-v4) * (-v1-v2+v3+v4); }

static inline MBCartVect tri_norm( const MBCartVect& v1,
                                   const MBCartVect& v2,
                                   const MBCartVect& v3 )
{ return (v2-v1) * (v3-v1); }


bool box_linear_elem_overlap( const MBCartVect *elem_corners,
                              MBEntityType type,
                              const MBCartVect& box_center,
                              const MBCartVect& box_halfdims )
{
    // Do Separating Axis Theorem:
    // If the element and the box overlap, then the 1D projections
    // onto at least one of the axes in the following three sets
    // must overlap (assuming convex polyhedral element).
    // 1) The normals of the faces of the box (the principal axes)
    // 2) The crossproduct of each element edge with each box edge
    //    (crossproduct of each edge with each principal axis)
    // 3) The normals of the faces of the element

  unsigned d, i, e, f;             // loop counters
  bool all_less, all_greater;      // track overlap (or lack thereof)
  double min, max, cross[2], tmp;
  MBCartVect norm;
  int indices[4]; // element edge/face vertex indices
    // calculate minimum and maximum corners of extents box
  const MBCartVect box_min( box_center - box_halfdims );
  const MBCartVect box_max( box_center + box_halfdims );
    // get element topology information
  const unsigned num_corner = MBCN::VerticesPerEntity( type );
  const unsigned num_edge = MBCN::NumSubEntities( type, 1 );
  const unsigned num_face = MBCN::NumSubEntities( type, 2 );
  
    // test box face normals (principal axes)
  for (d = 0; d < 3; ++d) {  // for each principal axis
    all_less = all_greater = true;
    for (i = 0; i < num_corner; ++i) { // for each element corner
      if (elem_corners[i][d] > box_min[d])
        all_less = false;
      if (elem_corners[i][d] < box_max[d])
        all_greater = false;
    }
    if (all_greater || all_less)
      return false;
  }
  
    // test edge-edge crossproducts
  for (d = 0; d < 3; ++d) {  // for each principal axis (box edge)
      // get indices of other two axes
    const int idx1 = (d+1)%3;
    const int idx2 = (d+2)%3;
    
    for (e = 0; e < num_edge; ++e) { // for each element edge
        // get which element vertices bound the edge
      MBCN::SubEntityVertexIndices( type, 1, e, indices );
        // calculate crossproduct: axis x (v1 - v0),
        // where v1 and v0 are edge vertices.
      cross[0] = elem_corners[indices[0]][idx2] - elem_corners[indices[1]][idx2];
      cross[1] = elem_corners[indices[1]][idx1] - elem_corners[indices[0]][idx1];
        // skip if orthogonal
      if ((cross[0]*cross[0] + cross[1]*cross[1]) < std::numeric_limits<double>::epsilon())
        continue;
      
        // first box vertex
      min = max = cross[0] * box_min[idx1] + cross[1] * box_min[idx2];
        // second box vertex
      tmp = cross[0] * box_min[idx1] + cross[1] * box_max[idx2];
      if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
        // third box vertex
      tmp = cross[0] * box_max[idx1] + cross[1] * box_max[idx2];
      if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
        // fourth box vertex
      tmp = cross[0] * box_max[idx1] + cross[1] * box_min[idx2];
      if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
    
      all_less = all_greater = true;
      for (i = (unsigned)(indices[0]+1)%num_corner; i != (unsigned)indices[0]; i = (i+1)%num_corner) { // for each element corner
        tmp = cross[0] * elem_corners[i][idx1] + cross[1] * elem_corners[i][idx2];
        if (tmp > min)
          all_less = false;
        if (tmp < max)
          all_greater = false;
      }
      
      if (all_less || all_greater)
        return false;
    }
  }
  
    // test element face normals
  for (f = 0; f < num_face; ++f) {
    MBCN::SubEntityVertexIndices( type, 2, f, indices );
    switch (MBCN::SubEntityType( type, 2, f )) {
      case MBTRI:
        norm = tri_norm( elem_corners[indices[0]], 
                         elem_corners[indices[1]], 
                         elem_corners[indices[2]] );
        break;
      case MBQUAD:
        norm = quad_norm( elem_corners[indices[0]], 
                          elem_corners[indices[1]], 
                          elem_corners[indices[2]], 
                          elem_corners[indices[3]] );
        break;
      default:
        assert(false);
        continue;
    }

      // first box corner
    min = max = norm % box_min;
      // second corner
    tmp = norm % MBCartVect( box_max[0], box_min[1], box_min[2] );
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
      // third corner
    tmp = norm % MBCartVect( box_max[0], box_max[1], box_min[2] );
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
      // fourth corner
    tmp = norm % MBCartVect( box_min[0], box_max[1], box_min[2] );
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
      // fifth corner
    tmp = norm % MBCartVect( box_min[0], box_min[1], box_max[2] );
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
      // sixth corner
    tmp = norm % MBCartVect( box_max[0], box_min[1], box_max[2] );
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
      // seventh corner
    tmp = norm % box_max;
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
      // eighth corner
    tmp = norm % MBCartVect( box_min[0], box_max[1], box_max[2] );
    if (tmp < min) min = tmp; else if (tmp > max) max = tmp;
    
    // for each element vertex
    all_less = all_greater = true;
    for (i = 0; i < num_corner; ++i) { 
      tmp = norm % elem_corners[i];
      if (tmp > min)
        all_less = false;
      if (tmp < max)
        all_greater = false;
    }

    if (all_less || all_greater)
      return false;
  }
  
    // Overlap on all tested axes.
  return true;
}
        
 
  

//from: http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf#search=%22closest%20point%20on%20triangle%22
/*       t
 *   \(2)^
 *    \  |
 *     \ |
 *      \|
 *       \
 *       |\
 *       | \
 *       |  \  (1)
 *  (3)  tv  \
 *       |    \
 *       | (0) \
 *       |      \
 *-------+---sv--\----> s
 *       |        \ (6)
 *  (4)  |   (5)   \
 */
// Worst case is either 61 flops and 5 compares or 53 flops and 6 compares,
// depending on relative costs.  For all paths that do not return one of the
// corner vertices, exactly one of the flops is a divide.
void closest_location_on_tri( const MBCartVect& location,
                              const MBCartVect* vertices,
                              MBCartVect& closest_out )
{                                                     // ops      comparisons
  const MBCartVect sv( vertices[1] - vertices[0] );   // +3 = 3
  const MBCartVect tv( vertices[2] - vertices[0] );   // +3 = 6
  const MBCartVect pv( vertices[0] - location );      // +3 = 9
  const double ss = sv % sv;                          // +5 = 14
  const double st = sv % tv;                          // +5 = 19
  const double tt = tv % tv;                          // +5 = 24
  const double sp = sv % pv;                          // +5 = 29
  const double tp = tv % pv;                          // +5 = 34
  const double det = ss*tt - st*st;                   // +3 = 37
  double s = st*tp - tt*sp;                           // +3 = 40
  double t = st*sp - ss*tp;                           // +3 = 43
  if (s+t < det) {                                    // +1 = 44, +1 = 1
    if (s < 0) {                                      //          +1 = 2
      if (t < 0) {                                    //          +1 = 3
        // region 4
        if (sp < 0) {                                 //          +1 = 4
          if (-sp > ss)                               //          +1 = 5
            closest_out = vertices[1];                //      44       5
          else
            closest_out = vertices[0] - (sp/ss) * sv; // +7 = 51,      5
        }
        else if (tp < 0) {                            //          +1 = 5
          if (-tp > tt)                               //          +1 = 6
            closest_out = vertices[2];                //      44,      6
          else
            closest_out = vertices[0] - (tp/tt) * tv; // +7 = 51,      6
        }
        else {
          closest_out = vertices[0];                  //      44,      5
        }
      }
      else {
        // region 3
        if (tp >= 0)                                  //          +1 = 4
          closest_out = vertices[0];                  //      44,      4
        else if (-tp >= tt)                           //          +1 = 5
          closest_out = vertices[2];                  //      44,      5
        else
          closest_out = vertices[0] - (tp/tt) * tv;   // +7 = 51,      5
      }
    }
    else if (t < 0) {                                 //          +1 = 3
      // region 5;
      if (sp >= 0.0)                                  //          +1 = 4
        closest_out = vertices[0];                    //      44,      4
      else if (-sp >= ss)                             //          +1 = 5
        closest_out = vertices[1];                    //      44       5
      else
        closest_out = vertices[0] - (sp/ss) * sv;     // +7 = 51,      5
    }
    else {
      // region 0
      const double inv_det = 1.0 / det;               // +1 = 45
      s *= inv_det;                                   // +1 = 46
      t *= inv_det;                                   // +1 = 47
      closest_out = vertices[0] + s*sv + t*tv;        //+12 = 59,      3  
    }
  }
  else {
    if (s < 0) {                                      //          +1 = 2
      // region 2
      s = st + sp;                                    // +1 = 45
      t = tt + tp;                                    // +1 = 46
      if (t > s) {                                    //          +1 = 3
        const double num = t - s;                     // +1 = 47
        const double den = ss - 2*st + tt;            // +3 = 50
        if (num > den)                                //          +1 = 4
          closest_out = vertices[1];                  //      50,      4
        else {
          s = num/den;                                // +1 = 51
          t = 1 - s;                                  // +1 = 52
          closest_out = s*vertices[1] + t*vertices[2];// +9 = 61,      4
        }
      }
      else if (t <= 0)                                //          +1 = 4
        closest_out = vertices[2];                    //      46,      4
      else if (tp >= 0)                               //          +1 = 5
        closest_out = vertices[0];                    //      46,      5
      else
        closest_out = vertices[0] - (tp/tt) * tv;     // +7 = 53,      5
    }
    else if (t < 0) {                                 //          +1 = 3
      // region 6
      t = st + tp;                                    // +1 = 45
      s = ss + sp;                                    // +1 = 46
      if (s > t) {                                    //          +1 = 4
        const double num = t - s;                     // +1 = 47
        const double den = tt - 2*st + ss;            // +3 = 50
        if (num > den)                                //          +1 = 5
          closest_out = vertices[2];                  //      50,      5
        else {
          t = num/den;                                // +1 = 51
          s = 1 - t;                                  // +1 = 52
          closest_out = s*vertices[1] + t*vertices[2];// +9 = 61,      5
        }
      }
      else if (s <= 0)                                //          +1 = 5
        closest_out = vertices[1];                    //      46,      5
      else if (sp >= 0)                               //          +1 = 6
        closest_out = vertices[0];                    //      46,      6
      else
        closest_out = vertices[0] - (sp/ss) * sv;     // +7 = 53,      6
    }
    else {
      // region 1
      const double num = tt + tp - st - sp;           // +3 = 47
      if (num <= 0) {                                 //          +1 = 4
        closest_out = vertices[2];                    //      47,      4
      }
      else {
        const double den = ss - 2*st + tt;            // +3 = 50
        if (num >= den)                               //          +1 = 5
          closest_out = vertices[1];                  //      50,      5
        else {
          s = num/den;                                // +1 = 51
          t = 1 - s;                                  // +1 = 52
          closest_out = s*vertices[1] + t*vertices[2];// +9 = 61,      5
        }
      }
    }
  }
}

void closest_location_on_tri( const MBCartVect& location,
                              const MBCartVect* vertices,
                              double tolerance,
                              MBCartVect& closest_out,
                              int& closest_topo )
{
  const double tsqr = tolerance*tolerance;
  int i;
  MBCartVect pv[3], ev, ep;
  double t;

  closest_location_on_tri( location, vertices, closest_out );
  
  for (i = 0; i < 3; ++i) {
    pv[i] = vertices[i] - closest_out;
    if ((pv[i] % pv[i]) <= tsqr) {
      closest_topo = i;
      return;
    }
  }
  
  for (i = 0; i < 3; ++i) {
    ev = vertices[(i+1)%3] - vertices[i];
    t = (ev % pv[i]) / (ev % ev);
    ep = closest_out - (vertices[i] + t * ev);
    if ((ep % ep) <= tsqr) {
      closest_topo = i+3;
      return;
    }
  }
  
  closest_topo = 6;
}
 
    
// We assume polygon is *convex*, but *not* planar.
void closest_location_on_polygon( const MBCartVect& location,
                                  const MBCartVect* vertices,
                                  int num_vertices,
                                  MBCartVect& closest_out )
{
  const int n = num_vertices;
  MBCartVect d, p, v;
  double shortest_sqr, dist_sqr, t_closest, t;
  int i, e;
  
    // Find closest edge of polygon.
  e = n - 1;
  v = vertices[0] - vertices[e];
  t_closest = (v % (location - vertices[e])) / (v % v);
  if (t_closest < 0.0)
    d = location - vertices[e];
  else if (t_closest > 1.0)
    d = location - vertices[0];
  else 
    d = location - vertices[e] - t_closest * v;
  shortest_sqr = d % d;
  for (i = 0; i < n - 1; ++i) {
    v = vertices[i+1] - vertices[i];
    t = (v % (location - vertices[i])) / (v % v);
    if (t < 0.0)
      d = location - vertices[i];
    else if (t > 1.0)
      d = location - vertices[i+1];
    else
      d = location - vertices[i] - t * v;
    dist_sqr = d % d;
    if (dist_sqr < shortest_sqr) {
      e = i;
      shortest_sqr = dist_sqr;
      t_closest = t;
    }
  }
  
    // If we are beyond the bounds of the edge, then
    // the point is outside and closest to a vertex
  if (t_closest <= 0.0) {
    closest_out = vertices[e];
    return;
  }
  else if (t_closest >= 1.0) {
    closest_out = vertices[(e+1)%n];
    return;
  }
  
    // Now check which side of the edge we are one
  const MBCartVect v0 = vertices[e] - vertices[(e+n-1)%n];
  const MBCartVect v1 = vertices[(e+1)%n] - vertices[e];
  const MBCartVect v2 = vertices[(e+2)%n] - vertices[(e+1)%n];
  const MBCartVect norm = (1.0 - t_closest) * (v0 * v1) + t_closest * (v1 * v2);
    // if on outside of edge, result is closest point on edge
  if ((norm % ((vertices[e] - location) * v1)) <= 0.0) { 
    closest_out = vertices[e] + t_closest * v1;
    return;
  }
  
    // Inside.  Project to plane defined by point and normal at
    // closest edge
  const double D = -(norm % (vertices[e] + t_closest * v1));
  closest_out = (location - (norm % location + D) * norm)/(norm % norm);
}

void closest_location_on_box( const MBCartVect& min,
                              const MBCartVect& max,
                              const MBCartVect& point,
                              MBCartVect& closest )
{
  closest[0] = point[0] < min[0] ? min[0] : point[0] > max[0] ? max[0] : point[0];
  closest[1] = point[1] < min[1] ? min[1] : point[1] > max[1] ? max[1] : point[1];
  closest[2] = point[2] < min[2] ? min[2] : point[2] > max[2] ? max[2] : point[2];
}

bool box_point_overlap( const MBCartVect& box_min_corner,
                        const MBCartVect& box_max_corner,
                        const MBCartVect& point,
                        double tolerance )
{
  MBCartVect closest;
  closest_location_on_box( box_min_corner, box_max_corner, point, closest );
  closest -= point;
  return closest % closest < tolerance * tolerance;
}

bool box_hex_overlap( const MBCartVect hex_vertices[8],
                      const MBCartVect& box_center,
                      const MBCartVect& box_dims)
{

      // Mapping of vertices to each face on MBHEX
      const unsigned int facids[] = {0, 1, 5, 4,   1, 2, 6, 5,   2, 3, 7, 6,
                                     3, 0, 4, 7,   3, 2, 1, 0,   4, 5, 6, 7};
      // Mapping of vertices to each edge on MBHEX      
      const unsigned int edgids[] = {0, 1,   1, 2,   2, 3,   3, 0,  
                                     0, 4,   1, 5,   2, 6,   3, 7,
                                     4, 5,   5, 6,   6, 7,   7, 4}; 

      const double eps = 1.e-10;

      // Center the hex such that origin located at the box
      MBCartVect hexv[8];
      for (unsigned int i=0; i < 8; i++) {
            hexv[i] = hex_vertices[i] - box_center;
      }

      // Check hex against box normals...
      // Need to check both positive and negative faces
      double       t;
      unsigned int positive = 0, negative = 0;
      unsigned int out_plus = 0, out_minus = 0, inside = 0;
      bool         test;
      // Negative face, loop over all vertices in hex
      for (unsigned int i=0; i < 3; i++) {   
            // vxyz = box_center[i] + box_dims[i];
            out_plus = 0; out_minus = 0; inside = 0;
            // Test each vertex on hex...  
            for (unsigned int j=0; j < 8; j++) {

                  const double bxd = box_dims[i] + eps;

                  if (hexv[j][i] > bxd )  out_plus++;
                  else if (hexv[j][i] < -bxd )  out_minus++;
                  else inside++;

                  test = ((inside) || ( out_plus && out_minus ));
                  if (test) break;
            }
            if (!test) return false;
      }

      // Construct the vertices of the box...
      MBCartVect boxv[8];

      boxv[0][0] = -box_dims[0]; boxv[0][1] = -box_dims[1]; boxv[0][2] = -box_dims[2];
      boxv[1][0] =  box_dims[0]; boxv[1][1] = -box_dims[1]; boxv[1][2] = -box_dims[2];
      boxv[2][0] =  box_dims[0]; boxv[2][1] =  box_dims[1]; boxv[2][2] = -box_dims[2];
      boxv[3][0] = -box_dims[0]; boxv[3][1] =  box_dims[1]; boxv[3][2] = -box_dims[2];
      boxv[4][0] = -box_dims[0]; boxv[4][1] = -box_dims[1]; boxv[4][2] =  box_dims[2];
      boxv[5][0] =  box_dims[0]; boxv[5][1] = -box_dims[1]; boxv[5][2] =  box_dims[2];
      boxv[6][0] =  box_dims[0]; boxv[6][1] =  box_dims[1]; boxv[6][2] =  box_dims[2];
      boxv[7][0] = -box_dims[0]; boxv[7][1] =  box_dims[1]; boxv[7][2] =  box_dims[2];

      // Check box against hex normals ...
      // Loop over each face of hex
      for (unsigned int i = 0; i < 6; i++ ) {
            // Compute the normal
        const MBCartVect midpt1 = 0.5 * (( hexv[(facids[4*i+2])] + hexv[(facids[4*i+3])] )
                                    - ( hexv[(facids[4*i])]   + hexv[(facids[4*i+1])] ));
        const MBCartVect midpt2 = 0.5 * (( hexv[(facids[4*i+3])] + hexv[(facids[4*i])]   )
                                    - ( hexv[(facids[4*i+1])] + hexv[(facids[4*i+2])] ));
        const MBCartVect normal = midpt1 * midpt2;

        // Loop over each vertex in the box
        positive = 0; negative = 0;
        for (unsigned int j=0; j < 8; j++) { 
          // Take dot product of the vector and normal
          t = ( boxv[j] - hexv[(facids[4*i])] ) % normal;

          // Do the comparison
          if (t > eps )  positive++;
          else           negative++;
          if (positive && negative) break;
        }
        // std::cout << positive << "  " << negative << std::endl;
        test = ( (positive && !negative) ? false : true );
        if (!test) return false;
      }

      // Edge check ...
      // For box, only need to check three edges due to orthogonality
      int side1, side2;
      MBCartVect edge_cross;

      for (unsigned int d = 0; d < 3; d++ ) {
        // const unsigned int d0 = d % 3;
        const unsigned int d1 = (d + 1) % 3;
        const unsigned int d2 = (d + 2) % 3;

        for (unsigned int i = 0; i < 12; i++) {
          // edge_cross[d0] = 0.;
          edge_cross[d1] = hexv[(edgids[2*i])][d2] - hexv[(edgids[2*i+1])][d2];
          edge_cross[d2] = hexv[(edgids[2*i+1])][d1] - hexv[(edgids[2*i])][d1];

          if ((edge_cross[d1]*edge_cross[d1] + edge_cross[d2]*edge_cross[d2]) < eps) continue;

          side1 = side2 = 0;

          // Hex edge and box vertex tests
          positive = 0; negative = 0; test = false;
          for (unsigned int j = 0; j < 8; j++) {
            // t = edge_cross % ( boxv[j] - hexv[(edgids[2*i])] );
            t =  edge_cross[d1] * ( boxv[j][d1] - hexv[(edgids[2*i])][d1] )
               + edge_cross[d2] * ( boxv[j][d2] - hexv[(edgids[2*i])][d2] );
            if (t > eps) positive++;
            else if (t < -eps) negative++;

            if (positive && negative) {
              test = true;
              break;
            }
          }
          if (!test) side1 = ( (positive && !negative) ? +1 : -1 );
          if (side1 == 0) continue;   

          // Hex edge and hex vertex tests
          positive = 0; negative = 0; test = false;
          for (unsigned int j = 0; j < 8; j++) {
            // t = edge_cross % ( hexv[j] - hexv[(edgids[2*i])] );
            t =  edge_cross[d1] * ( boxv[j][d1] - hexv[(edgids[2*i])][d1] )
               + edge_cross[d2] * ( boxv[j][d2] - hexv[(edgids[2*i])][d2] );
            if (t > eps) positive++;
            else if (t < -eps) negative++;

            if (positive && negative) {
              test = true;
              break;
            }
          }
          if (!test) side2 = ( (positive && !negative) ? +1 : -1 );
          if (side2 == 0) continue;        

          // Test if on opposite sides
          if (side1 * side2 < 0 ) return false;

        }

      }

      // All possibilities exhausted, must overlap ...
      return true;
}


bool point_in_trilinear_hex(MBCartVect hex[8], 
                            MBCartVect xyz,
                            double etol) 
{

      const double one = 1.000001;

      MBCartVect  nat(0.);
      MBElemUtil::nat_coords_trilinear_hex(hex, xyz, nat, etol);
      
      for (unsigned int i = 0; i < 3; i++) {
            if ((nat[i] > one) || (nat[i] < -one)) return false;
      }

      return true;

}


bool point_in_trilinear_hex(MBCartVect hex[8], 
                            MBCartVect xyz, 
                            MBCartVect box_min, MBCartVect box_max,
                            double etol) 
{

      const double one = 1.000001;

      if ((xyz[0] < box_min[0]) || (xyz[0] > box_max[0])) return false;
      if ((xyz[1] < box_min[1]) || (xyz[1] > box_max[1])) return false;
      if ((xyz[2] < box_min[2]) || (xyz[2] > box_max[2])) return false;

      MBCartVect  nat(0.);
      MBElemUtil::nat_coords_trilinear_hex(hex, xyz, nat, etol);
      
      for (unsigned int i = 0; i < 3; i++) {
            if ((nat[i] > one) || (nat[i] < -one)) return false;
      }

      return true;

}

} // namespace MBGeoemtry
