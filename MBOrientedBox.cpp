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

/* 
 * The algorithms for the calculation of the oriented box from a
 * set of points or a set of cells was copied from the implemenation
 " in the "Visualization Toolkit".  J.K. - 2006-07-19
 *
 * Program:   Visualization Toolkit
 * Module:    $RCSfile$
 *
 * Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 * All rights reserved.
 * See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
 */

/**\file MBOrientedBox.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */
 
#include "MBOrientedBox.hpp"
#include "MBRange.hpp"
#include "MBMatrix3.hpp"
#include <ostream>
#include <assert.h>
 
std::ostream& operator<<( std::ostream& s, const MBOrientedBox& b )
{
  return s << b.center 
           << " + " 
           << b.axis[0] 
#if MB_ORIENTED_BOX_UNIT_VECTORS
           << ":" << b.length[0] 
#endif
           << " x " 
           << b.axis[1] 
#if MB_ORIENTED_BOX_UNIT_VECTORS
           << ":" << b.length[1] 
#endif
           << " x " 
           << b.axis[2]
#if MB_ORIENTED_BOX_UNIT_VECTORS
           << ":" << b.length[2] 
#endif
            ;
}

/**\brief Find closest point on line
 *
 * Find the point on the line for which a line trough the 
 * input point \a p and the result position is orthogonal to
 * the input line.
 * \param p  The point for which to find the perpendicular
 * \param b  A point on the line
 * \param m  The direction of the line
 * \return   The location on the line specified as 't' in the
 *           formula t * m + b
 */
static double point_perp( const MBCartVect& p,   // closest to this point
                          const MBCartVect& b,   // point on line
                          const MBCartVect& m )  // line direction
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  double t = (m % (p - b));
#else
  double t = (m % (p - b)) / (m % m);
#endif
  return finite(t) ? t : 0.0;
}

MBErrorCode MBOrientedBox::tag_handle( MBTag& handle_out,
                                       MBInterface* instance,
                                       const char* name,
                                       bool create )
{
    // We're going to assume this when mapping the MBOrientedBox
    // to tag data, so assert it.  
#if MB_ORIENTED_BOX_OUTER_RADIUS
  const size_t rad_size = sizeof(double);
#else
  const size_t rad_size = 0;
#endif
#if MB_ORIENTED_BOX_UNIT_VECTORS
  const size_t SIZE = rad_size + 15 * sizeof(double);
#else
  const size_t SIZE = rad_size + 12 * sizeof(double);
#endif
  assert( sizeof(MBOrientedBox) == SIZE );
  
  MBErrorCode rval = instance->tag_get_handle( name, handle_out );
  if (rval == MB_TAG_NOT_FOUND)
  {
    rval = instance->tag_create( name, 
                                 SIZE,
                                 MB_TAG_DENSE,
                                 MB_TYPE_DOUBLE,
                                 handle_out,
                                 0 );
  }
  else if (rval == MB_SUCCESS)
  {
    MBDataType type;
    rval = instance->tag_get_data_type( handle_out, type );
    if (MB_SUCCESS != rval)
      return rval;
    
    int size;
    rval = instance->tag_get_size( handle_out, size );
    if (MB_SUCCESS != rval)
      return rval;
    
    if (type != MB_TYPE_DOUBLE && type != MB_TYPE_OPAQUE)
      return MB_FAILURE;
    if ((unsigned)size != SIZE)
      return MB_FAILURE;
  }
  
  return rval;
}

/**\brief Common code for box calculation
 *
 * Given the orientation of the box and an approximate center,
 * calculate the exact center and extents of the box.
 * 
 *\param result.center  As input, the approximate center of the box.
 *                      As output, the exact center of the box.
 *\param result.axes    As input, directions of principal axes corresponding
 *                      to the orientation of the box.  As output, the 
 *                      direction of the principal axes of the box with a 
 *                      length equal to the distance from the center to the
 *                      corresponding side of the box.
 *\param points  The set of points the box should contain.
 */
static MBErrorCode box_from_axes( MBOrientedBox& result,
                                  MBInterface* instance,
                                  const MBRange& points )
{ 
  MBErrorCode rval;
  
    // project points onto axes to get box extents
  MBCartVect min(std::numeric_limits<double>::max()), 
             max(-std::numeric_limits<double>::max());
  for (MBRange::iterator i = points.begin(); i != points.end(); ++i)
  {
    MBCartVect coords;
    rval = instance->get_coords( &*i, 1, coords.array() );
    if (MB_SUCCESS != rval)
      return rval;
    
    for (int d = 0; d < 3; ++d)
    {
      double t = point_perp( coords, result.center, result.axis[d] );
      if (t < min[d])
        min[d] = t;
      if (t > max[d])
        max[d] = t;
    }
  }
  
    // We now have a box defined by three orthogonal line segments
    // that intersect at the center of the box.  Each line segment
    // is defined as result.center + t * result.axis[i], where the
    // range of t is [min[i], max[i]].
  
    // Calculate new center
  MBCartVect mid = 0.5 * (min + max);
  result.center += mid[0] * result.axis[0] +
                   mid[1] * result.axis[1] +
                   mid[2] * result.axis[2];
  
    // reorder axes by length
  MBCartVect range = 0.5 * (max - min);
  if (range[2] < range[1])
  {
    if (range[2] < range[0]) {
      std::swap( range[0], range[2] );
      std::swap( result.axis[0], result.axis[2] );
    }
  }
  else if (range[1] < range[0]) {
    std::swap( range[0], range[1] );
    std::swap( result.axis[0], result.axis[1] );
  }
  if (range[1] > range[2]) {
    std::swap( range[1], range[2] );
    std::swap( result.axis[1], result.axis[2] );
  }

    // scale axis to encompass all points, divide in half
#if MB_ORIENTED_BOX_UNIT_VECTORS
  result.length = range;
#else
  result.axis[0] *= range[0];
  result.axis[1] *= range[1];
  result.axis[2] *= range[2];
#endif

#if MB_ORIENTED_BOX_OUTER_RADIUS
  result.radius = range.length();
#endif

  return MB_SUCCESS;
}


MBErrorCode MBOrientedBox::compute_from_vertices( MBOrientedBox& result,
                                                  MBInterface* instance,
                                                  const MBRange& vertices )
{
  const MBRange::iterator begin = vertices.lower_bound( MBVERTEX );
  const MBRange::iterator end = vertices.upper_bound( MBVERTEX );
  size_t count = 0;
  
    // compute mean
  MBCartVect v;
  result.center = MBCartVect( 0, 0, 0 );
  for (MBRange::iterator i = begin; i != end; ++i)
  {
    MBErrorCode rval = instance->get_coords( &*i, 1, v.array() );
    if (MB_SUCCESS != rval)
      return rval;
    result.center += v;
    ++count;
  }
  result.center /= count;
  
    // compute covariance matrix
  MBMatrix3 a( 0.0 );
  for (MBRange::iterator i = begin; i != end; ++i)
  {
    MBErrorCode rval = instance->get_coords( &*i, 1, v.array() );
    if (MB_SUCCESS != rval)
      return rval;
  
    v -= result.center;
    a += outer_product( v, v );
  }
  a /= count;

    // Get axes (Eigenvectors) from covariance matrix
  double lambda[3];
  EigenDecomp( a, lambda, result.axis );
  
    // Calculate center and extents of box given orientation defined by axes
  return box_from_axes( result, instance, vertices );
}


MBErrorCode MBOrientedBox::compute_from_2d_cells( MBOrientedBox& result,
                                                  MBInterface* instance,
                                                  const MBRange& elements )
{
  MBErrorCode rval;
  const MBRange::iterator begin = elements.lower_bound( MBCN::TypeDimensionMap[2].first );
  const MBRange::iterator end = elements.lower_bound( MBCN::TypeDimensionMap[3].first );
  
    // compute mean and moments
  MBMatrix3 a(0.0);
  result.center = MBCartVect(0.0);
  double total_area = 0;
  for (MBRange::iterator i = begin; i != end; ++i)
  {
    const MBEntityHandle* conn;
    int conn_len;
    rval = instance->get_connectivity( *i, conn, conn_len );
    if (MB_SUCCESS != rval)
      return rval;
    
      // for each triangle in the 2-D cell
    for (int j = 2; j < conn_len; ++j)
    {
      MBEntityHandle vertices[3] = { conn[0], conn[j-1], conn[j] };
      MBCartVect coords[3];
      rval = instance->get_coords( vertices, 3, coords[0].array() );
      if (MB_SUCCESS != rval)
        return rval;
      
        // edge vectors
      const MBCartVect edge0 = coords[1] - coords[0];
      const MBCartVect edge1 = coords[2] - coords[0];
      const MBCartVect centroid = (coords[0] + coords[1] + coords[2]) / 3;
      const double tri_area2 = (edge0 * edge1).length();
      total_area += tri_area2;
      result.center += tri_area2 * centroid;
      
      a += tri_area2 * (9 * outer_product( centroid,  centroid  ) +
                            outer_product( coords[0], coords[0] ) +
                            outer_product( coords[1], coords[1] ) +
                            outer_product( coords[2], coords[2] ));
    } // for each triangle
  } // for each element
  result.center /= total_area;
  
    // get covariance matrix from moments
  a /= 12 * total_area;
  a -= outer_product( result.center, result.center );

    // get axes (Eigenvectors) from covariance matrix
  double lamda[3];
  EigenDecomp( a, lamda, result.axis );
  
    // Calculate center and extents of box given orientation defined by axes
  MBRange points;
  rval = instance->get_adjacencies( elements, 0, false, points, MBInterface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  return box_from_axes( result, instance, points );
}      

bool MBOrientedBox::contained( const MBCartVect& point, double tol ) const
{
  MBCartVect from_center = point - center;
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return fabs(from_center % axis[0]) - length[0] <= tol &&
         fabs(from_center % axis[1]) - length[1] <= tol &&
         fabs(from_center % axis[2]) - length[2] <= tol ;
#else
  for (int i = 0; i < 3; ++i) {
    double length = axis[i].length();
    if (fabs(from_center % axis[i]) - length*length > length*tol)
      return false;
  }
  return true;
#endif
}


//bool MBOrientedBox::contained( const MBOrientedBox& box, double tol ) const
//{
//  for (int i = -1; i < 2; i += 2) 
//  {
//    for (int j = -1; j < 2; j += 2) 
//    {
//      for (int k = -1; k < 2; k += 2) 
//      {
//        MBCartVect corner( center );
//#ifdef MB_ORIENTED_BOX_UNIT_VECTORS
//        corner += i * box.length[0] * box.axis[0];
//        corner += j * box.length[1] * box.axis[1];
//        corner += k * box.length[2] * box.axis[2];
//#else
//        corner += i * box.axis[0];
//        corner += j * box.axis[1];
//        corner += k * box.axis[2];
//#endif
//        if (!contained( corner, tol ))
//          return false;
//      }
//    }
//  }
//  return true;
//}


/* This implementation copied from cgmMC (overlap.C).
 * Original author:  Tim Taugtes?
 */
bool MBOrientedBox::intersect_ray( const MBCartVect& b,
                                   const MBCartVect& m,
                                   double reps,
                                   const double* len ) const
{
    // test distance from box center to line
  const MBCartVect cx = center - b;
  double dist_s = cx % m;
  double dist_sq = cx % cx - (dist_s*dist_s);
  double max_diagsq = outer_radius_squared();
  
    // if greater than the longest diagonal, we don't hit
  if (dist_sq > max_diagsq+reps)
    return false;
  if (len && dist_s - max_diagsq > *len)
    return false;
  
    // if smaller than shortest diagonal, we do hit
  if (dist_sq < inner_radius_squared() - reps && dist_s >= 0.0)
    return true;
    
    // get transpose of axes
    // Note: if axes were stored as a matrix, could skip
    // transpose and just switch order of operands in
    // matrix-vector multiplies below. - J.K.
  //MBMatrix3 B( axis[0][0], axis[1][0], axis[2][0],
  //             axis[0][1], axis[1][1], axis[2][1],
  //             axis[0][2], axis[1][2], axis[2][2] );
  MBMatrix3 B( axis[0][0], axis[0][1], axis[0][2],
               axis[1][0], axis[1][1], axis[1][2],
               axis[2][0], axis[2][1], axis[2][2] );
  //MBCartVect T = B * -center;
  
    // transform ray to box coordintae system
  //MBCartVect par_pos = T + B * b;
  MBCartVect par_pos = B * (b - center);
  MBCartVect par_dir = B * m;
  
    //fast rejection test
  const double half_x = length[0] + reps;
  if ((par_pos[0] >  half_x && par_dir[0] >= 0) ||
      (par_pos[0] < -half_x && par_dir[0] <= 0))
    return false;
  
  const double half_y = length[1] + reps;
  if ((par_pos[1] >  half_y && par_dir[1] >= 0) ||
      (par_pos[1] < -half_y && par_dir[1] <= 0))
    return false;
    
  const double half_z = length[2] + reps;
  if ((par_pos[2] >  half_z && par_dir[2] >= 0) ||
      (par_pos[2] < -half_z && par_dir[2] <= 0))
    return false;
  
    // test if point is inside
  if (par_pos[0] <= half_x && par_pos[0] >= -half_x &&
      par_pos[1] <= half_y && par_pos[1] >= -half_y &&
      par_pos[2] <= half_z && par_pos[2] >= -half_z)
    return true;
  
  // then outside case
//bool write = false;
//if (write) {
//  FILE* file = fopen("dump","w+");
//  fprintf(file,"create vertex %f %f %f\n", -half_x, -half_y, -half_z );
//  fprintf(file,"create vertex %f %f %f\n",  half_x, -half_y, -half_z );
//  fprintf(file,"create vertex %f %f %f\n",  half_x,  half_y, -half_z );
//  fprintf(file,"create vertex %f %f %f\n", -half_x,  half_y, -half_z );
//  fprintf(file,"create vertex %f %f %f\n", -half_x, -half_y, +half_z );
//  fprintf(file,"create vertex %f %f %f\n",  half_x, -half_y, +half_z );
//  fprintf(file,"create vertex %f %f %f\n",  half_x,  half_y, +half_z );
//  fprintf(file,"create vertex %f %f %f\n", -half_x,  half_y, +half_z );
//  fprintf(file,"create surface vertex 1 2 6 5\n");
//  fprintf(file,"create surface vertex 2 3 7 6\n");
//  fprintf(file,"create surface vertex 3 4 8 7\n");
//  fprintf(file,"create surface vertex 4 1 5 8\n");
//  fprintf(file,"create surface vertex 4 3 2 1\n");
//  fprintf(file,"create surface vertex 5 6 7 8\n");
//  fprintf(file,"create volume surface all\n");
//  fprintf(file,"delete vertex all\n");
//  fprintf(file,"compress ids\n");
//  fprintf(file,"create vertex %f %f %f\n", par_pos[0], par_pos[1], par_pos[2]);
//  fprintf(file,"create vertex %f %f %f\n", par_pos[0] + 1000 * par_dir[0],
//                                           par_pos[1] + 1000 * par_dir[1],
//                                           par_dir[2] + 1000 * par_dir[2] );
//  fprintf(file,"create curve vertex 9 10\n");
//  fclose(file);
//}

    //test two xy plane
  if ((half_z - par_pos[2]) * par_dir[2] >= 0 &&
      fabs(par_dir[0] * (half_z - par_pos[2]) + par_dir[2] * par_pos[0]) 
        <= fabs(par_dir[2] * half_x) && 
      fabs(par_dir[1] * (half_z - par_pos[2]) + par_dir[2] * par_pos[1]) 
        <= fabs(par_dir[2] * half_y)) 
    return true;
  if ((-half_z - par_pos[2]) * par_dir[2] >= 0 &&
      fabs(par_dir[0] * (-half_z - par_pos[2]) + par_dir[2] * par_pos[0]) 
        <= fabs(par_dir[2] * half_x) && 
      fabs(par_dir[1] * (-half_z - par_pos[2]) + par_dir[2] * par_pos[1]) 
        <= fabs(par_dir[2] * half_y))
    return true;

    //test two xz plane
  if ((half_y - par_pos[1]) * par_dir[1] >= 0 &&
      fabs(par_dir[0] * (half_y - par_pos[1]) + par_dir[1] * par_pos[0]) 
        <= fabs(par_dir[1] * half_x) && 
      fabs(par_dir[2] * (half_y - par_pos[1]) + par_dir[1] * par_pos[2]) 
        <= fabs(par_dir[1] * half_z))
    return true;
  if ((-half_y - par_pos[1]) * par_dir[1] >= 0 &&
      fabs(par_dir[0] * (-half_y - par_pos[1]) + par_dir[1] * par_pos[0]) 
        <= fabs(par_dir[1] * half_x)  && 
      fabs(par_dir[2] * (-half_y - par_pos[1]) + par_dir[1] * par_pos[2])
        <= fabs(par_dir[1] * half_z))
    return true;

    //test two yz plane
  if ((half_x - par_pos[0]) * par_dir[0] >= 0 &&
      fabs(par_dir[1] * (half_x - par_pos[0]) + par_dir[0] * par_pos[1]) 
        <= fabs(par_dir[0] * half_y) &&
      fabs(par_dir[2] * (half_x - par_pos[0]) + par_dir[0] * par_pos[2]) 
        <= fabs(par_dir[0] * half_z))
    return true;
  if ((-half_x - par_pos[0]) * par_dir[0] >= 0 &&
      fabs(par_dir[1] * (-half_x - par_pos[0]) + par_dir[0] * par_pos[1])
        <= fabs(par_dir[0] * half_y) &&
      fabs(par_dir[2] * (-half_x - par_pos[0]) + par_dir[0] * par_pos[2]) 
        <= fabs(par_dir[0] * half_z))
    return true;

  return false;
}

MBErrorCode MBOrientedBox::make_hex( MBEntityHandle& hex, MBInterface* instance )
{
  MBErrorCode rval;
  int signs[8][3] = { { -1, -1, -1 },
                      {  1, -1, -1 },
                      {  1,  1, -1 },
                      { -1,  1, -1 },
                      { -1, -1,  1 },
                      {  1, -1,  1 },
                      {  1,  1,  1 },
                      { -1,  1,  1 } };
                      
  std::vector<MBEntityHandle> vertices;
  for (int i = 0; i < 8; ++i)
  {
    MBCartVect coords(center);
    for (int j = 0; j < 3; ++j)
      coords += signs[i][j] * axis[j];
    MBEntityHandle handle;
    rval = instance->create_vertex( coords.array(), handle );
    if (MB_SUCCESS != rval) {
      instance->delete_entities( &vertices[0], vertices.size() );
      return rval;
    }
    vertices.push_back( handle );
  }
  
  rval = instance->create_element( MBHEX, &vertices[0], vertices.size(), hex );
  if (MB_SUCCESS != rval) {
    instance->delete_entities( &vertices[0], vertices.size() );
    return rval;
  }
  
  return MB_SUCCESS;
}
  
