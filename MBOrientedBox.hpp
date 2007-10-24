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

/**\file MBOrientedBox.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */

#ifndef MB_ORIENTED_BOX_HPP
#define MB_ORIENTED_BOX_HPP

#include "MBForward.hpp"
#include "MBCartVect.hpp"
#include "MBMatrix3.hpp"

#include <iosfwd>

#define MB_ORIENTED_BOX_UNIT_VECTORS 1
#define MB_ORIENTED_BOX_OUTER_RADIUS 1

class MBRange;


/**\brief Oriented bounding box
 */
class MBOrientedBox
{
public:
  MBCartVect center;  //!< Box center
  MBCartVect axis[3]; //!< Box axes, unit vectors sorted by extent of box along axis
#if MB_ORIENTED_BOX_UNIT_VECTORS
  MBCartVect length;  //!< distance from center to plane along each axis
#endif
#if MB_ORIENTED_BOX_OUTER_RADIUS
  double radius;      //!< outer radius (1/2 diagonal length) of box
#endif

  inline MBOrientedBox() {}

  MBOrientedBox( const MBCartVect axis[3], const MBCartVect& center );

  inline double inner_radius() const; //!< radius of inscribed sphere
  inline double outer_radius() const; //!< radius of circumscribed sphere
  inline double outer_radius_squared() const; //!< square of radius of circumsphere
  inline double inner_radius_squared() const; //!< square of radius if inscribed sphere
  inline double volume() const;               //!< volume of box
  inline MBCartVect dimensions() const;       //!< number of dimensions for which box is not flat
  inline double area() const;                 //!< largest side area
  inline MBCartVect scaled_axis( int index ) const; //!< get vector in direction of axis, from box center to face
  
  /** Test if point is contained in box */
  bool contained( const MBCartVect& point, double tolerance ) const;
  
  //bool contained( const MBOrientedBox& other, double tolerance ) const;
  
  /**\brief get tag handle for storing oriented box
   *
   * Get the handle for the tag with the specified name and
   * check that the tag is appropriate for storing instances
   * of MBOrientedBox.  The resulting tag may be used to store
   * instances of MBOrientedBox directly.
   *
   *\param handle_out  The TagHandle, passed back to caller
   *\param name        The tag name
   *\param create      If true, tag will be created if it does not exist
   */
  static MBErrorCode tag_handle( MBTag& handle_out,
                                 MBInterface* instance, 
                                 const char* name,
                                 bool create = true );

  /**\brief Calculate an oriented box from a set of vertices */
  static MBErrorCode compute_from_vertices( MBOrientedBox& result,
                                            MBInterface* instance,
                                            const MBRange& vertices );
                                  
  /**\brief Calculate an oriented box from a set of 2D elements */
  static MBErrorCode compute_from_2d_cells( MBOrientedBox& result,
                                            MBInterface* instance,
                                            const MBRange& elements );

    /** Structure to hold temporary accumulated triangle data for
     *  caculating box orietation.  See box_from_covariance_data
     *  to see how this is used to calculate the final covariance matrix
     *  and resulting box orientation.
     */
  struct CovarienceData {
    MBMatrix3 matrix;    //!< Running sum for covariance matrix
    MBCartVect center;   //!< Sum of triangle centroids weighted by 2*triangle area
    double area;         //!< 2x the sum of the triangle areas
  };
  
    /** Calculate a CovarienceData struct from a list of triangles */
  static MBErrorCode covariance_data_from_tris( CovarienceData& result,
                                                MBInterface* moab_instance,
                                                const MBRange& elements );
  
    /** Calculate an MBOrientedBox given an arrray of CovarienceData and 
     *  the list  of vertices the box is to bound.
     */
  static MBErrorCode compute_from_covariance_data( MBOrientedBox& result,
                                          MBInterface* moab_instance,
                                          const CovarienceData* orient_array,
                                          unsigned orient_array_length,
                                          const MBRange& vertices );
  
    /** Test for intersection of a ray (or line segment) with this box
     *\param ray_start_point The base point of the ray
     *\param ray_unit_direction The direction of the ray (must be unit length)
     *\param distance_tolerance Tolerance to use in intersection checks
     *\param segment_length Optional length of ray
     */
  bool intersect_ray( const MBCartVect& ray_start_point,
                      const MBCartVect& ray_unit_direction,
                      double distance_tolerance,
                      const double* segment_length = 0 ) const;
                      
    /**\brief Find closest position on/within box to input position.
     * 
     * Find the closest position in the solid box to the input position.
     * If the input position is on or within the box, then the output
     * position will be the same as the input position.  If the input
     * position is outside the box, the outside position will be the
     * closest point on the box boundary to the input position.
     */
  void closest_location_in_box( const MBCartVect& input_position,
                                MBCartVect& output_position ) const;
                      
    //! Construct a hexahedral element with the same shape as this box.
  MBErrorCode make_hex( MBEntityHandle& hex, MBInterface* instance );
                                    
  
    /** Calculate an MBOrientedBox given a CovarienceData struct and
     *  the list of points the box is to bound.
     */
  static MBErrorCode compute_from_covariance_data( MBOrientedBox& result,
                                          MBInterface* moab_instance,
                                          CovarienceData& orientation_data,
                                          const MBRange& vertices );
};

std::ostream& operator<<( std::ostream&, const MBOrientedBox& );

double MBOrientedBox::inner_radius() const
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return length[0];
#else
  return axis[0].length();
#endif
}

double MBOrientedBox::outer_radius() const
{
#if MB_ORIENTED_BOX_OUTER_RADIUS
  return radius;
#elif MB_ORIENTED_BOX_UNIT_VECTORS
  return length.length();
#else
  return (axis[0] + axis[1] + axis[2]).length();
#endif
}

double MBOrientedBox::outer_radius_squared() const
{
#if MB_ORIENTED_BOX_OUTER_RADIUS
  return radius * radius;
#elif MB_ORIENTED_BOX_UNIT_VECTORS
  return length % length;
#else
  const MBCartVect half_diag = axis[0] + axis[1] + axis[2];
  return half_diag % half_diag;
#endif
}

double MBOrientedBox::inner_radius_squared() const
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return length[0] * length[0];
#else
  return axis[0] % axis[0];
#endif
}

double MBOrientedBox::volume() const
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return 8 * length[0] * length[1] * length[2];
#else
  return fabs(8 * axis[0] % (axis[1] * axis[2]));
#endif
}

MBCartVect MBOrientedBox::dimensions() const
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return 2.0 * length;
#else
  return 2.0 * MBCartVect( axis[0].length(), axis[1].length(), axis[2].length() );
#endif
}

double MBOrientedBox::area() const
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return 4 * length[1] * length[2];
#else
  return 4 * (axis[1] * axis[2]).length();
#endif
}

MBCartVect MBOrientedBox::scaled_axis( int index ) const
{
#if MB_ORIENTED_BOX_UNIT_VECTORS
  return length[index] * axis[index];
#else
  return axis[index];
#endif
}


#endif
