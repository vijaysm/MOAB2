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

/**\file MBGeometry.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-27
 */

#ifndef MB_GEOM_UTIL_HPP
#define MB_GEOM_UTIL_HPP

#include "MBCartVect.hpp"
#include <cmath>

namespace MBGeomUtil {

bool ray_tri_intersect( const MBCartVect vertices[3],
                        const MBCartVect& ray_point,
                        const MBCartVect& ray_unit_direction,
                        double tolerance,
                        double& t_out,
                        const double* ray_length = 0 );

/**\brief find closest location on triangle
 *
 * Find closest location on linear triangle.
 *\param location  Input position to evaluate from
 *\param vertices  Array of three corner vertex coordinates.
 *\param closest_out Result position 
 */
void closest_location_on_tri( const MBCartVect& location,
                              const MBCartVect* vertices,
                              MBCartVect& closest_out );

/**\brief find closest location on polygon
 *
 * Find closest location on polygon
 *\param location  Input position to evaluate from
 *\param vertices  Array of corner vertex coordinates.
 *\param num_vertices Length of 'vertices' array.
 *\param closest_out Result position 
 */
bool closest_location_on_polygon( const MBCartVect& location,
                                  const MBCartVect* vertices,
                                  int num_vertices,
                                  MBCartVect& closest_out );

/**\brief find closest topological location on triangle
 *
 * Find closest location on linear triangle.
 *\param location  Input position to evaluate from
 *\param vertices  Array of three corner vertex coordinates.
 *\param tolerance Tolerance to use when comparing to corners and edges
 *\param closest_out Result position 
 *\param closest_topo Closest topological entity
 *                     0-2 : vertex index
 *                     3-5 : edge beginning at closest_topo - 3
 *                       6 : triangle interior
 */
void closest_location_on_tri( const MBCartVect& location,
                              const MBCartVect* vertices,
                              double tolerance,
                              MBCartVect& closest_out,
                              int& closest_topo );

} // namespace MBGeoemtry

#endif
