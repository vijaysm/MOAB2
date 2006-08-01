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

#ifndef MB_GEOMETRY_HPP
#define MB_GEOMETRY_HPP

#include "MBCartVect.hpp"
#include <cmath>

namespace MBGeometry {

bool ray_tri_intersect( const MBCartVect vertices[3],
                        const MBCartVect& ray_point,
                        const MBCartVect& ray_unit_direction,
                        double tolerance,
                        double& t_out,
                        const double* ray_length = 0 );

} // namespace MBGeoemtry

#endif
