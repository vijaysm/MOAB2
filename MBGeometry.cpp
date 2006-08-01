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
#include <cmath>

namespace MBGeometry {

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
  
  if (beta < 0 || gamma < 0)
    return false;
  if (beta + gamma > 1)
    return false;
  if (t < 0.0)
    return false;
  
  t_out = t;
  return true;
}

} // namespace MBGeoemtry
