/**
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

//-------------------------------------------------------------------------
// Filename      : MBUtil.hpp
//
// Purpose       : This file contains utility functions that can be used
//                 with MB
//
// Special Notes : This is a pure virtual class, to prevent instantiation.
//                 All functions are static, called like this:
//                 MBUtil::function_name();
//
// Creator       : Ray J. Meyers
//
// Date          : 09/01/02
//
// Owner         : Ray J. meyers
//-------------------------------------------------------------------------


#ifndef MB_UTIL_HPP
#define MB_UTIL_HPP

#include "MBForward.hpp"

struct  MBCoord
{
  double x;
  double y;
  double z;
};


//! utility function to be use throughout VERDE
class MBUtil
{
public:

   
  static void normal(MBInterface* MB, MBEntityHandle handle, double& x, double& y, double& z);

  static void centroid(MBInterface *MB, MBEntityHandle handle,MBCoord &coord);
 // static void edge_centers(MBInterface *MB, MBEntityHandle handle, std::vector<MBCoord> &coords_list);

  //static void face_centers(MBInterface *MB, MBEntityHandle handle, std::vector<MBCoord> &coords_list);

private:

  MBUtil(){};

};

#endif
