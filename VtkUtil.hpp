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

#ifndef VTK_UTIL
#define VTK_UTIL

//
// ExoIIUtil class: utility class for functions used by both reader
// and writer

#ifndef IS_BUILDING_MB
#error "ExoIIUtil.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"


class VtkUtil 
{

public:

    //! vtk element type numbers
  static const int vtkElemType[];

    //! vtk element type names
  static const char *vtkElemNames[];
};

#endif
