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

#include "VtkUtil.hpp"
#include "MBInternals.hpp"
#include "MBInterface.hpp"

const int VtkUtil::vtkElemType[] =
{
    1,  //  MBVERTEX
    3,  //  MBEDGE,     
    5,  //  MBTRI,      
    9,  //  MBQUAD,     
    7,  //  MBPOLYGON,  
    10, //  MBTET,      
    14, //  MBPYRAMID,  
    13, //  MBPRISM,    
    -1, //  MBKNIFE,    
    12, //  MBHEX,      
    -1, //  MBPOLYHEDRON
    -1, //  MBENTITYSET,
    -1  //  MBMAXTYPE    
};

const char *VtkUtil::vtkElemNames[] = 
{
  "VTK_VERTEX", //  MBVERTEX     
  "VTK_LINE", //  MBEDGE,      
  "VTK_TRIANGLE", //  MBTRI,       
  "VTK_QUAD", //  MBQUAD,      
  "VTK_POLYGON", //  MBPOLYGON,   
  "VTK_TETRA", //  MBTET,       
  "VTK_PYRAMID", //  MBPYRAMID,   
  "VTK_WEDGE", //  MBPRISM,     
  "VTK_NO_EQUIVALENT", //  MBKNIFE,     
  "VTK_HEXAHEDRON", //  MBHEX,       
  "VTK_NO_EQUIVALENT", //  MBPOLYHEDRON 
  "VTK_NO_EQUIVALENT", //  MBENTITYSET, 
  "VTK_NO_EQUIVALENT"  //  MBMAXTYPE    
};
