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


const char* VtkUtil::vtkTypeNames[] = {
 "unsigned_char", // MB_TYPE_OPAQUE
 "int",           // MB_TYPE_INTEGER
 "double",        // MB_TYPE_DOUBLE
 "bit",           // MB_TYPE_BIT
 "unsigned_long", // MB_TYPE_HANDLE
};

/*
const unsigned VtkUtil::typeSizes[] = {
 1,                      // MB_TYPE_OPAQUE
 sizeof(int),            // MB_TYPE_INTEGER
 sizeof(double),         // MB_TYPE_DOUBLE
 1,                      // MB_TYPE_BIT
 sizeof(MBEntityHandle), // MB_TYPE_HANDLE
};
*/

const unsigned pixel[] = { 0, 1, 3, 2 };
const unsigned voxel[] = { 0, 1, 3, 2, 4, 5, 7, 6 };
const unsigned wedge[] = { 0, 2, 1, 3, 5, 4 };
const unsigned  qhex[] = {  0,  1,  2,  3, 
                            4,  5,  6,  7, 
                            8,  9, 10, 11,
                           16, 17, 18, 19,
                           12, 13, 14, 15 };

const VtkElemType VtkUtil::vtkElemTypes[] = {
      { 0,                0, MBMAXTYPE, 0, 0 },
      { "vertex",         1, MBMAXTYPE, 1, 0 },
      { "polyvertex",     2, MBMAXTYPE, 0, 0 },
      { "line",           3, MBEDGE,    2, 0 },
      { "polyline",       4, MBMAXTYPE, 0, 0 },
      { "triangle",       5, MBTRI,     3, 0 },
      { "triangle strip", 6, MBMAXTYPE, 0, 0 },
      { "polygon",        7, MBPOLYGON, 0, 0 },
      { "pixel",          8, MBQUAD,    4, pixel },
      { "quadrilateral",  9, MBQUAD,    4, 0 }, 
      { "tetrahedron",   10, MBTET,     4, 0 }, 
      { "voxel",         11, MBHEX,     8, voxel }, 
      { "hexahedron",    12, MBHEX,     8, 0 }, 
      { "wedge",         13, MBPRISM,   6, wedge }, 
      { "pyramid",       14, MBPYRAMID, 5, 0 },
      { 0,               15, MBMAXTYPE, 0, 0 },
      { 0,               16, MBMAXTYPE, 0, 0 },
      { 0,               17, MBMAXTYPE, 0, 0 },
      { 0,               18, MBMAXTYPE, 0, 0 },
      { 0,               19, MBMAXTYPE, 0, 0 },
      { 0,               20, MBMAXTYPE, 0, 0 },
      { "quadratic edge",21, MBEDGE,    3, 0 },
      { "quadratic tri", 22, MBTRI,     6, 0 },
      { "quadratic quad",23, MBQUAD,    8, 0 },
      { "quadratic tet", 24, MBTET,    10, 0 },
      { "quadratic hex", 25, MBHEX,    20, qhex } ,
      { 0,               26, MBMAXTYPE, 0, 0 }
    };

const unsigned VtkUtil::numVtkElemType = sizeof(VtkUtil::vtkElemTypes) / sizeof(VtkUtil::vtkElemTypes[0]);

const VtkElemType* mb_type_offsets[] = {
  NULL,                       // MBVERTEX
  VtkUtil::vtkElemTypes +  3, // MBEDGE
  VtkUtil::vtkElemTypes +  5, // MBTRI
  VtkUtil::vtkElemTypes +  9, // MBQUAD
  VtkUtil::vtkElemTypes +  7, // MBPOLYGON
  VtkUtil::vtkElemTypes + 10, // MBTET
  VtkUtil::vtkElemTypes + 14, // MBPYRAMID
  VtkUtil::vtkElemTypes + 13, // MBWEDGE
  NULL,                       // MBKNIFE
  VtkUtil::vtkElemTypes + 12, // MBHEX
  NULL,                       // MBPOLYHEDRON
  NULL,                       // MBENTITYSET
  NULL,                       // MBMAXTYPE
};

const VtkElemType* VtkUtil::get_vtk_type( MBEntityType type, unsigned num_nodes )
{
  const VtkElemType* iter = mb_type_offsets[type];
  if (!iter || type == MBPOLYGON)
    return iter;
  
  do {
    if (iter->num_nodes == num_nodes)
      return iter;
    ++iter;
  } while (iter->mb_type == type);
  
  return 0;
}

