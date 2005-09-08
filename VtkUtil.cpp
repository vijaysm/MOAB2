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

// Define node ordering for those types for which
// the VTK and MOAB ordering doesn't match.
const unsigned pixel[] = { 0, 1, 3, 2 };
const unsigned voxel[] = { 0, 1, 3, 2, 4, 5, 7, 6 };
const unsigned wedge[] = { 0, 2, 1, 3, 5, 4 };
const unsigned  qhex[] = {  0,  1,  2,  3, 
                            4,  5,  6,  7, 
                            8,  9, 10, 11,
                           16, 17, 18, 19,
                           12, 13, 14, 15 };

// List of VtkElemType structs, indexed by the VTK type number.
const VtkElemType VtkUtil::vtkElemTypes[] = {
      { 0,                0, MBMAXTYPE, 0, 0 },
      { "vertex",         1, MBMAXTYPE, 1, 0 },
      { "polyvertex",     2, MBMAXTYPE, 0, 0 },
      { "line",           3, MBEDGE,    2, 0 },
      { "polyline",       4, MBMAXTYPE, 0, 0 },
      { "triangle",       5, MBTRI,     3, 0 },
      { "triangle strip", 6, MBMAXTYPE, 0, 0 },
//      { "polygon",        7, MBPOLYGON, 0, 0 },
      { 0,                7, MBMAXTYPE, 0, 0 },
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

// Define an array, indexed by MBEntityType containing the corresponding 
// VTK element type numbers for the linear and quadratic (mid-edge) elements.
// Zero is used to indicate an invalid type (not supported by VTK.)  The
// VTK element type number may be used as an index into vtkElemTypes[].
const int mb_to_vtk_type[][2] = {
  {  1,  0 },  // MBVERTEX
  {  3, 21 },  // MBEDGE
  {  5, 22 },  // MBTRI
  {  9, 23 },  // MBQUAD
  {  0,  0 },  // MBPOLYGON
  { 10, 24 },  // MBTET
  { 14,  0 },  // MBPYRAMID
  { 13,  0 },  // MBWEDGE
  {  0,  0 },  // MBKNIFE
  { 12, 25 },  // MBHEX
  {  0,  0 },  // MBPOLYHEDRON
  {  0,  0 },  // MBENTITYSET
  {  0,  0 } };// MBMAXTYPE

const VtkElemType* VtkUtil::get_vtk_type( MBEntityType type, unsigned num_nodes )
{
  const int i = mb_to_vtk_type[type][0]; // Index for linear type
  const int j = mb_to_vtk_type[type][1]; // Index for quadratic type
  if (i) // If element type is supported at all (if not linear then not quadratic either)
  {
      // If the linear type is requested (all polygons are linear
      // irrespective of the number of nodes), return that.
    if (type == MBPOLYGON || vtkElemTypes[i].num_nodes == num_nodes)
      return vtkElemTypes + i;
      // Otherwise if there is a quadratic type and the number of
      // nodes specified corresponds to the quadratic type, return that.
    else if (j && vtkElemTypes[j].num_nodes == num_nodes)
      return vtkElemTypes + j;
  }
  
  return 0;
}

