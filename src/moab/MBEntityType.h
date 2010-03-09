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

#ifndef MB_ENTITY_TYPE_HPP
#define MB_ENTITY_TYPE_HPP

#ifdef __cplusplus
extern "C" {
#endif

/*! Entity types defined in MOAB and MBCN
 *  The ordering here must ensure that all element types are 
 *  grouped together and all elements of similar dimension are
 *  grouped together.
 */
typedef enum 
{
  MBVERTEX = 0, /**< Mesh Vertex AKA node */
  MBEDGE,       /**< Mesh Edge */
  MBTRI,        /**< Triangular element (including shells) */
  MBQUAD,       /**< Quadrilateral element (including shells) */
  MBPOLYGON,    /**< Polygon */
  MBTET,        /**< Tetrahedral element */
  MBPYRAMID,    /**< Pyramid element (where are the face ids for this defined?) */
  MBPRISM,      /**< Wedge element (Exodus has one, Cubit doesn't. Does Mesh need it?) */
  MBKNIFE,      /**< Knife element */
  MBHEX,        /**< Hexahedral element */
  MBPOLYHEDRON, /**< Polyhedron */
  MBENTITYSET,    /**< MeshSet */
  MBMAXTYPE  /**< Just a place keeper - must be the # of entities, for array */
    /**< dimensioning purposes  */
} MBEntityType;

#ifdef __cplusplus
} /* extern "C" */


/** prefix increment operator for MBEntityType */
inline MBEntityType & operator++(MBEntityType &type)
{
  return type = static_cast<MBEntityType>(type+1);
}

/** postfix increment operator for MBEntityType */
inline MBEntityType operator++(MBEntityType &type, int)
{
  MBEntityType oldval = type;
  ++type;
  return oldval;
}

/** prefix increment operator for MBEntityType */
inline MBEntityType & operator--(MBEntityType &type)
{
  return type = static_cast<MBEntityType>(type-1);
}

/** postfix increment operator for MBEntityType */
inline MBEntityType operator--(MBEntityType &type, int)
{
  MBEntityType oldval = type;
  --type;
  return oldval;
}

#endif

#endif
