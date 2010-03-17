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

#ifndef MB_ENTITY_TYPE_H
#define MB_ENTITY_TYPE_H

/*! Entity types defined in MOAB and MBCN
 *  The ordering here must ensure that all element types are 
 *  grouped together and all elements of similar dimension are
 *  grouped together.
 */

#ifdef __cplusplus
namespace moab { 
# define MOAB_ENTITY_TYPE_NAME EntityType
# else /* __cplusplus */
# define MOAB_ENTITY_TYPE_NAME MBEntityType
#endif /* __cplusplus */
enum MOAB_ENTITY_TYPE_NAME
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
};

/** prefix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME & operator++(MOAB_ENTITY_TYPE_NAME &type)
{
  return type = static_cast<MOAB_ENTITY_TYPE_NAME>(type+1);
}

/** postfix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME operator++(MOAB_ENTITY_TYPE_NAME &type, int)
{
  MOAB_ENTITY_TYPE_NAME oldval = type;
  ++type;
  return oldval;
}

/** prefix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME & operator--(MOAB_ENTITY_TYPE_NAME &type)
{
  return type = static_cast<MOAB_ENTITY_TYPE_NAME>(type-1);
}

/** postfix increment operator for MBEntityType */
inline MOAB_ENTITY_TYPE_NAME operator--(MOAB_ENTITY_TYPE_NAME &type, int)
{
  MOAB_ENTITY_TYPE_NAME oldval = type;
  --type;
  return oldval;
}

#ifdef __cplusplus
} /* namespace moab*/
#endif /* __cplusplus */

#undef MOAB_ENTITY_TYPE_NAME
#endif /* MB_ENTITY_TYPE_H */

#ifdef __cplusplus
#  ifndef MOAB_ENTITY_TYPE_NS_ONLY
#    define MOAB_ENTITY_TYPE_NS_ONLY
     typedef moab::EntityType MBEntityType;
     using moab::MBVERTEX;
     using moab::MBEDGE;
     using moab::MBTRI;
     using moab::MBQUAD;
     using moab::MBPOLYGON;
     using moab::MBTET;
     using moab::MBPYRAMID;
     using moab::MBPRISM;
     using moab::MBKNIFE;
     using moab::MBHEX;
     using moab::MBPOLYHEDRON;
     using moab::MBENTITYSET;
     using moab::MBMAXTYPE;
#  endif /* MOAB_ENTITY_TYPE_NS_ONLY */
#else /* __cplusplus */
#  ifndef MOAB_ENTITY_TYPE_C
#    define MOAB_ENTITY_TYPE_C
     typedef enum MBEntityType MBEntityType;
#  endif /* MOAB_ENTITY_TYPE_C */
#endif /* __cplusplus */
