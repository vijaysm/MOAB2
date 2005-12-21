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

/** 
 * \class MBCN
 * \author Tim Tautges
 * \date April 2004
 *
 * \brief Canonical numbering data and functions
 * This class represents canonical ordering of finite-element meshes.
 * Elements in the finite element "zoo" are represented.  Canonical numbering
 * denotes the vertex, edge, and face numbers making up each kind of element,
 * and the vertex numbers defining those entities.  Functions for evaluating
 * adjacencies and other things based on vertex numbering are also provided.
 * By default, this class defines a zero-based numbering system.
 * For a complete description of this class, see the document "MOAB Canonical
 * Numbering Conventions", Timothy J. Tautges, Sandia National Laboratories
 * Report #SAND2004-xxxx.
 */
//
#ifndef MBCN_HPP
#define MBCN_HPP

#include <vector>
#include <algorithm>

// the maximum number n-1 dimension adjacencies a element may have
#define MB_MAX_SUB_ENTITIES  24

// the maximum number of nodes an n-1 dimensional element may have
#define MB_MAX_SUB_ENTITY_VERTICES 16

// the maximum possible number of vertices, including higher-order
#define MB_MAX_NUM_VERTICES 27

/*! Entity types defined in MOAB and MBCN
 *  The ordering here must ensure that all element types are 
 *  grouped together and all elements of similar dimension are
 *  grouped together.
 */
enum MBEntityType
{
  MBVERTEX = 0, //!< Mesh Vertex AKA node
  MBEDGE,       //!< Mesh Edge
  MBTRI,        //!< Triangular element (including shells)
  MBQUAD,       //!< Quadrilateral element (including shells)
  MBPOLYGON,    //!< Polygon
  MBTET,        //!< Tetrahedral element
  MBPYRAMID,    //!< Pyramid element (where are the face ids for this defined?)
  MBPRISM,      //!< Wedge element (Exodus has one, Cubit doesn't. Does Mesh need it?)
  MBKNIFE,      //!< Knife element
  MBHEX,        //!< Hexahedral element
  MBPOLYHEDRON, //!< Polyhedron
  MBENTITYSET,    //!< MeshSet
  MBMAXTYPE  //!< Just a place keeper - must be the # of entities, for array
    //!< dimensioning purposes 
};

typedef void* Handle;

//! postfix increment operator for MBEntityType
inline MBEntityType operator++(MBEntityType &type, int)
{
  return static_cast<MBEntityType>(reinterpret_cast<int&>(type)++);
}  

//! prefix increment operator for MBEntityType
inline MBEntityType & operator++(MBEntityType &type)
{
  return reinterpret_cast<MBEntityType&>(++reinterpret_cast<int&>(type));
}

typedef std::pair<MBEntityType, MBEntityType> MBDimensionPair;

class MBCN
{
private:

//! entity names
  static const char *entityTypeNames[];
  
//! declare private constructor, since we don't want to create any of these
  MBCN();

//! the basis of the numbering system (normally 0 or 1, 0 by default)
  static int numberBasis;

public:

    //! enum used to specify operation type
  enum {INTERSECT, UNION};

    // each entity type has two ConnMap objects, holding information about the bounding
    // edges and faces for each entity; see comment for mConnectivityMap
    // (this struct not documented with Doxygen)
  struct ConnMap
  {
      // Topological dimension of this entry
    int topo_dimension;

      // Number of sub-elements of this dimension
    int num_sub_elements;
    
      // Number of nodes in each sub-element of this dimension
    int num_nodes_per_sub_element[MB_MAX_SUB_ENTITIES];

      // Type of each sub-element
    MBEntityType target_type[MB_MAX_SUB_ENTITIES];

      // Connectivity of each of the sub-elements
    int conn[MB_MAX_SUB_ENTITIES][MB_MAX_SUB_ENTITY_VERTICES];
  };

    // mConnectivityMap[i=entity type][j=0,1,2]:
    //  num_sub_elements = # bounding edges(j=0) or faces(j=1) for entity type i, or self (j=2)
    //  num_nodes_per_sub_element[k] (k=0..num_sub_elements-1) = number of nodes in sub-facet k
    //    (can vary over sub-facets, e.g. faces bounding a pyramid) or self (j=2)
    //  target_type[k] = entity type of sub-facet k (e.g. MBTRI or MBQUAD bounding a pyramid) or self (j=2)
    //  conn[k][l] (l=0..MBCN::VerticesPerEntity[target_type[k]]) = vertex connectivity of sub-facet k,
    //    with respect to entity i's canonical vertex ordering, or self (j=2)
    // (not documented with Doxygen)
  static ConnMap mConnectivityMap[MBMAXTYPE][3];

    // structure used to define reverse canonical ordering information
    // (not documented with Doxygen)
  struct UpConnMap
  {
      // Number of higher-dimensional entities using each sub-entity
    int num_targets_per_source_element[MB_MAX_SUB_ENTITIES];

      // Higher-dimensional entities using each sub-entity
    int targets_per_source_element[MB_MAX_SUB_ENTITIES][MB_MAX_SUB_ENTITIES];
  };

    // Reverse canonical numbering, duplicates data in mConnectivityMap, but 
    // connectivity data in this table must be in ascending order (used for
    // efficient sorting)
    // (not documented with Doxygen)
  static UpConnMap mUpConnMap[MBMAXTYPE][4][4];

  //! this const vector defines the starting and ending MBEntityType for 
  //! each dimension, e.g. TypeDimensionMap[2] returns a pair of MBEntityTypes 
  //! bounding dimension 2.
  static const MBDimensionPair TypeDimensionMap[];

  //! get the basis of the numbering system
  static int GetBasis();
  
  //! set the basis of the numbering system
  static void SetBasis(const int in_basis);

  //! return the string type name for this type
  static const char *EntityTypeName(const MBEntityType this_type);
  
  //! given a name, find the corresponding entity type
  static MBEntityType EntityTypeFromName(const char *name);
  
  //! return the topological entity dimension
  static int Dimension(const MBEntityType t);

  //! return the number of (corner) vertices contained in the specified type.  
  static int VerticesPerEntity(const MBEntityType t);
  
  //! return the number of sub-entities bounding the entity.
  static int NumSubEntities(const MBEntityType t, const int d);

  //! return the type of a particular sub-entity.
  //! \param this_type Type of entity for which sub-entity type is being queried
  //! \param sub_dimension Topological dimension of sub-entity whose type is being queried
  //! \param index Index of sub-entity whose type is being queried
  //! \return type Entity type of sub-entity with specified dimension and index
  static MBEntityType SubEntityType(const MBEntityType this_type,
                                         const int sub_dimension,
                                         const int index);
  
  //! return the vertex indices of the specified sub-entity.
  //! \param this_type Type of entity for which sub-entity connectivity is being queried
  //! \param sub_dimension Dimension of sub-entity
  //! \param sub_index Index of sub-entity
  //! \param sub_entity_conn Connectivity of sub-entity (returned to calling function)
  static void SubEntityVertexIndices(const MBEntityType this_type, 
                                     const int sub_dimension,
                                     const int sub_index,
                                     int sub_entity_conn[]);

  //! return the vertices of the specified sub entity
  //! \param parent_conn Connectivity of parent entity
  //! \param parent_type Entity type of parent entity
  //! \param sub_dimension Dimension of sub-entity being queried
  //! \param sub_index Index of sub-entity being queried
  //! \param sub_entity_conn Connectivity of sub-entity, based on parent_conn and canonical
  //!           ordering for parent_type
  //! \param num_sub_vertices Number of vertices in sub-entity
  static void SubEntityConn(const void *parent_conn, const MBEntityType parent_type,
                            const int sub_dimension,
                            const int sub_index,
                            void *sub_entity_conn, int &num_sub_vertices);

  //! For a specified set of sides of given dimension, return the intersection 
  //! or union of all sides of specified target dimension adjacent to those sides.
  //! \param this_type Type of entity for which sub-entity connectivity is being queried
  //! \param source_indices Indices of sides being queried
  //! \param num_source_indices Number of entries in <em>source_indices</em>
  //! \param source_dim Dimension of source entity
  //! \param target_dim Dimension of target entity
  //! \param index_list Indices of target entities (returned)
  //! \param operation_type Specify either MBCN::INTERSECT or MBCN::UNION to get intersection
  //!        or union of target entity lists over source entities
  static int AdjacentSubEntities(const MBEntityType this_type,
                                 const int *source_indices,
                                 const int num_source_indices,
                                 const int source_dim,
                                 const int target_dim,
                                 std::vector<int> &index_list,
                                 const int operation_type = MBCN::INTERSECT);

  //! return the side index represented in the input sub-entity connectivity in the input 
  //! parent entity connectivity array.
  //! \param parent_conn Connectivity of parent entity being queried
  //! \param parent_type Entity type of parent entity
  //! \param child_conn Connectivity of child whose index is being queried
  //! \param child_num_verts Number of vertices in <em>child_conn</em>
  //! \param child_dim Dimension of child entity being queried
  //! \param side_number Side number of child entity (returned)
  //! \param sense Sense of child entity with respect to order in <em>child_conn</em> (returned)
  //! \param offset Offset of <em>child_conn</em> with respect to canonical ordering data (returned)
  //! \return status Returns zero if successful, -1 if not
  static int SideNumber(const void *parent_conn, const MBEntityType parent_type,
                        const void *child_conn, const int child_num_verts,
                        const int child_dim,
                        int &side_number, int &sense, int &offset);

  //! given two connectivity arrays, determine whether or not they represent the same entity.
  //! \param conn1 Connectivity array of first entity
  //! \param conn2 Connectivity array of second entity
  //! \param num_vertices Number of entries in <em>conn1</em> and <em>conn2</em>
  //! \param direct If positive, entities have the same sense (returned)
  //! \param offset Offset of <em>conn2</em>'s first vertex in <em>conn1</em>
  //! \return bool Returns true if <em>conn1</em> and <em>conn2</em> match
  static bool ConnectivityMatch(const void *conn1,
                                const void *conn2,
                                const int num_vertices,
                                int &direct, int &offset);

  //! true if entities of a given type and number of nodes indicates mid edge nodes are present.
  //! \param this_type Type of entity for which sub-entity connectivity is being queried
  //! \param num_verts Number of nodes defining entity
  //! \return bool Returns true if <em>this_type</em> combined with <em>num_nodes</em> indicates
  //!  mid-edge nodes are likely
  static bool HasMidEdgeNodes(const MBEntityType this_type, 
                              const int num_verts);

  //! true if entities of a given type and number of nodes indicates mid face nodes are present.
  //! \param this_type Type of entity for which sub-entity connectivity is being queried
  //! \param num_verts Number of nodes defining entity
  //! \return bool Returns true if <em>this_type</em> combined with <em>num_nodes</em> indicates
  //!  mid-face nodes are likely
  static bool HasMidFaceNodes(const MBEntityType this_type, 
                              const int num_verts);

  //! true if entities of a given type and number of nodes indicates mid region nodes are present.
  //! \param this_type Type of entity for which sub-entity connectivity is being queried
  //! \param num_verts Number of nodes defining entity
  //! \return bool Returns true if <em>this_type</em> combined with <em>num_nodes</em> indicates
  //!  mid-region nodes are likely
  static bool HasMidRegionNodes(const MBEntityType this_type, 
                                const int num_verts);

  //! true if entities of a given type and number of nodes indicates mid edge/face/region nodes 
  //! are present.
  //! \param this_type Type of entity for which sub-entity connectivity is being queried
  //! \param num_verts Number of nodes defining entity
  //! \param mid_nodes If <em>mid_nodes[i], i=0..2</em> is true, indicates that mid-edge 
  //!    (i=0), mid-face (i=1), and/or mid-region (i=2) nodes are likely
  static void HasMidNodes(const MBEntityType this_type, 
                          const int num_verts, 
                          bool mid_nodes[3]);

  //! given data about an element and a vertex in that element, return the dimension
  //! and index of the sub-entity that the vertex resolves.  If it does not resolve a
  //! sub-entity, either because it's a corner node or it's not in the element, -1 is
  //! returned in both return values
  //! \param elem_conn Connectivity of the entity being queried
  //! \param elem_type Type of entity being queried
  //! \param num_verts Number of vertices in <em>elem_conn</em>
  //! \param ho_node Handle of high-order node being queried
  //! \param parent_dim Dimension of sub-entity high-order node resolves (returned)
  //! \param parent_index Index of sub-entity high-order node resolves (returned)
  static void HONodeParent(const void *elem_conn, 
                           const MBEntityType elem_type,
                           const int num_verts, 
                           const void *ho_node,
                           int &parent_dim, 
                           int &parent_index);

  //! for an entity of this type with num_verts vertices, and a specified subfacet 
  //! (dimension and index), return the index of the higher order node for that entity 
  //! in this entity's connectivity array
  //! \param this_type Type of entity being queried
  //! \param num_verts Number of vertices for the entity being queried
  //! \param subfacet_dim Dimension of sub-entity being queried
  //! \param subfacet_index Index of sub-entity being queried
  //! \return index Index of sub-entity's higher-order node
  static int HONodeIndex(const MBEntityType this_type, const int num_verts,
                         const int subfacet_dim, const int subfacet_index);

    //! 
};

  //! get the basis of the numbering system
inline int MBCN::GetBasis() {return numberBasis;}
  
inline const char *MBCN::EntityTypeName(const MBEntityType this_type) 
{
  return entityTypeNames[this_type];
}

inline int MBCN::Dimension(const MBEntityType t) 
{
  return mConnectivityMap[t][0].topo_dimension;
}

inline int MBCN::VerticesPerEntity(const MBEntityType t) 
{
  return (MBVERTEX == t ? 1 : mConnectivityMap[t][mConnectivityMap[t][0].topo_dimension-1].num_nodes_per_sub_element[0]);
}

inline int MBCN::NumSubEntities(const MBEntityType t, const int d)
{
  return (t != MBVERTEX && d > 0 ? mConnectivityMap[t][d-1].num_sub_elements :
          (d ? -1 : VerticesPerEntity(t)));
}

  //! return the type of a particular sub-entity.
inline MBEntityType MBCN::SubEntityType(const MBEntityType this_type,
                                        const int sub_dimension,
                                        const int index) 
{
  
  return (!sub_dimension ? MBVERTEX : 
          (Dimension(this_type) == sub_dimension && 0 == index ? this_type :
          mConnectivityMap[this_type][sub_dimension-1].target_type[index]));
}
  
  //! return the connectivity of the specified sub-entity.
inline void MBCN::SubEntityVertexIndices(const MBEntityType this_type, 
                                         const int sub_dimension,
                                         const int index,
                                         int sub_entity_conn[]) 
{
  for (int i = 0; i < VerticesPerEntity(SubEntityType(this_type, sub_dimension, index)); i++)
    sub_entity_conn[i] = (MBVERTEX == this_type && 0 == sub_dimension && 
                          0 == index) ? 0 : 
      (0 == sub_dimension ? index :
      mConnectivityMap[this_type][sub_dimension-1].conn[index][i]);
}

inline bool MBCN::HasMidEdgeNodes(const MBEntityType this_type, 
                                     const int num_nodes)
{
    // poly elements never have mid nodes as far as canonical ordering is concerned
  if (MBPOLYGON == this_type || MBPOLYHEDRON == this_type) return false;

  if (num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1)) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 
                    NumSubEntities(this_type, 2)) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 
                    NumSubEntities(this_type, 2) + 1) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 1) )
    return true;
  
  else
    return false;
}

inline bool MBCN::HasMidFaceNodes(const MBEntityType this_type, 
                                       const int num_nodes)
{
    // poly elements never have mid nodes as far as canonical ordering is concerned
  if (MBPOLYGON == this_type || MBPOLYHEDRON == this_type) return false;

  if (num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 2)) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 
                    NumSubEntities(this_type, 2)) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 
                    NumSubEntities(this_type, 2) + 1) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 2) + 1) )
    return true;
  
  else
    return false;
}

inline bool MBCN::HasMidRegionNodes(const MBEntityType this_type, 
                                         const int num_nodes)
{
    // poly elements never have mid nodes as far as canonical ordering is concerned
  if (MBPOLYGON == this_type || MBPOLYHEDRON == this_type) return false;

  if (num_nodes == (VerticesPerEntity(this_type) + 1) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 1) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 1) + 
                    NumSubEntities(this_type, 2) + 1) ||
      num_nodes == (VerticesPerEntity(this_type) + NumSubEntities(this_type, 2) + 1) )
    return true;
  
  else
    return false;
}

inline void MBCN::HasMidNodes(const MBEntityType this_type, const int num_nodes,
                                  bool mid_nodes[3])
{
  if (num_nodes <= VerticesPerEntity(this_type) ||
    // poly elements never have mid nodes as far as canonical ordering is concerned
      MBPOLYGON == this_type || MBPOLYHEDRON == this_type) {
    mid_nodes[0] = false;
    mid_nodes[1] = false;
    mid_nodes[2] = false;
    return;
  }
  
  mid_nodes[0] = HasMidEdgeNodes(this_type, num_nodes);
  mid_nodes[1] = HasMidFaceNodes(this_type, num_nodes);
  mid_nodes[2] = HasMidRegionNodes(this_type, num_nodes);
}


#endif
