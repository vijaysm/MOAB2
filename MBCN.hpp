//
// MBCN: generic connectivity tools for meshes
//
#ifndef MBCN_HPP
#define MBCN_HPP

#include "assert.h"
#include <vector>
#include <algorithm>

// the maximum number n-1 dimension adjacencies a element may have
#define MB_MAX_SUB_ENTITIES  12

// the maximum number of nodes an n-1 dimensional element may have
#define MB_MAX_SUB_ENTITY_VERTICES 8


/*! The ordering here must ensure that all element types are 
 *  grouped together and all elements of similar dimension are
 *  grouped together.
 */
enum MBEntityType
{
  MBVERTEX = 0, //!< Mesh Vertex AKA node
  MBEDGE,       //!< Mesh Edge
  MBTRI,        //!< Triangular element (including shells)
  MBQUAD,       //!< Quadrilateral element (including shells)
  MBTET,        //!< Tetrahedral element
  MBPYRAMID,    //!< Pyramid element (where are the face ids for this defined?)
  MBPRISM,      //!< Wedge element (Exodus has one, Cubit doesn't. Does Mesh need it?)
  MBKNIFE,      //!< Knife element
  MBHEX,        //!< Hexahedral element
  MBENTITYSET,    //!< MeshSet
  MBMAXTYPE  //!< Just a place keeper - must be the # of entities, for array
    //!< dimensioning purposes 
};

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

//! switch the basis
  static void SwitchBasis(const int old_basis, const int new_basis);
  
public:

    //! enum used to specify operation type
  enum {INTERSECT, UNION};

    //! each entity type has two ConnMap objects, holding information about the bounding
    //! edges and faces for each entity; see comment for mConnectivityMap
  struct ConnMap
  {
    int topo_dimension;
    int num_sub_elements;
    int num_nodes_per_sub_element[MB_MAX_SUB_ENTITIES];
    MBEntityType target_type[MB_MAX_SUB_ENTITIES];
    int conn[MB_MAX_SUB_ENTITIES][MB_MAX_SUB_ENTITY_VERTICES];
  };

    //! mConnectivityMap[i=entity type][j=0,1,2]:
    //!  num_sub_elements = # bounding edges(j=0) or faces(j=1) for entity type i, or self (j=2)
    //!  num_nodes_per_sub_element[k] (k=0..num_sub_elements-1) = number of nodes in sub-facet k
    //!    (can vary over sub-facets, e.g. faces bounding a pyramid) or self (j=2)
    //!  target_type[k] = entity type of sub-facet k (e.g. MBTRI or MBQUAD bounding a pyramid) or self (j=2)
    //!  conn[k][l] (l=0..MBCN::VerticesPerEntity[target_type[k]]) = vertex connectivity of sub-facet k,
    //!    with respect to entity i's canonical vertex ordering, or self (j=2)
  static const ConnMap mConnectivityMap[MBMAXTYPE][3];

  struct UpConnMap
  {
    int num_targets_per_source_element[MB_MAX_SUB_ENTITIES];
    int targets_per_source_element[MB_MAX_SUB_ENTITIES][MB_MAX_SUB_ENTITIES];
  };

  static const UpConnMap mUpConnMap[MBMAXTYPE][4][4];

  //! this const vector defines the starting and ending MBEntityType for 
  //! each dimension,
  //! i.e. TypeDimensionMap[2] returns a pair of MBEntityTypes bounding dimension 2.
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
  static MBEntityType SubEntityType(const MBEntityType this_type,
                                         const int sub_dimension,
                                         const int index);
  
  //! return the connectivity of the specified sub-entity.
  static void SubEntityConn(const MBEntityType this_type, 
                            const int sub_dimension,
                            const int index,
                            int sub_entity_conn[]);

  //! For a specified set of sides of given dimension, return the intersection 
  //! or union of all sides of specified target dimension adjacent to those sides.
  static int AdjacentSubEntities(const MBEntityType this_type,
                                 const int *source_indices,
                                 const int num_source_indices,
                                 const int source_dim,
                                 const int target_dim,
                                 std::vector<int> &index_list,
                                 const int operation_type = MBCN::INTERSECT);

  //! return the side index represented in the input sub-entity connectivity in the input 
  //! parent entity connectivity array.
  static int SideNumber(const void *parent_conn, const MBEntityType parent_type,
                        const void *child_conn, const int child_num_verts,
                        const int child_dim,
                        int &side_number, int &sense, int &offset);

  //! given two connectivity arrays, determine whether or not they represent the same entity.
  static bool ConnectivityMatch(const void *conn1,
                                const void *conn2,
                                const int num_vertices,
                                int &direct, int &offset);

  //! true if entities of a given type and number of nodes indicates mid edge nodes are present.
  static bool HasMidEdgeNodes(const MBEntityType this_type, const int num_nodes);

  //! true if entities of a given type and number of nodes indicates mid face nodes are present.
  static bool HasMidFaceNodes(const MBEntityType this_type, const int num_verts);

  //! true if entities of a given type and number of nodes indicates mid region nodes are present.
  static bool HasMidRegionNodes(const MBEntityType this_type, const int num_verts);

  //! true if entities of a given type and number of nodes indicates mid edge/face/region nodes 
  //! are present.
  static void HasMidNodes(const MBEntityType this_type, const int num_nodes, 
                          bool mid_nodes[3]);

  //! given data about an element and a vertex in that element, return the dimension
  //! and index of the sub-entity that the vertex resolves.  If it does not resolve a
  //! sub-entity, either because it's a corner node or it's not in the element, -1 is
  //! returned in both return values
  static void HONodeParent(const void *elem_conn, const MBEntityType elem_type,
                           const int num_verts, const void *ho_node,
                           int &parent_dim, int &parent_index);

  //! for an entity of this type with num_verts vertices, and a specified subfacet 
  //! (dimension and index), return the index of the higher order node for that entity 
  //! in this entity's connectivity array
  static int HONodeIndex(const MBEntityType this_type, const int num_verts,
                         const int subfacet_dim, const int subfacet_index);
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
  return mConnectivityMap[t][mConnectivityMap[t][0].topo_dimension-1].num_nodes_per_sub_element[0];
}

inline int MBCN::NumSubEntities(const MBEntityType t, const int d)
{
  return mConnectivityMap[t][d-1].num_sub_elements;
}

  //! return the type of a particular sub-entity.
inline MBEntityType MBCN::SubEntityType(const MBEntityType this_type,
                                               const int sub_dimension,
                                               const int index) 
{
  return mConnectivityMap[this_type][sub_dimension-1].target_type[index];
}
  
  //! return the connectivity of the specified sub-entity.
inline void MBCN::SubEntityConn(const MBEntityType this_type, 
                                  const int sub_dimension,
                                  const int index,
                                  int sub_entity_conn[]) 
{
  for (int i = 0; i < VerticesPerEntity(SubEntityType(this_type, sub_dimension, index)); i++)
    sub_entity_conn[i] = mConnectivityMap[this_type][sub_dimension-1].conn[index][i];
}

inline bool MBCN::HasMidEdgeNodes(const MBEntityType this_type, 
                                     const int num_nodes)
{
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
  if (num_nodes <= VerticesPerEntity(this_type)) {
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
