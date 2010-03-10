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



#ifndef MB_SKINNER_HPP
#define MB_SKINNER_HPP

#include "MBForward.hpp"
#include <vector>

class MBSkinner 
{

  enum direction{FORWARD=1, REVERSE=-1};
protected:
  //! the MB instance that this works with
  MBInterface* thisMB;

  MBTag mDeletableMBTag;
  MBTag mAdjTag;
  int mTargetDim;

public:
  //! constructor, takes mdb instance
  MBSkinner(MBInterface* mdb) 
    : thisMB(mdb), mDeletableMBTag(0), mAdjTag(0){}

  //! destructor
  ~MBSkinner();

  MBErrorCode find_geometric_skin(MBRange &forward_target_entities);
  
  // will accept entities all of one dimension
  // and return entities of n-1 dimension
  MBErrorCode find_skin( const MBRange &entities,
                         bool get_vertices,
                         MBRange &output_handles,
                         MBRange *output_reverse_handles = 0,
                         bool create_vert_elem_adjs = false,
                         bool create_skin_elements = true);

    // get skin entities of prescribed dimension
  MBErrorCode find_skin(const MBRange &entities,
                        int dim,
                        MBRange &skin_entities,
                        bool create_vert_elem_adjs = false);

    /**\brief Find vertices on the skin of a set of mesh entities.
     *\param entities The elements for which to find the skin.  Range
     *                may NOT contain vertices, polyhedra, or entity sets.
     *                All elements in range must be of the same dimension.
     *\param skin_verts Output: the vertices on the skin.
     *\param skin_elems Optional output: elements representing sides of entities 
     *                    that are on the skin
     *\param create_if_missing If skin_elemts is non-null and this is true, 
     *                    create new elements representing the sides of 
     *                    entities on the skin.  If this is false, skin_elems
     *                    will contain only those skin elements that already
     *                    exist.
     */
  MBErrorCode find_skin_vertices( const MBRange& entities,
                                  MBRange& skin_verts,
                                  MBRange* skin_elems = 0,
                                  bool create_if_missing = true );

  MBErrorCode classify_2d_boundary( const MBRange &boundary,
                                     const MBRange &bar_elements,
                                     MBEntityHandle boundary_edges,
                                     MBEntityHandle inferred_edges,
                                     MBEntityHandle non_manifold_edges,
                                     MBEntityHandle other_edges,
                                     int &number_boundary_nodes);
  
  //!given a skin of dimension 2, will classify and return edges
  //! as boundary, inferred, and non-manifold, and the rest (other)
  MBErrorCode classify_2d_boundary( const MBRange  &boundary,
                                     const MBRange  &mesh_1d_elements,
                                     MBRange  &boundary_edges,
                                     MBRange  &inferred_edges,
                                     MBRange  &non_manifold_edges,
                                     MBRange  &other_edges,
                                     int &number_boundary_nodes);

protected:
  
  void initialize();
  
  void deinitialize();

  MBErrorCode find_skin_noadj( const MBRange &source_entities,
                               MBRange &forward_target_entities,
                               MBRange &reverse_target_entities );

  void add_adjacency(MBEntityHandle entity);
  
  void add_adjacency(MBEntityHandle entity, const MBEntityHandle *conn,
                     const int num_nodes);

  MBErrorCode remove_adjacency(MBEntityHandle entity);

  bool entity_deletable(MBEntityHandle entity);

  void find_match( MBEntityType type, 
                   const MBEntityHandle *conn, 
                   const int num_nodes,
                   MBEntityHandle& match,
                   MBSkinner::direction &direct);

  bool connectivity_match(const MBEntityHandle *conn1,
                          const MBEntityHandle *conn2,
                          const int num_verts,
                          MBSkinner::direction &direct);

  void find_inferred_edges(MBRange &skin_boundary,
                           MBRange &candidate_edges,
                           MBRange &inferred_edges,
                           double reference_angle_degrees);

  bool has_larger_angle(MBEntityHandle &entity1,
                       MBEntityHandle &entity2,
                       double reference_angle_cosine);


    /**\brief Find vertices on the skin of a set of mesh entities.
     *\param entities The elements for which to find the skin.  Range
     *                may NOT contain vertices, polyhedra, or entity sets.
     *                All elements in range must be of the same dimension.
     *\param skin_verts Output: the vertices on the skin.
     *\param skin_elems Optional output: elements representing sides of entities 
     *                    that are on the skin
     *\param create_if_missing If skin_elemts is non-null and this is true, 
     *                    create new elements representing the sides of 
     *                    entities on the skin.  If this is false, skin_elems
     *                    will contain only those skin elements that already
     *                    exist.
     */
  MBErrorCode find_skin_vertices( const MBRange& entities,
                                  MBRange* skin_verts = 0,
                                  MBRange* skin_elems = 0,
                                  MBRange* rev_elems = 0,
                                  bool create_if_missing = true,
                                  bool corners_only = false );

  /**\brief Skin edges
   *
   * Return any vertices adjacent to exactly one of the input edges.
   */
  MBErrorCode find_skin_vertices_1D( MBTag tag,
                                     const MBRange& edges,
                                     MBRange& skin_verts );
                                     
  /**\brief Skin faces
   *
   * For the set of face sides (logical edges), return 
   * vertices on such sides and/or edges equivalent to such sides.
   *\param faces  Set of toplogically 2D entities to skin.
   *\param skin_verts If non-NULL, skin vertices will be added to this container.
   *\param skin_edges If non-NULL, skin edges will be added to this container
   *\param reverse_edges If skin_edges is not NULL and this is not NULL, then
   *                  any existing skin edges that are reversed with respect
   *                  to the skin side will be placed in this range instead of
   *                  skin_edges.  Note: this argument is ignored if skin_edges
   *                  is NULL.
   *\param create_edges If true, edges equivalent to face sides on the skin
   *                  that don't already exist will be created.  Note: this
   *                  parameter is honored regardless of whether or not skin
   *                  edges or vertices are returned.
   *\param corners_only If true, only skin vertices that correspond to the
   *                  corners of sides will be returned (i.e. no higher-order
   *                  nodes.)  This argument is ignored if skin_verts is NULL.
   */
  MBErrorCode find_skin_vertices_2D( MBTag tag,
                                     const MBRange& faces,
                                     MBRange* skin_verts = 0,
                                     MBRange* skin_edges = 0,
                                     MBRange* reverse_edges = 0,
                                     bool create_edges = false,
                                     bool corners_only = false );
                                     
  /**\brief Skin volume mesh
   *
   * For the set of element sides (logical faces), return 
   * vertices on such sides and/or faces equivalent to such sides.
   *\param entities  Set of toplogically 3D entities to skin.
   *\param skin_verts If non-NULL, skin vertices will be added to this container.
   *\param skin_faces If non-NULL, skin faces will be added to this container
   *\param reverse_faces If skin_faces is not NULL and this is not NULL, then
   *                  any existing skin faces that are reversed with respect
   *                  to the skin side will be placed in this range instead of
   *                  skin_faces.  Note: this argument is ignored if skin_faces
   *                  is NULL.
   *\param create_faces If true, face equivalent to sides on the skin
   *                  that don't already exist will be created.  Note: this
   *                  parameter is honored regardless of whether or not skin
   *                  faces or vertices are returned.
   *\param corners_only If true, only skin vertices that correspond to the
   *                  corners of sides will be returned (i.e. no higher-order
   *                  nodes.)  This argument is ignored if skin_verts is NULL.
   */
  MBErrorCode find_skin_vertices_3D( MBTag tag,
                                     const MBRange& entities,
                                     MBRange* skin_verts = 0,
                                     MBRange* skin_faces = 0,
                                     MBRange* reverse_faces = 0,
                                     bool create_faces = false,
                                     bool corners_only = false );

  MBErrorCode create_side( MBEntityHandle element,
                           MBEntityType side_type,
                           const MBEntityHandle* side_corners,
                           MBEntityHandle& side_elem_handle_out );
                           
  bool edge_reversed( MBEntityHandle face, const MBEntityHandle edge_ends[2] );
  bool face_reversed( MBEntityHandle region, const MBEntityHandle* face_conn, 
                      MBEntityType face_type );
};


#endif

