/*
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

/**\file MBAdaptiveKDTree.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-04-1
 */

#ifndef MB_ADAPTIVE_KD_TREE_HPP
#define MB_ADAPTIVE_KD_TREE_HPP

#include "MBTypes.h"

#include <string>
#include <vector>
#include <math.h>

class MBAdaptiveKDTreeIter;
class MBInterface;
class MBRange;

class MBAdaptiveKDTree
{
private:

  MBInterface* mbInstance;
  MBTag planeTag, axisTag, rootTag;
  unsigned meshSetFlags;
  
public:

  MBAdaptiveKDTree( MBInterface* iface, 
                    const char* tagname = 0,
                    unsigned meshset_creation_flags = MESHSET_SET );

  //! Enumeriate split plane directions
  enum Axis { X = 0, Y = 1, Z = 2 };
  
  //! Split plane 
  struct Plane {
    double coord;  //!< Location of plane as coordinate on normal axis
    int norm; //!< The principal axis that is the normal of the plane;
    
      /** return true if point is below/to the left of the split plane */
    bool left_side( const double point[3] ) {
      return point[norm] < coord;
    }
      /** return true if point is abve/to the right of the split plane */
    bool right_side( const double point[3] ) {
      return point[norm] > coord;
    }
      /** return distance from point to plane */
    double distance( const double point[3] ) const {
      return fabs(point[norm] - coord);
    }
  };
  
  //! Get split plane for tree node
  MBErrorCode get_split_plane( MBEntityHandle node, Plane& plane );
  
  //! Set split plane for tree node
  MBErrorCode set_split_plane( MBEntityHandle node, const Plane& plane );
  
  //! Get bounding box for entire tree
  MBErrorCode get_tree_box( MBEntityHandle root_node,
                            double box_min_out[3], 
                            double box_max_out[3] );
  
  //! Set bounding box for entire tree
  MBErrorCode set_tree_box( MBEntityHandle root_node,
                            const double box_min[3], 
                            const double box_max[3] );
  
  //! Create tree root node
  MBErrorCode create_tree( const double box_min[3],
                           const double box_max[3],
                           MBEntityHandle& root_handle );

  //! Find all tree roots
  MBErrorCode find_all_trees( MBRange& results );

  //! Destroy a tree
  MBErrorCode delete_tree( MBEntityHandle root_handle );

  MBInterface* moab() { return mbInstance; }

  //! Get iterator for tree
  MBErrorCode get_tree_iterator( MBEntityHandle tree_root,
                                 MBAdaptiveKDTreeIter& result );
  
  //! Get iterator at right-most ('last') leaf.
  MBErrorCode get_last_iterator( MBEntityHandle tree_root,
                                 MBAdaptiveKDTreeIter& result );

  //! Get iterator for tree or subtree
  MBErrorCode get_sub_tree_iterator( MBEntityHandle tree_root,
                                     const double box_min[3], 
                                     const double box_max[3],
                                     MBAdaptiveKDTreeIter& result );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBAdaptiveKDTreeIter& leaf, Plane plane );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBAdaptiveKDTreeIter& leaf, 
                          Plane plane,
                          MBEntityHandle& left_child,
                          MBEntityHandle& right_child );
  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBAdaptiveKDTreeIter& leaf, 
                          Plane plane,
                          const MBRange& left_entities,
                          const MBRange& right_entities );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBAdaptiveKDTreeIter& leaf, 
                          Plane plane,
                          const std::vector<MBEntityHandle>& left_entities,
                          const std::vector<MBEntityHandle>& right_entities );
  
  //! Merge the leaf pointed to by the current iterator with it's
  //! sibling.  If the sibling is not a leaf, multiple merges may
  //! be done.
  MBErrorCode merge_leaf( MBAdaptiveKDTreeIter& iter );
  
    //! methods for selecting candidate split planes
  enum CandidatePlaneSet {
    //! Candidiate planes at evenly spaced intervals 
    SUBDIVISION,
    //! Like SUBDIVISION, except snap to closest vertex coordinate
    SUBDIVISION_SNAP,
    //! Median vertex coodinate values
    VERTEX_MEDIAN,
    //! Random sampling of vertex coordinate values
    VERTEX_SAMPLE,
  };
  
    //! Settings used for tree construction
  struct Settings {
    Settings(); //!< initialize to defaults
    unsigned maxEntPerLeaf; //!< split leafs with more entities than this
    unsigned maxTreeDepth;  //!< limit on the depth of the tree
    unsigned candidateSplitsPerDir; //!< number of cadiditate split planes to consider in each axial direction
    CandidatePlaneSet candidatePlaneSet;
    double minBoxWidth; //!< Tolerance
  };
  
  //! Build a tree
  MBErrorCode build_tree( MBRange& entities,
                          MBEntityHandle& root_set_out,
                          const Settings* settings = 0 );
  
  //! Find triangle closest to input position. 
  //!\param from_coords  The input position to test against
  //!\param closest_point_out  The closest point on the set of triangles in the tree
  //!\param triangle_out The triangle closest to the input position
  MBErrorCode closest_triangle( MBEntityHandle tree_root,
                                const double from_coords[3],
                                double closest_point_out[3],
                                MBEntityHandle& triangle_out );

  MBErrorCode sphere_intersect_triangles( MBEntityHandle tree_root,
                                          const double center[3],
                                          double radius,
                                          std::vector<MBEntityHandle>& triangles );

  MBErrorCode ray_intersect_triangles( MBEntityHandle tree_root,
                                       const double tolerance,
                                       const double ray_unit_dir[3],
                                       const double ray_base_pt[3],
                                       std::vector<MBEntityHandle>& triangles_out,
                                       std::vector<double>& distance_out,
                                       int result_count_limit = 0,
                                       double distance_limit = -1.0);

  //! Get leaf contianing input position.
  //!
  //! Does not take into account global bounding box of tree.
  //! - Therefore there is always one leaf containing the point.
  //! - If caller wants to account for global bounding box, then
  //!   caller can test against that box and not call this method
  //!   at all if the point is outside the box, as there is no leaf
  //!   containing the point in that case.
  MBErrorCode leaf_containing_point( MBEntityHandle tree_root,
                                     const double point[3],
                                     MBEntityHandle& leaf_out );

  //! Get iterator at leaf containing input position.
  //! 
  //! Returns MB_ENTITY_NOT_FOUND if point is not within
  //! bounding box of tree.
  MBErrorCode leaf_containing_point( MBEntityHandle tree_root,
                                     const double xyz[3],
                                     MBAdaptiveKDTreeIter& result );

  //! Find all leaves within a given distance from point in space.
  MBErrorCode leaves_within_distance( MBEntityHandle tree_root,
                                      const double from_point[3],
                                      const double distance,
                                      std::vector<MBEntityHandle>& leaves_out );

private:
  
  /**\brief find a triangle near the input point */
  MBErrorCode find_close_triangle( MBEntityHandle root,
                                   const double from_point[3],
                                   double pt[3],
                                   MBEntityHandle& triangle );
};
                    

//! Iterate over leaves of an adapative kD-tree
class MBAdaptiveKDTreeIter
{
public:

  enum Direction { LEFT = 0, RIGHT = 1 };

private:
  
  struct StackObj {
    StackObj( MBEntityHandle e, double c ) : entity(e), coord(c) {}
    StackObj() {}
    MBEntityHandle entity; //!< handle for tree node
    double coord;          //!< box coordinate of parent
  };
  
  enum { BMIN = 0, BMAX = 1 };  //!< indices into mBox and child list
  
  double mBox[2][3];                //!< min and max corners of bounding box
  MBAdaptiveKDTree* treeTool;       //!< tool for tree
  std::vector<StackObj> mStack;     //!< stack storing path through tree
  mutable std::vector<MBEntityHandle> childVect; //!< tempory storage of child handles
  
  //! Descend tree to left most leaf from current position
  //! No-op if at leaf.
  MBErrorCode step_to_first_leaf( Direction direction );

  friend class MBAdaptiveKDTree;
public:

  MBAdaptiveKDTreeIter() : treeTool(0), childVect(2) {}
  
  MBErrorCode initialize( MBAdaptiveKDTree* tool,
                          MBEntityHandle root,
                          const double box_min[3],
                          const double box_max[3],
                          Direction direction );

  MBAdaptiveKDTree* tool() const
    { return treeTool; }

    //! Get handle for current leaf
  MBEntityHandle handle() const
    { return mStack.back().entity; }
  
    //! Get min corner of axis-aligned box for current leaf
  const double* box_min() const 
    { return mBox[BMIN]; }
    
    //! Get max corner of axis-aligned box for current leaf
  const double* box_max() const 
    { return mBox[BMAX]; }
    
    //! Get depth in tree. root is at depth of 1.
  unsigned depth() const
    { return mStack.size(); }
  
  //! Advance the iterator either left or right in the tree
  //! Note:  stepping past the end of the tree will invalidate
  //!        the iterator.  It will *not* be work step the
  //!        other direction.
  MBErrorCode step( Direction direction );

    //! Advance to next leaf
    //! Returns MB_ENTITY_NOT_FOUND if at end.
    //! Note: steping past the end of the tree will invalidate
    //!       the iterator. Calling back() will not work.
  MBErrorCode step() { return step(RIGHT); }

    //! Move back to previous leaf
    //! Returns MB_ENTITY_NOT_FOUND if at beginning.
    //! Note: steping past the start of the tree will invalidate
    //!       the iterator. Calling step() will not work.
  MBErrorCode back() { return step(LEFT); }
  
  
    //! Return the side of the box bounding this tree node
    //! that is shared with the immediately adjacent sibling
    //! (the tree node that shares a common parent node with
    //! this node in the binary tree.)
    //!
    //!\param axis_out The principal axis orthogonal to the side of the box
    //!\param neg_out  true if the side of the box is toward the decreasing
    //!                direction of the principal axis indicated by axis_out,
    //!                false if it is toward the increasing direction.
    //!\return MB_ENTITY_NOT FOUND if root node.
    //!        MB_FAILURE if internal error.
    //!        MB_SUCCESS otherwise.
  MBErrorCode sibling_side( MBAdaptiveKDTree::Axis& axis_out, bool& neg_out ) const;

    //! Get adjacent leaf nodes on side indicated by norm and neg.
    //!
    //! E.g. if norm == X and neg == true, then get neighbor(s)
    //! adjacent to the side of the box contained in the plane
    //! with normal to the X axis and with the x coordinate equal 
    //! to the minimum x of the bounding box.
    //!
    //! E.g. if norm == Y and neg == false, then get neighbor(s)
    //! adjacent to the side of the box with y = maximum y of bounding box.
    //!
    //!\param norm  Normal vector for box side (X, Y, or Z)
    //!\param neg   Which of two planes with norm (true->smaller coord, 
    //!             false->larget coord)
    //!\param results List to which to append results.  This function does
    //!             *not* clear existing values in list.
    //!\param epsilon Tolerance on overlap.  A positive value E will
    //!              result in nodes that are separated by as much as E
    //!              to be considered touching.  A negative value -E will
    //!              cause leaves that do not overlap by at least E to be
    //!              considered non-overlapping.  Amongst other things, 
    //!              this value can be used to control whether or not
    //!              leaves adjacent at only their edges or corners are
    //!              returned.
  MBErrorCode get_neighbors( MBAdaptiveKDTree::Axis norm, bool neg,
                             std::vector<MBAdaptiveKDTreeIter>& results,
                             double epsilon = 0.0 ) const;
  
    //! Get split plane that separates this node from its immediate sibling.
  MBErrorCode get_parent_split_plane( MBAdaptiveKDTree::Plane& plane ) const;
};

#endif
