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
  };
  
  //! Get split plane for tree node
  MBErrorCode get_split_plane( MBEntityHandle entity, Plane& plane );
  
  //! Set split plane for tree node
  MBErrorCode set_split_plane( MBEntityHandle entity, const Plane& plane );
  
  //! Get bounding box for entire tree
  MBErrorCode get_tree_box( MBEntityHandle root_node,
                            double box_min_out[3], 
                            double box_max_out[3] );
  
  //! Set bounding box for entire tree
  MBErrorCode set_tree_box( MBEntityHandle root_node,
                            const double box_min_out[3], 
                            const double box_max_out[3] );
  
  //! Create tree root node
  MBErrorCode create_tree( const double box_min[3],
                           const double box_max[3],
                           MBEntityHandle& root_handle );

  //! Destroy a tree
  MBErrorCode delete_tree( MBEntityHandle root_handle );

  MBInterface* moab() { return mbInstance; }

  //! Get iterator for tree
  MBErrorCode get_tree_iterator( MBEntityHandle tree_root,
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
  
  
  struct Settings {
    Settings(); // initialize to defaults
    unsigned maxEntPerLeaf; //! split leafs with more entities than this
    unsigned maxTreeDepth;  //! limit on the depth of the tree
  };
  
  //! Build tree from triangles using simple bisection of nodes.
  //! Builds tree quicker, but results in less efficient tree.
  MBErrorCode build_tree_bisect_triangles( MBRange& triangles, 
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

  //! Find all leaves within a given distance from point in space.
  MBErrorCode leaves_within_distance( MBEntityHandle tree_root,
                                      const double from_point[3],
                                      const double distance,
                                      std::vector<MBEntityHandle>& leaves_out );
};
                    

//! Iterate over leaves of an adapative kD-tree
class MBAdaptiveKDTreeIter
{
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
  std::vector<MBEntityHandle> childVect; //!< tempory storage of child handles
  
  //! Descend tree to left most leaf from current position
  //! No-op if at leaf.
  MBErrorCode step_to_first_leaf();

  friend class MBAdaptiveKDTree;
public:

  MBAdaptiveKDTreeIter() : treeTool(0), childVect(2) {}
  
  MBErrorCode initialize( MBAdaptiveKDTree* tool,
                          MBEntityHandle root,
                          const double box_min[3],
                          const double box_max[3] );

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

    //! Advance to next leaf
    //! Returns MB_ENTITY_NOT_FOUND if at end.
  MBErrorCode step();

    //! Get adjacent leaf nodes on side indicated by norm and neg.
    //!
    //! E.g. if norm == X and neg == true, then get neighbor(s)
    //! adjacent to the size of the box contained in the plane
    //! with normal to the X axis and with the x coordinate equal 
    //! to the minimum x of the bounding box.
    //!
    //! E.g. if norm == Y and neg == false, then get neighbor(s)
    //! adjacent to the side of the box with y = maximum y of bounding box.
    //!
    //! Results are appended to list.  This function does not clear any
    //! existing values in the 'results' vector.
  MBErrorCode get_neighbors( MBAdaptiveKDTree::Axis norm, bool neg,
                             std::vector<MBAdaptiveKDTreeIter>& results ) const;
};

#endif
