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
 *\date 2008-05-12
 */

#ifndef MB_BSP_TREE_HPP
#define MB_BSP_TREE_HPP

#include "MBTypes.h"
#include "MBInterface.hpp"

#include <math.h>
#include <string>
#include <vector>

class MBBSPTreeBoxIter;
class MBBSPTreeIter;
class MBRange;
class BSPTreePoly;

class MBBSPTree
{
private:
  MBInterface* mbInstance;
  MBTag planeTag, rootTag;
  unsigned meshSetFlags;
  bool cleanUpTrees;
  std::vector<MBEntityHandle> createdTrees;

  MBErrorCode init_tags( const char* tagname = 0 );

public:
  
  static double epsilon() { return 1e-6; }

  MBBSPTree( MBInterface* iface,
             const char* tagname = 0,
             unsigned meshset_creation_flags = MESHSET_SET );

  MBBSPTree( MBInterface* iface,
             bool destroy_created_trees,
             const char* tagname = 0,
             unsigned meshset_creation_flags = MESHSET_SET );
  
  ~MBBSPTree();

  //! Enumeriate split plane directions
  enum Axis { X = 0, Y = 1, Z = 2 };
  
  /**\brief struct to store a plane 
   *
   * If plane is defined as Ax+By+Cz+D=0, then
   * norm={A,B,C} and coeff=-D.
   */
  struct Plane {
    Plane() {}
    Plane( const double n[3], double d ) : coeff(d)
      { norm[0] = n[0]; norm[1] = n[1]; norm[2] = n[2]; }
      /** a x + b y + c z + d = 0 */
    Plane( double a, double b, double c, double d ) : coeff(d)
      { norm[0] = a; norm[1] = b; norm[2] = c; }
      /** Create Y = 1 plane by doing Plane( Y, 1.0 ); */
    Plane( Axis normal, double point_on_axis ) : coeff(-point_on_axis)
      { norm[0] = norm[1] = norm[2] = 0; norm[normal] = 1.0; }
  
    double norm[3]; //!< Unit normal of plane
    double coeff;   //!< norm[0]*x + norm[1]*y + norm[2]*z + coeff = 0
    
      /** Signed distance from point to plane: 
       *    absolute value is distance from point to plane
       *    positive if 'above' plane, negative if 'below'
       *  Note: assumes unit-length normal.
       */
    double signed_distance( const double point[3] ) const
      { return point[0]*norm[0]+point[1]*norm[1]+point[2]*norm[2] + coeff; }
    
      /** return true if point is below the plane */
    bool below( const double point[3] ) const
      { return signed_distance(point) <= 0.0; }
    
      /** return true if point is above the plane */
    bool above( const double point[3] ) const
      { return signed_distance(point) >= 0.0; }
      
    double distance( const double point[3] ) const
      { return fabs( signed_distance(point) ); }
  
      /** reverse plane normal */
    void flip()
      { norm[0] = -norm[0]; norm[1] = -norm[1]; norm[2] = -norm[2]; coeff = -coeff; }
  
    void set( const double normal[3], const double point[3] )
      { 
         const double dot = normal[0]*point[0] + normal[1]*point[1] + normal[2]*point[2];
        *this = Plane( normal, -dot ); 
      }
      
    void set( const double pt1[3], const double pt2[3], const double pt3[3] );
    
    void set( double i, double j, double k, double coeff )
      { *this = Plane( i, j, k, coeff ); }
    
      /** Create Y = 1 plane by doing set( Y, 1.0 ); */
    void set( Axis normal, double point_on_axis )
      { 
        coeff = -point_on_axis;
        norm[0] = norm[1] = norm[2] = 0; 
        norm[normal] = 1.0; 
      }
  };
  
  //! Get split plane for tree node
  MBErrorCode get_split_plane( MBEntityHandle node, Plane& plane )
    { return moab()->tag_get_data( planeTag, &node, 1, &plane ); }
  
  //! Set split plane for tree node
  MBErrorCode set_split_plane( MBEntityHandle node, const Plane& plane );
  
  //! Get bounding box for entire tree
  MBErrorCode get_tree_box( MBEntityHandle root_node,
                            double corner_coords[8][3] );
  
  //! Get bounding box for entire tree
  MBErrorCode get_tree_box( MBEntityHandle root_node,
                            double corner_coords[24] );
  
  //! Set bounding box for entire tree
  MBErrorCode set_tree_box( MBEntityHandle root_node,
                            const double box_min[3], 
                            const double box_max[3] );
  MBErrorCode set_tree_box( MBEntityHandle root_node,
                            const double corner_coords[8][3] );
  
  //! Create tree root node
  MBErrorCode create_tree( const double box_min[3],
                           const double box_max[3],
                           MBEntityHandle& root_handle );
  MBErrorCode create_tree( const double corner_coords[8][3],
                           MBEntityHandle& root_handle );
  
  //! Create tree root node
  MBErrorCode create_tree( MBEntityHandle& root_handle );

  //! Find all tree roots
  MBErrorCode find_all_trees( MBRange& results );

  //! Destroy a tree
  MBErrorCode delete_tree( MBEntityHandle root_handle );

  MBInterface* moab() { return mbInstance; }

  //! Get iterator for tree
  MBErrorCode get_tree_iterator( MBEntityHandle tree_root,
                                 MBBSPTreeIter& result );
  
  //! Get iterator at right-most ('last') leaf.
  MBErrorCode get_tree_end_iterator( MBEntityHandle tree_root,
                                     MBBSPTreeIter& result );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBBSPTreeIter& leaf, Plane plane );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBBSPTreeIter& leaf, 
                          Plane plane,
                          MBEntityHandle& left_child,
                          MBEntityHandle& right_child );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBBSPTreeIter& leaf, 
                          Plane plane,
                          const MBRange& left_entities,
                          const MBRange& right_entities );

  //! Split leaf of tree
  //! Updates iterator location to point to first new leaf node.
  MBErrorCode split_leaf( MBBSPTreeIter& leaf, 
                          Plane plane,
                          const std::vector<MBEntityHandle>& left_entities,
                          const std::vector<MBEntityHandle>& right_entities );
  
  //! Merge the leaf pointed to by the current iterator with it's
  //! sibling.  If the sibling is not a leaf, multiple merges may
  //! be done.
  MBErrorCode merge_leaf( MBBSPTreeIter& iter );

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
                                     MBBSPTreeIter& result );
};

class MBBSPTreeIter
{
public:
  
  enum Direction { LEFT = 0, RIGHT = 1 };

private:
  friend class MBBSPTree;

  MBBSPTree* treeTool;

protected:
  std::vector<MBEntityHandle> mStack;
  mutable std::vector<MBEntityHandle> childVect;

  virtual MBErrorCode step_to_first_leaf( Direction direction );

  virtual MBErrorCode up();
  virtual MBErrorCode down( const MBBSPTree::Plane& plane, Direction direction );
  
  virtual MBErrorCode initialize( MBBSPTree* tool,
                                  MBEntityHandle root,
                                  const double* point = 0 );
  
public:
  
  MBBSPTreeIter() : treeTool(0), childVect(2) {}
  virtual ~MBBSPTreeIter() {}
  

  MBBSPTree* tool() const
    { return treeTool; }

    //! Get handle for current leaf
  MBEntityHandle handle() const
    { return mStack.back(); }
    
    //! Get depth in tree. root is at depth of 1.
  unsigned depth() const
    { return mStack.size(); }
  
  //! Advance the iterator either left or right in the tree
  //! Note:  stepping past the end of the tree will invalidate
  //!        the iterator.  It will *not* be work step the
  //!        other direction.
  virtual MBErrorCode step( Direction direction );

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
  
  
    //! Get split plane that separates this node from its immediate sibling.
  MBErrorCode get_parent_split_plane( MBBSPTree::Plane& plane ) const;
  
    //! Get volume of leaf polyhedron
  virtual double volume() const;
  
    //! Find range of overlap between ray and leaf.
    //!
    //!\param ray_point Coordinates of start point of ray
    //!\param ray_vect  Directionion vector for ray such that
    //!                 the ray is defined by r(t) = ray_point + t * ray_vect
    //!                 for t > 0.
    //!\param t_enter   Output: if return value is true, this value
    //!                 is the parameter location along the ray at which
    //!                 the ray entered the leaf.  If return value is false,
    //!                 then this value is undefined.
    //!\param t_exit    Output: if return value is true, this value
    //!                 is the parameter location along the ray at which
    //!                 the ray exited the leaf.  If return value is false,
    //!                 then this value is undefined.
    //!\return true if ray intersects leaf, false otherwise.
  virtual bool intersect_ray( const double ray_point[3],
                              const double ray_vect[3],
                              double& t_enter, double& t_exit ) const;
  
    //! Return true if thos node and the passed node share the
    //! same immediate parent.
  bool is_sibling( const MBBSPTreeIter& other_leaf ) const;
  
    //! Return true if thos node and the passed node share the
    //! same immediate parent.
  bool is_sibling( MBEntityHandle other_leaf ) const;
  
    //! Returns true if calling step() will advance to the
    //! immediate sibling of the current node.  Returns false
    //! if current node is root or back() will move to the 
    //! immediate sibling.
  bool sibling_is_forward( ) const;
  
    //! Calculate the convex polyhedron bounding this leaf.
  virtual MBErrorCode calculate_polyhedron( BSPTreePoly& polyhedron_out ) const;
};

class MBBSPTreeBoxIter : public MBBSPTreeIter
{
  private:
    
    double leafCoords[8][3];
    struct Corners { double coords[4][3]; };
    std::vector<Corners> stackData;
  
  protected:

    virtual MBErrorCode step_to_first_leaf( Direction direction );
    
    virtual MBErrorCode up();
    virtual MBErrorCode down( const MBBSPTree::Plane& plane, Direction direction );
  
    virtual MBErrorCode initialize( MBBSPTree* tool,
                                    MBEntityHandle root,
                                    const double* point = 0 );
  public:
  

    //! Faces of a hex : corner bitmap
  enum SideBits { B0154 = 0x33, //!< Face defined by corners {0,1,5,4}: 1st Exodus side
                  B1265 = 0x66, //!< Face defined by corners {1,2,6,5}: 2nd Exodus side
                  B2376 = 0xCC, //!< Face defined by corners {2,3,7,6}: 3rd Exodus side
                  B3047 = 0x99, //!< Face defined by corners {3,0,4,7}: 4th Exodus side
                  B3210 = 0x0F, //!< Face defined by corners {3,2,1,0}: 5th Exodus side
                  B4567 = 0xF0  //!< Face defined by corners {4,5,6,7}: 6th Exodus side
                 };
  
  static SideBits side_above_plane( const double hex_coords[8][3],
                                    const MBBSPTree::Plane& plane );
  
  static SideBits side_on_plane( const double hex_coords[8][3],
                                 const MBBSPTree::Plane& plane );
  
  static SideBits opposite_face( const SideBits& bits ) 
    { return (SideBits)((~bits) & 0xFF); }
  
  static MBErrorCode face_corners( const SideBits face, 
                                   const double hex_corners[8][3],
                                   double face_corners_out[4][3] );
  
  //! Advance the iterator either left or right in the tree
  //! Note:  stepping past the end of the tree will invalidate
  //!        the iterator.  It will *not* work to subsequently 
  //!        step the other direction.
  virtual MBErrorCode step( Direction direction );

    //! Advance to next leaf
    //! Returns MB_ENTITY_NOT_FOUND if at end.
    //! Note: steping past the end of the tree will invalidate
    //!       the iterator. Calling back() will not work.
  MBErrorCode step() { return MBBSPTreeIter::step(); }

    //! Move back to previous leaf
    //! Returns MB_ENTITY_NOT_FOUND if at beginning.
    //! Note: steping past the start of the tree will invalidate
    //!       the iterator. Calling step() will not work.
  MBErrorCode back() { return MBBSPTreeIter::back(); }
  
    //! Get coordinates of box corners, in Exodus II hexahedral ordering.
  MBErrorCode get_box_corners( double coords[8][3] ) const;
  
    //! Get volume of leaf box
  double volume() const;
  
    //! test if a plane intersects the leaf box
  enum XSect { MISS = 0, SPLIT = 1, NONHEX = -1 };
  XSect splits( const MBBSPTree::Plane& plane ) const;
  
    //! test if a plane intersects the leaf box
  bool intersects( const MBBSPTree::Plane& plane ) const;
  
    //! Find range of overlap between ray and leaf.
    //!
    //!\param ray_point Coordinates of start point of ray
    //!\param ray_vect  Directionion vector for ray such that
    //!                 the ray is defined by r(t) = ray_point + t * ray_vect
    //!                 for t > 0.
    //!\param t_enter   Output: if return value is true, this value
    //!                 is the parameter location along the ray at which
    //!                 the ray entered the leaf.  If return value is false,
    //!                 then this value is undefined.
    //!\param t_exit    Output: if return value is true, this value
    //!                 is the parameter location along the ray at which
    //!                 the ray exited the leaf.  If return value is false,
    //!                 then this value is undefined.
    //!\return true if ray intersects leaf, false otherwise.
  bool intersect_ray( const double ray_point[3],
                      const double ray_vect[3],
                      double& t_enter, double& t_exit ) const;
  
    //! Return the side of the box bounding this tree node
    //! that is shared with the immediately adjacent sibling
    //! (the tree node that shares a common parent node with
    //! this node in the binary tree.)
    //!
    //!\return MB_ENTITY_NOT FOUND if root node.
    //!        MB_FAILURE if internal error.
    //!        MB_SUCCESS otherwise.
  MBErrorCode sibling_side( SideBits& side_out ) const;
  
    //! Get adjacent leaf nodes on indicated side
    //!
    //!\param side   Face of box for which to retrieve neighbors
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
  MBErrorCode get_neighbors( SideBits side,
                             std::vector<MBBSPTreeBoxIter>& results,
                             double epsilon = 0.0 ) const;
  
    //! Calculate the convex polyhedron bounding this leaf.
  MBErrorCode calculate_polyhedron( BSPTreePoly& polyhedron_out ) const;
};

#endif
