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

/**\file MBAdaptiveKDTree.cpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-04-1
 */

#include "MBAdaptiveKDTree.hpp"
#include "MBInterface.hpp"
#include "MBGeomUtil.hpp"
#include "MBRange.hpp"
#include "MBInternals.hpp"

#include <assert.h>
#include <algorithm>
#include <limits>

#ifdef _MSC_VER
#  include <float.h>
#  define finite(A) _finite(A)
#endif

MBAdaptiveKDTree::Settings::Settings()
  : maxEntPerLeaf(6), 
    maxTreeDepth(30),
    candidateSplitsPerDir(5),
    candidatePlaneSet(SUBDIVISION_SNAP),
    minBoxWidth( std::numeric_limits<double>::epsilon() )
  {}


#define MB_AD_KD_TREE_DEFAULT_TAG_NAME "AKDTree"

// If defined, use single tag for both axis and location of split plane
#define MB_AD_KD_TREE_USE_SINGLE_TAG 

// No effect if MB_AD_KD_TREE_USE_SINGLE_TAG is not defined.
// If defined, store plane axis as double so tag has consistent
// type (doubles for both location and axis).  If not defined,
// store struct Plane as opaque.
#define MB_AD_KD_TREE_USE_TWO_DOUBLE_TAG

#if defined(MB_AD_KD_TREE_USE_SINGLE_TAG) && defined(HDF5_FILE)
# include <H5Tpublic.h>
#endif

#define MAKE_TAG( NAME, STORAGE, TYPE, COUNT, HANDLE, DEFAULT ) \
  if (MB_SUCCESS != make_tag( moab(), \
                              (NAME), \
                              (STORAGE), \
                              (TYPE), \
                              (COUNT), \
                              (DEFAULT), \
                              (HANDLE), \
                              ctl )) { \
    planeTag = axisTag = rootTag = (MBTag)-1; \
    return; \
  }

static MBErrorCode make_tag( MBInterface* iface,
                             std::string name,
                             MBTagType storage, 
                             MBDataType type,
                             int count,
                             void* default_val,
                             MBTag& tag_handle,
                             std::vector<MBTag>& created_tags )
{
  int size;
  switch (type) {
    case MB_TYPE_DOUBLE:  size = sizeof(double);         break;
    case MB_TYPE_INTEGER: size = sizeof(int);            break;
    case MB_TYPE_OPAQUE:  size = 1;                      break;
    case MB_TYPE_HANDLE:  size = sizeof(MBEntityHandle); break;
    default: return MB_FAILURE;
  }
  size *= count;
  
  MBErrorCode rval = iface->tag_create( name.c_str(),
                                        size,
                                        storage,
                                        type,
                                        tag_handle,
                                        default_val,
                                        false );
  if (MB_SUCCESS == rval) 
    created_tags.push_back( tag_handle );
  else if (MB_ALREADY_ALLOCATED == rval)
    rval = iface->tag_create( name.c_str(),
                              size,
                              MB_TAG_DENSE,
                              type,
                              tag_handle,
                              default_val,
                              true );

  if (MB_SUCCESS != rval)
    while( !created_tags.empty() ) {
      iface->tag_delete( created_tags.back() );
      created_tags.pop_back();
    }
  
  return rval;
}

MBAdaptiveKDTree::MBAdaptiveKDTree( MBInterface* mb, const char* tagname_in, unsigned set_flags )
  : mbInstance(mb), meshSetFlags(set_flags)
{
  const char* tagname = tagname_in ? tagname_in : MB_AD_KD_TREE_DEFAULT_TAG_NAME;
  std::vector<MBTag> ctl;

#ifndef MB_AD_KD_TREE_USE_SINGLE_TAG
    // create two tags, one for axis direction and one for axis coordinate
  std::string n1(tagname), n2(tagname);
  n1 += "_coord";
  n2 += "_norm";
  MAKE_TAG( n1, MB_TAG_DENSE, MB_TYPE_DOUBLE, 1, planeTag, 0 )
  MAKE_TAG( n2, MB_TAG_DENSE, MB_TYPE_INT,    1, axisTag,  0 )

#elif defined(MB_AD_KD_TREE_USE_TWO_DOUBLE_TAG)
    // create tag to hold two doubles, one for location and one for axis
  MAKE_TAG( tagname, MB_TAG_DENSE, MB_TYPE_DOUBLE, 2, planeTag, 0 )
#else
    // create opaque tag to hold struct Plane
  MAKE_TAG( tagname, MB_TAG_DENSE, MB_TYPE_OPAQUE, sizeof(Plane), planeTag, 0 )

#ifdef HDF5_FILE  
    // create a mesh tag holding the HDF5 type for a struct Plane
  MBTag type_tag;
  std::string type_tag_name = "__hdf5_tag_type_";
  type_tag_name += tagname;
  MAKE_TAG( type_tag_name, MB_TAG_MESH), MB_TYPE_OPAQUE, sizeof(hid_t), type_tag, 0 )
    // create HDF5 type object describing struct Plane
  Plane p;
  hid_t handle = H5Tcreate( H5T_COMPOUND, sizeof(Plane) );
  H5Tinsert( handle, "coord", &(p.coord) - &p, H5T_NATIVE_DOUBLE );
  H5Tinsert( handle, "norm", &(p.axis) - &p, H5T_NATIVE_INT );
  mbInstance->tag_set_data( type_tag, 0, 0, &handle );
#endif
#endif

  std::string root_name(tagname);
  root_name += "_box";
  MAKE_TAG( root_name, MB_TAG_SPARSE, MB_TYPE_DOUBLE, 6, rootTag, 0 )
}

MBErrorCode MBAdaptiveKDTree::get_split_plane( MBEntityHandle entity,
                                               Plane& plane )
{
#ifndef MB_AD_KD_TREE_USE_SINGLE_TAG
  MBErrorCode r1, r2;
  r1 = moab()->tag_get_data( planeTag, &entity, 1, &plane.coord );
  r2 = moab()->tag_get_data( axisTag , &entity, 1, &plane.norm  );
  return MB_SUCCESS == r1 ? r2 : r1;
#elif defined(MB_AD_KD_TREE_USE_TWO_DOUBLE_TAG)
  double values[2];
  MBErrorCode rval = moab()->tag_get_data( planeTag, &entity, 1, values );
  plane.coord = values[0];
  plane.norm = (int)values[1];
  return rval;
#else
  return moab()->tag_get_data( planeTag, &entity, 1, &plane );
#endif
}

MBErrorCode MBAdaptiveKDTree::set_split_plane( MBEntityHandle entity, 
                                               const Plane& plane )
{
#ifndef MB_AD_KD_TREE_USE_SINGLE_TAG
  MBErrorCode r1, r2;
  r1 = moab()->tag_set_data( planeTag, &entity, 1, &plane.coord );
  r2 = moab()->tag_set_data( axisTag , &entity, 1, &plane.norm  );
  return MB_SUCCESS == r1 ? r2 : r1;
#elif defined(MB_AD_KD_TREE_USE_TWO_DOUBLE_TAG)
  double values[2] = { plane.coord, plane.norm };
  return moab()->tag_set_data( planeTag, &entity, 1, values );
#else
  return moab()->tag_set_data( planeTag, &entity, 1, &plane );
#endif
}


MBErrorCode MBAdaptiveKDTree::set_tree_box( MBEntityHandle root_handle,
                                            const double box_min[3],
                                            const double box_max[3] )
{
  const double box[6] = { box_min[0], box_min[1], box_min[2],
                          box_max[0], box_max[1], box_max[2] };
  return moab()->tag_set_data( rootTag, &root_handle, 1, box );
}

MBErrorCode MBAdaptiveKDTree::get_tree_box( MBEntityHandle root_handle,
                                            double box_min_out[3],
                                            double box_max_out[3] )
{
  double box[6];
  MBErrorCode rval = moab()->tag_get_data( rootTag, &root_handle, 1, box );
  box_min_out[0] = box[0]; box_min_out[1] = box[1]; box_min_out[2] = box[2];
  box_max_out[0] = box[3]; box_max_out[1] = box[4]; box_max_out[2] = box[5];
  return rval;
}


MBErrorCode MBAdaptiveKDTree::create_tree( const double box_min[3],
                                           const double box_max[3],
                                           MBEntityHandle& root_handle )
{
  MBErrorCode rval = moab()->create_meshset( meshSetFlags, root_handle );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = set_tree_box( root_handle, box_min, box_max );
  if (MB_SUCCESS != rval) {
    moab()->delete_entities( &root_handle, 1 );
    root_handle = 0;
    return rval;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTree::delete_tree( MBEntityHandle root_handle )
{
  MBErrorCode rval;
  
  std::vector<MBEntityHandle> children, dead_sets, current_sets;
  current_sets.push_back( root_handle );
  while (!current_sets.empty()) {
    MBEntityHandle set = current_sets.back();
    current_sets.pop_back();
    dead_sets.push_back( set );
    rval = moab()->get_child_meshsets( set, children );
    if (MB_SUCCESS != rval)
      return rval;
    std::copy( children.begin(), children.end(), std::back_inserter(current_sets) );
    children.clear();
  }
  
  rval = moab()->tag_delete_data( rootTag, &root_handle, 1 );
  if (MB_SUCCESS != rval)
    return rval;
  
  return moab()->delete_entities( &dead_sets[0], dead_sets.size() );
}

MBErrorCode MBAdaptiveKDTree::find_all_trees( MBRange& results )
{
  return moab()->get_entities_by_type_and_tag( 0, MBENTITYSET, 
                                               &rootTag, 0, 1,
                                               results );
}

MBErrorCode MBAdaptiveKDTree::get_tree_iterator( MBEntityHandle root,
                                                 MBAdaptiveKDTreeIter& iter )
{
  double box[6];
  MBErrorCode rval = moab()->tag_get_data( rootTag, &root, 1, box );
  if (MB_SUCCESS != rval)
    return rval;
  
  return get_sub_tree_iterator( root, box, box+3, iter );
}

MBErrorCode MBAdaptiveKDTree::get_last_iterator( MBEntityHandle root,
                                                 MBAdaptiveKDTreeIter& iter )
{
  double box[6];
  MBErrorCode rval = moab()->tag_get_data( rootTag, &root, 1, box );
  if (MB_SUCCESS != rval)
    return rval;
  
  return iter.initialize( this, root, box, box+3, MBAdaptiveKDTreeIter::RIGHT );
}

MBErrorCode MBAdaptiveKDTree::get_sub_tree_iterator( MBEntityHandle root,
                                                     const double min[3], 
                                                     const double max[3],
                                                     MBAdaptiveKDTreeIter& result ) 
{
  return result.initialize( this, root, min, max, MBAdaptiveKDTreeIter::LEFT );
}

MBErrorCode MBAdaptiveKDTree::split_leaf( MBAdaptiveKDTreeIter& leaf,
                                          Plane plane,
                                          MBEntityHandle& left,
                                          MBEntityHandle& right )
{
  MBErrorCode rval;
  
  rval = moab()->create_meshset( meshSetFlags, left );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = moab()->create_meshset( meshSetFlags, right );
  if (MB_SUCCESS != rval) {
    moab()->delete_entities( &left, 1 );
    return rval;
  }
  
  if (MB_SUCCESS != set_split_plane( leaf.handle(), plane ) ||
      MB_SUCCESS != moab()->add_child_meshset( leaf.handle(), left ) ||
      MB_SUCCESS != moab()->add_child_meshset( leaf.handle(), right) ||
      MB_SUCCESS != leaf.step_to_first_leaf(MBAdaptiveKDTreeIter::LEFT)) {
    MBEntityHandle children[] = { left, right };
    moab()->delete_entities( children, 2 );
    return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTree::split_leaf( MBAdaptiveKDTreeIter& leaf,
                                          Plane plane )
{
  MBEntityHandle left, right;
  return split_leaf( leaf, plane, left, right );
}

MBErrorCode MBAdaptiveKDTree::split_leaf( MBAdaptiveKDTreeIter& leaf, 
                                          Plane plane,
                                          const MBRange& left_entities,
                                          const MBRange& right_entities )
{
  MBEntityHandle left, right, parent = leaf.handle();
  MBErrorCode rval = split_leaf( leaf, plane, left, right );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (MB_SUCCESS == moab()->add_entities( left, left_entities ) &&
      MB_SUCCESS == moab()->add_entities(right,right_entities ) &&
      MB_SUCCESS == moab()->clear_meshset( &parent, 1 ))
    return MB_SUCCESS;
  
  moab()->remove_child_meshset( parent, left );
  moab()->remove_child_meshset( parent, right );
  MBEntityHandle children[] = { left, right };
  moab()->delete_entities( children, 2 );
  return MB_FAILURE;
}

MBErrorCode MBAdaptiveKDTree::split_leaf( MBAdaptiveKDTreeIter& leaf, 
                                          Plane plane,
                                          const std::vector<MBEntityHandle>& left_entities,
                                          const std::vector<MBEntityHandle>& right_entities )
{
  MBEntityHandle left, right, parent = leaf.handle();
  MBErrorCode rval = split_leaf( leaf, plane, left, right );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (MB_SUCCESS == moab()->add_entities( left, &left_entities[0], left_entities.size() ) &&
      MB_SUCCESS == moab()->add_entities(right,&right_entities[0],right_entities.size() ) &&
      MB_SUCCESS == moab()->clear_meshset( &parent, 1 ))
    return MB_SUCCESS;
  
  moab()->remove_child_meshset( parent, left );
  moab()->remove_child_meshset( parent, right );
  MBEntityHandle children[] = { left, right };
  moab()->delete_entities( children, 2 );
  return MB_FAILURE;
}

MBErrorCode MBAdaptiveKDTree::merge_leaf( MBAdaptiveKDTreeIter& iter )
{
  MBErrorCode rval;
  if (iter.depth() == 1) // at root
    return MB_FAILURE;
  
    // Move iter to parent
  
  MBAdaptiveKDTreeIter::StackObj node = iter.mStack.back();
  iter.mStack.pop_back();
  
  iter.childVect.clear();
  rval = moab()->get_child_meshsets( iter.mStack.back().entity, iter.childVect );
  if (MB_SUCCESS != rval)
    return rval;
  Plane plane;
  rval = get_split_plane( iter.mStack.back().entity, plane );
  if (MB_SUCCESS != rval)
    return rval;
  
  int child_idx = iter.childVect[0] == node.entity ? 0 : 1;
  assert(iter.childVect[child_idx] == node.entity);
  iter.mBox[1-child_idx][plane.norm] = node.coord;
  

    // Get all entities from children and put them in parent
  MBEntityHandle parent = iter.handle();
  moab()->remove_child_meshset( parent, iter.childVect[0] );
  moab()->remove_child_meshset( parent, iter.childVect[1] );
  std::vector<MBEntityHandle> stack( iter.childVect );
  
  MBRange range;
  while (!stack.empty()) {
    MBEntityHandle h = stack.back();
    stack.pop_back();
    range.clear();
    rval = moab()->get_entities_by_handle( h, range );
    if (MB_SUCCESS != rval)
      return rval;
    rval = moab()->add_entities( parent, range );
    if (MB_SUCCESS != rval)
      return rval;
    
    iter.childVect.clear();
    moab()->get_child_meshsets( h, iter.childVect );
    if (!iter.childVect.empty()) {
     moab()->remove_child_meshset( h, iter.childVect[0] );
     moab()->remove_child_meshset( h, iter.childVect[1] );
     stack.push_back( iter.childVect[0] );
     stack.push_back( iter.childVect[1] );
    }
  
    rval = moab()->delete_entities( &h, 1 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  return MB_SUCCESS;
}

  

MBErrorCode MBAdaptiveKDTreeIter::initialize( MBAdaptiveKDTree* tool,
                                              MBEntityHandle root,
                                              const double box_min[3],
                                              const double box_max[3],
                                              Direction direction )
{
  mStack.clear();
  treeTool = tool;
  mBox[BMIN][0] = box_min[0];
  mBox[BMIN][1] = box_min[1];
  mBox[BMIN][2] = box_min[2];
  mBox[BMAX][0] = box_max[0];
  mBox[BMAX][1] = box_max[1];
  mBox[BMAX][2] = box_max[2];
  mStack.push_back( StackObj(root,0) );
  return step_to_first_leaf( direction );
}

MBErrorCode MBAdaptiveKDTreeIter::step_to_first_leaf( Direction direction )
{
  MBErrorCode rval;
  MBAdaptiveKDTree::Plane plane;
  const Direction opposite = static_cast<Direction>(1-direction);
  
  for (;;) {
    childVect.clear();
    rval = treeTool->moab()->get_child_meshsets( mStack.back().entity, childVect );
    if (MB_SUCCESS != rval)
      return rval;
    if (childVect.empty()) // leaf
      break;
  
    rval = treeTool->get_split_plane( mStack.back().entity, plane );
    if (MB_SUCCESS != rval)
      return rval;
  
    mStack.push_back( StackObj(childVect[direction],mBox[opposite][plane.norm]) );
    mBox[opposite][plane.norm] = plane.coord;
  }
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTreeIter::step( Direction direction )
{
  StackObj node, parent;
  MBErrorCode rval;
  MBAdaptiveKDTree::Plane plane;
  const Direction opposite = static_cast<Direction>(1-direction);
  
    // If stack is empty, then either this iterator is uninitialized
    // or we reached the end of the iteration (and return 
    // MB_ENTITY_NOT_FOUND) already.
  if (mStack.empty())
    return MB_FAILURE;
    
    // Pop the current node from the stack.
    // The stack should then contain the parent of the current node.
    // If the stack is empty after this pop, then we've reached the end.
  node = mStack.back();
  mStack.pop_back();
  
  while(!mStack.empty()) {
      // Get data for parent entity
    parent = mStack.back();
    childVect.clear();
    rval = treeTool->moab()->get_child_meshsets( parent.entity, childVect );
    if (MB_SUCCESS != rval)
      return rval;
    rval = treeTool->get_split_plane( parent.entity, plane );
    if (MB_SUCCESS != rval)
      return rval;
    
      // If we're at the left child
    if (childVect[opposite] == node.entity) {
        // change from box of left child to box of parent
      mBox[direction][plane.norm] = node.coord;
        // push right child on stack
      node.entity = childVect[direction];
      node.coord = mBox[opposite][plane.norm];
      mStack.push_back( node );
        // change from box of parent to box of right child
      mBox[opposite][plane.norm] = plane.coord;
        // descend to left-most leaf of the right child
      return step_to_first_leaf(opposite);
    }
    
      // The current node is the right child of the parent,
      // continue up the tree.
    assert( childVect[direction] == node.entity );
    mBox[opposite][plane.norm] = node.coord;
    node = parent;
    mStack.pop_back();
  }
  
  return MB_ENTITY_NOT_FOUND;
}

MBErrorCode MBAdaptiveKDTreeIter::get_neighbors( 
                      MBAdaptiveKDTree::Axis norm, bool neg,
                      std::vector<MBAdaptiveKDTreeIter>& results,
                      double epsilon ) const
{
  StackObj node, parent;
  MBErrorCode rval;
  MBAdaptiveKDTree::Plane plane;
  int child_idx;
  
    // Find tree node at which the specified side of the box
    // for this node was created.
  MBAdaptiveKDTreeIter iter( *this ); // temporary iterator (don't modifiy *this)
  node = iter.mStack.back();
  iter.mStack.pop_back();
  for (;;) {
      // reached the root - original node was on boundary (no neighbors)
    if (iter.mStack.empty())
      return MB_SUCCESS;
    
      // get parent node data
    parent = iter.mStack.back();
    iter.childVect.clear();
    rval = treeTool->moab()->get_child_meshsets( parent.entity, iter.childVect );
    if (MB_SUCCESS != rval)
      return rval;
    rval = treeTool->get_split_plane( parent.entity, plane );
    if (MB_SUCCESS != rval)
      return rval;
    
    child_idx = iter.childVect[0] == node.entity ? 0 : 1;
    assert(iter.childVect[child_idx] == node.entity);
    
      // if we found the split plane for the desired side
      // push neighbor on stack and stop
    if (plane.norm == norm && (int)neg == child_idx) {
        // change from box of previous child to box of parent
      iter.mBox[1-child_idx][plane.norm] = node.coord;
        // push other child of parent onto stack
      node.entity = iter.childVect[1-child_idx];
      node.coord = iter.mBox[child_idx][plane.norm];
      iter.mStack.push_back( node );
        // change from parent box to box of new child
      iter.mBox[child_idx][plane.norm] = plane.coord;
      break;
    }
    
      // continue up the tree
    iter.mBox[1-child_idx][plane.norm] = node.coord;
    node = parent;
    iter.mStack.pop_back();
  }

    // now move down tree, searching for adjacent boxes
  std::vector<MBAdaptiveKDTreeIter> list;
    // loop over all potential paths to neighbors (until list is empty)
  for (;;) {
      // follow a single path to a leaf, append any other potential
      // paths to neighbors to 'list'
    node = iter.mStack.back();
    for (;;) { 
      iter.childVect.clear();
      rval = treeTool->moab()->get_child_meshsets( node.entity, iter.childVect );
      if (MB_SUCCESS != rval)
        return rval;
        
        // if leaf
      if (iter.childVect.empty()) {
        results.push_back( iter );
        break; 
      }
      
      rval = treeTool->get_split_plane( node.entity, plane );
      if (MB_SUCCESS != rval)
        return rval;
     
        // if split parallel to side
      if (plane.norm == norm) {
          // continue with whichever child is on the correct side of the split
        node.entity = iter.childVect[neg];
        node.coord = iter.mBox[1-neg][plane.norm];
        iter.mStack.push_back( node );
        iter.mBox[1-neg][plane.norm] = plane.coord;
      }
        // if left child is adjacent
      else if (this->mBox[BMIN][plane.norm] - plane.coord <= epsilon) {
          // if right child is also adjacent, add to list
        if (plane.coord - this->mBox[BMAX][plane.norm] <= epsilon) {
          list.push_back( iter );
          list.back().mStack.push_back( StackObj( iter.childVect[1], iter.mBox[BMIN][plane.norm] ) );
          list.back().mBox[BMIN][plane.norm] = plane.coord;
        }
          // continue with left child
        node.entity = iter.childVect[0];
        node.coord = iter.mBox[BMAX][plane.norm];
        iter.mStack.push_back( node );
        iter.mBox[BMAX][plane.norm] = plane.coord;
      }
        // right child is adjacent
      else {
          // if left child is not adjacent, right must be or something
          // is really messed up.
        assert(plane.coord - this->mBox[BMAX][plane.norm] <= epsilon);
           // continue with left child
        node.entity = iter.childVect[1];
        node.coord = iter.mBox[BMIN][plane.norm];
        iter.mStack.push_back( node );
        iter.mBox[BMIN][plane.norm] = plane.coord;
      }
    }
    
    if (list.empty())
      break;
    
    iter = list.back();
    list.pop_back();
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTreeIter::sibling_side( 
                            MBAdaptiveKDTree::Axis& axis_out,
                            bool& neg_out ) const
{
  if (mStack.size() < 2) // at tree root
    return MB_ENTITY_NOT_FOUND;
  
  MBEntityHandle parent = mStack[mStack.size()-2].entity;
  MBAdaptiveKDTree::Plane plane;
  MBErrorCode rval = tool()->get_split_plane( parent, plane );
  if (MB_SUCCESS != rval)
    return MB_FAILURE;
    
  childVect.clear();
  rval = tool()->moab()->get_child_meshsets( parent, childVect );
  if (MB_SUCCESS != rval || childVect.size() != 2)
    return MB_FAILURE;
  
  axis_out = static_cast<MBAdaptiveKDTree::Axis>(plane.norm);
  neg_out = (childVect[1] == handle());
  assert(childVect[neg_out] == handle());
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTreeIter::get_parent_split_plane( MBAdaptiveKDTree::Plane& plane ) const
{
  if (mStack.size() < 2) // at tree root
    return MB_ENTITY_NOT_FOUND;
  
  MBEntityHandle parent = mStack[mStack.size()-2].entity;
  return tool()->get_split_plane( parent, plane );
}


static MBErrorCode intersect_children_with_elems( MBInterface* moab,
                                        const MBRange& elems,
                                        MBAdaptiveKDTree::Plane plane,
                                        MBCartVect box_min,
                                        MBCartVect box_max,
                                        MBRange& left_tris,
                                        MBRange& right_tris,
                                        MBRange& both_tris,
                                        double& metric_value )
{
  left_tris.clear();
  right_tris.clear();
  both_tris.clear();
  MBCartVect coords[8];
  
    // get extents of boxes for left and right sides
  MBCartVect right_min( box_min ), left_max( box_max );
  right_min[plane.norm] = left_max[plane.norm] = plane.coord;
  right_min *= 0.5;
  left_max *= 0.5;
  box_min *= 0.5;
  box_max *= 0.5;
  const MBCartVect left_cen = left_max + box_min;
  const MBCartVect left_dim = left_max - box_min;
  const MBCartVect right_cen = box_max + right_min;
  const MBCartVect right_dim = box_max - right_min;
  
  
    // test each triangle
  MBErrorCode rval;
  int count;
  const MBEntityHandle* conn;
  for (MBRange::reverse_iterator i = elems.rbegin(); i != elems.rend(); ++i) {
    rval = moab->get_connectivity( *i, conn, count, true );
    if (MB_SUCCESS != rval) return rval;
    if (count > (int)(sizeof(coords)/sizeof(coords[0])))
      return MB_FAILURE;
    rval = moab->get_coords( &conn[0], count, coords[0].array() );
    if (MB_SUCCESS != rval) return rval;
    
    bool lo = MBGeomUtil::box_elem_overlap( coords, TYPE_FROM_HANDLE(*i), left_cen, left_dim );
    bool ro = MBGeomUtil::box_elem_overlap( coords, TYPE_FROM_HANDLE(*i),right_cen,right_dim );
    
      // didn't intersect either - tolerance issue
    while (!lo && !ro) {
        // calculate a good tolerance
      MBCartVect dim = box_max - box_min;
      double max_dim;
      if (dim[0] > dim[1] && dim[1] > dim[2])
        max_dim = dim[0];
      else if (dim[1] > dim[2])
        max_dim = dim[1];
      else
        max_dim = dim[2];
        // loop with increasing tolerance until we intersect something
      double tol = std::numeric_limits<double>::epsilon();
      while (!lo && !ro) {
        lo = MBGeomUtil::box_elem_overlap( coords, TYPE_FROM_HANDLE(*i), left_cen,tol*max_dim* left_dim );
        ro = MBGeomUtil::box_elem_overlap( coords, TYPE_FROM_HANDLE(*i),right_cen,tol*max_dim*right_dim );
        tol *= 10.0;
        if (tol > 1e-3)
          return MB_FAILURE;
      }
    }
    if (lo && ro)
      both_tris.insert( *i );
    else if (lo)
      left_tris.insert( *i );
    else //if (ro)
      right_tris.insert( *i );
  }
  
  MBCartVect box_dim = box_max - box_min;
  double area_left = left_dim[0]*left_dim[1] + left_dim[1]*left_dim[2] + left_dim[2]*left_dim[0];
  double area_right = right_dim[0]*right_dim[1] + right_dim[1]*right_dim[2] + right_dim[2]*right_dim[0];
  double area_both = box_dim[0]*box_dim[1] + box_dim[1]*box_dim[2] + box_dim[2]*box_dim[0];
  metric_value = (area_left * left_tris.size() + area_right * right_tris.size()) / area_both + both_tris.size();
  return MB_SUCCESS;
}

static MBErrorCode best_subdivision_plane( int num_planes,
                                           const MBAdaptiveKDTreeIter& iter,
                                           MBRange& best_left,
                                           MBRange& best_right,
                                           MBRange& best_both,
                                           MBAdaptiveKDTree::Plane& best_plane,
                                           double eps )
{
  double metric_val = std::numeric_limits<unsigned>::max();
  
  MBErrorCode r;
  const MBCartVect box_min(iter.box_min());
  const MBCartVect box_max(iter.box_max());
  const MBCartVect diff(box_max - box_min);
  
  MBRange entities;
  r = iter.tool()->moab()->get_entities_by_handle( iter.handle(), entities );
  if (MB_SUCCESS != r)
    return r;
  const size_t p_count = entities.size();
  
  for (int axis = 0; axis < 3; ++axis) {
    int plane_count = num_planes;
    if ((num_planes+1)*eps >= diff[axis])
      plane_count = (int)(diff[axis] / eps) - 1;
  
    for (int p = 1; p <= plane_count; ++p) {
      MBAdaptiveKDTree::Plane plane = { box_min[axis] + (p/(1.0+plane_count)) * diff[axis], axis };
      MBRange left, right, both;
      double val;
      r = intersect_children_with_elems( iter.tool()->moab(),
                                         entities, plane,
                                         box_min, box_max,
                                         left, right, both, 
                                         val );
      if (MB_SUCCESS != r)
        return r;
      const size_t diff = p_count - both.size();
      if (left.size() == diff || right.size() == diff)
        continue;
      
      if (val >= metric_val)
        continue;
      
      metric_val = val;
      best_plane = plane;
      best_left.swap(left);
      best_right.swap(right);
      best_both.swap(both);
    }
  }
      
  return MB_SUCCESS;
}


static MBErrorCode best_subdivision_snap_plane( int num_planes,
                                           const MBAdaptiveKDTreeIter& iter,
                                           MBRange& best_left,
                                           MBRange& best_right,
                                           MBRange& best_both,
                                           MBAdaptiveKDTree::Plane& best_plane,
                                           std::vector<double>& tmp_data,
                                           double eps )
{
  double metric_val = std::numeric_limits<unsigned>::max();
  
  MBErrorCode r;
  const MBCartVect box_min(iter.box_min());
  const MBCartVect box_max(iter.box_max());
  const MBCartVect diff(box_max - box_min);
  
  MBRange entities, vertices;
  r = iter.tool()->moab()->get_entities_by_handle( iter.handle(), entities );
  if (MB_SUCCESS != r)
    return r;
  const size_t p_count = entities.size();
  r = iter.tool()->moab()->get_adjacencies( entities, 0, false, vertices, MBInterface::UNION );
  if (MB_SUCCESS != r)
    return r;

  tmp_data.resize( vertices.size() );
  for (int axis = 0; axis < 3; ++axis) {
    int plane_count = num_planes;
    if ((num_planes+1)*eps >= diff[axis])
      plane_count = (int)(diff[axis] / eps) - 1;

    double *ptrs[] = { 0, 0, 0 };
    ptrs[axis] = &tmp_data[0];
    r = iter.tool()->moab()->get_coords( vertices, ptrs[0], ptrs[1], ptrs[2] );
    if (MB_SUCCESS != r)
      return r;
  
    for (int p = 1; p <= plane_count; ++p) {
      double coord = box_min[axis] + (p/(1.0+plane_count)) * diff[axis];
      double closest_coord = tmp_data[0];
      for (unsigned i = 1; i < tmp_data.size(); ++i) 
        if (fabs(coord-tmp_data[i]) < fabs(coord-closest_coord))
          closest_coord = tmp_data[i];
      if (closest_coord <= box_min[axis] || closest_coord >= box_max[axis])
        continue;
      MBAdaptiveKDTree::Plane plane = { closest_coord, axis };
      MBRange left, right, both;
      double val;
      r = intersect_children_with_elems( iter.tool()->moab(),
                                         entities, plane,
                                         box_min, box_max,
                                         left, right, both, 
                                         val );
      if (MB_SUCCESS != r)
        return r;
      const size_t diff = p_count - both.size();
      if (left.size() == diff || right.size() == diff)
        continue;
      
      if (val >= metric_val)
        continue;
      
      metric_val = val;
      best_plane = plane;
      best_left.swap(left);
      best_right.swap(right);
      best_both.swap(both);
    }
  }
      
  return MB_SUCCESS;
}

static MBErrorCode best_vertex_median_plane( int num_planes,
                                           const MBAdaptiveKDTreeIter& iter,
                                           MBRange& best_left,
                                           MBRange& best_right,
                                           MBRange& best_both,
                                           MBAdaptiveKDTree::Plane& best_plane,
                                           std::vector<double>& coords,
                                           double eps)
{
  double metric_val = std::numeric_limits<unsigned>::max();
  
  MBErrorCode r;
  const MBCartVect box_min(iter.box_min());
  const MBCartVect box_max(iter.box_max());
  
  MBRange entities, vertices;
  r = iter.tool()->moab()->get_entities_by_handle( iter.handle(), entities );
  if (MB_SUCCESS != r)
    return r;
  const size_t p_count = entities.size();
  r = iter.tool()->moab()->get_adjacencies( entities, 0, false, vertices, MBInterface::UNION );
  if (MB_SUCCESS != r)
    return r;

  coords.resize( vertices.size() );
  for (int axis = 0; axis < 3; ++axis) {
    if (box_max[axis] - box_min[axis] <= 2*eps)
      continue;
  
    double *ptrs[] = { 0, 0, 0 };
    ptrs[axis] = &coords[0];
    r = iter.tool()->moab()->get_coords( vertices, ptrs[0], ptrs[1], ptrs[2] );
    if (MB_SUCCESS != r)
      return r;
  
    std::sort( coords.begin(), coords.end() );
    std::vector<double>::iterator citer;
    citer = std::upper_bound( coords.begin(), coords.end(), box_min[axis] + eps );
    const size_t count = std::upper_bound( citer, coords.end(), box_max[axis] - eps ) - citer;
    size_t step;
    if ((int)count < 2*num_planes) {
      step = 1; num_planes = count - 1;
    }
    else {
      step = count / (num_planes + 1);
    }
  
    for (int p = 1; p <= num_planes; ++p) {
      
      citer += step;
      MBAdaptiveKDTree::Plane plane = { *citer, axis };
      MBRange left, right, both;
      double val;
      r = intersect_children_with_elems( iter.tool()->moab(),
                                         entities, plane,
                                         box_min, box_max,
                                         left, right, both, 
                                         val );
      if (MB_SUCCESS != r)
        return r;
      const size_t diff = p_count - both.size();
      if (left.size() == diff || right.size() == diff)
        continue;
      
      if (val >= metric_val)
        continue;
      
      metric_val = val;
      best_plane = plane;
      best_left.swap(left);
      best_right.swap(right);
      best_both.swap(both);
    }
  }
      
  return MB_SUCCESS;
}


static MBErrorCode best_vertex_sample_plane( int num_planes,
                                           const MBAdaptiveKDTreeIter& iter,
                                           MBRange& best_left,
                                           MBRange& best_right,
                                           MBRange& best_both,
                                           MBAdaptiveKDTree::Plane& best_plane,
                                           std::vector<double>& coords,
                                           std::vector<size_t>& indices,
                                           double eps )
{
  double metric_val = std::numeric_limits<unsigned>::max();
  
  MBErrorCode r;
  const MBCartVect box_min(iter.box_min());
  const MBCartVect box_max(iter.box_max());
  
  MBRange entities, vertices;
  r = iter.tool()->moab()->get_entities_by_handle( iter.handle(), entities );
  if (MB_SUCCESS != r)
    return r;
  const size_t p_count = entities.size();
  r = iter.tool()->moab()->get_adjacencies( entities, 0, false, vertices, MBInterface::UNION );
  if (MB_SUCCESS != r)
    return r;

  coords.resize( vertices.size() );
  for (int axis = 0; axis < 3; ++axis) {
    if (box_max[axis] - box_min[axis] <= 2*eps)
      continue;
  
    double *ptrs[] = { 0, 0, 0 };
    ptrs[axis] = &coords[0];
    r = iter.tool()->moab()->get_coords( vertices, ptrs[0], ptrs[1], ptrs[2] );
    if (MB_SUCCESS != r)
      return r;
      
    size_t num_valid_coords = 0;
    for (size_t i = 0; i < coords.size(); ++i) 
      if (coords[i] > box_min[axis]+eps && coords[i] < box_max[axis]-eps)
        ++num_valid_coords;
      
    if (2*(size_t)num_planes > num_valid_coords) {
      indices.clear();
      for (size_t i = 0; i < coords.size(); ++i) 
        if (coords[i] > box_min[axis]+eps && coords[i] < box_max[axis]-eps)
          indices.push_back( i );
    }
    else {
      indices.resize( num_planes );
        // make sure random indices are sufficient to cover entire range
      const int num_rand = coords.size() / RAND_MAX + 1;
      for (int j = 0; j < num_planes; ++j)
      {
        size_t rnd;
        do { 
          size_t rnd = rand();
          for (int i = num_rand; i > 1; --i)
            rnd *= rand();
          rnd %= coords.size();
        } while (coords[rnd] <= box_min[axis]+eps || coords[rnd] >= box_max[axis]-eps);
        indices[j] = rnd;
      }
    }
  
    for (unsigned p = 0; p <= indices.size(); ++p) {
      
      MBAdaptiveKDTree::Plane plane = { coords[indices[p]], axis };
      MBRange left, right, both;
      double val;
      r = intersect_children_with_elems( iter.tool()->moab(),
                                         entities, plane,
                                         box_min, box_max,
                                         left, right, both, 
                                         val );
      if (MB_SUCCESS != r)
        return r;
      const size_t diff = p_count - both.size();
      if (left.size() == diff || right.size() == diff)
        continue;
      
      if (val >= metric_val)
        continue;
      
      metric_val = val;
      best_plane = plane;
      best_left.swap(left);
      best_right.swap(right);
      best_both.swap(both);
    }
  }
      
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTree::build_tree( MBRange& elems,
                                       MBEntityHandle& root_set_out,
                                       const Settings* settings_ptr )
{
  Settings settings;
  if (settings_ptr)
    settings = *settings_ptr;
  if (settings.maxEntPerLeaf < 1)
    settings.maxEntPerLeaf = 1;
  if (settings.maxTreeDepth < 1)
    settings.maxTreeDepth = std::numeric_limits<unsigned>::max();
  if (settings.candidateSplitsPerDir < 1)
    settings.candidateSplitsPerDir = 1;
  
    // calculate bounding box of elements
    
  std::vector<double> tmp_data;
  std::vector<size_t> tmp_data2;
  MBRange vertices;
  MBErrorCode rval = moab()->get_adjacencies( elems, 0, false, vertices, MBInterface::UNION );
  if (MB_SUCCESS != rval)
    return rval;
  
  MBCartVect bmin(HUGE_VAL), bmax(-HUGE_VAL), coords;
  for (MBRange::iterator i = vertices.begin(); i != vertices.end(); ++i) {
    rval = moab()->get_coords( &*i, 1, coords.array() );
    if (MB_SUCCESS != rval)
      return rval;
    for (unsigned j = 0; j < 3; ++j) {
      if (coords[j] < bmin[j])
        bmin[j] = coords[j];
      if (coords[j] > bmax[j])
        bmax[j] = coords[j];
    }
  }
  
    // create tree root
  rval = create_tree( bmin.array(), bmax.array(), root_set_out );
  if (MB_SUCCESS != rval)
    return rval;
  rval = moab()->add_entities( root_set_out, elems );
  if (MB_SUCCESS != rval)
    return rval;
  
  MBAdaptiveKDTreeIter iter;
  iter.initialize( this, root_set_out, bmin.array(), bmax.array(), MBAdaptiveKDTreeIter::LEFT );
  
  for (;;) {
  
    int pcount;
    rval = moab()->get_number_entities_by_handle( iter.handle(), pcount );
    if (MB_SUCCESS != rval)
      break;

    const size_t p_count = pcount;
    MBRange best_left, best_right, best_both;
    Plane best_plane = { HUGE_VAL, -1 };
    if (p_count > settings.maxEntPerLeaf && iter.depth() < settings.maxTreeDepth) {
      switch (settings.candidatePlaneSet) {
        case MBAdaptiveKDTree::SUBDIVISION:
          rval = best_subdivision_plane( settings.candidateSplitsPerDir, 
                                               iter, 
                                               best_left, 
                                               best_right, 
                                               best_both, 
                                               best_plane, 
                                               settings.minBoxWidth );
          break;
        case MBAdaptiveKDTree::SUBDIVISION_SNAP:
          rval = best_subdivision_snap_plane( settings.candidateSplitsPerDir, 
                                               iter, 
                                               best_left, 
                                               best_right, 
                                               best_both, 
                                               best_plane, 
                                               tmp_data, 
                                               settings.minBoxWidth );
          break;
        case MBAdaptiveKDTree::VERTEX_MEDIAN:
          rval = best_vertex_median_plane( settings.candidateSplitsPerDir, 
                                               iter, 
                                               best_left, 
                                               best_right, 
                                               best_both, 
                                               best_plane, 
                                               tmp_data, 
                                               settings.minBoxWidth );
          break;
        case MBAdaptiveKDTree::VERTEX_SAMPLE:
          rval = best_vertex_sample_plane( settings.candidateSplitsPerDir, 
                                               iter, 
                                               best_left, 
                                               best_right, 
                                               best_both, 
                                               best_plane, 
                                               tmp_data, 
                                               tmp_data2,
                                               settings.minBoxWidth );
          break;
        default:
          rval = MB_FAILURE;
      }
    
      if (MB_SUCCESS != rval)
        return rval;
    }
    
    if (best_plane.norm >= 0) {
      best_left.merge( best_both );
      best_right.merge( best_both );
      rval = split_leaf( iter, best_plane,  best_left, best_right );
      if (MB_SUCCESS != rval)
        return rval;
    }
    else {
      rval = iter.step();
      if (MB_ENTITY_NOT_FOUND == rval) 
        return MB_SUCCESS;  // at end
      else if (MB_SUCCESS != rval)
        break;
    }
  }
  
  delete_tree( root_set_out );
  return rval;
}

MBErrorCode MBAdaptiveKDTree::leaf_containing_point( MBEntityHandle tree_root,
                                                     const double point[3],
                                                     MBEntityHandle& leaf_out )
{
  std::vector<MBEntityHandle> children;
  Plane plane;
  MBEntityHandle node = tree_root;
  MBErrorCode rval = moab()->get_child_meshsets( node, children );
  if (MB_SUCCESS != rval)
    return rval;
  while (!children.empty()) {
    rval = get_split_plane( node, plane );
    if (MB_SUCCESS != rval)
      return rval;
      
    const double d = point[plane.norm] - plane.coord;
    node = children[(d > 0.0)];
    
    children.clear();
    rval = moab()->get_child_meshsets( node, children );
    if (MB_SUCCESS != rval)
      return rval;
  }
  leaf_out = node;
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTree::leaf_containing_point( MBEntityHandle root,
                                                     const double point[3],
                                                     MBAdaptiveKDTreeIter& result )
{
    // get bounding box of tree
  MBErrorCode rval = moab()->tag_get_data( rootTag, &root, 1, result.mBox );
  if (MB_SUCCESS != rval)
    return rval;
    
    // test that point is inside tree
  if (point[0] < result.box_min()[0] || point[0] > result.box_max()[0] ||
      point[1] < result.box_min()[1] || point[1] > result.box_max()[1] ||
      point[2] < result.box_min()[2] || point[2] > result.box_max()[2])
    return MB_ENTITY_NOT_FOUND;  

    // initialize iterator at tree root
  result.treeTool = this;
  result.mStack.clear();
  result.mStack.push_back( MBAdaptiveKDTreeIter::StackObj(root,0) );
    
    // loop until we reach a leaf
  MBAdaptiveKDTree::Plane plane;
  for(;;) {
      // get children
    result.childVect.clear();
    rval = moab()->get_child_meshsets( result.handle(), result.childVect );
    if (MB_SUCCESS != rval)
      return rval;
      
      // if no children, then at leaf (done)
    if (result.childVect.empty())
      break;

      // get split plane
    rval = get_split_plane( result.handle(), plane );
    if (MB_SUCCESS != rval) 
      return rval;
    
      // step iterator to appropriate child
      // idx: 0->left, 1->right
    const int idx = (point[plane.norm] > plane.coord);
    result.mStack.push_back( MBAdaptiveKDTreeIter::StackObj( result.childVect[idx], 
                                                             result.mBox[1-idx][plane.norm] ) );
    result.mBox[1-idx][plane.norm] = plane.coord;
  }
    
  return MB_SUCCESS;
}

struct NodeDistance {
  MBEntityHandle handle;
  MBCartVect dist; // from_point - closest_point_on_box
};

MBErrorCode MBAdaptiveKDTree::leaves_within_distance( MBEntityHandle tree_root,
                                                      const double from_point[3],
                                                      const double distance,
                                                      std::vector<MBEntityHandle>& result_list )
{
  const double dist_sqr = distance * distance;
  const MBCartVect from(from_point);
  std::vector<NodeDistance> list;     // list of subtrees to traverse
    // pre-allocate space for default max tree depth
  Settings tmp_settings;
  list.reserve( tmp_settings.maxTreeDepth );

    // misc temporary values
  Plane plane;
  NodeDistance node; 
  MBErrorCode rval;
  std::vector<MBEntityHandle> children;
  
    // Get distance from input position to bounding box of tree
    // (zero if inside box)
  double min[3], max[3];
  rval = get_tree_box( tree_root, min, max );
    // if bounding box is not available (e.g. not starting from true root)
    // just start with zero.  Less efficient, but will work.
  node.dist = MBCartVect(0.0);
  if (MB_SUCCESS == rval) {
    for (int i = 0; i < 3; ++i) {
      if (from_point[i] < min[i])
        node.dist[i] = min[i] - from_point[i];
      else if (from_point[i] > max[i])
        node.dist[i] = from_point[i] - max[i];
    }
    if (node.dist % node.dist > dist_sqr)
      return MB_SUCCESS;
  }
  
    // begin with root in list  
  node.handle = tree_root;
  list.push_back( node );
  
  while( !list.empty() ) {

    node = list.back();
    list.pop_back();
      
      // If leaf node, test contained triangles
    children.clear();
    rval = moab()->get_child_meshsets( node.handle, children );
    if (children.empty()) {
      result_list.push_back( node.handle );
      continue;
    }
      
      // If not leaf node, add children to working list
    rval = get_split_plane( node.handle, plane );
    if (MB_SUCCESS != rval)
      return rval;
    
    const double d = from[plane.norm] - plane.coord;
    
      // right of plane?
    if (d > 0) {
      node.handle = children[1];
      list.push_back( node );
        // if the split plane is close to the input point, add
        // the left child also (we'll check the exact distance
        /// when we pop it from the list.)
      if (d <= distance) {
        node.dist[plane.norm] = d;
        if (node.dist % node.dist <= dist_sqr) {
          node.handle = children[0];
          list.push_back( node );
        }
      }
    }
      // left of plane
    else {
      node.handle = children[0];
      list.push_back( node );
        // if the split plane is close to the input point, add
        // the right child also (we'll check the exact distance
        /// when we pop it from the list.)
      if (-d <= distance) {
        node.dist[plane.norm] = -d;
        if (node.dist % node.dist <= dist_sqr) {
          node.handle = children[1];
          list.push_back( node );
        }
      }
    }
  }

  return MB_SUCCESS;
}

/** Find the triangles in a set that are closer to the input
 *  position than any triangles in the 'closest_tris' list.
 *
 *  closest_tris is assumed to contain a list of triangles for 
 *  which the first is the closest known triangle to the input 
 *  position and the first entry in 'closest_pts' is the closest
 *  location on that triangle.  Any other values in the lists must
 *  be other triangles for which the closest point is within the
 *  input tolernace of the closest closest point.  This function
 *  will update the lists as appropriate if any closer triangles
 *  or triangles within the tolerance of the current closest location
 *  are found.  The fisrt entry is maintaned as the closest of the
 *  list of triangles.
 */
static MBErrorCode closest_to_triangles( MBInterface* moab,
                                         MBEntityHandle set_handle,
                                         double tolerance,
                                         const MBCartVect& from,
                                         std::vector<MBEntityHandle>& closest_tris,
                                         std::vector<MBCartVect>& closest_pts )
{
  MBErrorCode rval;
  MBRange tris;
  MBCartVect pos, diff, verts[3];
  const MBEntityHandle* conn;
  int len;
  double shortest_dist_sqr = HUGE_VAL;
  if (!closest_pts.empty()) {
    diff = from - closest_pts.front();
    shortest_dist_sqr = diff % diff;
  }
  
  rval = moab->get_entities_by_type( set_handle, MBTRI, tris );
  if (MB_SUCCESS != rval)
    return rval;
      
  for (MBRange::iterator i = tris.begin(); i != tris.end(); ++i) {
    rval = moab->get_connectivity( *i, conn, len );
    if (MB_SUCCESS != rval)
      return rval;

    rval = moab->get_coords( conn, 3, verts[0].array() );
    if (MB_SUCCESS != rval)
      return rval;

    MBGeomUtil::closest_location_on_tri( from, verts, pos );
    diff = pos - from;
    double dist_sqr = diff % diff;
    if (dist_sqr < shortest_dist_sqr) {
        // new closest location
      shortest_dist_sqr = dist_sqr;

      if (closest_pts.empty()) {
        closest_tris.push_back( *i );
        closest_pts.push_back( pos );
      }
        // if have a previous closest location
      else {
          // if previous closest is more than 2*tolerance away
          // from new closest, then nothing in the list can
          // be within tolerance of new closest point.
        diff = pos - closest_pts.front();
        dist_sqr = diff % diff;
        if (dist_sqr > 4.0 * tolerance * tolerance) {
          closest_tris.clear();
          closest_pts.clear();
          closest_tris.push_back( *i );
          closest_pts.push_back( pos );
        }
          // otherwise need to remove any triangles that are
          // not within tolerance of the new closest point.
        else {
          unsigned r = 0, w = 0;
          for (r = 0; r < closest_pts.size(); ++r) {
            diff = pos - closest_pts[r];
            if (diff % diff <= tolerance*tolerance) {
              closest_pts[w] = closest_pts[r];
              closest_tris[w] = closest_tris[r];
              ++w;
            }
          }
          closest_pts.resize( w + 1 );
          closest_tris.resize( w + 1 );
            // always put the closest one in the front
          if (w > 0) {
            closest_pts.back() = closest_pts.front();
            closest_tris.back() = closest_tris.front();
          }
          closest_pts.front() = pos;
          closest_tris.front() = *i;
        }
      }
    }
    else {
        // If within tolerance of old closest triangle,
        // add this one to the list.
      diff = closest_pts.front() - pos;
      if (diff % diff <= tolerance*tolerance) {
        closest_pts.push_back( pos );
        closest_tris.push_back( *i );
      }
    }
  }
  
  return MB_SUCCESS;
}

static MBErrorCode closest_to_triangles( MBInterface* moab,
                                         MBEntityHandle set_handle,
                                         const MBCartVect& from,
                                         double& shortest_dist_sqr,
                                         MBCartVect& closest_pt,
                                         MBEntityHandle& closest_tri )
{
  MBErrorCode rval;
  MBRange tris;
  MBCartVect pos, diff, verts[3];
  const MBEntityHandle* conn;
  int len;
  
  rval = moab->get_entities_by_type( set_handle, MBTRI, tris );
  if (MB_SUCCESS != rval)
    return rval;
      
  for (MBRange::iterator i = tris.begin(); i != tris.end(); ++i) {
    rval = moab->get_connectivity( *i, conn, len );
    if (MB_SUCCESS != rval)
      return rval;

    rval = moab->get_coords( conn, 3, verts[0].array() );
    if (MB_SUCCESS != rval)
      return rval;

    MBGeomUtil::closest_location_on_tri( from, verts, pos );
    diff = pos - from;
    double dist_sqr = diff % diff;
    if (dist_sqr < shortest_dist_sqr) {
        // new closest location
      shortest_dist_sqr = dist_sqr;
      closest_pt = pos;
      closest_tri = *i;
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTree::closest_triangle( MBEntityHandle tree_root,
                                 const double from_coords[3],
                                 double closest_point_out[3],
                                 MBEntityHandle& triangle_out )
{
  MBErrorCode rval;
  double shortest_dist_sqr = HUGE_VAL;
  std::vector<MBEntityHandle> leaves;
  const MBCartVect from(from_coords);
  MBCartVect closest_pt;
  
    // Find the leaf containing the input point
    // This search does not take into account any bounding box for the
    // tree, so it always returns one leaf.
  MBEntityHandle leaf;
  rval = leaf_containing_point( tree_root, from_coords, leaf );
  if (MB_SUCCESS != rval) return rval;
  
    // Find the closest triangle(s) in the leaf containing the point
  rval = closest_to_triangles( moab(), leaf, from, shortest_dist_sqr, 
                               closest_pt, triangle_out );
  if (MB_SUCCESS != rval) return rval;
  
    // Find any other leaves for which the bounding box is within
    // the same distance from the input point as the current closest
    // point is.
  MBCartVect diff = closest_pt - from;
  rval = leaves_within_distance( tree_root, from_coords, 
                                 sqrt(diff%diff), leaves );
  if (MB_SUCCESS != rval) return rval;

    // Check any close leaves to see if they contain triangles that
    // are as close to or closer than the current closest triangle(s).
  for (unsigned i = 0; i < leaves.size(); ++i) {
    rval = closest_to_triangles( moab(), leaves[i], from, shortest_dist_sqr, 
                                 closest_pt, triangle_out );
    if (MB_SUCCESS != rval) return rval;
  }
  
    // pass back resulting position
  closest_pt.get( closest_point_out );
  return MB_SUCCESS;
}

MBErrorCode MBAdaptiveKDTree::sphere_intersect_triangles( 
                                   MBEntityHandle tree_root,
                                   const double center[3],
                                   double radius,
                                   std::vector<MBEntityHandle>& triangles )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle> leaves;
  const MBCartVect from(center);
  MBCartVect closest_pt;
  const MBEntityHandle* conn;
  MBCartVect coords[3];
  int conn_len;

    // get leaves of tree that intersect sphere
  rval = leaves_within_distance( tree_root, center, radius, leaves );
  if (MB_SUCCESS != rval) return rval;
  
    // search each leaf for triangles intersecting sphere
  for (unsigned i = 0; i < leaves.size(); ++i) {
    MBRange tris;
    rval = moab()->get_entities_by_type( leaves[i], MBTRI, tris );
    if (MB_SUCCESS != rval) return rval;
    
    for (MBRange::iterator j = tris.begin(); j != tris.end(); ++j) {
      rval = moab()->get_connectivity( *j, conn, conn_len );
      if (MB_SUCCESS != rval) return rval;
      rval = moab()->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval) return rval;
      MBGeomUtil::closest_location_on_tri( from, coords, closest_pt );
      closest_pt -= from;
      if ((closest_pt % closest_pt) <= (radius*radius)) 
        triangles.push_back( *j );
    }
  }
  
    // remove duplicates from triangle list
  std::sort( triangles.begin(), triangles.end() );
  triangles.erase( std::unique( triangles.begin(), triangles.end() ), triangles.end() );
  return MB_SUCCESS;
}
  
      

struct NodeSeg {
  NodeSeg( MBEntityHandle h, double b, double e )
    : handle(h), beg(b), end(e) {}
  MBEntityHandle handle;
  double beg, end;
};

MBErrorCode MBAdaptiveKDTree::ray_intersect_triangles( MBEntityHandle root,
                                                 const double tol,
                                                 const double ray_dir_in[3],
                                                 const double ray_pt_in[3],
                                                 std::vector<MBEntityHandle>& tris_out,
                                                 std::vector<double>& dists_out,
                                                 int max_ints,
                                                 double ray_end )
{
  MBErrorCode rval;
  double ray_beg = 0.0;
  if (ray_end < 0.0)
    ray_end = HUGE_VAL;
  
    // if root has bounding box, trim ray to that box
  MBCartVect bmin, bmax, tvec(tol);
  const MBCartVect ray_pt( ray_pt_in ), ray_dir( ray_dir_in );
  rval = get_tree_box( root, bmin.array(), bmax.array() );
  if (MB_SUCCESS == rval) {
    if (!MBGeomUtil::segment_box_intersect( bmin-tvec, bmax+tvec, ray_pt, ray_dir, ray_beg, ray_end ))
      return MB_SUCCESS; // ray misses entire tree.
  }
  
  MBRange tris;
  MBRange::iterator iter;
  MBCartVect tri_coords[3];
  const MBEntityHandle* tri_conn;
  int conn_len;
  double tri_t;
  
  Plane plane;
  std::vector<MBEntityHandle> children;
  std::vector<NodeSeg> list;
  NodeSeg seg(root, ray_beg, ray_end);
  list.push_back( seg );
  
  while (!list.empty()) {
    seg = list.back();
    list.pop_back();
    
      // If we are limited to a certain number of intersections
      // (max_ints != 0), then ray_end will contain the distance
      // to the furthest intersection we have so far.  If the
      // tree node is further than that, skip it.
    if (seg.beg > ray_end) 
      continue;

      // Check if at a leaf 
    children.clear();
    rval = moab()->get_child_meshsets( seg.handle, children );
    if (MB_SUCCESS != rval)
      return rval;
    if (children.empty()) { // leaf

      tris.clear();
      rval = moab()->get_entities_by_type( seg.handle, MBTRI, tris );
      if (MB_SUCCESS != rval)
        return rval;
    
      for (iter = tris.begin(); iter != tris.end(); ++iter) {
        rval = moab()->get_connectivity( *iter, tri_conn, conn_len );
        if (MB_SUCCESS != rval) return rval;
        rval = moab()->get_coords( tri_conn, 3, tri_coords[0].array() );
        if (MB_SUCCESS != rval) return rval;
        
        if (MBGeomUtil::ray_tri_intersect( tri_coords, ray_pt, ray_dir, tol, tri_t, &ray_end )) {
          if (!max_ints) {
            if (std::find(tris_out.begin(),tris_out.end(),*iter) == tris_out.end()) {
              tris_out.push_back( *iter );
              dists_out.push_back( tri_t );
            }
          } 
          else if (tri_t < ray_end) {
            if (std::find(tris_out.begin(),tris_out.end(),*iter) == tris_out.end()) {
              if (tris_out.size() < (unsigned)max_ints) {
                tris_out.resize( tris_out.size() + 1 );
                dists_out.resize( dists_out.size() + 1 );
              }
              int w = tris_out.size() - 1;
              for (; w > 0 && tri_t < dists_out[w-1]; --w) {
                tris_out[w] = tris_out[w-1];
                dists_out[w] = dists_out[w-1];
              }
              tris_out[w] = *iter;
              dists_out[w] = tri_t;
              ray_end = dists_out.back();
            }
          }
        }
      }

      continue;
    }
    
    rval = get_split_plane( seg.handle, plane );
    if (MB_SUCCESS != rval)
      return rval;
    
    const double t = (plane.coord - ray_pt[plane.norm]) / ray_dir[plane.norm];
    if (!finite(t)) {         // ray parallel to plane
      if (ray_pt[plane.norm] - tol <= plane.coord)
        list.push_back( NodeSeg( children[0], seg.beg, seg.end ) );
      if (ray_pt[plane.norm] + tol >= plane.coord)
        list.push_back( NodeSeg( children[1], seg.beg, seg.end ) );
    }
    else if (ray_dir[plane.norm] < 0.0) {
      if (seg.beg > t) {      // segment left of plane
        list.push_back( NodeSeg( children[0], seg.beg, seg.end ) );
//        if (plane.coord - ray_pt[plane.norm] + ray_dir[plane.norm] * seg.beg < tol)
//          list.push_back( NodeSeg( children[1], seg.beg, seg.end ) );
      }
      else if (seg.end < t) { // segment right of plane
        list.push_back( NodeSeg( children[1], seg.beg, seg.end ) );
//        if (ray_pt[plane.norm] + ray_dir[plane.norm] * seg.end - plane.coord < tol)
//          list.push_back( NodeSeg( children[0], seg.beg, seg.end ) );
      }
      else {                  // segment crosses plane
        list.push_back( NodeSeg( children[1], seg.beg, t ) );
        list.push_back( NodeSeg( children[0], t, seg.end ) );
      }
    }
    else {
      if (seg.beg > t) {      // segment right of plane
        list.push_back( NodeSeg( children[1], seg.beg, seg.end ) );
//        if (ray_pt[plane.norm] + ray_dir[plane.norm] * seg.beg - plane.coord < tol)
//          list.push_back( NodeSeg( children[0], seg.beg, seg.end ) );
      }
      else if (seg.end < t) { // segment left of plane
        list.push_back( NodeSeg( children[0], seg.beg, seg.end ) );
//        if (plane.coord - ray_pt[plane.norm] + ray_dir[plane.norm] * seg.end < tol)
//          list.push_back( NodeSeg( children[1], seg.beg, seg.end ) );
      }
      else {                  // segment crosses plane
        list.push_back( NodeSeg( children[0], seg.beg, t ) );
        list.push_back( NodeSeg( children[1], t, seg.end ) );
      }
    }
  }
  
  return MB_SUCCESS;
}
          
