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

/**\file MBOrientedBox.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 */

#include "MBInterface.hpp"
#include "MBOrientedBoxTreeTool.hpp"
#include "MBOrientedBox.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"
#include "MBGeomUtil.hpp"
#include "MBTagConventions.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <assert.h>

#define MB_OBB_USE_VECTOR_QUERIES
//#define MB_OBB_USE_TYPE_QUERIES

#if defined(MB_OBB_USE_VECTOR_QUERIES) && defined(MB_OBB_USE_TYPE_QUERIES)
# undef MB_OBB_USE_TYPE_QUERIES
#endif

const char DEFAULT_TAG_NAME[] = "OBB";

MBOrientedBoxTreeTool::Op::~Op() {}

MBOrientedBoxTreeTool::MBOrientedBoxTreeTool( MBInterface* i,
                                              const char* tag_name )
  : instance( i )
{
  if (!tag_name)
    tag_name = DEFAULT_TAG_NAME;
  MBErrorCode rval = MBOrientedBox::tag_handle( tagHandle, instance, tag_name, true );
  if (MB_SUCCESS != rval)
    tagHandle = 0;
}

MBOrientedBoxTreeTool::Settings::Settings() 
  : max_leaf_entities( 8 ),
    max_depth( 0 ),
    worst_split_ratio( 0.95 ),
    best_split_ratio( 0.4 ),
    set_options( MESHSET_SET )
  {}

bool MBOrientedBoxTreeTool::Settings::valid() const
{
  return max_leaf_entities > 0 
      && max_depth >= 0
      && worst_split_ratio <= 1.0
      && best_split_ratio >= 0.0
      && worst_split_ratio >= best_split_ratio
      ;
}


/********************** Simple Tree Access Methods ****************************/


MBErrorCode MBOrientedBoxTreeTool::box( MBEntityHandle set, MBOrientedBox& obb )
{
  return instance->tag_get_data( tagHandle, &set, 1, &obb );
}

MBErrorCode MBOrientedBoxTreeTool::box( MBEntityHandle set,
                                        double center[3],
                                        double axis1[3],
                                        double axis2[3],
                                        double axis3[3] )
{
  MBOrientedBox obb;
  MBErrorCode rval = this->box( set, obb );
  obb.center.get( center );
  obb.scaled_axis(0).get( axis1 );
  obb.scaled_axis(1).get( axis2 );
  obb.scaled_axis(2).get( axis3 );
  return rval;
}


/********************** Tree Construction Code ****************************/


struct MBOrientedBoxTreeTool::SetData {
  MBEntityHandle handle;
  MBOrientedBox::CovarienceData box_data;
  //MBRange vertices;
};


MBErrorCode MBOrientedBoxTreeTool::build( const MBRange& entities,
                                          MBEntityHandle& set_handle_out,
                                          const Settings* settings )
{
  if (!entities.all_of_dimension(2))
    return MB_TYPE_OUT_OF_RANGE;
  if (settings && !settings->valid())
    return MB_FAILURE;
    
  return build_tree( entities, set_handle_out, 0, 
                     settings ? *settings : Settings() );
}

MBErrorCode MBOrientedBoxTreeTool::join_trees( const MBRange& sets,
                                               MBEntityHandle& set_handle_out,
                                               const Settings* settings )
{
  if (!sets.all_of_type(MBENTITYSET))
    return MB_TYPE_OUT_OF_RANGE;
  if (settings && !settings->valid())
    return MB_FAILURE;
  
    // Build initial set data list.
  std::list<SetData> data;
  for (MBRange::const_iterator i = sets.begin(); i != sets.end(); ++i) {
    MBRange elements;
    MBErrorCode rval = instance->get_entities_by_dimension( *i, 2, elements, true );
    if (MB_SUCCESS != rval)
      return rval;
    if (elements.empty())
      continue;
    
    data.push_back( SetData() );
    SetData& set_data = data.back();
    set_data.handle = *i;
    rval = MBOrientedBox::covarience_data_from_tris( set_data.box_data, instance, elements );
    if (MB_SUCCESS != rval)
      return rval;
  }

  return build_sets( data, set_handle_out, 0, 
                     settings ? *settings : Settings() );
}
  

/**\brief Split trianges by which side of a plane they are on
 *
 * Given a plane specified as a bisecting plane normal to one
 * of the axes of a box, split triangles based on which side
 * of the plane they are on.
 *\param instance   MOAB instance
 *\param box        The oriented box containing all the entities
 *\param axis       The axis for which the split plane is orthogonal
 *\param left_list  Output, entities to the left of the plane
 *\param right_list Output, entities to the right of the plane
 *\param num_intersecting Output, number entities intersecting plane
 */
static MBErrorCode split_box( MBInterface* instance, 
                              const MBOrientedBox& box, 
                              int axis, 
                              const MBRange& entities, 
                              MBRange& left_list, 
                              MBRange& right_list )
{
  MBErrorCode rval;
  left_list.clear();
  right_list.clear();

  std::vector<MBCartVect> coords;
  for (MBRange::reverse_iterator i = entities.rbegin(); i != entities.rend(); ++i) {
    const MBEntityHandle *conn;
    int conn_len;
    rval = instance->get_connectivity( *i, conn, conn_len );
    if (MB_SUCCESS != rval)
      return rval;
    
    coords.resize( conn_len );
    rval = instance->get_coords( conn, conn_len, coords[0].array() );
    if (MB_SUCCESS != rval)
      return rval;
    
    MBCartVect centroid(0.0);
    for (int j = 0; j < conn_len; ++j)
      centroid += coords[j];
    centroid /= conn_len;
    
    if ((box.axis[axis] % (centroid - box.center)) < 0.0)
      left_list.insert( *i );
    else
      right_list.insert( *i );
  }
  
  return MB_SUCCESS;
}


MBErrorCode MBOrientedBoxTreeTool::build_tree( const MBRange& entities,
                                               MBEntityHandle& set,
                                               int depth,
                                               const Settings& settings )
{
  MBOrientedBox box;
  MBErrorCode rval;
  
  rval = MBOrientedBox::compute_from_2d_cells( box, instance, entities );
  if (MB_SUCCESS != rval)
    return rval;
  
    // create an entity set for the tree node
  rval = instance->create_meshset( settings.set_options, set );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = instance->tag_set_data( tagHandle, &set, 1, &box );
  if (MB_SUCCESS != rval) 
    { delete_tree( set ); return rval; }
  
    // check if should create children
  bool leaf = true;
  ++depth;
  if ((!settings.max_depth || depth < settings.max_depth) && 
      entities.size() > (unsigned)settings.max_leaf_entities) {
      // try splitting with planes normal to each axis of the box
      // until we find an acceptable split
    double best_ratio = settings.worst_split_ratio; // worst case ratio
    MBRange best_left_list, best_right_list;
      // Axes are sorted from shortest to longest, so search backwards
    for (int axis = 2; best_ratio > settings.best_split_ratio && axis >= 0; --axis) {
      MBRange left_list, right_list;

      rval = split_box( instance, box, axis, entities, left_list, right_list );
      if (MB_SUCCESS != rval) 
        { delete_tree( set ); return rval; }
        
      double ratio = fabs((double)right_list.size() - left_list.size()) / entities.size();
      
      if (ratio < best_ratio) {
        best_ratio = ratio;
        best_left_list.swap( left_list );
        best_right_list.swap( right_list );
      }
    }
    
      // create children
    if (!best_left_list.empty())
    {
      MBEntityHandle child = 0;
      
      rval = build_tree( best_left_list, child, depth, settings );
      if (MB_SUCCESS != rval)
        { delete_tree( set ); return rval; }
      rval = instance->add_child_meshset( set, child );
      if (MB_SUCCESS != rval)
        { delete_tree( set ); delete_tree( child ); return rval; }
      
      rval = build_tree( best_right_list, child, depth, settings );
      if (MB_SUCCESS != rval)
        { delete_tree( set ); return rval; }
      rval = instance->add_child_meshset( set, child );
      if (MB_SUCCESS != rval)
        { delete_tree( set ); delete_tree( child ); return rval; }
      
      leaf = false;
    }
  }
  
  if (leaf)
  {
    rval = instance->add_entities( set, entities );
    if (MB_SUCCESS != rval) 
      { delete_tree( set ); return rval; }
  }
  
  return MB_SUCCESS;
}


static MBErrorCode split_sets( MBInterface* , 
                               const MBOrientedBox& box, 
                               int axis, 
                               const std::list<MBOrientedBoxTreeTool::SetData>& sets,
                               std::list<MBOrientedBoxTreeTool::SetData>& left,
                               std::list<MBOrientedBoxTreeTool::SetData>& right )
{
  left.clear();
  right.clear();

  std::list<MBOrientedBoxTreeTool::SetData>::const_iterator i;
  for (i = sets.begin(); i != sets.end(); ++i) {
    MBCartVect centroid(i->box_data.center / i->box_data.area);
    if ((box.axis[axis] % (centroid - box.center)) < 0.0)
      left.push_back( *i );
    else
      right.push_back( *i );
  }
  
  return MB_SUCCESS;
}


MBErrorCode MBOrientedBoxTreeTool::build_sets( std::list<SetData>& sets,
                                               MBEntityHandle& node_set,
                                               int depth,
                                               const Settings& settings )
{
  MBErrorCode rval;
  int count = sets.size();
  
    // calculate box
  MBOrientedBox box;

  // make vector go out of scope when done, so memory is released
  { 
    MBRange elems;
    std::vector<MBOrientedBox::CovarienceData> data(sets.size());
    data.clear();
    for (std::list<SetData>::iterator i = sets.begin(); i != sets.end(); ++i) {
      data.push_back( i->box_data );
      rval = instance->get_entities_by_dimension( i->handle, 2, elems, true );
      if (MB_SUCCESS != rval)
        return rval;
    }
    
    MBRange points;
    rval = instance->get_adjacencies( elems, 0, false, points, MBInterface::UNION );
    if (MB_SUCCESS != rval)
      return rval;
    
    rval = MBOrientedBox::compute_from_covarience_data( box, instance, &data[0], data.size(), points );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
    // If only one set in list...
  if (count == 1) {
    node_set = sets.front().handle;
    return instance->tag_set_data( tagHandle, &node_set, 1, &box );
  }
  
    // create an entity set for the tree node
  rval = instance->create_meshset( settings.set_options, node_set );
  if (MB_SUCCESS != rval)
    return rval;
  
  rval = instance->tag_set_data( tagHandle, &node_set, 1, &box );
  if (MB_SUCCESS != rval) 
    { delete_tree( node_set ); return rval; }
  
  double best_ratio = 2.0; 
  std::list<SetData> best_left_list, best_right_list;
  for (int axis = 0; axis < 2; ++axis) {
    std::list<SetData> left_list, right_list;
    rval = split_sets( instance, box, axis, sets, left_list, right_list );
    if (MB_SUCCESS != rval) 
      { delete_tree( node_set ); return rval; }

    double ratio = fabs((double)right_list.size() - left_list.size()) / sets.size();

    if (ratio < best_ratio) {
      best_ratio = ratio;
      best_left_list.swap( left_list );
      best_right_list.swap( right_list );
    }
  }
  
    // We must subdivide the list of sets, because we want to guarantee that
    // there is a node in the tree corresponding to each of the sets.  So if
    // we couldn't find a usable split plane, just split them in an arbitrary
    // manner.
  if (best_left_list.empty() || best_right_list.empty()) {
    best_left_list.clear();
    best_right_list.clear();
    std::list<SetData>* lists[2] = {&best_left_list,&best_right_list};
    int i = 0;
    while (!sets.empty()) {
      lists[i]->push_back( sets.front() );
      sets.pop_front();
      i = 1 - i;
    }
  }
  else {
    sets.clear(); // release memory before recursion
  }
  
    // Create child sets
    
  MBEntityHandle child = 0;
      
  rval = build_sets( best_left_list, child, depth+1, settings );
  if (MB_SUCCESS != rval)
    { delete_tree( node_set ); return rval; }
  rval = instance->add_child_meshset( node_set, child );
  if (MB_SUCCESS != rval)
    { delete_tree( node_set ); delete_tree( child ); return rval; }

  rval = build_sets( best_right_list, child, depth+1, settings );
  if (MB_SUCCESS != rval)
    { delete_tree( node_set ); return rval; }
  rval = instance->add_child_meshset( node_set, child );
  if (MB_SUCCESS != rval)
    { delete_tree( node_set ); delete_tree( child ); return rval; }
  
  return MB_SUCCESS;
}

MBErrorCode MBOrientedBoxTreeTool::delete_tree( MBEntityHandle set )
{
  std::vector<MBEntityHandle> children;
  MBErrorCode rval = instance->get_child_meshsets( set, children, 0 );
  if (MB_SUCCESS != rval)
    return rval;
  
  children.insert( children.begin(), set );
  return instance->delete_entities( &children[0], children.size() );
}


/********************** Generic Tree Traversal ****************************/


struct Data { MBEntityHandle set; int depth; };
MBErrorCode MBOrientedBoxTreeTool::preorder_traverse( MBEntityHandle set,
                                                      Op& operation )
{
  MBErrorCode rval;
  std::vector<MBEntityHandle> children;
  std::vector<Data> the_stack;
  Data data = { set, 0 };
  the_stack.push_back( data );
  
  while (!the_stack.empty())
  {
    data = the_stack.back();
    the_stack.pop_back();
    
    bool descend = true;
    rval = operation.visit( data.set, data.depth, descend );
    if (MB_SUCCESS != rval)
      return rval;
    
    if (!descend)
      continue;
    
    children.clear();
    rval = instance->get_child_meshsets( data.set, children );
    if (MB_SUCCESS != rval)
      return rval;
    if (children.empty()) {
      rval = operation.leaf( data.set );
      if (MB_SUCCESS != rval)
        return rval;
    }
    else if (children.size() == 2) {
      data.depth++;
      data.set = children[0];
      the_stack.push_back( data );
      data.set = children[1];
      the_stack.push_back( data );
    }
    else
      return MB_MULTIPLE_ENTITIES_FOUND;
  }
  
  return MB_SUCCESS;
}


/********************** General Ray/Tree and Ray/Triangel Intersection ***************/


class RayIntersector : public MBOrientedBoxTreeTool::Op
{
  private:
    MBOrientedBoxTreeTool* tool;
    const MBCartVect b, m;
    const double* len;
    const double tol;
    MBRange& boxes;
    
  public:
    RayIntersector( MBOrientedBoxTreeTool* tool_ptr,
                    const double* ray_point,
                    const double* unit_ray_dir,
                    const double *ray_length,
                    double tolerance,
                    MBRange& leaf_boxes )
      : tool(tool_ptr),
        b(ray_point), m(unit_ray_dir),
        len(ray_length), tol(tolerance),
        boxes(leaf_boxes) 
      { }
  
    virtual MBErrorCode visit( MBEntityHandle node,
                               int depth,
                               bool& descend );
    virtual MBErrorCode leaf( MBEntityHandle node );
};

//#include <stdio.h>
//inline void dump_fragmentation( const MBRange& range ) {
//  static FILE* file = fopen( "fragmentation", "w" );
//  unsigned ranges = 0, entities = 0;
//  for (MBRange::const_pair_iterator i = range.const_pair_begin(); i != range.const_pair_end(); ++i)
//  {
//    ++ranges;
//    entities += i->second - i->first + 1;
//  }
//  fprintf( file, "%u %u\n", ranges, entities );
//}

MBErrorCode MBOrientedBoxTreeTool::ray_intersect_triangles( 
                          std::vector<double>& intersection_distances_out,
                          const MBRange& boxes,
                          double tolerance,
                          const double ray_point[3],
                          const double unit_ray_dir[3],
                          const double* ray_length )
{
  MBErrorCode rval;
  intersection_distances_out.clear();
  std::vector<MBEntityHandle> tris;
    
  const MBCartVect point( ray_point );
  const MBCartVect dir( unit_ray_dir );
  
  for (MBRange::iterator b = boxes.begin(); b != boxes.end(); ++b)
  {
#ifndef MB_OBB_USE_VECTOR_QUERIES
    MBRange tris;
# ifdef MB_OBB_USE_TYPE_QUERIES
    rval = instance->get_entities_by_type( *b, MBTRI, tris );
# else
    rval = instance->get_entities_by_handle( *b, tris );
# endif
#else
    tris.clear();
    rval = instance->get_entities_by_handle( *b, tris );
#endif
    if (MB_SUCCESS != rval)
      return rval;
//dump_fragmentation( tris );
    
#ifndef MB_OBB_USE_VECTOR_QUERIES
    for (MBRange::iterator t = tris.begin(); t != tris.end(); ++t)
#else
    for (std::vector<MBEntityHandle>::iterator t = tris.begin(); t != tris.end(); ++t)
#endif
    {
#ifndef MB_OBB_USE_TYPE_QUERIES
      if (instance->type_from_handle(*t) != MBTRI)
        continue;
#endif
    
      const MBEntityHandle* conn;
      int len;
      rval = instance->get_connectivity( *t, conn, len, true );
      if (MB_SUCCESS != rval)
        return rval;
      
      MBCartVect coords[3];
      rval = instance->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval)
        return rval;
      
      double t;
      if (MBGeomUtil::ray_tri_intersect( coords, point, dir, tolerance, t, ray_length ))
        intersection_distances_out.push_back(t);
    }
  }
  
  return MB_SUCCESS;
}                    

MBErrorCode MBOrientedBoxTreeTool::ray_intersect_triangles( 
                          std::vector<double>& intersection_distances_out,
                          MBEntityHandle root_set,
                          double tolerance,
                          const double ray_point[3],
                          const double unit_ray_dir[3],
                          const double* ray_length )
{
  MBRange boxes;
  MBErrorCode rval;
  
  rval = ray_intersect_boxes( boxes, root_set, tolerance, ray_point, unit_ray_dir, ray_length );
  if (MB_SUCCESS != rval)
    return rval;
    
  return ray_intersect_triangles( intersection_distances_out, boxes, tolerance, ray_point, unit_ray_dir, ray_length );
}

MBErrorCode MBOrientedBoxTreeTool::ray_intersect_boxes( 
                          MBRange& boxes_out,
                          MBEntityHandle root_set,
                          double tolerance,
                          const double ray_point[3],
                          const double unit_ray_dir[3],
                          const double* ray_length )
{
  RayIntersector op( this, ray_point, unit_ray_dir, ray_length, tolerance, boxes_out );
  return preorder_traverse( root_set, op );
}

MBErrorCode RayIntersector::visit( MBEntityHandle node,
                                   int ,
                                   bool& descend ) 
{
  MBOrientedBox box;
  MBErrorCode rval = tool->box( node, box );
  if (MB_SUCCESS != rval)
    return rval;
  
  descend = box.intersect_ray( b, m, tol, len);
  return MB_SUCCESS;
}


MBErrorCode RayIntersector::leaf( MBEntityHandle node )
{
  boxes.insert(node);
  return MB_SUCCESS;
}


/********************** Ray/Set Intersection ****************************/


class RayIntersectSets : public MBOrientedBoxTreeTool::Op
{
  private:
    MBOrientedBoxTreeTool* tool;
    unsigned minTolInt;
    const MBCartVect b, m;
    const double* len;
    const double tol;
    std::vector<double>& intersections;
    std::vector<MBEntityHandle>& sets;
    
    MBEntityHandle lastSet;
    int lastSetDepth;
    
    void add_intersection( double t );
    
  public:
    RayIntersectSets( MBOrientedBoxTreeTool* tool_ptr,
                   const double* ray_point,
                   const double* unit_ray_dir,
                   const double *ray_length,
                   double tolerance,
                   unsigned min_tol_intersections,
                   std::vector<double>& intersections,
                   std::vector<MBEntityHandle>& surfaces )
      : tool(tool_ptr), minTolInt(min_tol_intersections),
        b(ray_point), m(unit_ray_dir),
        len(ray_length), tol(tolerance),
        intersections(intersections),
        sets(surfaces), lastSet(0)
      { }
  
    virtual MBErrorCode visit( MBEntityHandle node,
                               int depth,
                               bool& descend );
    virtual MBErrorCode leaf( MBEntityHandle node );
};

MBErrorCode RayIntersectSets::visit( MBEntityHandle node,
                                     int depth,
                                     bool& descend ) 
{
  MBOrientedBox box;
  MBErrorCode rval = tool->box( node, box );
  if (MB_SUCCESS != rval)
    return rval;
  
  descend = box.intersect_ray( b, m, tol, len);
  
  if (lastSet && depth <= lastSetDepth)
    lastSet = 0;

  if (descend && !lastSet) {
    MBRange sets;
    rval = tool->get_moab_instance()->get_entities_by_type( node, MBENTITYSET, sets );
    if (MB_SUCCESS != rval)
      return rval;

    if (!sets.empty()) {
      if (sets.size() > 1)
        return MB_FAILURE;
      lastSet = *sets.begin();
      lastSetDepth = depth;
    }
  }
    
  return MB_SUCCESS;
}

MBErrorCode RayIntersectSets::leaf( MBEntityHandle node )
{
  if (!lastSet) // if no surface has been visited yet, something's messed up.
    return MB_FAILURE;
  
#ifndef MB_OBB_USE_VECTOR_QUERIES
  MBRange tris;
# ifdef MB_OBB_USE_TYPE_QUERIES
  MBErrorCode rval = tool->get_moab_instance()->get_entities_by_type( node, MBTRI, tris );
# else
  MBErrorCode rval = tool->get_moab_instance()->get_entities_by_handle( node, tris );
# endif
#else
  std::vector<MBEntityHandle> tris;
  MBErrorCode rval = tool->get_moab_instance()->get_entities_by_handle( node, tris );
#endif
  if (MB_SUCCESS != rval)
    return rval;

#ifndef MB_OBB_USE_VECTOR_QUERIES
  for (MBRange::iterator t = tris.begin(); t != tris.end(); ++t)
#else
  for (std::vector<MBEntityHandle>::iterator t = tris.begin(); t != tris.end(); ++t)
#endif
  {
#ifndef MB_OBB_USE_TYPE_QUERIES
    if (tool->get_moab_instance()->type_from_handle(*t) != MBTRI)
      continue;
#endif
    
    const MBEntityHandle* conn;
    int num_conn;
    rval = tool->get_moab_instance()->get_connectivity( *t, conn, num_conn, true );
    if (MB_SUCCESS != rval)
      return rval;

    MBCartVect coords[3];
    rval = tool->get_moab_instance()->get_coords( conn, 3, coords[0].array() );
    if (MB_SUCCESS != rval)
      return rval;

    double t;
    if (MBGeomUtil::ray_tri_intersect( coords, b, m, tol, t, len )) 
        // NOTE: add_intersection may modify the 'len' member, which
        //       will affect subsequent calls to ray_tri_intersect in 
        //       this loop.
      add_intersection( t );
  }
  return MB_SUCCESS;
}

void RayIntersectSets::add_intersection( double t )
{
    // Check if the 'len' pointer is pointing into the intersection
    // list.  If this is the case, then the list contains, at that
    // location, an intersection greater than the tolerance away from
    // the base point of the ray.
  int len_idx = -1;
  if (len && len >= &intersections[0] && len < &intersections[0] + intersections.size())
    len_idx = len - &intersections[0];

    // If the intersection is within tol of the ray base point, we 
    // always add it to the list.
  if (t <= tol) {
      // If the list contains an intersection outside the tolerance...
    if (len_idx >= 0) {
        // If we no longer want an intersection outside the tolerance,
        // remove it.
      if (intersections.size() >= minTolInt) {
        intersections[len_idx] = t;
        sets[len_idx] = lastSet;
          // From now on, we want only intersections within the tolerance,
          // so update length accordingly
        len = &tol;
      }
        // Otherwise appended to the list and update pointer
      else {
        intersections.push_back(t);
        sets.push_back(lastSet);
        len = &intersections[len_idx];
      }
    }
      // Otherwise just append it
    else {
      intersections.push_back(t);
      sets.push_back(lastSet);
        // If we have all the intersections we want, set
        // length such that we will only find further intersections
        // within the tolerance
      if (intersections.size() >= minTolInt)
        len = &tol;
    }
  }
    // Otherwise the intersection is outside the tolerance
    // If we already have an intersection outside the tolerance and
    // this one is closer, replace the existing one with this one.
  else if (len_idx >= 0) {
    if (t <= *len) {
      intersections[len_idx] = t;
      sets[len_idx] = lastSet;
    }
  }
    // Otherwise if we want an intersection outside the tolerance
    // and don'thave one yet, add it.
  else if (intersections.size() < minTolInt) {
    intersections.push_back( t );
    sets.push_back( lastSet );
      // udpate length.  this is currently the closest intersection, so
      // only want further intersections that are closer than this one.
    len = &intersections.back();
  }
}
  

MBErrorCode MBOrientedBoxTreeTool::ray_intersect_sets( 
                                    std::vector<double>& distances_out,
                                    std::vector<MBEntityHandle>& sets_out,
                                    MBEntityHandle root_set,
                                    double tolerance,
                                    unsigned min_tolerace_intersections,
                                    const double ray_point[3],
                                    const double unit_ray_dir[3],
                                    const double* ray_length )
{
  RayIntersectSets op( this, ray_point, unit_ray_dir, ray_length, tolerance, 
                       min_tolerace_intersections, distances_out, sets_out );
  return preorder_traverse( root_set, op );
}



/********************** Closest Point code ***************/

MBErrorCode MBOrientedBoxTreeTool::closest_to_location( 
                                     const double* point,
                                     MBEntityHandle root,
                                     double tolerance,
                                     double* point_out,
                                     MBEntityHandle& facet_out,
                                     MBEntityHandle* set_out ) 
{
  std::vector<MBEntityHandle> facets(1), sets;
  MBErrorCode rval = closest_to_location( point, root, tolerance, facets, set_out ? &sets : 0 );
  if (MB_SUCCESS != rval)
    return rval;
  if (facets.empty())
    return MB_ENTITY_NOT_FOUND;
  facet_out = facets.front();
  if (set_out)
    *set_out = sets.front();
  
  const MBEntityHandle* conn;
  int len;
  rval = instance->get_connectivity( facet_out, conn, len, true );
  if (MB_SUCCESS != rval)
    return rval;
  
  std::vector<MBCartVect> coords( len );
  rval = instance->get_coords( conn, len, coords[0].array() );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (len == 3)
    MBGeomUtil::closest_location_on_tri( *reinterpret_cast<const MBCartVect*>(point), 
                             &coords[0], 
                             *reinterpret_cast<MBCartVect*>(point_out) );
  else
    MBGeomUtil::closest_location_on_polygon( *reinterpret_cast<const MBCartVect*>(point), 
                                 &coords[0], len,
                                 *reinterpret_cast<MBCartVect*>(point_out) );
  return MB_SUCCESS;
}

struct MBOBBTreeCPFrame {
  MBOBBTreeCPFrame( double d, MBEntityHandle n, MBEntityHandle s )
    : dist_sqr(d), node(n), mset(s) {}
  double dist_sqr;
  MBEntityHandle node;
  MBEntityHandle mset;
};
                                     
MBErrorCode MBOrientedBoxTreeTool::closest_to_location( const double* point,
                                     MBEntityHandle root,
                                     double tolerance,
                                     std::vector<MBEntityHandle>& facets_out,
                                     std::vector<MBEntityHandle>* sets_out )
{
  MBErrorCode rval;
  const MBCartVect loc( point );
  double smallest_dist_sqr = std::numeric_limits<double>::max();
  double smallest_dist = smallest_dist_sqr;
  
  MBEntityHandle current_set = 0;
  MBRange sets;
  std::vector<MBEntityHandle> children(2);
  std::vector<double> coords;
  std::vector<MBOBBTreeCPFrame> stack;
    
  stack.push_back( MBOBBTreeCPFrame(0.0, root, current_set) );
  
  while( !stack.empty() ) {

      // pop from top of stack
    MBEntityHandle node = stack.back().node;
    double dist_sqr = stack.back().dist_sqr;
    current_set = stack.back().mset;
    stack.pop_back();

      // If current best result is closer than the box, skip this tree node.
    if (dist_sqr > smallest_dist_sqr + tolerance)
      continue;

      // Check if this node has a set associated with it
    if (sets_out && !current_set) {
      sets.clear();
      rval = instance->get_entities_by_type( node, MBENTITYSET, sets );
      if (MB_SUCCESS != rval)
        return rval;
      if (!sets.empty()) {
        if (sets.size() != 1)
          return MB_MULTIPLE_ENTITIES_FOUND;
        current_set = *sets.begin();
      }
    }

      // Get child boxes
    children.clear();
    rval = instance->get_child_meshsets( node, children );
    if (MB_SUCCESS != rval)
      return rval;

      // if not a leaf node
    if (!children.empty()) {
      if (children.size() != 2)
        return MB_MULTIPLE_ENTITIES_FOUND;
    
        // get boxes
      MBOrientedBox box1, box2;
      rval = box( children[0], box1 );
      if (MB_SUCCESS != rval) return rval;
      rval = box( children[1], box2 );
      if (MB_SUCCESS != rval) return rval;
      
        // get distance from each box
      MBCartVect pt1, pt2;
      box1.closest_location_in_box( loc, pt1 );
      box2.closest_location_in_box( loc, pt2 );
      pt1 -= loc;
      pt2 -= loc;
      const double dsqr1 = pt1 % pt1;
      const double dsqr2 = pt2 % pt2;
      
        // push children on tree such that closer one is on top
      if (dsqr1 < dsqr2) {
        stack.push_back( MBOBBTreeCPFrame(dsqr2, children[1], current_set ) );
        stack.push_back( MBOBBTreeCPFrame(dsqr1, children[0], current_set ) );
      }
      else {
        stack.push_back( MBOBBTreeCPFrame(dsqr1, children[0], current_set ) );
        stack.push_back( MBOBBTreeCPFrame(dsqr2, children[1], current_set ) );
      }
    }
    else { // LEAF NODE
      MBRange facets;
      rval = instance->get_entities_by_dimension( node, 2, facets );
      if (MB_SUCCESS != rval)
        return rval;
      
      const MBEntityHandle* conn;
      int len;
      MBCartVect tmp, diff;
      for (MBRange::iterator i = facets.begin(); i != facets.end(); ++i) {
        rval = instance->get_connectivity( *i, conn, len, true );
        if (MB_SUCCESS != rval)
          return rval;
        
        coords.resize(3*len);
        rval = instance->get_coords( conn, len, &coords[0] );
        if (MB_SUCCESS != rval)
          return rval;
        
        if (len == 3) 
          MBGeomUtil::closest_location_on_tri( loc, (MBCartVect*)(&coords[0]), tmp );
        else
          MBGeomUtil::closest_location_on_polygon( loc, (MBCartVect*)(&coords[0]), len, tmp );
        
        diff = tmp - loc;
        dist_sqr = diff % diff;
        if (dist_sqr < smallest_dist_sqr) {
          if (dist_sqr < smallest_dist_sqr + tolerance*(tolerance - 2*smallest_dist)) {
            facets_out.clear();
            if (sets_out)
              sets_out->clear();
          }
          smallest_dist_sqr = dist_sqr;
          smallest_dist = sqrt(smallest_dist_sqr);
            // put closest result at start of lists
          facets_out.push_back( *i );
          std::swap( facets_out.front(), facets_out.back() );
          if (sets_out) {
            sets_out->push_back( current_set );
            std::swap( sets_out->front(), sets_out->back() );
          }
        }
        else if (dist_sqr <= smallest_dist_sqr + tolerance*(tolerance + 2*smallest_dist)) {
          facets_out.push_back(*i);
          if (sets_out)
            sets_out->push_back( current_set );
        }
      }
    } // LEAF NODE
  }
  
  return MB_SUCCESS;
}
    
    

/********************** Tree Printing Code ****************************/


class TreeLayoutPrinter : public MBOrientedBoxTreeTool::Op
{
  public:
    TreeLayoutPrinter( std::ostream& stream,
                       MBInterface* instance );
    
    virtual MBErrorCode visit( MBEntityHandle node, 
                               int depth,
                               bool& descend );
    virtual MBErrorCode leaf( MBEntityHandle node );
  private:

    MBInterface* instance;
    std::ostream& outputStream;
    std::vector<bool> path;
};

TreeLayoutPrinter::TreeLayoutPrinter( std::ostream& stream,
                                      MBInterface* interface )
  : instance(interface),
    outputStream(stream)
  {}

MBErrorCode TreeLayoutPrinter::visit( MBEntityHandle node, 
                                      int depth,
                                      bool& descend )
{
  descend = true;
  
  if ((unsigned)depth > path.size()) {
    //assert(depth+1 == path.size); // preorder traversal
    path.push_back(true);
  }
  else {
    path.resize( depth );
    if (depth)
      path.back() = false;
  }
  
  for (unsigned i = 0; i+1 < path.size(); ++i) {
    if (path[i])
      outputStream << "|   ";
    else
      outputStream << "    ";
  }
  if (depth) {
    if (path.back())
      outputStream << "+---";
    else
      outputStream << "\\---";
  }
  outputStream << instance->id_from_handle( node ) << std::endl;
  return MB_SUCCESS;
}

MBErrorCode TreeLayoutPrinter::leaf( MBEntityHandle ) { return MB_SUCCESS; }
    

class TreeNodePrinter : public MBOrientedBoxTreeTool::Op
{
  public:
    TreeNodePrinter( std::ostream& stream,
                     bool list_contents,
                     bool list_box,
                     const char* id_tag_name,
                     MBOrientedBoxTreeTool* tool_ptr );
    
    virtual MBErrorCode visit( MBEntityHandle node, 
                               int depth,
                               bool& descend );
                               
    virtual MBErrorCode leaf( MBEntityHandle ) { return MB_SUCCESS; }
  private:
  
    MBErrorCode print_geometry( MBEntityHandle node );
    MBErrorCode print_contents( MBEntityHandle node );
    MBErrorCode print_counts( MBEntityHandle node );
  
    bool printContents;
    bool printGeometry;
    bool haveTag;
    MBTag tag, gidTag, geomTag;
    MBInterface* instance;
    MBOrientedBoxTreeTool* tool;
    std::ostream& outputStream;
};


TreeNodePrinter::TreeNodePrinter( std::ostream& stream,
                                  bool list_contents,
                                  bool list_box,
                                  const char* id_tag_name,
                                  MBOrientedBoxTreeTool* tool_ptr )
  : printContents( list_contents ),
    printGeometry( list_box ),
    haveTag( false ),
    tag( 0 ),
    gidTag(0), geomTag(0),
    instance( tool_ptr->get_moab_instance() ),
    tool(tool_ptr),
    outputStream( stream )
{
  MBErrorCode rval;
  if (id_tag_name) {
    rval = instance->tag_get_handle( id_tag_name, tag );
    if (!rval) {
      std::cerr << "Could not get tag \"" << id_tag_name << "\"\n";
      stream << "Could not get tag \"" << id_tag_name << "\"\n";
    }
    else {
      int size;
      rval = instance->tag_get_size( tag, size );
      if (!rval) {
        std::cerr << "Could not get size for tag \"" << id_tag_name << "\"\n";
        stream << "Could not get size for tag \"" << id_tag_name << "\"\n";
      }
      else if (size != sizeof(int)) {
        std::cerr << "Tag \"" << id_tag_name << "\" is not an integer\n";
        stream << "Tag \"" << id_tag_name << "\" is not an integer\n";
      }
      else {
        haveTag = true;
      }
    }
  }
  
  rval = instance->tag_get_handle( GLOBAL_ID_TAG_NAME, gidTag );
  if (MB_SUCCESS != rval)
    gidTag = 0;
  rval = instance->tag_get_handle( GEOM_DIMENSION_TAG_NAME, geomTag );
  if (MB_SUCCESS != rval)
    geomTag = 0;
}   

MBErrorCode TreeNodePrinter::visit( MBEntityHandle node, int, bool& descend )
{
  descend = true;
  MBEntityHandle setid = instance->id_from_handle( node );
  outputStream << setid << ":" << std::endl;
  
  MBRange surfs;
  MBErrorCode r3 = MB_SUCCESS;
  if (geomTag) {
    const int two = 2;
    const void* tagdata[] = {&two};
    r3 = instance->get_entities_by_type_and_tag( node, MBENTITYSET, &geomTag, tagdata, 1, surfs );
  
    if (MB_SUCCESS == r3 && surfs.size() == 1) {
      MBEntityHandle surf = *surfs.begin();
      int id;
      if (gidTag && MB_SUCCESS == instance->tag_get_data( gidTag, &surf, 1, &id ))
        outputStream << "  Surface " << id << std::endl;
      else
        outputStream << "  Surface w/ unknown ID (" << surf << ")" << std::endl;
    }
  }
  
  MBErrorCode r1 = printGeometry ? print_geometry( node ) : MB_SUCCESS;
  MBErrorCode r2 = printContents ? print_contents( node ) : print_counts( node );
  outputStream << std::endl;
  
  if (MB_SUCCESS != r1)
    return r1;
  else if (MB_SUCCESS != r2)
    return r2;
  else
    return r3;
}

MBErrorCode TreeNodePrinter::print_geometry( MBEntityHandle node )
{
  MBOrientedBox box;
  MBErrorCode rval= tool->box( node, box );
  if (MB_SUCCESS != rval)
    return rval;
  
  MBCartVect length = box.dimensions();
  
  outputStream << box.center << "  Radius: " 
               << box.inner_radius() << " - " << box.outer_radius() << std::endl
               << '+' << box.axis[0] << " : " << length[0] << std::endl
               << 'x' << box.axis[1] << " : " << length[1] << std::endl
               << 'x' << box.axis[2] << " : " << length[2] << std::endl;
  return MB_SUCCESS;
}

MBErrorCode TreeNodePrinter::print_counts( MBEntityHandle node )
{
  for (MBEntityType type = MBVERTEX; type != MBMAXTYPE; ++type) {
    int count = 0;
    MBErrorCode rval = instance->get_number_entities_by_type( node, type, count );
    if (MB_SUCCESS != rval)
      return rval;
    if(count > 0)
      outputStream << " " << count << " " << MBCN::EntityTypeName(type) << std::endl;
  }
  return MB_SUCCESS;
}

MBErrorCode TreeNodePrinter::print_contents( MBEntityHandle node )
{
    // list contents
  for (MBEntityType type = MBVERTEX; type != MBMAXTYPE; ++type) {
    MBRange range;
    MBErrorCode rval = instance->get_entities_by_type( node, type, range );
    if (MB_SUCCESS != rval)
      return rval;
    if (range.empty())
      continue;
    outputStream << " " << MBCN::EntityTypeName(type) << " ";  
    std::vector<int> ids( range.size() );
    if (haveTag) {
      rval = instance->tag_get_data( tag, range, &ids[0] );
      std::sort( ids.begin(), ids.end() );
    }
    else {
      MBRange::iterator ri = range.begin();
      std::vector<int>::iterator vi = ids.begin();
      while (ri != range.end()) {
        *vi = instance->id_from_handle( *ri );
        ++ri;
        ++vi;
      }
    }
    
    unsigned i = 0;
    for(;;) {
      unsigned beg = i, end;
      do { end = i++; } while (i < ids.size() && ids[end]+1 == ids[i]);
      if (end == beg)
        outputStream << ids[end];
      else if (end == beg+1) 
        outputStream << ids[beg] << ", " << ids[end];
      else
        outputStream << ids[beg] << "-" << ids[end];
        
      if (i == ids.size()) {
        outputStream << std::endl;
        break;
      }
      else 
        outputStream << ", ";
    }
  }
  
  return MB_SUCCESS;
}

  
void MBOrientedBoxTreeTool::print( MBEntityHandle set, std::ostream& str, bool list, const char* tag )
{
  TreeLayoutPrinter op1( str, instance );
  TreeNodePrinter op2( str, list, true, tag, this );
  MBErrorCode r1 = preorder_traverse( set, op1 );
  str << std::endl;
  MBErrorCode r2 = preorder_traverse( set, op2 );
  if (r1 != MB_SUCCESS || r2 != MB_SUCCESS) {
    std::cerr << "Errors encountered while printing tree\n";
    str << "Errors encountered while printing tree\n";
  }
}


/********************** Tree Statistics Code ****************************/


struct StatData {
  struct Ratio {
    double min, max, sum, sqr;
    int hist[10];
    Ratio() 
      : min(std::numeric_limits<double>::max()), 
        max(-std::numeric_limits<double>::max()), 
        sum(0.0), sqr(0.0)
      { hist[0] = hist[1] = hist[2] = hist[3] = hist[4] = hist[5] =
        hist[6] = hist[7] = hist[8] = hist[9] = 0; }
    void accum( double v ) {
      if (v < min) min = v;
      if (v > max) max = v;
      sum += v;
      sqr += v*v;
      int i = (int)(10*v);
      if (i < 0) i = 0;
      else if (i > 9) i = 9;
      ++hist[i];
    }
  };
  
  template <typename T> struct Stat {
    T min, max;
    double sum, sqr;
    Stat() : sum(0.0), sqr(0.0) {
      std::numeric_limits<T> lim;
      min = lim.max();
      if (lim.is_integer)
        max = lim.min();
      else
        max = -lim.max();
    }
    void accum( T v ) {
      if (v < min) min = v;
      if (v > max) max = v;
      sum += v;
      sqr += (double)v * v;
    }
  };

  StatData() :
    count(0) 
    {}

  Ratio volume;
  Ratio entities;
  Ratio radius;
  Stat<unsigned> leaf_ent;
  Stat<double> vol;
  Stat<double> area;
  std::vector<unsigned> leaf_depth;
  unsigned count;
};

static int measure( const MBCartVect& v, double& result )
{
  const double tol = 1e-6;
  int dims = 0;
  result = 1;
  for (int i = 0; i < 3; ++i)
    if (v[i] > tol) {
      ++dims; 
      result *= v[i];
  }
  return dims;
}
  

static MBErrorCode recursive_stats( MBOrientedBoxTreeTool* tool,
                                    MBInterface* instance,
                                    MBEntityHandle set,
                                    int depth,
                                    StatData& data,
                                    unsigned& count_out,
                                    MBCartVect& dimensions_out )
{
  MBErrorCode rval;
  MBOrientedBox box;
  std::vector<MBEntityHandle> children(2);
  unsigned counts[2];
  bool isleaf;
  
  ++data.count;
  
  rval = tool->box( set, box );
  if (MB_SUCCESS != rval) return rval;
  children.clear();
  rval = instance->get_child_meshsets( set, children );
  if (MB_SUCCESS != rval) return rval;
  isleaf = children.empty();
  if (!isleaf && children.size() != 2)
    return MB_MULTIPLE_ENTITIES_FOUND;
  
  dimensions_out = box.dimensions();
  data.radius.accum( box.inner_radius() / box.outer_radius());
  data.vol.accum( box.volume() );
  data.area.accum( box.area() );
  
  if (isleaf) {
    if (data.leaf_depth.size() <= (unsigned)depth)
      data.leaf_depth.resize( depth+1, 0 );
    ++data.leaf_depth[depth];
    
    int count;
    rval = instance->get_number_entities_by_handle( set, count );
    if (MB_SUCCESS != rval) return rval;
    count_out = count;
    data.leaf_ent.accum( count_out );
  }
  else {
    for (int i = 0; i < 2; ++i) {
      MBCartVect dims;
      rval = recursive_stats( tool, instance, children[i], depth+1, data, counts[i], dims );
      if (MB_SUCCESS != rval) return rval;
      double this_measure, chld_measure;
      int this_dim = measure( dimensions_out, this_measure );
      int chld_dim = measure( dims, chld_measure );
      double ratio;
      if (chld_dim < this_dim)
        ratio = 0;
      else
        ratio = chld_measure / this_measure;
  
      data.volume.accum( ratio );
    }
    count_out = counts[0] + counts[1];
    data.entities.accum( (double)counts[0] / count_out );
    data.entities.accum( (double)counts[1] / count_out );
  }
  return MB_SUCCESS;
}

static inline double std_dev( double sqr, double sum, double count )
{
  sum /= count;
  sqr /= count;
  return sqrt( sqr - sum*sum );
}

//#define WW <<std::setw(10)<<std::fixed<<
#define WE <<std::setw(10)<<
#define WW WE
MBErrorCode MBOrientedBoxTreeTool::stats( MBEntityHandle set, std::ostream& s )
{
  StatData d;
  MBErrorCode rval;
  unsigned total_entities, i;
  MBCartVect total_dim;
  
  rval = recursive_stats( this, instance, set, 0, d, total_entities, total_dim );
  if (MB_SUCCESS != rval)
    return rval;
  
  unsigned tree_height = d.leaf_depth.size();
  unsigned min_leaf_depth = tree_height, num_leaves = 0;
  unsigned max_leaf_per_depth = 0;
  double sum_leaf_depth = 0, sqr_leaf_depth = 0;
  for (i = 0; i < d.leaf_depth.size(); ++i) {
    unsigned val = d.leaf_depth[i];
    num_leaves += val;
    sum_leaf_depth += (double)val*i;
    sqr_leaf_depth += (double)val*i*i;
    if (val && i < min_leaf_depth)
      min_leaf_depth = i;
    if (max_leaf_per_depth < val)
      max_leaf_per_depth = val;
  }
  unsigned num_non_leaf = d.count - num_leaves;
  
  double rv = total_dim[0]*total_dim[1]*total_dim[2];
  s << "entities in tree:  " << total_entities << std::endl
    << "root volume:       " << rv << std::endl
    << "total node volume: " << d.vol.sum << std::endl
    << "total/root volume: " << d.vol.sum/rv << std::endl
    << "tree height:       " << tree_height << std::endl
    << "node count:        " << d.count << std::endl
    << "leaf count:        " << num_leaves << std::endl
    << std::endl;
  
  double avg_leaf_depth = sum_leaf_depth / num_leaves;
  double rms_leaf_depth = sqrt( sqr_leaf_depth / num_leaves );
  double std_leaf_depth = std_dev( sqr_leaf_depth, sum_leaf_depth, num_leaves );

  double avg_leaf_ent = d.leaf_ent.sum / num_leaves;
  double rms_leaf_ent = sqrt( d.leaf_ent.sqr / num_leaves );
  double std_leaf_ent = std_dev( d.leaf_ent.sqr, d.leaf_ent.sum, num_leaves );

  unsigned num_child = 2 * num_non_leaf;

  double avg_vol_ratio = d.volume.sum / num_child;
  double rms_vol_ratio = sqrt( d.volume.sqr / num_child );
  double std_vol_ratio = std_dev( d.volume.sqr, d.volume.sum, num_child);

  double avg_ent_ratio = d.entities.sum / num_child;
  double rms_ent_ratio = sqrt( d.entities.sqr / num_child );
  double std_ent_ratio = std_dev( d.entities.sqr, d.entities.sum, num_child);

  double avg_rad_ratio = d.radius.sum / d.count;
  double rms_rad_ratio = sqrt( d.radius.sqr / d.count );
  double std_rad_ratio = std_dev( d.radius.sqr, d.radius.sum, d.count );
  
  double avg_vol = d.vol.sum / d.count;
  double rms_vol = sqrt( d.vol.sqr / d.count );
  double std_vol = std_dev( d.vol.sqr, d.vol.sum, d.count );
  
  double avg_area = d.area.sum / d.count;
  double rms_area = sqrt( d.area.sqr / d.count );
  double std_area = std_dev( d.area.sqr, d.area.sum, d.count );
      
  int prec = s.precision();
  s <<                         "                   " WW "Minimum"      WW "Average"      WW "RMS"          WW "Maximum"             WW "Std.Dev."     << std::endl;
  s << std::setprecision(1) << "Leaf Depth         " WW min_leaf_depth WW avg_leaf_depth WW rms_leaf_depth WW d.leaf_depth.size()-1 WW std_leaf_depth << std::endl; 
  s << std::setprecision(0) << "Entities/Leaf      " WW d.leaf_ent.min WW avg_leaf_ent   WW rms_leaf_ent   WW d.leaf_ent.max        WW std_leaf_ent   << std::endl;
  s << std::setprecision(3) << "Child Volume Ratio " WW d.volume.min   WW avg_vol_ratio  WW rms_vol_ratio  WW d.volume.max          WW std_vol_ratio  << std::endl;
  s << std::setprecision(3) << "Child Entity Ratio " WW d.entities.min WW avg_ent_ratio  WW rms_ent_ratio  WW d.entities.max        WW std_ent_ratio  << std::endl;
  s << std::setprecision(3) << "Box Radius Ratio   " WW d.radius.min   WW avg_rad_ratio  WW rms_rad_ratio  WW d.radius.max          WW std_rad_ratio  << std::endl;
  s << std::setprecision(0) << "Box volume         " WE d.vol.min      WE avg_vol        WE rms_vol        WE d.vol.max             WE std_vol        << std::endl;
  s << std::setprecision(0) << "Largest side area  " WE d.area.min     WE avg_area       WE rms_area       WE d.area.max            WE std_area       << std::endl;
  s << std::setprecision(prec) << std::endl;
  
  s << "Leaf Depth Histogram (Root depth is 0)" << std::endl;
  double f = 60.0 / max_leaf_per_depth;
  for (i = min_leaf_depth; i < d.leaf_depth.size(); ++i)
    s << std::setw(2) << i << " " << std::setw(5) << d.leaf_depth[i] << " |"
      << std::setfill('*') << std::setw((int)round(f*d.leaf_depth[i])) << "" 
      << std::setfill(' ') << std::endl;
  s <<std::endl;
  
  s << "Child/Parent Volume Ratio Histogram" << std::endl;
  f = 60.0 / *(std::max_element(d.volume.hist, d.volume.hist+10));
  for (i = 0; i < 10u; ++i)
    s << "0." << i << " " << std::setw(5) << d.volume.hist[i] << " |"
      << std::setfill('*') << std::setw((int)round(f*d.volume.hist[i])) << ""
      << std::setfill(' ') << std::endl;
  s <<std::endl;
  
  s << "Child/Parent Entity Count Ratio Histogram" << std::endl;
  f = 60.0 / *(std::max_element(d.entities.hist, d.entities.hist+10));
  for (i = 0; i < 10u; ++i)
    s << "0." << i << " " << std::setw(5) << d.entities.hist[i] << " |"
      << std::setfill('*') << std::setw((int)round(f*d.entities.hist[i])) << ""
      << std::setfill(' ') << std::endl;
  s <<std::endl;
  
  s << "Inner/Outer Radius Ratio Histogram (~0.70 for cube)" << std::endl;
    // max radius ratio for a box is about 0.7071.  Move any boxes
    // in the .7 bucket into .6 and print .0 to .6.
  d.radius.hist[6] += d.radius.hist[7]; 
  f = 60.0 / *(std::max_element(d.entities.hist, d.entities.hist+7));
  for (i = 0; i < 7u; ++i)
    s << "0." << i << " " << std::setw(5) << d.entities.hist[i] << " |"
      << std::setfill('*') << std::setw((int)round(f*d.entities.hist[i])) << ""
      << std::setfill(' ') << std::endl;
  s <<std::endl;
  
  return MB_SUCCESS;
}
