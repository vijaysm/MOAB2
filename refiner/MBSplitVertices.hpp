/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2007 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**\class MBSplitVertices
  *\brief A dictionary of new vertices.
  *
  * An array of existing vertex handles used as a key in a dictionary of new vertices.
  */
#ifndef MB_SPLITVERTICES_HPP
#define MB_SPLITVERTICES_HPP

#include "MBTypes.h"
#include "MBProcessSet.hpp"
#include "MBTagConventions.hpp"

#include <map>
#include <vector>
#include <algorithm>

class MBRefinerTagManager;

template< int _n >
class MBSplitVertexIndex
{
public:
  MBSplitVertexIndex() { }
  MBSplitVertexIndex( const int* src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src[i]; std::sort( this->handles, this->handles + _n ); }
  MBSplitVertexIndex( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; this->process_set = src.process_set; }
  MBSplitVertexIndex& operator = ( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; this->process_set = src.process_set; return *this; }

  void set_common_processes( const MBProcessSet& procs )
    { this->process_set = procs; }
  MBProcessSet& common_processes()
    { return this->process_set; }
  const MBProcessSet& common_processes() const
    { return this->process_set; }

  bool operator < ( const MBSplitVertexIndex<_n>& other ) const
    {
    // Ignore the process set. Only program errors lead to mismatched process sets with identical handles.
    for ( int i = 0; i < _n; ++ i )
      if ( this->handles[i] < other.handles[i] )
        return true;
      else if ( this->handles[i] > other.handles[i] )
        return false;
    return false;
    }

  int handles[_n + 1];
  MBProcessSet process_set;
};

template< int _n >
std::ostream& operator << ( std::ostream& os, const MBSplitVertexIndex<_n>& idx )
{
  for ( int i = 0; i < _n; ++ i )
    {
    os << idx.handles[i] << " ";
    }
  os << "(" << idx.process_set << ")";
  return os;
}

/** A non-templated base class that the template subclasses all share.
  *
  * All methods that need to be accessed by other classes should be
  * declared by the base class so that no knowledge of template parameters
  * is required.
  */
class MBSplitVerticesBase
{
public:
  MBSplitVerticesBase( MBRefinerTagManager* tag_mgr );
  virtual ~MBSplitVerticesBase();

  virtual bool find_or_create(
    const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
    std::map<MBProcessSet,int>& proc_partition_counts ) = 0;

  virtual bool create_element(
    const MBEntityHandle* split_src, MBEntityHandle& elem_handle,
    std::map<MBProcessSet,int>& proc_partition_counts ) = 0;

  virtual void assign_global_ids( std::map<MBProcessSet,int>& gids ) = 0;

  /// Determine which processes will contain an output vertex given the split vertices defining it.
  void update_partition_counts( int num, const MBEntityHandle* split_src, std::map<MBProcessSet,int>& proc_partition_counts );

  /// Prepare to compute the processes on which a new split-vertex will live.
  void begin_vertex_procs();

  /// Call this for each existing corner vertex used to define a split-vertex.
  void add_vertex_procs( MBEntityHandle vert_in );

  /// Call this once after all the add_vertex_procs() calls for a split-vertex to prepare queues for the second stage MPI send. 
  void end_vertex_procs();

  /// Set the tags which indicate sharing process(es) for an entity.
  void set_sharing( MBEntityHandle vert_handle, MBProcessSet& procs );

  MBInterface* mesh_in; // Input mesh. Needed to determine tag values on split_src verts
  MBInterface* mesh_out; // Output mesh. Needed for new vertex set in vert_handle
  MBRefinerTagManager* tag_manager;
  std::vector<int> shared_procs_in; // Used to hold procs sharing an input vert.
  std::vector<int> shared_procs_out; // Used to hold procs sharing an output vert.
  std::vector<int> split_gids; // Used to hold global IDs of split vertices
  MBProcessSet current_shared_procs; // Holds process list as it is being accumulated
  MBProcessSet common_shared_procs; // Holds intersection of several shared_procs_ins.
  int rank; // This process' rank.
  bool first_vertex; // True just after begin_vertex_procs() is called.
  MBTag tag_gid;
};

/** A map from a set of pre-existing entities to a new mesh entity.
  *
  * This is used as a dictionary to determine whether a new vertex should be
  * created on the given n-simplex (n being the template parameter) or whether
  * it has already been created as part of the refinement of a neighboring entity.
  */
template< int _n >
class MBSplitVertices : public std::map<MBSplitVertexIndex<_n>,MBEntityHandle>, public MBSplitVerticesBase
{
public:
  typedef std::map<MBSplitVertexIndex<_n>,MBEntityHandle> MapType;
  typedef typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator MapIteratorType;

  MBSplitVertices( MBRefinerTagManager* tag_mgr );
  virtual ~MBSplitVertices();
  virtual bool find_or_create(
    const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
    std::map<MBProcessSet,int>& proc_partition_counts );
  virtual bool create_element(
    const MBEntityHandle* split_src, MBEntityHandle& elem_handle,
    std::map<MBProcessSet,int>& proc_partition_counts );

  virtual void assign_global_ids( std::map<MBProcessSet,int>& gids );
};

// ------------------------- Template member definitions ----------------------
template< int _n >
MBSplitVertices<_n>::MBSplitVertices( MBRefinerTagManager* tag_mgr )
  : MBSplitVerticesBase( tag_mgr )
{
  this->shared_procs_in.resize( _n * MAX_SHARING_PROCS );
  this->split_gids.resize( _n );
}

template< int _n >
MBSplitVertices<_n>::~MBSplitVertices()
{
}

template< int _n >
bool MBSplitVertices<_n>::find_or_create(
  const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
  std::map<MBProcessSet,int>& proc_partition_counts )
{
  // Get the global IDs of the input vertices
  int stat;
  for ( int i = 0; i < _n; ++ i )
    {
    int gid = -1;
    stat = this->mesh_in->tag_get_data( this->tag_gid, split_src + i, 1, &gid );
    this->split_gids[i] = gid;
    }
  MBSplitVertexIndex<_n> key( &this->split_gids[0] );
  MapIteratorType it = this->find( key );
  if ( it == this->end() )
    {
    this->update_partition_counts( _n, split_src, proc_partition_counts );
    key.set_common_processes( this->common_shared_procs );
    if ( this->mesh_out->create_vertex( coords, vert_handle ) != MB_SUCCESS )
      {
      return false;
      }
    (*this)[key] = vert_handle;
    this->set_sharing( vert_handle, this->common_shared_procs );
    return true;
    }
  vert_handle = it->second;
  return false;
}

template< int _n >
bool MBSplitVertices<_n>::create_element(
  const MBEntityHandle* split_src, MBEntityHandle& elem_handle,
  std::map<MBProcessSet,int>& proc_partition_counts )
{
}

template< int _n >
void MBSplitVertices<_n>::assign_global_ids( std::map<MBProcessSet,int>& gids )
{
  typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator it;
  for ( it = this->begin(); it != this->end(); ++ it )
    {
    int gid = gids[it->first.process_set] ++;
    this->mesh_out->tag_set_data( this->tag_gid, &it->second, 1, &gid );
    std::cout << "Assigning " << it->first << " -> " << gid << "\n";
    }
}

#endif /* MB_SPLITVERTICES_HPP */
