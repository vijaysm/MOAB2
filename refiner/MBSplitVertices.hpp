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
  * An array of existing vertex ids used as a key in a dictionary of new vertices.
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
    { for ( int i = 0; i < _n; ++ i ) this->ids[i] = src[i]; std::sort( this->ids, this->ids + _n ); }
  MBSplitVertexIndex( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->ids[i] = src.ids[i]; this->process_set = src.process_set; }
  MBSplitVertexIndex& operator = ( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->ids[i] = src.ids[i]; this->process_set = src.process_set; return *this; }

  void set_common_processes( const MBProcessSet& procs )
    { this->process_set = procs; }
  MBProcessSet& common_processes()
    { return this->process_set; }
  const MBProcessSet& common_processes() const
    { return this->process_set; }

  bool operator < ( const MBSplitVertexIndex<_n>& other ) const
    {
    // Ignore the process set. Only program errors lead to mismatched process sets with identical ids.
    for ( int i = 0; i < _n; ++ i )
      if ( this->ids[i] < other.ids[i] )
        return true;
      else if ( this->ids[i] > other.ids[i] )
        return false;
    return false;
    }

  int ids[_n + 1];
  MBProcessSet process_set;
};

template< int _n >
std::ostream& operator << ( std::ostream& os, const MBSplitVertexIndex<_n>& idx )
{
  for ( int i = 0; i < _n; ++ i )
    {
    os << idx.ids[i] << " ";
    }
  os << "(" << idx.process_set << ")";
  return os;
}

class MBEntitySourceRecord
{
public:
  MBEntitySourceRecord() { }
  MBEntitySourceRecord( int nc, MBEntityHandle ent, const MBProcessSet& procs )
    { this->ids.resize( nc ); this->handle = ent; this->process_set = procs; }
  MBEntitySourceRecord( const MBEntitySourceRecord& src )
    { this->handle = src.handle; this->process_set = src.process_set; this->ids = src.ids; }
  MBEntitySourceRecord& operator = ( const MBEntitySourceRecord& src )
    { this->handle = src.handle; this->process_set = src.process_set; this->ids = src.ids; return *this; }

  void set_common_processes( const MBProcessSet& procs )
    { this->process_set = procs; }
  MBProcessSet& common_processes()
    { return this->process_set; }
  const MBProcessSet& common_processes() const
    { return this->process_set; }

  bool operator < ( const MBEntitySourceRecord& other ) const
    {
    //assert( this->ids.size() == other.ids.size() );
    std::vector<int>::size_type N = this->ids.size();
    std::vector<int>::size_type i;
    // Ignore the process set. Only program errors lead to mismatched process sets with identical ids.
    for ( i = 0; i < N; ++ i )
      if ( this->ids[i] < other.ids[i] )
        return true;
      else if ( this->ids[i] > other.ids[i] )
        return false;
    return false;
    }

  std::vector<int> ids;
  MBProcessSet process_set;
  MBEntityHandle handle;
};


/** A non-templated base class that the MBSplitVertices template subclasses all share.
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
    std::map<MBProcessSet,int>& proc_partition_counts, bool handles_on_output_mesh ) = 0;

  virtual void assign_global_ids( std::map<MBProcessSet,int>& gids ) = 0;

  MBInterface* mesh_out; // Output mesh. Needed for new vertex set in vert_handle
  MBRefinerTagManager* tag_manager;
  std::vector<int> split_gids; // Used to hold global IDs of split vertices
  MBProcessSet common_shared_procs; // Holds intersection of several shared_procs_ins.
};

/** A vector of pre-existing entities to a new mesh entity.
  *
  * This is used as a dictionary to determine whether a new vertex should be
  * created on the given n-simplex (n being the template parameter) or whether
  * it has already been created as part of the refinement of a neighboring entity.
  */
class MBEntitySource : public std::vector<MBEntitySourceRecord>
{
public:
  typedef std::vector<MBEntitySourceRecord> VecType;
  typedef std::vector<MBEntitySourceRecord>::iterator VecIteratorType;

  MBEntitySource( int num_corners, MBRefinerTagManager* tag_mgr );
  ~MBEntitySource();
  bool create_element(
    MBEntityType etyp, int nconn, const MBEntityHandle* split_src, MBEntityHandle& elem_handle,
    std::map<MBProcessSet,int>& proc_partition_counts );

  void assign_global_ids( std::map<MBProcessSet,int>& gids );

  MBInterface* mesh_out; // Output mesh. Needed for new vertex set in vert_handle
  MBRefinerTagManager* tag_manager;
  MBProcessSet common_shared_procs; // Holds intersection of several shared_procs_ins.
  int num_corners;
};


/** A map from a set of pre-existing vertices to a new mesh vertex.
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
    std::map<MBProcessSet,int>& proc_partition_counts, bool handles_on_output_mesh );

  virtual void assign_global_ids( std::map<MBProcessSet,int>& gids );
};

// ------------------------- Template member definitions ----------------------
template< int _n >
MBSplitVertices<_n>::MBSplitVertices( MBRefinerTagManager* tag_mgr )
  : MBSplitVerticesBase( tag_mgr )
{
  this->split_gids.resize( _n );
}

template< int _n >
MBSplitVertices<_n>::~MBSplitVertices()
{
}

template< int _n >
bool MBSplitVertices<_n>::find_or_create(
  const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
  std::map<MBProcessSet,int>& proc_partition_counts, bool handles_on_output_mesh )
{
  // Get the global IDs of the input vertices
  int stat;
  if ( handles_on_output_mesh )
    {
    stat = this->tag_manager->get_output_gids( _n, split_src, this->split_gids );
    }
  else
    {
    stat = this->tag_manager->get_input_gids( _n, split_src, this->split_gids );
    }
  MBSplitVertexIndex<_n> key( &this->split_gids[0] );
  MapIteratorType it = this->find( key );
  if ( it == this->end() )
    {
    std::cout << " wrt output: " << handles_on_output_mesh << " ";
    this->tag_manager->get_common_processes( _n, split_src, this->common_shared_procs, handles_on_output_mesh );
    proc_partition_counts[this->common_shared_procs]++;
    key.set_common_processes( this->common_shared_procs );
    if ( this->mesh_out->create_vertex( coords + 3, vert_handle ) != MB_SUCCESS )
      {
      return false;
      }
    (*this)[key] = vert_handle;
    this->tag_manager->set_sharing( vert_handle, this->common_shared_procs );
    return true;
    }
  vert_handle = it->second;
  return false;
}

template< int _n >
void MBSplitVertices<_n>::assign_global_ids( std::map<MBProcessSet,int>& gids )
{
  typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator it;
  for ( it = this->begin(); it != this->end(); ++ it )
    {
    int gid = gids[it->first.process_set] ++;
    this->tag_manager->set_gid( it->second, gid );
    std::cout << "Assigning entity: " << it->first << " GID: " << gid << "\n";
    }
}

#endif /* MB_SPLITVERTICES_HPP */
