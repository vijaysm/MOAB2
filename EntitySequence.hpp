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

/*!
 *  \class   MBEntitySequence
 *  \authors Karl Merkley & Corey Ernst
 *  \date    3/27/02
 *  \brief   The MBEntitySequence represents a contiguous range of 
 *           mesh entities of a single type.  MBEntitySequence manages
 *           the internal ids of those entities using a starting id and 
 *           the number of entities in the range.  All MBEntitySequences 
 *           for a given MBEntity type are stored in an stl Map.  When a 
 *           MBHandle references an entity in mdb, mdb finds the Entity-
 *           Sequence associated with the MBHandle from the stl Map.  
 *          
 */ 

#ifndef ENTITY_SEQUENCE_HPP
#define ENTITY_SEQUENCE_HPP

#ifndef IS_BUILDING_MB
#error "EntitySequence.hpp isn't supposed to be included into an application"
#endif

#include "MBForward.hpp"
#include "MBInternals.hpp"
#include "MBCN.hpp"
#include <set>
#include <assert.h>

class EntitySequenceManager;
class MBCore;

class MBEntitySequence
{
public:

  //! Constructor-- takes a start handle, and number of entities
  //!  if all_handles_used is true, then there is no room to create
  //!  new entities
  MBEntitySequence(EntitySequenceManager* manager, MBEntityHandle start_handle, int num_entities);

  //! Destructor
  virtual ~MBEntitySequence();

  //! gets starting MBEntityHandle (mStartEntityHandle) 
  //! note: this doesn't account for unused entity handles
  MBEntityHandle get_start_handle() const { return mStartEntityHandle; }
 
  //! gets the end of allocated entities, does not account for holes 
  MBEntityHandle get_end_handle() const { return mStartEntityHandle + mNumAllocated - 1; }

  //! get type of elements in EntitySequence
  MBEntityType get_type() const { return TYPE_FROM_HANDLE(mStartEntityHandle); }

  //! returns number of MBEntityHandles in MBEntitySequence
  int number_entities() const {return mNumEntities;}

  //! returns the number of entities space was allocated for
  int number_allocated() const { return mNumAllocated; }
  
  inline bool is_valid_entity(MBEntityHandle entity) const;

  //! get an unused handle
  virtual MBEntityHandle get_unused_handle() = 0;

  //! free an entity handle
  virtual void free_handle(MBEntityHandle handle) = 0;

  //! get entities in this range, will not add unused entities
  virtual void get_entities(MBRange& entities) const = 0;

  //! split this sequence range [begin, end] into
  //! [begin, split_location-1] and [split_location, end]
  //! returns the new_entitiy sequence which holds [split_location, end]
  //! returns success
  virtual MBErrorCode split(MBEntityHandle split_location, MBEntitySequence*& new_sequence)
  {
    split_location = split_location;
    new_sequence = NULL;
    return MB_FAILURE;
  }
  
  virtual int get_next_free_index( int prev_free_index ) const = 0;

protected:

  EntitySequenceManager* mSequenceManager;

  //!id to 1st element in EntitySequence
  MBEntityHandle       mStartEntityHandle;

  //!number of Entities in EntitySequence 
  int                   mNumEntities;
  
  //!the number of entities that space has been allocated for
  int                   mNumAllocated;

  //! the head to the list of unused handles, -1 is no free handles
  int                   mFirstFreeIndex;

  //! the last entity deleted, for speed when deleting entities sequentially
  //! -1 is means no last index
  int                   mLastDeletedIndex;

  //! a list of whether entities are free or not
  std::vector<bool>     mFreeEntities;

};


//! entity sequence for vertex type entities
class VertexEntitySequence : public MBEntitySequence
{
public:
  //! constructor
  VertexEntitySequence(EntitySequenceManager* seq_manager,
        MBEntityHandle start_handle, int num_entities, bool all_handles_used);
  //! destructor
  virtual ~VertexEntitySequence();

  MBEntityHandle get_unused_handle();
  virtual void free_handle(MBEntityHandle handle);

  MBErrorCode get_coordinates(MBEntityHandle entity,
                               double& x, double& y, double& z) const;
  
  MBErrorCode get_coordinates(MBEntityHandle entity,
                               double* xyz) const;
  
  MBErrorCode get_coordinates_ref(MBEntityHandle entity,
                               const double*& x, const double*& y, const double*& z) const;
  
  MBErrorCode set_coordinates(MBEntityHandle entity,
                               const double& x,
                               const double& y,
                               const double& z);
  
  MBErrorCode get_coordinate_arrays(double*& x, double*& y, double*& z);

  virtual void get_entities(MBRange& entities) const;
  
  virtual int get_next_free_index( int prev_free_index ) const;

private:

  // coordinate arrays x,y,z
  double** mCoords;
  
};


class ElementEntitySequence: public MBEntitySequence
{
public:
  ElementEntitySequence(EntitySequenceManager* seq_manager, MBEntityHandle start_handle, int num_entities,
                        int nodes_per_element, bool all_handles_used,
                        bool allocate_connect = true);
  virtual ~ElementEntitySequence();

  virtual MBEntityHandle get_unused_handle();
  virtual void free_handle(MBEntityHandle handle);

  virtual unsigned int nodes_per_element() const { return mNodesPerElement; }

  virtual MBErrorCode get_connectivity(MBEntityHandle entity, 
                                       std::vector<MBEntityHandle>& connectivity,
                                       const bool topological_connectivity = false) const;
  virtual MBErrorCode get_connectivity(MBEntityHandle entity, 
                                        const MBEntityHandle*& connectivity,
                                        int &num_vertices,
                                       const bool topological_connectivity = false) const;

  virtual MBErrorCode set_connectivity(MBEntityHandle entity, const MBEntityHandle *conn,
                                const int num_vertices);

  virtual MBErrorCode get_connectivity_array(MBEntityHandle*& conn_array);
  
  virtual void get_entities(MBRange& entities) const;
  
  virtual int get_next_free_index( int prev_free_index ) const;
  
  virtual MBErrorCode split(MBEntityHandle split_location, MBEntitySequence*& new_sequence);

  // reallocated the sequence to hold extra/less nodes, pass in what you want, and will return whether it needed
  // reallocate space for those nodes
  virtual MBErrorCode convert_realloc(bool& mid_edge_nodes, bool& mid_face_nodes, bool& mid_volume_nodes, 
      MBCore* MB, MBTag bit_delete_mark );
  
  virtual bool has_mid_edge_nodes() const;
  virtual bool has_mid_face_nodes() const;
  virtual bool has_mid_volume_nodes() const;

protected:
  
  unsigned short mNodesPerElement;

  MBEntityHandle* mElements;
  bool tag_for_deletion( MBEntityHandle &node, MBCore *MB );

};


inline bool MBEntitySequence::is_valid_entity(MBEntityHandle entity) const
{
  return mLastDeletedIndex == -1 || !mFreeEntities[entity-mStartEntityHandle];
}



inline MBErrorCode VertexEntitySequence::get_coordinates(MBEntityHandle handle, 
    double& x, double& y, double& z) const
{
  unsigned int index = handle - mStartEntityHandle;
  x = mCoords[0][index];
  y = mCoords[1][index];
  z = mCoords[2][index];
  return MB_SUCCESS;
}

inline MBErrorCode VertexEntitySequence::get_coordinates(MBEntityHandle handle, 
    double* xyz) const
{
  unsigned int index = handle - mStartEntityHandle;
  xyz[0] = mCoords[0][index];
  xyz[1] = mCoords[1][index];
  xyz[2] = mCoords[2][index];
  return MB_SUCCESS;
}

inline MBErrorCode VertexEntitySequence::get_coordinates_ref(MBEntityHandle handle, 
    const double*& x, const double*& y, const double*& z) const
{
  unsigned int index = handle - mStartEntityHandle;
  x = &mCoords[0][index];
  y = &mCoords[1][index];
  z = &mCoords[2][index];
  return MB_SUCCESS;
}

inline MBErrorCode VertexEntitySequence::set_coordinates(MBEntityHandle handle, 
    const double& x, const double& y, const double& z)
{
  unsigned int index = handle - mStartEntityHandle;
  mCoords[0][index] = x;
  mCoords[1][index] = y;
  mCoords[2][index] = z;
  return MB_SUCCESS;
}

inline MBErrorCode VertexEntitySequence::get_coordinate_arrays(double*& x, double*& y, double*& z)
{
  x = mCoords[0];
  y = mCoords[1];
  z = mCoords[2]; 
  return MB_SUCCESS;
}

inline MBErrorCode ElementEntitySequence::get_connectivity(MBEntityHandle entity,
                                                           const MBEntityHandle*& conn,
                                                           int &num_vertices,
                                                           const bool topological_connectivity) const
{
  num_vertices = mNodesPerElement;
  int index = entity - mStartEntityHandle;
  conn = mElements+index*mNodesPerElement;
  num_vertices = (topological_connectivity ? MBCN::VerticesPerEntity(TYPE_FROM_HANDLE(entity))
              : mNodesPerElement);
  return MB_SUCCESS;
}

inline MBErrorCode ElementEntitySequence::get_connectivity(MBEntityHandle entity,
                                                           std::vector<MBEntityHandle> &conn,
                                                           const bool topological_connectivity) const
{
  int numv = (topological_connectivity ? MBCN::VerticesPerEntity(TYPE_FROM_HANDLE(entity))
              : mNodesPerElement);
  conn.reserve(numv);
  int index = entity - mStartEntityHandle;
  MBErrorCode result = MB_SUCCESS;
  if (!is_valid_entity(entity)) result = MB_FAILURE;
  else
    conn.insert(conn.end(), mElements+index*mNodesPerElement, 
                mElements+index*mNodesPerElement+numv);
  return result;
}

inline MBErrorCode ElementEntitySequence::set_connectivity(MBEntityHandle entity, 
                                                            const MBEntityHandle *conn,
                                                            const int num_vertices)
{
  if(num_vertices != mNodesPerElement)
    return MB_FAILURE;

  MBEntityHandle* iter = &mElements[(entity - mStartEntityHandle) * mNodesPerElement];
  std::copy(conn, (conn+num_vertices), iter);
  return MB_SUCCESS;
}


inline MBErrorCode ElementEntitySequence::get_connectivity_array(MBEntityHandle*& conn_array)
{
  conn_array = mElements;
  return MB_SUCCESS;
}



#endif

