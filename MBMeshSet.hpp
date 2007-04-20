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

//
//-------------------------------------------------------------------------
// Filename      : MB_MeshSet.hpp 
// Creator       : Tim Tautges
//
// Date          : 02/01/02
//
// Owner         : Tim Tautges
//
// Description   : MB_MeshSet is the MB implementation of MeshSet
//-------------------------------------------------------------------------

#ifndef MB_MESHSET_HPP
#define MB_MESHSET_HPP

#ifndef IS_BUILDING_MB
#error "MB_MeshSet.hpp isn't supposed to be included into an application"
#endif

#include "MBTypes.h"
#include "MBRange.hpp"

#include <vector>


#define MB_MESH_SET_COMPACT_PARENT_CHILD_LISTS


class AEntityFactory;

class MBMeshSet
{
public:

  //! create an empty meshset
  MBMeshSet( bool track_ownership );

  //! virtual destructor
  virtual ~MBMeshSet();

#ifndef MB_MESH_SET_COMPACT_PARENT_CHILD_LISTS
    //! get all children pointed to by this meshset
  const MBEntityHandle* get_children( int& count_out ) const 
    { count_out = childMeshSets.size(); return &childMeshSets[0]; }

    //! get all parents pointed to by this meshset
  const MBEntityHandle* get_parents( int& count_out ) const
    { count_out = parentMeshSets.size(); return &parentMeshSets[0]; }

    //! return the number of children pointed to by this meshset
  int num_children() const { return childMeshSets.size(); }
    
    //! return the number of parents pointed to by this meshset
  int num_parents() const { return parentMeshSets.size(); }

#else
    //! get all children pointed to by this meshset
  const MBEntityHandle* get_children( int& count_out ) const 
    { 
      count_out = mChildCount;
      if (count_out < MANY)
        return childMeshSets.hnd;
      
      count_out = childMeshSets.ptr[1] - childMeshSets.ptr[0];
      return childMeshSets.ptr[0];
    }

    //! get all parents pointed to by this meshset
  const MBEntityHandle* get_parents( int& count_out ) const
    { 
      count_out = mParentCount;
      if (count_out < MANY)
        return parentMeshSets.hnd;
        
      count_out = parentMeshSets.ptr[1] - parentMeshSets.ptr[0];
      return parentMeshSets.ptr[0];
    }

    //! return the number of children pointed to by this meshset
  int num_children() const 
    {
      if (mChildCount < MANY)
        return mChildCount;
      else
        return childMeshSets.ptr[1] - childMeshSets.ptr[0];
    }
    
    //! return the number of parents pointed to by this meshset
  int num_parents() const
    {
      if (mParentCount < MANY)
        return mParentCount;
      else
        return parentMeshSets.ptr[1] - parentMeshSets.ptr[0];
    }
#endif
  

    //! add a parent to this meshset; returns true if parent was added, 0 if it was
    //! already a parent of this meshset
  int add_parent(MBEntityHandle parent);
    
    //! add a child to this meshset; returns true if child was added, 0 if it was
    //! already a child of this meshset
  int add_child(MBEntityHandle child);
    
    //! remove a parent from this meshset; returns true if parent was removed, 0 if it was
    //! not a parent of this meshset
  int remove_parent(MBEntityHandle parent);
    
    //! remove a child from this meshset; returns true if child was removed, 0 if it was
    //! not a child of this meshset
  int remove_child(MBEntityHandle child);
   
  //! returns whether entities of meshsets know this meshset 
  bool tracking() { return mTracking; }

   //!  PURE VIRTUAL FUNCTIONS overwritten by derived classes  *******************

  virtual MBErrorCode clear( MBEntityHandle myhandle, AEntityFactory* adjacencies ) = 0;

  virtual MBErrorCode get_entities(std::vector<MBEntityHandle>& entities) const = 0;
  
  virtual MBErrorCode get_entities(MBRange& entities) const = 0;
  
    //! get all entities in this MeshSet with the specified type
  virtual MBErrorCode get_entities_by_type(MBEntityType entity_type,
      std::vector<MBEntityHandle> &entity_list) const = 0;
  
  virtual MBErrorCode get_entities_by_type(MBEntityType type,
      MBRange& entity_list) const = 0;
      
  virtual MBErrorCode get_entities_by_dimension( int dimension,
      std::vector<MBEntityHandle> &entity_list) const = 0;

  virtual MBErrorCode get_entities_by_dimension( int dimension,
      MBRange& entity_list) const = 0;
      
  virtual MBErrorCode get_non_set_entities( MBRange& range ) const = 0;
    
  //! subtract/intersect/unite meshset_2 from/with/into meshset_1; modifies meshset_1
  virtual MBErrorCode subtract(const MBMeshSet *meshset_2,
                               MBEntityHandle my_handle,
                               AEntityFactory* adjacencies) = 0;

  virtual MBErrorCode intersect(const MBMeshSet *meshset_2,
                                MBEntityHandle my_handle,
                                AEntityFactory* adjacencies) = 0;

  virtual MBErrorCode unite(const MBMeshSet *meshset_2,
                            MBEntityHandle my_handle,
                            AEntityFactory* adjacencies)=0;

  //! add these entities to this meshset
  virtual MBErrorCode add_entities(const MBEntityHandle *entity_handles,
                                   const int num_entities,
                                   MBEntityHandle my_handle,
                                   AEntityFactory* adjacencies) = 0;
    
    //! add these entities to this meshset
  virtual MBErrorCode add_entities(const MBRange &entities,
                                   MBEntityHandle my_handle,
                                   AEntityFactory* adjacencies)=0;
    
    //! add these entities to this meshset
  virtual MBErrorCode remove_entities(const MBRange& entities,
                                      MBEntityHandle my_handle,
                                      AEntityFactory* adjacencies)=0;
   
    //! remove these entities from this meshset
  virtual MBErrorCode remove_entities(const MBEntityHandle *entities,
                                      const int num_entities,
                                      MBEntityHandle my_handle,
                                      AEntityFactory* adjacencies)=0;

    //! return the number of entities contained in this meshset
  virtual unsigned int num_entities() const = 0;
    
    //! return the number of entities with the given type contained in this meshset
  virtual unsigned int num_entities_by_type(MBEntityType entity_type) const = 0;
    
    //! return the number of entities with the given type contained in this meshset
  virtual unsigned int num_entities_by_dimension(int dimension) const = 0;

  typedef std::vector<MBEntityHandle> LinkSet;
  
  virtual unsigned long get_memory_use() const = 0;

protected:

  unsigned long parent_child_memory_use() const;

#ifndef MB_MESH_SET_COMPACT_PARENT_CHILD_LISTS
    //! links to parents/children
  LinkSet parentMeshSets, childMeshSets;
  
  //!flag to indicate whether 'tracking' is occuring on this meshset
  //ie. all entities of the meshset know they belong to this meshset 
  bool mTracking;
#else

public:  
    //! Possible values of mParentCount and mChildCount
  enum Count { ZERO=0, ONE=1, TWO=2, MANY=3 };
    //! If the number of entities is less than 3, store
    //! the handles directly in the hnd member.  Otherwise
    //! use the ptr member to hold the beginning and end
    //! of a dynamically allocated array.
  union CompactList {
    MBEntityHandle hnd[2];  //!< Two handles
    MBEntityHandle* ptr[2]; //!< begin and end pointers for array
  };
protected:
  //!flag to indicate whether 'tracking' is occuring on this meshset
  //!ie. all entities of the meshset know they belong to this meshset 
  bool mTracking;
  //! If less than MANY, the number of parents stored inline in
  //! parentMeshSets.hnd.  If MANY, then parentMeshSets.ptr contains
  //! array begin and end pointers for a dynamically allocated array
  //! of parent handles.
  Count mParentCount : 2;
  //! If less than MANY, the number of children stored inline in
  //! childMeshSets.hnd.  If MANY, then childMeshSets.ptr contains
  //! array begin and end pointers for a dynamically allocated array
  //! of child handles.
  Count mChildCount : 2;
private:
  //! Storage for parent and child lists
  CompactList parentMeshSets, childMeshSets;
#endif
};

#define MESH_SET_VIRTUAL_FUNCTIONS      \
  virtual MBErrorCode clear( MBEntityHandle my_handle, AEntityFactory* adjacencies);                                                          \
  \
  virtual MBErrorCode get_entities(std::vector<MBEntityHandle>& entities) const;       \
  \
  virtual MBErrorCode get_entities(MBRange& entities) const;                           \
  \
  virtual MBErrorCode get_entities_by_type(MBEntityType type,                          \
      std::vector<MBEntityHandle>& entity_list) const;                                  \
  \
  virtual MBErrorCode get_entities_by_type(MBEntityType type,                          \
      MBRange& entity_list) const;                                                      \
  \
  virtual MBErrorCode get_entities_by_dimension(int dimension,                          \
      std::vector<MBEntityHandle>& entity_list) const;                                  \
  \
  virtual MBErrorCode get_entities_by_dimension(int dimension,                          \
      MBRange& entity_list) const;                                                      \
  \
  virtual MBErrorCode get_non_set_entities( MBRange& range ) const; \
  \
  virtual MBErrorCode subtract(const MBMeshSet*,                          \
                               MBEntityHandle my_handle,                  \
                               AEntityFactory* adjacencies);              \
  \
  virtual MBErrorCode intersect(const MBMeshSet *meshset_2,               \
                                MBEntityHandle my_handle,                 \
                                AEntityFactory* adjacencies);             \
  \
  virtual MBErrorCode unite(const MBMeshSet *meshset_2,                   \
                            MBEntityHandle my_handle,                     \
                            AEntityFactory* adjacencies);                 \
  \
  virtual MBErrorCode add_entities(const MBEntityHandle *entity_handles,  \
                                   const int num_entities,                \
                                   MBEntityHandle my_handle,              \
                                   AEntityFactory* adjacencies);          \
  \
  virtual MBErrorCode add_entities(const MBRange &entities,               \
                                   MBEntityHandle my_handle,              \
                                   AEntityFactory* adjacencies);          \
  \
  virtual MBErrorCode remove_entities(const MBRange& entities,            \
                                      MBEntityHandle my_handle,           \
                                      AEntityFactory* adjacencies);       \
  \
  virtual MBErrorCode remove_entities(const MBEntityHandle *entities,     \
                                      const int num_entities,             \
                                      MBEntityHandle my_handle,           \
                                      AEntityFactory* adjacencies);       \
  \
  virtual unsigned int num_entities() const;                                        \
  \
  virtual unsigned int num_entities_by_type(MBEntityType entity_type) const;            \
  \
  virtual unsigned int num_entities_by_dimension(int dimesion) const; \
  \
  virtual unsigned long get_memory_use() const;



class MBMeshSet_MBRange : public MBMeshSet
{
public:

  MBMeshSet_MBRange(bool track_ownership) 
    : MBMeshSet(track_ownership) {}
  virtual ~MBMeshSet_MBRange();

  MESH_SET_VIRTUAL_FUNCTIONS

private:
  MBRange mRange;

};

class MBMeshSet_Vector : public MBMeshSet
{
public:

  MBMeshSet_Vector(bool track_ownership) 
    : MBMeshSet(track_ownership) {}

  virtual ~MBMeshSet_Vector();

  MESH_SET_VIRTUAL_FUNCTIONS

private:
    std::vector<MBEntityHandle> mVector;
};


#endif
