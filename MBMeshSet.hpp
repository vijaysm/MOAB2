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

#include "MBInterface.hpp"
#include "MBRange.hpp"

class AEntityFactory;

class MBMeshSet
{
public:

  //! create an empty meshset
  MBMeshSet( MBEntityHandle entity_handle, AEntityFactory* a_entity_factory, bool track_ownership );

  //! virtual destructor
  virtual ~MBMeshSet();

    //! get all children pointed to by this meshset
  MBErrorCode get_children(const int num_hops, std::vector<MBEntityHandle> &children) const;

    //! get all parents pointed to by this meshset
  MBErrorCode get_parents(const int num_hops, std::vector<MBEntityHandle> &parents) const;
    
  //! add a parent/child link between the meshsets; returns error if entities are already
  //! related or if child is already a parent of parent
  static MBErrorCode add_parent_child(MBMeshSet *parent_meshset, 
                             MBMeshSet *child_meshset);

  //! remove a parent/child link between the meshsets; returns error if entities
  //! are not related
  static MBErrorCode remove_parent_child(MBMeshSet *parent_meshset, 
                                MBMeshSet *child_meshset);

    //! return the number of children pointed to by this meshset
  int num_children(int *, const int num_hops) const;
    
    //! return the number of parents pointed to by this meshset
  int num_parents(int *, const int num_hops) const;

    //! add a parent to this meshset; returns true if parent was added, 0 if it was
    //! already a parent of this meshset
  int add_parent(MBMeshSet *parent);
    
    //! add a child to this meshset; returns true if child was added, 0 if it was
    //! already a child of this meshset
  int add_child(MBMeshSet *child);
    
    //! remove a parent from this meshset; returns true if parent was removed, 0 if it was
    //! not a parent of this meshset
  int remove_parent(MBMeshSet *parent);
    
    //! remove a child from this meshset; returns true if child was removed, 0 if it was
    //! not a child of this meshset
  int remove_child(MBMeshSet *child);
   
  //! returns whether entities of meshsets know this meshset 
  bool tracking() { return mTracking; }

  void set_adj_factory(AEntityFactory* adj_fact) { mAdjFact = adj_fact; }

   //!  PURE VIRTUAL FUNCTIONS overwritten by derived classes  *******************

  virtual MBErrorCode clear() = 0;

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
  virtual MBErrorCode subtract(const MBMeshSet *meshset_2) = 0;

  virtual MBErrorCode intersect(const MBMeshSet *meshset_2) = 0;

  virtual MBErrorCode unite(const MBMeshSet *meshset_2)=0;

  //! add these entities to this meshset
  virtual MBErrorCode add_entities(const MBEntityHandle *entity_handles,
                                    const int num_entities) = 0;
    
    //! add these entities to this meshset
  virtual MBErrorCode add_entities(const MBRange &entities)=0;
    
    //! add these entities to this meshset
  virtual MBErrorCode remove_entities(const MBRange& entities)=0;
   
    //! remove these entities from this meshset
  virtual MBErrorCode remove_entities(const MBEntityHandle *entities,
                                       const int num_entities)=0;

    //! return the number of entities contained in this meshset
  virtual unsigned int num_entities() const = 0;
    
    //! return the number of entities with the given type contained in this meshset
  virtual unsigned int num_entities_by_type(MBEntityType entity_type) const = 0;
    
    //! return the number of entities with the given type contained in this meshset
  virtual unsigned int num_entities_by_dimension(int dimension) const = 0;

    //! rebuild parent-child relations for geometry
  static MBErrorCode rebuild_geometry_relations();
  
protected:

    //! links to parents/children
  std::vector<MBMeshSet*> parentMeshSets, childMeshSets;
   
  MBEntityHandle mEntityHandle;
  
  //!flag to indicate whether 'tracking' is occuring on this meshset
  //ie. all entities of the meshset know they belong to this meshset 
  bool mTracking;
  
  //! adjacency factory to handle adjacencies
  AEntityFactory* mAdjFact;
};

#define MESH_SET_VIRTUAL_FUNCTIONS      \
  virtual MBErrorCode clear();                                                          \
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
  virtual MBErrorCode subtract(const MBMeshSet*);                                     \
  \
  virtual MBErrorCode intersect(const MBMeshSet *meshset_2);                          \
  \
  virtual MBErrorCode unite(const MBMeshSet *meshset_2);                              \
  \
  virtual MBErrorCode add_entities(const MBEntityHandle *entity_handles,  \
                                    const int num_entities); \
  \
  virtual MBErrorCode add_entities(const MBRange &entities);                           \
  \
  virtual MBErrorCode remove_entities(const MBRange& entities);                        \
  \
  virtual MBErrorCode remove_entities(const MBEntityHandle *entities, \
                                       const int num_entities);    \
  \
  virtual unsigned int num_entities() const;                                        \
  \
  virtual unsigned int num_entities_by_type(MBEntityType entity_type) const;            \
  \
  virtual unsigned int num_entities_by_dimension(int dimesion) const;



class MBMeshSet_MBRange : public MBMeshSet
{
public:

  MBMeshSet_MBRange(MBEntityHandle handle, AEntityFactory* adj_fact, bool track_ownership) 
    : MBMeshSet(handle, adj_fact, track_ownership) {}
  virtual ~MBMeshSet_MBRange();

  MESH_SET_VIRTUAL_FUNCTIONS

private:
  MBRange mRange;

};

class MBMeshSet_Vector : public MBMeshSet
{
public:

  MBMeshSet_Vector(MBEntityHandle handle, AEntityFactory* adj_fact, bool track_ownership) 
    : MBMeshSet(handle, adj_fact, track_ownership) {}

  virtual ~MBMeshSet_Vector();

  MESH_SET_VIRTUAL_FUNCTIONS

private:
    std::vector<MBEntityHandle> mVector;
};


#endif
