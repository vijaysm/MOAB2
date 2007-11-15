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
#include "MBCN.hpp"

#include <vector>
#include <algorithm>

#define MB_MESH_SET_COMPACT_PARENT_CHILD_LISTS


class AEntityFactory;

class MBMeshSet
{
protected:

  //! create an empty meshset
  MBMeshSet( unsigned flags );

  //! destructor
  ~MBMeshSet();

public:

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
  
  unsigned flags() const { return mFlags; }
  //! returns whether entities of meshsets know this meshset 
  int tracking()     const { return mFlags & MESHSET_TRACK_OWNER; }
  int set()          const { return mFlags & MESHSET_SET; }
  int ordered()      const { return mFlags & MESHSET_ORDERED; }
  int vector_based() const { return ordered(); }

   //!  PURE VIRTUAL FUNCTIONS overwritten by derived classes  *******************

  inline MBErrorCode clear( MBEntityHandle myhandle, AEntityFactory* adjacencies );

  inline MBErrorCode get_entities(std::vector<MBEntityHandle>& entities) const;
  
  inline MBErrorCode get_entities(MBRange& entities) const;
  
    //! get all entities in this MeshSet with the specified type
  inline MBErrorCode get_entities_by_type(MBEntityType entity_type,
      std::vector<MBEntityHandle> &entity_list) const;
  
  inline MBErrorCode get_entities_by_type(MBEntityType type,
      MBRange& entity_list) const;
      
  inline MBErrorCode get_entities_by_dimension( int dimension,
      std::vector<MBEntityHandle> &entity_list) const;

  inline MBErrorCode get_entities_by_dimension( int dimension,
      MBRange& entity_list) const;
      
  inline MBErrorCode get_non_set_entities( MBRange& range ) const;
    
  //! subtract/intersect/unite meshset_2 from/with/into meshset_1; modifies meshset_1
  inline MBErrorCode subtract(const MBMeshSet *meshset_2,
                               MBEntityHandle my_handle,
                               AEntityFactory* adjacencies);

  inline MBErrorCode intersect(const MBMeshSet *meshset_2,
                                MBEntityHandle my_handle,
                                AEntityFactory* adjacencies);

  inline MBErrorCode unite(const MBMeshSet *meshset_2,
                            MBEntityHandle my_handle,
                            AEntityFactory* adjacencies);

  //! add these entities to this meshset
  inline MBErrorCode add_entities(const MBEntityHandle *entity_handles,
                                   const int num_entities,
                                   MBEntityHandle my_handle,
                                   AEntityFactory* adjacencies);
    
    //! add these entities to this meshset
  inline MBErrorCode add_entities(const MBRange &entities,
                                   MBEntityHandle my_handle,
                                   AEntityFactory* adjacencies);
    
    //! add these entities to this meshset
  inline MBErrorCode remove_entities(const MBRange& entities,
                                      MBEntityHandle my_handle,
                                      AEntityFactory* adjacencies);
   
    //! remove these entities from this meshset
  inline MBErrorCode remove_entities(const MBEntityHandle *entities,
                                      const int num_entities,
                                      MBEntityHandle my_handle,
                                      AEntityFactory* adjacencies);

    //! return the number of entities contained in this meshset
  inline unsigned int num_entities() const;
    
    //! return the number of entities with the given type contained in this meshset
  inline unsigned int num_entities_by_type(MBEntityType entity_type) const;
    
    //! return the number of entities with the given type contained in this meshset
  inline unsigned int num_entities_by_dimension(int dimension) const;

  typedef std::vector<MBEntityHandle> LinkSet;
  
  inline unsigned long get_memory_use() const;

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
  //!Meshset propery flags
  unsigned char mFlags;
  //! If less than MANY, the number of parents stored inline in
  //! parentMeshSets.hnd.  If MANY, then parentMeshSets.ptr contains
  //! array begin and end pointers for a dynamically allocated array
  //! of parent handles.
  unsigned mParentCount : 2;
  //! If less than MANY, the number of children stored inline in
  //! childMeshSets.hnd.  If MANY, then childMeshSets.ptr contains
  //! array begin and end pointers for a dynamically allocated array
  //! of child handles.
  unsigned mChildCount : 2;
private:
  //! Storage for parent and child lists
  CompactList parentMeshSets, childMeshSets;
#endif
};

#define MESH_SET_VIRTUAL_FUNCTIONS      \
  MBErrorCode clear( MBEntityHandle my_handle, AEntityFactory* adjacencies);                                                          \
  \
  inline MBErrorCode get_entities(std::vector<MBEntityHandle>& entities) const;       \
  \
  inline MBErrorCode get_entities(MBRange& entities) const;                           \
  \
  inline MBErrorCode get_entities_by_type(MBEntityType type,                          \
      std::vector<MBEntityHandle>& entity_list) const;                                  \
  \
  inline MBErrorCode get_entities_by_type(MBEntityType type,                          \
      MBRange& entity_list) const;                                                      \
  \
  inline MBErrorCode get_entities_by_dimension(int dimension,                          \
      std::vector<MBEntityHandle>& entity_list) const;                                  \
  \
  inline MBErrorCode get_entities_by_dimension(int dimension,                          \
      MBRange& entity_list) const;                                                      \
  \
  inline MBErrorCode get_non_set_entities( MBRange& range ) const; \
  \
  MBErrorCode subtract(const MBMeshSet*,                          \
                               MBEntityHandle my_handle,                  \
                               AEntityFactory* adjacencies);              \
  \
  MBErrorCode intersect(const MBMeshSet *meshset_2,               \
                                MBEntityHandle my_handle,                 \
                                AEntityFactory* adjacencies);             \
  \
  MBErrorCode unite(const MBMeshSet *meshset_2,                   \
                            MBEntityHandle my_handle,                     \
                            AEntityFactory* adjacencies);                 \
  \
  MBErrorCode add_entities(const MBEntityHandle *entity_handles,  \
                                   const int num_entities,                \
                                   MBEntityHandle my_handle,              \
                                   AEntityFactory* adjacencies);          \
  \
  MBErrorCode add_entities(const MBRange &entities,               \
                                   MBEntityHandle my_handle,              \
                                   AEntityFactory* adjacencies);          \
  \
  MBErrorCode remove_entities(const MBRange& entities,            \
                                      MBEntityHandle my_handle,           \
                                      AEntityFactory* adjacencies);       \
  \
  MBErrorCode remove_entities(const MBEntityHandle *entities,     \
                                      const int num_entities,             \
                                      MBEntityHandle my_handle,           \
                                      AEntityFactory* adjacencies);       \
  \
  inline unsigned int num_entities() const;                                        \
  \
  inline unsigned int num_entities_by_type(MBEntityType entity_type) const;            \
  \
  inline unsigned int num_entities_by_dimension(int dimesion) const; \
  \
  unsigned long get_memory_use() const;



class MBMeshSet_MBRange : public MBMeshSet
{
public:

  MBMeshSet_MBRange(bool track_ownership) 
    : MBMeshSet(track_ownership ? 
                MESHSET_SET|MESHSET_TRACK_OWNER :
                MESHSET_SET) {}

  MESH_SET_VIRTUAL_FUNCTIONS

private:
  MBRange mRange;

};

class MBMeshSet_Vector : public MBMeshSet
{
public:

  MBMeshSet_Vector(bool track_ownership) 
    : MBMeshSet(track_ownership ?
                MESHSET_ORDERED|MESHSET_TRACK_OWNER :
                MESHSET_ORDERED) {}

  MESH_SET_VIRTUAL_FUNCTIONS

private:
  static void vector_to_range( std::vector<MBEntityHandle>& vect, MBRange& range );

  struct not_type_test {
      inline not_type_test( MBEntityType type ) : mType(type) {}
      inline bool operator()( MBEntityHandle handle )
        { return TYPE_FROM_HANDLE(handle) != mType; }
      MBEntityType mType;
  };

  struct type_test {
      inline type_test( MBEntityType type ) : mType(type) {}
      inline bool operator()( MBEntityHandle handle )
        { return TYPE_FROM_HANDLE(handle) == mType; }
      MBEntityType mType;
  };

  struct not_dim_test {
      inline not_dim_test( int dimension ) : mDim(dimension) {}
      inline bool operator()( MBEntityHandle handle ) const
        { return MBCN::Dimension(TYPE_FROM_HANDLE(handle)) != mDim; }
      int mDim;
  };

  struct dim_test {
      inline dim_test( int dimension ) : mDim(dimension) {}
      inline bool operator()( MBEntityHandle handle ) const
        { return MBCN::Dimension(TYPE_FROM_HANDLE(handle)) == mDim; }
      int mDim;
  };


    std::vector<MBEntityHandle> mVector;
};

inline MBErrorCode 
MBMeshSet::clear( MBEntityHandle myhandle, AEntityFactory* adjacencies )
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->clear( myhandle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->clear( myhandle, adjacencies ) ;
}

inline MBErrorCode
MBMeshSet::get_entities(std::vector<MBEntityHandle>& entities) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_entities( entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_entities( entities ) ;
}

inline MBErrorCode
MBMeshSet::get_entities(MBRange& entities) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_entities( entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_entities( entities ) ;
}

  //! get all entities in this MeshSet with the specified type
inline MBErrorCode
MBMeshSet::get_entities_by_type(MBEntityType type, std::vector<MBEntityHandle> &entities) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_entities_by_type( type, entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_entities_by_type( type, entities ) ;
}

inline MBErrorCode
MBMeshSet::get_entities_by_type(MBEntityType type, MBRange& entities) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_entities_by_type( type, entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_entities_by_type( type, entities ) ;
}
    
inline MBErrorCode
MBMeshSet::get_entities_by_dimension( int dim, std::vector<MBEntityHandle> &entities) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_entities_by_dimension( dim, entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_entities_by_dimension( dim, entities ) ;
}

inline MBErrorCode
MBMeshSet::get_entities_by_dimension( int dim, MBRange& entities) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_entities_by_dimension( dim, entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_entities_by_dimension( dim, entities ) ;
}
    
inline MBErrorCode
MBMeshSet::get_non_set_entities( MBRange& entities ) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_non_set_entities( entities ) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_non_set_entities( entities ) ;
}
  
//! subtract/intersect/unite meshset_2 from/with/into meshset_1; modifies meshset_1
inline MBErrorCode
MBMeshSet::subtract( const MBMeshSet *meshset_2,
                     MBEntityHandle my_handle,
                     AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->subtract( meshset_2, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->subtract( meshset_2, my_handle, adjacencies ) ;
}

inline MBErrorCode
MBMeshSet::intersect( const MBMeshSet *meshset_2,
                      MBEntityHandle my_handle,
                      AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->intersect( meshset_2, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->intersect( meshset_2, my_handle, adjacencies ) ;
}

inline MBErrorCode
MBMeshSet::unite( const MBMeshSet *meshset_2,
                  MBEntityHandle my_handle,
                  AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->unite( meshset_2, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->unite( meshset_2, my_handle, adjacencies ) ;
}

//! add these entities to this meshset
inline MBErrorCode
MBMeshSet::add_entities( const MBEntityHandle *entities,
                         const int num_entities,
                         MBEntityHandle my_handle,
                         AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->add_entities( entities, num_entities, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->add_entities( entities, num_entities, my_handle, adjacencies ) ;
}
  
  //! add these entities to this meshset
inline MBErrorCode
MBMeshSet::add_entities( const MBRange &entities,
                         MBEntityHandle my_handle,
                         AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->add_entities( entities, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->add_entities( entities, my_handle, adjacencies ) ;
}
  
  //! add these entities to this meshset
inline MBErrorCode
MBMeshSet::remove_entities( const MBRange& entities,
                            MBEntityHandle my_handle,
                            AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->remove_entities( entities, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->remove_entities( entities, my_handle, adjacencies ) ;
}
 
  //! remove these entities from this meshset
inline MBErrorCode
MBMeshSet::remove_entities( const MBEntityHandle *entities,
                            const int num_entities,
                            MBEntityHandle my_handle,
                            AEntityFactory* adjacencies)
{ 
  return vector_based() ?
    reinterpret_cast<MBMeshSet_Vector* >(this)->remove_entities( entities, num_entities, my_handle, adjacencies ) :
    reinterpret_cast<MBMeshSet_MBRange*>(this)->remove_entities( entities, num_entities, my_handle, adjacencies ) ;
}

  //! return the number of entities contained in this meshset
inline unsigned int MBMeshSet::num_entities() const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->num_entities() :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->num_entities() ;
}
  
  //! return the number of entities with the given type contained in this meshset
inline unsigned int
MBMeshSet::num_entities_by_type(MBEntityType type) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->num_entities_by_type(type) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->num_entities_by_type(type) ;
}
  
  //! return the number of entities with the given type contained in this meshset
inline unsigned int
MBMeshSet::num_entities_by_dimension(int dim) const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->num_entities_by_dimension(dim) :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->num_entities_by_dimension(dim) ;
}

inline unsigned long MBMeshSet::get_memory_use() const
{ 
  return vector_based() ?
    reinterpret_cast<const MBMeshSet_Vector* >(this)->get_memory_use() :
    reinterpret_cast<const MBMeshSet_MBRange*>(this)->get_memory_use() ;
}

inline MBErrorCode MBMeshSet_MBRange::get_entities(std::vector<MBEntityHandle>& entity_list) const
{
  size_t old_size = entity_list.size();
  entity_list.resize(old_size + mRange.size());
  std::copy( mRange.begin(), mRange.end(), entity_list.begin()+old_size );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_MBRange::get_entities(MBRange& entity_list) const
{
  entity_list.merge( mRange );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_MBRange::get_entities_by_type(MBEntityType type,
    std::vector<MBEntityHandle>& entity_list) const
{
  std::pair<MBRange::const_iterator,MBRange::const_iterator> its;
  its = mRange.equal_range( type );
  std::copy( its.first, its.second, std::back_inserter( entity_list ) );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_MBRange::get_entities_by_type(MBEntityType type,
    MBRange& entity_list) const
{
  std::pair<MBRange::const_iterator,MBRange::const_iterator> its;
  its = mRange.equal_range( type );
  entity_list.merge( its.first, its.second );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_MBRange::get_entities_by_dimension(int dim,
    std::vector<MBEntityHandle>& entity_list) const
{
  MBRange::const_iterator beg = mRange.lower_bound(MBCN::TypeDimensionMap[dim].first);
  MBRange::const_iterator end = mRange.lower_bound(MBCN::TypeDimensionMap[dim+1].first);
  std::copy( beg, end, std::back_inserter( entity_list ) );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_MBRange::get_entities_by_dimension(int dim,
    MBRange& entity_list) const
{
  MBRange::const_iterator beg = mRange.lower_bound(MBCN::TypeDimensionMap[dim].first);
  MBRange::const_iterator end = mRange.lower_bound(MBCN::TypeDimensionMap[dim+1].first);
  entity_list.merge( beg, end );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_MBRange::get_non_set_entities( MBRange& range ) const
{
  range.merge( mRange.begin(), mRange.lower_bound( MBENTITYSET ) );
  return MB_SUCCESS;
}

inline unsigned int MBMeshSet_MBRange::num_entities() const
{
  return mRange.size();
}

inline unsigned int MBMeshSet_MBRange::num_entities_by_type(MBEntityType type) const
{
  return mRange.num_of_type( type );
}

inline unsigned int MBMeshSet_MBRange::num_entities_by_dimension(int dimension) const
{
  return mRange.num_of_dimension( dimension );
}

inline MBErrorCode MBMeshSet_Vector::get_entities(std::vector<MBEntityHandle>& entity_list) const
{
  size_t old_size = entity_list.size();
  entity_list.resize( mVector.size() + old_size );
  std::copy( mVector.begin(), mVector.end(), entity_list.begin()+old_size );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_Vector::get_entities(MBRange& entity_list) const
{
  std::vector<MBEntityHandle> tmp_vect( mVector );
  vector_to_range( tmp_vect, entity_list );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_Vector::get_entities_by_type(MBEntityType type,
    std::vector<MBEntityHandle>& entity_list) const
{
  std::remove_copy_if( mVector.begin(), mVector.end(), 
                       std::back_inserter( entity_list ),
                       not_type_test(type) );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_Vector::get_entities_by_type(MBEntityType type,
    MBRange& entity_list) const
{
  std::vector<MBEntityHandle> tmp_vect;
  get_entities_by_type( type, tmp_vect );
  vector_to_range( tmp_vect, entity_list );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_Vector::get_entities_by_dimension( int dimension, 
    std::vector<MBEntityHandle>& entity_list) const
{
  std::remove_copy_if( mVector.begin(), mVector.end(), 
                       std::back_inserter( entity_list ),
                       not_dim_test(dimension) );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_Vector::get_entities_by_dimension( int dimension, 
    MBRange& entity_list) const
{
  std::vector<MBEntityHandle> tmp_vect;
  get_entities_by_dimension( dimension, tmp_vect );
  vector_to_range( tmp_vect, entity_list );
  return MB_SUCCESS;
}

inline MBErrorCode MBMeshSet_Vector::get_non_set_entities( MBRange& entities ) const
{
  std::vector<MBEntityHandle> tmp_vect;
  std::remove_copy_if( mVector.begin(), mVector.end(), 
                       std::back_inserter( tmp_vect ),
                       type_test(MBENTITYSET) );
  vector_to_range( tmp_vect, entities );
  return MB_SUCCESS;
}



inline unsigned int MBMeshSet_Vector::num_entities() const
{
  return mVector.size();
}

inline unsigned int MBMeshSet_Vector::num_entities_by_type(MBEntityType type) const
{
#ifndef __SUNPRO_CC
  return std::count_if( mVector.begin(), mVector.end(), type_test(type) );
#else
  unsigned int result = 0;
  std::count_if( mVector.begin(), mVector.end(), type_test(type), result );
  return result;
#endif
}

inline unsigned int MBMeshSet_Vector::num_entities_by_dimension( int dim ) const
{
#ifndef __SUNPRO_CC
  return std::count_if( mVector.begin(), mVector.end(), dim_test(dim) );
#else
  unsigned int result = 0;
  std::count_if( mVector.begin(), mVector.end(), dim_test(dim), result );
  return result;
#endif
}



#endif
