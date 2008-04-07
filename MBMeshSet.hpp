#ifndef MB_MESHSET_HPP
#define MB_MESHSET_HPP

#ifndef IS_BUILDING_MB
#error "MB_MeshSet.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBRange.hpp"
#include "MBCN.hpp"

#include <assert.h>
#include <vector>
#include <algorithm>
#include <iterator>

class AEntityFactory;

/** \brief Class to implement entity set functionality 
  * \author Jason Kraftcheck <kraftche@cae.wisc.edu>
  */
class MBMeshSet
{
public:

  //! create an empty meshset
  inline MBMeshSet();
  inline MBMeshSet(unsigned flags);

  //! destructor
  inline ~MBMeshSet();
  
  inline MBErrorCode set_flags( unsigned flags, MBEntityHandle my_handle, AEntityFactory* adjacencies );
    

    //! get all children pointed to by this meshset
  inline const MBEntityHandle* get_children( int& count_out ) const ;

    //! get all parents pointed to by this meshset
  inline const MBEntityHandle* get_parents( int& count_out ) const;

    //! return the number of children pointed to by this meshset
  inline int num_children() const ;
    
    //! return the number of parents pointed to by this meshset
  inline int num_parents() const;

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

    //! replace one entity with another in the set (contents and parent/child
    //! lists); returns whether it was replaced or not
  bool replace_entities(MBEntityHandle my_handle, 
                        const MBEntityHandle *entities, 
                        size_t num_entities,
                        AEntityFactory* mAdjFact);
  
    /** Clear *contents* of set (not parents or children) */
  inline MBErrorCode clear( MBEntityHandle myhandle, AEntityFactory* adjacencies );
  
    /** Clear all set lists (contents, parents, and children) */
  inline MBErrorCode clear_all( MBEntityHandle myhandle, AEntityFactory* adjacencies );

    /** Get contents data array.  NOTE: this may not contain what you expect if not vector_based */
  inline const MBEntityHandle* get_contents( size_t& count_out ) const;
    /** Get contents data array.  NOTE: this may not contain what you expect if not vector_based */
  inline MBEntityHandle* get_contents( size_t& count_out );

    /** Get entities contained in set */
  inline MBErrorCode get_entities(std::vector<MBEntityHandle>& entities) const;

    /** Get entities contained in set */
  inline MBErrorCode get_entities( MBRange& entities ) const;
  
    //! get all entities in this MeshSet with the specified type
  inline MBErrorCode get_entities_by_type(MBEntityType entity_type, std::vector<MBEntityHandle> &entity_list) const;
  
  inline MBErrorCode get_entities_by_type( MBEntityType type, MBRange& entity_list) const;
    
    //! return the number of entities with the given type contained in this meshset
  inline unsigned int num_entities_by_type(MBEntityType type) const;
      
  inline MBErrorCode get_entities_by_dimension( int dimension, std::vector<MBEntityHandle> &entity_list) const;

  inline MBErrorCode get_entities_by_dimension( int dimension, MBRange& entity_list) const;
  
    //! return the number of entities with the given type contained in this meshset
  inline unsigned int num_entities_by_dimension(int dimension) const;
      
  inline MBErrorCode get_non_set_entities( MBRange& range ) const;

  /** Test of meshset contains some or all of passed entities
   *
   *\param entities Array of entities
   *\param num_entities Length of array of entities.
   *\param op - MBInterface::UNION     : Test if set contains any of the input entities
   *          - MBInterface::INTERSECT : Test if set contains all of the input entities
   */
  inline bool contains_entities(const MBEntityHandle *entities, int num_entities, const int op) const;

  
  //! subtract/intersect/unite meshset_2 from/with/into meshset_1; modifies meshset_1
  inline MBErrorCode subtract(const MBMeshSet *meshset_2,
                               MBEntityHandle my_handle,
                               AEntityFactory* adjacencies);

  MBErrorCode intersect(const MBMeshSet *meshset_2,
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
  
  inline bool empty() const { return mContentCount == ZERO; }
  
  unsigned long get_memory_use() const;

protected:
  
  /** Convert for changing flag values */
  MBErrorCode convert( unsigned flags, MBEntityHandle my_handle, AEntityFactory* adj );

  /** Add explicit adjacencies from all contained entities to this (i.e. convert to tracking) */
  MBErrorCode create_adjacencies( MBEntityHandle myhandle, AEntityFactory* adjacencies );

  /** Remvoe explicit adjacencies from all contained entities to this (i.e. convert from tracking) */
  MBErrorCode remove_adjacencies( MBEntityHandle myhandle, AEntityFactory* adjacencies );

  /** Insert vector of handles into MBMeshSet */
  MBErrorCode insert_entity_vector( const MBEntityHandle* vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj );
 
  /** Insert vector of handle range pairs into MBMeshSet */
  MBErrorCode insert_entity_ranges( const MBEntityHandle* range_vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj );

  /** Insert MBRange of handles into MBMeshSet */
  MBErrorCode insert_entity_ranges( const MBRange& range, MBEntityHandle my_h, AEntityFactory* adj );

  /** Remove vector of handles from MBMeshSet */
  MBErrorCode remove_entity_vector( const MBEntityHandle* vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj );

  /** Remove vector of handle range pairs from MBMeshSet */
  MBErrorCode remove_entity_ranges( const MBEntityHandle* range_vect, size_t len, MBEntityHandle my_h, AEntityFactory* adj );

  /** Remove MBRange of handles from MBMeshSet */
  MBErrorCode remove_entity_ranges( const MBRange& range, MBEntityHandle my_h, AEntityFactory* adj );

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

private:
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
  //! If less than MANY, the number of children stored inline in
  //! contentList.hnd.  If MANY, then contentList.ptr contains
  //! array begin and end pointers for a dynamically allocated array..
  unsigned mContentCount : 2;
  //! Storage for data lists
  CompactList parentMeshSets, childMeshSets, contentList;

public:
    /** get dimension of enity */
  static inline int DIM_FROM_HANDLE( MBEntityHandle h )
    { return MBCN::Dimension( TYPE_FROM_HANDLE( h ) ); }
  
    /** Get smallest possible handle with specified dimension (first handle for first type of dimension) */
  static inline MBEntityHandle FIRST_OF_DIM( int dim )
    { return FIRST_HANDLE( MBCN::TypeDimensionMap[dim].first ); }
  
    /** Get largest possible handle with specified dimension (largest handle for last type of dimension) */
  static inline MBEntityHandle LAST_OF_DIM( int dim )
    { return LAST_HANDLE( MBCN::TypeDimensionMap[dim].second ); }

    /** functor: test if handle is not of type */
  struct not_type_test {
      inline not_type_test( MBEntityType type ) : mType(type) {}
      inline bool operator()( MBEntityHandle handle )
        { return TYPE_FROM_HANDLE(handle) != mType; }
      MBEntityType mType;
  };

    /** functor: test if handle is of type */
  struct type_test {
      inline type_test( MBEntityType type ) : mType(type) {}
      inline bool operator()( MBEntityHandle handle )
        { return TYPE_FROM_HANDLE(handle) == mType; }
      MBEntityType mType;
  };

    /** functor: test if handle is not of dimension */
  struct not_dim_test {
      inline not_dim_test( int dimension ) : mDim(dimension) {}
      inline bool operator()( MBEntityHandle handle ) const
        { return DIM_FROM_HANDLE(handle) != mDim; }
      int mDim;
  };

    /** functor: test if handle is of dimension */
  struct dim_test {
      inline dim_test( int dimension ) : mDim(dimension) {}
      inline bool operator()( MBEntityHandle handle ) const
        { return DIM_FROM_HANDLE(handle) == mDim; }
      int mDim;
  };
  
    /** Iterate over range of handles.  That is, given [first_handle,last_handle],
     *  step through all contained values.
     */
  struct hdl_iter {
    MBEntityHandle h;
    hdl_iter( MBEntityHandle val ) : h(val) {}
    hdl_iter& operator++() { ++h; return *this; }
    hdl_iter& operator--() { --h; return *this; }
    hdl_iter operator++(int) { return hdl_iter(h++); }
    hdl_iter operator--(int) { return hdl_iter(h--); }
    hdl_iter& operator+=(size_t s) { h += s; return *this; }
    hdl_iter& operator-=(size_t s) { h -= s; return *this; }
    MBEntityHandle operator*() const { return h; }
    bool operator==(hdl_iter other) const { return h == other.h; }
    bool operator!=(hdl_iter other) const { return h != other.h; }
    bool operator< (hdl_iter other) const { return h <  other.h; }
    bool operator> (hdl_iter other) const { return h >  other.h; }
    bool operator<=(hdl_iter other) const { return h <= other.h; }
    bool operator>=(hdl_iter other) const { return h >= other.h; }
    
    struct iterator_category : public std::random_access_iterator_tag { };
    typedef MBEntityID difference_type;
    typedef MBEntityHandle value_type;
    typedef MBEntityHandle* pointer;
    typedef MBEntityHandle& reference;
  };
  
};

inline MBMeshSet::hdl_iter::difference_type
 operator-( const MBMeshSet::hdl_iter a, const MBMeshSet::hdl_iter b )
  { return (MBMeshSet::hdl_iter::difference_type)a.h 
         - (MBMeshSet::hdl_iter::difference_type)b.h; }


  //! create an empty meshset
MBMeshSet::MBMeshSet()
  : mFlags(0), mParentCount(ZERO), mChildCount(ZERO), mContentCount(ZERO)
{ }

  //! create an empty meshset
MBMeshSet::MBMeshSet(unsigned flags)
  : mFlags(flags), mParentCount(ZERO), mChildCount(ZERO), mContentCount(ZERO)
{ }

  //! destructor
MBMeshSet::~MBMeshSet() 
{
  if (mChildCount == MANY)
    free( childMeshSets.ptr[0] );
  if (mParentCount == MANY)
    free( parentMeshSets.ptr[0] );
  if (mContentCount == MANY)
    free( contentList.ptr[0] );
  mChildCount = mParentCount = mContentCount = ZERO;
}
  
MBErrorCode MBMeshSet::set_flags( unsigned flags, MBEntityHandle my_handle, AEntityFactory* adjacencies ) 
{
  if (mContentCount == 0) {
    mFlags = flags;
    return MB_SUCCESS;
  }
  else 
    return convert( flags, my_handle, adjacencies );
}
    

//! get all children pointed to by this meshset
const MBEntityHandle* MBMeshSet::get_children( int& count_out ) const 
{ 
  count_out = mChildCount;
  if (count_out < MANY)
    return childMeshSets.hnd;

  count_out = childMeshSets.ptr[1] - childMeshSets.ptr[0];
  return childMeshSets.ptr[0];
}

//! get all parents pointed to by this meshset
const MBEntityHandle* MBMeshSet::get_parents( int& count_out ) const
{ 
  count_out = mParentCount;
  if (count_out < MANY)
    return parentMeshSets.hnd;

  count_out = parentMeshSets.ptr[1] - parentMeshSets.ptr[0];
  return parentMeshSets.ptr[0];
}

//! return the number of children pointed to by this meshset
int MBMeshSet::num_children() const 
{
  if (mChildCount < MANY)
    return mChildCount;
  else
    return childMeshSets.ptr[1] - childMeshSets.ptr[0];
}

//! return the number of parents pointed to by this meshset
int MBMeshSet::num_parents() const
{
  if (mParentCount < MANY)
    return mParentCount;
  else
    return parentMeshSets.ptr[1] - parentMeshSets.ptr[0];
}
  
inline MBErrorCode MBMeshSet::clear( MBEntityHandle myhandle, AEntityFactory* adjacencies )
{ 
  if (tracking())
    remove_adjacencies( myhandle, adjacencies );
  if (mContentCount == MANY)
    free( contentList.ptr[0] );
  mContentCount = ZERO;
  return MB_SUCCESS;
}
  
inline MBErrorCode MBMeshSet::clear_all( MBEntityHandle myhandle, AEntityFactory* adjacencies )
{ 
  MBErrorCode rval = clear( myhandle, adjacencies );
  if (mChildCount == MANY)
    free( childMeshSets.ptr[0] );
  mChildCount = ZERO;
  if (mParentCount == MANY)
    free( parentMeshSets.ptr[0] );
  mParentCount = ZERO;
  return rval;
}

inline const MBEntityHandle* MBMeshSet::get_contents( size_t& count_out ) const
{
  if (mContentCount == MANY) {
    count_out = contentList.ptr[1] - contentList.ptr[0];
    return contentList.ptr[0];
  }
  else {
    count_out = mContentCount;
    return contentList.hnd;
  }
}

inline MBEntityHandle* MBMeshSet::get_contents( size_t& count_out )
{
  if (mContentCount == MANY) {
    count_out = contentList.ptr[1] - contentList.ptr[0];
    return contentList.ptr[0];
  }
  else {
    count_out = mContentCount;
    return contentList.hnd;
  }
}

inline MBErrorCode MBMeshSet::get_entities(std::vector<MBEntityHandle>& entities) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    size_t old_size = entities.size();
    entities.resize( count+old_size );
    std::copy( ptr, ptr+count, entities.begin()+old_size );
  }
  else {
    assert(count%2 == 0);
    for (size_t i = 0; i < count; i += 2) 
      std::copy( hdl_iter(ptr[i]), hdl_iter(ptr[i+1]+1), std::back_inserter(entities) );
  }
  return MB_SUCCESS;
} 

inline MBErrorCode MBMeshSet::get_entities(MBRange& entities) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    std::copy( ptr, ptr+count, mb_range_inserter(entities) );
  }
  else {
    assert(count%2 == 0);
    entities.clear();
    MBRange::iterator in = entities.begin();
    for (size_t i = 0; i < count; i += 2) 
      in = entities.insert( in, ptr[i], ptr[i+1] );
  }
  return MB_SUCCESS;
} 


  //! get all entities in this MeshSet with the specified type
inline MBErrorCode MBMeshSet::get_entities_by_type(MBEntityType type,
                                                   std::vector<MBEntityHandle> &entity_list
                                                   ) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    std::remove_copy_if( ptr, ptr+count, 
                         std::back_inserter( entity_list ),
                         not_type_test(type) );
  }
  else {
    size_t idx = std::lower_bound( ptr, ptr+count, FIRST_HANDLE(type) ) - ptr;
    if (idx < count && TYPE_FROM_HANDLE(ptr[idx]) == type) {
      if (idx % 2) { // only part of first block is of type
        std::copy( hdl_iter(FIRST_HANDLE(type)), hdl_iter(ptr[idx]+1), std::back_inserter( entity_list ) );
        ++idx;
      }
      for (; idx < count; idx += 2) {
        if (TYPE_FROM_HANDLE(ptr[idx+1]) == type) // whole block is of type
          std::copy( hdl_iter(ptr[idx]), hdl_iter(ptr[idx+1]+1), std::back_inserter( entity_list ) );
        else {
          if (TYPE_FROM_HANDLE(ptr[idx]) == type) // part of last block is of type
            std::copy( hdl_iter(ptr[idx]), hdl_iter(LAST_HANDLE(type)), std::back_inserter( entity_list ) );
          break;
        }
      }
    }
  }

  return MB_SUCCESS;
}


inline MBErrorCode MBMeshSet::get_entities_by_type( MBEntityType type,
                                                    MBRange& entity_list) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    std::remove_copy_if( ptr, ptr+count, 
                         mb_range_inserter( entity_list ),
                         not_type_test(type) );
  }
  else {
    size_t idx = std::lower_bound( ptr, ptr+count, FIRST_HANDLE(type) ) - ptr;
    MBRange::iterator in = entity_list.begin();
    if (idx < count && TYPE_FROM_HANDLE(ptr[idx]) == type) {
      if (idx % 2) { // only part of first block is of type
        in = entity_list.insert( in, FIRST_HANDLE(type), ptr[idx] );
        ++idx;
      }
      for (; idx < count; idx += 2) {
        if (TYPE_FROM_HANDLE(ptr[idx+1]) == type) // whole block is of type
          in = entity_list.insert( in, ptr[idx], ptr[idx+1] );
        else {
          if (TYPE_FROM_HANDLE(ptr[idx]) == type) // part of last block is of type
            entity_list.insert( in, ptr[idx], LAST_HANDLE(type) );
          break;
        }
      }
    }
  }

  return MB_SUCCESS;
}

  //! return the number of entities with the given type contained in this meshset
inline unsigned int MBMeshSet::num_entities_by_type(MBEntityType type) const
{
  unsigned int result;
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    #ifndef __SUNPRO_CC
      result = std::count_if( ptr, ptr+count, type_test(type) );
    #else
      std::count_if( ptr, ptr+count, type_test(type), result );
    #endif
  }
  else {
    result = 0;
    size_t idx = std::lower_bound( ptr, ptr+count, FIRST_HANDLE(type) ) - ptr;
    if (idx < count && TYPE_FROM_HANDLE(ptr[idx]) == type) {
      if (idx % 2) { // only part of first block is of type
        result += ptr[idx] - FIRST_HANDLE(type) + 1;
        ++idx;
      }
      for (; idx < count; idx += 2) {
        if (TYPE_FROM_HANDLE(ptr[idx+1]) == type) // whole block is of type
          result += ptr[idx+1] - ptr[idx] + 1;
        else {
          if (TYPE_FROM_HANDLE(ptr[idx]) == type) // part of last block is of type
            result += LAST_HANDLE(type) - ptr[idx] + 1;
          break;
        }
      }
    }
  }

  return result;
}

inline MBErrorCode MBMeshSet::get_entities_by_dimension( int dimension,
                                                         std::vector<MBEntityHandle> &entity_list
                                                         ) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    std::remove_copy_if( ptr, ptr+count, 
                         std::back_inserter( entity_list ),
                         not_dim_test(dimension) );
  }
  else {
    size_t idx = std::lower_bound( ptr, ptr+count, FIRST_OF_DIM(dimension) ) - ptr;
    if (idx < count && DIM_FROM_HANDLE(ptr[idx]) == dimension) {
      if (idx % 2) { // only part of first block is of type
        std::copy( hdl_iter(FIRST_OF_DIM(dimension)), hdl_iter(ptr[idx]+1), std::back_inserter( entity_list ) );
        ++idx;
      }
      for (; idx < count; idx += 2) {
        if (DIM_FROM_HANDLE(ptr[idx+1]) == dimension) // whole block is of type
          std::copy( hdl_iter(ptr[idx]), hdl_iter(ptr[idx+1]+1), std::back_inserter( entity_list ) );
        else {
          if (DIM_FROM_HANDLE(ptr[idx]) == dimension) // part of last block is of type
            std::copy( hdl_iter(ptr[idx]), hdl_iter(LAST_OF_DIM(dimension)), std::back_inserter( entity_list ) );
          break;
        }
      }
    }
  }

  return MB_SUCCESS;
}


inline MBErrorCode MBMeshSet::get_entities_by_dimension( int dimension,
                                                         MBRange& entity_list) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    std::remove_copy_if( ptr, ptr+count, 
                         mb_range_inserter( entity_list ),
                         not_dim_test(dimension) );
  }
  else {
    size_t idx = std::lower_bound( ptr, ptr+count, FIRST_OF_DIM(dimension) ) - ptr;
    MBRange::iterator in = entity_list.begin();
    if (idx < count && DIM_FROM_HANDLE(ptr[idx]) == dimension) {
      if (idx % 2) { // only part of first block is of type
        in = entity_list.insert( in, FIRST_OF_DIM(dimension), ptr[idx] );
        ++idx;
      }
      for (; idx < count; idx += 2) {
        if (DIM_FROM_HANDLE(ptr[idx+1]) == dimension) // whole block is of type
          in = entity_list.insert( in, ptr[idx], ptr[idx+1] );
        else {
          if (DIM_FROM_HANDLE(ptr[idx]) == dimension) // part of last block is of type
            entity_list.insert( in, ptr[idx], LAST_OF_DIM(dimension) );
          break;
        }
      }
    }
  }

  return MB_SUCCESS;
}

  //! return the number of entities with the given type contained in this meshset
inline unsigned int MBMeshSet::num_entities_by_dimension(int dimension) const
{
  unsigned int result;
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    #ifndef __SUNPRO_CC
      result = std::count_if( ptr, ptr+count, dim_test(dimension) );
    #else
      std::count_if( ptr, ptr+count, dim_test(dimension), result );
    #endif
  }
  else {
    result = 0;
    size_t idx = std::lower_bound( ptr, ptr+count, FIRST_OF_DIM(dimension) ) - ptr;
    if (idx < count && DIM_FROM_HANDLE(ptr[idx]) == dimension) {
      if (idx % 2) { // only part of first block is of type
        result += ptr[idx] - FIRST_OF_DIM(dimension) + 1;
        ++idx;
      }
      for (; idx < count; idx += 2) {
        if (DIM_FROM_HANDLE(ptr[idx+1]) == dimension) // whole block is of type
          result += ptr[idx+1] - ptr[idx] + 1;
        else {
          if (DIM_FROM_HANDLE(ptr[idx]) == dimension) // part of last block is of type
            result += LAST_OF_DIM(dimension) - ptr[idx] + 1;
          break;
        }
      }
    }
  }

  return result;
}

inline MBErrorCode MBMeshSet::get_non_set_entities( MBRange& range ) const
{
  size_t count;
  const MBEntityHandle* ptr = get_contents( count );
  if (vector_based()) {
    std::remove_copy_if( ptr, ptr+count, 
                         mb_range_inserter( range ),
                         type_test(MBENTITYSET) );
  }
  else {
    MBRange::iterator in = range.begin();
    for (size_t idx = 0; idx < count; idx += 2) {
      if (TYPE_FROM_HANDLE(ptr[idx+1]) != MBENTITYSET)
        in = range.insert( ptr[idx], ptr[idx+1] );
      else {
        if (TYPE_FROM_HANDLE(ptr[idx]) != MBENTITYSET)
          range.insert( ptr[idx], LAST_HANDLE( MBENTITYSET - 1 ) );
        break;
      }
    }
  }

  return MB_SUCCESS;
}

inline bool MBMeshSet::contains_entities(const MBEntityHandle *entities, 
                                         int num_entities,
                                         const int op) const
{
  size_t count;
  const MBEntityHandle* const ptr = get_contents( count );
  const MBEntityHandle* const end = ptr + count;
  size_t found_count = 0;
  if (vector_based()) {
    for (int i = 0; i < num_entities; ++i)
      if (std::find( ptr, end, entities[i] ) < end)
        ++found_count; 
  }
  else {
    assert(0 == count % 2);
    for (int i = 0; i < num_entities; ++i) {
      const unsigned long idx = std::lower_bound( ptr, end, entities[i] ) - ptr;
      if (idx < count && (idx%2 != 0 || ptr[idx] == entities[i]))
        ++found_count;
    }
  }

  return found_count >= ((MBInterface::INTERSECT == op) ? (unsigned)num_entities : 1u);
}



//! subtract/intersect/unite meshset_2 from/with/into meshset_1; modifies meshset_1
inline MBErrorCode MBMeshSet::subtract(const MBMeshSet *meshset_2,
                                       MBEntityHandle my_handle,
                                       AEntityFactory* adjacencies)
{
  size_t count;
  const MBEntityHandle* const ptr = meshset_2->get_contents( count );
  if (meshset_2->vector_based()) 
    return remove_entity_vector( ptr, count, my_handle, adjacencies );
  else
    return remove_entity_ranges( ptr, count, my_handle, adjacencies );
}

inline MBErrorCode MBMeshSet::unite(const MBMeshSet *meshset_2,
                                    MBEntityHandle my_handle,
                                    AEntityFactory* adjacencies)
{
  size_t count;
  const MBEntityHandle* const ptr = meshset_2->get_contents( count );
  if (meshset_2->vector_based()) 
    return insert_entity_vector( ptr, count, my_handle, adjacencies );
  else
    return insert_entity_ranges( ptr, count, my_handle, adjacencies );
}

//! add these entities to this meshset
inline MBErrorCode MBMeshSet::add_entities(const MBEntityHandle *entity_handles,
                                           const int num_entities,
                                           MBEntityHandle my_handle,
                                           AEntityFactory* adjacencies)
{
  return insert_entity_vector( entity_handles, num_entities, my_handle, adjacencies );
}

  //! add these entities to this meshset
inline MBErrorCode MBMeshSet::add_entities(const MBRange &entities,
                                           MBEntityHandle my_handle,
                                           AEntityFactory* adjacencies)
{
  return insert_entity_ranges( entities, my_handle, adjacencies );
}

  //! add these entities to this meshset
inline MBErrorCode MBMeshSet::remove_entities(const MBRange& entities,
                                              MBEntityHandle my_handle,
                                              AEntityFactory* adjacencies)
{
  return remove_entity_ranges(  entities, my_handle, adjacencies );
}


  //! remove these entities from this meshset
inline MBErrorCode MBMeshSet::remove_entities(const MBEntityHandle *entities,
                                              const int num_entities,
                                              MBEntityHandle my_handle,
                                              AEntityFactory* adjacencies)
{
  return remove_entity_vector( entities, num_entities, my_handle, adjacencies );
}

  //! return the number of entities contained in this meshset
unsigned int MBMeshSet::num_entities() const
{
  size_t count;
  const MBEntityHandle* list = get_contents( count );
  if (vector_based())
    return count;
  
  int result = 0;
  const MBEntityHandle *const end = list + count;
  for (; list < end; list += 2)
    result += list[1] - list[0] + 1;
  return result;
}

#endif
