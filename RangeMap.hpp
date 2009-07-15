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

/**\file RangeMap.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2007-04-25
 */

#ifndef RANGE_MAP_HPP
#define RANGE_MAP_HPP

#include <vector>
#include <algorithm>

/**\brief Map ranges of values 
 *
 * This class provides a map between ranges of values, such as
 * a map between file IDs and MBEntityHandles.  It is intended
 * for use in situations where there are relatively few insertions
 * of large contiguous ranges of values.
 */
template <typename KeyType, typename ValType, ValType NullVal = 0>
class RangeMap 
{
public:
  typedef KeyType key_type;
  typedef ValType value_type;

  struct Range {
    KeyType begin, count;
    ValType value;
    bool operator<( const Range& other ) const
      { return begin + count <= other.begin; } // equal if overlapping!
  };
  typedef typename std::vector<Range> RangeList;
  typedef typename RangeList::const_iterator iterator;
  typedef typename RangeList::const_iterator const_iterator;
  
  inline bool empty() const
    { return data.empty(); }
  
  /**\brief Insert mapping between range of keys and range of values
   * 
   * Insert mapping from [first_key, first_key+count) to [first_val, first_val+count) 
   *
   * Input range of keys many not overlap any other input range.  If it does overlap
   * an existing range, end() is returned and the internal map is not changed.
   */
  inline iterator insert( KeyType first_key, ValType first_val, KeyType count );
  
  /** Find the value corresponding to the specified key.  Returns NullVal if not found */
  inline ValType find( KeyType key ) const;
  
  /** Check if range contains key */
  inline bool exists( KeyType key ) const;
  
  /** Check if range contains key */
  inline bool intersects( KeyType start, KeyType count ) const;

  /** Remove a block of values */
  inline iterator erase( KeyType beg, KeyType count );

  inline unsigned num_ranges() const { return data.size(); }
  
  inline iterator begin() const { return data.begin(); }
  inline iterator end() const { return data.end(); }
  inline iterator lower_bound( KeyType key ) const 
    { 
      Range b = { key, 1, NullVal }; 
      return std::lower_bound( begin(), end(), b );
    }
  static inline iterator lower_bound( iterator s, iterator e, KeyType key )
    { 
      Range b = { key, 1, NullVal }; 
      return std::lower_bound( s, e, b );
    }
  inline iterator upper_bound( KeyType key ) const
    { 
      Range b = { key, 1, NullVal }; 
      return std::upper_bound( begin(), end(), b );
    }
  static inline iterator upper_bound( iterator s, iterator e, KeyType key )
    { 
      Range b = { key, 1, NullVal }; 
      return std::upper_bound( s, e, b );
    }
  inline std::pair<iterator,iterator>
    equal_range( KeyType key ) const
    { 
      Range b = { key, 1, NullVal }; 
      return std::equal_range( begin(), end(), b );
    }
    
  inline void clear() { data.clear(); }

protected:
  RangeList data;
};

template <typename KeyType, typename ValType, ValType NullVal> 
inline typename RangeMap<KeyType,ValType,NullVal>::iterator
RangeMap<KeyType,ValType,NullVal>::insert( KeyType first_key, ValType first_val, KeyType count )
{
  Range block = { first_key, count, first_val };
  typename RangeList::iterator i = std::lower_bound( data.begin(), data.end(), block );
  
  if (i == data.end()) {
    if (i != data.begin()) {
      --i;
      if (i->begin + i->count == first_key && 
          i->value + i->count == first_val) {
        i->count += count;
        return i;
      }
    }
    data.push_back( block );
    return data.end() - 1;
  }
  
  if (i->begin < first_key + count)
    return data.end();
  
  if (i->begin == first_key + count &&
      i->value == first_val + count) {
    i->begin = first_key;
    i->value = first_val;
    i->count += count;
    if (i != data.begin()) {
      count = i->count;
      --i;
      if (i->begin + i->count == first_key &&
          i->value + i->count == first_val) {
        i->count += count;
        ++i;
        i = data.erase(i);
        --i;
      }
    }
    return i;
  }
  
  if (i != data.begin()) {
    --i;
    if (i->begin + i->count == first_key &&
        i->value + i->count == first_val) {
      i->count += count;
      return i;
    }
    ++i;
  }
  
  return data.insert( i, block );
}



template <typename KeyType, typename ValType, ValType NullVal> inline
ValType RangeMap<KeyType,ValType,NullVal>::find( KeyType key ) const
{
  Range search = { key, 1, NullVal };
  typename RangeList::const_iterator i = std::lower_bound( data.begin(), data.end(), search );
  if (i == data.end() || i->begin > key)
    return NullVal;
  
  return i->value + key - i->begin;
}

template <typename KeyType, typename ValType, ValType NullVal> inline
bool RangeMap<KeyType,ValType,NullVal>::exists( KeyType key ) const
{
  Range search = { key, 1, NullVal };
  typename RangeList::const_iterator i = std::lower_bound( data.begin(), data.end(), search );
  return i != data.end() && key >= i->begin;
}

template <typename KeyType, typename ValType, ValType NullVal> inline 
bool RangeMap<KeyType,ValType,NullVal>::intersects( KeyType start, KeyType count ) const
{
  Range search = { start, count, NullVal };
  typename RangeList::const_iterator i = std::lower_bound( data.begin(), data.end(), search );
  return i != data.end() && start + count > i->begin && i->begin+i->count > start;
}

template <typename KeyType, typename ValType, ValType NullVal>
inline typename RangeMap<KeyType,ValType,NullVal>::iterator
RangeMap<KeyType,ValType,NullVal>::erase( KeyType key, KeyType count )
{
  Range search = { key, 1, NullVal };
  typename RangeList::iterator i, j;
  i = std::lower_bound( data.begin(), data.end(), search );
  
  if (i == data.end())
    return i;
  
  if (key > i->begin) {
    KeyType offset = key - i->begin;
      // special case - split existing entry
    if((offset + count) < i->count) {
      Range ins = { i->begin, offset, i->value };
      offset += count;
      i->begin += offset;
      i->value += offset;
      i->count -= offset;
      return data.insert( i, ins ) + 1;
    }
      // otherwise remove the latter part of i
    i->count = offset;
    ++i;
  }
  
    // remove any entries entirely convered by the input range
  for (j = i; j != data.end() && (j->begin + j->count) <= (key + count); ++j);
  i = data.erase(i,j);
  
    // remove beginning of last block
  if (i != data.end() && (key + count) >= i->begin) {
    KeyType offset = key + count - i->begin;
    i->begin += offset;
    i->value += offset;
    i->count -= offset;
  }
  
  return i;
}

#endif
