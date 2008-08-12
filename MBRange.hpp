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

/**
 * \class MBRange
 *
 * Stores contiguous or partially contiguous values in an optimized
 * fashion.  Partially contiguous accessing patterns is also optimized.
 *
 * \author Clinton Stimpson
 *
 * \date 15 April 2002
 *
 */


/*
 *************  MBRange FAQ and tips ********************

 The purpose of this FAQ is to familiarize a user with
 the appropriate use of the MBRange template.

   ******* A few points about MBRange: *******
 1.  MBRange is not the be all of generic containers.
 2.  MBRange has its strengths and weakneses as any other
     STL container has.
 3.  Strengths:
     a. For contiguous values, storage is extremely minimal.
     b. Searching through contiguous values, at best, is
        a constant time operation.
     b. Fairly compatible with most STL algorithms.
     c. Insertions of data from high value to low value
        is a linear operation (constant for each insertion).
     d. Removal of a value using an iterator is constant time.

 4.  Weaknesses:
     a. For non-contiguous values, storage is not minimal and is
        on the order of 4x the storage space as using a vector.
     b. Searching through non-contiguous values is linear
        time operation.
     c. Insertions of random data is VERY slow.

   Given the above characteristics of MBRanges, you can now
   decide between MBRange and another STL container for your
   particular needs.


   ******* Tips *******
 1.  Be aware of what iterators are available.  Currently, there are
     three.  MBRange<T>::iterator, MBRange<T>::const_iterator,
     and MBRange<T>::RangeListIterator.
     iterator is derived from const_iterator.  const_iterator
     is derived from RangeListIterator.  RangeListIterator is a
     std::list<std::pair<T,T> >::const_iterator.
     If a particular algorithm could be more efficient by using
     RangeListIterator, do so.

     ie.
     
     MBRange<char> range1;
     ... put some stuff in range1
     MBRange<char> range2;
     
     // the SLOW way.
     std::copy(range1.begin(), range1.end(), mb_range_inserter<...>(range2));

     // the FAST way.
     for(MBRange<char>::RangeListIterator iter = range1.begin(),
         iter != range1.end(); ++iter)
     {
       range2.insert(iter->first, iter->second);
     }

 2.  Prefer insert(val1, val2) over insert(val) where possible.
 
 3.  insert(val) and insert(val1, val2) have to perform searching
     to find out where to insert an item.  Searches are started
     from the beginning of the values.  Inserting larger values 
     before smaller values will increase efficiency.

     ie.
     std::set<int> my_set;
     MBRange<int> my_range;
     .. perform some operations which set does efficiently.

     // now copy the results from the set into the range.
     // copy from the end of the set to the beginning of the set
     std::copy(my_set.rbegin(), my_set.rend(), 
         mb_range_inserter< MBRange<int> > ( my_range );

 4.  Use empty() instead of size() if you only need to find out
     if there is anything in the list.

 5.  Know what swap() does.  Sometimes it is useful.  It'll replace
     the contents of one list with another.

     void compute_and_get_some_set( 
         MBRange<char> range1, 
         MBRange<char> range2,
         MBRange<char>& results
         );
     {
       MBRange<char> tmp_results;
       .. perform some computation on range1 and range2
          and put results in tmp_results;
       .. filter tmp_results out for some special type.
       .. etc....
       // return results
       results.swap(tmp_results);
     }


   ******* FAQ *******
 1. Why won't this code compile?
    ------------------------
    class SomeClass
    {
    public:
      MBRange<int> range;
    };
    ....
    void some_function( const SomeClass& some_class )
    {
      MBRange<int>::iterator = some_class.range.begin();
    }
    ---------------------
    Solution:  you need to use
    MBRange<int>::const_iterator instead.

 2. Why doesn't this work right when I try to change the
    contents of an MBRange?

    // make a range that has the letters A,B,C in it.
    MBRange<char> my_chars('A', 'C');
    // make an iterator that points to 'A'.
    MBRange<char>::iterator iter = my_chars.begin();
    // change A to D
    *iter = 'D';
    // print contents of my_chars to stdout
    std::copy(my_chars.begin(), my_chars.end(),
              std::ostream_iterator(std::cout, " "));

    result is 
      A B C
    instead of
      B C D

    When one calls *iter, which returns 'A', the actual storage of the value 'A' 
    which is returned is in the iterator and not in the MBRange.  This allows
    for multiple iterators on a single MBRange and the MBRange does not have
    to keep track of which iterator is referencing which value.



*/


#ifndef MB_RANGE_HPP
#define MB_RANGE_HPP

#include <iterator>
#include <iosfwd>
#include "MBTypes.h"

struct range_iter_tag : public std::bidirectional_iterator_tag {};

struct range_base_iter
{
  typedef range_iter_tag iterator_category;
  typedef MBEntityID difference_type;
  typedef MBEntityHandle value_type;
  typedef MBEntityHandle* pointer;
  typedef MBEntityHandle& reference;
};


//! the class MBRange
class MBRange
{
public:

    // forward declare the iterators
  class const_iterator;
  class const_reverse_iterator;
  typedef const_iterator iterator;
  typedef const_reverse_iterator reverse_iterator;
 
    //! intersect two ranges, placing the results in the return range
  MBRange intersect(const MBRange &range2) const;

    //! subtract range2 from this, placing the results in the return range
  MBRange subtract(const MBRange &range2) const;

  //! for short hand notation, lets typedef the 
  //! container class that holds the ranges
  typedef MBEntityHandle value_type;

  //! default constructor
  MBRange();

    //! copy constructor
  MBRange(const MBRange& copy);

  //! another constructor that takes an initial range
  MBRange( MBEntityHandle val1, MBEntityHandle val2 );

    //! operator=
  MBRange& operator=(const MBRange& copy);
  
  //! destructor
  ~MBRange();

  //! return the beginning const iterator of this range
  const_iterator begin() const;
  
  //! return the beginning const reverse iterator of this range
  const_reverse_iterator rbegin() const;
 
  //! return the ending const iterator for this range
  const_iterator end() const;
  
  //! return the ending const reverse iterator for this range
  const_reverse_iterator rend() const;

  //! return the number of values this Ranges represents
  MBEntityHandle size() const;
  
  //! return whether empty or not 
  //! always use "if(!Ranges::empty())" instead of "if(Ranges::size())"
  bool empty() const;

  //! insert an item into the list and return the iterator for the inserted item
  iterator insert(MBEntityHandle val);
  
  //! insert a range of items into this list and return the iterator for the first
  //! inserted item
  iterator insert(MBEntityHandle val1, MBEntityHandle val2);
  
    //! remove an item from this list and return an iterator to the next item
  iterator erase(iterator iter);

  //! remove a range of items from the list
  iterator erase( iterator iter1, iterator iter2);

  //! erases a value from this container
  iterator erase(MBEntityHandle val);
  
  //! get first entity in range
  const MBEntityHandle& front() const;
  //! get last entity in range
  const MBEntityHandle& back() const;
  //! remove first entity from range
  MBEntityHandle pop_front();
  //! remove last entity from range
  MBEntityHandle pop_back();
  
  //! find an item int the list and return an iterator at that value
  const_iterator find(MBEntityHandle val) const;

  //! return an iterator to the first value >= val
  static const_iterator lower_bound(const_iterator first,
                                    const_iterator last,
                                    MBEntityHandle val);
  static const_iterator upper_bound(const_iterator first,
                                    const_iterator last,
                                    MBEntityHandle val);
  
  const_iterator lower_bound( MBEntityType type ) const;
  const_iterator upper_bound( MBEntityType type ) const;
  std::pair<const_iterator, const_iterator> equal_range( MBEntityType type ) const;
  const_iterator lower_bound( MBEntityType type, const_iterator first ) const;
  const_iterator upper_bound( MBEntityType type, const_iterator first ) const;
  
  //! True if all entities in range are of passed type 
  //! (also true if range is empty)
  bool all_of_type( MBEntityType type ) const;
  //! True if all entities in range are of passed dimension 
  //! (also true if range is empty)
  bool all_of_dimension( int dimension ) const;
  
  unsigned num_of_type( MBEntityType type ) const;
  unsigned num_of_dimension( int dim ) const;
  
  //! clears the contents of the list 
  void clear();
  
  //! for debugging
  void print(const char *indent_prefix = NULL) const;
  void print(std::ostream& s, const char *indent_prefix = NULL) const;
  
  unsigned long get_memory_use() const;

  //! merges this MBRange with another range
  void merge( const MBRange& range );
  
  //! merge a subset of some other range
  void merge( MBRange::const_iterator begin,
              MBRange::const_iterator end );

  //! swap the contents of this range with another one
  void swap( MBRange &range );

    //! check for internal consistency
  void sanity_check() const;

    //! return a subset of this range, by type
  MBRange subset_by_type(MBEntityType t) const;
  
  struct PairNode : public std::pair<MBEntityHandle,MBEntityHandle>
  {

    PairNode() : std::pair<MBEntityHandle,MBEntityHandle>(0, 0), mNext(NULL), mPrev(NULL) {}
    PairNode(PairNode* next, PairNode* prev, 
             MBEntityHandle _first, MBEntityHandle _second)
      : std::pair<MBEntityHandle,MBEntityHandle>(_first,_second), mNext(next), mPrev(prev) {}

    PairNode* mNext;
    PairNode* mPrev;
  };

  
  MBEntityHandle operator[](MBEntityID index) const;

  int index(MBEntityHandle handle) const;
  
  MBRange::iterator insert( MBRange::iterator prev,
                            MBEntityHandle first,
                            MBEntityHandle last );
  
protected:

  //! the head of the list that contains pairs that represent the ranges 
  //! this list is sorted and unique at all times
  PairNode mHead;
  
  //! if dead_node is not mHead, remove it from the list and free it's memory.
  void delete_pair_node( PairNode* dead_node );

public:

    //! used to iterate over sub-ranges of a range
  class pair_iterator : public range_base_iter
  {
    friend class MBRange;
  public:
    pair_iterator() : mNode(NULL) {}
    pair_iterator(PairNode *node) : mNode(node) {}
    pair_iterator(const pair_iterator& copy)
      : mNode(copy.mNode) {}
    pair_iterator(const const_iterator& copy)
      : mNode(copy.mNode) {}

    std::pair<MBEntityHandle,MBEntityHandle>* operator->() { return mNode; }
    
    pair_iterator& operator++()
    {
      mNode = mNode->mNext;
      return *this;
    }
    pair_iterator operator++(int)
    {
      pair_iterator tmp(*this);
      this->operator ++();
      return tmp;
    }

    pair_iterator& operator--()
    {
      mNode = mNode->mPrev;
      return *this;
    }
    pair_iterator operator--(int)
    {
      pair_iterator tmp(*this);
      this->operator--();
      return tmp;
    }
    bool operator==(const pair_iterator& other) const
    {
      return mNode == other.mNode;
    }

    bool operator!=(const pair_iterator& other) const
    {
      return mNode != other.mNode;
    }

  private:
    
    PairNode* mNode;
  };

  class const_pair_iterator;

  //! a const iterator which iterates over an MBRange
  class const_iterator : public range_base_iter
  {
    friend class MBRange;
    friend class pair_iterator;
    friend class const_pair_iterator;
  public:
    //! default constructor - intialize base default constructor
    const_iterator() : mNode(NULL), mValue(0) {}

    //! constructor used by MBRange
    const_iterator( const PairNode* iter, const MBEntityHandle val) 
      : mNode(const_cast<PairNode*>(iter)), mValue(val)  {} 

    //! dereference that value this iterator points to
    //! returns a const reference
    const MBEntityHandle& operator*() const { return  mValue; }

    //! prefix incrementer
    const_iterator& operator++()
    {
      // see if we need to increment the base iterator
      if(mValue == mNode->second)
      {
        mNode = mNode->mNext;
        mValue = mNode->first;
      }
      // if not, just increment the value in the range
      else
        ++mValue;
      return *this;
    }

    //! postfix incrementer
    const_iterator operator++(int)
    {
      // make a temporary copy
      const_iterator tmp(*this);
      // increment self
      this->operator ++();
      // return the copy
      return tmp;
    }

    //! prefix decrementer
    const_iterator& operator--()
    {
      // see if we need to decrement the base iterator
      if(mValue == mNode->first)
      {
        mNode = mNode->mPrev;;
        mValue = mNode->second;
      }
      // if not, just decrement the value
      else
        --mValue;
      return *this;
    }

    //! postfix decrementer
    const_iterator operator--(int)
    {
      // make a copy of this
      const_iterator tmp(*this);
      // decrement self
      this->operator --();
      // return the copy
      return tmp;
    }
    
    //! Advance iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator++ step times.
    const_iterator& operator+=( MBEntityID step );
    
    //! Regress iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator-- step times.
    const_iterator& operator-=( MBEntityID step );

    //! equals operator
    bool operator==( const const_iterator& other ) const
    {
      // see if the base iterator is the same and the
      // value of this iterator is the same
      return (mNode == other.mNode) && (mValue == other.mValue);
    }

    //! not equals operator
    bool operator!=( const const_iterator& other ) const
    {
      // call == operator and not it.
      return (mNode != other.mNode) || (mValue != other.mValue);
    }
    
    /**\brief get an iterator at the end of the block
     *
     * Get an iterator at the end of the block of consecutive
     * handles that this iterator is currently contained in.
     * That is, if the range contains blocks of consecutive 
     * handles of the form { [1,5], [7,100], ... } and this
     * iterator is at any handle in the range [7,100], return
     * an iterator at the '100' handle.
     *
     * Never returns begin() or end() unless this iterator is
     * at begin() or end().  May return the same location as
     * this iterator.
     */
    inline const_iterator end_of_block() const;

  protected:

    //! the node we are pointing at
    PairNode* mNode;
    //! the value in the range
    MBEntityHandle mValue;
  };

  //! a const reverse iterator which iterates over an MBRange
  class const_reverse_iterator : public range_base_iter
  {
    friend class MBRange;
    friend class pair_iterator;
  public:
    //! default constructor - intialize base default constructor
    const_reverse_iterator() {}
    
    const_reverse_iterator( const_iterator fwd_iter ) : myIter(fwd_iter) {}

    //! constructor used by MBRange
    const_reverse_iterator( const PairNode* iter, const MBEntityHandle val) 
      : myIter(iter, val)  {} 

    //! dereference that value this iterator points to
    //! returns a const reference
    const MBEntityHandle& operator*() const { return  *myIter; }

    //! prefix incrementer
    const_reverse_iterator& operator++()
    {
      --myIter;
      return *this;
    }

    //! postfix incrementer
    const_reverse_iterator operator++(int)
    {
      return const_reverse_iterator( myIter-- );
    }

    //! prefix decrementer
    const_reverse_iterator& operator--()
    {
      ++myIter;
      return *this;
    }

    //! postfix decrementer
    const_reverse_iterator operator--(int)
    {
      return const_reverse_iterator( myIter++ );
    }
    
    //! Advance iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator++ step times.
    const_reverse_iterator& operator+=( MBEntityID step )
    {
      myIter -= step;
      return *this;
    }

    //! Regress iterator specified amount.
    //! Potentially O(n), but typically better.  Always
    //! more efficient than calling operator-- step times.
    const_reverse_iterator& operator-=( MBEntityID step )
    {
      myIter += step;
      return *this;
    }

    //! equals operator
    bool operator==( const const_reverse_iterator& other ) const
    {
      return myIter == other.myIter;
    }

    //! not equals operator
    bool operator!=( const const_reverse_iterator& other ) const
    {
      return myIter != other.myIter;
    }

  protected:

    //! the node we are pointing at
    const_iterator myIter;
  };

public:

  class const_pair_iterator {
    public:
      const_pair_iterator() : myNode(NULL) {}
      const_pair_iterator( const PairNode* node ) : myNode(node) {}
      const_pair_iterator( const const_iterator& i ) : myNode(i.mNode) {}
      
      const std::pair<MBEntityHandle, MBEntityHandle>& operator*() const
        { return *myNode; }
      
      const std::pair<MBEntityHandle, MBEntityHandle>* operator->() const
        { return myNode; }
      
      const_pair_iterator& operator--()
        { myNode = myNode->mPrev; return *this; }
      
      const_pair_iterator& operator++()
        { myNode = myNode->mNext; return *this; }
      
      const_pair_iterator operator--(int)
        { const_pair_iterator rval(*this); this->operator--(); return rval; }
      
      const_pair_iterator operator++(int)
        { const_pair_iterator rval(*this); this->operator++(); return rval; }
        
      bool operator==( const const_pair_iterator& other ) const
        { return other.myNode == myNode; }
    
      bool operator!=( const const_pair_iterator& other ) const
        { return other.myNode != myNode; }
        
    private:
      const PairNode* myNode;
  };
  
  pair_iterator pair_begin() { return pair_iterator(mHead.mNext); }
  pair_iterator pair_end() { return pair_iterator(&mHead); }

  const_pair_iterator const_pair_begin() const { return const_pair_iterator( mHead.mNext ); }
  const_pair_iterator const_pair_end() const { return const_pair_iterator( &mHead ); }
};

//! Use as you would an STL back_inserter
/**
 *  e.g. std::copy(list.begin(), list.end(), mb_range_inserter(my_range);
 * Also, see comments/instructions at the top of this class declaration
 */
class mb_range_inserter 
{
  
protected:
  MBRange* container;
 
public:
  //constructor
  explicit mb_range_inserter(MBRange& __x) : container(&__x) {}
  mb_range_inserter&
  operator=(const MBRange::value_type& __value) 
  {
    container->insert(__value);
    return *this;
  }

  mb_range_inserter& operator*() { return *this; }
  mb_range_inserter& operator++() { return *this; }
  mb_range_inserter& operator++(int) { return *this; }

  typedef MBEntityHandle            value_type;
  typedef MBEntityID                difference_type;
  typedef std::output_iterator_tag  iterator_category;
  typedef MBEntityHandle*           pointer;
  typedef MBEntityHandle&           reference;
};


inline MBRange::MBRange()
{
    // set the head node to point to itself
  mHead.mNext = mHead.mPrev = &mHead;
  mHead.first = mHead.second = 0;
}
  
  //! destructor
inline MBRange::~MBRange()
{
  clear();
}

  //! return the beginning const iterator of this range
inline MBRange::const_iterator MBRange::begin() const
{
  return const_iterator(mHead.mNext, mHead.mNext->first);
}
  
  //! return the beginning const reverse iterator of this range
inline MBRange::const_reverse_iterator MBRange::rbegin() const
{
  return const_reverse_iterator(mHead.mPrev, mHead.mPrev->second);
}
 
  //! return the ending const iterator for this range
inline MBRange::const_iterator MBRange::end() const
{
  return const_iterator(&mHead, mHead.first);
}
  
  //! return the ending const reverse iterator for this range
inline MBRange::const_reverse_iterator MBRange::rend() const
{
  return const_reverse_iterator(&mHead, mHead.second);
}

  //! return whether empty or not 
  //! always use "if(!Ranges::empty())" instead of "if(Ranges::size())"
inline bool MBRange::empty() const
{
  return (mHead.mNext == &mHead);
}

  //! remove a range of items from the list
inline MBRange::iterator MBRange::erase( iterator iter1, iterator iter2)
{
  while( iter1 != iter2 )
    erase( iter1++ );
  return iter1; 
}

  //! erases a value from this container
inline MBRange::iterator MBRange::erase(MBEntityHandle val) 
{ 
  return erase(find(val)); 
}
  
inline MBRange::const_iterator MBRange::const_iterator::end_of_block() const
  { return MBRange::const_iterator( mNode, mNode->second ); }

  //! get first entity in range
inline const MBEntityHandle& MBRange::front() const
  { return mHead.mNext->first; }
  //! get last entity in range
inline const MBEntityHandle& MBRange::back() const
  { return mHead.mPrev->second; }

inline std::ostream& operator<<( std::ostream& s, const MBRange& r )
  { r.print(s); return s; }
  
bool operator==( const MBRange& r1, const MBRange& r2 );
inline bool operator!=( const MBRange& r1, const MBRange& r2 )
  { return !(r1 == r2); }

inline MBEntityHandle MBRange::operator[](MBEntityID index) const
{
  MBRange::const_iterator i = begin();
  i += index;
  return *i;
}

#endif // MB_RANGE_HPP



