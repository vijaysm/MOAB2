/****************************************************
 * File     :      MBRange.cpp
 *
 * Purpose  :      Stores contiguous or partially
 *                 contiguous values in an optimized
 *                 fashion.  Partially contiguous
 *                 accessing patterns is also optimized.
 *
 * Creator  :      Clinton Stimpson
 *
 * Date     :      15 April 2002
 *
 *******************************************************/


#include <assert.h>
#include "MBRange.hpp"
#include <stdio.h>

/*! 
  returns the number of values this list represents
 */
unsigned int MBRange::size() const
{
  // go through each pair and add up the number of values
  // we have.
  unsigned int size=0;
  for(PairNode* iter = mHead.mNext; iter != &mHead; iter = iter->mNext)
  {
    size += ((iter->second - iter->first) + 1);
  }
  return size;
}

/*!
  inserts a single value into this range
*/

MBRange::iterator MBRange::insert(MBEntityHandle val)
{

  // if this is empty, just add it and return an iterator to it
  if(val == 0)
    return end();

  if(&mHead == mHead.mNext)
  {
    mHead.mNext = mHead.mPrev = new PairNode(&mHead, &mHead, val, val);
    return iterator(mHead.mNext, val);
  }
  
  // find the location in the list where we can safely insert
  // new items and keep it ordered
  PairNode* jter = mHead.mNext;
  for( ; (jter != &mHead) && (jter->second < val); jter=jter->mNext);
  PairNode* iter = jter;
  jter = jter->mPrev;

  // if this val is already in the list
  if( (iter->first <= val && iter->second >= val) && (iter != &mHead) )
  {
    // return an iterator pointing to this location
    return iterator( iter, val );
  }

  // one of a few things can happen at this point
  // 1. this range needs to be backwardly extended
  // 2. the previous range needs to be forwardly extended
  // 3. a new range needs to be added

  
  // extend this range back a bit
  else if( (iter->first == (val+1)) && (iter != &mHead) )
  {
    iter->first = val;
    // see if we need to merge two ranges
    if( (iter != mHead.mNext) && (jter->second == (val-1)))
    {
      jter->second = iter->second;
      iter->mPrev->mNext = iter->mNext;
      iter->mNext->mPrev = iter->mPrev;
      delete iter;
      return iterator( jter, val );
    }
    else
    {
      return iterator( iter, val );
    }

  }
  // extend the previous range forward a bit
  else if( (jter->second == (val-1)) && (iter != mHead.mNext) )
  {
    jter->second = val;
    return iterator(jter, val);
  }
  // make a new range
  else
  {
    PairNode* new_node = new PairNode(iter, iter->mPrev, val, val);
    iter->mPrev = new_node->mPrev->mNext = new_node;
    return iterator(new_node, val);
  }

}


/*!
  inserts a range of values
*/
MBRange::iterator MBRange::insert(MBEntityHandle val1, MBEntityHandle val2)
{

  if(val1 == 0 || val1 > val2)
    return end();

  // if this is empty, just add it and return an iterator to it
  if(&mHead == mHead.mNext)
  {
    mHead.mNext = mHead.mPrev = new PairNode(&mHead, &mHead, val1, val2);
    return iterator(mHead.mNext, val1);
  }

  // find the location in the list where we can safely begin insertions
  // and keep our list in order
  PairNode* jter = mHead.mNext;
  for( ; (jter != &mHead) && (jter->first < val1); jter = jter->mNext);

  PairNode* kter=NULL;

  if((jter->mPrev->second+1) >= val1)
  {
    jter->mPrev->second = val2;
    kter = jter->mPrev;
  }
  else
  {
    // just insert this range
    kter = new PairNode(jter, jter->mPrev, val1, val2);
    jter->mPrev->mNext = kter;
    jter->mPrev = kter;
  }

  jter = kter->mPrev;
  PairNode* iter = kter->mNext;
  // at this point, kter is the newly inserted node, 
  // jter is the previous after and iter is the one after
  
  // merge forward as many times as needed
  while( (kter->second >= iter->first - 1) && (iter != &mHead))
  {
    if(kter->second < iter->second)
      kter->second = iter->second;
    kter->mNext = iter->mNext;
    iter->mNext->mPrev = kter;
    delete iter;
    iter = kter->mNext;
  }

  // be sure something didn't go wrong
  assert( (kter->first <= val1) && (val1 <= kter->second) );
  return iterator( kter, val1 );
  
}



/*!
  erases an item from this list and returns an iterator to the next item
*/

MBRange::iterator MBRange::erase(iterator iter)
{
  // one of a few things could happen
  // 1. shrink a range
  // 2. split a range
  // 3. remove a range

  if(iter == end())
    return end();

  // the iterator most likely to be returned
  iterator new_iter = iter;
  ++new_iter;

  PairNode* kter = iter.mNode;
  
  // just remove the range
  if(kter->first == kter->second)
  {
    kter->mNext->mPrev = kter->mPrev;
    kter->mPrev->mNext = kter->mNext;
    delete kter;
    return new_iter;
  }
  // shrink it
  else if(kter->first == iter.mValue)
  {
    kter->first++;
    return new_iter;
  }
  // shrink it the other way
  else if(kter->second == iter.mValue)
  {
    kter->second--;
    return new_iter;
  }
  // split the range
  else
  {
    PairNode* new_node = new PairNode(iter.mNode, iter.mNode->mPrev, kter->first, iter.mValue-1);
    new_node->mPrev->mNext = new_node->mNext->mPrev = new_node;
    iter.mNode->first = iter.mValue+1;
    return new_iter;
  }

}

/*!
  finds a value in the list.
  this method is preferred over other algorithms because
  it can be found faster this way.
*/
MBRange::iterator MBRange::find(MBEntityHandle val)
{
  // iterator through the list
  PairNode* iter = mHead.mNext;
  for( ; iter != &mHead && (val > iter->second); iter=iter->mNext );
  return ((iter->second >= val) && (iter->first <= val)) ? iterator(iter,val) : end();
}

/*!
  finds a value in the list.
  this method is preferred over other algorithms because
  it can be found faster this way.
*/
MBRange::const_iterator MBRange::find(MBEntityHandle val) const
{
  // iterator through the list
  PairNode* iter = mHead.mNext;
  for( ; iter != &mHead && (val > iter->second); iter=iter->mNext );
  return ((iter->second >= val) && (iter->first <= val)) ? const_iterator(iter,val) : end();
}

/*!
  merges another MBRange with this one
*/


void MBRange::merge( const MBRange& range )
{

  if(range.empty())
    return;

  PairNode* iter = range.mHead.mPrev;
  
  //insert all the base ranges
  for(; iter != &(range.mHead); iter=iter->mPrev )
    insert(iter->first, iter->second);

}

#include <algorithm>


// checks the range to make sure everything is A-Ok.
void MBRange::sanity_check() const
{
  if(empty())
    return;

  const PairNode* node = mHead.mNext;
  std::vector<const PairNode*> seen_before;
  bool stop_it = false;
  
  for(; stop_it == false; node = node->mNext)
  {
    // have we seen this node before?
    assert(std::find(seen_before.begin(), seen_before.end(), node) == seen_before.end());
    seen_before.push_back(node);

    // is the connection correct?
    assert(node->mNext->mPrev == node);

    // are the values right?
    assert(node->first <= node->second);
    if(node != &mHead && node->mPrev != &mHead)
      assert(node->mPrev->second < node->first);

    if(node == &mHead)
      stop_it = true;

  }

}

// for debugging
void MBRange::print() const
{
  printf("List contents:\n");
  for(PairNode* iter = mHead.mNext; iter != &mHead; iter = iter->mNext)
  {
    printf("%u   %u\n", iter->first, iter->second);
  }
}

  // intersect two ranges, placing the results in the return range
#define MAX(a,b) (a < b ? b : a)
#define MIN(a,b) (a > b ? b : a)
MBRange MBRange::intersect(const MBRange &range2) const 
{
  pair_iterator r_it[2] = {pair_iterator(begin()), pair_iterator(range2.begin())};
  MBEntityHandle low_it, high_it;
  
  MBRange lhs;
  
    // terminate the while loop when at least one "start" iterator is at the
    // end of the list
  while (r_it[0] != end() && r_it[1] != range2.end()) {
    
    if (r_it[0]->second < r_it[1]->first)
        // 1st subrange completely below 2nd subrange
      r_it[0]++;
    else if (r_it[1]->second < r_it[0]->first) 
        // 2nd subrange completely below 1st subrange
      r_it[1]++;
    
    else {
        // else ranges overlap; first find greater start and lesser end
      low_it = MAX(r_it[0]->first, r_it[1]->first);
      high_it = MIN(r_it[0]->second, r_it[1]->second);
      
        // insert into result
      lhs.insert(low_it, high_it);
      
        // now find bounds of this insertion and increment corresponding iterator
      if (high_it == r_it[0]->second) r_it[0]++;
      if (high_it == r_it[1]->second) r_it[1]++;
    }
  }
  
  return lhs;
}

MBRange::const_iterator MBRange::lower_bound(MBRange::const_iterator first,
                                             MBRange::const_iterator last,
                                             MBEntityHandle val)
{
    // Find the first pair whose end is >= val
  for (PairNode* iter = first.mNode; iter != last.mNode; iter = iter->mNext)
  {
    if (iter->second >= val)
    {
        // This is the correct pair.  Either 'val' is in the range, or
        // the range starts before 'val' and iter->first IS the lower_bound.
      if (iter->first > val)
        return const_iterator(iter, iter->first);
      return const_iterator(iter, val);
    }
  }
  
  return last;
}

MBRange::iterator MBRange::lower_bound(MBRange::iterator first,
                                       MBRange::iterator last,
                                       MBEntityHandle val)
{
    // Find the first pair whose end is >= val
  for (PairNode* iter = first.mNode; iter != last.mNode; iter = iter->mNext)
  {
    if (iter->second >= val)
    {
        // This is the correct pair.  Either 'val' is in the range, or
        // the range starts before 'val' and iter->first IS the lower_bound.
      if (iter->first > val)
        return iterator(iter, iter->first);
      return iterator(iter, val);
    }
  }
  
  return last;
}

