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
  advance iterator
*/
MBRange::const_iterator& MBRange::const_iterator::operator+=( long sstep )
{
    // Check negative now to avoid infinite loop below.
  if (sstep < 0)
  {
    return operator-=( -sstep );
  }
  unsigned long step = sstep;
  
    // Handle current PairNode.  Either step is within the current
    // node or need to remove the remainder of the current node
    // from step.
  MBEntityHandle this_node_rem = mNode->second - mValue;
  if (this_node_rem >= step)
  {
    mValue += step;
    return *this;
  }
  step -= this_node_rem + 1;

    // For each node we are stepping past, decrement step
    // by the size of the node.
  PairNode* node = mNode->mNext;
  MBEntityHandle node_size = node->second - node->first + 1;
  while (step >= node_size)
  {
    step -= node_size;
    node = node->mNext;
    node_size = node->second - node->first + 1;
  }
  
    // Advance into the resulting node by whatever is
    // left in step.
  mNode = node;
  mValue = mNode->first + step;
  return *this;
}
  
    
 
/*!
  regress iterator
*/
MBRange::const_iterator& MBRange::const_iterator::operator-=( long sstep )
{
    // Check negative now to avoid infinite loop below.
  if (sstep < 0)
  {
    return operator+=( -sstep );
  }
  unsigned long step = sstep;
  
    // Handle current PairNode.  Either step is within the current
    // node or need to remove the remainder of the current node
    // from step.
  MBEntityHandle this_node_rem = mValue - mNode->first;
  if (this_node_rem >= step)
  {
    mValue -= step;
    return *this;
  }
  step -= this_node_rem + 1;

    // For each node we are stepping past, decrement step
    // by the size of the node.
  PairNode* node = mNode->mPrev;
  MBEntityHandle node_size = node->second - node->first + 1;
  while (step >= node_size)
  {
    step -= node_size;
    node = node->mPrev;
    node_size = node->second - node->first + 1;
  }
  
    // Advance into the resulting node by whatever is
    // left in step.
  mNode = node;
  mValue = mNode->second - step;
  return *this;
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

  return insert( begin(), val1, val2 );
}

MBRange::iterator MBRange::insert( MBRange::iterator prev,
                                   MBEntityHandle val1, 
                                   MBEntityHandle val2 )
{
  assert( val1 <= val2  );

  // Empty 
  if (mHead.mNext == &mHead)
  {
    assert( prev == end() );
    PairNode* new_node = new PairNode( &mHead, &mHead, val1, val2 );
    mHead.mNext = mHead.mPrev = new_node;
    return iterator( mHead.mNext, val1 );
  }
  
  PairNode* iter = prev.mNode;
  assert( iter != &mHead );
  assert( iter->mPrev == &mHead || iter->mPrev->second+1 < val1 );
  
  // Input range is before beginning?
  if (iter->mPrev == &mHead && val2 < iter->first - 1)
  {
    PairNode* new_node = new PairNode( iter, &mHead,  val1, val2 );
    mHead.mNext = iter->mPrev = new_node;
    return iterator( mHead.mNext, val1 );
  }
  
  // Find first intersecting list entry, or the next entry
  // if none intersects.
  while (iter != &mHead && iter->second+1 < val1)
    iter = iter->mNext;
  
  // Need to insert new pair (don't intersect any existing pair)?
  if (iter == &mHead || iter->first-1 > val1)
  {
    PairNode* new_node = new PairNode( iter, iter->mPrev, val1, val2 );
    iter->mPrev = iter->mPrev->mNext = new_node;
    return iterator( iter->mPrev, val1 );
  }
  
  // Make the first intersecting pair the union of itself with [val1,val2]
  if (iter->first > val1)
    iter->first = val1;
  if (iter->second >= val2)  
    return iterator( iter, val1 );
  iter->second = val2;
  
  // Merge any remaining pairs that intersect [val1,val2]
  while (iter->mNext != &mHead && iter->mNext->first <= val2 + 1)
  {
    PairNode* dead = iter->mNext;
    iter->mNext = dead->mNext;
    dead->mNext->mPrev = iter;
    
    if (dead->second > val2)
      iter->second = dead->second;
    delete dead;
  }
  
  return iterator( iter, val1 );
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
  merge( range.begin(), range.end() );
}

void MBRange::merge( MBRange::const_iterator begin,
                     MBRange::const_iterator end )
{
  if (begin == end)
    return;
  
  PairNode* node = begin.mNode;
  if (end.mNode == node)
  {
    insert( *begin, (*end)-1 );
    return;
  }
  
  MBRange::iterator hint = insert( *begin, node->second );
  node = node->mNext;
  while (node != end.mNode)
  {
    hint = insert( hint, node->first, node->second );
    node = node->mNext;
  }
  
  if (*end > node->first)
  {
    if (*end <= node->second)
      insert( hint, node->first, *(end) - 1 );
    else
      insert( hint, node->first, node->second );
  }
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
  PairNode* iter;
  for (iter = first.mNode; iter != last.mNode; iter = iter->mNext)
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
  
  if (iter->first >= val)
    return const_iterator( iter, iter->first );
  else if(*last > val)
    return const_iterator( iter, val );
  else
    return last;
}
