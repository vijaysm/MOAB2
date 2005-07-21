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


#if !defined(RANGE_TREE_RECURSIVE_BLOCK) && (!defined(RANGE_TREE_HPP) || defined(LINUX) || defined(WIN32))
#define RANGE_TREE_RECURSIVE_BLOCK

#ifndef RANGE_TREE_HPP
#  include "RangeTree.hpp"
#endif
#include <assert.h>

template <class T> typename RangeTree<T>::const_iterator& 
RangeTree<T>::const_iterator::operator++()
{
  if (span_ == end_)
    return *this;
  
  if (value_ == span_->last)
  {
    ++span_;
    value_ = (span_ == end_) ? 0 : span_->first;
  }
  else
  { 
    ++value_;
  }
  
  return *this;
}


template <class T> typename RangeTree<T>::const_iterator& 
RangeTree<T>::const_iterator::operator--()
{
  if (span_ == end_)
  {
    --span_;
    value_ = span_->last;
  }
  
  else if (value_ == span_->first && span_ != beg_)
  {
    --span_;
    value_ = span_->first;
  }
  
  else
  { 
    --value_;
  }
  
  return *this;
}

template <class T> typename RangeTree<T>::const_iterator
RangeTree<T>::find( T value ) const
{
  Span span = { value, value };
  typename SpanTree::iterator iter = set_.find( span );
  if (iter == set_.end())
    return end();
  else
    return const_iterator( iter, set_.end(), value );  
}

template <class T> typename RangeTree<T>::const_iterator
RangeTree<T>::insert( T value )
{
  Span span = { value, value };
  typename SpanTree::iterator iter = set_.find( span );
  if (iter != set_.end()) 
    return const_iterator(iter, set_.end(), value);
  
  --span.first;
  iter = set_.find( span );
  if (iter != set_.end())
  {
    typename SpanTree::iterator next = iter;
    ++next;
    if (next != set_.end() && next->first - 1 == value )
    {
      T tmp = next->last;
      set_.erase( next );
      iter->last = tmp;
    }
    else
    {
      ++iter->last;
    }
    return const_iterator(iter, set_.end(), value);
  }
  ++span.first;
  
  ++span.last;
  iter = set_.find( span );
  if (iter != set_.end())
  {
    iter->first--;
    return const_iterator(iter, set_.end(), value);
  }
  --span.last;
  
  iter = set_.insert( span ).first;
  return const_iterator(iter, set_.end(), value);
}

template <class T> T RangeTree<T>::size( ) const
{
  T count = 0;
  for (typename SpanTree::iterator iter = set_.begin(); iter != set_.end(); ++iter)
    count += iter->last - iter->first + 1;
  return count;
}

template <class T> bool RangeTree<T>::empty() const
{
  return set_.empty();
}





#ifdef DEBUG_TEST
#include <iostream>
#include <set>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void print(const RangeTree<int>& tree, bool brief)
{
  int rcount = 0;
  int vcount = 0;
  RangeTree<int>::span_iterator iter = tree.span_begin();
  const RangeTree<int>::span_iterator end = tree.span_end();
  
  for ( ; iter != end ; ++iter)
  {
    if (!brief)
      std::cout << (*iter).first << '-' << (*iter).second << ' ';
    rcount++;
    vcount += (*iter).second - (*iter).first + 1;
  }
  if (!brief)
    std::cout << std::endl;
  std::cout << rcount << " ranges, " << vcount << " values" << std::endl;
}

void validate(const RangeTree<int>& tree)
{
  RangeTree<int>::span_iterator iter = tree.span_begin();
  const RangeTree<int>::span_iterator end = tree.span_end();
  if (iter == end)
    return;
    
  for(;;)
  {
    long prev = (*iter).second;
    ++iter;
    if (iter == end)
      break;
    long next = (*iter).first;
    assert(next - prev > 1);
  }
}


int main( int argc, char* argv[] )
{
  RangeTree<int> tree;
  int size = 0;
  
    /* Random mode - "-r" optionally followed by the number
                     of random ranges to generate.
                   - Inserts random ranges of random length [1,64]
                   - Inserts values if each range in random order.
                   - Keeps std::set to compare to rangetree when done.
    */
  if (argc > 1 && strcmp( argv[1], "-r") == 0)
  {
    RangeTree<int>::iterator t_iter;
    std::set<int>::iterator s_iter;
    std::vector<int>::iterator v_iter;
    char marks[64];
    int count = argc > 2 ? atoi(argv[2]) : 100;
    if (count < 1)
      count = 100;
    int vallimit = sizeof(marks) * count;
    std::vector<int> vallist(0);
    
    std::set<int> values;

    clock_t clk = clock();
    std::cout << "Generating random values... ";
    std::cout.flush();
    while (count--)
    {
      int start = rand() % vallimit;  // first value in range to insert
      int length = rand() % sizeof(marks) + 1;  // length of range to insert
      memset( marks, 0, sizeof(marks) );  // values inserted
      while (length)  // insert range values in random order
      {
        int k = rand() % length + 1;
        int n = 0;
        
        while( k-- )
          while (marks[n]) 
            n++;
        
        marks[n] = 1;
        length--;
        n += start;
        vallist.push_back(n);
      }
    }
    double t = clock() - clk;
    t /= CLOCKS_PER_SEC;
    std::cout << t << " seconds." << std::endl;
    
    std::cout << "Building tree... ";
    std::cout.flush();
    clk = clock();
    for (v_iter = vallist.begin(); v_iter != vallist.end(); ++v_iter)
      tree.insert( *v_iter );
    t = clock() - clk;
    t /= CLOCKS_PER_SEC;
    std::cout << t << " seconds." << std::endl;
    
    std::cout << "Validating tree... ";
    std::cout.flush();
    clk = clock();
    validate(tree);
    t = clock() - clk;
    t /= CLOCKS_PER_SEC;
    std::cout << t << " seconds." << std::endl;
    print(tree,true);

    std::cout << "Building set... ";
    std::cout.flush();
    clk = clock();
    for (v_iter = vallist.begin(); v_iter != vallist.end(); ++v_iter)
      values.insert( *v_iter );
    t = clock() - clk;
    t /= CLOCKS_PER_SEC;
    std::cout << t << " seconds." << std::endl;
    
    std::cout << "Comparing tree to set... ";
    std::cout.flush();
    clk = clock();
    assert (values.size() == tree.size());
    t_iter = tree.begin();
    s_iter = values.begin();
    while (s_iter != values.end())
    {
      assert( t_iter != tree.end() );    
      assert( *t_iter == *s_iter );
      ++t_iter;
      ++s_iter;
    }
    t = clock() - clk;
    t /= CLOCKS_PER_SEC;
    std::cout << t << " seconds." << std::endl;
  }
    /* Literal mode - expects args to be individual values or
                      ranges of the form n-m.
                    - prints resulting tree
    */
  else
  {
    for (int i = 1; i < argc; i++ )
    {
      int start, end, count;
      count = sscanf( argv[i], "%d-%d", &start, &end );
      if (count == 1)
      {
        tree.insert( start );
        size++;
      }
      else if (count == 2 && start < end)
      {
        for ( ; start <= end; start++ )
        {
          tree.insert( start );
          size++;
        }
      }
      else
      {
        fprintf(stderr, "Ingoring invalid argument: \"%s\"\n", argv[i]);
      }
    }

    validate(tree);
    print(tree,false);
  }
  return 0;
}

#endif
#endif
