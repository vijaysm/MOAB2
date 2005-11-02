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


#ifdef DEBUG_TEST
#include "RangeTree.hpp"
#include <assert.h>
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
