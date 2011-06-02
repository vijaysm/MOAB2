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

// MOAB performance tests building mapped mesh with nodes and
// hexes created one at a time.  This also creates the node to hex adjacencies.

#include <math.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <assert.h>
#include <list>
#include "moab/Core.hpp"
#include "moab/Skinner.hpp"
#include "moab/ReadUtilIface.hpp"

using namespace moab;

double LENGTH = 1.0;
const int DEFAULT_INTERVALS = 50;

void create_regular_mesh( int interval, int dimension );
void skin_common( int interval, int dim, int blocks, bool use_adj );
void skin( int intervals, int dim, int num ) 
  { std::cout << "Skinning w/out adjacencies:" << std::endl;
    skin_common( intervals, dim, num, false ); }
void skin_adj( int intervals, int dim, int num )
  { std::cout << "Skinning with adjacencies:" << std::endl;
    skin_common( intervals, dim, num, true ); }

void tag_time( TagType storage, bool direct, int intervals, int dim, int blocks );

void dense_tag( int intervals, int dim, int blocks ) 
  { std::cout << "Dense Tag Time:"; 
    tag_time( MB_TAG_DENSE, false, intervals, dim, blocks ); }
void sparse_tag( int intervals, int dim, int blocks )
  { std::cout << "Sparse Tag Time:"; 
    tag_time( MB_TAG_SPARSE, false, intervals, dim, blocks ); }
void direct_tag( int intervals, int dim, int blocks )
  { std::cout << "Direct Tag Time:"; 
    tag_time( MB_TAG_DENSE, true, intervals, dim, blocks ); }

typedef void (*test_func_t)( int, int, int );
const struct {
  std::string testName;
  test_func_t testFunc;
  std::string testDesc;
} TestList[] = {
 { "skin",      &skin,      "Test time to get skin mesh w/out adjacencies" },
 { "skin_adj",  &skin_adj,  "Test time to get skin mesh with adjacencies" },
 { "sparse",    &sparse_tag,"Sparse tag data manipulation" },
 { "dense",     &dense_tag, "Dense tag data manipulation" },
 { "direct",    &direct_tag,"Dense tag data manipulation using direct data access" },
};
const int TestListSize = sizeof(TestList)/sizeof(TestList[0]); 

void usage( const char* argv0, bool error = true )
{
  std::ostream& str = error ? std::cerr : std::cout;
  str << "Usage: " << argv0 << " [-i <ints_per_side>] [-d <dimension>] [-n <test_specifc_int>] <test_name> [<test_name2> ...]" << std::endl;
  str << "       " << argv0 << " [-h|-l]" << std::endl;
  if (error)
    return;
  str << "  -i  : specify interverals per side (num hex = ints^3, default: " << DEFAULT_INTERVALS << std::endl;
  str << "  -d  : specify element dimension, default: 3" << std::endl;
  str << "  -n  : specify an integer value that for which the meaning is test-specific" << std::endl;
  str << "  -h  : print this help text." << std::endl;
  str << "  -l  : list available tests"  << std::endl;
}

void list_tests( ) 
{
  unsigned max_test_name = 0, max_test_desc = 0;
  for (int i = 0; i < TestListSize; ++i) {
    if (TestList[i].testName.size() > max_test_name)
      max_test_name = TestList[i].testName.size();
    if (TestList[i].testDesc.size() > max_test_desc)
      max_test_desc = TestList[i].testDesc.size();
  }
  std::cout << std::setw(max_test_name) << "NAME" << "   " 
            << std::setw(max_test_desc) << std::left << "DESCRIPTION" << std::endl
            << std::setfill('-') << std::setw(max_test_name) << "" << "   " 
            << std::setfill('-') << std::setw(max_test_desc) << "" 
            << std::setfill(' ') << std::endl;
  for (int i = 0; i < TestListSize; ++i) 
    std::cout << std::setw(max_test_name) << TestList[i].testName << " : "
              << std::setw(max_test_desc) << std::left << TestList[i].testDesc << std::endl;
}  
  
int main(int argc, char* argv[])
{
  int intervals = DEFAULT_INTERVALS;
  int dimension = 3;
  int number = 0;
  std::vector<test_func_t> test_list;
  std::list<int*> expected_list;
  bool did_help = false;
  for (int i = 1; i < argc; ++i) {
    if (!expected_list.empty()) {
      int* ptr = expected_list.front();
      expected_list.pop_front();
      char* endptr;
      *ptr = strtol( argv[i], &endptr, 0 );
      if (!endptr) {
        usage(argv[0]);
        std::cerr << "Expected integer value, got \"" << argv[i] << '"' << std::endl;
        return 1;
      }
    }
    else if (*argv[i] == '-') { // flag
      for (int j = 1; argv[i][j]; ++j) {
        switch (argv[i][j]) {
          case 'i': expected_list.push_back( &intervals ); break;
          case 'd': expected_list.push_back( &dimension ); break;
          case 'n': expected_list.push_back( &number    ); break;
          case 'h': did_help = true; usage(argv[0],false); break;
          case 'l': did_help = true; list_tests(); break;
          default:
            usage(argv[0]);
            std::cerr << "Invalid option: -" << argv[i][j] << std::endl;
            return 1;
        }
      }
    }
    else {
      int j = -1;
      do {
        ++j;
        if (j >= TestListSize) {
          usage(argv[0]);
          std::cerr << "Invalid test name: " << argv[i] << std::endl;
          return 1;
        }
      } while (TestList[j].testName != argv[i]);
      test_list.push_back( TestList[j].testFunc );
    }
  }
  
  if (!expected_list.empty()) {
    usage(argv[0]);
    std::cerr << "Missing final argument" <<std::endl;
    return 1;
  }
  
    // error if no input
  if (test_list.empty() && !did_help) {
    usage(argv[0]);
    std::cerr << "No tests specified" << std::endl;
    return 1;
  }

  if (intervals < 1) {
    std::cerr << "Invalid interval count: " << intervals << std::endl;
    return 1;
  }
  
  if (dimension < 1 || dimension > 3) {
    std::cerr << "Invalid dimension: " << dimension << std::endl;
    return 1;
  }
  
    // now run the tests
  for (std::vector<test_func_t>::iterator i = test_list.begin(); i != test_list.end(); ++i) {
    test_func_t fptr = *i;
    fptr(intervals,dimension,number);
  }
  
  return 0;
}

void create_regular_mesh( Interface* gMB, int interval, int dim )
{
  if (dim < 1 || dim > 3 || interval < 1) {
    std::cerr << "Invalid arguments" << std::endl;
    exit(1);
  }
  
  const int nvi = interval+1;
  const int dims[3] = { nvi, dim > 1 ? nvi : 1, dim > 2 ? nvi : 1 };
  int num_vert = dims[0] * dims[1] * dims[2];

  ReadUtilIface* readMeshIface;
  gMB->query_interface(readMeshIface);
  
  EntityHandle vstart;
  std::vector<double*> arrays;
  ErrorCode rval = readMeshIface->get_node_coords(3, num_vert, 1, vstart, arrays);
  if (MB_SUCCESS != rval || arrays.size() < 3) {
    std::cerr << "Vertex creation failed" << std::endl;
    exit(2);
  }
  double *x = arrays[0], *y = arrays[1], *z = arrays[2];
  
    // Calculate vertex coordinates
  for (int k = 0; k < dims[2]; ++k)
    for (int j = 0; j < dims[1]; ++j)
      for (int i = 0; i < dims[0]; ++i)
      {
        *x = i; ++x;
        *y = j; ++y;
        *z = k; ++z;
      }
  
  const long vert_per_elem = 1 << dim; // 2^dim
  const long intervals[3] = { interval, dim>1?interval:1, dim>2?interval:1 };
  const long num_elem = intervals[0]*intervals[1]*intervals[2];
  const EntityType type = (dim == 1) ? MBEDGE : (dim == 2) ? MBQUAD : MBHEX;
  
  EntityHandle estart, *conn = 0;
  rval = readMeshIface->get_element_connect( num_elem, vert_per_elem, type, 0, estart, conn );
  if (MB_SUCCESS != rval || !conn) {
    std::cerr << "Element creation failed" << std::endl;
    exit(2);
  }
  
  
    // Offsets of element vertices in grid relative to corner closest to origin 
  long c = dims[0]*dims[1];
  const long corners[8] = { 0, 1, 1+dims[0], dims[0], c, c+1, c+1+dims[0], c+dims[0] };
                             
    // Populate element list
  EntityHandle* iter = conn;
  for (long z = 0; z < intervals[2]; ++z)
    for (long y = 0; y < intervals[1]; ++y)
      for (long x = 0; x < intervals[0]; ++x)
      {
        const long index = x + y*dims[0] + z*(dims[0]*dims[1]);
        for (long j = 0; j < vert_per_elem; ++j, ++iter)
          *iter = index + corners[j] + vstart;
      }
  
    // notify MOAB of the new elements
  rval = readMeshIface->update_adjacencies(estart, num_elem, vert_per_elem, conn);
  if (MB_SUCCESS != rval) {
    std::cerr << "Element update failed" << std::endl;
    exit(2);
  }
}  

void skin_common( int interval, int dim, int num, bool use_adj ) 
{
  Core moab;
  Interface* gMB = &moab;
  ErrorCode rval;
  double d;
  clock_t t, tt;

  create_regular_mesh( gMB, interval, dim );

  Range skin, verts, elems;
  rval = gMB->get_entities_by_dimension( 0, dim, elems );
  assert(MB_SUCCESS == rval); assert(!elems.empty());

  Skinner tool(gMB);
  
  t = clock();
  rval = tool.find_skin( elems, true, verts, 0, use_adj, false );
  t = clock() - t;
  if (MB_SUCCESS != rval) {
    std::cerr << "Search for skin vertices failed" << std::endl;
    exit(2);
  }
  d = ((double)t)/CLOCKS_PER_SEC;
  std::cout << "Got " << verts.size() << " skin vertices in " << d << " seconds." << std::endl;
  
  t = 0;
  if (num < 1) num = 1000;
  long blocksize = elems.size() / num;
  if (!blocksize) blocksize = 1;
  long numblocks = elems.size()/blocksize;
  Range::iterator it = elems.begin();
  for (long i = 0; i < numblocks; ++i) {
    verts.clear();
    Range::iterator end = it + blocksize;
    Range blockelems;
    blockelems.merge( it, end );
    it = end;
    tt = clock();
    rval = tool.find_skin( blockelems, true, verts, 0, use_adj, false );
    t += clock() - tt;
    if (MB_SUCCESS != rval) {
      std::cerr << "Search for skin vertices failed" << std::endl;
      exit(2);
    }
  }
  d = ((double)t)/CLOCKS_PER_SEC;
  std::cout << "Got skin vertices for " << numblocks << " blocks of " 
            << blocksize << " elements in " << d << " seconds." << std::endl;
 
  for (int e = 0; e < 2; ++e) { // do this twice 
    if (e == 1) {
        // create all interior faces
      skin.clear();
      t = clock();
      gMB->get_adjacencies( elems, dim-1, true, skin, Interface::UNION );
      t = clock() - t;
      d = ((double)t)/CLOCKS_PER_SEC;
      std::cout << "Created " << skin.size() << " entities of dimension-1 in " << d << " seconds" << std::endl;
    }
  
    skin.clear();
    t = clock();
    rval = tool.find_skin( elems, false, skin, 0, use_adj, true );
    t = clock() - t;
    if (MB_SUCCESS != rval) {
      std::cerr << "Search for skin vertices failed" << std::endl;
      exit(2);
    }
    d = ((double)t)/CLOCKS_PER_SEC;
    std::cout << "Got " << skin.size() << " skin elements in " << d << " seconds." << std::endl;

    t = 0;
    it = elems.begin();
    for (long i = 0; i < numblocks; ++i) {
      skin.clear();
      Range::iterator end = it + blocksize;
      Range blockelems;
      blockelems.merge( it, end );
      it = end;
      tt = clock();
      rval = tool.find_skin( blockelems, false, skin, 0, use_adj, true );
      t += clock() - tt;
      if (MB_SUCCESS != rval) {
        std::cerr << "Search for skin elements failed" << std::endl;
        exit(2);
      }
    }
    d = ((double)t)/CLOCKS_PER_SEC;
    std::cout << "Got skin elements for " << numblocks << " blocks of " 
              << blocksize << " elements in " << d << " seconds." << std::endl;
  }
}

void tag_time( TagType storage, bool direct, int intervals, int dim, int blocks )
{
  Core moab;
  Interface& mb = moab;
  create_regular_mesh( &mb, intervals, dim );
  
    // Create tag in which to store data
  Tag tag;
  mb.tag_get_handle( "data", 1, MB_TYPE_DOUBLE, tag, storage|MB_TAG_CREAT );
  
    // Make up some arbitrary iterative calculation for timing purposes:
    // set each value v_n = (V + v_n)/2 until all values are within
    // epsilon of V.
  std::vector<double> data;
  Range verts;
  mb.get_entities_by_type( 0, MBVERTEX, verts );
  
  clock_t t = clock();
  
    // initialize
  if (direct) {
    Range::iterator i, j = verts.begin();
    void* ptr;
    int count;
    while (j != verts.end()) {
      mb.tag_iterate( tag, i, verts.end(), count, ptr );
      double* arr = reinterpret_cast<double*>(ptr);
      for (j = i + count; i != j; ++i, ++arr) 
        *arr = (11.0 * *i + 7.0)/(*i);
    }
  }
  else {
    data.resize( verts.size() );
    double* arr = &data[0];
    for (Range::iterator i = verts.begin(); i != verts.end(); ++i, ++arr)
      *arr = (11.0 * *i + 7.0)/(*i);
    mb.tag_set_data( tag, verts, &data[0] );
  }
  
    // iterate
  const double v0 = acos(-1.0); // pi
  size_t iter_count = 0;
  double max_diff;
  const size_t num_verts = verts.size();
  do {
    if (direct) {
      max_diff = 0.0;
      Range::iterator i, j = verts.begin();
      void* ptr;
      while (j != verts.end()) {
        int count;
        mb.tag_iterate( tag, i, verts.end(), count, ptr );
        double* arr = reinterpret_cast<double*>(ptr);
        
        for (j = i+count; i != j; ++i, ++arr) {
          *arr = 0.5 * (*arr + v0);
          double diff = fabs(*arr - v0);
          if (diff > max_diff)
            max_diff = diff;
        }
      }
    }
    else if (blocks < 1) {
      max_diff = 0.0;
      mb.tag_get_data( tag, verts, &data[0] );
      for (size_t i = 0; i < data.size(); ++i) { 
        data[i] = 0.5 * (v0+data[i]);
        double diff = fabs( data[i] - v0 );
        if (diff > max_diff)
          max_diff = diff;
      }
      mb.tag_set_data( tag, verts, &data[0] );
    }
    else {
      max_diff = 0.0;
      Range r;
      Range::iterator it = verts.begin();
      size_t step = num_verts / blocks;
      for (int j = 0; j < blocks-1; ++j) {
        Range::iterator nx = it;
        nx += step;
        r.clear();
        r.merge( it, nx );
        mb.tag_get_data( tag, r, &data[0] );
        it = nx;

        for (size_t i = 0; i < step; ++i) { 
          data[i] = 0.5 * (v0+data[i]);
          double diff = fabs( data[i] - v0 );
          if (diff > max_diff)
            max_diff = diff;
        }
        mb.tag_set_data( tag, r, &data[0] );
      }
      
      r.clear();
      r.merge( it, verts.end() );
      mb.tag_get_data( tag, r, &data[0] );
      for (size_t i = 0; i < (num_verts - (blocks-1)*step); ++i) { 
        data[i] = 0.5 * (v0+data[i]);
        double diff = fabs( data[i] - v0 );
        if (diff > max_diff)
          max_diff = diff;
      }
      mb.tag_set_data( tag, r, &data[0] );
    }
    ++iter_count;
//    std::cout << iter_count << " " << max_diff << std::endl;
  } while (max_diff > 1e-6);
  
  double secs = (clock() - t) / (double)CLOCKS_PER_SEC;
  std::cout << " " << iter_count << " iterations in " << secs << " seconds" << std::endl;
}
