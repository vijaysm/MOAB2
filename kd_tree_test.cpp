#include "MBCore.hpp"
#include "MBAdaptiveKDTree.hpp"

#include <math.h>
#include <assert.h>
#include <float.h>

const unsigned INTERVALS = 4;
const unsigned DEPTH = 7; // 3*log2(INTERVALS)+1
const char* TAG_NAME = "TEST_DATA";

int main()
{
  // Initialize MOAB & create tree tool
  MBCore moab;
  MBAdaptiveKDTree tool(&moab);
  
  // Create tree root
  MBErrorCode err;
  MBEntityHandle root, leaf;
  const double tree_box_min_corner[] = { 0, 0, 0 };  
    // Make each leaf box be 1x1x1.
  const double tree_box_max_corner[] = { INTERVALS, INTERVALS, INTERVALS };
  err = tool.create_tree( tree_box_min_corner, tree_box_max_corner, root );
  assert(!err);
  
  // Use iterator to create tree to fixed depth of DEPTH
  MBAdaptiveKDTreeIter iter;
  err = tool.get_tree_iterator( root, iter );
  assert(!err);
  while(err == MB_SUCCESS) { 
    if (iter.depth() < DEPTH) {
        // bisect leaves along alternating axes
      MBAdaptiveKDTree::Plane split;
      split.norm = iter.depth() % 3;  // alternate split axes;
      split.coord = 0.5 * (iter.box_min()[split.norm] + iter.box_max()[split.norm]);
      err = tool.split_leaf( iter, split ); // advances iter to first new leaf
      assert(!err);
    }
      // if current leaf is at desired depth, advance to next one
    else {
      err = iter.step();
    }
  }
  assert(MB_ENTITY_NOT_FOUND == err);
  
  // define a tag to use to store integer values on tree leaves
  MBTag data;
  err = moab.tag_create( TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, data, 0, false );
  assert(!err);
  
  // iterate over tree, verifying leaves and setting data
  int counter = 0;
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // store integer value on leaf
    ++counter;
    leaf = iter.handle();
    err = moab.tag_set_data( data, &leaf, 1, &counter );
    assert(!err);
      
    // check size of leaf
    const double* min = iter.box_min();
    const double* max = iter.box_max();
    double dims[] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] };
    double volume = dims[0] * dims[1] * dims[2];
    assert( fabs(volume - 1.0) <= DBL_EPSILON );  
    
    // check depth of leaf
    assert( iter.depth() == DEPTH );
  }
    // check number of leaves
  const int num_leaves = INTERVALS*INTERVALS*INTERVALS;
  assert( num_leaves == counter );
  
  // reduce tree depth to DEPTH-1 by merging adjacent leaf pairs, 
  // make new "leaf" have smaller of two data values on original pair
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // get data for first leaf
    int data1;
    leaf = iter.handle();
    err = moab.tag_get_data( data, &leaf, 1, &data1 );
    assert(!err);
    // tree traversal is always such that two leaves with same parent are consective
    err = iter.step();
    assert(!err);
    // get data for sibling
    int data2;
    leaf = iter.handle();
    err = moab.tag_get_data( data, &leaf, 1, &data2 );
    assert(!err);
    // as we stored increasing values, these had better be increasing
    assert( data2 - data1 == 1 );
    // merge leaf pair (iter can be at either one)
    err = tool.merge_leaf( iter );  // changes iter to be new "merged" leaf
    assert(!err);
    // store smaller of two values on new leaf
    leaf = iter.handle();
    err = moab.tag_set_data( data, &leaf, 1, &data1 );
    assert(!err);
  }
  
  // write to file
#if MB_VERSION_MAJOR <= 3
  err = moab.write_mesh( "tree.h5m" );
#else
  err = moab.write_file( "tree.h5m" );
#endif
  assert(!err);

  // clear everything
  moab.delete_mesh();
  
  // read tree from file
#if MB_VERSION_MAJOR <= 3
  err = moab.load_mesh( "tree.h5m" );
#else
  MBEntityHandle file;
  err = moab.load_file( "tree.h5m", file );
#endif
  assert(!err);
 
  // get tag handle by name, because the handle may have changed
  err = moab.tag_get_handle( TAG_NAME, data );
  assert(!err);

  // get root handle for tree
  MBRange range;
  err = tool.find_all_trees( range );
  assert(!err);
  assert(range.size() == 1);
  root = range.front(); // first (only) handle
  
  // Iterate over tree, verifying leaves and checking data
  // Initial leaves had volume of 1 : merged pairs of leaves so volume should now be 2.
  // Initial leaves were enumerated in order : merged pairs so new leaves should
  //   have data incrementing in steps of 2.
  counter = 1;
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // store integer value on leaf
    int data1;
    leaf = iter.handle();
    err = moab.tag_get_data( data, &leaf, 1, &data1 );
    assert(!err);
    assert( counter == data1 );
    counter += 2;
      
    // check size of leaf
    const double* min = iter.box_min();
    const double* max = iter.box_max();
    double dims[] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] };
    double volume = dims[0] * dims[1] * dims[2];
    assert( fabs(volume - 2.0) <= DBL_EPSILON );  
    
    // check depth of leaf
    assert( iter.depth() == DEPTH-1 );
  }
    // check number of leaves
    // (num_leaves is original number of leaves, twice current number,
    //  but counter is incremented by 2 for each iteration, so just
    //  compare them.)
  assert( counter-1 == num_leaves );
  
  
    // Delete data from tree
  err = moab.tag_delete( data );
  assert( !err );

  
  // write to file
#if MB_VERSION_MAJOR <= 3
  err = moab.write_mesh( "tree.h5m" );
#else
  err = moab.write_file( "tree.h5m" );
#endif
  assert(!err);

  // clear everything
  moab.delete_mesh();
  
  // read tree from file
#if MB_VERSION_MAJOR <= 3
  err = moab.load_mesh( "tree.h5m" );
#else
  MBEntityHandle file;
  err = moab.load_file( "tree.h5m", file );
#endif
  assert(!err);
  
  // check that tag doesn't exist
  err = moab.tag_get_handle( TAG_NAME, data );
  assert( MB_TAG_NOT_FOUND == err );

  return 0;
}

