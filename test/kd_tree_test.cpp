#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include <math.h>
#include <assert.h>
#include <float.h>
#include <cstdio>

#include "TestUtil.hpp"

using namespace moab;

const unsigned INTERVALS = 4;
const unsigned DEPTH = 7; // 3*log2(INTERVALS)+1
const char* TAG_NAME = "TEST_DATA";

EntityHandle create_tree( AdaptiveKDTree& tool, unsigned depth, int intervals, Tag* tag_handle = 0 );
void validate_tree( AdaptiveKDTree& tool, EntityHandle root, int depth, double intervals );

void test_tree_create();
void test_leaf_merge();
void test_tree_readwrite();
void test_tree_delete();
void test_iterator_back();
void test_point_search();

int main(int argc, char **argv)
{
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif

  int err = RUN_TEST(test_tree_create);
  if (err)  // can't run other tests if can't create tree
    return 1;
  err += RUN_TEST(test_leaf_merge);
#ifdef HDF5_FILE
  err += RUN_TEST(test_tree_readwrite);
#endif
  err += RUN_TEST(test_tree_delete);
  err += RUN_TEST(test_iterator_back);
  err += RUN_TEST(test_point_search);
  
#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return err;
}

EntityHandle create_tree( AdaptiveKDTree& tool, unsigned depth, int intervals, Tag* tag_handle )
{
  // Create tree root
  ErrorCode err;
  EntityHandle root, leaf;
  const double tree_box_min_corner[] = { 0, 0, 0 };  
    // Make each leaf box be 1x1x1.
  const double tree_box_max_corner[] = { intervals, intervals, intervals };
  err = tool.create_tree( tree_box_min_corner, tree_box_max_corner, root );
  assert(!err);
  
  // Use iterator to create tree to fixed depth of DEPTH
  AdaptiveKDTreeIter iter;
  err = tool.get_tree_iterator( root, iter );
  CHECK_ERR(err);
  while(err == MB_SUCCESS) { 
    if (iter.depth() < depth) {
        // bisect leaves along alternating axes
      AdaptiveKDTree::Plane split;
      split.norm = iter.depth() % 3;  // alternate split axes;
      split.coord = 0.5 * (iter.box_min()[split.norm] + iter.box_max()[split.norm]);
      err = tool.split_leaf( iter, split ); // advances iter to first new leaf
      CHECK_ERR(err);
    }
      // if current leaf is at desired depth, advance to next one
    else {
      err = iter.step();
    }
  }
  CHECK(MB_ENTITY_NOT_FOUND == err);
  
  if (!tag_handle)
    return root;
  
  // define a tag to use to store integer values on tree leaves
  err = tool.moab()->tag_get_handle( TAG_NAME, 1, MB_TYPE_INTEGER, *tag_handle, MB_TAG_DENSE|MB_TAG_EXCL );
  CHECK_ERR(err);

  // iterate over tree setting data
  int counter = 0;
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // store integer value on leaf
    ++counter;
    leaf = iter.handle();
    err = tool.moab()->tag_set_data( *tag_handle, &leaf, 1, &counter );
    CHECK_ERR(err);
  }
  
  return root;
}

void validate_tree( AdaptiveKDTree& tool, EntityHandle root, unsigned depth, int intervals, Tag data )
{
  ErrorCode err;
  const double VOL = 1.0; // all leaves should be 1x1x1 boxes
  int val;

  // iterate over tree, verifying leaves 
  AdaptiveKDTreeIter iter;
  int counter = 0;
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // store integer value on leaf
    ++counter;
    EntityHandle leaf = iter.handle();
    CHECK(leaf != 0);
    CHECK_EQUAL(MBENTITYSET, TYPE_FROM_HANDLE(leaf));
     
    // check size of leaf
    const double* min = iter.box_min();
    const double* max = iter.box_max();
    double dims[] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] };
    double volume = dims[0] * dims[1] * dims[2];
    CHECK_REAL_EQUAL( VOL, volume, DBL_EPSILON );  
   
    // check depth of leaf
    CHECK_EQUAL( depth, iter.depth() );
    
    // check tag value on leaf
    err = tool.moab()->tag_get_data( data, &leaf, 1, &val );
    CHECK_ERR(err);
    CHECK_EQUAL(counter, val);
  }
    // check number of leaves
  const int num_leaves = intervals*intervals*intervals;
  CHECK_EQUAL( num_leaves, counter );
}

void test_tree_create()
{
  Tag tag;
  Core mb;
  AdaptiveKDTree tool(&mb);
  const EntityHandle root = create_tree( tool, DEPTH, INTERVALS, &tag );
  validate_tree( tool, root, DEPTH, INTERVALS, tag );
}

void test_leaf_merge()
{
  ErrorCode err;
  Core mb;
  AdaptiveKDTree tool(&mb);
  Tag data;
  const EntityHandle root = create_tree( tool, DEPTH, INTERVALS, &data );
  
  // reduce tree depth to DEPTH-1 by merging adjacent leaf pairs, 
  // make new "leaf" have smaller of two data values on original pair
  AdaptiveKDTreeIter iter;
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // get data for first leaf
    int data1;
    EntityHandle leaf = iter.handle();
    err = mb.tag_get_data( data, &leaf, 1, &data1 );
    CHECK_ERR(err);
    // tree traversal is always such that two leaves with same parent are consective
    err = iter.step();
    CHECK_ERR(err);
    // get data for sibling
    int data2;
    leaf = iter.handle();
    err = mb.tag_get_data( data, &leaf, 1, &data2 );
    CHECK_ERR(err);
    // as we stored increasing values, these had better be increasing
    CHECK_EQUAL( 1, data2 - data1 );
    // merge leaf pair (iter can be at either one)
    err = tool.merge_leaf( iter );  // changes iter to be new "merged" leaf
    CHECK_ERR(err);
    // store smaller of two values on new leaf
    leaf = iter.handle();
    err = mb.tag_set_data( data, &leaf, 1, &data1 );
    CHECK_ERR(err);
  }

  
  // Iterate over tree, verifying leaves and checking data
  // Initial leaves had volume of 1 : merged pairs of leaves so volume should now be 2.
  // Initial leaves were enumerated in order : merged pairs so new leaves should
  //   have data incrementing in steps of 2.
  int counter = 1;
  for (err = tool.get_tree_iterator( root, iter ); !err; err = iter.step()) {
    // store integer value on leaf
    int data1;
    EntityHandle leaf = iter.handle();
    err = mb.tag_get_data( data, &leaf, 1, &data1 );
    CHECK_ERR(err);
    CHECK_EQUAL( counter, data1 );
    counter += 2;
      
    // check size of leaf
    const double* min = iter.box_min();
    const double* max = iter.box_max();
    double dims[] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] };
    double volume = dims[0] * dims[1] * dims[2];
    CHECK_REAL_EQUAL( 2.0, volume, DBL_EPSILON );  
    
    // check depth of leaf
    CHECK_EQUAL( DEPTH-1, iter.depth() );
  }
}

void test_tree_readwrite()
{
  ErrorCode err;
  Tag tag;
  Core mb;
  AdaptiveKDTree tool(&mb);
  EntityHandle root = create_tree( tool, DEPTH, INTERVALS, &tag );
  
  // write to file
  err = mb.write_file( "tree.h5m" );
  CHECK_ERR(err);

  // clear everything
  mb.delete_mesh();
  
  // read tree from file
  err = mb.load_file( "tree.h5m" );
  remove("tree.h5m");
  CHECK_ERR(err);
 
  // get tag handle by name, because the handle may have changed
  err = mb.tag_get_handle( TAG_NAME, 1, MB_TYPE_INTEGER, tag );
  CHECK_ERR(err);

  // get root handle for tree
  Range range;
  err = tool.find_all_trees( range );
  assert(!err);
  assert(range.size() == 1);
  root = range.front(); // first (only) handle
  
  validate_tree( tool, root, DEPTH, INTERVALS, tag );
}

void test_tree_delete()
{
  ErrorCode err;
  Core mb;
  AdaptiveKDTree tool(&mb);
  Tag data;
  const EntityHandle root = create_tree( tool, DEPTH, INTERVALS, &data );
  
  err = tool.delete_tree( root );
  CHECK_ERR(err);
  
  Range ents;
  err = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &data, 0, 1, ents );
  CHECK_ERR(err);
  CHECK(ents.empty());
}

void test_iterator_back( )
{
  Core mb;
  AdaptiveKDTree tool(&mb);
  const EntityHandle root = create_tree( tool, DEPTH, INTERVALS );
  
  AdaptiveKDTreeIter iter;
  ErrorCode rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
  CartVect min( iter.box_min() );
  CartVect max( iter.box_max() );
  EntityHandle leaf = iter.handle();
  
    // going back from first location should fail.
  rval = iter.back();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // make sure iterator is valid
  CHECK_REAL_EQUAL( min[0], iter.box_min()[0], DBL_EPSILON );
  CHECK_REAL_EQUAL( min[1], iter.box_min()[1], DBL_EPSILON );
  CHECK_REAL_EQUAL( min[2], iter.box_min()[2], DBL_EPSILON );
  CHECK_REAL_EQUAL( max[0], iter.box_max()[0], DBL_EPSILON );
  CHECK_REAL_EQUAL( max[1], iter.box_max()[1], DBL_EPSILON );
  CHECK_REAL_EQUAL( max[2], iter.box_max()[2], DBL_EPSILON );
  CHECK_EQUAL( leaf, iter.handle() );
  
  while (MB_SUCCESS == iter.step()) {
      // Get values at current iterator location
    CartVect next_min( iter.box_min() );
    CartVect next_max( iter.box_max() );
    EntityHandle next_leaf = iter.handle();
  
      // step back to previous location
    rval = iter.back();
    CHECK_ERR(rval);
    
      // check expected values for previous location
    CHECK_REAL_EQUAL( min[0], iter.box_min()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( min[1], iter.box_min()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( min[2], iter.box_min()[2], DBL_EPSILON );
    CHECK_REAL_EQUAL( max[0], iter.box_max()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( max[1], iter.box_max()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( max[2], iter.box_max()[2], DBL_EPSILON );
    CHECK_EQUAL( leaf, iter.handle() );
    
      // advance iterator to 'current' location
    rval = iter.step();
    CHECK_ERR(rval);
    
      // check that iterator values are correct
    CHECK_REAL_EQUAL( next_min[0], iter.box_min()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( next_min[1], iter.box_min()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( next_min[2], iter.box_min()[2], DBL_EPSILON );
    CHECK_REAL_EQUAL( next_max[0], iter.box_max()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( next_max[1], iter.box_max()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( next_max[2], iter.box_max()[2], DBL_EPSILON );
    CHECK_EQUAL( next_leaf, iter.handle() );
   
      // store values for next iteration
    min = next_min;
    max = next_max;
    leaf = next_leaf;
  }
}

void test_point_search()
{
  Core mb;
  AdaptiveKDTree tool(&mb);
  const EntityHandle root = create_tree( tool, DEPTH, INTERVALS );
  
  ErrorCode rval;
  EntityHandle leaf;
  AdaptiveKDTreeIter iter, iter2;
  
    // points to test
  CartVect left( 0.5 );
  CartVect right( CartVect(INTERVALS) - left );
 
    // compare leaf search to iterator search
  rval = tool.leaf_containing_point( root, left.array(), leaf );
  CHECK_ERR(rval);
  rval = tool.leaf_containing_point( root, left.array(), iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( leaf, iter.handle() );
  
    // iterator should be at 'first' leaf 
  rval = tool.get_tree_iterator( root, iter2 );
  CHECK_ERR(rval);
  for (;;) {
    CHECK_EQUAL( iter.handle(), iter2.handle() );
    CHECK_EQUAL( iter.depth(), iter2.depth() );
    CHECK_REAL_EQUAL( iter.box_min()[0], iter2.box_min()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_min()[1], iter2.box_min()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_min()[2], iter2.box_min()[2], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_max()[0], iter2.box_max()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_max()[1], iter2.box_max()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_max()[2], iter2.box_max()[2], DBL_EPSILON );
    
    rval = iter2.step();
    if (MB_ENTITY_NOT_FOUND == rval)
      break;
    CHECK_ERR(rval);
    rval = iter.step();
    CHECK_ERR(rval);
  }
  
    // compare leaf search to iterator search
  rval = tool.leaf_containing_point( root, right.array(), leaf );
  CHECK_ERR(rval);
  rval = tool.leaf_containing_point( root, right.array(), iter );
  CHECK_ERR(rval);
  assert( iter.handle() == leaf );
  
    // iterator should be at 'last' leaf 
  rval = tool.get_last_iterator( root, iter2 );
  CHECK_ERR(rval);
  for (;;) {
    CHECK_EQUAL( iter.handle(), iter2.handle() );
    CHECK_EQUAL( iter.depth(), iter2.depth() );
    CHECK_REAL_EQUAL( iter.box_min()[0], iter2.box_min()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_min()[1], iter2.box_min()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_min()[2], iter2.box_min()[2], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_max()[0], iter2.box_max()[0], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_max()[1], iter2.box_max()[1], DBL_EPSILON );
    CHECK_REAL_EQUAL( iter.box_max()[2], iter2.box_max()[2], DBL_EPSILON );
    
    rval = iter2.back();
    if (MB_ENTITY_NOT_FOUND == rval)
      break;
    CHECK_ERR(rval);
    rval = iter.back();
    CHECK_ERR(rval);
  }
}

  

