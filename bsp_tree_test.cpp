#include "MBCore.hpp"
#include "TestUtil.hpp"
#include "MBBSPTree.hpp"

void test_set_plane();
void test_iterator();
void test_box_iterator();
void test_tree_create();
void test_box_tree_create();
void test_leaf_containing_point_bounded_tree();
void test_leaf_containing_point_unbounded_tree();
void test_merge_leaf();

int main()
{
  int failures = 0;
  
  failures += RUN_TEST( test_set_plane );
  failures += RUN_TEST( test_iterator );
  failures += RUN_TEST( test_box_iterator );
  failures += RUN_TEST( test_tree_create );
  failures += RUN_TEST( test_box_tree_create );
  failures += RUN_TEST( test_leaf_containing_point_bounded_tree );
  failures += RUN_TEST( test_leaf_containing_point_unbounded_tree );
  failures += RUN_TEST( test_merge_leaf );

  return failures;
}


// Make CHECK_EQUAL macro work for planes
void check_equal( const MBBSPTree::Plane& p1, 
                  const MBBSPTree::Plane& p2,
                  const char* exp1, 
                  const char* exp2,
                  int line,
                  const char* file )
{
  if (fabs(p1.norm[0]-p2.norm[0]) > 1e-6 ||
      fabs(p1.norm[1]-p2.norm[1]) > 1e-6 ||
      fabs(p1.norm[2]-p2.norm[2]) > 1e-6 ||
      fabs(p1.coeff  -p2.coeff  ) > 1e-6) {
    printf( "Equality Test Failed: %s == %s\n", exp1, exp2 );
    printf( "  at line %d of '%s'\n", line, file );
    printf( "  Expected: %f x + %f y + %f z + %f = 0\n", 
            p1.norm[0], p1.norm[1], p1.norm[2], p1.coeff );
    printf( "  Actual  : %f x + %f y + %f z + %f = 0\n", 
            p2.norm[0], p2.norm[1], p2.norm[2], p2.coeff );
    printf( "\n" );
    FLAG_ERROR;
  }    
}         

void test_set_plane()
{
  MBBSPTree::Plane p;
  const double points[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
  p.set( points[0], points[1], points[2] );
  CHECK_REAL_EQUAL( 0.0, p.distance( points[0] ), 1e-10 );
  CHECK_REAL_EQUAL( 0.0, p.distance( points[1] ), 1e-10 );
  CHECK_REAL_EQUAL( 0.0, p.distance( points[2] ), 1e-10 );
}

void test_iterator()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeIter iter;
    
    // create a depth-1 tree
  rval = tool.create_tree( root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // check initial state of iterator
  CHECK_EQUAL( &tool, iter.tool() );
  CHECK_EQUAL(  root, iter.handle() );
  CHECK_EQUAL(  1u,   iter.depth() );
  
    // should fail if at root
  MBBSPTree::Plane p;
  rval = iter.get_parent_split_plane( p );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // check that step past end returns correct value
  rval = iter.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // reset iterator
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // check that step past start returns correct value
  rval = iter.back();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // reset iterator
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // insert a single split plane
  rval = tool.split_leaf( iter, MBBSPTree::Plane(2,0,0,0) );
  CHECK_ERR(rval);
  
    // check initial location
  CHECK_EQUAL( 2u, iter.depth() );
  CHECK( iter.handle() != root );
  
    // create new iterators at left and right ends
  MBBSPTreeIter left, right;
  rval = tool.get_tree_iterator( root, left );
  CHECK_ERR(rval);
  rval = tool.get_tree_end_iterator( root, right );
  CHECK_ERR(rval);
  
    // compare post-split iterator to left
  CHECK_EQUAL( left.depth(), iter.depth() );
  CHECK_EQUAL( left.handle(), iter.handle() );
  
    // step to other leaf
  rval = iter.step();
  CHECK_ERR(rval);
  
    // check location
  CHECK_EQUAL( 2u, iter.depth() );
  CHECK( iter.handle() != root );
  
    // compare stepped iterator to right
  CHECK_EQUAL( right.depth(), iter.depth() );
  CHECK_EQUAL( right.handle(), iter.handle() );
  
    // step to back to first leaf
  rval = iter.back();
  CHECK_ERR(rval);
  
    // compare to left
  CHECK_EQUAL( left.depth(), iter.depth() );
  CHECK_EQUAL( left.handle(), iter.handle() );
  
    // check plane
    // should have unit normal
  left.get_parent_split_plane( p );
  CHECK_EQUAL( MBBSPTree::Plane(1,0,0,0), p );
  p.norm[0] = 11;
  right.get_parent_split_plane( p );
  CHECK_EQUAL( MBBSPTree::Plane(1,0,0,0), p );
  
    // check step past ends
  rval = left.back();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = right.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}

bool compare_hexes( const double expected[8][3], 
                    const double actual[8][3],
                    double epsilon )
{
    // for each of three possible rotations
  const int rotation_maps[3][8] = { { 0, 1, 2, 3, 4, 5, 6, 7 },
                                    { 3, 2, 6, 7, 0, 1, 5, 4 },
                                    { 7, 3, 2, 6, 4, 0, 1, 5 } };
  for (int i = 0; i < 3; ++i) {
      // compare for rotationts about axis from face 4 to 5
    for (int j = 0; j < 4; ++j) {
      bool match = true;
        // for each face vertex
      for (int k = 0; k < 4 && match; ++k) {
          // for each coordinate
        for (int d = 0; d < 3; ++d) {
          if (fabs(expected[(j+k)%4  ][d] - actual[rotation_maps[i][k  ]][d]) > epsilon
           || fabs(expected[(j+k)%4+4][d] - actual[rotation_maps[i][k+4]][d]) > epsilon) {
            match = false; 
            break;
          }
        }
      }
      
      if (match)
        return true;
    }
  }
  
  printf("Hex vertex copmarison failed.\n");
  printf("           Exected                         Actual\n");
  for (int v = 0; v < 8; ++v) {
    printf("% 9f % 9f % 9f   % 9f % 9f % 9f\n",
      expected[v][0], expected[v][1], expected[v][2],
      actual[v][0], actual[v][1], actual[v][2]);
  }
  
  
  return false;
}

static void aabox_corners( const double min[3], 
                           const double max[3],
                           double corners[8][3] )
{
  const double expt[8][3] = { { min[0], min[1], min[2] }, 
                              { max[0], min[1], min[2] }, 
                              { max[0], max[1], min[2] }, 
                              { min[0], max[1], min[2] },
                              { min[0], min[1], max[2] }, 
                              { max[0], min[1], max[2] }, 
                              { max[0], max[1], max[2] }, 
                              { min[0], max[1], max[2] } };
  memcpy( corners, expt, sizeof(expt) );
}

static void aabox_corners( double min_x, double min_y, double min_z,
                           double max_x, double max_y, double max_z,
                           double corners[8][3] )
{
  double min[] = { min_x, min_y, min_z };
  double max[] = { max_x, max_y, max_z };
  aabox_corners( min, max, corners );
}
  

void test_box_iterator()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeBoxIter iter;
  const double min[3] = { -1, -2, -3 };
  const double max[3] = {  3,  2,  1 };
    
    // create a depth-1 tree
  rval = tool.create_tree( min, max, root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // check initial state of iterator
  CHECK_EQUAL( &tool, iter.tool() );
  CHECK_EQUAL(  root, iter.handle() );
  CHECK_EQUAL(  1u,   iter.depth() );
  
    // check initial corner values
  double corners[8][3], expt[8][3];
  aabox_corners( min, max, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
    // should fail if at root
  MBBSPTree::Plane p;
  rval = iter.get_parent_split_plane( p );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // check that step past end returns correct value
  rval = iter.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // reset iterator
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // check that step past start returns correct value
  rval = iter.back();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // reset iterator
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // insert a single split plane
  rval = tool.split_leaf( iter, MBBSPTree::Plane(2,0,0,0) );
  CHECK_ERR(rval);
  
    // check initial location
  CHECK_EQUAL( 2u, iter.depth() );
  CHECK( iter.handle() != root );
  
    // create new iterators at left and right ends
  MBBSPTreeIter left, right;
  rval = tool.get_tree_iterator( root, left );
  CHECK_ERR(rval);
  rval = tool.get_tree_end_iterator( root, right );
  CHECK_ERR(rval);
  
    // compare post-split iterator to left
  CHECK_EQUAL( left.depth(), iter.depth() );
  CHECK_EQUAL( left.handle(), iter.handle() );
  
    // check box
  aabox_corners( min, max, expt );
  for (int i = 0; i < 8; ++i)
    if (expt[i][0] > 0)
      expt[i][0] = 0;
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
    // step to other leaf
  rval = iter.step();
  CHECK_ERR(rval);
  
    // check location
  CHECK_EQUAL( 2u, iter.depth() );
  CHECK( iter.handle() != root );
  
    // compare stepped iterator to right
  CHECK_EQUAL( right.depth(), iter.depth() );
  CHECK_EQUAL( right.handle(), iter.handle() );
  
    // check box
  aabox_corners( min, max, expt );
  for (int i = 0; i < 8; ++i)
    if (expt[i][0] < 0)
      expt[i][0] = 0;
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
    // step to back to first leaf
  rval = iter.back();
  CHECK_ERR(rval);
  
    // compare to left
  CHECK_EQUAL( left.depth(), iter.depth() );
  CHECK_EQUAL( left.handle(), iter.handle() );
  
    // check box
  aabox_corners( min, max, expt );
  for (int i = 0; i < 8; ++i)
    if (expt[i][0] > 0)
      expt[i][0] = 0;
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
    // check plane
    // should have unit normal
  left.get_parent_split_plane( p );
  CHECK_EQUAL( MBBSPTree::Plane(1,0,0,0), p );
  p.norm[0] = 11;
  right.get_parent_split_plane( p );
  CHECK_EQUAL( MBBSPTree::Plane(1,0,0,0), p );
  
    // check step past ends
  rval = left.back();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = right.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}

void test_tree_create()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeIter iter;
  MBBSPTree::Plane p;
    
    // create a depth-1 tree
  rval = tool.create_tree( root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // check initial state of iterator
  CHECK_EQUAL( &tool, iter.tool() );
  CHECK_EQUAL(  root, iter.handle() );
  CHECK_EQUAL(  1u,   iter.depth() );

    // split with z=0
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,0,0.5,0) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,0,1,0), p );

    // split lower leaf with diagonal plane
  rval = tool.split_leaf( iter, MBBSPTree::Plane(1,1,0,0) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3u, iter.depth() );
  
  const double r = 0.5*sqrt(2.0);
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(r,r,0,0), p );
  
    // step to upper leaf
  rval = iter.step();
  CHECK_ERR(rval);
  
    // split upper leaf with diagonal plane
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,1,0,0) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(-r,r,0,0), p );
  
    // iterate over four leaves
  
    // first leaf
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(r,r,0,0), p );
  
    // second leaf
  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 4u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(-r,r,0,0), p );
   
    // third leaf
  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 4u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(-r,r,0,0), p );
   
    // fourth leaf
  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 2u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,0,1,0), p );
 
    // no more leaves
  rval = iter.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}

void test_box_tree_create()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeBoxIter iter;
  MBBSPTree::Plane p;
  const double min[3] = { -2, -2, -2 };
  const double max[3] = {  2,  2,  2 };
    
    // create a depth-1 tree
  rval = tool.create_tree( min, max, root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // check initial state of iterator
  CHECK_EQUAL( &tool, iter.tool() );
  CHECK_EQUAL(  root, iter.handle() );
  CHECK_EQUAL(  1u,   iter.depth() );
  
    // check initial corner values
  double corners[8][3], expt[8][3];
  aabox_corners( min, max, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

 
    // Try splits that should fail because they
    // do not intersect the leaf at all
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,0,1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,0,1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,1,0, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,1,0,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(1,0,0, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(1,0,0,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1,-1, 7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1,-1,-7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1,-1, 7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1,-1,-7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1,-1, 7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1,-1,-7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1,-1, 7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1,-1,-7) );
  CHECK_EQUAL( MB_FAILURE, rval );
  
 
    // Try a split that should fail because the
    // resulting leaf would not be a hexahedron.
    // Clip each corner twice using planes with opposing normals
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1,-1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1, 1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1,-1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1, 1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1,-1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1, 1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1,-1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1, 1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1, 1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1,-1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1, 1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1,-1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1, 1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1,-1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1, 1,-4) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1,-1, 4) );
  CHECK_EQUAL( MB_FAILURE, rval );
    // Clip each edge
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1,-1, 0,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1,-1, 0,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 1, 0,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 1, 0,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 0,-1,-1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 0, 1,-1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 0, 1, 1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 0,-1, 1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 0,-1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 0,-1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane( 1, 0, 1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-1, 0, 1,-2) );
  CHECK_EQUAL( MB_FAILURE, rval );
 
    // verify that iterator is still valid after many failed splits
  CHECK_EQUAL( &tool, iter.tool() );
  CHECK_EQUAL(  root, iter.handle() );
  CHECK_EQUAL(  1u,   iter.depth() );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );


    // split with z=0
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,0,0.5,0) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,0,1,0), p );
  
    // check that box corners are correct
  aabox_corners( min, max, expt );
  for (unsigned i = 0; i < 8; ++i)
    if (expt[i][2] > 0.0)
      expt[i][2] = 0.0;
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

  
    // split with z=-1 and normal in opposite direction
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,0,-2,-2) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3u, iter.depth() );
  for (unsigned i = 0; i < 8; ++i)
    if (expt[i][2] < -1.0)
      expt[i][2] = -1.0;
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
    // step to next leaf (z < -1)
  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 3u, iter.depth() );
  aabox_corners( min, max, expt );
  for (unsigned i = 0; i < 8; ++i)
    if (expt[i][2] > -1.0)
      expt[i][2] = -1.0;
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

    
    // split at x = -1
  rval = tool.split_leaf( iter, MBBSPTree::Plane(-0.1,0,0,-0.1) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 4u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(-1,0,0,-1), p );
  
    // check that leaf box is correct
  aabox_corners( -1, -2, -2, 2, 2, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );


    // split at x = 1
  rval = tool.split_leaf( iter, MBBSPTree::Plane(5,0,0,-5) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 5u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(1,0,0,-1), p );
  
    // check that leaf box is correct
  aabox_corners( -1, -2, -2, 1, 2, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

    
    // split at y = -1
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,-1,0,-1) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 6u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,-1,0,-1), p );
  
    // check that leaf box is correct
  aabox_corners( -1, -1, -2, 1, 2, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );


    // split at y = 1
  rval = tool.split_leaf( iter, MBBSPTree::Plane(0,1,0,-1) );
  CHECK_ERR(rval);
  CHECK_EQUAL( 7u, iter.depth() );
  
    // check that parent split plane is correct
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,1,0,-1), p );
  
    // check that leaf box is correct
  aabox_corners( -1, -1, -2, 1, 1, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
  
  
    // iterate over tree, verifying 
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,0,-1,-1), p );
  aabox_corners( -2, -2, -1, 2, 2, 0, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 7u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,1,0,-1), p );
  aabox_corners( -1, -1, -2, 1, 1, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 7u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,1,0,-1), p );
  aabox_corners( -1, 1, -2, 1, 2, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 6u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,-1,0,-1), p );
  aabox_corners( -1, -2, -2, 1, -1, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 5u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(1,0,0,-1), p );
  aabox_corners( 1, -2, -2, 2, 2, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );

  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 4u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(-1,0,0,-1), p );
  aabox_corners( -2, -2, -2, -1, 2, -1, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
      
  rval = iter.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( 2u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBBSPTree::Plane(0,0,1,0), p );
  aabox_corners( -2, -2, 0, 2, 2, 2, expt );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expt, corners, 1e-6 ) );
  
    // no more leaves
  rval = iter.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}

void test_leaf_containing_point_bounded_tree()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeIter iter;
  MBBSPTreeBoxIter b_iter;
  MBBSPTree::Plane p;
  MBEntityHandle h;
  double expected[8][3], corners[8][3];


/*  Build Tree

  1.0 +---------+--------------+
      |    A    |              |
      |         |              |
  0.7 +---------+      C       |
      |         |              |
      |         |              |
      |    B    |              |
      |         +--------------+ 0.3
      |         |      D       |
      |         |              |
  0.0 +---------+--------------+
      0.0       0.4            1.0  */

  const MBBSPTree::Plane  X_split(1.0, 0.0, 0.0,-0.4);
  const MBBSPTree::Plane AB_split(0.0,-1.0, 0.0, 0.7);
  const MBBSPTree::Plane CD_split(0.0,-1.0, 0.0, 0.3);
  
  
  const double min[3] = { 0, 0, 0 };
  const double max[3] = { 1, 1, 1 };
  rval = tool.create_tree( min, max, root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
  rval = tool.split_leaf( iter, X_split );
  CHECK_ERR(rval);
  
  rval = tool.split_leaf( iter, AB_split );
  CHECK_ERR(rval);
  const MBEntityHandle A = iter.handle();
  
  rval = iter.step();
  CHECK_ERR(rval);
  const MBEntityHandle B = iter.handle();
  
  rval = iter.step();
  CHECK_ERR(rval);
  rval = tool.split_leaf( iter, CD_split );
  CHECK_ERR(rval);
  const MBEntityHandle C = iter.handle();
  
  rval = iter.step();
  CHECK_ERR(rval);
  const MBEntityHandle D = iter.handle();

  
    // Test queries inside tree

  const double A_point[] = { 0.1, 0.8, 0.5 };
  rval = tool.leaf_containing_point( root, A_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( A, h );
  rval = tool.leaf_containing_point( root, A_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( A, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( AB_split, p );
  rval = tool.leaf_containing_point( root, A_point, b_iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( A, b_iter.handle() );
  aabox_corners( 0.0, 0.7, 0.0, 0.4, 1.0, 1.0, expected );
  rval = b_iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expected, corners, 1e-6 ) );

  const double B_point[] = { 0.3, 0.1, 0.6 };
  rval = tool.leaf_containing_point( root, B_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( B, h );
  rval = tool.leaf_containing_point( root, B_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( B, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( AB_split, p );
  rval = tool.leaf_containing_point( root, B_point, b_iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( B, b_iter.handle() );
  aabox_corners( 0.0, 0.0, 0.0, 0.4, 0.7, 1.0, expected );
  rval = b_iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expected, corners, 1e-6 ) );

  const double C_point[] = { 0.9, 0.9, 0.1 };
  rval = tool.leaf_containing_point( root, C_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( C, h );
  rval = tool.leaf_containing_point( root, C_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( C, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( CD_split, p );
  rval = tool.leaf_containing_point( root, C_point, b_iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( C, b_iter.handle() );
  aabox_corners( 0.4, 0.3, 0.0, 1.0, 1.0, 1.0, expected );
  rval = b_iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expected, corners, 1e-6 ) );

  const double D_point[] = { 0.5, 0.2, 0.9 };
  rval = tool.leaf_containing_point( root, D_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( D, h );
  rval = tool.leaf_containing_point( root, D_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( D, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( CD_split, p );
  rval = tool.leaf_containing_point( root, D_point, b_iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( D, b_iter.handle() );
  aabox_corners( 0.4, 0.0, 0.0, 1.0, 0.3, 1.0, expected );
  rval = b_iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expected, corners, 1e-6 ) );
  
  
    // Try a couple points outside of the tree 

  const double above_pt[] = { 0.5, 0.5, 2.0 };
  rval = tool.leaf_containing_point( root, above_pt, b_iter );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );

  const double left_pt[] = { -1.0, 0.5, 0.5 };
  rval = tool.leaf_containing_point( root, left_pt, b_iter );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}

void test_leaf_containing_point_unbounded_tree()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeIter iter;
  MBBSPTree::Plane p;
  MBEntityHandle h;


/*  Build Tree

    \      |
     \  C  o (0,4,0)
      \    |
       \   |
        \  |
     B   \ |         D
          \|
   ________o (0,0,0)
            \
             \
       A      \
               o (1,-2,0)
                \
 */
  const MBBSPTree::Plane  X_split( 2.0, 1.0, 0.0, 0.0);
  const MBBSPTree::Plane AB_split( 0.0, 1.0, 0.0, 0.0);
  const MBBSPTree::Plane CD_split( 1.0, 0.0, 0.0, 0.0);
  
  
  rval = tool.create_tree( root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
  rval = tool.split_leaf( iter, X_split );
  CHECK_ERR(rval);
  
  rval = tool.split_leaf( iter, AB_split );
  CHECK_ERR(rval);
  const MBEntityHandle A = iter.handle();
  
  rval = iter.step();
  CHECK_ERR(rval);
  const MBEntityHandle B = iter.handle();
  
  rval = iter.step();
  CHECK_ERR(rval);
  rval = tool.split_leaf( iter, CD_split );
  CHECK_ERR(rval);
  const MBEntityHandle C = iter.handle();
  
  rval = iter.step();
  CHECK_ERR(rval);
  const MBEntityHandle D = iter.handle();

  
    // Test queries inside tree

  const double A_point[] = { -1000, -1000, -1000 };
  rval = tool.leaf_containing_point( root, A_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( A, h );
  rval = tool.leaf_containing_point( root, A_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( A, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( AB_split, p );

  const double B_point[] = { -3, 1, 100 };
  rval = tool.leaf_containing_point( root, B_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( B, h );
  rval = tool.leaf_containing_point( root, B_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( B, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( AB_split, p );

  const double C_point[] = { -10, 500, 0 };
  rval = tool.leaf_containing_point( root, C_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( C, h );
  rval = tool.leaf_containing_point( root, C_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( C, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( CD_split, p );

  const double D_point[] = { 10, 500, 0 };
  rval = tool.leaf_containing_point( root, D_point, h );
  CHECK_ERR(rval);
  CHECK_EQUAL( D, h );
  rval = tool.leaf_containing_point( root, D_point, iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( D, iter.handle() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( CD_split, p );
}

void test_merge_leaf()
{
  MBCore moab;
  MBBSPTree tool( &moab );
  MBErrorCode rval;
  MBEntityHandle root;
  MBBSPTreeBoxIter iter;
  MBBSPTree::Plane p;
  double expected[8][3], corners[8][3];


/*  Build Tree

  1.0 +---------+--------------+
      |    A    |              |
      |         |              |
  0.7 +---------+      C       |
      |         |              |
      |         |              |
      |    B    |              |
      |         +--------------+ 0.3
      |         |      D       |
      |         |              |
  0.0 +---------+--------------+
      0.0       0.4            1.0  */

  const MBBSPTree::Plane  X_split(1.0, 0.0, 0.0,-0.4);
  const MBBSPTree::Plane AB_split(0.0,-1.0, 0.0, 0.7);
  const MBBSPTree::Plane CD_split(0.0,-1.0, 0.0, 0.3);
  
  const double min[3] = { 0, 0, 0 };
  const double max[3] = { 1, 1, 1 };
  rval = tool.create_tree( min, max, root );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  rval = tool.split_leaf( iter, X_split );
  CHECK_ERR(rval);
  const MBEntityHandle AB = iter.handle();
  rval = tool.split_leaf( iter, AB_split );
  CHECK_ERR(rval);
  rval = iter.step();
  CHECK_ERR(rval);
  rval = iter.step();
  CHECK_ERR(rval);
  const MBEntityHandle CD = iter.handle();
  rval = tool.split_leaf( iter, CD_split );
  CHECK_ERR(rval);
  rval = iter.step();
  CHECK_ERR(rval);

  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  rval = tool.merge_leaf( iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( AB, iter.handle() );
  CHECK_EQUAL( 2u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( X_split, p );
  aabox_corners( 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, expected );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expected, corners, 1e-6 ) );
  
  rval = iter.step();
  CHECK_ERR(rval);
  rval = iter.step();
  CHECK_ERR(rval);
  rval = tool.merge_leaf( iter );
  CHECK_ERR(rval);
  CHECK_EQUAL( CD, iter.handle() );
  CHECK_EQUAL( 2u, iter.depth() );
  rval = iter.get_parent_split_plane( p );
  CHECK_ERR(rval);
  CHECK_EQUAL( X_split, p );
  aabox_corners( 0.4, 0.0, 0.0, 1.0, 1.0, 1.0, expected );
  rval = iter.get_box_corners( corners );
  CHECK_ERR( rval );
  CHECK( compare_hexes( expected, corners, 1e-6 ) );
  
  rval = iter.step();
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}

