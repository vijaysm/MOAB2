#include "MBCore.hpp"
#include "testdir.h"
#include "TestUtil.hpp"
#include "MBRange.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

const char filename[] = "sets.h5m";
bool keep_file = false;

void read_write_file( MBInterface& output, MBInterface& input, MBEntityHandle* input_set = 0 )
{
  MBEntityHandle file;
  MBErrorCode rval;
  rval = output.write_file( filename );
  CHECK_ERR(rval);
  rval = input.load_file( filename, file );
  if (!keep_file)
    remove(filename);
  CHECK_ERR(rval);
  if (input_set)
    *input_set = file;
}

void test_ranged_set_with_holes()
{
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  MBRange verts;
  
  const int num_vtx = 40;
  std::vector<double> coords( 3*num_vtx, 0.0 );
  rval = mb.create_vertices( &coords[0], num_vtx, verts );
  CHECK_ERR(rval);
  CHECK_EQUAL(num_vtx, (int)verts.size());
  
  MBEntityHandle set;
  rval = mb.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  rval = mb.add_entities( set, verts );
  
  std::vector<MBEntityHandle> dead_verts;
  for (int i = num_vtx/4; i < num_vtx; i += num_vtx/4 ) {
    MBRange::iterator j = verts.begin();
    j += i;
    dead_verts.push_back( *j );
  }
  rval = mb.delete_entities( &dead_verts[0], dead_verts.size() );
  CHECK_ERR(rval);
  
  MBCore moab2;
  MBInterface& mb2 = moab2;
  MBEntityHandle file_set;
  read_write_file( mb, mb2, &file_set );
  MBRange sets;
  mb2.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_EQUAL( 2, (int)sets.size() );
  MBEntityHandle other_set = sets.front() == file_set ? sets.back() : sets.front();
  
  int num_vtx2 = -5;
  rval = mb2.get_number_entities_by_type( other_set, MBVERTEX, num_vtx2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( (int)(num_vtx - dead_verts.size()), num_vtx2 );
}

void test_file_set()
{
  MBErrorCode rval;
  MBCore moab;
  double vtxcoords[] = { 0.0, 0.0, 0.0, 
                         1.0, 0.0, 0.0, 
                         0.0, 1.0, 0.0 };
  MBRange verts;
  rval = moab.create_vertices( vtxcoords, 3, verts );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, (int)verts.size() );
  
  MBEntityHandle tri;
  MBEntityHandle conn[3];
  std::copy( verts.begin(), verts.end(), conn );
  rval = moab.create_element( MBTRI, conn, 3, tri );
  CHECK_ERR(rval);
  
  MBEntityHandle set;
  rval = moab.create_meshset( MESHSET_ORDERED, set );
  CHECK_ERR(rval);
  rval = moab.add_entities( set, &tri, 1 );
  CHECK_ERR(rval);
  
  MBEntityHandle file;
  read_write_file( moab, moab, &file );
  
  int count;
  rval = moab.get_number_entities_by_type( 0, MBVERTEX, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 6, count );
  rval = moab.get_number_entities_by_type( file, MBVERTEX, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, count );
  
  rval = moab.get_number_entities_by_type( 0, MBTRI, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, count );
  rval = moab.get_number_entities_by_type( file, MBTRI, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, count );
  
  rval = moab.get_number_entities_by_type( 0, MBENTITYSET, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, count );
  rval = moab.get_number_entities_by_type( file, MBENTITYSET, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, count );
}


int coords_by_idx( int idx, double coords[][3] )
{
  coords[0][0] = idx;
  coords[0][1] = 0;
  coords[0][2] = 0;
  coords[1][0] = 0;
  coords[1][1] = idx;
  coords[1][2] = 0;
  coords[2][0] = 0;
  coords[2][1] = 0;
  coords[2][2] = idx;
  coords[3][0] = 3.14*idx;
  coords[3][1] = 1;
  coords[3][2] = 1;
  coords[4][0] = 1;
  coords[4][1] = 3.14*idx;
  coords[4][2] = 1;
  coords[5][0] = 1;
  coords[5][1] = 1;
  coords[5][2] = 3.14*idx;
  return idx % 5  + 1;
}

void recursive_build_tree( int max_depth,
                           MBInterface& mb,
                           MBTag tag,
                           MBEntityHandle p,
                           int depth,
                           int& idx )
{
  MBErrorCode rval = mb.tag_set_data( tag, &p, 1, &idx ); CHECK_ERR(rval);
  
  MBRange verts;
  double coords[6][3];
  int num_vtx = coords_by_idx( idx, coords );
  rval = mb.create_vertices( &coords[0][0], num_vtx, verts );
  rval = mb.add_entities( p, verts );
  ++idx;
  if (depth == max_depth)
    return;

  MBEntityHandle l, r;
  rval = mb.create_meshset( MESHSET_SET, l ); CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, r ); CHECK_ERR(rval);
  rval = mb.add_parent_child( p, l ); CHECK_ERR(rval);
  rval = mb.add_parent_child( p, r ); CHECK_ERR(rval);
  
  recursive_build_tree( max_depth, mb, tag, l, depth+1, idx );
  recursive_build_tree( max_depth, mb, tag, r, depth+1, idx );
}
 
void recursive_check_tree( int max_depth,
                           MBInterface& mb,
                           MBTag tag,
                           MBEntityHandle p,
                           int depth,
                           int& idx )
{
  int id;
  MBErrorCode rval = mb.tag_get_data( tag, &p, 1, &id); CHECK_ERR(rval);
  CHECK_EQUAL( idx, id );
  
  MBRange verts;
  double coords[6][3];
  int num_vtx = coords_by_idx( idx, coords );
  rval = mb.get_entities_by_handle( p, verts );
  CHECK( verts.all_of_type( MBVERTEX ) );
  CHECK_EQUAL( num_vtx, (int)verts.size() );
  double coords2[6][3];
  rval = mb.get_coords( verts, &coords2[0][0] );
  std::vector<bool> match(6,true);
  for (int i = 0; i < num_vtx; ++i) {
    match[i] = false;
    for (int j = 0; j < num_vtx; ++j) {
      if (!match[j]) {
        double d[3] = { coords[i][0] - coords2[j][0],
                        coords[i][1] - coords2[j][1],
                        coords[i][2] - coords2[j][2] };
        double ds = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        if (ds < 1e-12) {
          match[j] = true;
          break;
        }
      }
    }
  }
  CHECK( match[0] );
  CHECK( match[1] );
  CHECK( match[2] );
  CHECK( match[3] );
  CHECK( match[4] );
  CHECK( match[5] );
 
  ++idx;
  
  std::vector<MBEntityHandle> children, parents;

  rval = mb.get_child_meshsets( p, children ); CHECK_ERR(rval);
  if (depth == max_depth) {
    CHECK_EQUAL( (size_t)0, children.size() );
    return;
  }
  
  CHECK_EQUAL( (size_t)2, children.size() );
  MBEntityHandle l = children.front();
  MBEntityHandle r = children.back();
  
  parents.clear();
  rval = mb.get_parent_meshsets( l, parents ); CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)1, parents.size() );
  CHECK_EQUAL( p, parents.front() );
  parents.clear();
  rval = mb.get_parent_meshsets( r, parents ); CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)1, parents.size() );
  CHECK_EQUAL( p, parents.front() );
  
  recursive_check_tree( max_depth, mb, tag, l, depth+1, idx );
  recursive_check_tree( max_depth, mb, tag, r, depth+1, idx );
}
 

void test_tree( int max_depth ) 
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  MBEntityHandle root;
  
  // create tag in which to store number for each tree node,
  // in depth-first in-order search order.
  MBTag tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", tag ); CHECK_ERR(rval);
  
  // create a binary tree to a depth of 20 (about 1 million nodes)
  rval = mb.create_meshset( MESHSET_SET, root ); CHECK_ERR(rval);
  int idx = 1;
  recursive_build_tree( max_depth, mb, tag, root, 1, idx );
  const int last_idx = idx;
  std::cerr << "Created binary tree containing " << last_idx << " nodes." << std::endl;
  
  std::ostringstream str;
  str << "tree-" << max_depth << ".h5m";
  
  // write file and read back in
  rval = mb.write_file( str.str().c_str(), 0, "BUFFER_SIZE=1024" ); CHECK_ERR(rval);
  mb.delete_mesh();
  rval = mb.load_file( str.str().c_str(), root );
  if (!keep_file)
    remove( str.str().c_str() );
  CHECK_ERR(rval);
  
  // get tree root
  rval = mb.tag_get_handle( "GLOBAL_ID", tag ); CHECK_ERR(rval);
  MBRange roots;
  idx = 1;
  const void* vals[] = {&idx};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, vals, 1, roots );
  CHECK_EQUAL( (MBEntityHandle)1, roots.size() );
  root = roots.front();
  
  // check that tree is as we expect it
  idx = 1;
  recursive_check_tree( max_depth, mb, tag, root, 1, idx );
  CHECK_EQUAL( last_idx, idx );
}

void test_small_tree() 
{ 
  int max_depth = 5;
  const char* str = getenv("MAX_DEPTH");
  if (str) {
    max_depth = atoi(str);
    CHECK(max_depth > 0);
  }
  test_tree( max_depth );
}

void test_big_tree()
  { test_tree( 20 ); }
  
int main(int argc, char* argv[])
{
  bool do_big_tree_test = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "-k")
      keep_file = true;
    else if (std::string(argv[i]) == "-b")
      do_big_tree_test = true;
    else {
      std::cerr << "Usage: " << argv[0] << " [-k] [-b]" << std::endl;
      return 1;
    }
  }

  // only one test so far... should probably add second test
  // for really-old-format  entityset parent/child links
  int exitval = 0;
  exitval += RUN_TEST( test_ranged_set_with_holes );
  exitval += RUN_TEST( test_file_set );
  exitval += RUN_TEST( test_small_tree );
  if (do_big_tree_test) {
    exitval += RUN_TEST( test_big_tree );
  }
  return exitval;
}

  
