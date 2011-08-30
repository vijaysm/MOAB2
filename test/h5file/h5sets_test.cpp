#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "TestUtil.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

using namespace moab;

#ifndef SRCDIR
#  define SRCDIR .
#endif

const char filename[] = "sets.h5m";
bool keep_file = false;

void read_write_file( Interface& output, Interface& input, EntityHandle* input_set = 0 )
{
  ErrorCode rval;
  rval = output.write_file( filename, 0, "DEBUG_BINIO" );
  CHECK_ERR(rval);
  if (input_set) {
    rval = input.create_meshset( MESHSET_SET, *input_set );
    CHECK_ERR(rval);
  }
  rval = input.load_file( filename, input_set );
  if (!keep_file)
    remove(filename);
  CHECK_ERR(rval);
}

void test_ranged_set_with_stale_handles()
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  Range verts;
  
  const int num_vtx = 40;
  std::vector<double> coords( 3*num_vtx, 0.0 );
  rval = mb.create_vertices( &coords[0], num_vtx, verts );
  CHECK_ERR(rval);
  CHECK_EQUAL(num_vtx, (int)verts.size());
  
  EntityHandle set;
  rval = mb.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  rval = mb.add_entities( set, verts );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> dead_verts;
  for (int i = num_vtx/4; i < num_vtx; i += num_vtx/4 ) {
    Range::iterator j = verts.begin();
    j += i;
    dead_verts.push_back( *j );
  }
  rval = mb.delete_entities( &dead_verts[0], dead_verts.size() );
  CHECK_ERR(rval);
  
  Core moab2;
  Interface& mb2 = moab2;
  EntityHandle file_set;
  read_write_file( mb, mb2, &file_set );
  Range sets;
  rval = mb2.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)sets.size() );
  EntityHandle other_set = sets.front() == file_set ? sets.back() : sets.front();
  
  int num_vtx2 = -5;
  rval = mb2.get_number_entities_by_type( other_set, MBVERTEX, num_vtx2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( (int)(num_vtx - dead_verts.size()), num_vtx2 );
}

void test_list_set_with_stale_handles()
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  Range verts;
  
  const int num_vtx = 40;
  std::vector<double> coords( 3*num_vtx, 0.0 );
  rval = mb.create_vertices( &coords[0], num_vtx, verts );
  CHECK_ERR(rval);
  CHECK_EQUAL(num_vtx, (int)verts.size());
  
  EntityHandle set;
  rval = mb.create_meshset( MESHSET_ORDERED, set );
  CHECK_ERR(rval);
  rval = mb.add_entities( set, verts );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> dead_verts;
  for (int i = num_vtx/4; i < num_vtx; i += num_vtx/4 ) {
    Range::iterator j = verts.begin();
    j += i;
    dead_verts.push_back( *j );
  }
  rval = mb.delete_entities( &dead_verts[0], dead_verts.size() );
  CHECK_ERR(rval);
  
  Core moab2;
  Interface& mb2 = moab2;
  EntityHandle file_set;
  read_write_file( mb, mb2, &file_set );
  Range sets;
  rval = mb2.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)sets.size() );
  EntityHandle other_set = sets.front() == file_set ? sets.back() : sets.front();
  
  std::vector<EntityHandle> list;
  rval = mb2.get_entities_by_handle( other_set, list );
  CHECK_ERR(rval);
  CHECK_EQUAL( verts.size() - dead_verts.size(), list.size() );
}

void test_file_set()
{
  ErrorCode rval;
  Core moab;
  double vtxcoords[] = { 0.0, 0.0, 0.0, 
                         1.0, 0.0, 0.0, 
                         0.0, 1.0, 0.0 };
  Range verts;
  rval = moab.create_vertices( vtxcoords, 3, verts );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, (int)verts.size() );
  
  EntityHandle tri;
  EntityHandle conn[3];
  std::copy( verts.begin(), verts.end(), conn );
  rval = moab.create_element( MBTRI, conn, 3, tri );
  CHECK_ERR(rval);
  
  EntityHandle set;
  rval = moab.create_meshset( MESHSET_ORDERED, set );
  CHECK_ERR(rval);
  rval = moab.add_entities( set, &tri, 1 );
  CHECK_ERR(rval);
  
  EntityHandle file;
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
                           Interface& mb,
                           Tag tag,
                           EntityHandle p,
                           int depth,
                           int& idx )
{
  ErrorCode rval = mb.tag_set_data( tag, &p, 1, &idx ); CHECK_ERR(rval);
  
  Range verts;
  double coords[6][3];
  int num_vtx = coords_by_idx( idx, coords );
  rval = mb.create_vertices( &coords[0][0], num_vtx, verts );
  rval = mb.add_entities( p, verts );
  ++idx;
  if (depth == max_depth)
    return;

  EntityHandle l, r;
  rval = mb.create_meshset( MESHSET_SET, l ); CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, r ); CHECK_ERR(rval);
  rval = mb.add_parent_child( p, l ); CHECK_ERR(rval);
  rval = mb.add_parent_child( p, r ); CHECK_ERR(rval);
  
  recursive_build_tree( max_depth, mb, tag, l, depth+1, idx );
  recursive_build_tree( max_depth, mb, tag, r, depth+1, idx );
}
 
void recursive_check_tree( int max_depth,
                           Interface& mb,
                           Tag tag,
                           EntityHandle p,
                           int depth,
                           int& idx )
{
  int id;
  ErrorCode rval = mb.tag_get_data( tag, &p, 1, &id); CHECK_ERR(rval);
  CHECK_EQUAL( idx, id );
  
  Range verts;
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
  
  std::vector<EntityHandle> children, parents;

  rval = mb.get_child_meshsets( p, children ); CHECK_ERR(rval);
  if (depth == max_depth) {
    CHECK_EQUAL( (size_t)0, children.size() );
    return;
  }
  
  CHECK_EQUAL( (size_t)2, children.size() );
  EntityHandle l = children.front();
  EntityHandle r = children.back();
  
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
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  EntityHandle root;
  
  // create tag in which to store number for each tree node,
  // in depth-first in-order search order.
  Tag tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, tag ); CHECK_ERR(rval);
  
  // create a binary tree to a depth of 20 (about 1 million nodes)
  rval = mb.create_meshset( MESHSET_SET, root ); CHECK_ERR(rval);
  int idx = 1;
  recursive_build_tree( max_depth, mb, tag, root, 1, idx );
  const int last_idx = idx;
  std::cerr << "Created binary tree containing " << last_idx << " nodes." << std::endl;
  
  std::ostringstream str;
  str << "tree-" << max_depth << ".h5m";
  
  // write file and read back in
  rval = mb.write_file( str.str().c_str(), 0, "BUFFER_SIZE=1024;DEBUG_BINIO" ); CHECK_ERR(rval);
  mb.delete_mesh();
  rval = mb.load_file( str.str().c_str() );
  if (!keep_file)
    remove( str.str().c_str() );
  CHECK_ERR(rval);
  
  // get tree root
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, tag ); CHECK_ERR(rval);
  Range roots;
  idx = 1;
  const void* vals[] = {&idx};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, vals, 1, roots );
  CHECK_EQUAL( (size_t)1, roots.size() );
  root = roots.front();
  
  // check that tree is as we expect it
  idx = 1;
  recursive_check_tree( max_depth, mb, tag, root, 1, idx );
  CHECK_EQUAL( last_idx, idx );
}

void test_small_tree() 
{ 
  int max_depth = 8;
  const char* str = getenv("MAX_DEPTH");
  if (str) {
    max_depth = atoi(str);
    CHECK(max_depth > 0);
  }
  test_tree( max_depth );
}

void test_big_tree()
  { test_tree( 20 ); }

 
void regression_mmiller_8_2010();

void test_set_flags();
  
int main(int argc, char* argv[])
{
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif

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

  int exitval = 0;
  exitval += RUN_TEST( test_ranged_set_with_stale_handles );
  exitval += RUN_TEST( test_list_set_with_stale_handles );
  exitval += RUN_TEST( test_file_set );
  exitval += RUN_TEST( test_small_tree );
  exitval += RUN_TEST( test_set_flags );
  exitval += RUN_TEST( regression_mmiller_8_2010 );
  if (do_big_tree_test) {
    exitval += RUN_TEST( test_big_tree );
  }

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return exitval;
}

// NOTE: this test makes some assuptions about handles:
//       mainly that they will be assigned sequentially 
//       in the same order as defined in the file and
//       beginning with ID 1
void regression_mmiller_8_2010()
{
  Core moab;
  Interface& mb = moab;
  
  const size_t num_vtx = 171;
  const size_t num_pri = 12;
  const size_t num_pyr = 8;
  const size_t num_hex = 100;
  const size_t num_set = 25;

  mb.load_file( STRINGIFY(SRCDIR) "/rocket_ents_in_assm.h5m" );

/* Dump of set contents from input file:
 1r: 172, 4, 
 2r: 192, 4, 204, 4, 216, 4, 228, 4, 240, 4, 252, 4, 264, 4, 276, 4, 
 3r: 288, 4, 
 4 : 181, 183, 185, 187, 
 5r: 176, 5, 182, 1, 184, 1, 186, 1, 188, 4, 196, 8, 208, 8, 220, 8, 232, 8, 244, 8, 256, 8, 268, 8, 280, 8, 
 6r: 172, 4, 192, 4, 204, 4, 216, 4, 228, 4, 240, 4, 252, 4, 264, 4, 276, 4, 288, 4, 
 7r: 176, 4, 188, 4, 196, 8, 208, 8, 220, 8, 232, 8, 244, 8,
 8r: 180, 8, 256, 8, 268, 8, 280, 8,
 9r: 172, 120, 301, 1, 309, 1, 
10r: 176, 4, 188, 100, 302, 1, 308, 1, 
11r: 176, 4, 188, 52, 303, 1, 
12r: 176, 4, 188, 4, 304, 4, 
13 : 177, 189, 
14 : 178, 190, 
15 : 179, 191, 
16 : 
17r: 240, 48,
18r: 172, 4, 180, 8, 288, 4, 310, 1, 312, 1, 
19r: 180, 8, 288, 4, 311, 1, 
20r: 180, 8, 
21r: 172, 4, 313, 4, 
22 : 173, 
23 : 174, 
24 : 175, 
25 : 176, 188
*/


    // check expected handles
    
  const EntityHandle VTX1 = CREATE_HANDLE(MBVERTEX,1);
  Range range, expected;
  mb.get_entities_by_type( 0, MBVERTEX, range );
  CHECK_EQUAL( num_vtx, range.size() );
  expected.insert( VTX1, VTX1+num_vtx-1 );
  CHECK_EQUAL( expected, range );
    
  const EntityHandle PRI1 = CREATE_HANDLE(MBPRISM,1);
  range.clear();
  expected.clear();
  mb.get_entities_by_type( 0, MBPRISM, range );
  CHECK_EQUAL( num_pri, range.size() );
  expected.insert( PRI1, PRI1+num_pri-1 );
  CHECK_EQUAL( expected, range );
    
  const EntityHandle PYR1 = CREATE_HANDLE(MBPYRAMID,1);
  range.clear();
  expected.clear();
  mb.get_entities_by_type( 0, MBPYRAMID, range );
  CHECK_EQUAL( num_pyr, range.size() );
  expected.insert( PYR1, PYR1+num_pyr-1 );
  CHECK_EQUAL( expected, range );
    
  const EntityHandle HEX1 = CREATE_HANDLE(MBHEX,1);
  range.clear();
  expected.clear();
  mb.get_entities_by_type( 0, MBHEX, range );
  CHECK_EQUAL( num_hex, range.size() );
  expected.insert( HEX1, HEX1+num_hex-1 );
  CHECK_EQUAL( expected, range );
  
  const EntityHandle SET1 = CREATE_HANDLE(MBENTITYSET,1);
  range.clear();
  expected.clear();
  mb.get_entities_by_type( 0, MBENTITYSET, range );
  CHECK_EQUAL( num_set, range.size() );
  expected.insert( SET1, SET1+num_set-1 );
  CHECK_EQUAL( expected, range );
  
  // Check set contents
  
  // Set 1: Pyramids 1 to 4
  range.clear();
  mb.get_entities_by_handle( SET1, range );
  expected.clear();
  expected.insert( PYR1+0, PYR1+3 );
  CHECK_EQUAL( expected, range );

  // Skip sets 2 through 8 because they're long and complicated and
  // I doubt I could code up the content lists explicitly from the
  // dump of the HDF5 file w/out many mistakes
  
  // Set 9: Pyramids 1 to 8, Prism 1 to 12, Hex 1 to 100, and Sets 10 and 18
  range.clear();
  mb.get_entities_by_handle( SET1+8, range );
  expected.clear();
  expected.insert( PYR1+0, PYR1+7 );
  expected.insert( PRI1+0, PRI1+11 );
  expected.insert( HEX1+0, HEX1+99 );
  expected.insert( SET1+9 );
  expected.insert( SET1+17 );
  CHECK_EQUAL( expected, range );
  
  // Set 10: Pyramids 5 to 8, Prism 9 to 12, Hex 1 to 96, and Sets 11 and 17
  range.clear();
  mb.get_entities_by_handle( SET1+9, range );
  expected.clear();
  expected.insert( PYR1+4, PYR1+7 );
  expected.insert( PRI1+8, PRI1+11 );
  expected.insert( HEX1+0, HEX1+95 );
  expected.insert( SET1+10 );
  expected.insert( SET1+16 );
  CHECK_EQUAL( expected, range );
  
  // Set 11: Pyramids 5 to 8, Prism 9 to 12, Hex 1 to 48, and Set 12
  range.clear();
  mb.get_entities_by_handle( SET1+10, range );
  expected.clear();
  expected.insert( PYR1+4, PYR1+7 );
  expected.insert( PRI1+8, PRI1+11 );
  expected.insert( HEX1+0, HEX1+47 );
  expected.insert( SET1+11 );
  CHECK_EQUAL( expected, range );
  
  // Set 12: Pyramids 5 to 8, Prism 9 to 12, and Sets 13 to 16
  range.clear();
  mb.get_entities_by_handle( SET1+11, range );
  expected.clear();
  expected.insert( PYR1+4, PYR1+7 );
  expected.insert( PRI1+8, PRI1+11 );
  expected.insert( SET1+12, SET1+15 );
  CHECK_EQUAL( expected, range );
  
  // Set 13: Pyramids 6 and Prism 10
  range.clear();
  mb.get_entities_by_handle( SET1+12, range );
  expected.clear();
  expected.insert( PYR1+5 );
  expected.insert( PRI1+9 );
  CHECK_EQUAL( expected, range );
  
  // Set 14: Pyramids 7 and Prism 11
  range.clear();
  mb.get_entities_by_handle( SET1+13, range );
  expected.clear();
  expected.insert( PYR1+6 );
  expected.insert( PRI1+10 );
  CHECK_EQUAL( expected, range );
  
  // Set 15: Pyramids 8 and Prism 12
  range.clear();
  mb.get_entities_by_handle( SET1+14, range );
  expected.clear();
  expected.insert( PYR1+7 );
  expected.insert( PRI1+11 );
  CHECK_EQUAL( expected, range );
  
  // Set 16: Empty
  range.clear();
  mb.get_entities_by_handle( SET1+15, range );
  expected.clear();
  CHECK_EQUAL( expected, range );
  
  // Set 17: Hex 49 to 96
  range.clear();
  mb.get_entities_by_handle( SET1+16, range );
  expected.clear();
  expected.insert( HEX1+48, HEX1+95 );
  CHECK_EQUAL( expected, range );
   
  // Set 18: Pyramids 1 to 4, Prism 1 to 8, Hex 97 to 100, and Sets 19 and 21
  range.clear();
  mb.get_entities_by_handle( SET1+17, range );
  expected.clear();
  expected.insert( PYR1+0, PYR1+3 );
  expected.insert( PRI1+0, PRI1+7 );
  expected.insert( HEX1+96, HEX1+99 );
  expected.insert( SET1+18 );
  expected.insert( SET1+20 );
  CHECK_EQUAL( expected, range );
   
  // Set 19: Prism 1 to 8, Hex 97 to 100, and Set 20  
  range.clear();
  mb.get_entities_by_handle( SET1+18, range );
  expected.clear();
  expected.insert( PRI1+0, PRI1+7 );
  expected.insert( HEX1+96, HEX1+99 );
  expected.insert( SET1+19 );
  CHECK_EQUAL( expected, range );
   
  // Set 20: Prism 1 to 8
  range.clear();
  mb.get_entities_by_handle( SET1+19, range );
  expected.clear();
  expected.insert( PRI1+0, PRI1+7 );
  CHECK_EQUAL( expected, range );
   
  // Set 21: Pyramids 1 to 4, and Sets 22 to 25
  range.clear();
  mb.get_entities_by_handle( SET1+20, range );
  expected.clear();
  expected.insert( PYR1+0, PYR1+3 );
  expected.insert( SET1+21, SET1+24 );
  CHECK_EQUAL( expected, range );
   
  // Set 22: Pyramid 2
  range.clear();
  mb.get_entities_by_handle( SET1+21, range );
  expected.clear();
  expected.insert( PYR1+1 );
  CHECK_EQUAL( expected, range );
   
  // Set 23: Pyramid 3
  range.clear();
  mb.get_entities_by_handle( SET1+22, range );
  expected.clear();
  expected.insert( PYR1+2 );
  CHECK_EQUAL( expected, range );
   
  // Set 24: Pyramid 4
  range.clear();
  mb.get_entities_by_handle( SET1+23, range );
  expected.clear();
  expected.insert( PYR1+3 );
  CHECK_EQUAL( expected, range );
   
  // Set 25: Pyramid 5 and Prism 9
  range.clear();
  mb.get_entities_by_handle( SET1+24, range );
  expected.clear();
  expected.insert( PYR1+4 );
  expected.insert( PRI1+8 );
  CHECK_EQUAL( expected, range );
}

// Test to reproduce bug reported by brandom smith on 2011-3-7
// and test other possible issues with the somewhat inconsistant
// meshset creation flags.  Bug was fixed in SVN revision 4548.
void test_set_flags()
{
  const char filename[] = "test_set_flags.h5m";
  ErrorCode rval;
  Core core;
  Interface& mb = core;

  // create a bunch of vertices so we have something to put in sets
  const int nverts = 20;
  double coords[3*nverts] = {0.0};
  Range verts;
  rval = mb.create_vertices( coords, nverts, verts );
  CHECK_ERR(rval);
  
  // Assign IDs to things so that we can identify them in the
  // data we read back in.
  Tag tag;
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, tag ); CHECK_ERR(rval);
  int ids[nverts];
  for (int i = 0; i < nverts; ++i)
    ids[i] = i+1;
  rval = mb.tag_set_data( tag, verts, ids ); CHECK_ERR(rval);
  
  // define two lists of vertex ids corresponding to the
  // vertices that we are going to put into different sets
  const int set_verts1[] = { 1, 2, 3, 4, 8, 13, 14, 15 };
  const int set_verts2[] = { 3, 9, 10, 11, 12, 13, 14, 15, 16, 17 };
  const int num_verts1 = sizeof(set_verts1)/sizeof(set_verts1[0]);
  const int num_verts2 = sizeof(set_verts1)/sizeof(set_verts1[0]);
  
  // convert to handle lists
  EntityHandle set_handles1[num_verts1], set_handles2[num_verts2];
  for (int i = 0; i < num_verts1; ++i)
    set_handles1[i] = *(verts.begin() + set_verts1[i] - 1);
  for (int i = 0; i < num_verts2; ++i)
    set_handles2[i] = *(verts.begin() + set_verts2[i] - 1);
  
  // now create some sets with different flag combinations
  EntityHandle sets[6];
  rval = mb.create_meshset( 0, sets[0] );
  rval = mb.create_meshset( MESHSET_TRACK_OWNER, sets[1] );
  rval = mb.create_meshset( MESHSET_SET, sets[2] );
  rval = mb.create_meshset( MESHSET_SET|MESHSET_TRACK_OWNER, sets[3] );
  rval = mb.create_meshset( MESHSET_ORDERED, sets[4] );
  rval = mb.create_meshset( MESHSET_ORDERED|MESHSET_TRACK_OWNER, sets[5] );
  
  // assign IDs to sets so that we can identify them later
  rval = mb.tag_set_data( tag, sets, 6, ids ); CHECK_ERR(rval);
  // add entities to sets
  rval = mb.add_entities( sets[0], set_handles1, num_verts1 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], set_handles2, num_verts2 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[2], set_handles1, num_verts1 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[3], set_handles2, num_verts2 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[4], set_handles1, num_verts1 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[5], set_handles2, num_verts2 ); CHECK_ERR(rval);
  
  // now write the file and read it back in
  rval = mb.write_file( filename, 0, "BUFFER_SIZE=1024;DEBUG_BINIO" ); CHECK_ERR(rval);
  mb.delete_mesh();
  rval = mb.load_file( filename );
  if (!keep_file)
    remove( filename );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( "GLOBAL_ID", 1, MB_TYPE_INTEGER, tag ); CHECK_ERR(rval);
  
  // find our sets
  Range tmp;
  for (int i = 0; i < 6; ++i) {
    int id = i+1;
    tmp.clear();
    const void* vals[] = {&id};
    rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, vals, 1, tmp ); CHECK_ERR(rval);
    CHECK_EQUAL( 1u, (unsigned)tmp.size() );
    sets[i] = tmp.front();
  }
  
  // check that sets have correct flags
  unsigned opts;
  rval = mb.get_meshset_options( sets[0], opts ); CHECK_ERR(rval);
  CHECK_EQUAL( 0u, opts );
  rval = mb.get_meshset_options( sets[1], opts ); CHECK_ERR(rval);
  CHECK_EQUAL( (unsigned)MESHSET_TRACK_OWNER, opts );
  rval = mb.get_meshset_options( sets[2], opts ); CHECK_ERR(rval);
  CHECK_EQUAL( (unsigned)MESHSET_SET, opts );
  rval = mb.get_meshset_options( sets[3], opts ); CHECK_ERR(rval);
  CHECK_EQUAL( (unsigned)(MESHSET_SET|MESHSET_TRACK_OWNER), opts );
  rval = mb.get_meshset_options( sets[4], opts ); CHECK_ERR(rval);
  CHECK_EQUAL( (unsigned)MESHSET_ORDERED, opts );
  rval = mb.get_meshset_options( sets[5], opts ); CHECK_ERR(rval);
  CHECK_EQUAL( (unsigned)(MESHSET_ORDERED|MESHSET_TRACK_OWNER), opts );
  
  // check that sets have correct contents
  int set_ids1[num_verts1], set_ids2[num_verts2];
  
  tmp.clear();
  rval = mb.get_entities_by_handle( sets[0], tmp ); CHECK_ERR(rval);
  CHECK_EQUAL( num_verts1, (int)tmp.size() );
  rval = mb.tag_get_data( tag, tmp, set_ids1 ); CHECK_ERR(rval);
  std::sort( set_ids1, set_ids1+num_verts1 );
  CHECK_ARRAYS_EQUAL( set_verts1, num_verts1, set_ids1, num_verts1 );
  
  tmp.clear();
  rval = mb.get_entities_by_handle( sets[1], tmp ); CHECK_ERR(rval);
  CHECK_EQUAL( num_verts2, (int)tmp.size() );
  rval = mb.tag_get_data( tag, tmp, set_ids2 ); CHECK_ERR(rval);
  std::sort( set_ids2, set_ids2+num_verts2 );
  CHECK_ARRAYS_EQUAL( set_verts2, num_verts2, set_ids2, num_verts2 );
  
  tmp.clear();
  rval = mb.get_entities_by_handle( sets[2], tmp ); CHECK_ERR(rval);
  CHECK_EQUAL( num_verts1, (int)tmp.size() );
  rval = mb.tag_get_data( tag, tmp, set_ids1 ); CHECK_ERR(rval);
  std::sort( set_ids1, set_ids1+num_verts1 );
  CHECK_ARRAYS_EQUAL( set_verts1, num_verts1, set_ids1, num_verts1 );
  
  tmp.clear();
  rval = mb.get_entities_by_handle( sets[3], tmp ); CHECK_ERR(rval);
  CHECK_EQUAL( num_verts2, (int)tmp.size() );
  rval = mb.tag_get_data( tag, tmp, set_ids2 ); CHECK_ERR(rval);
  std::sort( set_ids2, set_ids2+num_verts2 );
  CHECK_ARRAYS_EQUAL( set_verts2, num_verts2, set_ids2, num_verts2 );
  
  tmp.clear();
  rval = mb.get_entities_by_handle( sets[4], tmp ); CHECK_ERR(rval);
  CHECK_EQUAL( num_verts1, (int)tmp.size() );
  rval = mb.tag_get_data( tag, tmp, set_ids1 ); CHECK_ERR(rval);
  std::sort( set_ids1, set_ids1+num_verts1 );
  CHECK_ARRAYS_EQUAL( set_verts1, num_verts1, set_ids1, num_verts1 );
  
  tmp.clear();
  rval = mb.get_entities_by_handle( sets[5], tmp ); CHECK_ERR(rval);
  CHECK_EQUAL( num_verts2, (int)tmp.size() );
  rval = mb.tag_get_data( tag, tmp, set_ids2 ); CHECK_ERR(rval);
  std::sort( set_ids2, set_ids2+num_verts2 );
  CHECK_ARRAYS_EQUAL( set_verts2, num_verts2, set_ids2, num_verts2 );
}
