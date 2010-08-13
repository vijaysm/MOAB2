#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "TestUtil.hpp"

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

void test_ranged_set_with_holes()
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
  mb2.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_EQUAL( 2, (int)sets.size() );
  EntityHandle other_set = sets.front() == file_set ? sets.back() : sets.front();
  
  int num_vtx2 = -5;
  rval = mb2.get_number_entities_by_type( other_set, MBVERTEX, num_vtx2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( (int)(num_vtx - dead_verts.size()), num_vtx2 );
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
  rval = mb.write_file( str.str().c_str(), 0, "BUFFER_SIZE=1024;DEBUG_BINIO" ); CHECK_ERR(rval);
  mb.delete_mesh();
  rval = mb.load_file( str.str().c_str() );
  if (!keep_file)
    remove( str.str().c_str() );
  CHECK_ERR(rval);
  
  // get tree root
  rval = mb.tag_get_handle( "GLOBAL_ID", tag ); CHECK_ERR(rval);
  Range roots;
  idx = 1;
  const void* vals[] = {&idx};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, vals, 1, roots );
  CHECK_EQUAL( (EntityHandle)1, roots.size() );
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

 
void regression_mmiller_8_2010();
  
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
  //exitval += RUN_TEST( test_ranged_set_with_holes );
  //exitval += RUN_TEST( test_file_set );
  //exitval += RUN_TEST( test_small_tree );
  exitval += RUN_TEST( regression_mmiller_8_2010 );
  if (do_big_tree_test) {
    exitval += RUN_TEST( test_big_tree );
  }
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
  
  const EntityHandle num_vtx = 171;
  const EntityHandle num_pri = 12;
  const EntityHandle num_pyr = 8;
  const EntityHandle num_hex = 100;
  const EntityHandle num_set = 25;

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
