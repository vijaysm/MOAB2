#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "TestRunner.hpp"
#include "ReadHDF5.hpp"
#include "MBTagConventions.hpp"
#include "FileOptions.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include <vector>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace moab;

const char TEST_FILE[] = "partial.h5m";
#define READ_OPTS "BUFFER_SIZE=256"
const char ID_TAG_NAME[] = "test_id_tag";


static void test_read_nothing_common( bool non_existant );
static void test_read_nodes_common( int num_read_sets, bool blocked_coordinate_io );
static void test_read_handle_tag_common( bool var_len );

const int MBQUAD_INT = 20; 
const int NUM_SETS = 10;
const int SET_WIDTH = (MBQUAD_INT + NUM_SETS - 1) / NUM_SETS; // ceil(MBQUAD_INT/NUM_SETS)
const char LOGICAL_NAME[] = "logical";   // tag storing logical (i,j) coordinates
const char CENTROID_NAME[] = "centroid"; // tag storing position of centroid (x,y,0)
//! Create a regular MBQUAD_INT^2 element quad mesh with regularly
//! spaced coordinates in the range [1,100].  Group elements
//! into 10 vertical strips MBQUAD_INT/10 elements wide.  Tag elements,
//! vertices and/or sets with ID in [1,10] stored in ID_TAG_NAME
//! tag.  Write new mesh to TEST_FILE.
void create_mesh( bool create_element_sets,
                  bool create_vertex_sets,
                  bool tag_elements_with_id,
                  bool tag_vertices_with_id,
                  const char* adj_elem_tag_name = 0,
                  bool var_len_adj_elems = false );
// Given a list of vertices adjacent to a quad strip, identify it as one of the 
// NUM_SETS strips of quads written by create_mesh.
int identify_set( Interface& mb, const Range& verts );
int identify_set( Interface& mb, EntityHandle set );

static Tag check_tag( Interface& mb, 
                        const char* name,
                        TagType storage,
                        DataType type,
                        int size );

enum GatherTestMode { GATHER_SETS, GATHER_CONTENTS, GATHER_NONE };
void test_gather_sets_common( bool contained_sets, GatherTestMode mode, bool no_parent_containing_sets = false );
void test_gather_sets_ranged( bool contained_sets, GatherTestMode mode, bool no_parent_containing_sets = false );


//! Read a set containing no entities
void test_read_empty_set()
  { test_read_nothing_common( false ); }
  
//! Specify ID that doesn't exist in file
void test_read_non_existant_set()
  { test_read_nothing_common( true ); }

//! Read in the nodes contained in a set.
void test_read_one_set_nodes()
  { test_read_nodes_common(1, false); }

//! Read in the nodes contained in a set.
void test_read_one_set_nodes_blocked()
  { test_read_nodes_common(1, true); }

//! Read in the elems contained in a set
void test_read_one_set_elems();

//! Read in the polyhedra contained in a set
void test_read_one_set_polyhedra();

//! Read in the sets contained in a set.  
//! Should read all sets containing read elements or nodes
//! and all sets that are contained the the specified "read"
//! set.  Test the later here.
void test_read_set_sets();

//! Read in the nodes contained in a sets.
void test_read_two_sets_nodes()
  { test_read_nodes_common(2,false); }

//! Read in the elems contained in a sets
void test_read_two_sets_elems();

//! For any set selected to be read by either explicit designation,
//! containing read entities, or contained in an explcitly designated
//! set, any child sets are also read.  Check that here.
void test_read_child_sets_only()
{ test_gather_sets_common( false, GATHER_SETS);
  test_gather_sets_ranged( false, GATHER_SETS); }
void test_read_child_set_contents()
{ test_gather_sets_common( false, GATHER_CONTENTS); 
  test_gather_sets_ranged( false, GATHER_CONTENTS); }
void test_read_no_child_sets()
{ test_gather_sets_common( false, GATHER_NONE);
  test_gather_sets_ranged( false, GATHER_NONE); }

//! For any set selected to be read by either explicit designation,
//! containing read entities, or contained in an explcitly designated
//! set, any contained sets are also read.  Check that here.
void test_read_contained_sets_only()
{ test_gather_sets_common( true, GATHER_SETS, true);
  test_gather_sets_ranged( true, GATHER_SETS); }
void test_read_contained_set_contents()
{ test_gather_sets_common( true, GATHER_CONTENTS, true );
  test_gather_sets_ranged( true, GATHER_CONTENTS); }
void test_read_no_contained_sets()
{ test_gather_sets_common( true, GATHER_NONE, true );
  test_gather_sets_ranged( true, GATHER_NONE); }

//! Read in the sets contained in a set.  
//! Should read all sets containing read elements or nodes
//! and all sets that are contained the the specified "read"
//! set.  Test the former here.
void test_read_containing_sets();

//! Test reading of explicit adjacencies
void test_read_adjacencies();

//! Test reading of sparse double tag data
void test_read_double_tag();

//! Test reading of sparse opaque tag data
void test_read_opaque_tag();

//! Test reading of sparse handle tag data
void test_read_handle_tag()
  { test_read_handle_tag_common(false); }

//! Test reading of variable-length tag data
void test_var_len_tag()
  { test_read_handle_tag_common(true); }

void test_read_tagged_elems();

void test_read_tagged_nodes();

void test_read_sides();

void test_read_ids();

void test_read_partial_ids();

int main( int argc, char* argv[] )
{
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif

  REGISTER_TEST(test_read_empty_set);
  REGISTER_TEST(test_read_non_existant_set);
  REGISTER_TEST(test_read_one_set_nodes);
  REGISTER_TEST(test_read_one_set_nodes_blocked);
  REGISTER_TEST(test_read_one_set_elems);
  REGISTER_TEST(test_read_one_set_polyhedra);
  REGISTER_TEST(test_read_set_sets);
  REGISTER_TEST(test_read_two_sets_nodes);
  REGISTER_TEST(test_read_two_sets_elems);
  REGISTER_TEST(test_read_child_sets_only);
  REGISTER_TEST(test_read_child_set_contents);
  REGISTER_TEST(test_read_no_child_sets);
  REGISTER_TEST(test_read_contained_sets_only);
  REGISTER_TEST(test_read_contained_set_contents);
  REGISTER_TEST(test_read_no_contained_sets);
  REGISTER_TEST(test_read_containing_sets);
  REGISTER_TEST(test_read_double_tag);
  REGISTER_TEST(test_read_opaque_tag);
  REGISTER_TEST(test_read_handle_tag);
  REGISTER_TEST(test_var_len_tag);
  REGISTER_TEST(test_read_adjacencies);
  REGISTER_TEST(test_read_tagged_elems);
  REGISTER_TEST(test_read_tagged_nodes);
  REGISTER_TEST(test_read_sides);
  REGISTER_TEST(test_read_ids);
  REGISTER_TEST(test_read_partial_ids);
  int result = RUN_TESTS( argc, argv );

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return result;
}

void test_read_nothing_common( bool non_existant )
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;

  // create a few nodes to write to file
  std::vector<double> coords( 3000 );
  Range verts;
  rval= mb.create_vertices( &coords[0], coords.size()/3, verts );
  CHECK_ERR(rval);
  
  // create three entity sets
  EntityHandle sets[3];
  rval = mb.create_meshset( MESHSET_SET, sets[0] );
  CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, sets[1] );
  CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_ORDERED, sets[2] );
  CHECK_ERR(rval);
  
  // put all vertices into two of the sets
  rval = mb.add_entities( sets[0], verts );
  CHECK_ERR(rval);
  rval = mb.add_entities( sets[2], verts );
  CHECK_ERR(rval);
  
  // tag all three sets
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  int ids[3] = { 5, 7, 9 };
  rval = mb.tag_set_data( id_tag, sets, 3, ids );
  CHECK_ERR(rval);
  
  // write mesh
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  rval = mb.delete_mesh();
  CHECK_ERR(rval);
  
  // now read back in only the empty set
  EntityHandle file_set;
  int id = non_existant ? 8 : 7;
  rval = mb.create_meshset( MESHSET_SET, file_set );
  CHECK_ERR(rval);
  rval = mb.load_file( TEST_FILE, &file_set, READ_OPTS, ID_TAG_NAME, &id, 1 );
  if (non_existant) {
    CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
    return;
  }
  else
    CHECK_ERR( rval );
  
  // the file should contain exactly two sets (the specified one and the new
  // file set, and nothing else.)
  for (EntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
    int count = -1;
    rval = mb.get_number_entities_by_type( 0, t, count );
    CHECK_ERR(rval);
    CHECK_EQUAL( 0, count );
  }
  Range setrange;
  rval = mb.get_entities_by_type( 0, MBENTITYSET, setrange );
  CHECK_ERR(rval);
  CHECK_EQUAL( (non_existant ? 1 : 2), (int)setrange.size() );
  CHECK( setrange.find(file_set) != setrange.end() );
}


static void vtx_coords( int set_id, int j, int num_sets, double coords[3] )
{
  int i = num_sets*j + set_id;
  coords[0] = i;
  coords[1] = i+0.25;
  coords[2] = i+0.5;
}

void test_read_nodes_common( int num_read_sets, bool blocked )
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
    
    // create 1000 nodes
  const int num_sets = 2*num_read_sets;
  std::vector<EntityHandle> verts(1000);
  std::vector< std::vector<EntityHandle> > set_verts(num_sets);
  for (size_t i = 0; i < verts.size(); ++i) {
    double coords[3];
    int j = i % num_sets;
    vtx_coords( j+1, set_verts[j].size(), num_sets, coords );
    rval = mb.create_vertex( coords, verts[i] );
    set_verts[ j ].push_back( verts[i] );
    CHECK_ERR(rval);
  }
  
    // create two sets, each containing half of the nodes
  std::vector<EntityHandle> sets(num_sets);
  for (int i = 0; i < num_sets; ++i) {
    rval = mb.create_meshset( MESHSET_ORDERED, sets[i] );
    CHECK_ERR(rval);
    rval = mb.add_entities( sets[i], &set_verts[i][0], set_verts[i].size() );
    CHECK_ERR(rval);
  }
 
    // tag both sets
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  std::vector<int> values( num_sets );
  for (int i = 0; i < num_sets; ++i) values[i] = i+1;
  rval = mb.tag_set_data( id_tag, &sets[0], num_sets, &values[0] );
  CHECK_ERR(rval);
  
    // write file
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  rval = mb.delete_mesh();
  CHECK_ERR(rval);
  
    // now read back in only the specified number of sets
  std::string opts(READ_OPTS);
  if (!opts.empty())
    opts += ';';
  if (blocked)
    opts += "BLOCKED_COORDINATE_IO=yes";
  else
    opts += "BLOCKED_COORDINATE_IO=no";
    
  values.resize( num_read_sets );
  for (int i = 0; i < num_read_sets; ++i) values[i] = 2*(i+1);
  EntityHandle file_set;
  rval = mb.create_meshset( MESHSET_SET, file_set );
  CHECK_ERR(rval);
  rval = mb.load_file( TEST_FILE, &file_set, opts.c_str(), ID_TAG_NAME, &values[0], num_read_sets );
  CHECK_ERR(rval);
  
  int count, expected = 0;
  rval = mb.get_number_entities_by_dimension( 0, 0, count );
  CHECK_ERR(rval);
  for (int i = 0; i < num_sets; ++i)
    if (i % 2)
      expected += set_verts[i].size();
  CHECK_EQUAL( expected, count );
  
  Range sets2;
  rval = mb.get_entities_by_type( 0, MBENTITYSET, sets2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1+num_read_sets, (int)sets2.size() );
  Range::iterator it = sets2.find( file_set );
  CHECK( it != sets2.end() );
  sets2.erase( it );
  
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
  CHECK_ERR(rval);  
  while (!sets2.empty()) {
    EntityHandle set = sets2.pop_front();
    int id;
    rval = mb.tag_get_data( id_tag, &set, 1, &id );
    CHECK_ERR(rval);
    CHECK( std::find(values.begin(), values.end(), id) != values.end() );
    CHECK( id > 0 );
    CHECK( (unsigned)id <= set_verts.size() );
    
    std::vector<EntityHandle> verts;
    rval = mb.get_entities_by_handle( set, verts );
    CHECK_ERR(rval);
    CHECK_EQUAL( set_verts[id-1].size(), verts.size() );
    
    for (size_t i = 0; i < verts.size(); ++i) {
      double exp_coords[3], coords[3];
      vtx_coords( id, i, num_sets, exp_coords );
      rval = mb.get_coords( &verts[i], 1, coords );
      CHECK_ERR(rval);
      CHECK_REAL_EQUAL( exp_coords[0], coords[0], 1e-12 );
      CHECK_REAL_EQUAL( exp_coords[1], coords[1], 1e-12 );
      CHECK_REAL_EQUAL( exp_coords[2], coords[2], 1e-12 );
    }
  }
}

//! Create a regular MBQUAD_INT^2 element quad mesh with regularly
//! spaced coordinates in the range [1,100].  Group elements
//! into 10 vertical strips MBQUAD_INT/10 elements wide.  Tag elements,
//! vertices and/or sets with ID in [1,10] stored in ID_TAG_NAME
//! tag.  Write new mesh to TEST_FILE.
void create_mesh( bool create_element_sets,
                  bool create_vertex_sets,
                  bool tag_elements_with_id,
                  bool tag_vertices_with_id,
                  const char* adj_elem_tag_name,
                  bool var_len_adj_elems )
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create tags
  Tag logical_tag, centroid_tag, id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( LOGICAL_NAME, 2*sizeof(int), MB_TYPE_OPAQUE, logical_tag, MB_TAG_DENSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( CENTROID_NAME, 3, MB_TYPE_DOUBLE, centroid_tag, MB_TAG_DENSE|MB_TAG_EXCL );
  CHECK_ERR(rval);

  EntityHandle sets[NUM_SETS];
  if (create_element_sets || create_vertex_sets) {
    for (int i = 0; i < NUM_SETS; ++i) {
      rval = mb.create_meshset( MESHSET_ORDERED, sets[i] );
      CHECK_ERR(rval);
      int id = i+1;
      rval = mb.tag_set_data( id_tag, &sets[i], 1, &id );
      CHECK_ERR(rval);
    }
  }

    // create elements
  EntityHandle verts[MBQUAD_INT+1][MBQUAD_INT+1], quads[MBQUAD_INT][MBQUAD_INT];
  for (int i = 0; i <= MBQUAD_INT; ++i) for(int j = 0; j <= MBQUAD_INT; ++j) {
    double coords[3] = { i, j, 0 };
    rval = mb.create_vertex( coords, verts[j][i] );
    CHECK_ERR(rval);
    int logical[2] = { i, j };
    rval = mb.tag_set_data( logical_tag, &verts[j][i], 1, logical );
    CHECK_ERR(rval);
    rval = mb.tag_set_data( centroid_tag, &verts[j][i], 1, coords );
    CHECK_ERR(rval);
    int id = (i-1)/SET_WIDTH + 1; // Note: assumes SET_WIDTH > 1
    if (tag_vertices_with_id) {
      rval = mb.tag_set_data( id_tag, &verts[j][i], 1, &id );
      CHECK_ERR(rval);
    }
    if (create_vertex_sets) {
      rval = mb.add_entities( sets[id-1], &verts[j][i], 1 );
      CHECK_ERR(rval);
        // Some vertices are shared by quads in different sets.
        // put such vertices in both sets.
      int id2 = i/SET_WIDTH + 1;
      if (id2 != id && id2 <= NUM_SETS) {
        rval = mb.add_entities( sets[id2-1], &verts[j][i], 1 );
        CHECK_ERR(rval);
      }
    }
  }
  for (int i = 0; i < MBQUAD_INT; ++i) for(int j = 0; j < MBQUAD_INT; ++j) {
    EntityHandle conn[4] = { verts[j  ][i  ], 
                               verts[j  ][i+1],
                               verts[j+1][i+1],
                               verts[j+1][i  ] };
    rval = mb.create_element( MBQUAD, conn, 4, quads[j][i] );
    CHECK_ERR(rval);
    int logical[2] = { i, j };
    rval = mb.tag_set_data( logical_tag, &quads[j][i], 1, logical );
    CHECK_ERR(rval);
    double centroid[3] = { i + 0.5, j + 0.5, 0 };
    rval = mb.tag_set_data( centroid_tag, &quads[j][i], 1, centroid );
    CHECK_ERR(rval);
    int id = i/SET_WIDTH + 1;
    if (tag_elements_with_id) {
      rval = mb.tag_set_data( id_tag, &quads[j][i], 1, &id );
      CHECK_ERR(rval);
    }
    if (create_element_sets) {
      rval = mb.add_entities( sets[id-1], &quads[j][i], 1 );
      CHECK_ERR(rval);
    }
  }
  
  if (adj_elem_tag_name && !var_len_adj_elems) {
    Tag handle_tag;
    rval = mb.tag_get_handle( adj_elem_tag_name, 
                              4,
                              MB_TYPE_HANDLE,
                              handle_tag,
                              MB_TAG_DENSE|MB_TAG_EXCL );
    CHECK_ERR(rval);
    for (int i = 0; i <= MBQUAD_INT; ++i) for(int j = 0; j <= MBQUAD_INT; ++j) {
      EntityHandle val[4] = { (i > 0        && j > 0       ) ? quads[j-1][i-1] : 0,
                                (i > 0        && j < MBQUAD_INT) ? quads[j  ][i-1] : 0,
                                (i < MBQUAD_INT && j < MBQUAD_INT) ? quads[j  ][i  ] : 0,
                                (i < MBQUAD_INT && j > 0       ) ? quads[j-1][i  ] : 0 };
      rval = mb.tag_set_data( handle_tag, &verts[j][i], 1, val );
      CHECK_ERR(rval);
    }
  }
  else if (adj_elem_tag_name && var_len_adj_elems) {
    Tag handle_tag;
    rval = mb.tag_get_handle( adj_elem_tag_name,
                              0, MB_TYPE_HANDLE,
                              handle_tag,
                              MB_TAG_DENSE|MB_TAG_VARLEN|MB_TAG_EXCL );
    CHECK_ERR(rval);
    for (int i = 0; i <= MBQUAD_INT; ++i) for(int j = 0; j <= MBQUAD_INT; ++j) {
      EntityHandle val[4];
      int num = 0;
      if (i > 0        && j > 0       ) val[num++] = quads[j-1][i-1];
      if (i > 0        && j < MBQUAD_INT) val[num++] = quads[j  ][i-1];
      if (i < MBQUAD_INT && j < MBQUAD_INT) val[num++] = quads[j  ][i  ];
      if (i < MBQUAD_INT && j > 0       ) val[num++] = quads[j-1][i  ];
      const void* ptr = val;
      rval = mb.tag_set_by_ptr( handle_tag, &verts[j][i], 1, &ptr, &num );
      CHECK_ERR(rval);
    }
  }
  
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
}

// Given a list of vertices adjacent to a quad strip, identify it as one of the 
// NUM_SETS strips of quads written by create_mesh.
int identify_set( Interface& mb, const Range& verts )
{
  const int COL = SET_WIDTH+1;
  CHECK_EQUAL( (1+MBQUAD_INT)*COL, (int)verts.size() );

    // Get X range of vertices
  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  for (Range::const_iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    ErrorCode rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
      // Expect whole-valued coorindates
    int int_x = (int)coords[0];
    CHECK( fabs( coords[0] - (double)int_x ) < 1e-12 );
    
    if (int_x < min_x) min_x = int_x;
    if (int_x > max_x) max_x = int_x;
  }
  CHECK( max_x - min_x < COL );
  
    // Calculate ID (return value) from coordinate range)
  const int ID = min_x / SET_WIDTH + 1;
  
    // Now verify that all vertices correctly form a grid
  EntityHandle grid[MBQUAD_INT+1][COL];
  memset( grid, 0, sizeof(grid) );
  for (Range::const_iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    ErrorCode rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
      // Expect whole-valued coorindates
    int x = (int)coords[0] - (ID-1)*SET_WIDTH, y = (int)coords[1];
    CHECK( fabs( coords[1] - (double)y ) < 1e-12 );
    CHECK( fabs( coords[2] ) < 1e-12 );
    CHECK( y >= 0 && y <= MBQUAD_INT );
    CHECK_EQUAL( (EntityHandle)0, grid[y][x] );
    grid[y][x] = *i;
  }
  
  return ID;
}
int identify_set( Interface& mb, EntityHandle set )
{
  ErrorCode rval;
  Range verts, elems;
  rval = mb.get_entities_by_handle( set, elems );
  CHECK_ERR(rval);
  Range::iterator it = elems.upper_bound( MBVERTEX );
  verts.merge( elems.begin(), it );
  elems.erase( elems.begin(), it );
  it = elems.lower_bound( MBENTITYSET );
  elems.erase( it, elems.end() );
  rval = mb.get_adjacencies( elems, 0, false, verts, Interface::UNION );
  CHECK_ERR(rval);
  return identify_set( mb, verts );
}

//! Read in the elems contained in a set
void test_read_one_set_elems()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  create_mesh( true, false, false, false );
  
  for (int id = 1; id <= NUM_SETS; ++id) {
    rval = mb.delete_mesh();
    CHECK_ERR(rval);
    rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, &id, 1 );
    CHECK_ERR(rval);
    Range verts;
    rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
    int act_id = identify_set( mb, verts );
    CHECK_EQUAL( id, act_id );
  }
}

//! Read in the elems contained in a sets
void test_read_two_sets_elems()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  create_mesh( true, false, false, false );
  int ids[2] = { 2, 8 };
  EntityHandle file_set;
  rval = mb.create_meshset( MESHSET_SET, file_set );
  CHECK_ERR(rval);
  rval = mb.load_file( TEST_FILE, &file_set, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);
  
  Range sets;
  rval = mb.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, (int)sets.size() );
  Range::iterator it = sets.find( file_set );
  CHECK( it != sets.end() );
  sets.erase( it );
  
  int id1 = identify_set( mb, sets.front() );
  int id2 = identify_set( mb, sets.back() );
  if (id1 == ids[0]) {
    CHECK_EQUAL( ids[1], id2 );
  }
  else {
    CHECK_EQUAL( ids[1], id1 );
    CHECK_EQUAL( ids[0], id2 );
  }
}

Tag check_tag( Interface& mb, 
                 const char* name,
                 TagType storage,
                 DataType type,
                 int size )
{
  
  Tag tag;
  ErrorCode rval = mb.tag_get_handle( name, size, type, tag );
  CHECK_ERR(rval);

  TagType storage1;
  rval = mb.tag_get_type( tag, storage1 );
  CHECK_ERR(rval);
  CHECK_EQUAL( storage, storage1 );
  
  DataType type1;
  rval = mb.tag_get_data_type( tag, type1 );
  CHECK_ERR(rval);
  CHECK_EQUAL( type, type1 );
  
  int size1;
  rval = mb.tag_get_length( tag, size1 );
  if (size <= 0) { // variable-length tag
    CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, rval );
  }
  else {
    CHECK_ERR(rval);
    CHECK_EQUAL( size, size1 );
  }
  
  return tag;
}

//! Test reading of sparse double tag data
void test_read_double_tag()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  create_mesh( true, false, false, false );
  int ids[2] = { 1, 4 };
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);
  
  Tag tag = check_tag( mb, CENTROID_NAME, MB_TAG_DENSE, MB_TYPE_DOUBLE, 3 );
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3], data[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    rval = mb.tag_get_data( tag, &*i, 1, data );
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL( coords[0], data[0], 1e-12 );
    CHECK_REAL_EQUAL( coords[1], data[1], 1e-12 );
    CHECK_REAL_EQUAL( coords[2], data[2], 1e-12 );
  }
}

//! Test reading of sparse opaque tag data
void test_read_opaque_tag()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  create_mesh( true, false, false, false );
  int ids[2] = { 1, 4 };
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);

  Tag tag = check_tag( mb, LOGICAL_NAME, MB_TAG_DENSE, MB_TYPE_OPAQUE, 2*sizeof(int) );
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    int data[2];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    rval = mb.tag_get_data( tag, &*i, 1, data );
    CHECK_ERR(rval);
    CHECK_REAL_EQUAL( coords[0], (double)data[0], 1e-12 );
    CHECK_REAL_EQUAL( coords[1], (double)data[1], 1e-12 );
  }
}

static void test_read_handle_tag_common( bool var_len )
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  const char tag_name[] = "VTX_ADJ";
  create_mesh( true, false, false, false, tag_name, var_len );
  int ids[2] = { 7, 10 };
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);
  
  Tag tag = check_tag( mb, tag_name, MB_TAG_DENSE, MB_TYPE_HANDLE, 
                         var_len ? 0 : 4 );
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    std::vector<EntityHandle> adj, val;
    rval = mb.get_adjacencies( &*i, 1, 2, false, adj, Interface::UNION );
    CHECK_ERR(rval);
    CHECK(!adj.empty());
    
    int num;
    const void* ptr;
    rval = mb.tag_get_by_ptr( tag, &*i, 1, &ptr, &num );
    CHECK_ERR(rval);
    
    if (var_len) {
      CHECK( num > 0 );
      CHECK( num < 5 );
    }
    else {
      CHECK_EQUAL( 4, num );
    }
    
    val.clear();
    const EntityHandle* dat = (const EntityHandle*)ptr;
    for (const EntityHandle* end = dat+num; dat != end; ++dat)
      if (*dat)
        val.push_back(*dat);

    CHECK_EQUAL( adj.size(), val.size() );
    std::sort( adj.begin(), adj.end() );
    std::sort( val.begin(), val.end() );
    CHECK( adj == val );
  }
}


void test_read_tagged_elems()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  create_mesh( false, false, true, false );
  int id = 5;
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, &id, 1 );
  CHECK_ERR(rval);
  
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  int id2 = identify_set( mb, verts );
  CHECK_EQUAL( id, id2 );
  
  int elems;
  rval = mb.get_number_entities_by_type( 0, MBQUAD, elems );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBQUAD_INT*MBQUAD_INT/NUM_SETS, elems );
}

void test_read_tagged_nodes()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
  create_mesh( false, false, false, true );
  int id = 1; // NOTE: this test will only succeed for ID == 1 
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, &id, 1 );
  CHECK_ERR(rval);
  
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  int id2 = identify_set( mb, verts );
  CHECK_EQUAL( id, id2 );
  
  int elems;
  rval = mb.get_number_entities_by_type( 0, MBQUAD, elems );
  CHECK_ERR(rval);
  CHECK_EQUAL( MBQUAD_INT*MBQUAD_INT/NUM_SETS, elems );
}


//! Read in the polyhedra contained in a set
void test_read_one_set_polyhedra()
{
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
    // create a 2x2x1 block of hexes, splitting each hex face
    // into two triangles to form an 12-sided polyhedron
  EntityHandle verts[18], hexes[4];
  double coords[18][3] = { {0, 0, 0},
                           {1, 0, 0},
                           {2, 0, 0},
                           {0, 1, 0},
                           {1, 1, 0},
                           {2, 1, 0},
                           {0, 0, 1},
                           {1, 0, 1},
                           {2, 0, 1},
                           {0, 1, 1},
                           {1, 1, 1},
                           {2, 1, 1},
                           {0, 0, 2},
                           {1, 0, 2},
                           {2, 0, 2},
                           {0, 1, 2},
                           {1, 1, 2},
                           {2, 1, 2} };
  int hexconn[4][8] = {  { 0, 1, 4, 3, 6, 7, 10, 9 },
                         { 1, 2, 5, 4, 7, 8, 11, 10 },
                         { 6, 7, 10, 9, 12, 13, 16, 15 },
                         { 7, 8, 11, 10, 13, 14, 17, 16 } };
  for (int i = 0; i < 18; ++i) {
    rval = mb.create_vertex( coords[i], verts[i] );
    CHECK_ERR(rval);
  }
  for (int i = 0; i < 4; ++i) {
    EntityHandle conn[8];
    for (int j = 0; j < 8; ++j)
      conn[j] = verts[hexconn[i][j]];
    rval = mb.create_element( MBHEX, conn, 8, hexes[i] );
    CHECK_ERR(rval);
  }
  
  Tag tri_tag;
  rval = mb.tag_get_handle( "tris", 2, MB_TYPE_HANDLE, tri_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> quads;
  EntityHandle tris[12], poly[4];
  for (int i = 0; i < 4; ++i) {
    quads.clear();
    rval = mb.get_adjacencies( &hexes[i], 1, 2, true, quads );
    CHECK_ERR(rval);
    CHECK_EQUAL( 6, (int)quads.size() );
    
    for (int j = 0; j < 6; ++j) {
      rval = mb.tag_get_data( tri_tag, &quads[j], 1, tris + 2*j );
      if (MB_SUCCESS == rval)
        continue;
      CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
      const EntityHandle* conn;
      int len;
      rval = mb.get_connectivity( quads[j], conn, len );
      CHECK_ERR(rval);
      CHECK_EQUAL( 4, len );
      EntityHandle tri_conn[2][3] = {  { conn[0], conn[1], conn[2] },
                                         { conn[2], conn[3], conn[0] } };
      rval = mb.create_element( MBTRI, tri_conn[0], 3, tris[2*j  ] ); CHECK_ERR(rval);
      rval = mb.create_element( MBTRI, tri_conn[1], 3, tris[2*j+1] ); CHECK_ERR(rval);
      rval = mb.tag_set_data( tri_tag, &quads[j], 1, tris + 2*j ); CHECK_ERR(rval);
    }
    
    rval = mb.create_element( MBPOLYHEDRON, tris, 12, poly[i] );
    CHECK_ERR(rval);
  }
  
  Range all_tri;
  rval = mb.get_entities_by_type( 0, MBTRI, all_tri );
  CHECK_ERR(rval);
  CHECK_EQUAL( 40, (int)all_tri.size() );
  
  rval = mb.delete_entities( hexes, 4 ); CHECK_ERR(rval);
  rval = mb.delete_entities( &quads[0], quads.size() ); CHECK_ERR(rval);
  
  EntityHandle sets[2];
  rval = mb.create_meshset( 0, sets[0] ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[0], poly, 2 ); CHECK_ERR(rval);
  rval = mb.create_meshset( 0, sets[1] ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], poly+2, 2 ); CHECK_ERR(rval);
  
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  int ids[2] = { 2, 3 };
  rval = mb.tag_set_data( id_tag, sets, 2, ids ); CHECK_ERR(rval);
  
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  rval = mb.delete_mesh();
  CHECK_ERR(rval);
  
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, 1 );
  CHECK_ERR(rval);
  
  Range rpoly;
  rval = mb.get_entities_by_type( 0, MBPOLYHEDRON, rpoly );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)rpoly.size() );
  
  Range polyverts;
  rval = mb.get_adjacencies( rpoly, 0, false, polyverts, Interface::UNION );
  CHECK_ERR(rval);
  CHECK_EQUAL( 12, (int)polyverts.size() );
  
  for (Range::iterator it = polyverts.begin(); it != polyverts.end(); ++it) {
    double coords2[3];
    rval = mb.get_coords( &*it, 1, coords2 );
    CHECK_ERR(rval);
    CHECK( coords2[0] > -1e-12 && coords2[0]-2 < 1e-12 );
    CHECK( coords2[1] > -1e-12 && coords2[1]-1 < 1e-12 );
    CHECK( coords2[2] > -1e-12 && coords2[2]-1 < 1e-12 );
  }
}

//! Read in the sets contained in a set.  
//! Should read all sets containing read elements or nodes
//! and all sets that are contained the the specified "read"
//! set.  Test the later here.
void test_read_set_sets()
{
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
    // create sets and assign an ID to each
  const int len = 5;
  EntityHandle set[2*len+2];
  for (int i = 0; i < 2*len+2; ++i) {
    rval = mb.create_meshset( MESHSET_SET, set[i] );
    CHECK_ERR(rval);
    int id = i + 1;
    rval = mb.tag_set_data( id_tag, set + i, 1, &id );
    CHECK_ERR(rval);
  }
  
    // make set containment as follows (values are assigned IDs):
  int cont_ids[2][len] = { { 3, 4, 5, 9, 10 }, { 6, 7, 8, 11, 12 } };
  for (int i = 0; i < 2; ++i) {
    EntityHandle contents[len] = { set[cont_ids[i][0] - 1],
                                     set[cont_ids[i][1] - 1],
                                     set[cont_ids[i][2] - 1],
                                     set[cont_ids[i][3] - 1],
                                     set[cont_ids[i][4] - 1] };
    rval = mb.add_entities( set[i], contents, len );
    CHECK_ERR(rval);
  }
  
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  
  for (int i = 0; i < 2; ++i) {
    rval = mb.delete_mesh();
    CHECK_ERR(rval);
    
    EntityHandle file;
    rval = mb.create_meshset( MESHSET_SET, file );
    CHECK_ERR(rval);
    int id = i+1;
    rval = mb.load_file( TEST_FILE, &file, READ_OPTS ";SETS=NONE", ID_TAG_NAME, &id, 1 );
    CHECK_ERR(rval);
    
      // check that the total number of sets read is as expected
    Range sets;
    rval = mb.get_entities_by_type( 0, MBENTITYSET, sets );
    CHECK_ERR(rval);
    Range::iterator it = sets.find( file );
    if (it != sets.end())
      sets.erase( it );
    CHECK_EQUAL( len+1, (int)sets.size() );
    
      // check that we read in the set specified by ID to the reader
    rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
    CHECK_ERR(rval);
    sets.clear();
    const void* data[] = { &id };
    rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &id_tag, data, 1, sets );
    CHECK_ERR(rval);
    CHECK_EQUAL( 1, (int)sets.size() );
      
      // check that it contains the expected sets
    EntityHandle owner = sets.front();
    sets.clear();
    rval = mb.get_entities_by_type( owner, MBENTITYSET, sets );
    CHECK_ERR(rval);
    CHECK_EQUAL( len, (int)sets.size() );
    
    std::vector<int> expected( cont_ids[i], cont_ids[i] + len );
    std::vector<int> actual( len );
    rval = mb.tag_get_data( id_tag, sets, &actual[0] );
    CHECK_ERR(rval);
    std::sort( expected.begin(), expected.end() );
    std::sort( actual.begin(), actual.end() );
    CHECK( expected == actual );
  }
}

static void check_children( bool contents, GatherTestMode mode, Interface& mb, int id, Tag id_tag, EntityHandle file )
{
    // Increase number of expected sets by one if contents is true because
    // we always read immediately contained (depth 1) sets.
  const int exp_num_sets = (mode == GATHER_NONE) ? 1+contents : id;
  const int exp_num_edges = (mode == GATHER_CONTENTS) ? id : 1;
  
  ErrorCode rval;
  Range range;
  rval = mb.get_entities_by_type( 0, MBEDGE , range );
  CHECK_ERR(rval);
  CHECK_EQUAL( exp_num_edges, (int)range.size() );
  range.clear();
  rval = mb.get_entities_by_type( 0, MBENTITYSET , range );
  CHECK_ERR(rval);
  Range::iterator it = range.find( file );
  CHECK( it != range.end() );
  range.erase( it );
  CHECK_EQUAL( exp_num_sets, (int)range.size() );
  
  EntityHandle set;
  const void* val[] = {&id};
  range.clear();
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &id_tag, val, 1, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() );
  set = range.front();
  
  if (mode == GATHER_NONE) {
    range.clear();
    rval = mb.get_entities_by_type( set, MBEDGE , range );
    CHECK_ERR(rval);
    CHECK_EQUAL( 1, (int)range.size() );
    return;
  }
  
  for (int i = id; i > 0; --i) {
    int act_id;
    rval = mb.tag_get_data( id_tag, &set, 1, &act_id );
    CHECK_ERR(rval);
    CHECK_EQUAL( i, act_id );

    range.clear();
    rval = mb.get_entities_by_type( set, MBEDGE, range );
    CHECK_ERR(rval);
    if (mode == GATHER_CONTENTS || i == id) {
      CHECK_EQUAL( 1, (int)range.size() );
      const EntityHandle* conn;
      int len;
      rval = mb.get_connectivity( range.front(), conn, len );
      CHECK_ERR(rval);
      CHECK_EQUAL( 2, len );
      double coords[3];
      rval = mb.get_coords( conn + 1, 1, coords );
      CHECK_ERR(rval);
      CHECK_EQUAL( i, (int)coords[0] );
    }
    else {
      CHECK( range.empty() );
    }
    
    std::vector<EntityHandle> children;
    if (contents)
      rval = mb.get_entities_by_type( set, MBENTITYSET, children );
    else
      rval = mb.get_child_meshsets( set, children );
    CHECK_ERR(rval);
    if (i == 1) {
      CHECK( children.empty() );
    }
    else {
      CHECK_EQUAL( 1, (int)children.size() );
      set = children[0];
    }
  } 
}


const char* set_read_opts[] = { "SETS", "CONTENTS", "NONE" };
void test_gather_sets_common( bool contents, GatherTestMode mode, bool no_parent_containing_sets )
{
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
    // Create a string of edges from [0,INT] along the X axis each 1 unit in length.
    // Create a set for edge edge containing the edge and make it the parent of the 
    // set containing the previous (closer to origin) edge.  Assign each set an
    // ID that is the X coordinate of the larger of the two vertices of the edge
    // contained in the set.
  const int INT = 64;
  EntityHandle verts[INT+1], edges[INT], sets[INT];
  double coords[] = { 0, 0, 0 };
  rval = mb.create_vertex( coords, verts[0] );
  CHECK_ERR(rval);
  for (int i = 0; i < INT; ++i) {
    const int id = i + 1;
    coords[0] = id;
    rval = mb.create_vertex( coords, verts[id] );
    CHECK_ERR(rval);
    rval = mb.create_element( MBEDGE, verts + i, 2, edges[i] );
    CHECK_ERR(rval);
    rval = mb.create_meshset( MESHSET_SET, sets[i] );
    CHECK_ERR(rval);
    rval = mb.add_entities( sets[i], edges + i, 1 );
    CHECK_ERR(rval);
    rval = mb.tag_set_data( id_tag, sets + i, 1, &id );
    CHECK_ERR(rval);
    if (i > 0) {
      if (contents) 
        rval = mb.add_entities( sets[i], sets + (i-1), 1 );
      else
        rval = mb.add_child_meshset( sets[i], sets[i-1] );
      CHECK_ERR(rval);
    }
  }
  
    // Write the data
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
 
  EntityHandle file;
  std::string opt( READ_OPTS );
  if (contents)
    opt += ";CHILDREN=NONE;SETS=";
  else
    opt += ";SETS=NONE;CHILDREN=";
  opt += set_read_opts[mode];

  if (no_parent_containing_sets) opt += ";NO_SET_CONTAINING_PARENTS";

  const int test_ids[] = { 2, 7, INT/3-1, INT/2+1, INT-3 };
  const int num_test_ids = sizeof(test_ids)/sizeof(int);
  for (int i = 0; i < num_test_ids; ++i) {
    CHECK (test_ids[i] <= INT);
    
    rval = mb.delete_mesh();
    CHECK_ERR(rval);
    
    rval = mb.create_meshset( MESHSET_SET, file );
    CHECK_ERR(rval);
    
    rval = mb.load_file( TEST_FILE, 0, opt.c_str(), ID_TAG_NAME, test_ids+i, 1 );
    CHECK_ERR(rval);
    rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
    CHECK_ERR(rval);
    
    check_children( contents, mode, mb, test_ids[i], id_tag, file );
  }
}


void test_gather_sets_ranged( bool contents, GatherTestMode mode, bool no_parent_containing_sets )
{
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
  Range verts;
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
    // create four groups of vertices, where all vertices in the same group
    // have the same x-coordinate
  const int NUM_GRP_VTX = 20;
  const int NUM_GRP = 4;
  EntityHandle sets[NUM_GRP];
  for (int i = 0; i < NUM_GRP; ++i) {
    double coords[3*NUM_GRP_VTX];
    for (int j = 0; j < NUM_GRP_VTX; ++j) {
      coords[3*j  ] = i;
      coords[3*j+1] = j;
      coords[3*j+2] = 0;
    }
    rval = mb.create_vertices( coords, NUM_GRP_VTX, verts );
    CHECK_ERR(rval);
    
    rval = mb.create_meshset( MESHSET_SET, sets[i] );
    CHECK_ERR(rval);
    rval = mb.add_entities( sets[i], verts );
    CHECK_ERR(rval);
    int id = i + 1;
    rval = mb.tag_set_data( id_tag, sets+i, 1, &id );
    CHECK_ERR(rval);
  }
  
    // place two of the sets inside the others
  if (contents) {
    rval = mb.add_entities( sets[0], &sets[1], 1 ); CHECK_ERR(rval);
    rval = mb.add_entities( sets[2], &sets[3], 1 ); CHECK_ERR(rval);
  }
  else {
    rval = mb.add_child_meshset( sets[0], sets[1] ); CHECK_ERR(rval);
    rval = mb.add_child_meshset( sets[2], sets[3] ); CHECK_ERR(rval);
  }
    
    // Write the data
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
 
    // Read the data
  std::string opt( READ_OPTS );
  if (contents)
    opt += ";CHILDREN=NONE;SETS=";
  else
    opt += ";SETS=NONE;CHILDREN=";
  opt += set_read_opts[mode];

  if (no_parent_containing_sets) opt += ";NO_PARENT_CONTAINING_SETS";

  EntityHandle file;
  const int read_id = 3;
  rval = mb.delete_mesh(); CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, file ); CHECK_ERR(rval);
  rval = mb.load_file( TEST_FILE, &file, opt.c_str(), ID_TAG_NAME, &read_id, 1 );
  CHECK_ERR(rval);
  
    // get any sets that were read it
  Range read_sets;
  rval = mb.get_entities_by_type( file, MBENTITYSET, read_sets );
  CHECK_ERR(rval);
  
    // count number of vertices in each group
  int counts[NUM_GRP];
  memset( counts, 0, sizeof(counts) );
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  for (Range::iterator it = verts.begin(); it != verts.end(); ++it) {
    double coords[3];
    rval = mb.get_coords( &*it, 1, coords );
    CHECK_ERR(rval);
    int i = (int)(coords[0]+1e-12);
    CHECK( i >= 0 && i < NUM_GRP );
    counts[i]++;
  }
  
    // check expected counts
  CHECK_EQUAL( 0, counts[0] );
  CHECK_EQUAL( 0, counts[1] );
  CHECK_EQUAL( NUM_GRP_VTX, counts[2] );
  switch (mode) {
    case GATHER_NONE:
      CHECK_EQUAL( 0, counts[3] );
      CHECK_EQUAL( 1+contents, (int)read_sets.size() );
      break;
    case GATHER_SETS:
      CHECK_EQUAL( 0, counts[3] );
      CHECK_EQUAL( 2, (int)read_sets.size() );
      break;
    case GATHER_CONTENTS:
      CHECK_EQUAL( NUM_GRP_VTX, counts[3] );
      CHECK_EQUAL( 2, (int)read_sets.size() );
      break;
  }
}

static void check_num_verts( Interface& mb, Tag tag, int id, int num_vtx )
{
  ErrorCode rval;
  const void* val[] = {&id};
  Range range;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, val, 1, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() );

  EntityHandle set = range.front();
  range.clear();
  rval = mb.get_entities_by_type( set, MBVERTEX, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_vtx, (int)range.size() );
}


//! Read in the sets contained in a set.  
//! Should read all sets containing read elements or nodes
//! and all sets that are contained the the specified "read"
//! set.  Test the former here.
void test_read_containing_sets()
{
    // create mesh decomposed by elements but create
    // sets containing all vertices of decomposed elements
    // such that adjacent sets share vertices.
  create_mesh( false, true, false, false );
  
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
    // read some sets
  const int ids[] = { 1, 5, 9 };
  const int num_sets = sizeof(ids)/sizeof(int);
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, num_sets );
  CHECK_ERR(rval);

  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag );
  CHECK_ERR(rval);
  
    // expect all sets adjacent to the specified sets because
    // they share vertices.
  Range verts;
  for (int i = 0; i < num_sets; ++i) {
    if (ids[i] > 1)
      check_num_verts( mb, id_tag, ids[i]-1, MBQUAD_INT+1 );
    check_num_verts( mb, id_tag, ids[i], (MBQUAD_INT+1)*(SET_WIDTH+1) );
    if (ids[i] < NUM_SETS)
      check_num_verts( mb, id_tag, ids[i]+1, MBQUAD_INT+1 );
  }
}

//! Test reading of explicit adjacencies
void test_read_adjacencies()
{
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
    // create four hexes sharing an edge
  EntityHandle verts[3][3][2], hexes[2][2];
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double coords[] = { i, j, k };
        rval = mb.create_vertex( coords, verts[i][j][k] );
        CHECK_ERR(rval);
      }
    }
  }
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      EntityHandle conn[] = { verts[i  ][j  ][0],
                                verts[i+1][j  ][0],
                                verts[i+1][j+1][0],
                                verts[i  ][j+1][0],
                                verts[i  ][j  ][1],
                                verts[i+1][j  ][1],
                                verts[i+1][j+1][1],
                                verts[i  ][j+1][1] };
      rval = mb.create_element( MBHEX, conn, 8, hexes[i][j] );
      CHECK_ERR(rval);
    }
  }
  
    // create two duplicate edges that connect the vertices common to all four hexes
  EntityHandle edge_conn[2] = { verts[1][1][0], verts[1][1][1] };
  EntityHandle edges[2];
  rval = mb.create_element( MBEDGE, edge_conn, 2, edges[0] ); CHECK_ERR(rval);
  rval = mb.create_element( MBEDGE, edge_conn, 2, edges[1] ); CHECK_ERR(rval);
    // mark one edge as adjacent to the left two hexes and the
    // other as adjacent to the right two
  rval = mb.add_adjacencies( edges[0], hexes[0], 2, true ); CHECK_ERR(rval);
  rval = mb.add_adjacencies( edges[1], hexes[1], 2, true ); CHECK_ERR(rval);
    // create two sets containing the front two and the rear two
    // hexes, respectively.
  EntityHandle sets[2];
  rval = mb.create_meshset( MESHSET_SET, sets[0] ); CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, sets[1] ); CHECK_ERR(rval);
  EntityHandle set1[4] = { hexes[0][0], hexes[1][0], edges[0], edges[1] };
  EntityHandle set2[4] = { hexes[0][1], hexes[1][1], edges[0], edges[1] };
  rval = mb.add_entities( sets[0], set1, 4 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], set2, 4 ); CHECK_ERR(rval);
  
    // assign IDs to sets
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  int ids[2] = { 1, 2 };
  rval = mb.tag_set_data( id_tag, sets, 2, ids );
  CHECK_ERR(rval);
  
    // write mesh
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  
    // read mesh
  rval = mb.delete_mesh(); CHECK_ERR(rval);
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, 1 );
  CHECK_ERR(rval);
  
    // expect two hexes and two edges
  Range range;
  rval = mb.get_entities_by_type( 0, MBHEX, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)range.size() );
  EntityHandle h1 = range.front(), h2 = range.back();
  range.clear();
  rval = mb.get_entities_by_type( 0, MBEDGE, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)range.size() );
 
    // expecte each edge to have one of the hexes
  range.clear();
  rval = mb.get_adjacencies( &h1, 1, 1, false, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() ); 
  EntityHandle e1 = range.front();
  range.clear();
  rval = mb.get_adjacencies( &h2, 1, 1, false, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() ); 
  EntityHandle e2 = range.front();
  
  CHECK( e1 != e2 );
}


void test_read_sides()
{
  ErrorCode rval;
  Core instance;
  Interface& mb = instance;
  
    // create 4x4 grid of quads with edges
  const int INT = 4;
  EntityHandle verts[INT+1][INT+1];
  for (int j = 0; j <= INT; ++j) {
    for (int i = 0; i <= INT; ++i) {
      double coords[3] = { i, j, 0 };
      rval = mb.create_vertex( coords, verts[INT-j][i] );
      CHECK_ERR(rval);
    }
  }
  EntityHandle quads[INT][INT];
  for (int j = 0; j < INT; ++j) {
    for (int i = 0; i < INT; ++i) {
      EntityHandle conn[4] = { verts[INT-j][i],
                                 verts[INT-j][i+1],
                                 verts[INT-j-1][i+1],
                                 verts[INT-j-1][i] };
      rval = mb.create_element( MBQUAD, conn, 4, quads[INT-j-1][i] );
      CHECK_ERR(rval);
    }
  }
  Range edges;
  rval = mb.get_adjacencies( &quads[0][0], INT*INT, 1, true, edges, Interface::UNION );
  CHECK_ERR(rval);
  CHECK_EQUAL( 40, (int)edges.size() );
  
    // group quads into two sets
  EntityHandle sets[2];
  rval = mb.create_meshset( MESHSET_SET, sets[0] ); CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, sets[1] ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[0], quads[0], INT ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], quads[1], INT ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[0], quads[2], INT ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], quads[3], INT ); CHECK_ERR(rval);
  
    // assign IDS
  Tag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  int ids[2] = { 4, 5 };
  rval = mb.tag_set_data( id_tag, sets, 2, ids );
  CHECK_ERR(rval);
  
    // write mesh
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  
    // read first set back in
  rval = mb.delete_mesh(); CHECK_ERR(rval);
  rval = mb.load_file( TEST_FILE, 0, READ_OPTS, ID_TAG_NAME, ids, 1 );
  CHECK_ERR(rval);
  
    // check expected counts
  int count;
  rval = mb.get_number_entities_by_type( 0, MBVERTEX, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( (INT+1)*INT, count );
  rval = mb.get_number_entities_by_type( 0, MBQUAD, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( INT*INT/2, count );
  rval = mb.get_number_entities_by_type( 0, MBEDGE, count );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2*(INT+1) + INT*INT, count );

    // check edges adjacent to each quad
  Range elems;
  rval = mb.get_entities_by_type( 0, MBQUAD, elems );
  CHECK_ERR(rval);
  for (Range::iterator it = elems.begin(); it != elems.end(); ++it) {
    edges.clear();
    rval = mb.get_adjacencies( &*it, 1, 1, false, edges );
    CHECK_ERR(rval);
    CHECK_EQUAL( 4, (int)edges.size() );
  }
}

const int expected_ids[] = { 2, 4, 6, 8, 10, 12, 14, 16, 18 };
const int expected_vols[] = { 3, 7, 10 };

void write_id_test_file()
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // create 12 entity sets
  EntityHandle sets[12];
  for (int i = 0; i < 12; ++i) {
    rval = mb.create_meshset( MESHSET_SET, sets[i] );
    CHECK_ERR(rval);
  }
  
    // create tag handles
  Tag id = 0, gid = 0, dim = 0;
  mb.tag_get_handle( ID_TAG_NAME, 1, MB_TYPE_INTEGER, id, MB_TAG_SPARSE|MB_TAG_EXCL );
  mb.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, dim, MB_TAG_SPARSE|MB_TAG_EXCL );
  mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE|MB_TAG_EXCL );
  
    // set ID tag on first 10 sets
  rval = mb.tag_set_data( id, sets, sizeof(expected_ids)/sizeof(int), expected_ids );
  CHECK_ERR(rval);
    // set geom dim on all sets, only three of them have dim == 3
  int num_vol = sizeof(expected_vols)/sizeof(int);
  int dims[12], ids[12];
  int v = 0;
  for (int i = 0; i < 12; ++i) {
    dims[i] = i % 3 + 1;
    if (dims[i] == 3) {
      if (v < num_vol) 
        ids[i] = expected_vols[v++];
      else
        ids[i] = expected_vols[0];
    }
    else
      ids[i] = 100;
  }
  rval = mb.tag_set_data( gid, sets, 12, ids );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( dim, sets, 12, dims );
  CHECK_ERR(rval);
  
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
}

void test_read_ids()
{
  write_id_test_file();
  
  Core moab;
  ReadHDF5 reader(&moab);
  FileOptions opts("");
  ErrorCode rval;
  std::vector<int> values;
  rval = reader.read_tag_values( TEST_FILE, ID_TAG_NAME, opts, values );
  remove( TEST_FILE );
  CHECK_ERR(rval);
  
  std::sort( values.begin(), values.end() );
  std::vector<int> expected( expected_ids, expected_ids+sizeof(expected_ids)/sizeof(int) );
  CHECK_EQUAL( expected, values );
}

void test_read_partial_ids()
{
  write_id_test_file();
  
  const int three = 3;
  ReaderIface::IDTag vols = { GEOM_DIMENSION_TAG_NAME, &three, 1 };
  ReaderIface::SubsetList subset = { &vols, 1, 0, 0 };
  
  Core moab;
  ReadHDF5 reader(&moab);
  FileOptions opts("");
  ErrorCode rval;
  std::vector<int> values;
  rval = reader.read_tag_values( TEST_FILE, GLOBAL_ID_TAG_NAME, opts, values, &subset );
  remove( TEST_FILE );
  CHECK_ERR(rval);
  
  std::sort( values.begin(), values.end() );
  std::vector<int> expected( expected_ids, expected_ids+sizeof(expected_ids)/sizeof(int) );
}
