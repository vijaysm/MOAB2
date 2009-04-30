#include "MBCore.hpp"
#include "MBRange.hpp"
#include "TestUtil.hpp"
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <limits>

const char TEST_FILE[] = "partial.h5m";
const char READ_OPTS[] = "BUFFER_SIZE=256";
const char ID_TAG_NAME[] = "test_id_tag";


static void test_read_nothing_common( bool non_existant );
static void test_read_nodes_common( int num_read_sets );
static void test_read_handle_tag_common( bool var_len );

const int QUAD_INT = 20; 
const int NUM_SETS = 10;
const int SET_WIDTH = (QUAD_INT + NUM_SETS - 1) / NUM_SETS; // ceil(QUAD_INT/NUM_SETS)
const char LOGICAL_NAME[] = "logical";   // tag storing logical (i,j) coordinates
const char CENTROID_NAME[] = "centroid"; // tag storing position of centroid (x,y,0)
//! Create a regular QUAD_INT^2 element quad mesh with regularly
//! spaced coordinates in the range [1,100].  Group elements
//! into 10 vertical strips QUAD_INT/10 elements wide.  Tag elements,
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
int identify_set( MBInterface& mb, const MBRange& verts );
int identify_set( MBInterface& mb, MBEntityHandle set );

static MBTag check_tag( MBInterface& mb, 
                        const char* name,
                        MBTagType storage,
                        MBDataType type,
                        int size );

enum ChildTestMode { CHILD_SETS, CHILD_CONTENTS, CHILD_NONE };
void test_read_children_common( ChildTestMode mode );


//! Read a set containing no entities
void test_read_empty_set()
  { test_read_nothing_common( false ); }
  
//! Specify ID that doesn't exist in file
void test_read_non_existant_set()
  { test_read_nothing_common( true ); }

//! Read in the nodes contained in a set.
void test_read_one_set_nodes()
  { test_read_nodes_common(1); }

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
  { test_read_nodes_common(2); }

//! Read in the elems contained in a sets
void test_read_two_sets_elems();

//! For any set selected to be read by either explicit designation,
//! containing read entities, or contained in an explcitly designated
//! set, any child sets are also read.  Check that here.
void test_read_child_sets_only()
{ test_read_children_common( CHILD_SETS ); }
void test_read_child_set_contents()
{ test_read_children_common( CHILD_CONTENTS ); }
void test_read_no_child_sets()
{ test_read_children_common( CHILD_NONE ); }

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



int main( int argc, char* argv[] )
{
  if (argc > 1) {
    if (argc != 2 || strcmp(argv[1],"-k")) {
      std::cerr << "Usage: " << argv[0] << " [-k]" << std::endl;
      return 1;
    }
  }
  
  int result = 0;
  
  result += RUN_TEST(test_read_empty_set);
  result += RUN_TEST(test_read_non_existant_set);
  result += RUN_TEST(test_read_one_set_nodes);
  result += RUN_TEST(test_read_one_set_elems);
  result += RUN_TEST(test_read_one_set_polyhedra);
  result += RUN_TEST(test_read_set_sets);
  result += RUN_TEST(test_read_two_sets_nodes);
  result += RUN_TEST(test_read_two_sets_elems);
  result += RUN_TEST(test_read_child_sets_only);
  result += RUN_TEST(test_read_child_set_contents);
  result += RUN_TEST(test_read_no_child_sets);
  result += RUN_TEST(test_read_containing_sets);
  result += RUN_TEST(test_read_double_tag);
  result += RUN_TEST(test_read_opaque_tag);
  result += RUN_TEST(test_read_handle_tag);
  result += RUN_TEST(test_var_len_tag);
  result += RUN_TEST(test_read_adjacencies);
  result += RUN_TEST(test_read_tagged_elems);
  result += RUN_TEST(test_read_tagged_nodes);

  if (argc == 1)
    remove( TEST_FILE );
  return result;
}

void test_read_nothing_common( bool non_existant )
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;

  // create a few nodes to write to file
  std::vector<double> coords( 3000 );
  MBRange verts;
  rval= mb.create_vertices( &coords[0], coords.size()/3, verts );
  CHECK_ERR(rval);
  
  // create three entity sets
  MBEntityHandle sets[3];
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
  MBTag id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
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
  MBEntityHandle file_set;
  int id = non_existant ? 8 : 7;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, &id, 1 );
  CHECK_ERR( rval );
  
  // the file should contain exactly two sets (the specified one and the new
  // file set, and nothing else.)
  for (MBEntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
    int count = -1;
    rval = mb.get_number_entities_by_type( 0, t, count );
    CHECK_ERR(rval);
    CHECK_EQUAL( 0, count );
  }
  MBRange setrange;
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

void test_read_nodes_common( int num_read_sets )
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
    
    // create 1000 nodes
  const int num_sets = 2*num_read_sets;
  std::vector<MBEntityHandle> verts(1000);
  std::vector< std::vector<MBEntityHandle> > set_verts(num_sets);
  for (size_t i = 0; i < verts.size(); ++i) {
    double coords[3];
    int j = i % num_sets;
    vtx_coords( j+1, set_verts[j].size(), num_sets, coords );
    rval = mb.create_vertex( coords, verts[i] );
    set_verts[ j ].push_back( verts[i] );
    CHECK_ERR(rval);
  }
  
    // create two sets, each containing half of the nodes
  std::vector<MBEntityHandle> sets(num_sets);
  for (int i = 0; i < num_sets; ++i) {
    rval = mb.create_meshset( MESHSET_ORDERED, sets[i] );
    CHECK_ERR(rval);
    rval = mb.add_entities( sets[i], &set_verts[i][0], set_verts[i].size() );
    CHECK_ERR(rval);
  }
 
    // tag both sets
  MBTag id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
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
  values.resize( num_read_sets );
  for (int i = 0; i < num_read_sets; ++i) values[i] = 2*(i+1);
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, &values[0], num_read_sets );
  CHECK_ERR(rval);
  
  int count, expected = 0;
  rval = mb.get_number_entities_by_dimension( 0, 0, count );
  CHECK_ERR(rval);
  for (int i = 0; i < num_sets; ++i)
    if (i % 2)
      expected += set_verts[i].size();
  CHECK_EQUAL( expected, count );
  
  MBRange sets2;
  rval = mb.get_entities_by_type( 0, MBENTITYSET, sets2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1+num_read_sets, (int)sets2.size() );
  MBRange::iterator it = sets2.find( file_set );
  CHECK( it != sets2.end() );
  sets2.erase( it );
  
  rval = mb.tag_get_handle( ID_TAG_NAME, id_tag );
  CHECK_ERR(rval);  
  while (!sets2.empty()) {
    MBEntityHandle set = sets2.pop_front();
    int id;
    rval = mb.tag_get_data( id_tag, &set, 1, &id );
    CHECK_ERR(rval);
    CHECK( std::find(values.begin(), values.end(), id) != values.end() );
    CHECK( id > 0 );
    CHECK( (unsigned)id <= set_verts.size() );
    
    std::vector<MBEntityHandle> verts;
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

//! Create a regular QUAD_INT^2 element quad mesh with regularly
//! spaced coordinates in the range [1,100].  Group elements
//! into 10 vertical strips QUAD_INT/10 elements wide.  Tag elements,
//! vertices and/or sets with ID in [1,10] stored in ID_TAG_NAME
//! tag.  Write new mesh to TEST_FILE.
void create_mesh( bool create_element_sets,
                  bool create_vertex_sets,
                  bool tag_elements_with_id,
                  bool tag_vertices_with_id,
                  const char* adj_elem_tag_name,
                  bool var_len_adj_elems )
{
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
    // create tags
  MBTag logical_tag, centroid_tag, id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
  CHECK_ERR(rval);
  rval = mb.tag_create( LOGICAL_NAME, 2*sizeof(int), MB_TAG_DENSE, MB_TYPE_OPAQUE, logical_tag, 0 );
  CHECK_ERR(rval);
  rval = mb.tag_create( CENTROID_NAME, 3*sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, centroid_tag, 0 );
  CHECK_ERR(rval);

  MBEntityHandle sets[NUM_SETS];
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
  MBEntityHandle verts[QUAD_INT+1][QUAD_INT+1], quads[QUAD_INT][QUAD_INT];
  for (int i = 0; i <= QUAD_INT; ++i) for(int j = 0; j <= QUAD_INT; ++j) {
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
  for (int i = 0; i < QUAD_INT; ++i) for(int j = 0; j < QUAD_INT; ++j) {
    MBEntityHandle conn[4] = { verts[j  ][i  ], 
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
    MBTag handle_tag;
    rval = mb.tag_create( adj_elem_tag_name, 
                          4 * sizeof(MBEntityHandle), 
                          MB_TAG_DENSE,
                          MB_TYPE_HANDLE,
                          handle_tag,
                          0 );
    CHECK_ERR(rval);
    for (int i = 0; i <= QUAD_INT; ++i) for(int j = 0; j <= QUAD_INT; ++j) {
      MBEntityHandle val[4] = { (i > 0        && j > 0       ) ? quads[j-1][i-1] : 0,
                                (i > 0        && j < QUAD_INT) ? quads[j  ][i-1] : 0,
                                (i < QUAD_INT && j < QUAD_INT) ? quads[j  ][i  ] : 0,
                                (i < QUAD_INT && j > 0       ) ? quads[j-1][i  ] : 0 };
      rval = mb.tag_set_data( handle_tag, &verts[j][i], 1, val );
      CHECK_ERR(rval);
    }
  }
  else if (adj_elem_tag_name && var_len_adj_elems) {
    MBTag handle_tag;
    rval = mb.tag_create_variable_length( adj_elem_tag_name,
                                          MB_TAG_DENSE,
                                          MB_TYPE_HANDLE,
                                          handle_tag );
    CHECK_ERR(rval);
    for (int i = 0; i <= QUAD_INT; ++i) for(int j = 0; j <= QUAD_INT; ++j) {
      MBEntityHandle val[4];
      int num = 0;
      if (i > 0        && j > 0       ) val[num++] = quads[j-1][i-1];
      if (i > 0        && j < QUAD_INT) val[num++] = quads[j  ][i-1];
      if (i < QUAD_INT && j < QUAD_INT) val[num++] = quads[j  ][i  ];
      if (i < QUAD_INT && j > 0       ) val[num++] = quads[j-1][i  ];
      const void* ptr = val;
      num *= sizeof(MBEntityHandle);
      rval = mb.tag_set_data( handle_tag, &verts[j][i], 1, &ptr, &num );
      CHECK_ERR(rval);
    }
  }
  
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
}

// Given a list of vertices adjacent to a quad strip, identify it as one of the 
// NUM_SETS strips of quads written by create_mesh.
int identify_set( MBInterface& mb, const MBRange& verts )
{
  const int COL = SET_WIDTH+1;
  CHECK_EQUAL( (1+QUAD_INT)*COL, (int)verts.size() );

    // Get X range of vertices
  int min_x = std::numeric_limits<int>::max();
  int max_x = std::numeric_limits<int>::min();
  for (MBRange::const_iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    MBErrorCode rval = mb.get_coords( &*i, 1, coords );
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
  MBEntityHandle grid[QUAD_INT+1][COL];
  memset( grid, 0, sizeof(grid) );
  for (MBRange::const_iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    MBErrorCode rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
      // Expect whole-valued coorindates
    int x = (int)coords[0] - (ID-1)*SET_WIDTH, y = (int)coords[1];
    CHECK( fabs( coords[1] - (double)y ) < 1e-12 );
    CHECK( fabs( coords[2] ) < 1e-12 );
    CHECK( y >= 0 && y <= QUAD_INT );
    CHECK_EQUAL( (MBEntityHandle)0, grid[y][x] );
    grid[y][x] = *i;
  }
  
  return ID;
}
int identify_set( MBInterface& mb, MBEntityHandle set )
{
  MBErrorCode rval;
  MBRange verts, elems;
  rval = mb.get_entities_by_handle( set, elems );
  CHECK_ERR(rval);
  MBRange::iterator it = elems.upper_bound( MBVERTEX );
  verts.merge( elems.begin(), it );
  elems.erase( elems.begin(), it );
  it = elems.lower_bound( MBENTITYSET );
  elems.erase( it, elems.end() );
  rval = mb.get_adjacencies( elems, 0, false, verts, MBInterface::UNION );
  CHECK_ERR(rval);
  return identify_set( mb, verts );
}

//! Read in the elems contained in a set
void test_read_one_set_elems()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  create_mesh( true, false, false, false );
  
  for (int id = 1; id <= NUM_SETS; ++id) {
    rval = mb.delete_mesh();
    CHECK_ERR(rval);
    MBEntityHandle file_set;
    rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, &id, 1 );
    CHECK_ERR(rval);
    MBRange verts;
    rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
    int act_id = identify_set( mb, verts );
    CHECK_EQUAL( id, act_id );
  }
}

//! Read in the elems contained in a sets
void test_read_two_sets_elems()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  create_mesh( true, false, false, false );
  int ids[2] = { 2, 8 };
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);
  
  MBRange sets;
  rval = mb.get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 3, (int)sets.size() );
  MBRange::iterator it = sets.find( file_set );
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

MBTag check_tag( MBInterface& mb, 
                 const char* name,
                 MBTagType storage,
                 MBDataType type,
                 int size )
{
  
  MBTag tag;
  MBErrorCode rval = mb.tag_get_handle( name, tag );
  CHECK_ERR(rval);

  MBTagType storage1;
  rval = mb.tag_get_type( tag, storage1 );
  CHECK_ERR(rval);
  CHECK_EQUAL( storage, storage1 );
  
  MBDataType type1;
  rval = mb.tag_get_data_type( tag, type1 );
  CHECK_ERR(rval);
  CHECK_EQUAL( type, type1 );
  
  int size1;
  rval = mb.tag_get_size( tag, size1 );
  if (size < 0) { // variable-length tag
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
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  create_mesh( true, false, false, false );
  int ids[2] = { 1, 4 };
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);
  
  MBTag tag = check_tag( mb, CENTROID_NAME, MB_TAG_DENSE, MB_TYPE_DOUBLE, 3*sizeof(double) );
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  for (MBRange::iterator i = verts.begin(); i != verts.end(); ++i) {
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
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  create_mesh( true, false, false, false );
  int ids[2] = { 1, 4 };
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);

  MBTag tag = check_tag( mb, LOGICAL_NAME, MB_TAG_DENSE, MB_TYPE_OPAQUE, 2*sizeof(int) );
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  for (MBRange::iterator i = verts.begin(); i != verts.end(); ++i) {
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
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  const char tag_name[] = "VTX_ADJ";
  create_mesh( true, false, false, false, tag_name, var_len );
  int ids[2] = { 7, 10 };
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, ids, 2 );
  CHECK_ERR(rval);
  
  MBTag tag = check_tag( mb, tag_name, MB_TAG_DENSE, MB_TYPE_HANDLE, 
                         var_len ? -1 : 4*sizeof(MBEntityHandle) );
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  CHECK( !verts.empty() );
  for (MBRange::iterator i = verts.begin(); i != verts.end(); ++i) {
    std::vector<MBEntityHandle> adj, val;
    rval = mb.get_adjacencies( &*i, 1, 2, false, adj, MBInterface::UNION );
    CHECK_ERR(rval);
    CHECK(!adj.empty());
    
    int num;
    const void* ptr;
    rval = mb.tag_get_data( tag, &*i, 1, &ptr, &num );
    CHECK_ERR(rval);
    CHECK_EQUAL( 0, num % (int)sizeof(MBEntityHandle) );
    num /= sizeof(MBEntityHandle);
    
    if (var_len) {
      CHECK( num > 0 );
      CHECK( num < 5 );
    }
    else {
      CHECK_EQUAL( 4, num );
    }
    
    val.clear();
    const MBEntityHandle* dat = (const MBEntityHandle*)ptr;
    for (const MBEntityHandle* end = dat+num; dat != end; ++dat)
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
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  create_mesh( false, false, true, false );
  int id = 5;
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, &id, 1 );
  CHECK_ERR(rval);
  
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  int id2 = identify_set( mb, verts );
  CHECK_EQUAL( id, id2 );
  
  int elems;
  rval = mb.get_number_entities_by_type( 0, MBQUAD, elems );
  CHECK_ERR(rval);
  CHECK_EQUAL( QUAD_INT*QUAD_INT/NUM_SETS, elems );
}

void test_read_tagged_nodes()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface& mb = moab;
  
  create_mesh( false, false, false, true );
  int id = 1; // NOTE: this test will only succeed for ID == 1 
  MBEntityHandle file_set;
  rval = mb.load_file( TEST_FILE, file_set, READ_OPTS, ID_TAG_NAME, &id, 1 );
  CHECK_ERR(rval);
  
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  int id2 = identify_set( mb, verts );
  CHECK_EQUAL( id, id2 );
  
  int elems;
  rval = mb.get_number_entities_by_type( 0, MBQUAD, elems );
  CHECK_ERR(rval);
  CHECK_EQUAL( QUAD_INT*QUAD_INT/NUM_SETS, elems );
}


//! Read in the polyhedra contained in a set
void test_read_one_set_polyhedra()
{
  MBErrorCode rval;
  MBCore instance;
  MBInterface& mb = instance;
  
    // create a 2x2x1 block of hexes, splitting each hex face
    // into two triangles to form an 12-sided polyhedron
  MBEntityHandle verts[18], hexes[4];
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
    MBEntityHandle conn[8];
    for (int j = 0; j < 8; ++j)
      conn[j] = verts[hexconn[i][j]];
    rval = mb.create_element( MBHEX, conn, 8, hexes[i] );
    CHECK_ERR(rval);
  }
  
  MBTag tri_tag;
  rval = mb.tag_create( "tris", 2*sizeof(MBEntityHandle), MB_TAG_SPARSE, MB_TYPE_HANDLE, tri_tag, 0 );
  CHECK_ERR(rval);
  
  std::vector<MBEntityHandle> quads;
  MBEntityHandle tris[12], poly[4];
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
      const MBEntityHandle* conn;
      int len;
      rval = mb.get_connectivity( quads[j], conn, len );
      CHECK_ERR(rval);
      CHECK_EQUAL( 4, len );
      MBEntityHandle tri_conn[2][3] = {  { conn[0], conn[1], conn[2] },
                                         { conn[2], conn[3], conn[0] } };
      rval = mb.create_element( MBTRI, tri_conn[0], 3, tris[2*j  ] ); CHECK_ERR(rval);
      rval = mb.create_element( MBTRI, tri_conn[1], 3, tris[2*j+1] ); CHECK_ERR(rval);
      rval = mb.tag_set_data( tri_tag, &quads[j], 1, tris + 2*j ); CHECK_ERR(rval);
    }
    
    rval = mb.create_element( MBPOLYHEDRON, tris, 12, poly[i] );
    CHECK_ERR(rval);
  }
  
  MBRange all_tri;
  rval = mb.get_entities_by_type( 0, MBTRI, all_tri );
  CHECK_ERR(rval);
  CHECK_EQUAL( 40, (int)all_tri.size() );
  
  rval = mb.delete_entities( hexes, 4 ); CHECK_ERR(rval);
  rval = mb.delete_entities( &quads[0], quads.size() ); CHECK_ERR(rval);
  
  MBEntityHandle sets[2];
  rval = mb.create_meshset( 0, sets[0] ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[0], poly, 2 ); CHECK_ERR(rval);
  rval = mb.create_meshset( 0, sets[1] ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], poly+2, 2 ); CHECK_ERR(rval);
  
  MBTag id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
  CHECK_ERR(rval);
  int ids[2] = { 2, 3 };
  rval = mb.tag_set_data( id_tag, sets, 2, ids ); CHECK_ERR(rval);
  
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  rval = mb.delete_mesh();
  CHECK_ERR(rval);
  
  MBEntityHandle file;
  rval = mb.load_file( TEST_FILE, file, READ_OPTS, ID_TAG_NAME, ids, 1 );
  CHECK_ERR(rval);
  
  MBRange rpoly;
  rval = mb.get_entities_by_type( 0, MBPOLYHEDRON, rpoly );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)rpoly.size() );
  
  MBRange polyverts;
  rval = mb.get_adjacencies( rpoly, 0, false, polyverts, MBInterface::UNION );
  CHECK_ERR(rval);
  CHECK_EQUAL( 12, (int)polyverts.size() );
  
  for (MBRange::iterator it = polyverts.begin(); it != polyverts.end(); ++it) {
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
  MBErrorCode rval;
  MBCore instance;
  MBInterface& mb = instance;
  
  MBTag id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
  CHECK_ERR(rval);
  
    // create sets and assign an ID to each
  const int len = 5;
  MBEntityHandle set[2*len+2];
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
    MBEntityHandle contents[len] = { set[cont_ids[i][0] - 1],
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
    
    MBEntityHandle file;
    int id = i+1;
    rval = mb.load_file( TEST_FILE, file, READ_OPTS, ID_TAG_NAME, &id, 1 );
    CHECK_ERR(rval);
    
      // check that the total number of sets read is as expected
    MBRange sets;
    rval = mb.get_entities_by_type( 0, MBENTITYSET, sets );
    CHECK_ERR(rval);
    MBRange::iterator it = sets.find( file );
    if (it != sets.end())
      sets.erase( it );
    CHECK_EQUAL( len+1, (int)sets.size() );
    
      // check that we read in the set specified by ID to the reader
    rval = mb.tag_get_handle( ID_TAG_NAME, id_tag );
    CHECK_ERR(rval);
    sets.clear();
    const void* data[] = { &id };
    rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &id_tag, data, 1, sets );
    CHECK_ERR(rval);
    CHECK_EQUAL( 1, (int)sets.size() );
      
      // check that it contains the expected sets
    MBEntityHandle owner = sets.front();
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

static void check_children( ChildTestMode mode, MBInterface& mb, int id, MBTag id_tag, MBEntityHandle file )
{
  const int exp_num_sets = (mode == CHILD_NONE) ? 1 : id;
  const int exp_num_edges = (mode == CHILD_CONTENTS) ? id : 1;
  
  MBErrorCode rval;
  MBRange range;
  rval = mb.get_entities_by_type( 0, MBEDGE , range );
  CHECK_ERR(rval);
  CHECK_EQUAL( exp_num_edges, (int)range.size() );
  range.clear();
  rval = mb.get_entities_by_type( 0, MBENTITYSET , range );
  CHECK_ERR(rval);
  MBRange::iterator it = range.find( file );
  CHECK( it != range.end() );
  range.erase( it );
  CHECK_EQUAL( exp_num_sets, (int)range.size() );
  
  MBEntityHandle set;
  const void* val[] = {&id};
  range.clear();
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &id_tag, val, 1, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() );
  set = range.front();
  
  if (mode == CHILD_NONE) {
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
    if (mode == CHILD_CONTENTS || i == id) {
      CHECK_EQUAL( 1, (int)range.size() );
      const MBEntityHandle* conn;
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
      
    std::vector<MBEntityHandle> children;
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
void test_read_children_common( ChildTestMode mode )
{
  MBErrorCode rval;
  MBCore instance;
  MBInterface& mb = instance;
  
  MBTag id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
  CHECK_ERR(rval);
  
    // Create a string of edges from [0,INT] along the X axis each 1 unit in length.
    // Create a set for edge edge, containing the edge and the parent of the 
    // set containing the previous (closer to origin) edge.  Assign each set an
    // ID that is the X coordinate of the larger of the two vertices of the edge
    // contained in the set.
  const int INT = 64;
  MBEntityHandle verts[INT+1], edges[INT], sets[INT];
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
      rval = mb.add_child_meshset( sets[i], sets[i-1] );
      CHECK_ERR(rval);
    }
  }
  
    // Write the data
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
 
  MBEntityHandle file;
  std::string opt( READ_OPTS );
  opt += ";CHILDREN=";
  opt += set_read_opts[mode];

  const int test_ids[] = { 2, 7, INT/3-1, INT/2+1, INT-3 };
  const int num_test_ids = sizeof(test_ids)/sizeof(int);
  for (int i = 0; i < num_test_ids; ++i) {
    CHECK (test_ids[i] <= INT);
    
    rval = mb.delete_mesh();
    CHECK_ERR(rval);
    
    rval = mb.load_file( TEST_FILE, file, opt.c_str(), ID_TAG_NAME, test_ids+i, 1 );
    CHECK_ERR(rval);
    rval = mb.tag_get_handle( ID_TAG_NAME, id_tag );
    CHECK_ERR(rval);
    
    check_children( mode, mb, test_ids[i], id_tag, file );
  }
}

static void check_num_verts( MBInterface& mb, MBTag tag, int id, int num_vtx )
{
  MBErrorCode rval;
  const void* val[] = {&id};
  MBRange range;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, val, 1, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() );

  MBEntityHandle set = range.front();
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
  
  MBErrorCode rval;
  MBCore instance;
  MBInterface& mb = instance;
  
    // read some sets
  MBEntityHandle file;
  const int ids[] = { 1, 5, 9 };
  const int num_sets = sizeof(ids)/sizeof(int);
  rval = mb.load_file( TEST_FILE, file, READ_OPTS, ID_TAG_NAME, ids, num_sets );
  CHECK_ERR(rval);

  MBTag id_tag;
  rval = mb.tag_get_handle( ID_TAG_NAME, id_tag );
  CHECK_ERR(rval);
  
    // expect all sets adjacent to the specified sets because
    // they share vertices.
  MBRange verts;
  for (int i = 0; i < num_sets; ++i) {
    if (ids[i] > 1)
      check_num_verts( mb, id_tag, ids[i]-1, QUAD_INT+1 );
    check_num_verts( mb, id_tag, ids[i], (QUAD_INT+1)*(SET_WIDTH+1) );
    if (ids[i] < NUM_SETS)
      check_num_verts( mb, id_tag, ids[i]+1, QUAD_INT+1 );
  }
}

//! Test reading of explicit adjacencies
void test_read_adjacencies()
{
  MBErrorCode rval;
  MBCore instance;
  MBInterface& mb = instance;
  
    // create four hexes sharing an edge
  MBEntityHandle verts[3][3][2], hexes[2][2];
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
      MBEntityHandle conn[] = { verts[i  ][j  ][0],
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
  MBEntityHandle edge_conn[2] = { verts[1][1][0], verts[1][1][1] };
  MBEntityHandle edges[2];
  rval = mb.create_element( MBEDGE, edge_conn, 2, edges[0] ); CHECK_ERR(rval);
  rval = mb.create_element( MBEDGE, edge_conn, 2, edges[1] ); CHECK_ERR(rval);
    // mark one edge as adjacent to the left two hexes and the
    // other as adjacent to the right two
  rval = mb.add_adjacencies( edges[0], hexes[0], 2, true ); CHECK_ERR(rval);
  rval = mb.add_adjacencies( edges[1], hexes[1], 2, true ); CHECK_ERR(rval);
    // create two sets containing the front two and the rear two
    // hexes, respectively.
  MBEntityHandle sets[2];
  rval = mb.create_meshset( MESHSET_SET, sets[0] ); CHECK_ERR(rval);
  rval = mb.create_meshset( MESHSET_SET, sets[1] ); CHECK_ERR(rval);
  MBEntityHandle set1[4] = { hexes[0][0], hexes[1][0], edges[0], edges[1] };
  MBEntityHandle set2[4] = { hexes[0][1], hexes[1][1], edges[0], edges[1] };
  rval = mb.add_entities( sets[0], set1, 4 ); CHECK_ERR(rval);
  rval = mb.add_entities( sets[1], set2, 4 ); CHECK_ERR(rval);
  
    // assign IDs to sets
  MBTag id_tag;
  rval = mb.tag_create( ID_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, id_tag, 0 );
  CHECK_ERR(rval);
  int ids[2] = { 1, 2 };
  rval = mb.tag_set_data( id_tag, sets, 2, ids );
  CHECK_ERR(rval);
  
    // write mesh
  rval = mb.write_file( TEST_FILE, "MOAB" );
  CHECK_ERR(rval);
  
    // read mesh
  rval = mb.delete_mesh(); CHECK_ERR(rval);
  MBEntityHandle file;
  rval = mb.load_file( TEST_FILE, file, READ_OPTS, ID_TAG_NAME, ids, 1 );
  CHECK_ERR(rval);
  
    // expect two hexes and two edges
  MBRange range;
  rval = mb.get_entities_by_type( 0, MBHEX, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)range.size() );
  MBEntityHandle h1 = range.front(), h2 = range.back();
  range.clear();
  rval = mb.get_entities_by_type( 0, MBEDGE, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 2, (int)range.size() );
 
    // expecte each edge to have one of the hexes
  range.clear();
  rval = mb.get_adjacencies( &h1, 1, 1, false, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() ); 
  MBEntityHandle e1 = range.front();
  range.clear();
  rval = mb.get_adjacencies( &h2, 1, 1, false, range );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)range.size() ); 
  MBEntityHandle e2 = range.front();
  
  CHECK( e1 != e2 );
}


