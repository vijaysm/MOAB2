#include "moab/Core.hpp"
#include "TestUtil.hpp"
#include "moab/Range.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include <stdlib.h>
#include <stdio.h>

using namespace moab;

static bool keep_files = false;

void test_var_length_no_data();

void test_var_length_data();

void test_var_length_data_big();

void test_var_length_big_data();

void test_var_length_opaque();

void test_var_length_mesh_data();

void test_var_length_default_data();

void test_var_length_mesh_opaque();

void test_var_length_default_opaque();

void test_var_length_handle_tag();

void create_mesh( Interface& mb );

void create_big_mesh( Interface& mb );

void compare_tags( const char* name, Interface& mb1, Interface& mb2 );

void read_write( const char* filename, Interface& write, Interface& reader );

#define CHECK_ERR_FILE( ERRCODE, FILENAME ) \
do { \
  if (MB_SUCCESS != ERRCODE && !keep_files) \
    remove( FILENAME ); \
  CHECK_ERR(ERRCODE); \
while (false)

int main(int argc, char* argv[])
{
#ifdef USE_MPI
  int fail = MPI_Init(&argc, &argv);
  if (fail) return fail;
#endif

  if (argc != 1) {
    if (argc != 2 || strcmp(argv[1],"-k")) {
      fprintf( stderr, "Usage: %s [-k]\n", argv[0] );
      abort();
    }
    keep_files = true;
  } 

  int err_count = 0;
  err_count += RUN_TEST( test_var_length_no_data );
  err_count += RUN_TEST( test_var_length_data );
  err_count += RUN_TEST( test_var_length_big_data );
  err_count += RUN_TEST( test_var_length_opaque );
  err_count += RUN_TEST( test_var_length_mesh_data );
  err_count += RUN_TEST( test_var_length_default_data );
  err_count += RUN_TEST( test_var_length_mesh_opaque );
  err_count += RUN_TEST( test_var_length_default_opaque );
  err_count += RUN_TEST( test_var_length_handle_tag );
  err_count += RUN_TEST( test_var_length_data_big );

#ifdef USE_MPI
  fail = MPI_Finalize();
  if (fail) return fail;
#endif

  return err_count;
}

void test_var_length_no_data()
{
  ErrorCode rval;
  Core moab1, moab2;
  Interface &mb1 = moab1, &mb2 = moab2;
  Tag tag;
  
  create_mesh( mb1 );
  rval = mb1.tag_get_handle( "test_tag", 0, MB_TYPE_DOUBLE, tag, MB_TAG_EXCL|MB_TAG_VARLEN|MB_TAG_DENSE );
  CHECK_ERR( rval );
  
  read_write( "test_var_length_no_data.h5m", mb1, mb2 );
  compare_tags( "test_tag", mb1, mb2 );
}

void test_var_length_data_common( const char* filename, Interface& mb1, bool opaque = false )
{
    // create tag
  ErrorCode rval;
  Tag tag;
  DataType type = opaque ? MB_TYPE_OPAQUE : MB_TYPE_INTEGER;
  rval = mb1.tag_get_handle( "test_tag", 0, type, tag, MB_TAG_EXCL|MB_TAG_VARLEN|MB_TAG_SPARSE );
  CHECK_ERR( rval );
  
    // get all entities
  Range entities;
  rval = mb1.get_entities_by_handle( 0, entities );
  CHECK_ERR(rval);
  
    // Set tag data.
    // Tag data will be list of integer data as follows:
    //   number of values (counting this value)
    //   step, 2*step, 3*step, ...
  for (Range::const_iterator i = entities.begin(); i != entities.end(); ++i) {
    EntityHandle h = *i;
      // generate some data to write
    int num_values = h % 6 + 1;
    EntityType etype = mb1.type_from_handle(h);
    int step = (h%2) ? 1+(int)etype : -1-(int)etype;
    std::vector<int> tag_data( num_values, num_values );
    for (int j = 1; j < num_values; ++j)
      tag_data[j] = j*step;
      // set tag data
    const void* ptrarr[]= { &tag_data[0] };
    if (opaque) num_values *= sizeof(int);
    rval = mb1.tag_set_by_ptr(tag, &h,1, ptrarr, &num_values );
    CHECK_ERR(rval);
  }
  
    // write and read tag data
  Core moab;
  Interface &mb2 = moab;
  read_write( filename, mb1, mb2 );
  compare_tags( "test_tag", mb1, mb2 );
  
    // get new tag handle
  tag = 0;
  rval = mb2.tag_get_handle( "test_tag", 0, type, tag );
  CHECK_ERR(rval);
  
    // check consistency of tag values
  entities.clear();
  mb2.get_entities_by_handle( 0, entities );
    // remove sets created during read/write process
  Range sets;
  mb2.get_entities_by_type( 0, MBENTITYSET, sets );
  entities = subtract( entities, sets);
  for (Range::const_iterator i = entities.begin(); i != entities.end(); ++i) {
      // get data
    const void* ptrarr[] = { NULL };
    int size;
    rval = mb2.tag_get_by_ptr( tag, &*i, 1, ptrarr, &size );
    CHECK_ERR(rval);
    if (opaque) {
      CHECK_EQUAL( (size_t)0, size % sizeof(int) );
      size /= sizeof(int);
    }
    const int* dataptr = reinterpret_cast<const int*>(ptrarr[0]);
    CHECK( NULL != dataptr );
      // check size 
    CHECK( size > 0 );
    CHECK_EQUAL( size, dataptr[0] );
      // check other values
    if (size > 2) {
      int step = dataptr[1];
      for (int j = 2; j < size; ++j)
        CHECK_EQUAL( j*step, dataptr[j] );
    }
  }
}

void test_var_length_data()
{
  Core moab;
  create_mesh( moab );
  test_var_length_data_common( "test_var_length_data.h5m", moab );
}

void test_var_length_data_big()
{
  Core moab;
  create_big_mesh( moab );
  test_var_length_data_common( "test_var_length_data_big.h5m", moab );
}


void calculate_big_value( Interface& moab, EntityHandle vert, size_t size, double* data )
{
    // Make values like Fibonacci numbers, except use X and Y coords
    // rather than 0 and 1 as first two values.

  CHECK( size >= 3 );
  ErrorCode rval = moab.get_coords( &vert, 1, data );
  CHECK_ERR(rval);

  for (size_t j = 2; j < size; ++j)
    data[j] = data[j-2] + data[j-1];
  CHECK_ERR(rval);
}


void test_var_length_big_data()
{
  ErrorCode rval;
  Core moab1, moab2;
  Interface &mb1 = moab1, &mb2 = moab2;
  Tag tag;
  
  create_mesh( mb1 );
  rval = mb1.tag_get_handle( "test_tag", 0, MB_TYPE_DOUBLE, tag, MB_TAG_SPARSE|MB_TAG_VARLEN|MB_TAG_EXCL );
  CHECK_ERR( rval );
  
    // choose 3 vertices upon which to set data
  Range range;
  rval = mb1.get_entities_by_type( 0, MBVERTEX, range );
  CHECK_ERR(rval);
  EntityHandle verts[3] = { range.front(), 
                              *(range.begin() += range.size()/3), 
                              *(range.begin() += 2*range.size()/3) };
  
    // set 1-millon value tag data on three vertices
  std::vector<double> data(1000000);
  for (int i = 0; i < 3; ++i) {
    calculate_big_value( mb1, verts[i], data.size(), &data[0] );
    const void* ptr = &data[0];
    const int size = data.size();
    rval = mb1.tag_set_by_ptr( tag, verts + i, 1, &ptr, &size );
    CHECK_ERR(rval);
  }
  
  read_write( "test_var_length_big_data.h5m", mb1, mb2 );
  compare_tags( "test_tag", mb1, mb2 );
  
    // check 3 tagged vertices
  rval = mb2.tag_get_handle( "test_tag", 0, MB_TYPE_DOUBLE, tag );
  CHECK_ERR(rval);
  range.clear();
  rval = mb2.get_entities_by_type_and_tag( 0, MBVERTEX, &tag, 0, 1, range, Interface::UNION );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)3, range.size() );
  
    // check tag values
  for (Range::const_iterator i = range.begin(); i != range.end(); ++i) {
      // calculate expected value
    const EntityHandle h = *i;
    calculate_big_value( mb2, h, data.size(), &data[0] );
    
      // get actual value
    const void* ptr;
    int size;
    rval = mb2.tag_get_by_ptr( tag, &h, 1, &ptr, &size );
    CHECK_ERR(rval);
    CHECK_EQUAL( data.size(), (size_t)size );
    
      // compare values
    const double* act_data = reinterpret_cast<const double*>(ptr);
    int wrong_count = 0;
    for (size_t j = 0; j < data.size(); ++j)
      if (act_data[j] != data[j])
        ++wrong_count;
    CHECK_EQUAL( 0, wrong_count );
  }    
}  
  


void test_var_length_opaque()
{
  Core moab;
  create_mesh( moab );
  test_var_length_data_common( "test_var_length_opaque.h5m", moab, true );  
}

void test_global_value_common( bool mesh_value )
{
  Core moab;
  Interface &mb = moab;
  create_mesh( mb );
  
  // get three vertices
  Range vertices;
  ErrorCode rval = mb.get_entities_by_type( 0, MBVERTEX, vertices );
  CHECK_ERR(rval);
  CHECK( vertices.size() >= 3 );
  EntityHandle handles[3];
  Range::const_iterator i = vertices.begin();
  handles[0] = *i; ++i;
  handles[1] = *i; ++i;
  handles[2] = *i; ++i;

  // get vertex coordinates
  double coords[9];
  rval = mb.get_coords( handles, 3, coords );
  CHECK_ERR( rval );
  
  // create tag to hold vertex data
  Tag handle_tag = 0;
  void* default_val = mesh_value ? 0 : handles;
  int default_val_size = mesh_value ? 0 : 3;
  rval = mb.tag_get_handle( "handle_tag", default_val_size, MB_TYPE_HANDLE, 
                            handle_tag, MB_TAG_DENSE|MB_TAG_VARLEN|MB_TAG_EXCL,
                            default_val );
  CHECK_ERR( rval );
  
  // create tag to hold vertex coordinates
  Tag coord_tag = 0;
  default_val = mesh_value ? 0 : coords;
  default_val_size = mesh_value ? 0 : 9;
  rval = mb.tag_get_handle( "coord_tag", default_val_size, MB_TYPE_DOUBLE, 
                             coord_tag, MB_TAG_SPARSE|MB_TAG_VARLEN|MB_TAG_EXCL,
                            default_val );
  CHECK_ERR( rval );
  
  
  // if doing mesh tag, set it
  if (mesh_value) {
    int size = 3;
    const void* ptrarr[] = { handles };
    const EntityHandle root = 0;
    rval = mb.tag_set_by_ptr( handle_tag, &root, 1, ptrarr, &size );
    CHECK_ERR(rval);
    
    size = 9;
    ptrarr[0] = coords;
    rval = mb.tag_set_by_ptr( coord_tag, &root, 1, ptrarr, &size );
    CHECK_ERR(rval);
  }
  
    // write and read file
  Core moab2;
  Interface &mb2 = moab2;
  read_write( mesh_value ? "test_var_length_mesh_data.h5m" : "test_var_length_default_data.h5m", mb, mb2 );
  compare_tags( "handle_tag", mb, mb2 );
  compare_tags( "coord_tag", mb, mb2 );
  
    // get tag handles
  handle_tag = coord_tag = 0;
  rval = mb2.tag_get_handle( "handle_tag", 0, MB_TYPE_HANDLE, handle_tag );
  CHECK_ERR(rval);
  rval = mb2.tag_get_handle( "coord_tag", 0, MB_TYPE_DOUBLE, coord_tag );
  CHECK_ERR(rval);
  
    // get tag data
  int handle_tag_size = 0, coord_tag_size = 0;
  const void* ptrs[2];
  if (mesh_value) {
    const EntityHandle root = 0;
    rval = mb2.tag_get_by_ptr( handle_tag, &root, 1, ptrs, &handle_tag_size );
    CHECK_ERR(rval);
    rval = mb2.tag_get_by_ptr( coord_tag, &root, 1, ptrs+1, &coord_tag_size );
    CHECK_ERR(rval);
  }
  else {
    rval = mb2.tag_get_default_value( handle_tag, ptrs[0], handle_tag_size );
    CHECK_ERR(rval);
    rval = mb2.tag_get_default_value( coord_tag, ptrs[1], coord_tag_size );
    CHECK_ERR(rval);
  }
  
    // check expected sizes
  CHECK_EQUAL( 3, handle_tag_size );
  CHECK_EQUAL( 9, coord_tag_size );
  CHECK( ptrs[0] != NULL );
  CHECK( ptrs[1] != NULL );
  
    // check valid handles
  const EntityHandle* handle_vals = reinterpret_cast<const EntityHandle*>(ptrs[0]);
  CHECK( handle_vals[0] != 0 );
  CHECK( handle_vals[1] != 0 );
  CHECK( handle_vals[2] != 0 );
  CHECK_EQUAL( MBVERTEX, mb2.type_from_handle( handle_vals[0] ) );
  CHECK_EQUAL( MBVERTEX, mb2.type_from_handle( handle_vals[1] ) );
  CHECK_EQUAL( MBVERTEX, mb2.type_from_handle( handle_vals[2] ) );
  
    // check correct coordinate values
  const double* coord_vals = reinterpret_cast<const double*>(ptrs[1]);
  rval = mb2.get_coords( handle_vals, 3, coords );
  CHECK_ERR( rval );
  CHECK_REAL_EQUAL( coords[0], coord_vals[0], 1e-12 );
  CHECK_REAL_EQUAL( coords[1], coord_vals[1], 1e-12 );
  CHECK_REAL_EQUAL( coords[2], coord_vals[2], 1e-12 );
  CHECK_REAL_EQUAL( coords[3], coord_vals[3], 1e-12 );
  CHECK_REAL_EQUAL( coords[4], coord_vals[4], 1e-12 );
  CHECK_REAL_EQUAL( coords[5], coord_vals[5], 1e-12 );
  CHECK_REAL_EQUAL( coords[6], coord_vals[6], 1e-12 );
  CHECK_REAL_EQUAL( coords[7], coord_vals[7], 1e-12 );
  CHECK_REAL_EQUAL( coords[8], coord_vals[8], 1e-12 );
}

void test_var_length_mesh_data()
{
  test_global_value_common( true );
}

void test_var_length_default_data()
{
  test_global_value_common( false );
}

void test_global_opaque_common( bool mesh_value )
{
  ErrorCode rval;
  Core moab;
  Interface &mb = moab;
  create_mesh( mb );
  
  const char data[] = { 'J', 'A', 'S', 'O', 'N' };
  const int datalen = sizeof(data);
  CHECK_EQUAL( 5, datalen );
  
  // create tag 
  Tag tag = 0;
  const void* default_val = mesh_value ? 0 : data;
  int default_val_size = mesh_value ? 0 : datalen;
  rval = mb.tag_get_handle( "opaque_tag", default_val_size, MB_TYPE_OPAQUE, tag, 
                            MB_TAG_DENSE|MB_TAG_VARLEN|MB_TAG_EXCL, default_val );
  CHECK_ERR( rval );
  
  // if doing mesh tag, set it
  if (mesh_value) {
    const void* ptrarr[] = { data };
    const EntityHandle root = 0;
    rval = mb.tag_set_by_ptr( tag, &root, 1, ptrarr, &datalen );
    CHECK_ERR(rval);
  }
  
    // write and read file
  Core moab2;
  Interface &mb2 = moab2;
  read_write( mesh_value ? "test_var_length_mesh_opaque.h5m" : "test_var_length_default_opaque.h5m", mb, mb2 );
  compare_tags( "opaque_tag", mb, mb2 );
  
    // get tag handles
  tag = 0;
  rval = mb2.tag_get_handle( "opaque_tag", 0, MB_TYPE_OPAQUE, tag );
  CHECK_ERR(rval);
  
    // get tag data
  int tag_size = 0;
  const void* ptrs[1];
  if (mesh_value) {
    const EntityHandle root = 0;
    rval = mb2.tag_get_by_ptr( tag, &root, 1, ptrs, &tag_size );
    CHECK_ERR(rval);
  }
  else {
    rval = mb2.tag_get_default_value( tag, ptrs[0], tag_size );
    CHECK_ERR(rval);
  }
  
    // check size
  CHECK_EQUAL( datalen, tag_size );
  CHECK( ptrs[0] != NULL );
  
    // check values
  const char* tag_data = reinterpret_cast<const char*>(ptrs[0]);
  for (int i = 0; i < datalen; ++i)
    CHECK_EQUAL( data[i], tag_data[i] );
}

void test_var_length_mesh_opaque()
{
  test_global_opaque_common( true );
}

void test_var_length_default_opaque()
{
  test_global_opaque_common( false );
}

  
void test_var_length_handle_tag()
{
  ErrorCode rval;
  Core moab1, moab2;
  Interface &mb1 = moab1, &mb2 = moab2;
  Tag tag;
  Range::const_iterator  i;
  
  create_mesh( mb1 );
  rval = mb1.tag_get_handle( "test_tag", 0, MB_TYPE_HANDLE, tag, 
                             MB_TAG_SPARSE|MB_TAG_VARLEN|MB_TAG_EXCL );
  CHECK_ERR( rval );
  
    // Get all entities
  Range range;
  rval = mb1.get_entities_by_handle( 0, range );
  CHECK_ERR(rval);

    // For each entity, if it is a vertex store its own handle
    // in its tag.  Otherwise store the element connectivity list
    // in the tag.  Skip entity sets.
  size_t num_tagged_entities = 0;
  for (i = range.begin(); i != range.end(); ++i) {
    EntityHandle h = *i;
    EntityType type = mb1.type_from_handle( h );
    if (type == MBVERTEX) {
      const int size = 1;
      const void* ptr = &h;
      rval = mb1.tag_set_by_ptr( tag, &h, 1, &ptr, &size );
      CHECK_ERR(rval);
      ++num_tagged_entities;
    }
    else if (type != MBENTITYSET) {
      int size = 0;
      const EntityHandle* conn = 0;
      rval = mb1.get_connectivity( h, conn, size );
      CHECK_ERR(rval);
      const void* ptr = conn;
      rval = mb1.tag_set_by_ptr( tag, &h, 1, &ptr, &size );
      CHECK_ERR(rval);
      ++num_tagged_entities;
   }
  }
  
  read_write( "test_var_length_handle_tag.h5m", mb1, mb2 );
  compare_tags( "test_tag", mb1, mb2 );
  
    // check number of tagged entities
  rval = mb2.tag_get_handle( "test_tag", 0, MB_TYPE_HANDLE, tag );
  CHECK_ERR(rval);
  range.clear();
  for (EntityType t = MBVERTEX; t != MBENTITYSET; ++t) {
    rval = mb2.get_entities_by_type_and_tag( 0, t, &tag, 0, 1, range, Interface::UNION );
    CHECK_ERR(rval);
  }
  CHECK_EQUAL( num_tagged_entities, range.size() );
  
    // check tag values
  for (i = range.begin(); i != range.end(); ++i) {
    EntityHandle h = *i;
    
    const void* ptr;
    int size;
    rval = mb2.tag_get_by_ptr( tag, &h, 1, &ptr, &size );
    CHECK_ERR(rval);
    
    const EntityHandle* handles = reinterpret_cast<const EntityHandle*>(ptr);
    
    if (mb2.type_from_handle(h) == MBVERTEX) {
      CHECK_EQUAL( 1, size );
      CHECK_EQUAL( h, *handles );
    }
    else {
      int len;
      const EntityHandle* conn;
      rval = mb2.get_connectivity( h, conn, len );
      CHECK_ERR(rval);
      CHECK_EQUAL( len, size );
      for (int j = 0; j < len; ++j)
        CHECK_EQUAL( conn[j], handles[j] );
    }
  }
}

void create_structured_quad_mesh( Interface& mb, int x, int y )
{
  ErrorCode rval;
  
  const double z = 2.1;
  std::vector<EntityHandle> verts((x+1)*(y+1));
  for (int i = 0; i <= x; ++i) {
    for (int j = 0; j<= y; ++j) {
      double coords[3] = { i, j, z };
      rval = mb.create_vertex( coords, verts[i + (x+1)*j] );
      CHECK_ERR( rval );
    }
  }
  
  std::vector<EntityHandle> elems(x*y);
  for (int i = 0; i < x; ++i) {
    for (int j = 0; j < y; ++j ) {
      EntityHandle conn[4] = { verts[i +     (x+1)* j   ],
                                 verts[i + 1 + (x+1)* j   ],
                                 verts[i + 1 + (x+1)*(j+1)],
                                 verts[i +     (x+1)*(j+1)] };
      rval = mb.create_element( MBQUAD, conn, 4, elems[i + x*j] );
      CHECK_ERR( rval );
    }
  }
}

void create_mesh( Interface& mb )
{
  create_structured_quad_mesh( mb, 2, 2 );
}

void create_big_mesh( Interface& mb )
{
  create_structured_quad_mesh( mb, 300, 300 );
}

void compare_tags( const char* name, Interface& mb1, Interface& mb2 )
{
  ErrorCode rval;
  Tag tag1, tag2;
  rval = mb1.tag_get_handle( name, 0, MB_TYPE_OPAQUE, tag1, MB_TAG_ANY );
  CHECK_ERR(rval);
  rval = mb2.tag_get_handle( name, 0, MB_TYPE_OPAQUE, tag2, MB_TAG_ANY );
  CHECK_ERR(rval);
  
  int size;
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, mb1.tag_get_length( tag1, size ) );
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, mb2.tag_get_length( tag2, size ) );
  
  TagType storage1, storage2;
  rval = mb1.tag_get_type( tag1, storage1 );
  CHECK_ERR(rval);
  rval = mb2.tag_get_type( tag2, storage2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( storage1, storage2 );
  
  DataType type1, type2;
  rval = mb1.tag_get_data_type( tag1, type1 );
  CHECK_ERR(rval);
  rval = mb2.tag_get_data_type( tag2, type2 );
  CHECK_ERR(rval);
  
  const void *defval1, *defval2;
  int defsize1, defsize2;
  ErrorCode rval1 = mb1.tag_get_default_value( tag1, defval1, defsize1 );
  ErrorCode rval2 = mb2.tag_get_default_value( tag2, defval2, defsize2 );
  if (MB_SUCCESS == rval1) {
    CHECK_ERR(rval2);
 
    CHECK_EQUAL( defsize1, defsize2 );
    CHECK( !memcmp( defval1, defval2, defsize1 ) );
  }
  else if (MB_ENTITY_NOT_FOUND == rval1 || MB_TAG_NOT_FOUND == rval1) 
    CHECK_EQUAL( rval1, rval2 );
  else 
    CHECK_ERR( rval1 );
  
}

void read_write( const char* filename, Interface& writer, Interface& reader )
{
  ErrorCode rval = writer.write_mesh( filename );
  if (!keep_files && MB_SUCCESS != rval)
    remove( filename );
  CHECK_ERR(rval);
  rval = reader.load_mesh( filename );
  if (!keep_files)
    remove( filename );
  CHECK_ERR(rval);
}
