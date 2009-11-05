#include "MBCore.hpp"
#include "MBRange.hpp"
#include "TestUtil.hpp"
#include <stdlib.h>
#include <algorithm>

void test_create_tag();
void test_invalid_tag_size();
void test_get_set_sparse_int();
void test_get_set_dense_int();
void test_get_set_dense_double();
void test_get_set_bit();
void test_get_by_tag( );
void test_get_by_tag_value( );
void test_mesh_value();
void test_get_pointers_sparse();
void test_get_pointers_dense();
void test_set_pointers_sparse();
void test_set_pointers_dense();
void test_get_entity_tags();
void test_delete_sparse_tag();
void test_delete_dense_tag();
void test_delete_mesh_tag();
void test_delete_sparse_data();
void test_delete_dense_data();
void test_delete_bit_data();
void test_create_variable_length_tag();
void test_get_set_variable_length_sparse();
void test_get_set_variable_length_dense();
void test_get_set_variable_length_mesh();
void test_get_ents_with_default_value();

void regression_one_entity_by_var_tag();
void regression_tag_on_nonexistent_entity();

int main()
{
  int failures = 0;
  
  failures += RUN_TEST( test_create_tag );
  failures += RUN_TEST( test_invalid_tag_size );
  failures += RUN_TEST( test_get_set_sparse_int );
  failures += RUN_TEST( test_get_set_dense_int );
  failures += RUN_TEST( test_get_set_dense_double );
  failures += RUN_TEST( test_get_set_bit );
  failures += RUN_TEST( test_get_by_tag );
  failures += RUN_TEST( test_get_by_tag_value );
  failures += RUN_TEST( test_mesh_value );
  failures += RUN_TEST( test_get_pointers_sparse );
  failures += RUN_TEST( test_get_pointers_dense );
  failures += RUN_TEST( test_set_pointers_sparse );
  failures += RUN_TEST( test_set_pointers_dense );
  failures += RUN_TEST( test_get_entity_tags );
  failures += RUN_TEST( test_delete_sparse_tag );
  failures += RUN_TEST( test_delete_dense_tag );
  failures += RUN_TEST( test_delete_mesh_tag );
  failures += RUN_TEST( test_delete_sparse_data );
  failures += RUN_TEST( test_delete_dense_data );
  failures += RUN_TEST( test_delete_bit_data );
  failures += RUN_TEST( test_create_variable_length_tag );
  failures += RUN_TEST( test_get_set_variable_length_sparse );
  failures += RUN_TEST( test_get_set_variable_length_dense );
  failures += RUN_TEST( test_get_set_variable_length_mesh );  
  failures += RUN_TEST( test_get_ents_with_default_value );  
  failures += RUN_TEST( regression_one_entity_by_var_tag );
  failures += RUN_TEST( regression_tag_on_nonexistent_entity );
  
  return failures;
}
  

void setup_mesh( MBInterface& mesh );

MBTag test_create_tag( MBInterface& mb,
                       const char* name,
                       unsigned bytes,
                       MBTagType storage,
                       MBDataType type,
                       const void* defval,
                       MBErrorCode expect = MB_SUCCESS );

MBTag test_create_var_len_tag( MBInterface& mb,
                               const char* name,
                               MBTagType storage,
                               MBDataType type,
                               const void* defval,
                               int defval_size,
                               MBErrorCode expect = MB_SUCCESS );

void test_get_set( const char* name,
                   unsigned bytes,
                   MBTagType storage, 
                   MBDataType type, 
                   const void* some_values, 
                   int num_values,
                   const void* default_value,
                   bool set_by_pointer = false,
                   bool get_by_pointer = false );

void test_get_set_variable_length( const char* name,
                                   MBTagType storage, 
                                   MBDataType type, 
                                   const void** values,
                                   const int* lengths,
                                   int num_values,
                                   const void* default_value,
                                   int default_value_length );

void test_mesh_value( MBInterface& mb,
                      const char* tag_name,
                      unsigned tag_size,
                      MBTagType tag_storage,
                      MBDataType tag_type,
                      const void* value );

MBTag test_create_tag( MBInterface& mb,
                       const char* name,
                       unsigned bytes,
                       MBTagType storage,
                       MBDataType type,
                       const void* defval,
                       MBErrorCode expect  )
{
  MBErrorCode rval;
  MBTag tag;
  
  if (type == MB_TYPE_OPAQUE) 
    rval = mb.tag_create( name, bytes, storage, tag, defval );
  else
    rval = mb.tag_create( name, bytes, storage, type, tag, defval );
    
  if (expect != MB_SUCCESS) {
    CHECK_EQUAL( expect, rval );
    return 0;
  }
  CHECK_ERR(rval);
  
  std::string n;
  rval = mb.tag_get_name( tag, n );
  CHECK_ERR(rval);
  CHECK( n == name );
  
  MBTag tag2;
  rval = mb.tag_get_handle( name, tag2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( tag, tag2 );
  
  int s;
  rval = mb.tag_get_size( tag, s );
  CHECK_ERR(rval);
  CHECK_EQUAL( (int)bytes, s );
  
  MBTagType t;
  rval = mb.tag_get_type( tag, t );
  CHECK_ERR(rval);
  CHECK_EQUAL( storage, t );
  
  MBDataType d;
  rval = mb.tag_get_data_type( tag, d );
  CHECK_ERR(rval);
  CHECK_EQUAL( type, d );
  
  std::vector<unsigned char> defv( bytes );
  rval = mb.tag_get_default_value( tag, &defv[0] );
  if (defval) {
    CHECK_ERR(rval);
    unsigned real_bytes = (storage == MB_TAG_BIT) ? 1 : bytes;
    CHECK(!memcmp( defval, &defv[0], real_bytes ));
  }
  else {
    CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  }
  
    // make sure we can't create a second tag w/ the same name
  rval = mb.tag_create( name, bytes, storage, type, tag2, defval );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, rval );
    // we should get back the handle of the existing tag
  CHECK_EQUAL( tag, tag2 );
  
  return tag;
}
  

MBTag test_create_var_len_tag( MBInterface& mb,
                               const char* name,
                               MBTagType storage,
                               MBDataType type,
                               const void* defval,
                               int defval_size,
                               MBErrorCode expect  )
{
  MBErrorCode rval;
  MBTag tag;
  
  rval = mb.tag_create_variable_length( name, storage, type, tag, defval, defval_size );
  if (expect != MB_SUCCESS) {
    CHECK_EQUAL( expect, rval );
    return 0;
  }
  CHECK_ERR(rval);
  
  std::string n;
  rval = mb.tag_get_name( tag, n );
  CHECK_ERR(rval);
  CHECK( n == name );
  
  MBTag tag2;
  rval = mb.tag_get_handle( name, tag2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( tag, tag2 );
  
  int s;
  rval = mb.tag_get_size( tag, s );
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, rval );
  //CHECK_ERR(rval);
  //CHECK_EQUAL( MB_VARIABLE_LENGTH, s );
  
  MBTagType t;
  rval = mb.tag_get_type( tag, t );
  CHECK_ERR(rval);
  CHECK_EQUAL( storage, t );
  
  MBDataType d;
  rval = mb.tag_get_data_type( tag, d );
  CHECK_ERR(rval);
  CHECK_EQUAL( type, d );
  
  int size;
  const void* defv;
  rval = mb.tag_get_default_value( tag, defv, size );
  if (defval) {
    CHECK_ERR(rval);
    CHECK_EQUAL( defval_size, size );
    CHECK(!memcmp( defval, defv, size ));
  }
  else {
    CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  }
  
    // make sure we can't create a second tag w/ the same name
  rval = mb.tag_create_variable_length( name, storage, type, tag2 );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, rval );
    // we should get back the handle of the existing tag
  CHECK_EQUAL( tag, tag2 );
  
  return tag;
}
  

void test_create_tag()
{ 
  MBCore mb;
  unsigned char defval[] = { 1, 2, 3, 4, 0, 0, 0, 0 };
   
    // opaque tags
  test_create_tag( mb, "opaque_tag_sparse", 8, MB_TAG_SPARSE, MB_TYPE_OPAQUE, defval );
  test_create_tag( mb, "opaque_tag_dense" , 4, MB_TAG_DENSE,  MB_TYPE_OPAQUE, defval );
  test_create_tag( mb, "opaque_tag_dense2", 1, MB_TAG_DENSE,  MB_TYPE_OPAQUE, 0 );
  
    // integer tags
  test_create_tag( mb, "int_tag_sparse", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, defval );
  test_create_tag( mb, "int_tag_dense", 8*sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, 0 );
  
    // double tags
  double defval2[] = { 3.14159, 2.71828 };
  test_create_tag( mb, "dbl_tag_sparse", sizeof(double), MB_TAG_SPARSE, MB_TYPE_DOUBLE, 0 );
  test_create_tag( mb, "dbl_tag_dense", 2*sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, defval2 );
    
    // handle tags
  test_create_tag( mb, "h_tag_dense", sizeof(MBEntityHandle), MB_TAG_DENSE, MB_TYPE_HANDLE, defval );
  
    // bit tags
  unsigned char def_bit_val = 0xBF;
  test_create_tag( mb, "bit_tag_1", 2, MB_TAG_BIT, MB_TYPE_BIT, &def_bit_val );
}


void test_create_variable_length_tag()
{ 
  MBCore mb;
  unsigned char defval[] = { 1, 2, 3, 4, 0, 0, 0, 0 };
   
    // opaque tags
  const void* defptr = defval;
  test_create_var_len_tag( mb, "opaque_tag_sparse", MB_TAG_SPARSE, MB_TYPE_OPAQUE, defptr, 8 );
  test_create_var_len_tag( mb, "opaque_tag_dense" , MB_TAG_DENSE,  MB_TYPE_OPAQUE, defptr, 4 );
  test_create_var_len_tag( mb, "opaque_tag_dense2", MB_TAG_DENSE,  MB_TYPE_OPAQUE, 0, 0 );
  
    // integer tags
  test_create_var_len_tag( mb, "int_tag_sparse", MB_TAG_SPARSE, MB_TYPE_INTEGER, defptr, sizeof(int) );
  test_create_var_len_tag( mb, "int_tag_dense", MB_TAG_DENSE, MB_TYPE_INTEGER, 0, 0 );
  
    // double tags
  double defval2[] = { 3.14159, 2.71828 };
  defptr = defval2;
  test_create_var_len_tag( mb, "dbl_tag_sparse", MB_TAG_SPARSE, MB_TYPE_DOUBLE, 0, 0 );
  test_create_var_len_tag( mb, "dbl_tag_dense", MB_TAG_DENSE, MB_TYPE_DOUBLE, defptr, 2*sizeof(double) );
    
    // handle tags
  test_create_var_len_tag( mb, "h_tag_dense", MB_TAG_DENSE, MB_TYPE_HANDLE, 0, 0 );
}


  // test that MOAB enforces the rule that the the size must be multiple of the type
void test_invalid_tag_size()
{
  const MBErrorCode err = MB_INVALID_SIZE;
  MBCore mb;

    // double
  test_create_tag( mb, "std ",   0,                MB_TAG_DENSE , MB_TYPE_DOUBLE, 0, err );
  test_create_tag( mb, "std0",   1,                MB_TAG_SPARSE, MB_TYPE_DOUBLE, 0, err );
  test_create_tag( mb, "std1",   sizeof(double)-1, MB_TAG_DENSE , MB_TYPE_DOUBLE, 0, err );
  test_create_tag( mb, "std2",   sizeof(double)-2, MB_TAG_SPARSE, MB_TYPE_DOUBLE, 0, err );
  test_create_tag( mb, "std3",   sizeof(double)/2, MB_TAG_SPARSE, MB_TYPE_DOUBLE, 0, err );
  test_create_tag( mb, "std4", 2*sizeof(double)-1, MB_TAG_DENSE,  MB_TYPE_DOUBLE, 0, err );
  
    // integer
  test_create_tag( mb, "sti ",   0,             MB_TAG_DENSE , MB_TYPE_INTEGER, 0, err );
  test_create_tag( mb, "sti0",   1,             MB_TAG_SPARSE, MB_TYPE_INTEGER, 0, err );
  test_create_tag( mb, "sti1",   sizeof(int)-1, MB_TAG_DENSE , MB_TYPE_INTEGER, 0, err );
  test_create_tag( mb, "sti2",   sizeof(int)/2, MB_TAG_SPARSE, MB_TYPE_INTEGER, 0, err );
  test_create_tag( mb, "sti3", 2*sizeof(int)-1, MB_TAG_DENSE,  MB_TYPE_INTEGER, 0, err );
    
    // handle
  test_create_tag( mb, "sth ",   0,                        MB_TAG_DENSE , MB_TYPE_HANDLE, 0, err );
  test_create_tag( mb, "sth0",   1,                        MB_TAG_SPARSE, MB_TYPE_HANDLE, 0, err );
  test_create_tag( mb, "sth1",   sizeof(MBEntityHandle)-1, MB_TAG_DENSE , MB_TYPE_HANDLE, 0, err );
  test_create_tag( mb, "sth2",   sizeof(MBEntityHandle)-2, MB_TAG_SPARSE, MB_TYPE_HANDLE, 0, err );
  test_create_tag( mb, "sth3",   sizeof(MBEntityHandle)/2, MB_TAG_SPARSE, MB_TYPE_HANDLE, 0, err );
  test_create_tag( mb, "sth4", 2*sizeof(MBEntityHandle)-1, MB_TAG_DENSE,  MB_TYPE_HANDLE, 0, err );
}

// Given a list of sequential values in memory (pointed to by 'concat'), 
// populate a list of pointers to to each value.  The number of values
// is assumed to be the size of 'list'.  Individual values have a size of 'bytes'
static void concat_to_list( const void* concat, std::vector<const void*>& list, size_t bytes )
{
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(concat);
  std::vector<const void*>::iterator iter = list.begin();
  for (; iter != list.end(); ++iter, ptr += bytes)
    *iter = ptr;
}

  // test get/set of tag values
void test_get_set( const char* name,
                   unsigned bytes,
                   MBTagType storage, 
                   MBDataType type, 
                   const void* some_values, 
                   int num_values,
                   const void* default_value,
                   bool set_by_pointer,
                   bool get_by_pointer )
{
  std::vector<unsigned char> data;
  
    // create mesh and tag
  MBCore moab;
  MBInterface& mb = moab;
  setup_mesh( mb );
  MBTag tag = test_create_tag( mb, name, bytes, storage, type, default_value );
  
    // get some handles to work with
  MBRange entities;
  MBErrorCode rval = mb.get_entities_by_handle( 0, entities );
  CHECK_ERR(rval);
  CHECK( !entities.empty() );
  
    // split handles into four groups
    // a) a single handle
    // b) some non-consecutive handles in an array
    // c) some handles in an MBRange
    // d) remaining handles (remaining in 'entities');
  MBEntityHandle one_handle;
  MBRange::iterator i = entities.begin() += entities.size()/2;
  one_handle = *i;
  entities.erase( i );
  
  MBRange handle_range;
  std::vector<MBEntityHandle> handle_list;
  for (MBRange::const_pair_iterator i =  entities.const_pair_begin(); 
                                    i != entities.const_pair_end(); ++i) {
    if (i->first == i->second || i->second - i->first == 1) {
      MBEntityHandle h1 = i->first, h2 = i->second;
      ++i;
      handle_range.insert( h1, h2 );
    }
    else {
      MBEntityHandle mid = (MBEntityHandle)(i->first + (i->second - i->first + 1) / 2);
      handle_list.push_back( mid );
      handle_range.insert( mid+1, i->second );
    }
  }
  entities = subtract( entities,  handle_range );
  for (unsigned i = 0; i < handle_list.size(); ++i)
    entities.erase( handle_list[i] );
  
    // try getting/setting single handle value
    
  std::vector<const void*> list(1);
  if (set_by_pointer) {
    list[0] = some_values;
    rval = mb.tag_set_data( tag, &one_handle, 1, &list[0] );
  }
  else {
    rval = mb.tag_set_data( tag, &one_handle, 1, some_values );
  }
  CHECK_ERR( rval );
  data.resize( bytes );
  if (get_by_pointer) {
      // test that correct size is returned
    int rsize;
    rval = mb.tag_get_data( tag, &one_handle, 1, &list[0], &rsize );
    CHECK_ERR( rval );
    CHECK_EQUAL( (int)bytes, rsize );
      // try again with NULL size pointer
    list[0] = 0;
    rval = mb.tag_get_data( tag, &one_handle, 1, &list[0] );
  }
  else {
    rval = mb.tag_get_data( tag, &one_handle, 1, &data[0] );
    list[0] = &data[0];
  }
  CHECK_ERR( rval );
  CHECK( !memcmp( some_values, list[0], bytes ) );

  
    // try getting/setting for arrays of handles
    
  const int step = std::min( (int)handle_list.size(), num_values );
  data.resize( step * bytes );
  for (int i = 0; i < (int)handle_list.size(); i += step) {
    const int n = std::min( (int)handle_list.size() - i, step );
    list.resize(n);
    if (set_by_pointer) {
      concat_to_list( some_values, list, bytes );
      rval = mb.tag_set_data( tag, &handle_list[i], n, &list[0] );
    }
    else
      rval = mb.tag_set_data( tag, &handle_list[i], n, some_values );
    CHECK_ERR( rval );
    
    if (get_by_pointer) {
        // check that valid sizes are returned if requested
      std::vector<int> rsizes(n,0);
      rval = mb.tag_get_data( tag, &handle_list[i], n, &list[0], &rsizes[0] );
      CHECK_ERR( rval );
      for (int j = 0; j < n; ++j)
        CHECK_EQUAL( (int)bytes, rsizes[j] );
        // query a second time to verify that it works w/ NULL size array
      list.clear();
      list.resize( n, 0 );
      rval = mb.tag_get_data( tag, &handle_list[i], n, &list[0] );
    }
    else {
      rval = mb.tag_get_data( tag, &handle_list[i], n, &data[0] );
      concat_to_list( &data[0], list, bytes );
    }
    CHECK_ERR( rval );
    
    const unsigned char* ptr = reinterpret_cast<const unsigned char*>(some_values);
    for (int j = 0; j < n; ++j, ptr += bytes)
      CHECK( !memcmp( ptr, list[j], bytes ) );
  }
  
    // try getting/setting for MBRange of handles
    
    // set data for range
  if (set_by_pointer) {
    list.resize( num_values );
    concat_to_list( some_values, list, bytes );
    while (list.size() < handle_range.size()) {
      size_t s = list.size();
      list.resize( 2*s );
      std::copy( list.begin(), list.begin()+s, list.begin()+s );
    }
    rval = mb.tag_set_data( tag, handle_range, &list[0] );
  }
  else {
    std::vector<unsigned char> input_data( handle_range.size() * bytes );
    for (int i = 0; i < (int)input_data.size(); i += num_values*bytes)
      memcpy( &input_data[i], some_values, std::min((unsigned)input_data.size()-i,num_values*bytes) );
    rval = mb.tag_set_data( tag, handle_range, &input_data[0] );
  }
  CHECK_ERR( rval );
  
    // get data for range
  list.clear();
  list.resize( handle_range.size(), 0 );
  if (get_by_pointer) {
      // check that valid sizes are returned if requested
    std::vector<int> rsizes( handle_range.size(), 0 );
    rval = mb.tag_get_data( tag, handle_range, &list[0], &rsizes[0] );
    CHECK_ERR(rval);
    for (size_t j = 0; j <handle_range.size(); ++j)
      CHECK_EQUAL( (int)bytes, rsizes[j] );
      // query w/ NULL size array to make sure that works also
    list.clear();
    list.resize( handle_range.size(), 0 );
    rval = mb.tag_get_data( tag, handle_range, &list[0] );
  }
  else {
    data.resize( handle_range.size() * bytes );
    concat_to_list( &data[0], list, bytes );
    rval = mb.tag_get_data( tag, handle_range, &data[0] );
  }
  CHECK_ERR( rval );
  
    // compare values
  for (size_t i = 0; i < list.size(); ++i) {
    CHECK( list[i] != NULL );
    const void* ptr = reinterpret_cast<const char*>(some_values) + (i % num_values) * bytes;
    CHECK( !memcmp( list[i], ptr, bytes ) );
  }
  
  
    // try getting unset values
 
  list.resize( entities.size() );
  if (get_by_pointer) {
    rval = mb.tag_get_data( tag, entities, &list[0] );
  }
  else {
    data.clear(),
    data.resize( entities.size() * bytes, '\001' );
    concat_to_list( &data[0], list, bytes );
    rval = mb.tag_get_data( tag, entities, &data[0] );
  }
    // if there was a default value, we should have gotten it for all unset entities
  if (default_value) {
    CHECK_ERR( rval );
    for (unsigned i = 0; i < entities.size(); ++i)
      CHECK( !memcmp( default_value, list[i], bytes ) );
  }
    // otherwise we should get MB_TAG_NOT_FOUND, *unless* the tag
    // is dense, in which case we /might/ get all zero bytes instead.
  else if (MB_TAG_NOT_FOUND != rval) {
    CHECK_EQUAL( MB_TAG_DENSE, storage );
    std::vector<unsigned char> zeros( bytes, 0 );
    for (unsigned i = 0; i < entities.size(); ++i) 
      CHECK( !memcmp( &zeros[0], list[i], bytes ) );
  }
  
    // Check that handles for other entities didn't change.
    
    // Ignore get_by_pointer/set_by_pointer flags from here on.
    // We've estabilished (hopefully) that both methods work in
    // the above code.  Now we just want to verify correct values,
    // regarless of the API used to get them.
  
    // one handle
  data.resize( bytes );
  rval = mb.tag_get_data( tag, &one_handle, 1, &data[0] );
  CHECK_ERR( rval );
  CHECK( !memcmp( some_values, &data[0], bytes ) );
  
    // array values
  data.resize( step * bytes );
  for (int i = 0; i < (int)handle_list.size(); i += step) {
    const int n = std::min( (int)handle_list.size() - i, step );
    rval = mb.tag_get_data( tag, &handle_list[i], n, &data[0] );
    CHECK_ERR( rval );
    CHECK( !memcmp( some_values, &data[0], step*bytes ) );
    rval = mb.tag_set_data( tag, &handle_list[i], n, some_values );
    CHECK_ERR( rval );
  }
  
    // range values
  list.clear();
  list.resize( handle_range.size(), 0 );
  rval = mb.tag_get_data( tag, handle_range, &list[0] );
  CHECK_ERR( rval );
  for (size_t i = 0; i< handle_range.size(); ++i) {
    const void* ptr = reinterpret_cast<const char*>(some_values) + (i % num_values) * bytes;
    CHECK( !memcmp( ptr, list[i], bytes ) );
  }
}


void test_get_set_sparse_int()
{ 
  const int data[] = { 21, 00, 46, 30, 26, 63, 05, 49, 
                       31, 39, 86, 77, 24, 37, 25, 98, 
                       26, 20, 01, 54, 16, 28, 55, 49, 
                       96, 18, 28, 18, 53, 00, 80, 48 };
  const int num_val = sizeof(data)/sizeof(data[0]);
  
  test_get_set( "sparse_int", sizeof(int)*2, 
                MB_TAG_SPARSE, MB_TYPE_INTEGER, 
                data, num_val/2, 
                0 );

  const int defaultval = 19740508;
  test_get_set( "sparse_int_def", sizeof(int), 
                MB_TAG_SPARSE, MB_TYPE_INTEGER, 
                data, num_val, 
                &defaultval );
}

void test_get_set_dense_int()
{
  const int data[] = { 231, 416, 294, 504, 318, 558, 494, 006, 
                       464, 648, 737, 045, 179, 852, 944, 336,
                       773, 248, 434, 615, 677, 667, 521, 748,
                       820, 533, 955, 300, 108, 726, 747, 597 };
  const int num_val = sizeof(data)/sizeof(data[0]);

  
  test_get_set( "dense_int", sizeof(int), 
                MB_TAG_DENSE, MB_TYPE_INTEGER, 
                data, num_val, 
                0 );

  const int defaultval[] = { 5, 8, 1974 };
  test_get_set( "dense_int_def", sizeof(int)*3, 
                MB_TAG_DENSE, MB_TYPE_INTEGER, 
                data, num_val/3, 
                defaultval );
}

void test_get_set_dense_double()
{
  const double pi = 3.1415926535897931;
  const double e  = 2.7182818284590451;
  const double data[] = { 1, 1.,     pi, e,  
                          2, 1./2, 2*pi, e/2,
                          3, 1./3, 3*pi, e/3,
                          4, 1./4, 4*pi, e/4,
                          5, 1./5, 5*pi, e/5,
                          6, 1./6, 6*pi, e/6,
                          0, 100, 1000, 10000 };
  const int num_val = sizeof(data)/sizeof(data[0]);
  
  
  test_get_set( "dense_dbl", sizeof(double), 
                MB_TAG_DENSE, MB_TYPE_DOUBLE, 
                data, num_val, 
                0 );

  const double defaultval[] = { 0.11, 0.22 };
  test_get_set( "dense_dbl_def", sizeof(double)*2, 
                MB_TAG_DENSE, MB_TYPE_DOUBLE, 
                data, num_val/2, 
                defaultval );
}


void test_get_pointers_sparse()
{
  const double data[] = { 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16 };
  const int num_val = sizeof(data)/sizeof(data[0]);
  
  test_get_set( "sparse_dbl_ptr", 2*sizeof(double), 
                 MB_TAG_SPARSE, MB_TYPE_DOUBLE,
                 data, num_val/2, 0, 
                 false, true );
  
  const double defaultval[] = { -1, -2 };
  test_get_set( "sparse_dbl_ptr_def", 2*sizeof(double), 
                 MB_TAG_SPARSE, MB_TYPE_DOUBLE,
                 data, num_val/2, defaultval, 
                 false, true );
}
  

void test_get_pointers_dense()
{
  const unsigned char data[] = "a few aribtrary bytes entered as a string";
  const int num_val = sizeof(data)/sizeof(data[0]);
  
  test_get_set( "dense_byte_ptr", 2, 
                 MB_TAG_DENSE, MB_TYPE_OPAQUE,
                 data, num_val/2, 0, 
                 false, true );
  
  const unsigned char defaultval[] = "XY";
  test_get_set( "dense_byte_ptr_def", 2, 
                 MB_TAG_DENSE, MB_TYPE_OPAQUE,
                 data, num_val/2, defaultval, 
                 false, true );
}

void test_set_pointers_sparse()
{
  const double data[] = { 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16 };
  const int num_val = sizeof(data)/sizeof(data[0]);
  
  test_get_set( "sparse_dbl_ptr", 2*sizeof(double), 
                 MB_TAG_SPARSE, MB_TYPE_DOUBLE,
                 data, num_val/2, 0, 
                 true, false );
  
  const double defaultval[] = { -1, -2 };
  test_get_set( "sparse_dbl_ptr_def", 2*sizeof(double), 
                 MB_TAG_SPARSE, MB_TYPE_DOUBLE,
                 data, num_val/2, defaultval, 
                 true, false );
}

void test_set_pointers_dense()
{
  const unsigned char data[] = "a few aribtrary bytes entered as a string";
  const int num_val = sizeof(data)/sizeof(data[0]);
  
  test_get_set( "dense_byte_ptr", 2, 
                 MB_TAG_DENSE, MB_TYPE_OPAQUE,
                 data, num_val/2, 0, 
                 true, false );
  
  const unsigned char defaultval[] = "XY";
  test_get_set( "dense_byte_ptr_def", 2, 
                 MB_TAG_DENSE, MB_TYPE_OPAQUE,
                 data, num_val/2, defaultval, 
                 true, false );
}

void test_get_set_bit()
{
    // create mesh and tag
  MBCore moab;
  MBInterface& mb = moab;
  setup_mesh( mb );
  MBTag tag = test_create_tag( mb, "bit_val", 2, MB_TAG_BIT, MB_TYPE_BIT, 0 );
  
    // get some handles to work with
  MBRange entities;
  MBErrorCode rval = mb.get_entities_by_handle( 0, entities );
  CHECK_ERR(rval);
  CHECK( !entities.empty() );
  
    // set bits on every entity
  unsigned counter = 0;
  for (MBRange::iterator i = entities.begin(); i != entities.end(); ++i) {
    srand( counter++ );
    unsigned char bits = (unsigned char)(rand() & 3);
    rval = mb.tag_set_data( tag, &*i, 1, &bits );
    unsigned char bits_out = 0;
    rval = mb.tag_get_data( tag, &*i, 1, &bits_out );
    CHECK_EQUAL( bits, bits_out );
  }
  
    // test default value
  unsigned char defval = '\003';
  unsigned char zero = '\0';
  MBTag tag2 = test_create_tag( mb, "bit_val2", 3, MB_TAG_BIT, MB_TYPE_BIT, &defval );
  CHECK( entities.size() >= 3 );
  MBRange::iterator j = entities.begin();
  MBEntityHandle h1 = *j; ++j;
  MBEntityHandle h2 = *j; ++j;
  MBEntityHandle h3 = *j; ++j;
  rval = mb.tag_set_data( tag2, &h1, 1, &zero );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( tag2, &h3, 1, &zero );
  CHECK_ERR(rval);
  unsigned char byte;
  rval = mb.tag_get_data( tag2, &h2, 1, &byte );
  CHECK_ERR( rval );
  CHECK_EQUAL( defval, byte );
  rval = mb.tag_get_data( tag2, &h1, 1, &byte );
  CHECK_ERR( rval );
  CHECK_EQUAL( zero, byte );
  rval = mb.tag_get_data( tag2, &h3, 1, &byte );
  CHECK_ERR( rval );
  CHECK_EQUAL( zero, byte );
  
    // test default value for uninitialized data (tag not set for any entity)
  defval = '\002';
  MBTag tag3 = test_create_tag( mb, "bit_val3", 2, MB_TAG_BIT, MB_TYPE_BIT, &defval );
  rval = mb.tag_get_data( tag3, &h2, 1, &byte );
  CHECK_ERR( rval );
  CHECK_EQUAL( defval, byte );
}

void test_get_by_tag( )
{
    // create mesh and tag
  MBCore moab;
  MBInterface& mb = moab;
  setup_mesh( mb );
  MBTag tag = test_create_tag( mb, "sparse_count", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, 0 );
  
    // get some handles to work with
  MBRange entities;
  MBErrorCode rval = mb.get_entities_by_type( 0, MBVERTEX, entities );
  CHECK_ERR(rval);
  CHECK( !entities.empty() );

    // get five handles
  CHECK( entities.size() > 6 );
  MBEntityHandle arr[5] = { *(entities.begin() +=   entities.size()/6),
                            *(entities.begin() += 2*entities.size()/6),
                            *(entities.begin() += 3*entities.size()/6),
                            *(entities.begin() += 4*entities.size()/6),
                            *(entities.begin() += 5*entities.size()/6) };
  int values[5] = { 1, 2, 3, 4, 5 };
  rval = mb.tag_set_data( tag, arr, 5, values );
  CHECK_ERR( rval );
  const void* const valarr[1] = { 0 };
  
    // put some in a mesh set
  MBEntityHandle set;
  const int num_in_set = 3;
  rval = mb.create_meshset( 0, set );
  CHECK_ERR(rval);
  rval = mb.add_entities( set, arr, num_in_set );
  CHECK_ERR(rval);
  
    // try for whole mesh will null tag value array
  int count = -1;
  rval = mb.get_number_entities_by_type_and_tag( 0, MBVERTEX, &tag, 0, 1, count );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5, count );
  
    // try for whole mesh will null tag value, but non-null array
  count = -1;
  rval = mb.get_number_entities_by_type_and_tag( 0, MBVERTEX, &tag, valarr, 1, count );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5, count );
  
    // try for mesh set
  rval = mb.get_number_entities_by_type_and_tag( set, MBVERTEX, &tag, 0, 1, count );
  CHECK_ERR( rval );
  CHECK_EQUAL( num_in_set, count );
  
  
    // try for whole mesh will null tag value array
  MBRange found;
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag, 0, 1, found );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5u, (unsigned)found.size() );
  MBRange::iterator i = found.begin();
  CHECK_EQUAL( arr[0], *i ); ++i;
  CHECK_EQUAL( arr[1], *i ); ++i;
  CHECK_EQUAL( arr[2], *i ); ++i;
  CHECK_EQUAL( arr[3], *i ); ++i;
  CHECK_EQUAL( arr[4], *i ); ++i;
  
    // try for whole mesh will null tag value, but non-null array
  found.clear();
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag, valarr, 1, found );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5u, (unsigned)found.size() );
  i = found.begin();
  CHECK_EQUAL( arr[0], *i ); ++i;
  CHECK_EQUAL( arr[1], *i ); ++i;
  CHECK_EQUAL( arr[2], *i ); ++i;
  CHECK_EQUAL( arr[3], *i ); ++i;
  CHECK_EQUAL( arr[4], *i ); ++i;
  
    // try for mesh set
  found.clear();
  rval = mb.get_entities_by_type_and_tag( set, MBVERTEX, &tag, 0, 1, found );
  CHECK_ERR( rval );
  CHECK_EQUAL( 3u, (unsigned)found.size() );
  i = found.begin();
  CHECK_EQUAL( arr[0], *i ); ++i;
  CHECK_EQUAL( arr[1], *i ); ++i;
  CHECK_EQUAL( arr[2], *i ); ++i;
}

void test_get_by_tag_value( )
{
    // create mesh and tag
  MBCore moab;
  MBInterface& mb = moab;
  setup_mesh( mb );
  MBTag tag = test_create_tag( mb, "sparse_count", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, 0 );
  
    // get some handles to work with
  MBRange entities;
  MBErrorCode rval = mb.get_entities_by_type( 0, MBVERTEX, entities );
  CHECK_ERR(rval);
  CHECK( !entities.empty() );

    // get five handles
  CHECK( entities.size() > 6 );
  MBEntityHandle arr[5] = { *(entities.begin() +=   entities.size()/6),
                            *(entities.begin() += 2*entities.size()/6),
                            *(entities.begin() += 3*entities.size()/6),
                            *(entities.begin() += 4*entities.size()/6),
                            *(entities.begin() += 5*entities.size()/6) };
  int values[5] = { 0xBEEF, 0xBEEF, 0xBEEF, 0xBEEF, 0xBEEF };
  rval = mb.tag_set_data( tag, arr, 5, values );
  CHECK_ERR( rval );
  const void* const valarr[1] = { values };
  
    // put some in a mesh set
  MBEntityHandle set;
  const int num_in_set = 3;
  rval = mb.create_meshset( 0, set );
  CHECK_ERR(rval);
  rval = mb.add_entities( set, arr, num_in_set );
  CHECK_ERR(rval);
  
    // try for whole mesh
  int count = -1;
  rval = mb.get_number_entities_by_type_and_tag( 0, MBVERTEX, &tag, valarr, 1, count );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5, count );
  
    // try for mesh set
  rval = mb.get_number_entities_by_type_and_tag( set, MBVERTEX, &tag, valarr, 1, count );
  CHECK_ERR( rval );
  CHECK_EQUAL( num_in_set, count );
  
  
    // try for whole mesh 
  MBRange found;
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag, valarr, 1, found );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5u, (unsigned)found.size() );
  MBRange::iterator i = found.begin();
  CHECK_EQUAL( arr[0], *i ); ++i;
  CHECK_EQUAL( arr[1], *i ); ++i;
  CHECK_EQUAL( arr[2], *i ); ++i;
  CHECK_EQUAL( arr[3], *i ); ++i;
  CHECK_EQUAL( arr[4], *i ); ++i;
  
    // try for mesh set
  found.clear();
  rval = mb.get_entities_by_type_and_tag( set, MBVERTEX, &tag, valarr, 1, found );
  CHECK_ERR( rval );
  CHECK_EQUAL( 3u, (unsigned)found.size() );
  i = found.begin();
  CHECK_EQUAL( arr[0], *i ); ++i;
  CHECK_EQUAL( arr[1], *i ); ++i;
  CHECK_EQUAL( arr[2], *i ); ++i;
}

void test_mesh_value( MBInterface& mb,
                      const char* tag_name,
                      unsigned tag_size,
                      MBTagType tag_storage,
                      MBDataType tag_type,
                      const void* value )
{
    // create  tag
  MBTag tag = test_create_tag( mb, tag_name, tag_size, tag_storage, tag_type, 0 );
  
  unsigned memcmp_size = tag_size;
  if (tag_storage == MB_TAG_BIT || tag_type == MB_TYPE_BIT)
    memcmp_size = 1;
  
  MBErrorCode rval = mb.tag_set_data( tag, 0, 0, value );
  CHECK_ERR(rval);
  std::vector<unsigned char> bytes(memcmp_size, 0);
  rval = mb.tag_get_data( tag, 0, 0, &bytes[0] );
  CHECK_ERR(rval);
  CHECK( !memcmp( value, &bytes[0], memcmp_size ) );
  
    // test again, this time for default value
  std::string name2( tag_name );
  name2 += "_DEF";
  MBTag tag2 = test_create_tag( mb, name2.c_str(), tag_size, tag_storage, tag_type, value );
  bytes.clear();
  bytes.resize(memcmp_size, 0);
  rval = mb.tag_get_data( tag2, 0, 0, &bytes[0] );
  CHECK_ERR(rval);
  CHECK( !memcmp( value, &bytes[0], memcmp_size ) );
}

void test_mesh_value()
{
  MBCore moab;
  
  double dval = -0.5;
  test_mesh_value( moab, "mvd", sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, &dval );

  int sval = 42;
  test_mesh_value( moab, "mvs", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, &sval );

  MBEntityHandle mval = 0;
  test_mesh_value( moab, "mvm", sizeof(MBEntityHandle), MB_TAG_MESH, MB_TYPE_HANDLE, &mval );

  unsigned char bits = '\002';
  test_mesh_value( moab, "mvb", 2, MB_TAG_BIT, MB_TYPE_BIT, &bits );
}

static void test_delete_type_tag( MBTagType storage )
{
  MBCore moab;
  MBInterface &mb = moab;
  MBErrorCode rval;

  setup_mesh( mb );
  
    // create tag
  int default_val = 42;
  const char* tagname = "dead_tag";
  MBTag tag = test_create_tag( mb, tagname, sizeof(int), storage, MB_TYPE_INTEGER, &default_val );
  
    // get an entity handle to work with
  MBRange verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR( rval );
  CHECK( !verts.empty() );
  MBEntityHandle handle = verts.front();
  
    // set tag value on entity
  int value = -5;
  if (storage != MB_TAG_MESH) {
    rval = mb.tag_set_data( tag, &handle, 1, &value );
    CHECK_ERR( rval );
  }
  
    // set tag value on mesh
  value = 2;
  rval = mb.tag_set_data( tag, 0, 0, &value );
  
    // delete tag
  rval = mb.tag_delete( tag );
  CHECK_ERR( rval );
  
    // make sure all basic queries fail with MB_TAG_NOT_FOUND
  std::string name;
  rval = mb.tag_get_name( tag, name );
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  MBTag tag2;
  rval = mb.tag_get_handle( tagname, tag2 );
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  int size;
  rval = mb.tag_get_size( tag, size );
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
    // get get the type from the handle, so this still succeeds
  //MBTagType storage2;
  //rval = mb.tag_get_type( tag, storage2 );
  //CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  MBDataType type;
  rval = mb.tag_get_data_type( tag, type );
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  rval = mb.tag_get_default_value( tag, &value );
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  
    // check global list of tags
  std::vector<MBTag> tags;
  rval = mb.tag_get_tags( tags );
  CHECK_ERR( rval );
  CHECK( std::find( tags.begin(), tags.end(), tag ) == tags.end() );
  
    // check tags on entity
  tags.clear();
  rval = mb.tag_get_tags_on_entity( handle, tags );
  CHECK_ERR( rval );
  CHECK( std::find( tags.begin(), tags.end(), tag ) == tags.end() );
  
    // check get_data for entity
  rval = mb.tag_get_data( tag, &handle, 1, &value );  
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  
    // check get_data for mesh
  rval = mb.tag_get_data( tag, 0, 0, &value );  
  CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  
    // check that a new tag w/ the same name can be created
  tag = test_create_tag( mb, tagname, sizeof(double), storage, MB_TYPE_DOUBLE, 0 );
  rval = mb.tag_delete( tag );
  CHECK_ERR( rval );
}


template <typename Container>  
MBErrorCode set_bit_data( MBInterface& mb, 
                          MBTag tag,
                          const Container& handles,
                          const std::vector<unsigned char>& data )
{
  MBErrorCode rval;
  data.resize( handles.size() );
  std::vector<unsigned char>::const_iterator j = data.begin();
  for (typename Container::const_iterator i = handles.begin(); i != handles.end(); ++i, ++j) {
    rval = mb.tag_set_data( tag, &*i, 1, &*j );
    if (MB_SUCCESS != rval)
      return rval;
  }
}

static bool contains_tag( MBTag tag, const std::vector<MBTag>& list )
{
  return std::find( list.begin(), list.end(), tag ) != list.end();
}

void test_get_entity_tags()
{
  MBCore moab;
  MBInterface &mb = moab;
  MBErrorCode rval;

    // get 8 handles to work with
  setup_mesh( mb );
  MBRange entities;
  rval = mb.get_entities_by_handle( 0, entities );
  CHECK_ERR( rval );
  CHECK( entities.size() >= 8 );
  MBRange::iterator i = entities.begin();
  MBEntityHandle sparse_ent       = *i; ++i;
  MBEntityHandle dense_ent        = *i; ++i;
  MBEntityHandle bit_ent          = *i; ++i;
  MBEntityHandle sparse_dense_ent = *i; ++i;
  MBEntityHandle sparse_bit_ent   = *i; ++i;
  MBEntityHandle dense_bit_ent    = *i; ++i;
  MBEntityHandle all_tag_ent      = *i; ++i;
  MBEntityHandle no_tag_ent       = *i; ++i;
  
    // create three tags to work with
  MBTag sparse, dense, bit;
  sparse = test_create_tag( mb, "sparse", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, 0 );
  dense  = test_create_tag( mb, "dense_", sizeof(int), MB_TAG_DENSE , MB_TYPE_INTEGER, 0 );
  bit    = test_create_tag( mb, "bit___", sizeof(int), MB_TAG_BIT   , MB_TYPE_BIT    , 0 );
  
    // set tags on handles
  MBEntityHandle sparse_ents[4] = { sparse_ent, sparse_dense_ent, sparse_bit_ent, all_tag_ent };
  MBEntityHandle  dense_ents[4] = { dense_ent, sparse_dense_ent, dense_bit_ent, all_tag_ent };
  MBEntityHandle    bit_ents[4] = { bit_ent, sparse_bit_ent, dense_bit_ent, all_tag_ent };
  int values[4] = { -1, -2, -3, -4 };
  rval = mb.tag_set_data(  sparse, sparse_ents, 4, &values );
  rval = mb.tag_set_data(   dense,  dense_ents, 4, &values );
  for (int j = 0; j < 4; ++j) {
    unsigned char bitval = 0xF;
    rval = mb.tag_set_data( bit, bit_ents+j, 1, &bitval );
    CHECK_ERR(rval);
  }
  
    // get tags on each entity
  std::vector<MBTag> sparse_ent_tags, dense_ent_tags, bit_ent_tags, 
                     sparse_dense_ent_tags, sparse_bit_ent_tags, dense_bit_ent_tags,
                     all_tag_ent_tags, no_tag_ent_tags;
  rval = mb.tag_get_tags_on_entity( sparse_ent, sparse_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( dense_ent, dense_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( bit_ent, bit_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( sparse_dense_ent, sparse_dense_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( sparse_bit_ent, sparse_bit_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( dense_bit_ent, dense_bit_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( all_tag_ent, all_tag_ent_tags );
  CHECK_ERR(rval);
  rval = mb.tag_get_tags_on_entity( no_tag_ent, no_tag_ent_tags );
  CHECK_ERR(rval);
  
    // check expected values
    // NOTE:  could potentially get back bit and dense tags for any entity
    //        depending on the storage layout.  Only sparse tags guarantee
    //        no false positives.  False negatives should never happen for any type.
    //        Also, there could be other dense tags already defined.

    // verify sparse tag in all expected lists
  CHECK( contains_tag( sparse, sparse_ent_tags ) );
  CHECK( contains_tag( sparse, sparse_dense_ent_tags ) );
  CHECK( contains_tag( sparse, sparse_bit_ent_tags ) );
  CHECK( contains_tag( sparse, all_tag_ent_tags ) );
  
    // verify sparse tag not in any other lists
  CHECK( !contains_tag( sparse, dense_ent_tags ) );
  CHECK( !contains_tag( sparse, bit_ent_tags ) );
  CHECK( !contains_tag( sparse, dense_bit_ent_tags ) );
  CHECK( !contains_tag( sparse, no_tag_ent_tags ) );

    // verify dense tag in all expected lists
  CHECK( contains_tag( dense, dense_ent_tags ) );
  CHECK( contains_tag( dense, sparse_dense_ent_tags ) );
  CHECK( contains_tag( dense, dense_bit_ent_tags ) );
  CHECK( contains_tag( dense, all_tag_ent_tags ) );

    // verify bit tag in all expected lists
  CHECK( contains_tag( bit, bit_ent_tags ) );
  CHECK( contains_tag( bit, sparse_bit_ent_tags ) );
  CHECK( contains_tag( bit, dense_bit_ent_tags ) );
  CHECK( contains_tag( bit, all_tag_ent_tags ) );
}

void test_delete_sparse_tag()
{
  test_delete_type_tag( MB_TAG_SPARSE );
}

void test_delete_dense_tag()
{
  test_delete_type_tag( MB_TAG_DENSE );
}

void test_delete_mesh_tag()
{
  test_delete_type_tag( MB_TAG_MESH );
}
  

void test_delete_tag_data( MBTagType storage, bool with_default_value )
{
  MBCore moab;
  MBInterface &mb = moab;
  MBErrorCode rval;

  setup_mesh( mb );
  
    // subdivide entity handles into three groups:
    // 1) entities for which the tag data will be deleted using the array-based function
    // 2) entities for which the tag data will be deleted using the range-based function
    // 3) entities for which the tag data will not be deleted
  MBRange all_entities, del1_range, keep_range;
  std::vector<MBEntityHandle> del1_list, del2_list, keep_list;
  rval = mb.get_entities_by_handle( 0, all_entities );
  CHECK_ERR( rval );
  int c = 0;
  for (MBRange::iterator i = all_entities.begin(); i != all_entities.end(); ++i, ++c)  {
    switch (c%3) {
      case 0: del1_range.insert( *i ); break;
      case 1: keep_range.insert( *i ); break;
      case 2: del2_list.push_back( *i ); break;
    }
  }
  del1_list.resize( del1_range.size() );
  std::copy( del1_range.begin(), del1_range.end(), del1_list.begin() );
  keep_list.resize( keep_range.size() );
  std::copy( keep_range.begin(), keep_range.end(), keep_list.begin() );
 
    // create tag
  MBEntityHandle first = all_entities.front();
  MBEntityHandle* defval = with_default_value ? &first : 0;
  const char* tagname = "dead_tag";
  MBTag tag = test_create_tag( mb, tagname, sizeof(MBEntityHandle), storage, 
                               MB_TYPE_HANDLE, defval );
                               
    // set value for each entity to its handle
  rval = mb.tag_set_data( tag, del1_range, &del1_list[0] );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( tag, keep_range, &keep_list[0] );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( tag, &del2_list[0], del2_list.size(), &del2_list[0] );
  CHECK_ERR(rval);
  
    // delete tag data
  rval = mb.tag_delete_data( tag, del1_range );
  CHECK_ERR(rval);
  rval = mb.tag_delete_data( tag, &del2_list[0], del2_list.size() );
  CHECK_ERR(rval);
  
    // test that keep list is unaffected
  std::vector<MBEntityHandle> tag_data( keep_range.size(), 0 );
  rval = mb.tag_get_data( tag, keep_range, &tag_data[0] );
  CHECK_ERR(rval);
  CHECK( tag_data == keep_list );
  
    // try to get data for deleted range
  tag_data.clear();
  tag_data.resize( del1_range.size(), (MBEntityHandle)-1 );
  rval = mb.tag_get_data( tag, del1_range, &tag_data[0] );
    // if we have a default value, should get that for deleted entities
  if (with_default_value) {
    CHECK_ERR(rval);
    std::vector<MBEntityHandle> expected( del1_range.size(), *defval );
    CHECK( expected == tag_data );
  }
  else if (rval != MB_TAG_NOT_FOUND) {
      // dense and bit tags might return either MB_TAG_NOT_FOUND or zero bytes.
      // sparse tags should always return MB_TAG_NOT_FOUND
    CHECK( MB_TAG_SPARSE != storage );
    std::vector<MBEntityHandle> expected( del1_range.size(), 0 );
    CHECK( expected == tag_data );
  }
  
    // try to get data for deleted list
  tag_data.clear();
  tag_data.resize( del1_range.size(), (MBEntityHandle)-1 );
  rval = mb.tag_get_data( tag, del1_range, &tag_data[0] );
    // if we have a default value, should get that for deleted entities
  if (with_default_value) {
    CHECK_ERR(rval);
    std::vector<MBEntityHandle> expected( del1_range.size(), *defval );
    CHECK( expected == tag_data );
  }
  else if (rval != MB_TAG_NOT_FOUND) {
      // dense and bit tags might return either MB_TAG_NOT_FOUND or zero bytes.
      // sparse tags should always return MB_TAG_NOT_FOUND
    CHECK( MB_TAG_SPARSE != storage );
    std::vector<MBEntityHandle> expected( del1_range.size(), 0 );
    CHECK( expected == tag_data );
  }
}

void test_delete_sparse_data()
{ 
  test_delete_tag_data( MB_TAG_DENSE, false );
  test_delete_tag_data( MB_TAG_DENSE, true  );
}

void test_delete_dense_data()
{ 
  test_delete_tag_data( MB_TAG_SPARSE, false );
  test_delete_tag_data( MB_TAG_SPARSE, true  );
}

void test_delete_bit_data()
{ 
  MBCore moab;
  MBInterface &mb = moab;
  MBErrorCode rval;

    // get an entity to set data on
  setup_mesh( mb );
  MBRange entities;
  rval = mb.get_entities_by_handle( 0, entities );
  CHECK_ERR( rval );
  CHECK( !entities.empty() );
  MBEntityHandle handle = entities.front();
  
    // create two tags, one with a default value and one without
  unsigned char defval = '\006';  // 110
  MBTag tag1, tag2;
  tag1 = test_create_tag( mb, "tag1", 2, MB_TAG_BIT, MB_TYPE_BIT, 0 );
  tag2 = test_create_tag( mb, "tag2", 3, MB_TAG_BIT, MB_TYPE_BIT, &defval );
  
    // set value for each tag
  unsigned char val = '\001';
  rval = mb.tag_set_data( tag1, &handle, 1, &val );
  CHECK_ERR( rval );
  rval = mb.tag_set_data( tag2, &handle, 1, &val );
  CHECK_ERR( rval );
  
    // delete both tag values
  rval = mb.tag_delete_data( tag1, &handle, 1 );
  CHECK_ERR( rval );
  rval = mb.tag_delete_data( tag2, &handle, 1 );
  CHECK_ERR( rval );
  
    // try to get value for tag w/out default.
    // should get back either not found, or a zero value.
  val = 0xFF;
  rval = mb.tag_get_data( tag1, &handle, 1, &val );
  if (MB_SUCCESS == rval) 
    CHECK_EQUAL( val, '\0' );
  else
    CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  
    // test that we get back default value for second tag
  val = 0xFF;
  rval = mb.tag_get_data( tag2, &handle, 1, &val );
  CHECK_ERR( rval );
  CHECK_EQUAL( defval, val );
}

void test_get_set_variable_length( const char* name,
                                   MBTagType storage, 
                                   MBDataType type, 
                                   const void** values,
                                   const int* lengths,
                                   int num_values,
                                   const void* default_value,
                                   int default_value_length )
{
  std::vector<const void*> data;
  std::vector<int> data_lens;
  
    // create mesh and tag
  MBCore moab;
  MBInterface& mb = moab;
  setup_mesh( mb );
  MBTag tag = test_create_var_len_tag( mb, name,  storage, type, 
                           default_value, default_value_length );
  
    // get some handles to work with
  MBRange entities;
  MBErrorCode rval = mb.get_entities_by_handle( 0, entities );
  CHECK_ERR(rval);
  CHECK( !entities.empty() );
  
    // split handles into four groups
    // a) a single handle
    // b) some non-consecutive handles in an array
    // c) some handles in an MBRange
    // d) remaining handles (remaining in 'entities');
  MBEntityHandle one_handle;
  MBRange::iterator i = entities.begin() += entities.size()/2;
  one_handle = *i;
  entities.erase( i );
  
  MBRange handle_range;
  std::vector<MBEntityHandle> handle_list;
  for (MBRange::const_pair_iterator i =  entities.const_pair_begin(); 
                                    i != entities.const_pair_end(); ++i) {
    if (i->first == i->second || i->second - i->first == 1) {
      MBEntityHandle h1 = i->first, h2 = i->second;
      ++i;
      handle_range.insert( h1, h2 );
    }
    else {
      MBEntityHandle mid = (MBEntityHandle)(i->first + (i->second - i->first + 1) / 2);
      handle_list.push_back( mid );
      handle_range.insert( mid+1, i->second );
    }
  }
  entities = subtract( entities,  handle_range );
  for (unsigned i = 0; i < handle_list.size(); ++i)
    entities.erase( handle_list[i] );
  
    // try getting/setting single handle value
  rval = mb.tag_set_data( tag, &one_handle, 1, values, lengths );
  CHECK_ERR(rval);
  const void* data_ptr;
  int data_len;
  rval = mb.tag_get_data( tag, &one_handle, 1, &data_ptr, &data_len );
  CHECK_ERR(rval);
  CHECK_EQUAL( lengths[0], data_len );
  CHECK( !memcmp( values[0], data_ptr, data_len ) );

  
    // try getting/setting for arrays of handles
    
  int count = std::min( (int)handle_list.size(), num_values );
  rval = mb.tag_set_data( tag, &handle_list[0], count, values, lengths );
  CHECK_ERR( rval );
  data.clear();
  data.resize( count, 0 );
  data_lens.clear();
  data_lens.resize( count, 0 );
  rval = mb.tag_get_data( tag, &handle_list[0], count, &data[0], &data_lens[0] );
  CHECK_ERR( rval );
  for (int i = 0; i < count; ++i) {
    CHECK_EQUAL( lengths[i], data_lens[i] );
    CHECK( NULL != data[i] );
    CHECK( !memcmp( values[i], data[i], lengths[i] ) );
  }
  
    // try getting/setting for MBRange of handles
  
  data.resize( num_values );
  data_lens.resize( num_values );
  std::copy( values, values + num_values, data.begin() );
  std::copy( lengths, lengths + num_values, data_lens.begin() );
  while (data.size() < handle_range.size()) {
    size_t s = data.size();
    data.resize( 2*s );
    std::copy( data.begin(), data.begin() + s, data.begin() + s );
    data_lens.resize( 2*s );
    std::copy( data_lens.begin(), data_lens.begin() + s, data_lens.begin() + s );
  }
  rval = mb.tag_set_data( tag, handle_range, &data[0], &data_lens[0] );
  CHECK_ERR( rval );
  
  data.clear();
  data.resize( handle_range.size(), 0 );
  data_lens.clear();
  data_lens.resize( handle_range.size(), 0 );
  rval = mb.tag_get_data( tag, handle_range, &data[0], &data_lens[0] );
  CHECK_ERR( rval );
  
  for (size_t i = 0; i < data.size(); ++i) {
    const void* expect = values[i%num_values];
    int expect_len = lengths[i%num_values];
    CHECK_EQUAL( expect_len, data_lens[i] );
    CHECK( NULL != data[i] );
    CHECK( !memcmp( expect, data[i], expect_len ) );
  }

  
    // try getting unset values
 
  data.resize( entities.size() );
  data_lens.resize( entities.size() );
  rval = mb.tag_get_data( tag, entities, &data[0], &data_lens[0] );
    // if there was a default value, we should have gotten it for all unset entities
  if (default_value) {
    CHECK_ERR( rval );
    for (unsigned i = 0; i < entities.size(); ++i) {
      CHECK_EQUAL( default_value_length, data_lens[i] );
      CHECK( NULL != data[i] );
      CHECK( !memcmp( default_value, data[i], default_value_length ) );
    }
  }
    // otherwise we should get MB_TAG_NOT_FOUND
  else {
    CHECK_EQUAL( MB_TAG_NOT_FOUND, rval );
  }
  
    // Check that handles for other entities didn't change.
  
    // one handle
  rval = mb.tag_get_data( tag, &one_handle, 1, &data_ptr, &data_len );
  CHECK_ERR(rval);
  CHECK_EQUAL( lengths[0], data_len );
  CHECK( !memcmp( values[0], data_ptr, data_len ) );

    // array values
  count = std::min( (int)handle_list.size(), num_values );
  data.clear();
  data.resize( count, 0 );
  data_lens.clear();
  data_lens.resize( count, 0 );
  rval = mb.tag_get_data( tag, &handle_list[0], count, &data[0], &data_lens[0] );
  CHECK_ERR( rval );
  for (int i = 0; i < count; ++i) {
    CHECK_EQUAL( lengths[i], data_lens[i] );
    CHECK( NULL != data[i] );
    CHECK( !memcmp( values[i], data[i], lengths[i] ) );
  }
  
    // range values
  data.clear();
  data.resize( handle_range.size(), 0 );
  data_lens.clear();
  data_lens.resize( handle_range.size(), 0 );
  rval = mb.tag_get_data( tag, handle_range, &data[0], &data_lens[0] );
  CHECK_ERR( rval );
  
  for (size_t i = 0; i < data.size(); ++i) {
    const void* expect = values[i%num_values];
    int expect_len = lengths[i%num_values];
    CHECK_EQUAL( expect_len, data_lens[i] );
    CHECK( NULL != data[i] );
    CHECK( !memcmp( expect, data[i], expect_len ) );
  }
  
}

void test_get_set_variable_length_sparse()
{
  const double doubles[14] = { 1, 2, 3, 4, -4, -3, -2, -1, 42, 0, 1974, -0.5, 1./3, -1e-10 };
  const void* dvals[5] = { doubles, doubles + 3, doubles + 4, doubles + 8, doubles + 10 };
  const int ds = sizeof(double);
  const int dlens[5] = { 3*ds, 1*ds, 4*ds, 2*ds, 4*ds };
  test_get_set_variable_length( "vnodef", MB_TAG_SPARSE, MB_TYPE_DOUBLE, dvals, dlens, 5, 0, 0 );

  const int ints[32] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                        -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15 };
  const void* ivals[9] = { ints, ints+1, ints+3, ints+12, ints+17, ints+21, ints+28, ints+29, ints+31 };
  const int is = sizeof(double);
  const int ilens[9] = { 1*is, 2*is, 9*is, 5*is, 4*is, 7*is, 1*is, 2*is, 1*is };
  const int defvals[] = { 42, 5, 8, 74 };
  test_get_set_variable_length( "vdef", MB_TAG_SPARSE, MB_TYPE_INTEGER, ivals, ilens, 9, defvals, sizeof(defvals) );
}

void test_get_set_variable_length_dense()
{
  const double doubles[14] = { 1, 2, 3, 4, -4, -3, -2, -1, 42, 0, 1974, -0.5, 1./3, -1e-10 };
  const void* dvals[5] = { doubles, doubles + 3, doubles + 4, doubles + 8, doubles + 10 };
  const int ds = sizeof(double);
  const int dlens[5] = { 3*ds, 1*ds, 4*ds, 2*ds, 4*ds };
  test_get_set_variable_length( "vnodef", MB_TAG_DENSE, MB_TYPE_DOUBLE, dvals, dlens, 5, 0, 0 );

  const int ints[32] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                        -1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15 };
  const void* ivals[9] = { ints, ints+1, ints+3, ints+12, ints+17, ints+21, ints+28, ints+29, ints+31 };
  const int is = sizeof(double);
  const int ilens[9] = { 1*is, 2*is, 9*is, 5*is, 4*is, 7*is, 1*is, 2*is, 1*is };
  const int defvals[] = { 42, 5, 8, 74 };
  test_get_set_variable_length( "vdef", MB_TAG_DENSE, MB_TYPE_INTEGER, ivals, ilens, 9, defvals, sizeof(defvals) );
}

void test_get_set_variable_length_mesh()
{
  MBCore moab;
  MBInterface &mb = moab;
  MBErrorCode rval;

  MBTag tag = test_create_var_len_tag( mb, "vmesh", MB_TAG_MESH, MB_TYPE_INTEGER, 0, 0 );
  int values1[] = { 6 };
  int values5[] = { 1, 2, 3, 4, 5 };
  
  int one = sizeof(int);
  const void* data[1];
  data[0] = values1;
  rval = mb.tag_set_data( tag, 0, 0, data, &one );
  CHECK_ERR( rval );
  
  int len;
  rval = mb.tag_get_data( tag, 0, 0, data, &len );
  CHECK_ERR( rval );
  CHECK_EQUAL( (int)sizeof(int), len );
  CHECK_EQUAL( values1[0], *reinterpret_cast<const int*>(data[0]) );
  
  int five = 5*sizeof(int);
  data[0] = values5;
  rval = mb.tag_set_data( tag, 0, 0, data, &five );
  CHECK_ERR( rval );
  
  rval = mb.tag_get_data( tag, 0, 0, data, &len );
  CHECK_ERR( rval );
  CHECK_EQUAL( 5*(int)sizeof(int), len );
  CHECK_EQUAL( values5[0], reinterpret_cast<const int*>(data[0])[0] );
  CHECK_EQUAL( values5[1], reinterpret_cast<const int*>(data[0])[1] );
  CHECK_EQUAL( values5[2], reinterpret_cast<const int*>(data[0])[2] );
  CHECK_EQUAL( values5[3], reinterpret_cast<const int*>(data[0])[3] );
  CHECK_EQUAL( values5[4], reinterpret_cast<const int*>(data[0])[4] );
}

void test_get_ents_with_default_value()
{
  MBCore moab;
  MBInterface &mb = moab;
  MBErrorCode rval;
  MBRange result;

    // create a bunch of vertices
  std::vector<double> coords(90,0.0);
  MBRange verts;
  rval = mb.create_vertices( &coords[0], coords.size()/3, verts );
  CHECK_ERR( rval );
  CHECK_EQUAL( coords.size()/3, verts.size() );
    // create one edge, which we should never get back from 
    // our queries with type == MBVERTEX
  MBEntityHandle edge, ends[] = { verts.front(), verts.back() };
  rval = mb.create_element( MBEDGE, ends, 2, edge );
  CHECK_ERR(rval);
  
    // split vertices into four groups
  MBRange sets[4];
  size_t s = 0;
  for (MBRange::iterator i = verts.begin(); i != verts.end(); ++i) {
    sets[s].insert(*i);
    s = (s+1)%4;
  }


    // create a sparse tag and set some verts to non-default value
  int default_sparse = 5;
  MBTag tag_sparse = test_create_tag( mb, "int", sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER, &default_sparse );
  std::vector<int> sparse_vals(sets[0].size(), -1);
  rval = mb.tag_set_data( tag_sparse, sets[0], &sparse_vals[0] );
  CHECK_ERR(rval);
  
    // get all entities with default value for sparse tag
  result.clear();
  const void* ptrs[] = { &default_sparse };
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag_sparse, ptrs, 1, result );
  CHECK_ERR(rval);
  CHECK_EQUAL( subtract(verts, sets[0]), result );


    // create a dense tag and set some verts to non-default value
  double default_dense = -1.0;
  MBTag tag_dense = test_create_tag( mb, "double", sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, &default_dense );
  std::vector<double> dense_vals(sets[1].size(), 3.14159);
  rval = mb.tag_set_data( tag_dense, sets[1], &dense_vals[0] );
  CHECK_ERR(rval);
  
    // get all entities with default value for dense tag
  result.clear();
  ptrs[0] = &default_dense;
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag_dense, ptrs, 1, result );
  CHECK_ERR(rval);
  CHECK_EQUAL( subtract(verts, sets[1]), result );
  
  
    // create a variable-length tag and set some verts to non-default value
  // SKIP THIS: NO API FOR QUERYING ENTITIES WITH VARIABLE-LENGTH VALUE
  //int default_vlen[] = { 1, 2, 3 };
  //MBTag tag_vlen = test_create_var_len_tag( mb, "vlen", MB_TAG_SPARSE, MB_TYPE_INTEGER, default_vlen, sizeof(default_vlen) );
  //int other_vlen[] = { 4, 5, 6, 7 };
  //std::vector<const void*> vlen_ptrs( sets[2].size(), other_vlen );
  //std::vector<int> vlen_sizes( sets[2].size)(), sizeof(other_vlen) );
  //rval = mb.tag_set_data( tag_vlen, sets[2], &vlen_ptrs[0], &vlen_sizes[0] );
  //CHECK_ERR(rval);
  
  
    // check that INTERSECT option works as expected
  result.clear();
  result.insert( sets[1].front() );
  ptrs[0] = &default_sparse;
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag_sparse, ptrs, 1, result, MBInterface::INTERSECT );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)1, result.size() );
  CHECK_EQUAL( sets[1].front(), result.front() );
  
  
    // check that UNITE option works as expected
  result.clear();
  result.insert( edge );
  ptrs[0] = &default_sparse;
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag_sparse, ptrs, 1, result, MBInterface::UNION );
  CHECK_ERR(rval);
  CHECK_EQUAL( edge, result.back() );
}

void setup_mesh( MBInterface& mb )
{
  MBRange vertex_handles;
  const double vertex_coords[] = { 0, 0, 0,  1, 0, 0,  2, 0, 0,
                                   0, 1, 0,  1, 1, 0,  2, 1, 0,
                                   0, 2, 0,  1, 2, 0,  2, 2, 0,
                                   
                                   0, 0, 1,  1, 0, 1,  2, 0, 1,
                                   0, 1, 1,  1, 1, 1,  2, 1, 1,
                                   0, 2, 1,  1, 2, 1,  2, 2, 1,
                                   
                                   0, 0, 2,  1, 0, 2,  2, 0, 2,
                                   0, 1, 2,  1, 1, 2,  2, 1, 2,
                                   0, 2, 2,  1, 2, 2,  2, 2, 2 };
  const unsigned num_vtx = sizeof(vertex_coords)/(3*sizeof(double));
  MBErrorCode rval = mb.create_vertices( vertex_coords, num_vtx, vertex_handles );
  CHECK_ERR(rval);
  CHECK_EQUAL( num_vtx, (unsigned)vertex_handles.size() );
  
  CHECK_EQUAL( 27u, num_vtx );
  MBEntityHandle elements[8];
  MBEntityHandle conn[8][8] = { {  0,  1,  4,  3,   9, 10, 13, 12 },
                                {  1,  2,  5,  4,  10, 11, 14, 13 }, 
                                {  3,  4,  7,  6,  12, 13, 16, 15 },
                                {  4,  5,  8,  7,  13, 14, 17, 16 },
                                {  9, 10, 13, 12,  18, 19, 22, 21 },
                                { 10, 11, 14, 13,  19, 20, 23, 22 },
                                { 12, 13, 16, 15,  21, 22, 25, 24 },
                                { 13, 14, 17, 16,  22, 23, 26, 25 } };
  for (unsigned i = 0; i < 8; ++i) {
    rval = mb.create_element( MBHEX, conn[i], 8, elements[i] );
    CHECK_ERR(rval);
  }
  
    // delete some stuff so there are multiple sequences
  rval = mb.delete_entities( elements + 2, 2 );
  CHECK_ERR(rval);
  MBRange::iterator i = vertex_handles.begin();
  i += 16;
  rval = mb.delete_entities( &*i, 1 );
  CHECK_ERR(rval);
  i += 2;
  rval = mb.delete_entities( &*i, 1 );
  CHECK_ERR(rval);
}

/* Found bug where last entity in sequence is not
   returned for get entity by tag for a variable-
   length tag.  Test to make sure that it remains
   fixed.
 */
void regression_one_entity_by_var_tag()
{
  MBCore moab;
  MBErrorCode rval;
  
  MBEntityHandle vertex;
  const double coords[] = { 0, 0, 0 };
  rval = moab.create_vertex( coords, vertex );
  CHECK_ERR(rval);

  MBTag tag;
  rval = moab.tag_create_variable_length( "testtag", MB_TAG_DENSE, MB_TYPE_INTEGER, tag );
  CHECK_ERR(rval);
  
  int taglen = sizeof(int);
  const void* ptrarr[1] = { &taglen };
  rval = moab.tag_set_data( tag, &vertex, 1, ptrarr, &taglen );
  CHECK_ERR(rval);
  
  MBRange ents;
  rval = moab.get_entities_by_type_and_tag( 0, MBVERTEX, &tag, 0, 1, ents );
  CHECK_ERR(rval);
  
  CHECK_EQUAL( (MBEntityHandle)1, ents.size() );
  CHECK_EQUAL( vertex, ents.front() );
}

/* Return MB_ENTITY_NOT_FOUND if asked to set a tag on an 
   entity that doesn't exist.
 */
void regression_tag_on_nonexistent_entity()
{
  MBCore moab;
  MBErrorCode rval;
  const int tagval = 0xdeadbeef;
  const void* valarr[1] = { &tagval };
  const int numval = sizeof(int);
  
    // create all three types of tags
  MBTag dense, sparse, bit;
  rval = moab.tag_create( "test_dense", numval, MB_TAG_DENSE, MB_TYPE_INTEGER, dense, 0, false );
  CHECK_ERR(rval);
  rval = moab.tag_create( "test_sparse", numval, MB_TAG_SPARSE, MB_TYPE_INTEGER, sparse, 0, false );
  CHECK_ERR(rval);
  rval = moab.tag_create( "test_bit", numval, MB_TAG_BIT, MB_TYPE_BIT, bit, 0, false );
  CHECK_ERR(rval);
  
    // for each tag type, check all four mechanisms for setting tag data
    // (fixed and variable length given array or range).
  MBEntityHandle handle = (MBEntityHandle)1;
  MBRange handles;
  handles.insert( handle );
  
  rval = moab.tag_set_data( dense,  &handle, 1, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, &handle, 1, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( bit,    &handle, 1, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
  rval = moab.tag_set_data( dense,  handles, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, handles, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( bit,    handles, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
  rval = moab.tag_set_data( dense,  &handle, 1, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, &handle, 1, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
  rval = moab.tag_set_data( dense,  handles, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, handles, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
    // now add create an entity and try an adjacent handle
  MBEntityHandle set;
  rval = moab.create_meshset( 0, set );
  CHECK_ERR(rval);
  
  handle = (MBEntityHandle)(set+1);
  handles.clear();
  handles.insert( handle );
  
  rval = moab.tag_set_data( dense,  &handle, 1, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, &handle, 1, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( bit,    &handle, 1, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
  rval = moab.tag_set_data( dense,  handles, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, handles, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( bit,    handles, &tagval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
  rval = moab.tag_set_data( dense,  &handle, 1, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, &handle, 1, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  
  rval = moab.tag_set_data( dense,  handles, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
  rval = moab.tag_set_data( sparse, handles, valarr, &numval );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, rval );
}
