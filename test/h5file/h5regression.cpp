#include "moab/Core.hpp"
#include "TestUtil.hpp"
#include "moab/Range.hpp"
#include "moab/ReadUtilIface.hpp"
#include "WriteHDF5.hpp"
#include "FileOptions.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

using namespace moab;

const char filename[] = "bad.h5m";
  
void test_write_invalid_elem();
void test_write_read_many_tags();

int main(int argc, char* argv[])
{
  int exitval = 0;
  exitval += RUN_TEST( test_write_invalid_elem );
  exitval += RUN_TEST( test_write_read_many_tags );
  return exitval;
}

  
void test_write_invalid_elem()
{
  Core mbcore;
  Interface& moab = mbcore;
  ReadUtilIface* readtool = 0;
  ErrorCode rval;
  
  void* ptr = 0;
  rval = moab.query_interface( "ReadUtilIface", &ptr );
  CHECK_ERR(rval);
  CHECK( ptr != 0 );
  readtool = reinterpret_cast<ReadUtilIface*>(ptr);
  
    // create two nodes
  EntityHandle first_node;
  std::vector<double*> coords;
  rval = readtool->get_node_coords( 3, 2, 1, first_node, coords );
  CHECK_ERR(rval);
  coords[0][0] = coords[0][1] = 0.0;
  coords[1][0] = coords[1][1] = 0.0;
  coords[2][0] = coords[2][1] = 0.0;
  
    // create a triangle with an invalid node handle for its
    // third vertex
  EntityHandle tri;
  EntityHandle* conn = 0;
  rval = readtool->get_element_connect( 1, 3, MBTRI, 1, tri, conn );
  CHECK_ERR(rval);
  conn[0] = first_node;   // valid
  conn[1] = first_node+1; // valid
  conn[2] = first_node+2; // invalid
  
    // try to write the file (should fail)
  WriteHDF5 writer( &moab );
  FileOptions opts(0);
  rval = writer.write_file( filename, true, opts, 0, 0, std::vector<std::string>() );
  CHECK(MB_SUCCESS != rval);
}

void test_write_read_many_tags()
{
  const int N = 200;
  Core mbcore;
  Interface& mb = mbcore;
  ErrorCode rval;
  
  double coords[3] = { 0, 0, 0 };
  EntityHandle node;
  rval = mb.create_vertex( coords, node );
  CHECK_ERR(rval);
  
    // create a lot of tags
  std::vector<Tag> tags;
  for (int i = 0; i < N; ++i) {
    Tag t;
    std::ostringstream name("IntTag");
    name << i;
    rval = mb.tag_create( name.str().c_str(),
                          sizeof(int),
                          i % 2 ? MB_TAG_SPARSE : MB_TAG_DENSE,
                          MB_TYPE_INTEGER,
                          t,
                          &i );
    CHECK_ERR(rval);
    tags.push_back(t);
  }
  
    // write the file
  rval = mb.write_file( filename, "MOAB" );
  CHECK_ERR(rval);
  
    // clear moab instance
  rval = mb.delete_mesh();
  CHECK_ERR(rval);
  for (int i = 0; i < N; ++i) {
    rval = mb.tag_delete( tags[i] );
    CHECK_ERR(rval);
  }
  
  
    // read the file
  rval = mb.load_file( filename );
  CHECK_ERR(rval);
  remove(filename);
  
    // check that we have the expected tags
  for (int i = 0; i < N; ++i) {
    Tag t;
    std::ostringstream name("IntTag");
    name << i;
    rval = mb.tag_get_handle( name.str().c_str(), t );
    CHECK_ERR(rval);
    
    TagType storage;
    rval = mb.tag_get_type( t, storage );
    CHECK_ERR(rval);
    CHECK_EQUAL( (i%2)?MB_TAG_SPARSE:MB_TAG_DENSE, storage );
    
    DataType type;
    rval = mb.tag_get_data_type( t, type );
    CHECK_ERR(rval);
    CHECK_EQUAL( MB_TYPE_INTEGER, type );
    
    int size;
    rval = mb.tag_get_size( t, size );
    CHECK_ERR(rval);
    CHECK_EQUAL( (int)sizeof(int), size );
    
    int def;
    rval = mb.tag_get_default_value( t, &def );
    CHECK_ERR(rval);
    CHECK_EQUAL( i, def );
  }
}  
