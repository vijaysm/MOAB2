#include "MBCore.hpp"
#include "testdir.h"
#include "TestUtil.hpp"
#include "MBRange.hpp"

#include <algorithm>
#include <iostream>
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
  

int main(int argc, char* argv[])
{
  if (argc == 2 &&  std::string(argv[1]) == "-k")
    keep_file = true;
  else if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " [-k]" << std::endl;
    return 1;
  }

  // only one test so far... should probably add second test
  // for really-old-format  entityset parent/child links
  int exitval = 0;
  exitval += RUN_TEST( test_ranged_set_with_holes );
  exitval += RUN_TEST( test_file_set );
  return exitval;
}

  
