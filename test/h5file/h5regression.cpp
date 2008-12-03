#include "MBCore.hpp"
#include "testdir.h"
#include "TestUtil.hpp"
#include "MBRange.hpp"
#include "MBReadUtilIface.hpp"
#include "WriteHDF5.hpp"
#include "FileOptions.hpp"

#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <math.h>

const char filename[] = "bad.h5m";
  
void test_write_invalid_elem();

int main(int argc, char* argv[])
{
  int exitval = 0;
  exitval += RUN_TEST( test_write_invalid_elem );
  return exitval;
}

  
void test_write_invalid_elem()
{
  MBCore mbcore;
  MBInterface& moab = mbcore;
  MBReadUtilIface* readtool = 0;
  MBErrorCode rval;
  
  void* ptr = 0;
  rval = moab.query_interface( "MBReadUtilIface", &ptr );
  CHECK_ERR(rval);
  CHECK( ptr != 0 );
  readtool = reinterpret_cast<MBReadUtilIface*>(ptr);
  
    // create two nodes
  MBEntityHandle first_node;
  std::vector<double*> coords;
  rval = readtool->get_node_arrays( 3, 2, 1, first_node, coords );
  CHECK_ERR(rval);
  coords[0][0] = coords[0][1] = 0.0;
  coords[1][0] = coords[1][1] = 0.0;
  coords[2][0] = coords[2][1] = 0.0;
  
    // create a triangle with an invalid node handle for its
    // third vertex
  MBEntityHandle tri;
  MBEntityHandle* conn = 0;
  rval = readtool->get_element_array( 1, 3, MBTRI, 1, tri, conn );
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

  
