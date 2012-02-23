#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/HomXform.hpp"
#include "moab/ReadUtilIface.hpp"
#include "TestUtil.hpp"
#include <stdlib.h>
#include <algorithm>

void test_coords_connect_iterate();
void test_scd_invalid();

using namespace moab;

int main()
{
  int failures = 0;
  
  failures += RUN_TEST(test_coords_connect_iterate);
  failures += RUN_TEST(test_scd_invalid);
  
  if (failures) 
    std::cerr << "<<<< " << failures << " TESTS FAILED >>>>" << std::endl;
    
  return failures;
}

void test_coords_connect_iterate()
{
  // create 1000 vertices
  const unsigned int NUM_VTX = 1000;
  Core moab;
  Interface& mb = moab;
  std::vector<double> coords(3*NUM_VTX);
  for (unsigned int i = 0; i < NUM_VTX; i++)
    coords[3*i] = coords[3*i+1] = coords[3*i+2] = i;
    
  Range verts, hexes, faces, edges, dead;
  ErrorCode rval = mb.create_vertices( &coords[0], NUM_VTX, verts );
  CHECK_ERR(rval);

    // create a bunch of hexes from those
  ReadUtilIface *rui;
  EntityHandle *orig_connect, start_hex;
  rval = mb.query_interface(rui);
  CHECK_ERR(rval);
  rval = rui->get_element_connect(NUM_VTX/8, 8, MBHEX, 1, start_hex, orig_connect);
  CHECK_ERR(rval);
  std::copy(verts.begin(), verts.end(), orig_connect);
  hexes.insert(start_hex, start_hex + NUM_VTX/8 - 1);

  // delete about 1% of vertices
  const int step = 100;
  int remaining = NUM_VTX;
  Range::iterator vit = verts.begin();
  EntityHandle entities[2];
  for (int j = 0; j < remaining; j += step) {
    entities[0] = start_hex + j/8;
    entities[1] = *vit;
    rval = mb.delete_entities( entities, 2);
    CHECK_ERR(rval);
    dead.insert(*vit); dead.insert(start_hex + j/8);
    vit = verts.erase(vit);
    vit += step-1;
    hexes.erase(start_hex + j/8);
  }
  
  // Remove some additional values from the range 
  // so that our handle blocks don't always align with
  // sequences
  verts.erase( verts.begin() + (step-5), verts.begin() + (step+5) );
  hexes.erase( hexes.begin() + (step/8-5), hexes.begin() + (step/8+5) );
  
  // Check that we get back expected values
  double *xcoord, *ycoord, *zcoord;
  vit = verts.begin();
  int count, total = 0;
  while (vit != verts.end()) {
    rval = mb.coords_iterate(vit, verts.end(), xcoord, ycoord, zcoord, count);
    if (MB_SUCCESS && (!xcoord || !ycoord || !zcoord)) rval = MB_FAILURE;
    CHECK_ERR(rval);
    
    assert(total + count <= (int)verts.size());
    for (int i = 0; i < count; i++) {
        // vertex handles start at 1, so need to subtract one
      double val = *vit + (double)i - 1.0;
      CHECK_REAL_EQUAL(val, xcoord[i], 1.0e-10);
      CHECK_REAL_EQUAL(val, ycoord[i], 1.0e-10);
      CHECK_REAL_EQUAL(val, zcoord[i], 1.0e-10);
    }
    
  // Check that we can set values and get the right values back
    for (int i = 0; i < count; i++) {
      xcoord[i] *= 2.0;
      ycoord[i] *= 2.0;
      zcoord[i] *= 2.0;
    }
    
    std::vector<double> dum(3*count);
    Range dum_verts(*vit, *vit + count - 1);
    rval = mb.get_coords(dum_verts, &dum[0]);
    CHECK_ERR(rval);
    for (int i = 0; i < count; i++) {
        // vertex handles start at 1, so need to subtract 1 from expected value
      double val = 2.0*(*vit + (double)i - 1);
      CHECK_REAL_EQUAL(val, xcoord[i], 1.0e-10);
      CHECK_REAL_EQUAL(val, ycoord[i], 1.0e-10);
      CHECK_REAL_EQUAL(val, zcoord[i], 1.0e-10);
    }

    vit += count;
    total += count;  
  }
  
    // now check connectivity
  Range::iterator hit = hexes.begin();
  EntityHandle *connect = NULL;
  EntityHandle dum_connect[8];
  int num_connect;
  while (hit != hexes.end()) {
    rval = mb.connect_iterate(hit, hexes.end(), connect, num_connect, count);
    if (MB_SUCCESS && !connect) rval = MB_FAILURE;
    CHECK_ERR(rval);
    CHECK_EQUAL(num_connect, 8);
    
      // should be equal to initial connectivity
    for (int i = 0; i < count; i++) {
      EntityHandle first = 8*(*hit - start_hex + i) + 1;
      for (unsigned int j = 0; j < 8; j++)
        dum_connect[j] = first + j;
      CHECK_ARRAYS_EQUAL(connect, 8, &dum_connect[0], 8);
      connect += 8;
    }
    
    hit += count;
    total += count;
  }

    // ok, done
}
  
void test_scd_invalid()
{
    // check that we get errors from structured mesh
  Core moab;
  Interface& mb = moab;
  ScdInterface *scdi;
  ErrorCode rval = mb.query_interface(scdi);
  CHECK_ERR(rval);
  
    // make an arbitrary structured mesh
  const int NUM_DIMS = 10;
  HomCoord low(0, 0, 0), high(NUM_DIMS, NUM_DIMS, NUM_DIMS);
  ScdBox *new_box = NULL;
  rval = scdi->construct_box(low, high, NULL, 0, new_box);
  CHECK_ERR(rval);
  CHECK(new_box != NULL);
  
  EntityHandle start_hex = new_box->start_element();
  Range hexes(start_hex, start_hex + NUM_DIMS*NUM_DIMS*NUM_DIMS - 1);
  
    // should be able to get vertices used by this box
  Range verts;
  rval = mb.get_adjacencies(hexes, 0, false, verts, Interface::UNION);
  CHECK_ERR(rval);
  CHECK_EQUAL((int)verts.size(), (int)((NUM_DIMS+1)*(NUM_DIMS+1)*(NUM_DIMS+1)));
  
    // should NOT be able to get connect iterator
  EntityHandle *connect;
  int count, num_connect;
  rval = mb.connect_iterate(hexes.begin(), hexes.end(), connect, num_connect, count);
  CHECK_EQUAL(rval, MB_FAILURE);
}
