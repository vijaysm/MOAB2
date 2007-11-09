#include "MBCore.hpp"
#include "testdir.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>

// define our own assert so it doesn't get compiled out during release builds
#define assert( A ) if (!(A)) { \
  std::cerr << "Assertion failed at line " << __LINE__ << " : " << #A << std::endl; \
  abort(); \
}


void calc_centroid( MBInterface* iface, MBEntityHandle pent, double result[3] )
{
  int len;
  const MBEntityHandle* conn;
  MBErrorCode rval;
  rval = iface->get_connectivity( pent, conn, len );
  assert(MB_SUCCESS == rval && 5 == len);
  
  double coords[15];
  rval = iface->get_coords( conn, len, coords );
  assert(MB_SUCCESS == rval);
  
  for (int d = 0; d < 3; ++d) {
    result[d] = 0;
    for (int i = 0; i < 5; ++i)
      result[d] += coords[3*i+d];
    result[d] /= 5;
  }
}

void test_moab_v3_poly_format()
{
  MBCore moab;
  MBInterface& mb = moab;
  MBErrorCode rval;
  
    // load file containing a dodecahedron
  rval = mb.load_mesh( "../../" TEST_DIR "/v3_dodec.h5m" );
  assert( MB_SUCCESS == rval );
  
    // get entities from file
  MBRange verts, faces, polyhedrons;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  assert( MB_SUCCESS == rval );
  rval = mb.get_entities_by_type( 0, MBPOLYGON, faces );
  assert( MB_SUCCESS == rval );
  rval = mb.get_entities_by_type( 0, MBPOLYHEDRON, polyhedrons );
  assert( MB_SUCCESS == rval );
  
    // check expected number of entities
  assert( 20 == verts.size() );
  assert( 12 == faces.size() );
  assert( 1 == polyhedrons.size() );
  const MBEntityHandle polyhedron = polyhedrons.front();
  
    // check the polyhedron connectivity list
  std::vector<MBEntityHandle> faces1, faces2;
  std::copy( faces.begin(), faces.end(), std::back_inserter(faces1) );
  rval = mb.get_connectivity( &polyhedron, 1, faces2 );
  assert( MB_SUCCESS == rval );
  std::sort( faces2.begin(), faces2.end() );
  assert( faces1 == faces2 );
  
    // each polygonshould have a tag value storing its centroid.
    // compare this value against the centroid calculated from
    // the vertex coords.
  
    // get tag for saved values
  MBTag centroid;
  rval = mb.tag_get_handle( "CENTROID", centroid );
  assert( MB_SUCCESS == rval );
  
    // for each face...
  for (MBRange::iterator i = faces.begin(); i != faces.end(); ++i) {
    double saved[3], calc[3];
    rval = mb.tag_get_data( centroid, &*i, 1, saved );
    assert( MB_SUCCESS == rval );
    calc_centroid( &mb, *i, calc );
    assert( fabs(saved[0] - calc[0]) < 1e-6 );
    assert( fabs(saved[1] - calc[1]) < 1e-6 );
    assert( fabs(saved[2] - calc[2]) < 1e-6 );
  }
}

  


int main()
{
  // only one test so far... should probably add second test
  // for really-old-format  entityset parent/child links
  test_moab_v3_poly_format();
  return 0;
}

  
