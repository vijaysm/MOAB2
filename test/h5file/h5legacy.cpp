#include "moab/Core.hpp"
#include "TestUtil.hpp"
#include "moab/Range.hpp"

#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace moab;

void calc_centroid( Interface* iface, EntityHandle pent, double result[3] )
{
  int len;
  const EntityHandle* conn;
  ErrorCode rval;
  rval = iface->get_connectivity( pent, conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL( 5, len );
  
  double coords[15];
  rval = iface->get_coords( conn, len, coords );
  CHECK_ERR(rval);
  
  for (int d = 0; d < 3; ++d) {
    result[d] = 0;
    for (int i = 0; i < 5; ++i)
      result[d] += coords[3*i+d];
    result[d] /= 5;
  }
}

void test_moab_v3_poly_format()
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
    // load file containing a dodecahedron
  rval = mb.load_mesh( STRINGIFY(SRCDIR) "/v3_dodec.h5m" );
  CHECK_ERR(rval);
  
    // get entities from file
  Range verts, faces, polyhedrons;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type( 0, MBPOLYGON, faces );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type( 0, MBPOLYHEDRON, polyhedrons );
  CHECK_ERR(rval);
  
    // check expected number of entities
  CHECK_EQUAL( (size_t)20, verts.size() );
  CHECK_EQUAL( (size_t)12, faces.size() );
  CHECK_EQUAL( (size_t)1,  polyhedrons.size() );
  const EntityHandle polyhedron = polyhedrons.front();
  
    // check the polyhedron connectivity list
  std::vector<EntityHandle> faces1, faces2;
  std::copy( faces.begin(), faces.end(), std::back_inserter(faces1) );
  rval = mb.get_connectivity( &polyhedron, 1, faces2 );
  CHECK_ERR(rval);
  std::sort( faces2.begin(), faces2.end() );
  CHECK( faces1 == faces2 );
  
    // each polygonshould have a tag value storing its centroid.
    // compare this value against the centroid calculated from
    // the vertex coords.
  
    // get tag for saved values
  Tag centroid;
  rval = mb.tag_get_handle( "CENTROID", 3, MB_TYPE_DOUBLE, centroid );
  CHECK_ERR(rval);
  
    // for each face...
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
    double saved[3], calc[3];
    rval = mb.tag_get_data( centroid, &*i, 1, saved );
    CHECK_ERR(rval);
    calc_centroid( &mb, *i, calc );
    CHECK_REAL_EQUAL( saved[0], calc[0], 1e-6 );
    CHECK_REAL_EQUAL( saved[1], calc[1], 1e-6 );
    CHECK_REAL_EQUAL( saved[2], calc[2], 1e-6 );
  }
}

  


int main()
{
  // only one test so far... should probably add second test
  // for really-old-format  entityset parent/child links
  return RUN_TEST( test_moab_v3_poly_format );
}

  
