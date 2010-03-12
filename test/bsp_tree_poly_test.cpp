#include "moab/BSPTreePoly.hpp"
#include "TestUtil.hpp"
#include "moab/CartVect.hpp"

using namespace moab;

void test_construct_from_hex();
void test_cut_with_plane();
void test_volume();

int main( )
{
  int error_count = 0;
  error_count += RUN_TEST( test_construct_from_hex );
  error_count += RUN_TEST( test_cut_with_plane );
  error_count += RUN_TEST( test_volume );
  return error_count;
}

  
const int hex_faces[6][4] = { { 0, 1, 5, 4 },
                              { 1, 2, 6, 5 },
                              { 2, 3, 7, 6 },
                              { 3, 0, 4, 7 },
                              { 3, 2, 1, 0 },
                              { 4, 5, 6, 7 } };

static void get_corners( CartVect corners[8] )
{
  corners[0] = CartVect( 1, 1, 0 );
  corners[1] = CartVect( 6, 1, 0 );
  corners[2] = CartVect( 6, 3, 0 );
  corners[3] = CartVect( 1, 3, 0 );
  corners[4] = CartVect( 1, 1, 2 );
  corners[5] = CartVect( 4, 1, 2 );
  corners[6] = CartVect( 4, 3, 2 );
  corners[7] = CartVect( 1, 3, 2 );
}

const BSPTreePoly::Face* find_face( const BSPTreePoly& poly,
                                     const CartVect* coords,
                                     int num_corners,
                                     const int* face_indices = 0 )
{
  std::vector<const BSPTreePoly::Face*>::iterator i;
  std::vector<const BSPTreePoly::Face*> faces;
  std::vector<CartVect> corners;
  poly.get_faces( faces );
  for (i = faces.begin(); i != faces.end(); ++i) {
    corners.clear();
    poly.get_vertices( *i, corners );
    if (corners.size() != (unsigned)num_corners)
      continue;    
    
    int j;
    for (j = 0; j < num_corners; ++j) {
      int corner = face_indices ? face_indices[j] : j;
      if ((coords[corner] - corners.front()).length_squared() < 1e-12)
        break;
    }
    if (j == num_corners)
      continue;
    
    int k;
    for (k = 1; k < num_corners; ++k) {
      int corner = face_indices ? face_indices[(j+k)%num_corners] : (j+k)%num_corners;
      if ((coords[corner] - corners[k]).length_squared() > 1e-12)
        break;
    }
    
    if (k == num_corners) 
      return *i;
  }    
  return 0;
}
                                     

void test_construct_from_hex()
{
  BSPTreePoly::reset_debug_ids();

  CartVect corners[8];
  get_corners( corners );
  BSPTreePoly poly( corners );
  CHECK( poly.is_valid() );
  
  std::vector<const BSPTreePoly::Face*> faces;
  poly.get_faces( faces );
  CHECK_EQUAL( (size_t)6, faces.size() );
  
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[0] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[1] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[2] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[3] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[4] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[5] ) );
}

void test_cut_with_plane()
{
    // create a hexahedron
  BSPTreePoly::reset_debug_ids();
  CartVect corners[8];
  get_corners( corners );
  BSPTreePoly poly( corners );
  CHECK( poly.is_valid() );
  
    // check that a plane entirely above the 
    // polyhedron (coincident with one face)
    // doesn't modify the polyhedron
  bool r = poly.cut_polyhedron( CartVect(0,0,1), -2 );
  CHECK(!r);
  CHECK( poly.is_valid() );
  
    // cut in half with Z=1 plane
  r = poly.cut_polyhedron( CartVect(0,0,1), -1 );
  CHECK(r);
  CHECK( poly.is_valid() );
  for (int i = 0; i < 8; ++i) {
    if (fabs(corners[i][2] - 2) < 1e-6)
      corners[i][2] = 1;
    if (fabs(corners[i][0] - 4) < 1e-6)
      corners[i][0] = 5;
  }
  
  std::vector<const BSPTreePoly::Face*> faces;
  poly.get_faces( faces );
  CHECK_EQUAL( (size_t)6, faces.size() );
  
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[0] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[1] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[2] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[3] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[4] ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[5] ) );
  
    // create a hexahedron
  BSPTreePoly::reset_debug_ids();
  get_corners( corners );
  poly.set( corners );
  CHECK( poly.is_valid() );
  
    // cut off two corners using X=5 plane
  r = poly.cut_polyhedron( CartVect(1,0,0), -5 );
  CHECK(r);
  CHECK( poly.is_valid() );
  
  faces.clear();
  poly.get_faces( faces );
  CHECK_EQUAL( (size_t)7, faces.size() );
  
  CartVect new_vtx1( 5, 1, 1 );
  CartVect new_vtx2( 5, 1, 0 );
  CartVect new_vtx3( 5, 3, 0 );
  CartVect new_vtx4( 5, 3, 1 );
  
  CartVect face1[5] = { corners[0], new_vtx2, new_vtx1, corners[5], corners[4] };
  CartVect face2[4] = { new_vtx1, new_vtx4, corners[6], corners[5] };
  CartVect face3[5] = { new_vtx4, new_vtx3, corners[3], corners[7], corners[6] };
  CartVect face5[4] = { corners[3], new_vtx3, new_vtx2, corners[0] };
  CartVect face7[4] = { new_vtx1, new_vtx2, new_vtx3, new_vtx4 };
  
  CHECK( 0 != find_face( poly, face1, 5 ) );
  CHECK( 0 != find_face( poly, face2, 4 ) );
  CHECK( 0 != find_face( poly, face3, 5 ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[3] ) );
  CHECK( 0 != find_face( poly, face5, 4 ) );
  CHECK( 0 != find_face( poly, corners, 4, hex_faces[5] ) );
  CHECK( 0 != find_face( poly, face7, 4 ) );
}

void test_volume()
{
  CartVect corners[8];
  get_corners( corners );
  BSPTreePoly poly( corners );
  CHECK( poly.is_valid() );
  
  CHECK_REAL_EQUAL( 16.0, poly.volume(), 1e-6 );
}

