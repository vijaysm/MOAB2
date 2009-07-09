#include "MBCore.hpp"
#include "MBAdaptiveKDTree.hpp"
#include "MBCartVect.hpp"
#include "MBGeomUtil.hpp"
#include "MBRange.hpp"

#include <iostream>
#include <algorithm>

/* Like assert, except returns false rather than call abort() */
#define CHECK( A ) \
  if (!(A)) return check_failure( #A, __LINE__ )

/* used in CHECK macro */
bool check_failure( const char* str, int line )
{
    std::cout << "FAILURE:" << std::endl 
              << "  Condtion: " << str << std::endl 
              << "  Line: " << line << std::endl; 
    return false; 
}

/* Utility method - compare two range boxes */
bool box_equal( const MBAdaptiveKDTreeIter& iter, 
                double x_min, double y_min, double z_min,
                double x_max, double y_max, double z_max )
{
  return iter.box_min()[0] == x_min &&
         iter.box_min()[1] == y_min &&
         iter.box_min()[2] == z_min &&
         iter.box_max()[0] == x_max &&
         iter.box_max()[1] == y_max &&
         iter.box_max()[2] == z_max;
}

void build_triangles( MBInterface* moab, MBRange& tris,
                      int num_vert, const double* coords,
                      int num_tri,  const unsigned* conn )
{
  std::vector<MBEntityHandle> verts(num_vert);
  for (int i = 0; i < num_vert; ++i)
    moab->create_vertex( coords + 3*i, verts[i] );
  
  MBEntityHandle tri, tri_verts[3];
  for (int i = 0; i < num_tri; ++i) {
    tri_verts[0] = verts[conn[3*i  ]];
    tri_verts[1] = verts[conn[3*i+1]];
    tri_verts[2] = verts[conn[3*i+2]];
    moab->create_element( MBTRI, tri_verts, 3, tri );
    tris.insert(tri);
  }

}  

/* Utility method - build a 2x2x2 box composed of four triagles on a side */
void build_triangle_box_small( MBInterface* moab, MBRange& tris )
{
  const double coords[] = { -1, -1, -1,
                             1, -1, -1,
                             1,  1, -1,
                            -1,  1, -1,
                            -1, -1,  1,
                             1, -1,  1,
                             1,  1,  1,
                            -1,  1,  1 };
  
  const unsigned conn[] = { 0, 1, 5,   0, 5, 4,
                            2, 6, 7,   2, 7, 3,
                            1, 2, 6,   1, 6, 5,
                            0, 3, 4,   3, 7, 4,
                            0, 1, 3,   1, 2, 3,
                            4, 5, 6,   4, 6, 7 };
  
  build_triangles( moab, tris, 8, coords, 12, conn );
}

/* build 6x6x6 box centered at origin and composed of 18 triangles per side,
   Each side of the box is composed of a 3x3 grid of quads, where each
   quad is split diagonally to form two triangles. */
void build_triangle_box_large( MBInterface* moab, MBRange& tris )
{
  const double coords[] = { // corners
                            -3, -3, -3,
                             3, -3, -3,
                             3,  3, -3,
                            -3,  3, -3,
                            -3, -3,  3,
                             3, -3,  3,
                             3,  3,  3,
                            -3,  3,  3,
                            
                            // edges
                            -1, -3, -3,
                             1, -3, -3,
                             3, -1, -3,
                             3,  1, -3,
                             1,  3, -3,
                            -1,  3, -3,
                            -3,  1, -3,
                            -3, -1, -3,

                            -1, -3,  3,
                             1, -3,  3,
                             3, -1,  3,
                             3,  1,  3,
                             1,  3,  3,
                            -1,  3,  3,
                            -3,  1,  3,
                            -3, -1,  3,
                            
                            -3, -3, -1,
                            -3, -3,  1,
                             3, -3, -1,
                             3, -3,  1,
                             3,  3, -1,
                             3,  3,  1,
                            -3,  3, -1,
                            -3,  3,  1,
                            
                            // faces
                            -1, -3, -1,
                             1, -3, -1,
                             1, -3,  1,
                            -1, -3,  1,
                            
                             3, -1, -1,
                             3,  1, -1,
                             3,  1,  1,
                             3, -1,  1,
                             
                             1,  3, -1,
                            -1,  3, -1,
                            -1,  3,  1,
                             1,  3,  1,
                             
                            -3,  1, -1,
                            -3, -1, -1,
                            -3, -1,  1,
                            -3,  1,  1,
                            
                            -1, -1, -3,
                             1, -1, -3,
                             1,  1, -3,
                            -1,  1, -3,
                            
                            -1, -1,  3,
                             1, -1,  3,
                             1,  1,  3,
                            -1,  1,  3 };
 const unsigned conn[] = { // face 0
                            0,  8, 24,
                            8, 32, 24,
                            8, 33, 32,
                            8,  9, 33,
                            9,  1, 33,
                            1, 26, 33,
                           24, 35, 25, 
                           24, 32, 35, 
                           32, 33, 35,
                           33, 34, 35,
                           33, 27, 34,
                           33, 26, 27,
                           35,  4, 25,
                           35, 16,  4,
                           35, 17, 16,
                           35, 34, 17,
                           27, 17, 34,
                           27,  5, 17,
                           // face 1
                           36, 26,  1,
                           36,  1, 10,
                           36, 10, 11,
                           36, 11, 37,
                           11, 28, 37,
                           11,  2, 28,
                           36, 27, 26,
                           36, 39, 27,
                           36, 38, 39,
                           36, 37, 38,
                           37, 28, 38,
                           28, 29, 38,
                           18,  5, 27,
                           18, 27, 39,
                           18, 39, 38,
                           18, 38, 19,
                            6, 19, 38,
                            6, 38, 29,
                           // face 2
                           12, 28,  2,
                           12, 40, 28,
                           12, 41, 40,
                           12, 13, 41,
                            3, 41, 13,
                            3, 30, 41,
                           43, 29, 28,
                           43, 28, 40,
                           43, 40, 41,
                           43, 41, 42,
                           41, 31, 42,
                           41, 30, 31,
                           43,  6, 29,
                           43, 20,  6,
                           43, 21, 20,
                           43, 42, 21,
                           21, 42, 31,
                           21, 31,  7,
                           // face 3
                           44, 30,  3,
                           44,  3, 14,
                           44, 14, 15,
                           44, 15, 45,
                           15, 24, 45,
                           15,  0, 24,
                           44, 31, 30,
                           44, 47, 31,
                           44, 46, 47,
                           44, 45, 46,
                           46, 45, 24,
                           46, 24, 25,
                           31, 22,  7,
                           31, 47, 22,
                           46, 22, 47,
                           46, 23, 22,
                           46,  4, 23,
                           46, 25,  4,
                           // face 4
                            8, 15,  0,
                            8, 48, 15,
                            8, 49, 48,
                            8,  9, 49,
                            1, 49,  9,
                            1, 10, 49,
                           51, 14, 15,
                           51, 15, 48,
                           51, 48, 49,
                           51, 49, 50,
                           11, 50, 49,
                           11, 49, 10,
                           51,  3, 14,
                           51, 13,  3,
                           51, 12, 13,
                           51, 50, 12,
                           11, 12, 50,
                           11,  2, 12,
                           // face 5
                            4, 52, 16,
                            4, 23, 52,
                           22, 52, 23,
                           22, 55, 52,
                           22, 21, 55,
                           22,  7, 21,
                           17, 16, 52,
                           17, 52, 53,
                           54, 53, 52,
                           54, 52, 55,
                           54, 55, 21,
                           54, 21, 20,
                           18,  5, 17,
                           18, 17, 53, 
                           18, 53, 54,
                           18, 54, 19,
                            6, 19, 54,
                            6, 54, 20 };
                            
  build_triangles( moab, tris, 56, coords, 108, conn );
}
                            
                            
                            

/* Utility method - build 2x2x2 octahedron (3D diamond)*/
void build_triangle_octahedron( MBInterface* moab, MBRange& tris )
{
  const double coords[] = {  1,  0,  0,
                             0,  1,  0,
                            -1,  0,  0,
                             0, -1,  0,
                             0,  0,  1,
                             0,  0, -1 };
  const unsigned conn[] = { 0, 1, 4,
                            1, 2, 4,
                            2, 3, 4,
                            3, 0, 4,
                            1, 0, 5,
                            2, 1, 5,
                            3, 2, 5,
                            0, 3, 5 };

  build_triangles( moab, tris, 6, coords, 8, conn );
}

bool test_valid_tree( MBAdaptiveKDTree* tool, MBEntityHandle root, 
                      MBAdaptiveKDTree::Settings& settings,
                      const MBRange& expected_tris )
{
  MBRange all_tris;
  MBErrorCode rval;
  MBAdaptiveKDTreeIter iter;
  CHECK(MB_SUCCESS == tool->get_tree_iterator( root, iter ));
  do {
    CHECK( !settings.maxTreeDepth || iter.depth() <= settings.maxTreeDepth );
    
    MBRange tris;
    CHECK( tool->moab()->get_entities_by_type( iter.handle(), MBTRI, tris ) == MB_SUCCESS );
    //CHECK( !tris.empty() );
    all_tris.merge( tris );

    const MBCartVect min(iter.box_min()), max(iter.box_max());
    const MBCartVect cen(0.5*(min+max)), hdim(0.5*(max-min));
    
    for (MBRange::iterator i = tris.begin(); i != tris.end(); ++i) {
      MBEntityHandle tri = *i;
      const MBEntityHandle* conn;
      int conn_len;
      CHECK( tool->moab()->get_connectivity( tri, conn, conn_len ) == MB_SUCCESS );
      CHECK( conn_len == 3 );
      MBCartVect coords[3];
      CHECK( tool->moab()->get_coords( conn, 3, coords[0].array() ) == MB_SUCCESS );
      CHECK( MBGeomUtil::box_tri_overlap( coords, cen, hdim ) );
    }
  } while (MB_SUCCESS == (rval = iter.step()));  
  
  CHECK( MB_ENTITY_NOT_FOUND == rval );
  CHECK( expected_tris == all_tris );
  return true;
}

/* utility method - check that all tris share a vertex */
bool check_common_vertex( MBInterface* moab,
                          const MBEntityHandle* tris,
                          unsigned num_tri,
                          MBCartVect point )
{
  for (unsigned i = 0; i < num_tri; ++i)
    CHECK( MBTRI == moab->type_from_handle( tris[i] ) );
  
  MBErrorCode rval;
  MBCartVect tri_coords[3];
  const MBEntityHandle* conn;
  int conn_len;
  rval = moab->get_connectivity( tris[0], conn, conn_len );
  CHECK( MB_SUCCESS == rval );
  rval = moab->get_coords( conn, 3, tri_coords[0].array() );
  CHECK( MB_SUCCESS == rval );
  tri_coords[0] -= point;
  tri_coords[1] -= point;
  tri_coords[2] -= point;
  MBEntityHandle vertex = 0;
  if (tri_coords[0] % tri_coords[0] < 1e-6)
    vertex = conn[0];
  else if (tri_coords[1] % tri_coords[1] < 1e-6)
    vertex = conn[1];
  else if (tri_coords[2] % tri_coords[2] < 1e-6)
    vertex = conn[2];
  else CHECK(false);
  for (unsigned j = 1; j < num_tri; ++j) {
    rval = moab->get_connectivity( tris[j], conn, conn_len );
    CHECK( conn[0] == vertex || conn[1] == vertex || conn[2] == vertex );
  }
  return true;
}                        
 
/* utility method - check that all tris share a vertex */
bool check_point_in_triangles( MBInterface* moab,
                               const MBEntityHandle* tris,
                               unsigned num_tris,
                               MBCartVect point )
{
  MBErrorCode rval;
  MBCartVect tri_coords[3], tript;
  const MBEntityHandle* conn;
  int conn_len;

  for (unsigned i = 0; i < num_tris; ++i) {
    CHECK( MBTRI == moab->type_from_handle( tris[i] ) );

    rval = moab->get_connectivity( tris[i], conn, conn_len );
    CHECK( MB_SUCCESS == rval );
    rval = moab->get_coords( conn, 3, tri_coords[0].array() );
    CHECK( MB_SUCCESS == rval );
    
    MBGeomUtil::closest_location_on_tri( point, tri_coords, tript );
    tript -= point;
    CHECK( fabs(tript[0]) < 1e-6 );
    CHECK( fabs(tript[1]) < 1e-6 );
    CHECK( fabs(tript[2]) < 1e-6 );
  }
  return true;
}

  /* Create the following 2D tree (no Z splits)
  
          (6) X = -3    (5) X = 0
          /             |
          /             |      [8]
          /             |---------------- (7) Y = 3
      [5] /    [6]      |
          /             |      [7]
          /             |
          /             |
    -------------------------------------- (1) Y = 0
                    |                |
          [2]       |                |
 (3)----------------|                |
 Y = -2             |       [3]      |   [4]
                    |                |
          [1]       |                |
                    |                |
                    (2) X = -1      (4) X = 4
   */
bool create_simple_2d_tree( MBAdaptiveKDTree& tool, 
                            MBEntityHandle& root,
                            MBEntityHandle leaves[9] )
{
  MBErrorCode rval;
  MBAdaptiveKDTree::Plane plane;

    // create a single-node tree
  double min[3] = { -5, -4, -1 };
  double max[3] = {  5,  4,  1 };
  rval = tool.create_tree( min, max, root );
  CHECK( MB_SUCCESS == rval );

    // get iterator for tree
  MBAdaptiveKDTreeIter iter;
  rval = tool.get_tree_iterator( root, iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal(iter, -5, -4, -1, 5, 4, 1 ) );
   
   // split plane (1)
  plane.norm = MBAdaptiveKDTree::Y;
  plane.coord = 0.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, 5, 0, 1 ) );
  
   // split plane (2)
  plane.norm = MBAdaptiveKDTree::X;
  plane.coord = -1.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, 0, 1 ) );
  
   // split plane (3), leaf [1]
  plane.norm = MBAdaptiveKDTree::Y;
  plane.coord = -2.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, -2, 1 ) );
  leaves[1] = iter.handle();
  
    // leaf [2]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -2, -1, -1, 0, 1 ) );
  leaves[2] = iter.handle();
  
    // non-leaf right of split plane (2)
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -1, -4, -1, 5, 0, 1 ) );
  
    // split plane (4) and leaf [3]
  plane.norm = MBAdaptiveKDTree::X;
  plane.coord = 4.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -1, -4, -1, 4, 0, 1 ) );
  leaves[3] = iter.handle();
  
    // leaf [4]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 4, -4, -1, 5, 0, 1 ) );
  leaves[4] = iter.handle();
  
    // non-leaf above split plane (1)
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, 5, 4, 1 ) );
  
    // split plane (5)
  plane.norm = MBAdaptiveKDTree::X;
  plane.coord = 0.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, 0, 4, 1 ) );
  
    // split plane (6) and leaf [5]
  plane.norm = MBAdaptiveKDTree::X;
  plane.coord = -3.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, -3, 4, 1 ) );
  leaves[5] = iter.handle();
  
    // leaf [6];
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -3, 0, -1, 0, 4, 1 ) );
  leaves[6] = iter.handle();
  
    // non-leaf right of split plane (5)
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 0, -1, 5, 4, 1 ) );
  
    // split plane (7) and leaf [7]
  plane.norm = MBAdaptiveKDTree::Y;
  plane.coord = 3.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 0, -1, 5, 3, 1 ) );
  leaves[7] = iter.handle();
  
    // leaf [8];
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 3, -1, 5, 4, 1 ) );
  leaves[8] = iter.handle();
  
    // the end
  rval = iter.step();
  CHECK( MB_ENTITY_NOT_FOUND == rval );

  return true;
}

bool leaf_iterator_test()
{
  MBErrorCode rval;
  MBCore moab;
  MBInterface* mb = &moab;
  MBAdaptiveKDTree tool(mb);


    // create a single-node tree
  MBEntityHandle root;
  double min[3] = { -3, -2, -1 };
  double max[3] = {  1,  2,  3 };
  rval = tool.create_tree( min, max, root );
  CHECK( MB_SUCCESS == rval );

    // get iterator for tree
  MBAdaptiveKDTreeIter iter;
  rval = tool.get_tree_iterator( root, iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal(iter, -3, -2, -1, 1, 2, 3 ) );

    // check that steping the iterator fails correctly.
  rval = iter.step();
  CHECK( MB_ENTITY_NOT_FOUND == rval );
  rval = iter.step();
  CHECK( MB_SUCCESS != rval );

  rval = tool.delete_tree( root );
  CHECK( MB_SUCCESS == rval );
  root = 0;
  
  
    // construct a simple 2D tree for remaining tests
  MBEntityHandle leaves[9];
  CHECK( create_simple_2d_tree( tool, root, leaves ) );

  
  /**** Now traverse tree again, and check neighbors of each leaf *****/
  std::vector<MBAdaptiveKDTreeIter> list;
  
    // leaf [1]
  
  rval = tool.get_tree_iterator( root, iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, -2, 1 ) );
  CHECK( iter.handle() == leaves[1] );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::X, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::X, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 1 );
    // should be leaf [3]
  CHECK( box_equal( list.front(), -1, -4, -1, 4, 0, 1 ) );
  CHECK( list.front().handle() == leaves[3] );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Y, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Y, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 1 );
    // should be leaf [2]
  CHECK( box_equal( list.front(), -5, -2, -1, -1, 0, 1 ) );
  CHECK( list.front().handle() == leaves[2] );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Z, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Z, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
    // leaf [2]
  
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -2, -1, -1, 0, 1 ) );
  CHECK( iter.handle() == leaves[2] );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::X, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::X, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 1 );
    // should be leaf [3]
  CHECK( box_equal( list.front(), -1, -4, -1, 4, 0, 1 ) );
  CHECK( list.front().handle() == leaves[3] );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Y, true, list );
  CHECK( MB_SUCCESS == rval );
    // should be leaf [1]
  CHECK( list.size() == 1 );
  CHECK( box_equal( list.front(), -5, -4, -1, -1, -2, 1 ) );
  CHECK( list.front().handle() == leaves[1] );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Y, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 2 );
    // should be leaf [5] and leaf[6]
  if (list[0].handle() == leaves[6])
    std::swap( list[0], list[1] );
  CHECK( list[0].handle() == leaves[5] );
  CHECK( list[1].handle() == leaves[6] );
  CHECK( box_equal( list[0], -5, 0, -1, -3, 4, 1 ) );
  CHECK( box_equal( list[1], -3, 0, -1, 0, 4, 1 ) );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Z, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( MBAdaptiveKDTree::Z, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
    // leaf [3]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -1, -4, -1, 4, 0, 1 ) );
  CHECK( iter.handle() == leaves[3] );
  
    // leaf [4]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 4, -4, -1, 5, 0, 1 ) );
  CHECK( iter.handle() == leaves[4] );
  
    // leaf [5]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, -3, 4, 1 ) );
  CHECK( iter.handle() == leaves[5] );
  
    // leaf [6]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -3, 0, -1, 0, 4, 1 ) );
  CHECK( iter.handle() == leaves[6] );
  
    // leaf [7]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 0, -1, 5, 3, 1 ) );
  CHECK( iter.handle() == leaves[7] );
  
    // leaf [8]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 3, -1, 5, 4, 1 ) );
  CHECK( iter.handle() == leaves[8] );
  
  return true;
}

bool test_tree_merge_nodes()
{
    // build simple tree for tests
  MBErrorCode rval;
  MBCore moab;
  MBInterface* mb = &moab;
  MBAdaptiveKDTree tool(mb);
  MBEntityHandle root;
  MBEntityHandle leaves[9];
  CHECK( create_simple_2d_tree( tool, root, leaves ) );

    // get iterator for tree
  MBAdaptiveKDTreeIter iter;
  rval = tool.get_tree_iterator( root, iter );
  CHECK( MB_SUCCESS == rval );

    // merge leaves 1 and 2
  rval = tool.merge_leaf( iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, 0, 1 ) );
  
    // merge leaf 1,2 with 3,4 (implicity merges 3 and 4)
  rval = tool.merge_leaf( iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, 5, 0, 1 ) );
  
    // make sure iterator remains valid
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
    // leaf 5
  CHECK( box_equal( iter, -5, 0, -1, -3, 4, 1 ) );
  CHECK( iter.handle() == leaves[5] );
  
  return true;  
}

bool test_build_tree_bisect_triangles( )
{
  MBCore moab;
  MBAdaptiveKDTree tool( &moab );
  MBRange box_tris;
  
  MBAdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION;
  MBEntityHandle root;
  MBErrorCode rval;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_octahedron( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );
  
  settings.maxTreeDepth = 2;
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );
  
  return true;
}
  
bool test_closest_triangle()
{
  MBCore moab;
  MBAdaptiveKDTree tool( &moab );
  MBRange box_tris;
  
  MBAdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION;
  MBEntityHandle root;
  MBErrorCode rval;
  MBEntityHandle tri;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );

    // test closest to each corner of the box.
  for (unsigned i = 0; i < 8; ++i) {
    int x = i & 1 ? 1 : -1;
    int y = i & 2 ? 1 : -1;
    int z = i & 4 ? 1 : -1;
    double coords[] = { 2.0*x, 2.0*y, 2.0*z };
    MBCartVect point;
    rval = tool.closest_triangle( root, coords, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    CHECK( fabs(x - point[0]) < 1e-6 );
    CHECK( fabs(y - point[1]) < 1e-6 );
    CHECK( fabs(z - point[2]) < 1e-6 );
      // check point is at triangle vertex
    CHECK( check_common_vertex( tool.moab(), &tri, 1, point ) );
  }
  
    // test location on each face
  const MBCartVect facepts[] = { MBCartVect( -1.0, 0.5, 0.5 ),
                                 MBCartVect(  1.0, 0.5, 0.5 ),
                                 MBCartVect( 0.5, -1.0, 0.5 ),
                                 MBCartVect( 0.5,  1.0, 0.5 ),
                                 MBCartVect( 0.5, 0.5, -1.0 ),
                                 MBCartVect( 0.5, 0.5,  1.0 ) };
  for (unsigned i = 0; i < 6; ++i) {
    MBCartVect point;
    rval = tool.closest_triangle( root, facepts[i].array(), point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    const MBCartVect diff = facepts[i] - point;
    CHECK( fabs(diff[0]) < 1e-6 );
    CHECK( fabs(diff[1]) < 1e-6 );
    CHECK( fabs(diff[2]) < 1e-6 );
       // check that point is contained within triangle
    CHECK( check_point_in_triangles( tool.moab(), &tri, 1, point ) );
  }
  
    // test a location really far from the tree
  {
    const double far[] = { 0.75, 0.75, 200 };
    MBCartVect point;
    rval = tool.closest_triangle( root, far, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
    CHECK( fabs(point[0] - 0.75) < 1e-6 );    
    CHECK( fabs(point[1] - 0.75) < 1e-6 );    
    CHECK( fabs(point[2] - 1.00) < 1e-6 );    
       // check that point is contained within triangle
    CHECK( check_point_in_triangles( tool.moab(), &tri, 1, point ) );
  }
  
    // now do it all again with a lot more triangles
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );

    // test closest to each corner of the box.
  for (unsigned i = 0; i < 8; ++i) {
    int x = i & 1 ? 1 : -1;
    int y = i & 2 ? 1 : -1;
    int z = i & 4 ? 1 : -1;
    double coords[] = { 4.0*x, 4.0*y, 4.0*z };
    MBCartVect point;
    rval = tool.closest_triangle( root, coords, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    CHECK( fabs(3.0*x - point[0]) < 1e-6 );
    CHECK( fabs(3.0*y - point[1]) < 1e-6 );
    CHECK( fabs(3.0*z - point[2]) < 1e-6 );
      // check point is at triangle vertex
    CHECK( check_common_vertex( tool.moab(), &tri, 1, point ) );
  }
  
    // test location on each face
  const MBCartVect facepts2[] = { MBCartVect( -3.0, 0.5, 0.5 ),
                                 MBCartVect(  3.0, 0.5, 0.5 ),
                                 MBCartVect( 0.5, -3.0, 0.5 ),
                                 MBCartVect( 0.5,  3.0, 0.5 ),
                                 MBCartVect( 0.5, 0.5, -3.0 ),
                                 MBCartVect( 0.5, 0.5,  3.0 ) };
  for (unsigned i = 0; i < 6; ++i) {
    MBCartVect point;
    rval = tool.closest_triangle( root, facepts2[i].array(), point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    const MBCartVect diff = facepts2[i] - point;
    CHECK( fabs(diff[0]) < 1e-6 );
    CHECK( fabs(diff[1]) < 1e-6 );
    CHECK( fabs(diff[2]) < 1e-6 );
       // check that point is contained within triangle
    CHECK( check_point_in_triangles( tool.moab(), &tri, 1, point ) );
  }
  
    // test a location really far from the tree
  {
    const double far[] = { 2.75, 2.75, 200 };
    MBCartVect point;
    rval = tool.closest_triangle( root, far, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
    CHECK( fabs(point[0] - 2.75) < 1e-6 );    
    CHECK( fabs(point[1] - 2.75) < 1e-6 );    
    CHECK( fabs(point[2] - 3.00) < 1e-6 );    
       // check that point is contained within triangle
    CHECK( check_point_in_triangles( tool.moab(), &tri, 1, point ) );
  }
  
  return true;
}

bool test_sphere_intersect_triangles()
{
  MBCore moab;
  MBAdaptiveKDTree tool( &moab );
  MBRange box_tris;
  
  MBAdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION;
  MBEntityHandle root;
  MBErrorCode rval;
  std::vector<MBEntityHandle> triangles;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );

    // test closest to each corner of the box.
  for (unsigned i = 0; i < 8; ++i) {
    int x = i & 1 ? 1 : -1;
    int y = i & 2 ? 1 : -1;
    int z = i & 4 ? 1 : -1;
    double center[] = { x, y, z };
    triangles.clear();
    rval = tool.sphere_intersect_triangles( root, center, 0.26, triangles );
    CHECK( MB_SUCCESS == rval );
      // expect 3 to 6 triangles, depending on the corner
    CHECK( triangles.size() >= 3 );
    CHECK( triangles.size() <= 6 );
      // check point is at the same vertex for each triangle
    CHECK( check_common_vertex( tool.moab(), &triangles[0], triangles.size(), MBCartVect(center) ) );
  }
  
    // now do it all again with a lot more triangles

  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );

    // test closest to each corner of the box.
  for (unsigned i = 0; i < 8; ++i) {
    int x = i & 1 ? 1 : -1;
    int y = i & 2 ? 1 : -1;
    int z = i & 4 ? 1 : -1;
    double center[] = { 3.0*x, 3.0*y, 3.0*z };
    triangles.clear();
    rval = tool.sphere_intersect_triangles( root, center, 0.26, triangles );
    CHECK( MB_SUCCESS == rval );
      // expect 3 to 6 triangles, depending on the corner
    CHECK( triangles.size() >= 3 );
    CHECK( triangles.size() <= 6 );
      // check point is at the same vertex for each triangle
    CHECK( check_common_vertex( tool.moab(), &triangles[0], triangles.size(), MBCartVect(center) ) );
  }
  
  return true;
}

bool test_ray_intersect_triangles()
{
  MBCore moab;
  MBAdaptiveKDTree tool( &moab );
  MBRange box_tris;
  
  MBAdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = MBAdaptiveKDTree::SUBDIVISION;
  MBEntityHandle root;
  MBErrorCode rval;
  std::vector<MBEntityHandle> tris;
  std::vector<double> dists;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  CHECK( test_valid_tree( &tool, root, settings, box_tris ) );
  
    // test ray through box parallel to X axis
  MBCartVect dir( 1, 0, 0 );
  MBCartVect pt( -2, 0.75, 0.75 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 3 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 || fabs(dists[i] - 3) < 1e-6 );
    MBCartVect tript = pt + dists[i] * dir;
    CHECK( check_point_in_triangles( &moab, &tris[i], 1, tript ) );
  }
   
    // test ray through box parallel to X axis, closest tri only
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists, 1 );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 1 );
  CHECK( dists.size() == tris.size() );
  CHECK( fabs(dists[0] - 1) < 1e-6  );
  CHECK( check_point_in_triangles( &moab, &tris[0], 1, pt + dists[0] * dir ) );

    // test ray ending within box, parallel to X axis
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists, 0, 2.0 );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 1 );
  CHECK( dists.size() == tris.size() );
  CHECK( fabs(dists[0] - 1) < 1e-6  );
  CHECK( check_point_in_triangles( &moab, &tris[0], 1, pt + dists[0] * dir ) );
  
    // test ray starting within box parallel to X axis
  tris.clear(); dists.clear();
  pt = MBCartVect( 0, .75, .75 );
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 2 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 );
    MBCartVect tript = pt + dists[i] * dir;
    CHECK( check_point_in_triangles( &moab, &tris[i], 1, tript ) );
  }
  
    // test skew ray through box
  dir = MBCartVect( 0.5*sqrt(2.0), 0.5*sqrt(2.0), 0.0 );
  pt = MBCartVect( 0, -1.5, 0 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 2 );
  CHECK( dists.size() == tris.size() );
  if (dists[0] < dists[1]) {
    CHECK( fabs(dists[0] - 0.5*sqrt(2.0)) < 1e-6 );
    CHECK( fabs(dists[1] -     sqrt(2.0)) < 1e-6 );
  }
  CHECK( check_point_in_triangles( &moab, &tris[0], 1, pt + dists[0] * dir ) );
  CHECK( check_point_in_triangles( &moab, &tris[1], 1, pt + dists[1] * dir ) );
  
  
    // test ray through box parallel to -Y axis
  dir = MBCartVect( 0, -1, 0 );
  pt = MBCartVect( 0, 2, 0 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 4 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 || fabs(dists[i] - 3) < 1e-6 );
    MBCartVect tript = pt + dists[i] * dir;
    CHECK( check_point_in_triangles( &moab, &tris[i], 1, tript ) );
  }
  
    // test ray through box parallel to Z axis
  dir = MBCartVect( 0, 0, 1 );
  pt = MBCartVect( 0, 0, -2 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 4 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 || fabs(dists[i] - 3) < 1e-6 );
    MBCartVect tript = pt + dists[i] * dir;
    CHECK( check_point_in_triangles( &moab, &tris[i], 1, tript ) );
  }
  
    // test ray through box parallel to Z axis, limit 2 first 2 intersections
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists, 2 );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 2 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 );
    MBCartVect tript = pt + dists[i] * dir;
    CHECK( check_point_in_triangles( &moab, &tris[i], 1, tript ) );
  }
  
  return true;
}
  
  
int main()
{
  int error_count = 0;
  
  error_count += !leaf_iterator_test();
  error_count += !test_build_tree_bisect_triangles();
  error_count += !test_closest_triangle();
  error_count += !test_sphere_intersect_triangles();
  error_count += !test_ray_intersect_triangles();
  error_count += !test_tree_merge_nodes();
  return error_count;
}
