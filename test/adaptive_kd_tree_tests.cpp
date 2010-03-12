#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/Range.hpp"
#include "TestUtil.hpp"

using namespace moab;

#include <iostream>
#include <algorithm>

/* Utility method - compare two range boxes */
bool box_equal( const AdaptiveKDTreeIter& iter, 
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

void build_triangles( Interface* moab, Range& tris,
                      int num_vert, const double* coords,
                      int num_tri,  const unsigned* conn )
{
  std::vector<EntityHandle> verts(num_vert);
  for (int i = 0; i < num_vert; ++i)
    moab->create_vertex( coords + 3*i, verts[i] );
  
  EntityHandle tri, tri_verts[3];
  for (int i = 0; i < num_tri; ++i) {
    tri_verts[0] = verts[conn[3*i  ]];
    tri_verts[1] = verts[conn[3*i+1]];
    tri_verts[2] = verts[conn[3*i+2]];
    moab->create_element( MBTRI, tri_verts, 3, tri );
    tris.insert(tri);
  }

}  

/* Utility method - build a 2x2x2 box composed of four triagles on a side */
void build_triangle_box_small( Interface* moab, Range& tris )
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
void build_triangle_box_large( Interface* moab, Range& tris )
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
void build_triangle_octahedron( Interface* moab, Range& tris )
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

void test_valid_tree( AdaptiveKDTree* tool, EntityHandle root, 
                      AdaptiveKDTree::Settings& settings,
                      const Range& expected_tris )
{
  Range all_tris;
  ErrorCode rval;
  AdaptiveKDTreeIter iter;
  CHECK(MB_SUCCESS == tool->get_tree_iterator( root, iter ));
  do {
    CHECK( !settings.maxTreeDepth || iter.depth() <= settings.maxTreeDepth );
    
    Range tris;
    CHECK( tool->moab()->get_entities_by_type( iter.handle(), MBTRI, tris ) == MB_SUCCESS );
    //CHECK( !tris.empty() );
    all_tris.merge( tris );

    const CartVect min(iter.box_min()), max(iter.box_max());
    const CartVect cen(0.5*(min+max)), hdim(0.5*(max-min));
    
    for (Range::iterator i = tris.begin(); i != tris.end(); ++i) {
      EntityHandle tri = *i;
      const EntityHandle* conn;
      int conn_len;
      CHECK( tool->moab()->get_connectivity( tri, conn, conn_len ) == MB_SUCCESS );
      CHECK( conn_len == 3 );
      CartVect coords[3];
      CHECK( tool->moab()->get_coords( conn, 3, coords[0].array() ) == MB_SUCCESS );
      CHECK( GeomUtil::box_tri_overlap( coords, cen, hdim ) );
    }
  } while (MB_SUCCESS == (rval = iter.step()));  
  
  CHECK( MB_ENTITY_NOT_FOUND == rval );
  CHECK( expected_tris == all_tris );
}

/* utility method - check that all tris share a vertex */
void check_common_vertex( Interface* moab,
                          const EntityHandle* tris,
                          unsigned num_tri,
                          CartVect point )
{
  for (unsigned i = 0; i < num_tri; ++i)
    CHECK( MBTRI == moab->type_from_handle( tris[i] ) );
  
  ErrorCode rval;
  CartVect tri_coords[3];
  const EntityHandle* conn;
  int conn_len;
  rval = moab->get_connectivity( tris[0], conn, conn_len );
  CHECK( MB_SUCCESS == rval );
  rval = moab->get_coords( conn, 3, tri_coords[0].array() );
  CHECK( MB_SUCCESS == rval );
  tri_coords[0] -= point;
  tri_coords[1] -= point;
  tri_coords[2] -= point;
  EntityHandle vertex = 0;
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
}                        
 
/* utility method - check that all tris share a vertex */
void check_point_in_triangles( Interface* moab,
                               const EntityHandle* tris,
                               unsigned num_tris,
                               CartVect point )
{
  ErrorCode rval;
  CartVect tri_coords[3], tript;
  const EntityHandle* conn;
  int conn_len;

  for (unsigned i = 0; i < num_tris; ++i) {
    CHECK( MBTRI == moab->type_from_handle( tris[i] ) );

    rval = moab->get_connectivity( tris[i], conn, conn_len );
    CHECK( MB_SUCCESS == rval );
    rval = moab->get_coords( conn, 3, tri_coords[0].array() );
    CHECK( MB_SUCCESS == rval );
    
    GeomUtil::closest_location_on_tri( point, tri_coords, tript );
    tript -= point;
    CHECK( fabs(tript[0]) < 1e-6 );
    CHECK( fabs(tript[1]) < 1e-6 );
    CHECK( fabs(tript[2]) < 1e-6 );
  }
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
void create_simple_2d_tree( AdaptiveKDTree& tool, 
                            EntityHandle& root,
                            EntityHandle leaves[9] = 0 )
{
  ErrorCode rval;
  AdaptiveKDTree::Plane plane;

    // create a single-node tree
  double min[3] = { -5, -4, -1 };
  double max[3] = {  5,  4,  1 };
  rval = tool.create_tree( min, max, root );
  CHECK( MB_SUCCESS == rval );

    // get iterator for tree
  AdaptiveKDTreeIter iter;
  rval = tool.get_tree_iterator( root, iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal(iter, -5, -4, -1, 5, 4, 1 ) );
   
   // split plane (1)
  plane.norm = AdaptiveKDTree::Y;
  plane.coord = 0.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, 5, 0, 1 ) );
  
   // split plane (2)
  plane.norm = AdaptiveKDTree::X;
  plane.coord = -1.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, 0, 1 ) );
  
   // split plane (3), leaf [1]
  plane.norm = AdaptiveKDTree::Y;
  plane.coord = -2.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, -2, 1 ) );
  if (leaves)
    leaves[1] = iter.handle();
  
    // leaf [2]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -2, -1, -1, 0, 1 ) );
  if (leaves)
    leaves[2] = iter.handle();
  
    // non-leaf right of split plane (2)
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -1, -4, -1, 5, 0, 1 ) );
  
    // split plane (4) and leaf [3]
  plane.norm = AdaptiveKDTree::X;
  plane.coord = 4.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -1, -4, -1, 4, 0, 1 ) );
  if (leaves)
    leaves[3] = iter.handle();
  
    // leaf [4]
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 4, -4, -1, 5, 0, 1 ) );
  if (leaves)
    leaves[4] = iter.handle();
  
    // non-leaf above split plane (1)
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, 5, 4, 1 ) );
  
    // split plane (5)
  plane.norm = AdaptiveKDTree::X;
  plane.coord = 0.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, 0, 4, 1 ) );
  
    // split plane (6) and leaf [5]
  plane.norm = AdaptiveKDTree::X;
  plane.coord = -3.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, 0, -1, -3, 4, 1 ) );
  if (leaves)
    leaves[5] = iter.handle();
  
    // leaf [6];
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -3, 0, -1, 0, 4, 1 ) );
  if (leaves)
    leaves[6] = iter.handle();
  
    // non-leaf right of split plane (5)
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 0, -1, 5, 4, 1 ) );
  
    // split plane (7) and leaf [7]
  plane.norm = AdaptiveKDTree::Y;
  plane.coord = 3.0;
  rval = tool.split_leaf( iter, plane );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 0, -1, 5, 3, 1 ) );
  if (leaves)
    leaves[7] = iter.handle();
  
    // leaf [8];
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, 0, 3, -1, 5, 4, 1 ) );
  if (leaves)
    leaves[8] = iter.handle();
  
    // the end
  rval = iter.step();
  CHECK( MB_ENTITY_NOT_FOUND == rval );
}

void leaf_iterator_test()
{
  ErrorCode rval;
  Core moab;
  Interface* mb = &moab;
  AdaptiveKDTree tool(mb);


    // create a single-node tree
  EntityHandle root;
  double min[3] = { -3, -2, -1 };
  double max[3] = {  1,  2,  3 };
  rval = tool.create_tree( min, max, root );
  CHECK( MB_SUCCESS == rval );

    // get iterator for tree
  AdaptiveKDTreeIter iter;
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
  EntityHandle leaves[9];
  create_simple_2d_tree( tool, root, leaves );

  
  /**** Now traverse tree again, and check neighbors of each leaf *****/
  std::vector<AdaptiveKDTreeIter> list;
  
    // leaf [1]
  
  rval = tool.get_tree_iterator( root, iter );
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -4, -1, -1, -2, 1 ) );
  CHECK( iter.handle() == leaves[1] );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::X, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::X, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 1 );
    // should be leaf [3]
  CHECK( box_equal( list.front(), -1, -4, -1, 4, 0, 1 ) );
  CHECK( list.front().handle() == leaves[3] );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Y, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Y, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 1 );
    // should be leaf [2]
  CHECK( box_equal( list.front(), -5, -2, -1, -1, 0, 1 ) );
  CHECK( list.front().handle() == leaves[2] );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Z, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Z, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
    // leaf [2]
  
  rval = iter.step();
  CHECK( MB_SUCCESS == rval );
  CHECK( box_equal( iter, -5, -2, -1, -1, 0, 1 ) );
  CHECK( iter.handle() == leaves[2] );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::X, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::X, false, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.size() == 1 );
    // should be leaf [3]
  CHECK( box_equal( list.front(), -1, -4, -1, 4, 0, 1 ) );
  CHECK( list.front().handle() == leaves[3] );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Y, true, list );
  CHECK( MB_SUCCESS == rval );
    // should be leaf [1]
  CHECK( list.size() == 1 );
  CHECK( box_equal( list.front(), -5, -4, -1, -1, -2, 1 ) );
  CHECK( list.front().handle() == leaves[1] );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Y, false, list );
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
  rval = iter.get_neighbors( AdaptiveKDTree::Z, true, list );
  CHECK( MB_SUCCESS == rval );
  CHECK( list.empty() );
  
  list.clear();
  rval = iter.get_neighbors( AdaptiveKDTree::Z, false, list );
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
}

void test_tree_merge_nodes()
{
    // build simple tree for tests
  ErrorCode rval;
  Core moab;
  Interface* mb = &moab;
  AdaptiveKDTree tool(mb);
  EntityHandle root;
  EntityHandle leaves[9];
  create_simple_2d_tree( tool, root, leaves );

    // get iterator for tree
  AdaptiveKDTreeIter iter;
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
}

void test_build_tree_bisect_triangles( )
{
  Core moab;
  AdaptiveKDTree tool( &moab );
  Range box_tris;
  
  AdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;
  EntityHandle root;
  ErrorCode rval;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_octahedron( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );
  
  settings.maxTreeDepth = 2;
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );
}
  
void test_closest_triangle()
{
  Core moab;
  AdaptiveKDTree tool( &moab );
  Range box_tris;
  
  AdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;
  EntityHandle root;
  ErrorCode rval;
  EntityHandle tri;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );

    // test closest to each corner of the box.
  for (unsigned i = 0; i < 8; ++i) {
    int x = i & 1 ? 1 : -1;
    int y = i & 2 ? 1 : -1;
    int z = i & 4 ? 1 : -1;
    double coords[] = { 2.0*x, 2.0*y, 2.0*z };
    CartVect point;
    rval = tool.closest_triangle( root, coords, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    CHECK( fabs(x - point[0]) < 1e-6 );
    CHECK( fabs(y - point[1]) < 1e-6 );
    CHECK( fabs(z - point[2]) < 1e-6 );
      // check point is at triangle vertex
    check_common_vertex( tool.moab(), &tri, 1, point );
  }
  
    // test location on each face
  const CartVect facepts[] = { CartVect( -1.0, 0.5, 0.5 ),
                                 CartVect(  1.0, 0.5, 0.5 ),
                                 CartVect( 0.5, -1.0, 0.5 ),
                                 CartVect( 0.5,  1.0, 0.5 ),
                                 CartVect( 0.5, 0.5, -1.0 ),
                                 CartVect( 0.5, 0.5,  1.0 ) };
  for (unsigned i = 0; i < 6; ++i) {
    CartVect point;
    rval = tool.closest_triangle( root, facepts[i].array(), point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    const CartVect diff = facepts[i] - point;
    CHECK( fabs(diff[0]) < 1e-6 );
    CHECK( fabs(diff[1]) < 1e-6 );
    CHECK( fabs(diff[2]) < 1e-6 );
       // check that point is contained within triangle
    check_point_in_triangles( tool.moab(), &tri, 1, point );
  }
  
    // test a location really far from the tree
  {
    const double far[] = { 0.75, 0.75, 200 };
    CartVect point;
    rval = tool.closest_triangle( root, far, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
    CHECK( fabs(point[0] - 0.75) < 1e-6 );    
    CHECK( fabs(point[1] - 0.75) < 1e-6 );    
    CHECK( fabs(point[2] - 1.00) < 1e-6 );    
       // check that point is contained within triangle
    check_point_in_triangles( tool.moab(), &tri, 1, point );
  }
  
    // now do it all again with a lot more triangles
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );

    // test closest to each corner of the box.
  for (unsigned i = 0; i < 8; ++i) {
    int x = i & 1 ? 1 : -1;
    int y = i & 2 ? 1 : -1;
    int z = i & 4 ? 1 : -1;
    double coords[] = { 4.0*x, 4.0*y, 4.0*z };
    CartVect point;
    rval = tool.closest_triangle( root, coords, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    CHECK( fabs(3.0*x - point[0]) < 1e-6 );
    CHECK( fabs(3.0*y - point[1]) < 1e-6 );
    CHECK( fabs(3.0*z - point[2]) < 1e-6 );
      // check point is at triangle vertex
    check_common_vertex( tool.moab(), &tri, 1, point );
  }
  
    // test location on each face
  const CartVect facepts2[] = { CartVect( -3.0, 0.5, 0.5 ),
                                 CartVect(  3.0, 0.5, 0.5 ),
                                 CartVect( 0.5, -3.0, 0.5 ),
                                 CartVect( 0.5,  3.0, 0.5 ),
                                 CartVect( 0.5, 0.5, -3.0 ),
                                 CartVect( 0.5, 0.5,  3.0 ) };
  for (unsigned i = 0; i < 6; ++i) {
    CartVect point;
    rval = tool.closest_triangle( root, facepts2[i].array(), point.array(), tri );
    CHECK( MB_SUCCESS == rval );
      // check result position
    const CartVect diff = facepts2[i] - point;
    CHECK( fabs(diff[0]) < 1e-6 );
    CHECK( fabs(diff[1]) < 1e-6 );
    CHECK( fabs(diff[2]) < 1e-6 );
       // check that point is contained within triangle
    check_point_in_triangles( tool.moab(), &tri, 1, point );
  }
  
    // test a location really far from the tree
  {
    const double far[] = { 2.75, 2.75, 200 };
    CartVect point;
    rval = tool.closest_triangle( root, far, point.array(), tri );
    CHECK( MB_SUCCESS == rval );
    CHECK( fabs(point[0] - 2.75) < 1e-6 );    
    CHECK( fabs(point[1] - 2.75) < 1e-6 );    
    CHECK( fabs(point[2] - 3.00) < 1e-6 );    
       // check that point is contained within triangle
    check_point_in_triangles( tool.moab(), &tri, 1, point );
  }
}

void test_sphere_intersect_triangles()
{
  Core moab;
  AdaptiveKDTree tool( &moab );
  Range box_tris;
  
  AdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;
  EntityHandle root;
  ErrorCode rval;
  std::vector<EntityHandle> triangles;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );

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
    check_common_vertex( tool.moab(), &triangles[0], triangles.size(), CartVect(center) );
  }
  
    // now do it all again with a lot more triangles

  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_large( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );

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
    check_common_vertex( tool.moab(), &triangles[0], triangles.size(), CartVect(center) );
  }
}

void test_ray_intersect_triangles()
{
  Core moab;
  AdaptiveKDTree tool( &moab );
  Range box_tris;
  
  AdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 1;
  settings.candidateSplitsPerDir = 1;
  settings.candidatePlaneSet = AdaptiveKDTree::SUBDIVISION;
  EntityHandle root;
  ErrorCode rval;
  std::vector<EntityHandle> tris;
  std::vector<double> dists;
  
  moab.delete_mesh(); box_tris.clear();
  build_triangle_box_small( &moab, box_tris );
  rval = tool.build_tree( box_tris, root, &settings );
  CHECK( MB_SUCCESS == rval );
  test_valid_tree( &tool, root, settings, box_tris );
  
    // test ray through box parallel to X axis
  CartVect dir( 1, 0, 0 );
  CartVect pt( -2, 0.75, 0.75 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 3 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 || fabs(dists[i] - 3) < 1e-6 );
    CartVect tript = pt + dists[i] * dir;
    check_point_in_triangles( &moab, &tris[i], 1, tript );
  }
   
    // test ray through box parallel to X axis, closest tri only
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists, 1 );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 1 );
  CHECK( dists.size() == tris.size() );
  CHECK( fabs(dists[0] - 1) < 1e-6  );
  check_point_in_triangles( &moab, &tris[0], 1, pt + dists[0] * dir );

    // test ray ending within box, parallel to X axis
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists, 0, 2.0 );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 1 );
  CHECK( dists.size() == tris.size() );
  CHECK( fabs(dists[0] - 1) < 1e-6  );
  check_point_in_triangles( &moab, &tris[0], 1, pt + dists[0] * dir );
  
    // test ray starting within box parallel to X axis
  tris.clear(); dists.clear();
  pt = CartVect( 0, .75, .75 );
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 2 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 );
    CartVect tript = pt + dists[i] * dir;
    check_point_in_triangles( &moab, &tris[i], 1, tript );
  }
  
    // test skew ray through box
  dir = CartVect( 0.5*sqrt(2.0), 0.5*sqrt(2.0), 0.0 );
  pt = CartVect( 0, -1.5, 0 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 2 );
  CHECK( dists.size() == tris.size() );
  if (dists[0] < dists[1]) {
    CHECK( fabs(dists[0] - 0.5*sqrt(2.0)) < 1e-6 );
    CHECK( fabs(dists[1] -     sqrt(2.0)) < 1e-6 );
  }
  check_point_in_triangles( &moab, &tris[0], 1, pt + dists[0] * dir );
  check_point_in_triangles( &moab, &tris[1], 1, pt + dists[1] * dir );
  
  
    // test ray through box parallel to -Y axis
  dir = CartVect( 0, -1, 0 );
  pt = CartVect( 0, 2, 0 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 4 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 || fabs(dists[i] - 3) < 1e-6 );
    CartVect tript = pt + dists[i] * dir;
    check_point_in_triangles( &moab, &tris[i], 1, tript );
  }
  
    // test ray through box parallel to Z axis
  dir = CartVect( 0, 0, 1 );
  pt = CartVect( 0, 0, -2 );
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 4 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 || fabs(dists[i] - 3) < 1e-6 );
    CartVect tript = pt + dists[i] * dir;
    check_point_in_triangles( &moab, &tris[i], 1, tript );
  }
  
    // test ray through box parallel to Z axis, limit 2 first 2 intersections
  tris.clear(); dists.clear();
  rval = tool.ray_intersect_triangles( root, 1e-6, dir.array(), pt.array(), tris, dists, 2 );
  CHECK( MB_SUCCESS == rval );
  CHECK( tris.size() == 2 );
  CHECK( dists.size() == tris.size() );
  for ( unsigned i = 0; i < dists.size(); ++i ) {
    CHECK( fabs(dists[i] - 1) < 1e-6 );
    CartVect tript = pt + dists[i] * dir;
    check_point_in_triangles( &moab, &tris[i], 1, tript );
  }
}

void test_leaf_volume()
{
  Core moab;
  AdaptiveKDTree tool( &moab );
  EntityHandle root;
  
  create_simple_2d_tree( tool, root );
  AdaptiveKDTreeIter iter;
  ErrorCode rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
  CHECK_REAL_EQUAL( 16.0, iter.volume(), 1e-12 ); // 1
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL( 16.0, iter.volume(), 1e-12 ); // 2
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL( 40.0, iter.volume(), 1e-12 ); // 3
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL(  8.0, iter.volume(), 1e-12 ); // 4
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL( 16.0, iter.volume(), 1e-12 ); // 5
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL( 24.0, iter.volume(), 1e-12 ); // 6
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL( 30.0, iter.volume(), 1e-12 ); // 7
  CHECK_ERR(iter.step());
  CHECK_REAL_EQUAL( 10.0, iter.volume(), 1e-12 ); // 8
}

void test_leaf_sibling()
{
  ErrorCode rval;
  Core moab;
  AdaptiveKDTree tool( &moab );
  EntityHandle root;
  
  create_simple_2d_tree( tool, root );
  AdaptiveKDTreeIter iter1, iter2;
  rval = tool.get_tree_iterator( root, iter1 );
  CHECK_ERR(rval);
  rval = tool.get_tree_iterator( root, iter2 );
  CHECK_ERR(rval);
  
    // iter1 == 1, iter2 == 1
  CHECK( !iter1.is_sibling( iter1 ) );
  CHECK( !iter1.is_sibling( iter1.handle() ) );
  CHECK( !iter1.is_sibling( iter2 ) );
  CHECK(  iter1.sibling_is_forward() );
  
    // iter1 == 1, iter2 == 2
  rval = iter2.step();
  CHECK_ERR(rval);
  CHECK(  iter1.is_sibling( iter2 ) );
  CHECK(  iter1.is_sibling( iter2.handle() ) );
  CHECK(  iter2.is_sibling( iter1 ) );
  CHECK(  iter2.is_sibling( iter1.handle() ) );
  CHECK( !iter2.sibling_is_forward() );
  
    // iter1 == 1, iter2 == 3
  rval = iter2.step();
  CHECK_ERR(rval);
  CHECK( !iter1.is_sibling( iter2 ) );
  CHECK( !iter1.is_sibling( iter2.handle() ) );
  CHECK( !iter2.is_sibling( iter1 ) );
  CHECK( !iter2.is_sibling( iter1.handle() ) );
  CHECK(  iter2.sibling_is_forward() );
  
    // iter1 == 2, iter2 == 3
  rval = iter1.step();
  CHECK_ERR(rval);
  CHECK( !iter1.is_sibling( iter2 ) );
  CHECK( !iter1.is_sibling( iter2.handle() ) );
  CHECK( !iter2.is_sibling( iter1 ) );
  CHECK( !iter2.is_sibling( iter1.handle() ) );
  CHECK( !iter1.sibling_is_forward() );

    // iter1 == 4, iter2 == 3
  rval = iter1.step();
  CHECK_ERR(rval);
  CHECK_EQUAL( iter1.handle(), iter2.handle() );
  rval = iter1.step();
  CHECK_ERR(rval);
  CHECK(  iter1.is_sibling( iter2 ) );
  CHECK(  iter1.is_sibling( iter2.handle() ) );
  CHECK(  iter2.is_sibling( iter1 ) );
  CHECK(  iter2.is_sibling( iter1.handle() ) );
  CHECK( !iter1.sibling_is_forward() );
}  
  

void test_leaf_intersects_plane()
{
  ErrorCode rval;
  Core moab;
  AdaptiveKDTree tool( &moab );
  
  EntityHandle root;
  const double min[3] = { -5, -4, -1 };
  const double max[3] = {  1,  2,  3 };
  rval = tool.create_tree( min, max, root );
  CHECK_ERR(rval);
  
  AdaptiveKDTreeIter iter;
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
  AdaptiveKDTree::Plane x1 = { min[0] - 1, AdaptiveKDTree::X };
  CHECK( !iter.intersects( x1 ) );
  AdaptiveKDTree::Plane x2 = { min[0]    , AdaptiveKDTree::X };
  CHECK(  iter.intersects( x2 ) );
  AdaptiveKDTree::Plane x3 = { min[0] + 1, AdaptiveKDTree::X };
  CHECK(  iter.intersects( x3 ) );
  AdaptiveKDTree::Plane x4 = { max[0] - 1, AdaptiveKDTree::X };
  CHECK(  iter.intersects( x4 ) );
  AdaptiveKDTree::Plane x5 = { max[0]    , AdaptiveKDTree::X };
  CHECK(  iter.intersects( x5 ) );
  AdaptiveKDTree::Plane x6 = { max[0] + 1, AdaptiveKDTree::X };
  CHECK( !iter.intersects( x6 ) );
  
  AdaptiveKDTree::Plane y1 = { min[1] - 1, AdaptiveKDTree::Y };
  CHECK( !iter.intersects( y1 ) );
  AdaptiveKDTree::Plane y2 = { min[1]    , AdaptiveKDTree::Y };
  CHECK(  iter.intersects( y2 ) );
  AdaptiveKDTree::Plane y3 = { min[1] + 1, AdaptiveKDTree::Y };
  CHECK(  iter.intersects( y3 ) );
  AdaptiveKDTree::Plane y4 = { max[1] - 1, AdaptiveKDTree::Y };
  CHECK(  iter.intersects( y4 ) );
  AdaptiveKDTree::Plane y5 = { max[1]    , AdaptiveKDTree::Y };
  CHECK(  iter.intersects( y5 ) );
  AdaptiveKDTree::Plane y6 = { max[1] + 1, AdaptiveKDTree::Y };
  CHECK( !iter.intersects( y6 ) );
  
  AdaptiveKDTree::Plane z1 = { min[2] - 1, AdaptiveKDTree::Z };
  CHECK( !iter.intersects( z1 ) );
  AdaptiveKDTree::Plane z2 = { min[2]    , AdaptiveKDTree::Z };
  CHECK(  iter.intersects( z2 ) );
  AdaptiveKDTree::Plane z3 = { min[2] + 1, AdaptiveKDTree::Z };
  CHECK(  iter.intersects( z3 ) );
  AdaptiveKDTree::Plane z4 = { max[2] - 1, AdaptiveKDTree::Z };
  CHECK(  iter.intersects( z4 ) );
  AdaptiveKDTree::Plane z5 = { max[2]    , AdaptiveKDTree::Z };
  CHECK(  iter.intersects( z5 ) );
  AdaptiveKDTree::Plane z6 = { max[2] + 1, AdaptiveKDTree::Z };
  CHECK( !iter.intersects( z6 ) );
}

#define CHECK_RAY_XSECTS( PT, DIR, T_IN, T_OUT ) do { \
  CHECK(iter.intersect_ray( (PT), (DIR), t_in, t_out )); \
  CHECK_REAL_EQUAL( (T_IN), t_in, 1e-6 ); \
  CHECK_REAL_EQUAL( (T_OUT), t_out, 1e-6 ); \
  } while(false)
    
void test_leaf_intersects_ray()
{
  ErrorCode rval;
  Core moab;
  AdaptiveKDTree tool( &moab );
  double t_in, t_out;
  
  EntityHandle root;
  const double min[3] = { -5, -4, -1 };
  const double max[3] = {  1,  2,  3 };
  rval = tool.create_tree( min, max, root );
  CHECK_ERR(rval);
  
  AdaptiveKDTreeIter iter;
  rval = tool.get_tree_iterator( root, iter );
  CHECK_ERR(rval);
  
    // start with point inside box
  const double pt1[] = { 0, 0, 0 };
  const double dir1[] = { 0.1, 0.1, 0.1 };
  CHECK_RAY_XSECTS( pt1, dir1, 0, 10 );
  const double dir2[] = { 5, 5, 5 };
  CHECK_RAY_XSECTS( pt1, dir2, 0, 0.2 );
  const double pxdir[] = { 1, 0, 0 };
  CHECK_RAY_XSECTS( pt1, pxdir, 0, 1 );
  const double nxdir[] = { -1, 0, 0 };
  CHECK_RAY_XSECTS( pt1, nxdir, 0, 5 );
  const double pydir[] = { 0, 1, 0 };
  CHECK_RAY_XSECTS( pt1, pydir, 0, 2 );
  const double nydir[] = { 0, -1, 0 };
  CHECK_RAY_XSECTS( pt1, nydir, 0, 4 );
  const double pzdir[] = { 0, 0, 1 };
  CHECK_RAY_XSECTS( pt1, pzdir, 0, 3 );
  const double nzdir[] = { 0, 0, -1 };
  CHECK_RAY_XSECTS( pt1, nzdir, 0, 1 );
  
    // point below box
  const double pt2[] = { 0, 0, -2 };
  CHECK_RAY_XSECTS( pt2, dir1, 10, 10 );
  CHECK_RAY_XSECTS( pt2, dir2, 0.2, 0.2 );
  CHECK(!iter.intersect_ray( pt2, pxdir, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt2, nxdir, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt2, pydir, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt2, nydir, t_in, t_out ));
  CHECK_RAY_XSECTS( pt2, pzdir, 1, 5 );
  CHECK(!iter.intersect_ray( pt2, nzdir, t_in, t_out ));
  
    // point right of box
  const double pt3[] = { 3, 0, 0 };
  CHECK(!iter.intersect_ray( pt3, dir1, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt3, dir2, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt3, pxdir, t_in, t_out ));
  CHECK_RAY_XSECTS( pt3, nxdir, 2, 8 );
  CHECK(!iter.intersect_ray( pt3, pydir, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt3, nydir, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt3, pzdir, t_in, t_out ));
  CHECK(!iter.intersect_ray( pt3, nzdir, t_in, t_out ));
  
    // a few more complex test cases
  const double dira[] = { -3, 0, 3 };
  CHECK_RAY_XSECTS( pt3, dira, 2./3., 1.0 );
  const double dirb[] = { -2, 0, 3.1 };
  CHECK(!iter.intersect_ray( pt3, dirb, t_in, t_out ));
}

int main()
{
  int error_count = 0;
  
  error_count += RUN_TEST(leaf_iterator_test);
  error_count += RUN_TEST(test_build_tree_bisect_triangles);
  error_count += RUN_TEST(test_closest_triangle);
  error_count += RUN_TEST(test_sphere_intersect_triangles);
  error_count += RUN_TEST(test_ray_intersect_triangles);
  error_count += RUN_TEST(test_tree_merge_nodes);
  error_count += RUN_TEST(test_leaf_volume);
  error_count += RUN_TEST(test_leaf_sibling);
  error_count += RUN_TEST(test_leaf_intersects_plane);
  error_count += RUN_TEST(test_leaf_intersects_ray);
  return error_count;
}
