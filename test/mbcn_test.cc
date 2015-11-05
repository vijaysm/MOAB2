#include "TestUtil.hpp"
#include "moab/CN.hpp"

using namespace moab;

void test_dimension_pair();
void test_type_names();
void test_dimension();
void test_vertices_per_entity();
void test_num_sub_entities();

void test_sub_entity_type_vtx();
void test_sub_entity_type_edge();
void test_sub_entity_type_tri();
void test_sub_entity_type_quad();
void test_sub_entity_type_tet();
void test_sub_entity_type_pyr();
void test_sub_entity_type_pri();
void test_sub_entity_type_knife();
void test_sub_entity_type_hex();

void test_sub_entity_indices_vtx();
void test_sub_entity_indices_edge();
void test_sub_entity_indices_tri();
void test_sub_entity_indices_quad();
void test_sub_entity_indices_tet();
void test_sub_entity_indices_pyr();
void test_sub_entity_indices_pri();
void test_sub_entity_indices_hex();

void test_side_number_tri();
void test_side_number_quad();
void test_side_number_tet();
void test_side_number_pyr();
void test_side_number_pri();
void test_side_number_hex();

void test_opposite_side_tri();
void test_opposite_side_quad();
void test_opposite_side_tet();
void test_opposite_side_hex();

void test_has_mid_nodes( EntityType type );
void test_has_mid_nodes_edge() { test_has_mid_nodes(MBEDGE); }
void test_has_mid_nodes_tri()  { test_has_mid_nodes(MBTRI); }
void test_has_mid_nodes_quad() { test_has_mid_nodes(MBQUAD); }
void test_has_mid_nodes_tet()  { test_has_mid_nodes(MBTET); }
void test_has_mid_nodes_pyr()  { test_has_mid_nodes(MBPYRAMID); }
void test_has_mid_nodes_pri()  { test_has_mid_nodes(MBPRISM); }
void test_has_mid_nodes_knife(){ test_has_mid_nodes(MBKNIFE); }
void test_has_mid_nodes_hex()  { test_has_mid_nodes(MBHEX); }

void test_ho_node_parent();
void test_ho_node_index();

void test_sub_entity_nodes( EntityType parent, int sub_dimension );
void test_sub_entity_nodes( EntityType parent, int num_nodes, int sub_dimension );
void test_sub_entity_nodes_tri_edges()  { test_sub_entity_nodes(MBTRI,     1 ); }
void test_sub_entity_nodes_quad_edges() { test_sub_entity_nodes(MBQUAD,    1 ); }
void test_sub_entity_nodes_tet_edges()  { test_sub_entity_nodes(MBTET,     1 ); }
void test_sub_entity_nodes_tet_faces()  { test_sub_entity_nodes(MBTET,     2 ); }
void test_sub_entity_nodes_pyr_edges()  { test_sub_entity_nodes(MBPYRAMID, 1 ); }
void test_sub_entity_nodes_pyr_faces()  { test_sub_entity_nodes(MBPYRAMID, 2 ); }
void test_sub_entity_nodes_pri_edges()  { test_sub_entity_nodes(MBPRISM,   1 ); }
void test_sub_entity_nodes_pri_faces()  { test_sub_entity_nodes(MBPRISM,   2 ); }
void test_sub_entity_nodes_kni_edges()  { test_sub_entity_nodes(MBKNIFE,   1 ); }
void test_sub_entity_nodes_kni_faces()  { test_sub_entity_nodes(MBKNIFE,   2 ); }
void test_sub_entity_nodes_hex_edges()  { test_sub_entity_nodes(MBHEX,     1 ); }
void test_sub_entity_nodes_hex_faces()  { test_sub_entity_nodes(MBHEX,     2 ); }

int main()
{
  int result = 0;
  result += RUN_TEST(test_dimension_pair);
  result += RUN_TEST(test_type_names);
  result += RUN_TEST(test_dimension);
  result += RUN_TEST(test_vertices_per_entity);
  result += RUN_TEST(test_num_sub_entities);
  
  result += RUN_TEST(test_sub_entity_type_vtx);
  result += RUN_TEST(test_sub_entity_type_edge);
  result += RUN_TEST(test_sub_entity_type_tri);
  result += RUN_TEST(test_sub_entity_type_quad);
  result += RUN_TEST(test_sub_entity_type_tet);
  result += RUN_TEST(test_sub_entity_type_pyr);
  result += RUN_TEST(test_sub_entity_type_pri);
  result += RUN_TEST(test_sub_entity_type_knife);
  result += RUN_TEST(test_sub_entity_type_hex);
  
  result += RUN_TEST(test_sub_entity_indices_vtx);
  result += RUN_TEST(test_sub_entity_indices_edge);
  result += RUN_TEST(test_sub_entity_indices_tri);
  result += RUN_TEST(test_sub_entity_indices_quad);
  result += RUN_TEST(test_sub_entity_indices_tet);
  result += RUN_TEST(test_sub_entity_indices_pyr);
  result += RUN_TEST(test_sub_entity_indices_pri);
  result += RUN_TEST(test_sub_entity_indices_hex);
  
  result += RUN_TEST(test_side_number_tri);
  result += RUN_TEST(test_side_number_quad);
  result += RUN_TEST(test_side_number_tet);
  result += RUN_TEST(test_side_number_pyr);
  result += RUN_TEST(test_side_number_pri);
  result += RUN_TEST(test_side_number_hex);
  
  result += RUN_TEST(test_opposite_side_tri);
  result += RUN_TEST(test_opposite_side_quad);
  result += RUN_TEST(test_opposite_side_tet);
  result += RUN_TEST(test_opposite_side_hex);
  
  result += RUN_TEST(test_has_mid_nodes_edge);
  result += RUN_TEST(test_has_mid_nodes_tri);
  result += RUN_TEST(test_has_mid_nodes_quad);
  result += RUN_TEST(test_has_mid_nodes_tet);
  result += RUN_TEST(test_has_mid_nodes_pyr);
  result += RUN_TEST(test_has_mid_nodes_pri);
  result += RUN_TEST(test_has_mid_nodes_knife);
  result += RUN_TEST(test_has_mid_nodes_hex);  

  result += RUN_TEST(test_sub_entity_nodes_tri_edges);  
  result += RUN_TEST(test_sub_entity_nodes_quad_edges);  
  result += RUN_TEST(test_sub_entity_nodes_tet_edges);  
  result += RUN_TEST(test_sub_entity_nodes_tet_faces);  
  result += RUN_TEST(test_sub_entity_nodes_pyr_edges);  
  result += RUN_TEST(test_sub_entity_nodes_pyr_faces);  
  result += RUN_TEST(test_sub_entity_nodes_pri_edges);  
  result += RUN_TEST(test_sub_entity_nodes_pri_faces);  
  result += RUN_TEST(test_sub_entity_nodes_kni_edges);  
  result += RUN_TEST(test_sub_entity_nodes_kni_faces);  
  result += RUN_TEST(test_sub_entity_nodes_hex_edges);  
  result += RUN_TEST(test_sub_entity_nodes_hex_faces);  

  result += RUN_TEST(test_ho_node_parent);
  result += RUN_TEST(test_ho_node_index);
  return result;
}

const EntityType elem_types[] = { MBTRI,
                                    MBQUAD,
                                    MBTET,
                                    MBPYRAMID,
                                    MBPRISM,
                                    MBHEX,
                                    MBMAXTYPE };
// const int num_elem_types = sizeof(elem_types)/sizeof(elem_types[0]) - 1;

void test_dimension_pair()
{
  DimensionPair dp;
  
  dp = CN::TypeDimensionMap[0];
  CHECK_EQUAL( MBVERTEX, dp.first );  
  CHECK_EQUAL( MBVERTEX, dp.second );
  
  dp = CN::TypeDimensionMap[1];
  CHECK_EQUAL( MBEDGE, dp.first );  
  CHECK_EQUAL( MBEDGE, dp.second );
  
  dp = CN::TypeDimensionMap[2];
  CHECK_EQUAL( MBTRI, dp.first );  
  CHECK_EQUAL( MBPOLYGON, dp.second );
  
  dp = CN::TypeDimensionMap[3];
  CHECK_EQUAL( MBTET, dp.first );  
  CHECK_EQUAL( MBPOLYHEDRON, dp.second );
}

void test_type_names()
{
  for (EntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
    const char* name = CN::EntityTypeName(t);
    CHECK_EQUAL( t, CN::EntityTypeFromName(name) );
  }
}

void test_dimension()
{
  CHECK_EQUAL( 0, CN::Dimension(MBVERTEX) );
  CHECK_EQUAL( 1, CN::Dimension(MBEDGE) );
  CHECK_EQUAL( 2, CN::Dimension(MBTRI) );
  CHECK_EQUAL( 2, CN::Dimension(MBQUAD) );
  CHECK_EQUAL( 2, CN::Dimension(MBPOLYGON) );
  CHECK_EQUAL( 3, CN::Dimension(MBTET) );
  CHECK_EQUAL( 3, CN::Dimension(MBPYRAMID) );
  CHECK_EQUAL( 3, CN::Dimension(MBPRISM) );
  CHECK_EQUAL( 3, CN::Dimension(MBKNIFE) );
  CHECK_EQUAL( 3, CN::Dimension(MBHEX) );
  CHECK_EQUAL( 3, CN::Dimension(MBPOLYHEDRON) );
}

void test_vertices_per_entity()
{
  CHECK_EQUAL( 1, CN::VerticesPerEntity(MBVERTEX) );
  CHECK_EQUAL( 2, CN::VerticesPerEntity(MBEDGE) );
  CHECK_EQUAL( 3, CN::VerticesPerEntity(MBTRI) );
  CHECK_EQUAL( 4, CN::VerticesPerEntity(MBQUAD) );
  CHECK_EQUAL( 4, CN::VerticesPerEntity(MBTET) );
  CHECK_EQUAL( 5, CN::VerticesPerEntity(MBPYRAMID) );
  CHECK_EQUAL( 6, CN::VerticesPerEntity(MBPRISM) );
  CHECK_EQUAL( 7, CN::VerticesPerEntity(MBKNIFE) );
  CHECK_EQUAL( 8, CN::VerticesPerEntity(MBHEX) );
}

void test_num_sub_entities()
{
  CHECK_EQUAL( 1, CN::NumSubEntities(MBVERTEX, 0));

  CHECK_EQUAL( 2, CN::NumSubEntities(MBEDGE, 0));
  CHECK_EQUAL( 1, CN::NumSubEntities(MBEDGE, 1));

  CHECK_EQUAL( 3, CN::NumSubEntities(MBTRI, 0));
  CHECK_EQUAL( 3, CN::NumSubEntities(MBTRI, 1));
  CHECK_EQUAL( 1, CN::NumSubEntities(MBTRI, 2));

  CHECK_EQUAL( 4, CN::NumSubEntities(MBQUAD, 0));
  CHECK_EQUAL( 4, CN::NumSubEntities(MBQUAD, 1));
  CHECK_EQUAL( 1, CN::NumSubEntities(MBQUAD, 2));

  CHECK_EQUAL( 4, CN::NumSubEntities(MBTET, 0));
  CHECK_EQUAL( 6, CN::NumSubEntities(MBTET, 1));
  CHECK_EQUAL( 4, CN::NumSubEntities(MBTET, 2));

  CHECK_EQUAL( 5, CN::NumSubEntities(MBPYRAMID, 0));
  CHECK_EQUAL( 8, CN::NumSubEntities(MBPYRAMID, 1));
  CHECK_EQUAL( 5, CN::NumSubEntities(MBPYRAMID, 2));

  CHECK_EQUAL( 6, CN::NumSubEntities(MBPRISM, 0));
  CHECK_EQUAL( 9, CN::NumSubEntities(MBPRISM, 1));
  CHECK_EQUAL( 5, CN::NumSubEntities(MBPRISM, 2));

  CHECK_EQUAL( 7, CN::NumSubEntities(MBKNIFE, 0));
  CHECK_EQUAL(10, CN::NumSubEntities(MBKNIFE, 1));
  CHECK_EQUAL( 5, CN::NumSubEntities(MBKNIFE, 2));

  CHECK_EQUAL( 8, CN::NumSubEntities(MBHEX, 0));
  CHECK_EQUAL( 12, CN::NumSubEntities(MBHEX, 1));
  CHECK_EQUAL( 6, CN::NumSubEntities(MBHEX, 2));
}

void do_test_sub_entity_type_2d( EntityType type )
{
  for (int j = 0; j < CN::VerticesPerEntity(type); ++j) {
    CHECK_EQUAL( MBVERTEX, CN::SubEntityType(type, 0, j ) );
    CHECK_EQUAL( MBEDGE,   CN::SubEntityType(type, 1, j ) );
  }
  CHECK_EQUAL( type, CN::SubEntityType(type, 2, 0) );
}

void do_test_sub_entity_type_3d( EntityType type,
                                 int num_faces,
                                 const EntityType* face_types )
{
  for (int j = 0; j < CN::VerticesPerEntity(type); ++j) {
    CHECK_EQUAL( MBVERTEX, CN::SubEntityType(type, 0, j ) );
  }

  for (int j = 0; j < CN::NumSubEntities(type,1); ++j) {
    CHECK_EQUAL( MBEDGE, CN::SubEntityType(type, 1, j ) );
  }

  for (int j = 0; j < num_faces; ++j) {
    EntityType sub_type = CN::SubEntityType( type, 2, j );
    CHECK_EQUAL( face_types[j], sub_type );
  }

  CHECK_EQUAL( type, CN::SubEntityType(type, 3, 0) );
}

void test_sub_entity_type_vtx()
{
  CHECK_EQUAL( MBVERTEX, CN::SubEntityType(MBVERTEX, 0, 0 ));
}

void test_sub_entity_type_edge()
{
  CHECK_EQUAL( MBVERTEX, CN::SubEntityType(MBEDGE, 0, 0 ));
  CHECK_EQUAL( MBVERTEX, CN::SubEntityType(MBEDGE, 0, 1 ));
  CHECK_EQUAL( MBEDGE  , CN::SubEntityType(MBEDGE, 1, 0 ));
}

void test_sub_entity_type_tri()
{
  do_test_sub_entity_type_2d(MBTRI);
}

void test_sub_entity_type_quad()
{
  do_test_sub_entity_type_2d(MBQUAD);
}

void test_sub_entity_type_tet()
{
  const EntityType types[] = { MBTRI, MBTRI, MBTRI, MBTRI };
  do_test_sub_entity_type_3d( MBTET, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_pyr()
{
  const EntityType types[] = { MBTRI, MBTRI, MBTRI, MBTRI, MBQUAD};
  do_test_sub_entity_type_3d( MBPYRAMID, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_pri()
{
  const EntityType types[] = { MBQUAD, MBQUAD, MBQUAD, MBTRI, MBTRI };
  do_test_sub_entity_type_3d( MBPRISM, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_knife()
{
  const EntityType types[] = { MBQUAD, MBQUAD, MBQUAD, MBQUAD, MBQUAD };
  do_test_sub_entity_type_3d( MBKNIFE, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_hex()
{
  const EntityType types[] = { MBQUAD, MBQUAD, MBQUAD, MBQUAD, MBQUAD, MBQUAD };
  do_test_sub_entity_type_3d( MBHEX, sizeof(types)/sizeof(types[0]), types );
}


void test_0d_sub_entity_indices( EntityType type, int num_vtx )
{
  for (int i = 0; i < num_vtx; ++i) {
    // zero input array
    int indices[2] = { 0, -100 };
    // check correct results
    CN::SubEntityVertexIndices( type, 0, i, indices );
    CHECK_EQUAL( i, indices[0] );
    // check didn't write past end of array
    CHECK_EQUAL( -100, indices[1] );
  }
}

void test_1d_sub_entity_indices( EntityType type, int num_edges, 
                                 const int (*edge_indices)[2] )
{
  for (int i = 0; i < num_edges; ++i) {
    // zero input array
    int indices[3] = { 0, 0, -99 };
    // check correct results
    CN::SubEntityVertexIndices( type, 1, i, indices );
    if (edge_indices[i][0] == indices[0]) {
      CHECK_EQUAL( edge_indices[i][1], indices[1] );
    }
    else {
      CHECK_EQUAL( edge_indices[i][0], indices[1] );
      CHECK_EQUAL( edge_indices[i][1], indices[0] );
    }
    // check didn't write past end of array
    CHECK_EQUAL( -99, indices[2] );
  }
}

void test_2d_sub_entity_indices( EntityType type, int num_faces,
                                 const int (*face_indices)[5] )
{
  for (int i = 0; i < num_faces; ++i) {
    // zero input array
    int indices[5] = { 0, 0, 0, -99, -99 };
    // check correct results
    CN::SubEntityVertexIndices( type, 2, i, indices );
    const int num_vtx = face_indices[i][0];
    if (num_vtx != 3) CHECK_EQUAL( 4, num_vtx );
    const int* exp_index = face_indices[i] + 1;
    int off = std::find( indices, indices+num_vtx, exp_index[0] ) - indices;
    CHECK( off < num_vtx );
      /* Expect faces to be ordered such that CCW normal is outwards 
       *
    bool reverse = indices[(off+num_vtx-1)%num_vtx] == exp_index[1];
    if (reverse) {
      CHECK_EQUAL( exp_index[1], indices[(off+num_vtx-1)%num_vtx] );
      CHECK_EQUAL( exp_index[2], indices[(off+num_vtx-2)%num_vtx] );
      if (num_vtx == 4)
        CHECK_EQUAL( exp_index[3], indices[(off+num_vtx-3)%num_vtx] );
    }
    else {
        */
      CHECK_EQUAL( exp_index[1], indices[(off+1)%num_vtx] );
      CHECK_EQUAL( exp_index[2], indices[(off+2)%num_vtx] );
      if (num_vtx == 4)
        CHECK_EQUAL( exp_index[3], indices[(off+3)%num_vtx] );
      /*
    }
      */
    
    // check didn't write past end of array
    if (num_vtx == 3)
      CHECK_EQUAL( -99, indices[3] );
    CHECK_EQUAL( -99, indices[4] );
  }
}

void test_elem_as_sub_entity( EntityType type, int dim, int num_vertices )
{ 
  int indices[9] = { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
  CN::SubEntityVertexIndices( type, dim, 0, indices );
  for (int i = 0; i < num_vertices; ++i)
    CHECK_EQUAL( i, indices[i] );
  // make sure didn't write past end
  CHECK_EQUAL( -2, indices[num_vertices] );
}

void test_sub_entity_indices_vtx()
{
  test_elem_as_sub_entity( MBVERTEX, 0, 1 );
}

void test_sub_entity_indices_edge()
{
  test_0d_sub_entity_indices( MBEDGE, 2 );
  test_elem_as_sub_entity( MBEDGE, 1, 2 );
}    

void test_sub_entity_indices_tri()
{
  const int edges[3][2] = { { 0, 1 }, 
                            { 1, 2 },
                            { 2, 0 } };
  test_0d_sub_entity_indices( MBTRI, 3 );
  test_1d_sub_entity_indices( MBTRI, 3, edges );
  test_elem_as_sub_entity( MBTRI, 2, 3 );
}

void test_sub_entity_indices_quad()
{
  const int edges[4][2] = { { 0, 1 }, 
                            { 1, 2 },
                            { 2, 3 },
                            { 3, 0 } };
  test_0d_sub_entity_indices( MBQUAD, 4 );
  test_1d_sub_entity_indices( MBQUAD, 4, edges );
  test_elem_as_sub_entity( MBQUAD, 2, 4 );
}

void test_sub_entity_indices_tet()
{
  const EntityType type = MBTET;
  const int num_vtx = 4;
  const int edges[][2] = { { 0, 1 },
                           { 1, 2 },
                           { 2, 0 },
                           { 0, 3 },
                           { 1, 3 },
                           { 2, 3 } };
  const int faces[][5] = { { 3, 0, 1, 3, 0 },
                           { 3, 1, 2, 3, 0 },
                           { 3, 2, 0, 3, 0 },
                           { 3, 2, 1, 0, 0 } };
  test_0d_sub_entity_indices( type, num_vtx );
  test_1d_sub_entity_indices( type, sizeof(edges)/sizeof(edges[0]), edges );
  test_2d_sub_entity_indices( type, sizeof(faces)/sizeof(faces[0]), faces );
  test_elem_as_sub_entity( type, 3, num_vtx );
}

void test_sub_entity_indices_pyr()
{
  const EntityType type = MBPYRAMID;
  const int num_vtx = 5;
  const int edges[][2] = { { 0, 1 },
                           { 1, 2 },
                           { 2, 3 },
                           { 3, 0 },
                           { 0, 4 },
                           { 1, 4 },
                           { 2, 4 },
                           { 3, 4 } };
  const int faces[][5] = { { 3, 0, 1, 4, 0 },
                           { 3, 1, 2, 4, 0 },
                           { 3, 2, 3, 4, 0 },
                           { 3, 3, 0, 4, 0 },
                           { 4, 3, 2, 1, 0 }};
  test_0d_sub_entity_indices( type, num_vtx );
  test_1d_sub_entity_indices( type, sizeof(edges)/sizeof(edges[0]), edges );
  test_2d_sub_entity_indices( type, sizeof(faces)/sizeof(faces[0]), faces );
  test_elem_as_sub_entity( type, 3, num_vtx );
}

void test_sub_entity_indices_pri()
{
  const EntityType type = MBPRISM;
  const int num_vtx = 6;
  const int edges[][2] = { { 0, 1 },
                           { 1, 2 },
                           { 2, 0 },
                           { 0, 3 },
                           { 1, 4 },
                           { 2, 5 },
                           { 3, 4 },
                           { 4, 5 },
                           { 5, 3 } };
  const int faces[][5] = { { 4, 0, 1, 4, 3 },
                           { 4, 1, 2, 5, 4 },
                           { 4, 2, 0, 3, 5 },
                           { 3, 2, 1, 0, 0 },
                           { 3, 3, 4, 5, 0 } };
  test_0d_sub_entity_indices( type, num_vtx );
  test_1d_sub_entity_indices( type, sizeof(edges)/sizeof(edges[0]), edges );
  test_2d_sub_entity_indices( type, sizeof(faces)/sizeof(faces[0]), faces );
  test_elem_as_sub_entity( type, 3, num_vtx );
}

void test_sub_entity_indices_hex()
{
  const EntityType type = MBHEX;
  const int num_vtx = 8;
  const int edges[][2] = { { 0, 1 },
                           { 1, 2 },
                           { 2, 3 },
                           { 3, 0 },
                           { 0, 4 },
                           { 1, 5 },
                           { 2, 6 },
                           { 3, 7 },
                           { 4, 5 },
                           { 5, 6 },
                           { 6, 7 },
                           { 7, 4 } };
  const int faces[][5] = { { 4, 0, 1, 5, 4 },
                           { 4, 1, 2, 6, 5 },
                           { 4, 2, 3, 7, 6 },
                           { 4, 3, 0, 4, 7 },
                           { 4, 3, 2, 1, 0 },
                           { 4, 4, 5, 6, 7 } };
  test_0d_sub_entity_indices( type, num_vtx );
  test_1d_sub_entity_indices( type, sizeof(edges)/sizeof(edges[0]), edges );
  test_2d_sub_entity_indices( type, sizeof(faces)/sizeof(faces[0]), faces );
  test_elem_as_sub_entity( type, 3, num_vtx );
}

static void do_test_side_number_1d( EntityType type, int idx )
{
  // define a random handle list
  const int elem_verts[] = { 7400, 6233, 3027, 0454, 6839, 5391, 7735, 3603 };
  // get side indices
  int side_idx[4] = { 0, 0 };
  CN::SubEntityVertexIndices( type, 1, idx, side_idx );

  // "reversed" and "offset" are the same thing for edges.
  int side_conn[2] = { elem_verts[side_idx[0]], elem_verts[side_idx[1]] };
  int rev_conn[2] = { elem_verts[side_idx[1]], elem_verts[side_idx[0]] };
  int result_side = -100, result_sense = -100, result_offset = -100;
  int err= CN::SideNumber( type, elem_verts, side_conn, 2, 1, 
                             result_side, result_sense, result_offset );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( idx, result_side );
  CHECK_EQUAL( 1, result_sense );
  CHECK_EQUAL( 0, result_offset);
  err= CN::SideNumber( type, elem_verts, rev_conn, 2, 1, 
                         result_side, result_sense, result_offset );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( idx, result_side );
  CHECK(result_offset == 1 || result_sense == -1);
}     

static void do_test_side_number_2d( EntityType type, int idx )
{
  // define a random handle list
  const int elem_verts[] = { 7400, 6233, 3027, 0454, 6839, 5391, 7735, 3603 };
  // get side indices
  const int side_size = CN::VerticesPerEntity( CN::SubEntityType( type, 2, idx ) );
  int side_idx[4] = { 0, 0, 0, 0 };
  CN::SubEntityVertexIndices( type, 2, idx, side_idx );

  // for each possible forward or reverse offset
  for (int rev = -1; rev < 2; rev += 2) {
    for (int off = 0; off < side_size; ++off) {
      int side_conn[4]; side_conn[3] = 0;
      for (int i = 0; i < side_size; ++i)
        side_conn[(side_size+rev*i)%side_size] = elem_verts[side_idx[(i+off)%side_size]];
      
      int result_side = -100, result_sense = -100, result_offset = -100;
      int err = CN::SideNumber( type, elem_verts, side_conn, side_size, 2, 
                                  result_side, result_sense, result_offset );
      CHECK_EQUAL( 0, err );
      CHECK_EQUAL( idx, result_side );
      CHECK_EQUAL( rev, result_sense );
      CHECK_EQUAL( off, result_offset );
    }
  }
}

void test_side_number_tri()
{
  for (int side = 0; side < 3; ++side)
    do_test_side_number_1d( MBTRI, side );
}

void test_side_number_quad()
{
  for (int side = 0; side < 4; ++side)
    do_test_side_number_1d( MBQUAD, side );
}

void test_side_number_tet()
{
  for (int edge = 0; edge < 6; ++edge)
    do_test_side_number_1d( MBTET, edge );
  for (int face = 0; face < 4; ++face)
    do_test_side_number_2d( MBTET, face );
}

void test_side_number_pyr()
{
  for (int edge = 0; edge < 8; ++edge)
    do_test_side_number_1d( MBPYRAMID, edge );
  for (int face = 0; face < 5; ++face)
    do_test_side_number_2d( MBPYRAMID, face );
}

void test_side_number_pri()
{
  for (int edge = 0; edge < 9; ++edge)
    do_test_side_number_1d( MBPRISM, edge );
  for (int face = 0; face < 5; ++face)
    do_test_side_number_2d( MBPRISM, face );
}

void test_side_number_hex()
{
  for (int edge = 0; edge < 12; ++edge)
    do_test_side_number_1d( MBHEX, edge );
  for (int face = 0; face < 6; ++face)
    do_test_side_number_2d( MBHEX, face );
}


void test_opposite_side_tri()
{
  int idx, dim, err;
  err = CN::OppositeSide( MBTRI, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
  err = CN::OppositeSide( MBTRI, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBTRI, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBTRI, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBTRI, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBTRI, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );
}

void test_opposite_side_quad()
{
  int idx, dim, err;
  err = CN::OppositeSide( MBQUAD, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBQUAD, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 3, idx );
  err = CN::OppositeSide( MBQUAD, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBQUAD, 3, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );

  err = CN::OppositeSide( MBQUAD, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBQUAD, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 3, idx );
  err = CN::OppositeSide( MBQUAD, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBQUAD, 3, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
}

void test_opposite_side_tet()
{
  int idx, dim, err;
  
  err = CN::OppositeSide( MBTET, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 1, idx );
  err = CN::OppositeSide( MBTET, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBTET, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBTET, 3, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 3, idx );
  
  err = CN::OppositeSide( MBTET, 0, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBTET, 1, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBTET, 2, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );
  err = CN::OppositeSide( MBTET, 3, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 3, idx );
  
  err = CN::OppositeSide( MBTET, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 5, idx );
  err = CN::OppositeSide( MBTET, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 3, idx );
  err = CN::OppositeSide( MBTET, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 4, idx );
  err = CN::OppositeSide( MBTET, 3, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
  err = CN::OppositeSide( MBTET, 4, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBTET, 5, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );

}

void test_opposite_side_hex()
{
  int idx, dim, err;
  
  err = CN::OppositeSide( MBHEX, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 6, idx );
  err = CN::OppositeSide( MBHEX, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 7, idx );
  err = CN::OppositeSide( MBHEX, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 4, idx );
  err = CN::OppositeSide( MBHEX, 3, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 5, idx );
  err = CN::OppositeSide( MBHEX, 4, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBHEX, 5, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 3, idx );
  err = CN::OppositeSide( MBHEX, 6, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBHEX, 7, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );

  err = CN::OppositeSide( MBHEX, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 10, idx );
  err = CN::OppositeSide( MBHEX, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 11, idx );
  err = CN::OppositeSide( MBHEX, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 8, idx );
  err = CN::OppositeSide( MBHEX, 3, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 9, idx );
  err = CN::OppositeSide( MBHEX, 4, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 6, idx );
  err = CN::OppositeSide( MBHEX, 5, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 7, idx );
  err = CN::OppositeSide( MBHEX, 6, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 4, idx );
  err = CN::OppositeSide( MBHEX, 7, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 5, idx );
  err = CN::OppositeSide( MBHEX, 8, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBHEX, 9, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 3, idx );
  err = CN::OppositeSide( MBHEX, 10, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBHEX, 11, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
  
  err = CN::OppositeSide( MBHEX, 0, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 2, idx );
  err = CN::OppositeSide( MBHEX, 1, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 3, idx );
  err = CN::OppositeSide( MBHEX, 2, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 0, idx );
  err = CN::OppositeSide( MBHEX, 3, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 1, idx );
  err = CN::OppositeSide( MBHEX, 4, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 5, idx );
  err = CN::OppositeSide( MBHEX, 5, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 4, idx );
}

void test_has_mid_nodes(EntityType type)
{
  const int combinations[][4] = { { 0, 0, 0, 0 },
                                  { 0, 1, 0, 0 },
                                  { 0, 0, 1, 0 },
                                  { 0, 1, 1, 0 },
                                  { 0, 0, 0, 1 },
                                  { 0, 1, 0, 1 },
                                  { 0, 0, 1, 1 },
                                  { 0, 1, 1, 1 } };
  
    const int dim = CN::Dimension(type);
      // calculate number of valid combinations of ho node flags
    int num_comb = 1;
    for (int i = 0; i < dim; ++i)
      num_comb *= 2;
      // for each valid combination
    for (int c = 0; c < num_comb; ++c) {
        // calculate corresponding number of vertices in element
      const int* ho_nodes = combinations[c];
      int num_vtx = CN::VerticesPerEntity(type);
      switch (dim) {
        case 3: if (ho_nodes[2]) num_vtx += CN::NumSubEntities(type,2);
        case 2: if (ho_nodes[1]) num_vtx += CN::NumSubEntities(type,1);
      }
      if (ho_nodes[dim]) ++num_vtx;
      
      CHECK_EQUAL( ho_nodes[1], (int)CN::HasMidEdgeNodes( type, num_vtx ) );
      CHECK_EQUAL( ho_nodes[2], (int)CN::HasMidFaceNodes( type, num_vtx ) );
      CHECK_EQUAL( ho_nodes[3], (int)CN::HasMidRegionNodes( type, num_vtx ) );
      
      int results[4] = { 0, -1, -1, -1 };
      CN::HasMidNodes( type, num_vtx, results );
      CHECK_EQUAL(           0, !!results[0] );
      CHECK_EQUAL( ho_nodes[1], !!results[1] );
      CHECK_EQUAL( ho_nodes[2], !!results[2] );
      CHECK_EQUAL( ho_nodes[3], !!results[3] );
    }
}

void test_ho_node_parent()
{
  const int combinations[][4] = { { 0, 0, 0, 0 },
                                  { 0, 1, 0, 0 },
                                  { 0, 0, 1, 0 },
                                  { 0, 1, 1, 0 },
                                  { 0, 0, 0, 1 },
                                  { 0, 1, 0, 1 },
                                  { 0, 0, 1, 1 },
                                  { 0, 1, 1, 1 } };
  
  for (const EntityType* t = elem_types; *t != MBMAXTYPE; ++t) {
    const EntityType type = *t;
    const int dim = CN::Dimension(type);
      // calculate number of valid combinations of ho node flags
    int num_comb = 1;
    for (int i = 0; i < dim; ++i)
      num_comb *= 2;
      // for each valid combination
    for (int c = 0; c < num_comb; ++c) {
        // calculate corresponding number of vertices in element
      const int* ho_nodes = combinations[c];
      int num_vtx = CN::VerticesPerEntity(type);
      switch (dim) {
        case 3: if (ho_nodes[2]) num_vtx += CN::NumSubEntities(type,2);
        case 2: if (ho_nodes[1]) num_vtx += CN::NumSubEntities(type,1);
      }
      if (ho_nodes[dim]) ++num_vtx;

        // start at first higher-order node
      int pos = CN::VerticesPerEntity(type);
      
        // check mid-edge
      if (dim > 1 && ho_nodes[1]) {
        for (int i = 0; i < CN::NumSubEntities(type,1); ++i) {
          int pdim = -1, pidx = -1;
          CN::HONodeParent( type, num_vtx, pos++, pdim, pidx );
          CHECK_EQUAL( 1, pdim );
          CHECK_EQUAL( i, pidx );
        }
      }
      
        // check mid-face
      if (dim > 2 && ho_nodes[2]) {
        for (int i = 0; i < CN::NumSubEntities(type,2); ++i) {
          int pdim = -1, pidx = -1;
          CN::HONodeParent( type, num_vtx, pos++, pdim, pidx );
          CHECK_EQUAL( 2, pdim );
          CHECK_EQUAL( i, pidx );
        }
      }
      
        // check mid-volume
      if (ho_nodes[dim]) {
        int pdim = -1, pidx = -1;
        CN::HONodeParent( type, num_vtx, pos++, pdim, pidx );
        CHECK_EQUAL( dim, pdim );
        CHECK_EQUAL( 0, pidx );
      }
    } // for ho_node combinatinos
  } // for each type
}

void test_ho_node_index()
{
  const int combinations[][4] = { { 0, 0, 0, 0 },
                                  { 0, 1, 0, 0 },
                                  { 0, 0, 1, 0 },
                                  { 0, 1, 1, 0 },
                                  { 0, 0, 0, 1 },
                                  { 0, 1, 0, 1 },
                                  { 0, 0, 1, 1 },
                                  { 0, 1, 1, 1 } };
  
  for (const EntityType* t = elem_types; *t != MBMAXTYPE; ++t) {
    const EntityType type = *t;
    const int dim = CN::Dimension(type);
      // calculate number of valid combinations of ho node flags
    int num_comb = 1;
    for (int i = 0; i < dim; ++i)
      num_comb *= 2;
      // for each valid combination
    for (int c = 0; c < num_comb; ++c) {
        // calculate corresponding number of vertices in element
      const int* ho_nodes = combinations[c];
      int num_vtx = CN::VerticesPerEntity(type);
      switch (dim) {
        case 3: if (ho_nodes[2]) num_vtx += CN::NumSubEntities(type,2);
        case 2: if (ho_nodes[1]) num_vtx += CN::NumSubEntities(type,1);
      }
      if (ho_nodes[dim]) ++num_vtx;

        // start at first higher-order node
      int pos = CN::VerticesPerEntity(type);
      
        // check mid-edge
      if (dim > 1 && ho_nodes[1]) {
        for (int i = 0; i < CN::NumSubEntities(type,1); ++i) {
          int idx = CN::HONodeIndex( type, num_vtx, 1, i );
          CHECK_EQUAL( pos++, idx );
        }
      }
      
        // check mid-face
      if (dim > 2 && ho_nodes[2]) {
        for (int i = 0; i < CN::NumSubEntities(type,2); ++i) {
          int idx = CN::HONodeIndex( type, num_vtx, 2, i );
          CHECK_EQUAL( pos++, idx );
        }
      }
      
        // check mid-volume
      if (ho_nodes[dim]) {
        int idx = CN::HONodeIndex( type, num_vtx, dim, 0 );
        CHECK_EQUAL( pos++, idx );
      }
    } // for ho_node combinatinos
  } // for each type
}

void test_sub_entity_nodes( EntityType parent, int sub_dimension )
{
  const int num_corner = CN::VerticesPerEntity( parent );
  const int num_edge   = CN::NumSubEntities( parent, 1 );
  const int num_face   = CN::NumSubEntities( parent, 2 );
 
  switch (CN::Dimension(parent)) {
    case 3:
      test_sub_entity_nodes( parent, num_corner+num_face, sub_dimension );
      test_sub_entity_nodes( parent, num_corner+num_edge+num_face, sub_dimension );
      test_sub_entity_nodes( parent, num_corner+num_face+1, sub_dimension );
      test_sub_entity_nodes( parent, num_corner+num_edge+num_face+1, sub_dimension );
    case 2:
      test_sub_entity_nodes( parent, num_corner+num_edge, sub_dimension );
      test_sub_entity_nodes( parent, num_corner+num_edge+1, sub_dimension );
    case 1:
      test_sub_entity_nodes( parent, num_corner, sub_dimension );
      test_sub_entity_nodes( parent, num_corner+1, sub_dimension );
      break;
    default:
      CHECK(false);
  }
}

void test_sub_entity_nodes( EntityType parent, int num_nodes, int sub_dimension )
{
  const int num_sub = CN::NumSubEntities( parent, sub_dimension );
  const int parent_ho = CN::HasMidNodes( parent, num_nodes );
  int child_ho = 0;
  for (int d = 1; d <= sub_dimension; ++d)
    child_ho |= (parent_ho & (1<<d));
   
    // first test the types
  for (int i = 0; i < num_sub; ++i) {
    int num, conn[moab::MAX_SUB_ENTITY_VERTICES];
    EntityType type;
    CN::SubEntityNodeIndices( parent, num_nodes, sub_dimension, i, type, num, conn );
    CHECK_EQUAL( CN::SubEntityType(parent, sub_dimension, i), type );
  }
 
    // now test that they have the correct number of higher-order node
  for (int i = 0; i < num_sub; ++i) {
    int num, conn[moab::MAX_SUB_ENTITY_VERTICES];
    EntityType type;
    CN::SubEntityNodeIndices( parent, num_nodes, sub_dimension, i, type, num, conn );
    const int ho = CN::HasMidNodes( type, num );
    CHECK_EQUAL( child_ho, ho );
  }
  
    // now test the actual indices
  for (int i = 0; i < num_sub; ++i) {
    int num, conn[moab::MAX_SUB_ENTITY_VERTICES], corners[moab::MAX_SUB_ENTITY_VERTICES];
    EntityType type;
    CN::SubEntityNodeIndices( parent, num_nodes, sub_dimension, i, type, num, conn );
    
      // check corner indices against SubEntityVertexIndices
    const int num_corner = CN::VerticesPerEntity(type);
    CHECK( num >= num_corner );
    CN::SubEntityVertexIndices( parent, sub_dimension, i, corners );
    for (int j = 0; j < num_corner; ++j)
      CHECK_EQUAL( corners[j], conn[j] );
    
      // check mid-edge indices, if present
    int idx = num_corner;
    if (child_ho & CN::MID_EDGE_BIT) {
        // for each edge in the sub-entity type
      const int num_edge = CN::NumSubEntities( type, 1 );
      for (int j = 0; j < num_edge; ++j) {
          // get edge indices for sub-entity connectivity
        int edge_ends[2];
        CN::SubEntityVertexIndices( type, 1, j, edge_ends );
          // convert to indices into parent type's connectivity
        CHECK( edge_ends[0] < num_corner );
        edge_ends[0] = corners[edge_ends[0]];
        CHECK( edge_ends[1] < num_corner );
        edge_ends[1] = corners[edge_ends[1]];
          // find edge index in parent element
        int side, sense, off;
        int result = CN::SideNumber( parent, edge_ends, 2, 1, side, sense, off );
        CHECK_EQUAL( 0, result );
          // get location in parent entity connectivity for mid-edge node
        int loc = CN::HONodeIndex( parent, num_nodes, 1, side );
        CHECK_EQUAL( loc, conn[idx++] );
      }
    }
    
      // check mid-face indices, if present
    if (child_ho & CN::MID_FACE_BIT) {
      CHECK_EQUAL( 2, CN::Dimension(type) );
      int loc = CN::HONodeIndex( parent, num_nodes, 2, i );
      CHECK_EQUAL( loc, conn[idx++] );
    }
    
      // make sure there were no extra node indices returned
    CHECK_EQUAL( idx, num );
  }
}
