#include "TestUtil.hpp"
#include "MBCN.hpp"

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

void test_has_mid_nodes();
void test_ho_node_parent();
void test_ho_node_index();

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
  
  result += RUN_TEST(test_has_mid_nodes);
  result += RUN_TEST(test_ho_node_parent);
  result += RUN_TEST(test_ho_node_index);
  return result;
}

const MBEntityType elem_types[] = { MBTRI,
                                    MBQUAD,
                                    MBTET,
                                    MBPYRAMID,
                                    MBPRISM,
                                    MBHEX,
                                    MBMAXTYPE };
const int num_elem_types = sizeof(elem_types)/sizeof(elem_types[0]) - 1;

void test_dimension_pair()
{
  MBDimensionPair dp;
  
  dp = MBCN::TypeDimensionMap[0];
  CHECK_EQUAL( MBVERTEX, dp.first );  
  CHECK_EQUAL( MBVERTEX, dp.second );
  
  dp = MBCN::TypeDimensionMap[1];
  CHECK_EQUAL( MBEDGE, dp.first );  
  CHECK_EQUAL( MBEDGE, dp.second );
  
  dp = MBCN::TypeDimensionMap[2];
  CHECK_EQUAL( MBTRI, dp.first );  
  CHECK_EQUAL( MBPOLYGON, dp.second );
  
  dp = MBCN::TypeDimensionMap[3];
  CHECK_EQUAL( MBTET, dp.first );  
  CHECK_EQUAL( MBPOLYHEDRON, dp.second );
}

void test_type_names()
{
  for (MBEntityType t = MBVERTEX; t != MBMAXTYPE; ++t) {
    const char* name = MBCN::EntityTypeName(t);
    CHECK_EQUAL( t, MBCN::EntityTypeFromName(name) );
  }
}

void test_dimension()
{
  CHECK_EQUAL( 0, MBCN::Dimension(MBVERTEX) );
  CHECK_EQUAL( 1, MBCN::Dimension(MBEDGE) );
  CHECK_EQUAL( 2, MBCN::Dimension(MBTRI) );
  CHECK_EQUAL( 2, MBCN::Dimension(MBQUAD) );
  CHECK_EQUAL( 2, MBCN::Dimension(MBPOLYGON) );
  CHECK_EQUAL( 3, MBCN::Dimension(MBTET) );
  CHECK_EQUAL( 3, MBCN::Dimension(MBPYRAMID) );
  CHECK_EQUAL( 3, MBCN::Dimension(MBPRISM) );
  CHECK_EQUAL( 3, MBCN::Dimension(MBKNIFE) );
  CHECK_EQUAL( 3, MBCN::Dimension(MBHEX) );
  CHECK_EQUAL( 3, MBCN::Dimension(MBPOLYHEDRON) );
}

void test_vertices_per_entity()
{
  CHECK_EQUAL( 1, MBCN::VerticesPerEntity(MBVERTEX) );
  CHECK_EQUAL( 2, MBCN::VerticesPerEntity(MBEDGE) );
  CHECK_EQUAL( 3, MBCN::VerticesPerEntity(MBTRI) );
  CHECK_EQUAL( 4, MBCN::VerticesPerEntity(MBQUAD) );
  CHECK_EQUAL( 4, MBCN::VerticesPerEntity(MBTET) );
  CHECK_EQUAL( 5, MBCN::VerticesPerEntity(MBPYRAMID) );
  CHECK_EQUAL( 6, MBCN::VerticesPerEntity(MBPRISM) );
  CHECK_EQUAL( 7, MBCN::VerticesPerEntity(MBKNIFE) );
  CHECK_EQUAL( 8, MBCN::VerticesPerEntity(MBHEX) );
}

void test_num_sub_entities()
{
  CHECK_EQUAL( 1, MBCN::NumSubEntities(MBVERTEX, 0));

  CHECK_EQUAL( 2, MBCN::NumSubEntities(MBEDGE, 0));
  CHECK_EQUAL( 1, MBCN::NumSubEntities(MBEDGE, 1));

  CHECK_EQUAL( 3, MBCN::NumSubEntities(MBTRI, 0));
  CHECK_EQUAL( 3, MBCN::NumSubEntities(MBTRI, 1));
  CHECK_EQUAL( 1, MBCN::NumSubEntities(MBTRI, 2));

  CHECK_EQUAL( 4, MBCN::NumSubEntities(MBQUAD, 0));
  CHECK_EQUAL( 4, MBCN::NumSubEntities(MBQUAD, 1));
  CHECK_EQUAL( 1, MBCN::NumSubEntities(MBQUAD, 2));

  CHECK_EQUAL( 4, MBCN::NumSubEntities(MBTET, 0));
  CHECK_EQUAL( 6, MBCN::NumSubEntities(MBTET, 1));
  CHECK_EQUAL( 4, MBCN::NumSubEntities(MBTET, 2));

  CHECK_EQUAL( 5, MBCN::NumSubEntities(MBPYRAMID, 0));
  CHECK_EQUAL( 8, MBCN::NumSubEntities(MBPYRAMID, 1));
  CHECK_EQUAL( 5, MBCN::NumSubEntities(MBPYRAMID, 2));

  CHECK_EQUAL( 6, MBCN::NumSubEntities(MBPRISM, 0));
  CHECK_EQUAL( 9, MBCN::NumSubEntities(MBPRISM, 1));
  CHECK_EQUAL( 5, MBCN::NumSubEntities(MBPRISM, 2));

  CHECK_EQUAL( 7, MBCN::NumSubEntities(MBKNIFE, 0));
  CHECK_EQUAL( 8, MBCN::NumSubEntities(MBKNIFE, 1));
  CHECK_EQUAL( 5, MBCN::NumSubEntities(MBKNIFE, 2));

  CHECK_EQUAL( 8, MBCN::NumSubEntities(MBHEX, 0));
  CHECK_EQUAL( 12, MBCN::NumSubEntities(MBHEX, 1));
  CHECK_EQUAL( 6, MBCN::NumSubEntities(MBHEX, 2));
}

void do_test_sub_entity_type_2d( MBEntityType type )
{
  for (int j = 0; j < MBCN::VerticesPerEntity(type); ++j) {
    CHECK_EQUAL( MBVERTEX, MBCN::SubEntityType(type, 0, j ) );
    CHECK_EQUAL( MBEDGE,   MBCN::SubEntityType(type, 1, j ) );
  }
  CHECK_EQUAL( type, MBCN::SubEntityType(type, 2, 0) );
}

void do_test_sub_entity_type_3d( MBEntityType type,
                                 int num_faces,
                                 const MBEntityType* face_types )
{
  for (int j = 0; j < MBCN::VerticesPerEntity(type); ++j) {
    CHECK_EQUAL( MBVERTEX, MBCN::SubEntityType(type, 0, j ) );
  }

  for (int j = 0; j < MBCN::NumSubEntities(type,1); ++j) {
    CHECK_EQUAL( MBEDGE, MBCN::SubEntityType(type, 1, j ) );
  }

  for (int j = 0; j < num_faces; ++j) {
    MBEntityType sub_type = MBCN::SubEntityType( type, 2, j );
    CHECK_EQUAL( face_types[j], sub_type );
  }

  CHECK_EQUAL( type, MBCN::SubEntityType(type, 3, 0) );
}

void test_sub_entity_type_vtx()
{
  CHECK_EQUAL( MBVERTEX, MBCN::SubEntityType(MBVERTEX, 0, 0 ));
}

void test_sub_entity_type_edge()
{
  CHECK_EQUAL( MBVERTEX, MBCN::SubEntityType(MBEDGE, 0, 0 ));
  CHECK_EQUAL( MBVERTEX, MBCN::SubEntityType(MBEDGE, 0, 1 ));
  CHECK_EQUAL( MBEDGE  , MBCN::SubEntityType(MBEDGE, 1, 0 ));
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
  const MBEntityType types[] = { MBTRI, MBTRI, MBTRI, MBTRI };
  do_test_sub_entity_type_3d( MBTET, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_pyr()
{
  const MBEntityType types[] = { MBQUAD, MBTRI, MBTRI, MBTRI, MBTRI };
  do_test_sub_entity_type_3d( MBPYRAMID, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_pri()
{
  const MBEntityType types[] = { MBQUAD, MBQUAD, MBQUAD, MBTRI, MBTRI };
  do_test_sub_entity_type_3d( MBPRISM, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_knife()
{
  const MBEntityType types[] = { MBQUAD, MBQUAD, MBQUAD, MBQUAD, MBQUAD };
  do_test_sub_entity_type_3d( MBKNIFE, sizeof(types)/sizeof(types[0]), types );
}

void test_sub_entity_type_hex()
{
  const MBEntityType types[] = { MBQUAD, MBQUAD, MBQUAD, MBQUAD, MBQUAD, MBQUAD };
  do_test_sub_entity_type_3d( MBHEX, sizeof(types)/sizeof(types[0]), types );
}


void test_0d_sub_entity_indices( MBEntityType type, int num_vtx )
{
  for (int i = 0; i < num_vtx; ++i) {
    // zero input array
    int indices[2] = { 0, -100 };
    // check correct results
    MBCN::SubEntityVertexIndices( type, 0, i, indices );
    CHECK_EQUAL( i, indices[0] );
    // check didn't write past end of array
    CHECK_EQUAL( -100, indices[1] );
  }
}

void test_1d_sub_entity_indices( MBEntityType type, int num_edges, 
                                 const int (*edge_indices)[2] )
{
  for (int i = 0; i < num_edges; ++i) {
    // zero input array
    int indices[3] = { 0, 0, -99 };
    // check correct results
    MBCN::SubEntityVertexIndices( type, 1, i, indices );
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

void test_2d_sub_entity_indices( MBEntityType type, int num_faces,
                                 const int (*face_indices)[5] )
{
  for (int i = 0; i < num_faces; ++i) {
    // zero input array
    int indices[5] = { 0, 0, 0, -99, -99 };
    // check correct results
    MBCN::SubEntityVertexIndices( type, 2, i, indices );
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

void test_elem_as_sub_entity( MBEntityType type, int dim, int num_vertices )
{ 
  int indices[9] = { -2, -2, -2, -2, -2, -2, -2, -2, -2 };
  MBCN::SubEntityVertexIndices( type, dim, 0, indices );
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
  const MBEntityType type = MBTET;
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
  const MBEntityType type = MBPYRAMID;
  const int num_vtx = 5;
  const int edges[][2] = { { 0, 1 },
                           { 1, 2 },
                           { 2, 3 },
                           { 3, 0 },
                           { 0, 4 },
                           { 1, 4 },
                           { 2, 4 },
                           { 3, 4 } };
  const int faces[][5] = { { 4, 3, 2, 1, 0 },
                           { 3, 0, 1, 4, 0 },
                           { 3, 1, 2, 4, 0 },
                           { 3, 2, 3, 4, 0 },
                           { 3, 3, 0, 4, 0 } };
  test_0d_sub_entity_indices( type, num_vtx );
  test_1d_sub_entity_indices( type, sizeof(edges)/sizeof(edges[0]), edges );
  test_2d_sub_entity_indices( type, sizeof(faces)/sizeof(faces[0]), faces );
  test_elem_as_sub_entity( type, 3, num_vtx );
}

void test_sub_entity_indices_pri()
{
  const MBEntityType type = MBPRISM;
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
  const MBEntityType type = MBHEX;
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

static void do_test_side_number_1d( MBEntityType type, int idx )
{
  // define a random handle list
  const int elem_verts[] = { 7400, 6233, 3027, 0454, 6839, 5391, 7735, 3603 };
  // get side indices
  int side_idx[4] = { 0, 0 };
  MBCN::SubEntityVertexIndices( type, 1, idx, side_idx );

  // "reversed" and "offset" are the same thing for edges.
  int side_conn[2] = { elem_verts[side_idx[0]], elem_verts[side_idx[1]] };
  int rev_conn[2] = { elem_verts[side_idx[1]], elem_verts[side_idx[0]] };
  int result_side = -100, result_sense = -100, result_offset = -100;
  int err= MBCN::SideNumber( type, elem_verts, side_conn, 2, 1, 
                             result_side, result_sense, result_offset );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( idx, result_side );
  CHECK_EQUAL( 1, result_sense );
  CHECK_EQUAL( 0, result_offset);
  err= MBCN::SideNumber( type, elem_verts, rev_conn, 2, 1, 
                         result_side, result_sense, result_offset );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( idx, result_side );
  CHECK(result_offset == 1 || result_sense == -1);
}     

static void do_test_side_number_2d( MBEntityType type, int idx )
{
  // define a random handle list
  const int elem_verts[] = { 7400, 6233, 3027, 0454, 6839, 5391, 7735, 3603 };
  // get side indices
  const int side_size = MBCN::VerticesPerEntity( MBCN::SubEntityType( type, 2, idx ) );
  int side_idx[4] = { 0, 0, 0, 0 };
  MBCN::SubEntityVertexIndices( type, 2, idx, side_idx );

  // for each possible forward or reverse offset
  for (int rev = -1; rev < 2; rev += 2) {
    for (int off = 0; off < side_size; ++off) {
      int side_conn[4]; side_conn[3] = 0;
      for (int i = 0; i < side_size; ++i)
        side_conn[(i+side_size-rev*off)%side_size] = elem_verts[side_idx[i]];
      
      int result_side = -100, result_sense = -100, result_offset = -100;
      int err = MBCN::SideNumber( type, elem_verts, side_conn, side_size, 2, 
                                  result_side, result_sense, result_offset );
      CHECK_EQUAL( 0, err );
      CHECK_EQUAL( idx, result_side );
      CHECK_EQUAL( 1-2*rev, result_sense );
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
    do_test_side_number_1d( MBTET, face );
}

void test_side_number_pyr()
{
  for (int edge = 0; edge < 8; ++edge)
    do_test_side_number_1d( MBPYRAMID, edge );
  for (int face = 0; face < 5; ++face)
    do_test_side_number_1d( MBPYRAMID, face );
}

void test_side_number_pri()
{
  for (int edge = 0; edge < 9; ++edge)
    do_test_side_number_1d( MBPRISM, edge );
  for (int face = 0; face < 5; ++face)
    do_test_side_number_1d( MBPRISM, face );
}

void test_side_number_hex()
{
  for (int edge = 0; edge < 12; ++edge)
    do_test_side_number_1d( MBHEX, edge );
  for (int face = 0; face < 8; ++face)
    do_test_side_number_1d( MBHEX, face );
}


void test_opposite_side_tri()
{
  int idx, dim, err;
  err = MBCN::OppositeSide( MBTRI, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
  err = MBCN::OppositeSide( MBTRI, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBTRI, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBTRI, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBTRI, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBTRI, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );
}

void test_opposite_side_quad()
{
  int idx, dim, err;
  err = MBCN::OppositeSide( MBQUAD, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBQUAD, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 3, idx );
  err = MBCN::OppositeSide( MBQUAD, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBQUAD, 3, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );

  err = MBCN::OppositeSide( MBQUAD, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBQUAD, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 3, idx );
  err = MBCN::OppositeSide( MBQUAD, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBQUAD, 3, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
}

void test_opposite_side_tet()
{
  int idx, dim, err;
  
  err = MBCN::OppositeSide( MBTET, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 1, idx );
  err = MBCN::OppositeSide( MBTET, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBTET, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBTET, 3, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 3, idx );
  
  err = MBCN::OppositeSide( MBTET, 0, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBTET, 1, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBTET, 2, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );
  err = MBCN::OppositeSide( MBTET, 3, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 3, idx );
  
  err = MBCN::OppositeSide( MBTET, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 5, idx );
  err = MBCN::OppositeSide( MBTET, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 3, idx );
  err = MBCN::OppositeSide( MBTET, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 4, idx );
  err = MBCN::OppositeSide( MBTET, 3, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
  err = MBCN::OppositeSide( MBTET, 4, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBTET, 5, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );

}

void test_opposite_side_hex()
{
  int idx, dim, err;
  
  err = MBCN::OppositeSide( MBHEX, 0, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 6, idx );
  err = MBCN::OppositeSide( MBHEX, 1, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 7, idx );
  err = MBCN::OppositeSide( MBHEX, 2, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 4, idx );
  err = MBCN::OppositeSide( MBHEX, 3, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 5, idx );
  err = MBCN::OppositeSide( MBHEX, 4, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBHEX, 5, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 3, idx );
  err = MBCN::OppositeSide( MBHEX, 6, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBHEX, 7, 0, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 0, dim );
  CHECK_EQUAL( 1, idx );

  err = MBCN::OppositeSide( MBHEX, 0, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 10, idx );
  err = MBCN::OppositeSide( MBHEX, 1, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 11, idx );
  err = MBCN::OppositeSide( MBHEX, 2, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 8, idx );
  err = MBCN::OppositeSide( MBHEX, 3, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 9, idx );
  err = MBCN::OppositeSide( MBHEX, 4, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 6, idx );
  err = MBCN::OppositeSide( MBHEX, 5, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 7, idx );
  err = MBCN::OppositeSide( MBHEX, 6, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 4, idx );
  err = MBCN::OppositeSide( MBHEX, 7, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 5, idx );
  err = MBCN::OppositeSide( MBHEX, 8, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBHEX, 9, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 3, idx );
  err = MBCN::OppositeSide( MBHEX, 10, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBHEX, 11, 1, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 1, dim );
  CHECK_EQUAL( 1, idx );
  
  err = MBCN::OppositeSide( MBHEX, 0, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 2, idx );
  err = MBCN::OppositeSide( MBHEX, 1, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 3, idx );
  err = MBCN::OppositeSide( MBHEX, 2, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 0, idx );
  err = MBCN::OppositeSide( MBHEX, 3, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 1, idx );
  err = MBCN::OppositeSide( MBHEX, 4, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 5, idx );
  err = MBCN::OppositeSide( MBHEX, 5, 2, idx, dim );
  CHECK_EQUAL( 0, err );
  CHECK_EQUAL( 2, dim );
  CHECK_EQUAL( 4, idx );
}

void test_has_mid_nodes()
{
  const int combinations[][4] = { { 0, 0, 0, 0 },
                                  { 0, 1, 0, 0 },
                                  { 0, 0, 1, 0 },
                                  { 0, 1, 1, 0 },
                                  { 0, 0, 0, 1 },
                                  { 0, 1, 0, 1 },
                                  { 0, 0, 1, 1 },
                                  { 0, 1, 1, 1 } };
  
  for (const MBEntityType* t = elem_types; *t != MBMAXTYPE; ++t) {
    const MBEntityType type = *t;
    const int dim = MBCN::Dimension(type);
      // calculate number of valid combinations of ho node flags
    int num_comb = 1;
    for (int i = 0; i < dim; ++i)
      num_comb *= 2;
      // for each valid combination
    for (int c = 0; c < num_comb; ++c) {
        // calculate corresponding number of vertices in element
      const int* ho_nodes = combinations[c];
      int num_vtx = MBCN::VerticesPerEntity(type);
      switch (dim) {
        case 3: if (ho_nodes[2]) num_vtx += MBCN::NumSubEntities(type,2);
        case 2: if (ho_nodes[1]) num_vtx += MBCN::NumSubEntities(type,1);
      }
      if (ho_nodes[dim]) ++num_vtx;
      
      CHECK_EQUAL( ho_nodes[1], (int)MBCN::HasMidEdgeNodes( type, num_vtx ) );
      CHECK_EQUAL( ho_nodes[2], (int)MBCN::HasMidFaceNodes( type, num_vtx ) );
      CHECK_EQUAL( ho_nodes[3], (int)MBCN::HasMidRegionNodes( type, num_vtx ) );
      
      int results[4] = { 0, -1, -1, -1 };
      MBCN::HasMidNodes( type, num_vtx, results );
      CHECK_EQUAL(           0, results[0] );
      CHECK_EQUAL( ho_nodes[1], results[1] );
      CHECK_EQUAL( ho_nodes[2], results[2] );
      CHECK_EQUAL( ho_nodes[3], results[3] );
    }
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
  
  for (const MBEntityType* t = elem_types; *t != MBMAXTYPE; ++t) {
    const MBEntityType type = *t;
    const int dim = MBCN::Dimension(type);
      // calculate number of valid combinations of ho node flags
    int num_comb = 1;
    for (int i = 0; i < dim; ++i)
      num_comb *= 2;
      // for each valid combination
    for (int c = 0; c < num_comb; ++c) {
        // calculate corresponding number of vertices in element
      const int* ho_nodes = combinations[c];
      int num_vtx = MBCN::VerticesPerEntity(type);
      switch (dim) {
        case 3: if (ho_nodes[2]) num_vtx += MBCN::NumSubEntities(type,2);
        case 2: if (ho_nodes[1]) num_vtx += MBCN::NumSubEntities(type,1);
      }
      if (ho_nodes[dim]) ++num_vtx;

        // start at first higher-order node
      int pos = MBCN::VerticesPerEntity(type);
      
        // check mid-edge
      if (dim > 1 && ho_nodes[1]) {
        for (int i = 0; i < MBCN::NumSubEntities(type,1); ++i) {
          int pdim = -1, pidx = -1;
          MBCN::HONodeParent( type, num_vtx, pos++, pdim, pidx );
          CHECK_EQUAL( 1, pdim );
          CHECK_EQUAL( i, pidx );
        }
      }
      
        // check mid-face
      if (dim > 2 && ho_nodes[2]) {
        for (int i = 0; i < MBCN::NumSubEntities(type,2); ++i) {
          int pdim = -1, pidx = -1;
          MBCN::HONodeParent( type, num_vtx, pos++, pdim, pidx );
          CHECK_EQUAL( 2, pdim );
          CHECK_EQUAL( i, pidx );
        }
      }
      
        // check mid-volume
      if (ho_nodes[dim]) {
        int pdim = -1, pidx = -1;
        MBCN::HONodeParent( type, num_vtx, pos++, pdim, pidx );
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
  
  for (const MBEntityType* t = elem_types; *t != MBMAXTYPE; ++t) {
    const MBEntityType type = *t;
    const int dim = MBCN::Dimension(type);
      // calculate number of valid combinations of ho node flags
    int num_comb = 1;
    for (int i = 0; i < dim; ++i)
      num_comb *= 2;
      // for each valid combination
    for (int c = 0; c < num_comb; ++c) {
        // calculate corresponding number of vertices in element
      const int* ho_nodes = combinations[c];
      int num_vtx = MBCN::VerticesPerEntity(type);
      switch (dim) {
        case 3: if (ho_nodes[2]) num_vtx += MBCN::NumSubEntities(type,2);
        case 2: if (ho_nodes[1]) num_vtx += MBCN::NumSubEntities(type,1);
      }
      if (ho_nodes[dim]) ++num_vtx;

        // start at first higher-order node
      int pos = MBCN::VerticesPerEntity(type);
      
        // check mid-edge
      if (dim > 1 && ho_nodes[1]) {
        for (int i = 0; i < MBCN::NumSubEntities(type,1); ++i) {
          int idx = MBCN::HONodeIndex( type, num_vtx, 1, i );
          CHECK_EQUAL( pos++, idx );
        }
      }
      
        // check mid-face
      if (dim > 2 && ho_nodes[2]) {
        for (int i = 0; i < MBCN::NumSubEntities(type,2); ++i) {
          int idx = MBCN::HONodeIndex( type, num_vtx, 2, i );
          CHECK_EQUAL( pos++, idx );
        }
      }
      
        // check mid-volume
      if (ho_nodes[dim]) {
        int idx = MBCN::HONodeIndex( type, num_vtx, dim, 0 );
        CHECK_EQUAL( pos++, idx );
      }
    } // for ho_node combinatinos
  } // for each type
}

