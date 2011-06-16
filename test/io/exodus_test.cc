#include "TestUtil.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include "moab/Range.hpp"
#include "ReadNCDF.hpp"
#include "FileOptions.hpp"
#define IS_BUILDING_MB
#include "ExoIIUtil.hpp"
#include <math.h>
#include <algorithm>

using namespace moab;

/* Input test file: ho_test.g
 * 
 * File is expected to contain at least one block for every
 * supported higher-order element type.  The coordinates of
 * every higher-order node are expected to be the mean of the
 * adjacent corner vertices of the element.
 */
#ifdef MESHDIR
static const char ho_file[] = STRINGIFY(MESHDIR) "/io/ho_test.g";
static const char file_one[] = STRINGIFY(MESHDIR) "/mbtest1.g";
static const char alt_file[] = STRINGIFY(MESHDIR) "/io/hex_2x2x2_ss.exo";
#else
static const char ho_file[] = "ho_test.g";
static const char file_one[] = "mbtest1.g";
static const char alt_one[] = "hex_2x2x2_ss.exo";
#endif

void read_file( Interface& moab, 
                const char* input_file );

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( Interface& moab, 
                       EntityHandle entity,
                       int mid_nodes[4] );

void test_read_side( int sideset_id,
                     EntityType sideset_type,
                     int sideset_nodes_per_elem,
                     bool shell_side = false );

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( EntityType type, int num_nodes );

// Tests originally in MBTest.cpp
void mb_vertex_coordinate_test();
void mb_bar_connectivity_test();
void mb_tri_connectivity_test();
void mb_quad_connectivity_test();
void mb_hex_connectivity_test();
void mb_tet_connectivity_test();
void mb_write_mesh_test();

void test_types();

void test_tri6 () { test_ho_elements(MBTRI, 6); }
void test_tri7 () { test_ho_elements(MBTRI, 7); }

void test_quad5() { test_ho_elements(MBQUAD, 5); }
void test_quad8() { test_ho_elements(MBQUAD, 8); }
void test_quad9() { test_ho_elements(MBQUAD, 9); }

void test_tet8 () { test_ho_elements(MBTET,  8); }
void test_tet10() { test_ho_elements(MBTET, 10); }
void test_tet14() { test_ho_elements(MBTET, 14); }

void test_hex9 () { test_ho_elements(MBHEX,  9); }
void test_hex20() { test_ho_elements(MBHEX, 20); }
void test_hex27() { test_ho_elements(MBHEX, 27); }

void test_read_tri6_side()    { test_read_side( 1, MBEDGE, 3 ); }  // sideset 1
void test_read_shell_side()   { test_read_side( 3, MBQUAD, 9, true ); } // sideset 3
void test_read_shell_edge()   { test_read_side( 4, MBEDGE, 3 ); } // sideset 4
void test_read_hex20_side()   { test_read_side( 2, MBQUAD, 8 ); }  // sideset 2

void test_read_block_ids();
void test_read_sideset_ids();
void test_read_nodeset_ids();

void test_read_alternate_coord_format();

int main()
{
  int result = 0;
  
  result += RUN_TEST(mb_vertex_coordinate_test);
  result += RUN_TEST(mb_bar_connectivity_test);
  result += RUN_TEST(mb_tri_connectivity_test);
  result += RUN_TEST(mb_quad_connectivity_test);
  result += RUN_TEST(mb_hex_connectivity_test);
  result += RUN_TEST(mb_tet_connectivity_test);
  result += RUN_TEST(mb_write_mesh_test);
  
  result += RUN_TEST(test_types);
  
  result += RUN_TEST(test_tri6 );
  result += RUN_TEST(test_tri7 );
  result += RUN_TEST(test_quad5);
  result += RUN_TEST(test_quad8);
  result += RUN_TEST(test_quad9);
  result += RUN_TEST(test_tet8 );
  result += RUN_TEST(test_tet10);
  result += RUN_TEST(test_tet14);
  result += RUN_TEST(test_hex9 );
  result += RUN_TEST(test_hex20);
  result += RUN_TEST(test_hex27);
  
  result += RUN_TEST(test_read_tri6_side );
  result += RUN_TEST(test_read_shell_side);
  result += RUN_TEST(test_read_shell_edge);
  result += RUN_TEST(test_read_hex20_side);
  
  result += RUN_TEST(test_read_block_ids );
  result += RUN_TEST(test_read_sideset_ids);
  result += RUN_TEST(test_read_nodeset_ids);

  result += RUN_TEST(test_read_alternate_coord_format);
  
  return result;
}

void load_file_one( Interface* iface )
{
  ErrorCode error = iface->load_mesh( file_one );
  if (MB_SUCCESS != error) {
    std::cout << "Failed to load input file: " << file_one << std::endl;
    std::string error_reason;
    iface->get_last_error(error_reason);
    std::cout << error_reason << std::endl;
  }
  CHECK_ERR(error);
}

  /*!
    @test 
    Vertex Coordinates
    @li Get coordinates of vertex 1 correctly
    @li Get coordinates of vertex 8 correctly
    @li Get coordinates of vertex 6 correctly
  */
void mb_vertex_coordinate_test()
{
  double coords[3];
  EntityHandle handle;
  ErrorCode error;
  int err;

  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );

    // coordinate 2 should be {1.5, -1.5, 3.5}

  handle = CREATE_HANDLE(MBVERTEX, 2, err);
  error = MB->get_coords(&handle, 1, coords );
  CHECK_ERR(error);
  const double exp2[] = { 1.5, -1.5, 3.5 };
  CHECK_ARRAYS_EQUAL( exp2, 3, coords, 3 );


    // coordinate 9 should be {1, -2, 3.5}
  handle = CREATE_HANDLE(MBVERTEX, 9, err);
  error = MB->get_coords(&handle, 1, coords );
  CHECK_ERR(error);
  const double exp9[] = { 1, -2, 3.5 };
  CHECK_ARRAYS_EQUAL( exp9, 3, coords, 3 );

    // coordinate 7 should be {0.5, -2, 3.5}
  handle = CREATE_HANDLE(MBVERTEX, 7, err);
  error = MB->get_coords(&handle, 1, coords);
  CHECK_ERR(error);
  const double exp7[] = {0.5, -2, 3.5};
  CHECK_ARRAYS_EQUAL( exp7, 3, coords, 3 );
  
  int node_count = 0;
  error = MB->get_number_entities_by_type( 0, MBVERTEX, node_count );
  CHECK_ERR(error);
    // Number of vertices (node_count) should be 83 assuming no gaps in the handle space
  CHECK_EQUAL( 47, node_count );
}


  /*!
    @test
    MB Bar Element Connectivity Test
    @li Get coordinates for 2 node bar elements
  */

void mb_bar_connectivity_test()
{
  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );

  std::vector<EntityHandle> conn;
  Range bars;

  ErrorCode error = MB->get_entities_by_type(0, MBEDGE, bars);
  CHECK_ERR(error);

    // get the connectivity of the second bar
  EntityHandle handle = *(++bars.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  CHECK_ERR(error);

  CHECK_EQUAL( (size_t)2, conn.size() );

    // from ncdump the connectivity of bar 2 (0 based) is
    //  14, 13 
  CHECK_EQUAL( (EntityHandle)14, conn[0] );
  CHECK_EQUAL( (EntityHandle)13, conn[1] );
  
    // Now try getting the connectivity of one of the vertices for fun.
    // just return the vertex in the connectivity
  handle = conn[0];
  error = MB->get_connectivity(&handle, 1, conn);
  CHECK_EQUAL( MB_FAILURE, error );
}

void mb_tri_connectivity_test()
{
  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );

  std::vector<EntityHandle> conn; 
  Range tris;
  ErrorCode error = MB->get_entities_by_type(0, MBTRI, tris);
  CHECK_ERR(error);

    // get the connectivity of the second tri
  EntityHandle handle = *(++tris.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  CHECK_ERR(error);

  CHECK_EQUAL( (size_t)3, conn.size() );

    // from ncdump the connectivity of tri 2 (0 based) is
    //  45, 37, 38

  CHECK_EQUAL( (EntityHandle)45, conn[0] );
  CHECK_EQUAL( (EntityHandle)37, conn[1] );
  CHECK_EQUAL( (EntityHandle)38, conn[2] );
}

void mb_quad_connectivity_test()
{
  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );

  std::vector<EntityHandle> conn;
  Range quads;

  ErrorCode error = MB->get_entities_by_type(0, MBQUAD, quads);
  CHECK_ERR(error);

    // get the connectivity of the second quad
  EntityHandle handle = *(++quads.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  CHECK_ERR(error);

  CHECK_EQUAL( (size_t)4, conn.size() );

    // from ncdump the connectivity of quad 2 (0 based) is
    // 20, 11, 12, 26,

  CHECK_EQUAL( (EntityHandle)20, conn[0] );
  CHECK_EQUAL( (EntityHandle)11, conn[1] );
  CHECK_EQUAL( (EntityHandle)12, conn[2] );
  CHECK_EQUAL( (EntityHandle)26, conn[3] );
}

void mb_hex_connectivity_test()
{
  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );

  std::vector<EntityHandle> conn;
  Range hexes;

  ErrorCode error = MB->get_entities_by_type(0,  MBHEX, hexes);
  CHECK_ERR(error);

    // get the connectivity of the second hex
  EntityHandle handle = *(++hexes.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  CHECK_ERR(error);

  CHECK_EQUAL( (size_t)8, conn.size() );

    // from ncdump the connectivity of hex 1 (0 based) is
    //19, 13, 16, 23, 21, 14, 18, 27

  CHECK_EQUAL( (EntityHandle)19, conn[0] );
  CHECK_EQUAL( (EntityHandle)13, conn[1] );
  CHECK_EQUAL( (EntityHandle)16, conn[2] );
  CHECK_EQUAL( (EntityHandle)23, conn[3] );
  CHECK_EQUAL( (EntityHandle)21, conn[4] );
  CHECK_EQUAL( (EntityHandle)14, conn[5] );
  CHECK_EQUAL( (EntityHandle)18, conn[6] );
  CHECK_EQUAL( (EntityHandle)27, conn[7] );
}

void mb_tet_connectivity_test()
{
  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );

  std::vector<EntityHandle> conn; 
  Range tets;
  ErrorCode error = MB->get_entities_by_type(0, MBTET, tets);
  CHECK_ERR(error);

    // get the connectivity of the second tet
  EntityHandle handle = *(++tets.begin());

  error = MB->get_connectivity(&handle, 1, conn);
  CHECK_ERR(error);

  CHECK_EQUAL( (size_t)4, conn.size() );

    // from ncdump the connectivity of tet 2 (0 based) is: 
    // 35, 34, 32, 43 

  CHECK_EQUAL( (EntityHandle)35, conn[0] );
  CHECK_EQUAL( (EntityHandle)34, conn[1] );
  CHECK_EQUAL( (EntityHandle)32, conn[2] );
  CHECK_EQUAL( (EntityHandle)43, conn[3] );
}

void mb_write_mesh_test()
{
  Core moab;
  Interface* MB = &moab;
  load_file_one( MB );
  ErrorCode result;

  std::string file_name = "mb_write.g";

    // no need to get lists, write out the whole mesh
  result = MB->write_mesh(file_name.c_str());
  CHECK_ERR(result);

    //---------The following tests outputting meshsets that are in meshsets of blocks ---/

    //lets create a block meshset and put some entities and meshsets into it
  EntityHandle block_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, block_ms );
  CHECK_ERR(result);

    //make another meshset to put quads in, so SHELLs can be written out
  EntityHandle block_of_shells;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, block_of_shells); 
  CHECK_ERR(result);

    //tag the meshset so it's a block, with id 100
  int id = 100;
  Tag tag_handle;
  result = MB->tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, tag_handle ) ;
  CHECK_ERR(result);
  result = MB->tag_set_data( tag_handle, &block_ms, 1, &id ) ;
  CHECK_ERR(result);
  id = 101;
  result = MB->tag_set_data( tag_handle, &block_of_shells, 1, &id ) ;
  CHECK_ERR(result);

    // set dimension tag on this to ensure shells get output; reuse id variable
  result = MB->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, tag_handle) ;
  CHECK_ERR(result);
  id = 3;
  result = MB->tag_set_data( tag_handle, &block_of_shells, 1, &id ) ;
  CHECK_ERR(result);

    //get some entities (tets) 
  Range temp_range;
  result = MB->get_entities_by_type(0,  MBHEX, temp_range ) ;
  CHECK_ERR(result);

  Range::iterator iter, end_iter;
  iter = temp_range.begin();
  end_iter = temp_range.end();

    //add evens to 'block_ms'
  std::vector<EntityHandle> temp_vec; 
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 == 0 ) 
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( block_ms, &temp_vec[0], temp_vec.size()); 
  CHECK_ERR(result);


    //make another meshset
  EntityHandle ms_of_block_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ms_of_block_ms);
  CHECK_ERR(result);

    //add some entities to it
  temp_vec.clear();
  iter = temp_range.begin();
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 )  //add all odds
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( ms_of_block_ms, &temp_vec[0], temp_vec.size() ); 
  CHECK_ERR(result);

    //add the other meshset to the block's meshset
  result = MB->add_entities( block_ms, &ms_of_block_ms, 1);
  CHECK_ERR(result);


    //---------------testing sidesets----------------/

    //lets create a sideset meshset and put some entities and meshsets into it
  EntityHandle sideset_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, sideset_ms );
  CHECK_ERR(result);

    //tag the meshset so it's a sideset, with id 104
  id = 104;
  result = MB->tag_get_handle( NEUMANN_SET_TAG_NAME, 1, MB_TYPE_INTEGER, tag_handle ) ;
  CHECK_ERR(result);

  result = MB->tag_set_data( tag_handle, &sideset_ms, 1, &id ) ;
  CHECK_ERR(result);

    //get some entities (tris) 
  temp_range.clear();
  result = MB->get_entities_by_type(0,  MBQUAD, temp_range ) ;
  CHECK_ERR(result);

  iter = temp_range.begin();
  end_iter = temp_range.end();

    //add evens to 'sideset_ms'
  temp_vec.clear(); 
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 == 0 ) 
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( sideset_ms, &temp_vec[0], temp_vec.size() ); 
  CHECK_ERR(result);

    //make another meshset
  EntityHandle ms_of_sideset_ms;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ms_of_sideset_ms);
  CHECK_ERR(result);

    //add some entities to it
  temp_vec.clear();
  iter = temp_range.begin();
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 )  //add all odds
      temp_vec.push_back( *iter );
  }
  result = MB->add_entities( ms_of_sideset_ms, &temp_vec[0], temp_vec.size() ); 
  CHECK_ERR(result);

    //add the other meshset to the sideset's meshset
  result = MB->add_entities( sideset_ms, &ms_of_sideset_ms, 1);
  CHECK_ERR(result);

    //---------test sense on meshsets (reverse/foward)-------//

    //get all quads whose x-coord = 2.5 and put them into a meshset_a 
  EntityHandle meshset_a;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_a );
  CHECK_ERR(result);

  temp_range.clear();
  result = MB->get_entities_by_type(0,  MBQUAD, temp_range ) ;
  CHECK_ERR(result);

  std::vector<EntityHandle> nodes, entity_vec;
  std::copy(temp_range.begin(), temp_range.end(), std::back_inserter(entity_vec));
  result = MB->get_connectivity(&entity_vec[0], entity_vec.size(), nodes);
  CHECK_ERR(result);
  assert( nodes.size() == 4 * temp_range.size() );
  temp_vec.clear(); 
  std::vector<double> coords(3*nodes.size());
  result = MB->get_coords(&nodes[0], nodes.size(), &coords[0]);
  CHECK_ERR(result);
  
  unsigned int k = 0;
  for(Range::iterator it = temp_range.begin(); it != temp_range.end(); it++) {
    if( coords[12*k] == 2.5 && coords[12*k+3] == 2.5 &&
        coords[12*k+6] == 2.5 && coords[12*k+9] == 2.5 )
      temp_vec.push_back(*it);
    k++;
  }
  result = MB->add_entities( meshset_a, &temp_vec[0], temp_vec.size() );
  CHECK_ERR(result);
  result = MB->add_entities( block_of_shells, &temp_vec[0], temp_vec.size());
  CHECK_ERR(result);

    //put these quads into a different meshset_b and tag them with a reverse sense tag
  EntityHandle meshset_b;
  result = MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_b );
  CHECK_ERR(result);

  result = MB->add_entities( meshset_b, &meshset_a, 1);
  CHECK_ERR(result);


  result = MB->tag_get_handle( "SENSE", 1, MB_TYPE_INTEGER, tag_handle, MB_TAG_SPARSE|MB_TAG_CREAT );
  CHECK_ERR(result);

  int reverse_value = -1;
  result = MB->tag_set_data( tag_handle, &meshset_b, 1, &reverse_value ) ; 
  CHECK_ERR(result);


    //get some random quad, whose x-coord != 2.5, and put it into a different meshset_c
    //and tag it with a reverse sense tag

  iter = temp_range.begin();
  end_iter = temp_range.end();

  temp_vec.clear();
  for(; iter != end_iter; iter++ )
  {
    std::vector<EntityHandle> nodes;
    result = MB->get_connectivity( &(*iter), 1, nodes );
    CHECK_ERR(result);

    bool not_equal_2_5 = true; 
    for(unsigned int k=0; k<nodes.size(); k++ )
    {
      double coords[3] = {0};

      result = MB->get_coords( &(nodes[k]), 1, coords );
      CHECK_ERR(result);

      if( coords[0] == 2.5 )
      {
        not_equal_2_5 = false;
        break;
      }
    }

    if( not_equal_2_5 && nodes.size()> 0)
    {
      temp_vec.push_back( *iter );
      break;
    }
  }

  EntityHandle meshset_c;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_c );
    
  
  result = MB->tag_get_handle( "SENSE", 1, MB_TYPE_INTEGER, tag_handle ); 
  CHECK_ERR(result);

  reverse_value = -1;
  result = MB->tag_set_data( tag_handle, &meshset_c, 1, &reverse_value ) ; 
  CHECK_ERR(result);

  MB->add_entities( meshset_c, &temp_vec[0], temp_vec.size() );
  MB->add_entities( block_of_shells, &temp_vec[0], temp_vec.size());


    //create another meshset_abc, adding meshset_a, meshset_b, meshset_c 
  EntityHandle meshset_abc;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, meshset_abc );

  temp_vec.clear();
  temp_vec.push_back( meshset_a );
  temp_vec.push_back( meshset_b );
  temp_vec.push_back( meshset_c );

  MB->add_entities( meshset_abc, &temp_vec[0], temp_vec.size());


    //tag it so it's a sideset
  id = 444;
  result = MB->tag_get_handle( "NEUMANN_SET", 1, MB_TYPE_INTEGER, tag_handle ) ;
  CHECK_ERR(result);

  result = MB->tag_set_data( tag_handle, &meshset_abc, 1, &id ) ;
  CHECK_ERR(result);



    //---------------do nodesets now -----------------//


    //lets create a nodeset meshset and put some entities and meshsets into it
  EntityHandle nodeset_ms;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, nodeset_ms );

    //tag the meshset so it's a nodeset, with id 119
  id = 119;
  result = MB->tag_get_handle( DIRICHLET_SET_TAG_NAME, 1, MB_TYPE_INTEGER, tag_handle ) ;
  CHECK_ERR(result);

  result = MB->tag_set_data( tag_handle, &nodeset_ms, 1, &id ) ;
  CHECK_ERR(result);

    //get all Quads 
  temp_range.clear();
  result = MB->get_entities_by_type(0,  MBQUAD, temp_range ) ;
  CHECK_ERR(result);


    //get all the nodes of the tris
  Range nodes_of_quads;
  iter = temp_range.begin();
  end_iter = temp_range.end();


  for(; iter != end_iter; iter++ )
  {
    std::vector<EntityHandle> nodes;
    result = MB->get_connectivity( &(*iter), 1, nodes);
    CHECK_ERR(result);

    for(unsigned int k=0; k<nodes.size(); k++ )
      nodes_of_quads.insert( nodes[k] ); 

  }

  iter = nodes_of_quads.begin();
  end_iter = nodes_of_quads.end();

    //add evens to 'nodeset_ms'
  temp_vec.clear(); 
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 == 0 ) 
      temp_vec.push_back( *iter );
  }
  MB->add_entities( nodeset_ms, &temp_vec[0], temp_vec.size() ); 


    //make another meshset
  EntityHandle ms_of_nodeset_ms;
  MB->create_meshset(MESHSET_ORDERED | MESHSET_TRACK_OWNER, ms_of_nodeset_ms);

    //add some entities to it
  temp_vec.clear();
  iter = nodes_of_quads.begin();
  end_iter = nodes_of_quads.end();
  for(; iter != end_iter; iter++)
  {
    if( ID_FROM_HANDLE( *iter ) % 2 )  //add all odds
      temp_vec.push_back( *iter );
  }
  MB->add_entities( ms_of_nodeset_ms, &temp_vec[0], temp_vec.size() ); 

    //add the other meshset to the nodeset's meshset
  MB->add_entities( nodeset_ms, &ms_of_nodeset_ms, 1);


    // no need to get lists, write out the whole mesh
  file_name = "mb_write2.g";
  std::vector<EntityHandle> output_list;
  output_list.push_back( block_ms );
  output_list.push_back( sideset_ms );
  output_list.push_back( meshset_abc );
  output_list.push_back( nodeset_ms );
  output_list.push_back( block_of_shells );
  result = MB->write_mesh(file_name.c_str(), &output_list[0], output_list.size());
  CHECK_ERR(result);
}


struct TestType {
  EntityType moab_type;
  ExoIIElementType exo_type;
  int num_nodes;
  std::string name;
};

void check_type( const TestType& type )
{
  int has_mid_nodes[4];
  CN::HasMidNodes( type.moab_type, type.num_nodes, has_mid_nodes );
  
  CHECK_EQUAL( type.moab_type, ExoIIUtil::ExoIIElementMBEntity[type.exo_type] );
  CHECK_EQUAL( type.name, std::string( ExoIIUtil::ElementTypeNames[type.exo_type] ) );
  CHECK_EQUAL( type.num_nodes, ExoIIUtil::VerticesPerElement[type.exo_type] );
  switch (CN::Dimension(type.moab_type)) {
    case 3: CHECK_EQUAL( has_mid_nodes[3], ExoIIUtil::HasMidNodes[type.exo_type][3] );
    case 2: CHECK_EQUAL( has_mid_nodes[2], ExoIIUtil::HasMidNodes[type.exo_type][2] );
    case 1: CHECK_EQUAL( has_mid_nodes[1], ExoIIUtil::HasMidNodes[type.exo_type][1] );
  }
  
  Core moab;
  ExoIIUtil tool(&moab);
  CHECK_EQUAL( type.exo_type, tool.element_name_to_type( type.name.c_str() ) );
  CHECK_EQUAL( type.name, std::string(tool.element_type_name( type.exo_type ) ) );
}
  
void test_types()
{
  const TestType types[] = {
    { MBVERTEX,  EXOII_SPHERE,     1, "SPHERE" },
    { MBEDGE,    EXOII_SPRING,     1, "SPRING" },
    { MBEDGE,    EXOII_BAR,        2, "BAR" },
    { MBEDGE,    EXOII_BAR2,       2, "BAR2" },
    { MBEDGE,    EXOII_BAR3,       3, "BAR3" },
    { MBEDGE,    EXOII_BEAM,       2, "BEAM" },
    { MBEDGE,    EXOII_BEAM2,      2, "BEAM2" },
    { MBEDGE,    EXOII_BEAM3,      3, "BEAM3" },
    { MBEDGE,    EXOII_TRUSS,      2, "TRUSS" },
    { MBEDGE,    EXOII_TRUSS2,     2, "TRUSS2" },
    { MBEDGE,    EXOII_TRUSS3,     3, "TRUSS3" },
    { MBTRI,     EXOII_TRI,        3, "TRI" },
    { MBTRI,     EXOII_TRI3,       3, "TRI3" },
    { MBTRI,     EXOII_TRI6,       6, "TRI6" },
    { MBTRI,     EXOII_TRI7,       7, "TRI7" },
    { MBQUAD,    EXOII_QUAD,       4, "QUAD" },
    { MBQUAD,    EXOII_QUAD4,      4, "QUAD4" },
    { MBQUAD,    EXOII_QUAD5,      5, "QUAD5" },
    { MBQUAD,    EXOII_QUAD8,      8, "QUAD8" },
    { MBQUAD,    EXOII_QUAD9,      9, "QUAD9" },
    { MBQUAD,    EXOII_SHELL,      4, "SHELL" },
    { MBQUAD,    EXOII_SHELL4,     4, "SHELL4" },
    { MBQUAD,    EXOII_SHELL5,     5, "SHELL5" },
    { MBQUAD,    EXOII_SHELL8,     8, "SHELL8" },
    { MBQUAD,    EXOII_SHELL9,     9, "SHELL9" },
    { MBTET,     EXOII_TETRA,      4, "TETRA" },
    { MBTET,     EXOII_TETRA4,     4, "TETRA4" },
    { MBTET,     EXOII_TETRA8,     8, "TETRA8" },
    { MBTET,     EXOII_TETRA10,   10, "TETRA10" },
    { MBTET,     EXOII_TETRA14,   14, "TETRA14" },
    { MBPYRAMID, EXOII_PYRAMID,    5, "PYRAMID" },
    { MBPYRAMID, EXOII_PYRAMID5,   5, "PYRAMID5" },
    { MBPYRAMID, EXOII_PYRAMID10, 10, "PYRAMID10" },
    { MBPYRAMID, EXOII_PYRAMID13, 13, "PYRAMID13" },
    { MBPYRAMID, EXOII_PYRAMID18, 18, "PYRAMID18" },
    { MBPRISM,   EXOII_WEDGE,      6, "WEDGE" },
    { MBKNIFE,   EXOII_KNIFE,      7, "KNIFE" },
    { MBHEX,     EXOII_HEX,        8, "HEX" },
    { MBHEX,     EXOII_HEX8,       8, "HEX8" },
    { MBHEX,     EXOII_HEX9,       9, "HEX9" },
    { MBHEX,     EXOII_HEX20,     20, "HEX20" },
    { MBHEX,     EXOII_HEX27,     27, "HEX27" },
    { MBHEX,     EXOII_HEXSHELL,  12, "HEXSHELL" } };
  const int num_types = sizeof(types)/sizeof(types[0]);
  for (int i = 0; i < num_types; ++i) 
    check_type( types[i] );
}

void read_file( Interface& moab, 
                const char* input_file )
{
  ErrorCode rval = moab.load_file( input_file );
  CHECK_ERR(rval);
}

void write_and_read( Interface& write_mb,
                     Interface& read_mb,
                     EntityHandle block = 0 )
{
  const char* tmp_file = "exodus_test_tmp.g";
  ErrorCode rval;
  
  EntityHandle* write_set_list = &block;
  int write_set_list_len = 0;//(block != 0);
  rval = write_mb.write_file( tmp_file, "EXODUS", 0, 
                            write_set_list, write_set_list_len );
  if (MB_SUCCESS != rval) 
    remove(tmp_file);
  CHECK_ERR(rval);
  
  rval = read_mb.load_file( tmp_file );
  remove( tmp_file );
  CHECK_ERR(rval);
}

void check_ho_elements( Interface& moab, 
                        EntityHandle block,
                        EntityType type,
                        int mid_nodes[4] )
{
  ErrorCode rval;
  Range elems;
  rval = moab.get_entities_by_handle( block, elems );
  CHECK_ERR(rval);
  CHECK(!elems.empty());
  CHECK(elems.all_of_type(type));
  for (Range::const_iterator i = elems.begin(); i != elems.end(); ++i)
    check_ho_element( moab, *i, mid_nodes );
}

// Check that element has expected higher-order nodes
// and that each higher-order node is at the center
// of the sub-entity it is on.
void check_ho_element( Interface& moab, 
                       EntityHandle entity,
                       int mid_nodes[4] )
{
    // get element info
  const EntityType type = TYPE_FROM_HANDLE(entity);
  const EntityHandle* conn;
  int conn_len;
  ErrorCode rval = moab.get_connectivity( entity, conn, conn_len );
  CHECK_ERR(rval);
  std::vector<double> coords(3*conn_len);
  rval = moab.get_coords( conn, conn_len, &coords[0] );
  CHECK_ERR(rval);
  
    // calculate and verify expected number of mid nodes
  int num_nodes = CN::VerticesPerEntity(type);
  for (int d = 1; d <= CN::Dimension(type); ++d)
    if (mid_nodes[d])
      num_nodes += CN::NumSubEntities(type, d);
  CHECK_EQUAL( num_nodes, conn_len );
  
    // verify that each higher-order node is at the center
    // of its respective sub-entity.
  for (int i = CN::VerticesPerEntity(type); i < num_nodes; ++i) {
      // get sub-entity owning ho-node  
    int sub_dim, sub_num;
    CN::HONodeParent( type, num_nodes, i, sub_dim, sub_num );
      // get corner vertex indices
    int sub_conn[8], num_sub;
    if (sub_dim < CN::Dimension(type)) {
      CN::SubEntityVertexIndices( type, sub_dim, sub_num, sub_conn );
      EntityType sub_type = CN::SubEntityType( type, sub_dim, sub_num );
      num_sub = CN::VerticesPerEntity( sub_type );
    }
    else {
      num_sub = CN::VerticesPerEntity(type);
      for (int j = 0; j < num_sub; ++j)
        sub_conn[j] = j;
    }
      // calculate mean of corner vertices
    double mean[3] = {0,0,0};
    for (int j = 0; j < num_sub; ++j) {
      int co = 3*sub_conn[j];
      mean[0] += coords[co  ];
      mean[1] += coords[co+1];
      mean[2] += coords[co+2];
    }
    mean[0] /= num_sub;
    mean[1] /= num_sub;
    mean[2] /= num_sub;
      // verify that higher-order node is at expected location
    CHECK_REAL_EQUAL( mean[0], coords[3*i  ], 1e-6 );
    CHECK_REAL_EQUAL( mean[1], coords[3*i+1], 1e-6 );
    CHECK_REAL_EQUAL( mean[2], coords[3*i+2], 1e-6 );
  }
}


EntityHandle find_block( Interface& mb, EntityType type, const int has_mid_nodes[4] )
{

  ErrorCode rval;
  Tag ho_tag, block_tag;
  rval = mb.tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, block_tag );
  CHECK_ERR(rval);
  rval = mb.tag_get_handle( HAS_MID_NODES_TAG_NAME, 4, MB_TYPE_INTEGER, ho_tag );
  CHECK_ERR(rval);
  
  // get material sets with expected higher-order nodes
  Range blocks;
  Tag tags[2] = {ho_tag, block_tag};
  const void* vals[2] = {has_mid_nodes, NULL};
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, blocks );
  CHECK_ERR(rval);
  
  for (Range::iterator i = blocks.begin(); i != blocks.end(); ++i) {
    int n;
    rval = mb.get_number_entities_by_type( *i, type, n );
    CHECK_ERR(rval);
    if (n > 0)
      return *i;
  }
  
  CHECK(false); // no block matching element type description
  return 0;
}

EntityHandle find_sideset( Interface& mb, 
                             int sideset_id,
                             EntityType side_type )
{
  ErrorCode rval;
  Tag ss_tag;
  rval = mb.tag_get_handle( NEUMANN_SET_TAG_NAME, 1, MB_TYPE_INTEGER, ss_tag );
  CHECK_ERR(rval);
  
  const void* tag_vals[] = { &sideset_id };
  Range side_sets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &ss_tag, tag_vals, 1, side_sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)side_sets.size() );
  return side_sets.front();
}

// Validate elements of specified type.
// Looks for a block containing the specified entity type
// and with the specified mid-node flags set in its
// HAS_MID_NODES_TAG.
void test_ho_elements( EntityType type, int num_nodes )
{
  Core mb_impl1, mb_impl2;
  Interface &mb1 = mb_impl1, &mb2 = mb_impl2;
  int ho_flags[4];
  CN::HasMidNodes( type, num_nodes, ho_flags );

    // read file 
  read_file( mb1, ho_file );
    // test element connectivity order
  EntityHandle block = find_block( mb1, type, ho_flags );
  CHECK(block != 0);
  check_ho_elements( mb1, block, type, ho_flags );
  
    // write block and read it back in
  write_and_read( mb1, mb2, block );
    // test element connectivity order on re-read data
  block = find_block( mb2, type, ho_flags );
  CHECK(block != 0);
  check_ho_elements( mb2, block, type, ho_flags );
}

void test_read_side( int id,
                     EntityType sideset_type,
                     int sideset_nodes_per_elem,
                     bool shell_side )
{
  // read test file
  Core mb_impl;
  Interface& moab = mb_impl;
  read_file( moab, ho_file );
  
  // get side set 
  EntityHandle set = find_sideset( moab, id, sideset_type );
  CHECK(set != 0);
  
  // check expected element connectivity
  int ho_flags[4];
  CN::HasMidNodes( sideset_type, sideset_nodes_per_elem, ho_flags );
  check_ho_elements( moab, set, sideset_type, ho_flags );
  
  if (shell_side)
    return;
  
  // check that each element is on the boundary of the mesh
  Range elems;
  ErrorCode rval = mb_impl.get_entities_by_handle( set, elems );
  CHECK_ERR(rval);
  
  int dim = CN::Dimension( sideset_type );
  for (Range::iterator i= elems.begin(); i != elems.end(); ++i) {
    Range adj;
    rval = mb_impl.get_adjacencies( &*i, 1, dim+1, false, adj, Interface::UNION );
    CHECK_ERR(rval);
    CHECK_EQUAL( 1, (int)adj.size() );
  }
    
}

void test_read_ids_common( const char* file_name,
                           const char* tag_name,
                           const int* expected_vals,
                           int num_expected )
{
  Core mb;
  ReadNCDF reader(&mb);
  
  FileOptions opts("");
  std::vector<int> values;
  ErrorCode rval = reader.read_tag_values( file_name, tag_name, opts, values );
  CHECK_ERR(rval);
  
  std::vector<int> expected( expected_vals, expected_vals+num_expected );
  std::sort( values.begin(), values.end() );
  std::sort( expected.begin(), expected.end() );
  CHECK_EQUAL( expected, values );
}

void test_read_block_ids() {
  const int expected[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
  test_read_ids_common( ho_file, MATERIAL_SET_TAG_NAME, expected, sizeof(expected)/sizeof(expected[0]) );
}

void test_read_sideset_ids() {
  const int expected[] = { 1, 2, 3, 4 };
  test_read_ids_common( ho_file, NEUMANN_SET_TAG_NAME, expected, sizeof(expected)/sizeof(expected[0]) );
}

void test_read_nodeset_ids() {
  test_read_ids_common( ho_file, DIRICHLET_SET_TAG_NAME, 0, 0 );
}

void test_read_alternate_coord_format()
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval = mb.load_file( alt_file );
  CHECK_ERR(rval);
  
  const double exp[] = { 0, 0, 0,
                         1, 0, 0,
                         1, 1, 0,
                         0, 1, 0,
                         0, 0, 1,
                         1, 0, 1,
                         1, 1, 1,
                         0, 1, 1 };
  
  Range hexes;
  rval = mb.get_entities_by_type( 0, MBHEX, hexes );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)1, hexes.size() );
  EntityHandle hex = hexes.front();
  const EntityHandle* conn;
  int len;
  rval = mb.get_connectivity( hex, conn, len );
  CHECK_ERR(rval);
  CHECK_EQUAL(8, len);
  double act[3*8];
  rval = mb.get_coords( conn, len, act );
  CHECK_ERR(rval);
  CHECK_ARRAYS_EQUAL( exp, 3*8, act, 3*8 );
}
