
#include "moab/RefineSlabs.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iostream>
#include <assert.h>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>


namespace moab{


  // ======================================================== AHF functions needing to be filled in =============

  ErrorCode RefineSlabs::new_refinement_ahf( size_t /*num_hexes_memory_estimate*/ )
  {
    refinement_ahf = new HalfFacetRep(mbImpl /*, num_hexes_memory_estimate*/ ); // AHF todo, make use of the memory estimate
    if (!ahf)
      return MB_MEMORY_ALLOCATION_FAILED;
    return MB_SUCCESS;
  }


  ErrorCode RefineSlabs::create_hex( EntityHandle fine_nodes[8], EntityHandle & new_hex )
  {
    ErrorCode error = mbImpl->create_element(MBHEX, fine_nodes, 8, new_hex);    
    assert( error = MB_SUCCESS );
    // refinement_AHF
    // tell AHF about the new hex // AHF todo
    return error;
  } 

  ErrorCode RefineSlabs::create_node( EntityHandle node, EntityHandle &new_node )
  {
    // get coordinates from node
    const double *xp, *yp, *zp;
    ErrorCode error = mbImpl->get_coords( node, xp, yp, zp );
    assert( error = MB_SUCCESS );


    double coord[3];
    coord[0] = *xp; 
    coord[1] = *yp;
    coord[2] = *zp;
    error = mbImpl->create_vertex(coord, new_node);
    assert( error = MB_SUCCESS );


    // refinement_AHF
    // tell AHF about the new node // AHF todo
    return error;
  } 

  // replace the coarse hexes and quads with the fine hexes and quads in the global database mbImpl
  // return the new hexes for the caller. The coarse entities no longer exist.
  ErrorCode RefineSlabs::replace_mesh( Entities &, Entities &, Entities &, Entities& )
  //  ErrorCode RefineSlabs::replace_mesh( Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities &fine_quads )
  {
    // refinement_ahf elements -> moved into ahf elements
    // we don't use ahf after this function call, so maybe we don't need to build the ahf, just moab core?
    // AHF todo
    // return MB_FAILURE;

    return MB_SUCCESS;
  }
  void RefineSlabs::replace_node( EntityHandle, int, EntityHandle )
  //  void RefineSlabs::replace_node( EntityHandle chex, int node_lid, EntityHandle new_node)
  {
    // refinement_ahf
    // just replace the node in the individual hex, don't worry about updating the hex-to-hex connectivity until later when we call update_AHF_connectivity
    // AHF todo
    ;
  }
  void RefineSlabs::udpate_AHF_connectivity()
  {
    // since we've replaced a lot of nodes in hexes, update the hex-to-hex connectivity in refinement_ahf for traversals
    // refinement_ahf; 
    // AHF todo
  }
  
  void RefineSlabs::get_all_hexes( EntityHandle node, Entities &hexes, bool is_coarse )
  {
    HalfFacetRep *use_ahf = (is_coarse) ? ahf : refinement_ahf;
    use_ahf->get_up_adjacencies_vert_3d( node, hexes );
  }
  void RefineSlabs::get_all_quads( EntityHandle node, Entities &quads, bool is_coarse )
  {
    HalfFacetRep *use_ahf = (is_coarse) ? ahf : refinement_ahf;
    use_ahf->get_up_adjacencies_vert_2d( node, quads );
  }

    // given a SlabEdge (edge and vertex) of a hex, find the other two edges sharing that vertex
  void RefineSlabs::get_adj( const SlabEdge &edge, SlabEdge &adj1, SlabEdge &adj2 )
  {
    // todo AHF, move to AHF
    // this should go into AHF and not be hard-coded here

    adj1.hex = edge.hex;
    adj2.hex = edge.hex;
    adj1.head_node = adj2.head_node = edge.head_node;    

    int head_index = get_hex_node_index( edge.hex, edge.head_node );
    int tail_index = get_hex_node_index( edge.hex, edge.tail_node );

    EntityHandle hex_nodes[8];
    get_hex_nodes( edge.hex, hex_nodes);

    // todo AHF, move to AHF
    // 24 different edges and orientations, check for each one
    switch ( head_index )
    {
      case 0:
      switch( tail_index)
      {        
        case 1:
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 3;
        adj2.tail_node = hex_nodes[4];  adj1.edge_lid = 4;
        break;
        case 3:
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 0;
        adj2.tail_node = hex_nodes[4];  adj1.edge_lid = 4;
        break;
        case 4:
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 3;
        adj2.tail_node = hex_nodes[1];  adj1.edge_lid = 0;
        break;
      }
      break;

      case 1:
      switch( tail_index)
      {        
        case 0:
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[2];  adj1.edge_lid = 1;
        break;
        case 5:
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 0;
        adj2.tail_node = hex_nodes[2];  adj1.edge_lid = 1;
        break;
        case 2:
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[0];  adj1.edge_lid = 0;
        break;
      }
      break;

      case 2:
      switch( tail_index)
      {        
        case 1:
        adj1.tail_node = hex_nodes[6];  adj1.edge_lid = 6;
        adj2.tail_node = hex_nodes[3];  adj1.edge_lid = 2;
        break;
        case 3:
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 1;
        adj2.tail_node = hex_nodes[6];  adj1.edge_lid = 6;
        break;
        case 6:
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 1;
        adj2.tail_node = hex_nodes[3];  adj1.edge_lid = 2;
        break;
      }
      break;


      case 3:
      switch( tail_index)
      {        
        case 2:
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 3;
        adj2.tail_node = hex_nodes[7];  adj1.edge_lid = 7;
        break;
        case 0:
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 2;
        adj2.tail_node = hex_nodes[7];  adj1.edge_lid = 7;
        break;
        case 7:
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 2;
        adj2.tail_node = hex_nodes[0];  adj1.edge_lid = 3;
        break;
      }
      break;


      case 4:
      switch( tail_index)
      {        
        case 0:
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 8;
        adj2.tail_node = hex_nodes[7];  adj1.edge_lid = 11;
        break;
        case 5:
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 4;
        adj2.tail_node = hex_nodes[7];  adj1.edge_lid = 11;
        break;
        case 7:
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 4;
        adj2.tail_node = hex_nodes[5];  adj1.edge_lid = 8;
        break;
      }
      break;

      case 5:
      switch( tail_index)
      {        
        case 1:
        adj1.tail_node = hex_nodes[4];  adj1.edge_lid = 8;
        adj2.tail_node = hex_nodes[6];  adj1.edge_lid = 9;
        break;
        case 4:
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[6];  adj1.edge_lid = 9;
        break;
        case 6:
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[4];  adj1.edge_lid = 8;
        break;
      }
      break;



      case 6:
      switch( tail_index)
      {        
        case 2:
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 9;
        adj2.tail_node = hex_nodes[7];  adj1.edge_lid = 10;
        break;
        case 5:
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 6;
        adj2.tail_node = hex_nodes[7];  adj1.edge_lid = 10;
        break;
        case 7:
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 6;
        adj2.tail_node = hex_nodes[5];  adj1.edge_lid = 9;
        break;
      }
      break;



      case 7:
      switch( tail_index)
      {        
        case 3:
        adj1.tail_node = hex_nodes[4];  adj1.edge_lid = 11;
        adj2.tail_node = hex_nodes[6];  adj1.edge_lid = 10;
        break;
        case 4:
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 7;
        adj2.tail_node = hex_nodes[6];  adj1.edge_lid = 10;
        break;
        case 6:
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 7;
        adj2.tail_node = hex_nodes[4];  adj1.edge_lid = 11;
        break;
      }
      break;

      default: ;
    }
  }

  // given a SlabEdge (edge and vertex) of a hex, find the other two edges that are opposite the edge in a quad of the hex
  void RefineSlabs::get_opp( const SlabEdge &edge, SlabEdge &opp1, SlabEdge &opp2 ) 
  {
    // todo AHF, move to AHF
    // this should go into AHF and not be hard-coded here
    opp1.hex = edge.hex;
    opp2.hex = edge.hex;

    int head_index = get_hex_node_index( edge.hex, edge.head_node );
    int tail_index = get_hex_node_index( edge.hex, edge.tail_node );

    EntityHandle hex_nodes[8];
    get_hex_nodes( edge.hex, hex_nodes);

    // todo AHF, move to AHF
    // 24 different edges and orientations, check for each one
    switch ( head_index )
    {
      case 0:
      switch( tail_index)
      {        
        case 1:
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[5];  opp1.edge_lid = 8;
        break;
        case 3:
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[7];  opp1.edge_lid = 11;
        break;
        case 4:
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[5];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[3];  opp2.tail_node = hex_nodes[7];  opp1.edge_lid = 7;
        break;
      }
      break;

      case 1:
      switch( tail_index)
      {        
        case 0:
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[4];  opp1.edge_lid = 8;
        break;
        case 5:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[4];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[2];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 6;
        break;
        case 2:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 9;
        break;
      }
      break;

      case 2:
      switch( tail_index)
      {        
        case 1:
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[5];  opp1.edge_lid = 9;
        break;
        case 3:
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[7];  opp1.edge_lid = 10;
        break;
        case 6:
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[5];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[3];  opp2.tail_node = hex_nodes[7];  opp1.edge_lid = 7;
        break;
      }
      break;


      case 3:
      switch( tail_index)
      {        
        case 2:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 10;
        break;
        case 0:
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[4];  opp1.edge_lid = 11;
        break;
        case 7:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[4];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[2];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 6;
        break;
      }
      break;


      case 4:
      switch( tail_index)
      {        
        case 0:
        opp1.head_node = hex_nodes[5];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[3];  opp1.edge_lid = 7;
        break;
        case 5:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 10;
        break;
        case 7:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 9;
        break;
      }
      break;

      case 5:
      switch( tail_index)
      {        
        case 1:
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[4];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[2];  opp2.tail_node = hex_nodes[6];  opp1.edge_lid = 6;
        break;
        case 4:
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[7];  opp1.edge_lid = 10;
        break;
        case 6:
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[7];  opp1.edge_lid = 11;
        break;
      }
      break;



      case 6:
      switch( tail_index)
      {        
        case 2:
        opp1.head_node = hex_nodes[5];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[3];  opp1.edge_lid = 7; 
        break;
        case 5:
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[4];  opp1.edge_lid = 11;
        break;
        case 7:
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[4];  opp1.edge_lid = 8;
        break;
      }
      break;



      case 7:
      switch( tail_index)
      {        
        case 3:
        opp1.head_node = hex_nodes[4];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[2];  opp1.edge_lid = 6;
        break;
        case 4:
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[5];  opp1.edge_lid = 9;
        break;
        case 6:
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[5];  opp1.edge_lid = 8;
        break;
      }
      break;

      default: ;
    }
  }


  // ========================================================

  RefineSlabs::RefineSlabs(Core *impl)
  {
    assert(NULL != impl);
    mbImpl = impl;

    ErrorCode error;
    error = initialize();
    if (error != MB_SUCCESS)
    {
      std::cout<<"Error initializing RefineSlabs\n"<<std::endl;
      exit(1);
    }
  }

  RefineSlabs::~RefineSlabs()
  {
    delete ahf;
  }

  ErrorCode RefineSlabs::initialize()
  {
#ifdef USE_AHF
  std::cout << "macro USE_AHF is defined" << std::endl;
#else
  std::cout << "macro USE_AHF is undefined" << std::endl;
#endif

    ErrorCode error;
    ahf = new HalfFacetRep(mbImpl);
    if (!ahf)
      return MB_MEMORY_ALLOCATION_FAILED;

    //Check for mixed entity type
    bool chk_mixed = ahf->check_mixed_entity_type();
    if (chk_mixed)
      MB_SET_ERR(  MB_NOT_IMPLEMENTED, "Encountered a mesh with mixed entity types");

    error = ahf->initialize(); MB_CHK_ERR(error);
    Range _inverts, _inedges, _infaces, _incells;
    error = ahf->get_entity_ranges(_inverts, _inedges, _infaces, _incells);  MB_CHK_ERR(error);

    // Check for mixed dimensional mesh
    if ((!_inedges.empty() && !_infaces.empty()) ||(!_inedges.empty() &&  !_incells.empty()) || (!_infaces.empty() && !_incells.empty()))
      MB_SET_ERR(MB_NOT_IMPLEMENTED, "Encountered a mixed-dimensional mesh");

    // Only hexes are supported
    //Check for supported entity type
    size_t meshdim;
    if (!_inedges.empty())
      {
        MB_SET_ERR(MB_FAILURE, "Encountered a 1d mesh, but only hexes are supported");
        meshdim = 1;
      }
    else if (!_infaces.empty())
      {
        MB_SET_ERR(MB_FAILURE, "Encountered a 2D mesh, but only hexes are supported");
        meshdim = 2;
      }
    else if (!_incells.empty())
      {
        EntityType type = mbImpl->type_from_handle(_incells[0]);
        if(type != MBHEX)
          MB_SET_ERR(MB_FAILURE, "Encountered a 3D mesh that wasn't hexes, but only hexes are supported");
        meshdim = 3;
      }

    return MB_SUCCESS;
  }

  /************************************************************
   *     Interface Functions                                  *
   ************************************************************/

  ErrorCode RefineSlabs::refine_mesh(Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities &fine_quads)
  {
    
    // find boundary
    // define which vertices are on the boundary of the coarse_hexes refinement set?
    ErrorCode err = MB_SUCCESS;
    err = mark_hex_nodes(coarse_hexes);
    if (err == MB_SUCCESS)
      err = mark_surface_nodes(coarse_quads);

    // find slabs    
    // ideally two parallel sheets of hexes
    EntitiesVec slabs;
    if (err == MB_SUCCESS)
      err = find_slabs( coarse_hexes, coarse_quads,  slabs );

    // initialize refinement
    if (err == MB_SUCCESS)
      err = initialize_refinement( coarse_hexes, coarse_quads );

    // pillow slabs
    // refine the hexes of each slab, referencing each hex's current refinement
    if (err == MB_SUCCESS)
      err = pillow_slabs( slabs );

    // replace mesh
    if (err == MB_SUCCESS)
      err = replace_mesh( coarse_hexes, coarse_quads, fine_hexes, fine_quads );

    delete refinement_ahf;


    return err;
  }


  ErrorCode RefineSlabs::initialize_refinement( Entities &coarse_hexes, Entities &/*coarse_quads*/ )
  {

    // estimate memory size for AHF
    // this could be made tighter by considering the quads, and the actual slabs found.
    size_t memory_estimate = 8 * coarse_hexes.size() + 8 * coarse_hexes.size();
    //                        refine each hex into 8      transition == outer pillow

    // make a new AHF to store the refinement
    ErrorCode error = new_refinement_ahf( memory_estimate ); 
    assert( error = MB_SUCCESS );

    // copy nodes
    // each node and its copy point to each other
    copy_nodes( coarse_hexes );

    // initial refined hexes = original hexes
    // copy hexes
    // a hex points to its refinement (other way not needed?)
    copy_hexes( coarse_hexes );

    // ? anything needed for the quads at this point?
    // should we copy them?
    return error;

  }

  ErrorCode RefineSlabs::copy_nodes( Entities &coarse_hexes )
  {
    for (size_t i = 0; i < coarse_hexes.size(); ++i )
    {
      EntityHandle chex = coarse_hexes[i];
      EntityHandle coarse_nodes [8];
      // get 8 coarse nodes of hex
      get_hex_nodes( chex, coarse_nodes);

      for (size_t j = 0; j < 8; ++j)
      {
        if ( !get_hex_refinement( coarse_nodes[j] ))
        {
          // create a new node in the ahf with the same coordinates as the coarse one
          // associate them with each other
          copy_node( coarse_nodes[j] );          
        }
      }
    }
    return MB_SUCCESS;
  }

  void RefineSlabs::get_fine_nodes( EntityHandle *coarse_nodes, EntityHandle *fine_nodes, int num_nodes)
  {
    for (int i = 0; i < num_nodes; ++i)
    {
      get_fine_node( coarse_nodes[i], fine_nodes[i] );
      assert( get_fine_node( coarse_nodes[i], fine_nodes[i] ) );
    }
  }


  ErrorCode RefineSlabs::copy_hexes( Entities &coarse_hexes )
  {
    for (size_t i = 0; i < coarse_hexes.size(); ++i )
    {
      EntityHandle chex = coarse_hexes[i];
      EntityHandle coarse_nodes [8], fine_nodes[8]; 
      // todo next level of sophistication
      // generalize to 2nd order hexes with more nodes, etc
      get_hex_nodes( chex, coarse_nodes);
      get_fine_nodes( coarse_nodes, fine_nodes, 8);
      EntityHandle fhex;
      create_hex( fine_nodes, fhex );
      add_refined_hex( chex, fhex );
    }
    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::find_slabs( Entities &coarse_hexes, Entities &/*coarse_quads*/, EntitiesVec &slabs )
  {
    for ( size_t c = 0; c < coarse_hexes.size(); ++c )
    {
      EntityHandle hex = coarse_hexes[c];
      SlabEdge slab_edge;
      for (int edge_lid = 0; edge_lid < 12; ++edge_lid)
      {
        // find a seed edge, then
        // recursively extend the slab
        // find_seed_edge may increment edge_lid until a good one is found, or 12 is reached
        if ( find_seed_edge( hex, edge_lid, slab_edge ) )
        {
          Entities slab;
          reset_in_slab();
          add_edge( slab_edge, slab );
          extend_slab( slab_edge, slab);

          // uniquify the slab, removing redundant hexes
          std::sort( slab.begin(), slab.end() );
          slab.erase( std::unique( slab.begin(), slab.end() ), slab.end() );

          slabs.push_back( slab );
        }
      }
    }
    return MB_SUCCESS;
  }

  void RefineSlabs::extend_slab( SlabEdge slab_edge, Entities &slab )
  {
    // enqueue the passed in slab_edge
    // while the queue is not empty
    //   pop edge
    //   if the edge is good, then
    //      mark it as being in the set 
    //      enqueue the parallel edges 

    // enqueue passed-in edge
    std::vector< SlabEdge > edge_queue;
    edge_queue.push_back(slab_edge);

    std::vector< SlabEdge > ortho_slab_edges, upper_slab_edges, parallel_slab_edges;

    // while the queue is not empty
    // note: edge_queue.size() grows inside the loop
    for( size_t i = 0; i < edge_queue.size(); ++i )
    {
      slab_edge = edge_queue[i];

      if ( is_good_slab_edge(slab_edge, ortho_slab_edges, upper_slab_edges, parallel_slab_edges) )
      {
        add_edge( slab_edge, slab );

        // enqueue the parallel edges
        // grow in bfs order, not dfs, in order to get better geometrically shaped slabs
        edge_queue.insert( edge_queue.end(), parallel_slab_edges.begin(), parallel_slab_edges.end() );
      }
    }
  }

  void RefineSlabs::add_edge( SlabEdge &slab_edge, Entities &slab )
  {
    // find the distinguished vertex
    // get all the (coarse) hexes containing the vertex
    // add them to the slab

    EntityHandle star_node = slab_edge.head_node;
    std::vector<EntityHandle> star_hexes;
    get_all_hexes( star_node, star_hexes );
    slab.insert( slab.end(), star_hexes.begin(), star_hexes.end() );

    set_in_slab( slab_edge, true);
  }


  bool RefineSlabs::find_seed_edge( EntityHandle hex, int &edge_lid, SlabEdge &slab_edge )
  {
    // caller sets initial value of edge_lid
    for ( /*caller initializes edge_lid*/; edge_lid < 12; ++edge_lid)
    {
      bool good_edge = 
        is_good_slab_edge( hex, edge_lid, 0, slab_edge ) ||
        is_good_slab_edge( hex, edge_lid, 1, slab_edge );
      if( good_edge )
        return true;
    }
    return false;
  }

  bool RefineSlabs::is_good_slab_edge( const SlabEdge &slab_edge )
  {
    std::vector< SlabEdge > ortho_slab_edges, upper_slab_edges, parallel_slab_edges;
    return is_good_slab_edge( slab_edge, ortho_slab_edges, upper_slab_edges, parallel_slab_edges );
  }
  
  bool RefineSlabs::is_good_slab_edge( const SlabEdge &slab_edge, std::vector<SlabEdge> &ortho, std::vector<SlabEdge> &upper, std::vector<SlabEdge> &parallel )
  {
    ortho.clear(); 
    upper.clear();  
    parallel.clear(); 

    // skip if edge has already been refined, or in the current set
    if (get_edge_refinement_level(slab_edge.hex, slab_edge.edge_lid) > 0 || get_in_slab( slab_edge) )
      return false;

    // keep going if a vertex is interior to the volume
    if ( get_geometry_dimension( slab_edge.head_node ) < 3)
      return false;

    get_adjacent_slab_edges( slab_edge, ortho, upper, parallel );
    // ensure no non-ortho edge is already refined
    if ( !none_refined( upper ) )
    {
      return false;
    }

    // ensure upper and ortho edges are not already in the current set!
    if ( any_in_slab(ortho) || any_in_slab(upper) )
      return false;

    return true;
  }

  bool RefineSlabs::is_good_slab_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge )
  {
    // make SlabEdge from one of the two vertices
    get_edge( hex, edge_lid, node_01, slab_edge );

    return is_good_slab_edge( slab_edge );
  }

  bool RefineSlabs::none_refined( std::vector<SlabEdge> &slab_edges )
  {
    for (size_t i = 0; i < slab_edges.size(); ++i )
    {
      SlabEdge &slab_edge = slab_edges[i];
      if (get_edge_refinement_level(slab_edge.hex, slab_edge.edge_lid) > 0)
        return false;
    }
    return true;
  }

  void RefineSlabs::get_adjacent_slab_edges( const SlabEdge &slab_edge, std::vector< SlabEdge > &ortho_slab_edges, 
        std::vector< SlabEdge > &upper_slab_edges, std::vector< SlabEdge > parallel_slab_edges )
  {
    ortho_slab_edges.clear(); 
    upper_slab_edges.clear();  
    parallel_slab_edges.clear(); 

    Entities hexes;
    get_all_hexes( slab_edge.head_node, hexes );
    Entities non_sheet_hexes;
    for (size_t h = 0; h < hexes.size(); ++h )
    {
      EntityHandle hex = hexes[h];      
      SlabEdge match;
      if( get_matching_edge( hex, slab_edge, match ) )
      {
        SlabEdge adj1, adj2, opp1, opp2;
        get_adj( match, adj1, adj2 );
        add_unique( ortho_slab_edges, adj1 );
        add_unique( ortho_slab_edges, adj2 );
        get_opp( match, opp1, opp2 );
        add_unique( parallel_slab_edges, opp1 );
        add_unique( parallel_slab_edges, opp2 );
      }
      else
        non_sheet_hexes.push_back( hex );
    }
    for (size_t h = 0; h < non_sheet_hexes.size(); ++h )
    {
      EntityHandle hex = non_sheet_hexes[h];
      //get_matching_node( slab_edge.head_node, hex, node_lid );
      SlabEdge star[3];
      get_star_edges( hex, slab_edge.head_node, star );
      for ( size_t e = 0; e < 3; ++e )
        if ( unique( ortho_slab_edges, star[e] ) )
          add_unique( upper_slab_edges, star[e] );
    }
  }

  int RefineSlabs::get_hex_node_index( EntityHandle hex, EntityHandle node )
  {
    EntityHandle hex_nodes[8];
    get_hex_nodes( hex, hex_nodes);
    for ( int n = 0; n < 12; ++n )
      if ( hex_nodes[n] == node )
        return n;
    return -1;
  }

  void RefineSlabs::get_star_edges( EntityHandle hex, EntityHandle node, SlabEdge star[3] )
  {
    int i = 0;
    for (size_t e = 0; e < 12; ++e )
    {
      get_edge( hex, e, 0, star[i] );
      if ( star[i].head_node == node )
      {
        // match
        assert( i < 3 );
        ++i;
      }
      else if ( star[i].tail_node == node )
      {
        // match flipped
        assert( i < 3 );
        star[i].flip();
        ++i;
      }
    }
  }


  bool RefineSlabs::get_matching_edge( EntityHandle hex, const SlabEdge &slab_edge, SlabEdge &match )
  {
    for (int e = 0; e < 12; ++e)
    {
      get_edge( hex, e, 0, match );
      if ( slab_edge.nodes_match(match) )
      { 
        if( !slab_edge.directions_match(match) )
        {
          match.flip();
        }
        return true;
      }
    }
    return false;
  }
  bool RefineSlabs::unique( std::vector< SlabEdge > edges, SlabEdge edge )
  {
    for (size_t i = 0; i < edges.size(); ++i )
    {
      if ( edges[i].nodes_match( edge ))
        return false;
    }
   return true;
  }
  bool RefineSlabs::add_unique( std::vector< SlabEdge > edges, SlabEdge edge )
  {
    if (unique( edges, edge ))
    {
      edges.push_back( edge );
      return true;
    }
    return false;
  }
  
  ErrorCode RefineSlabs::pillow_slabs( EntitiesVec &slabs )
  {
    for (size_t i = 0; i < slabs.size(); ++i)
    {
      ErrorCode err_i = pillow_slab( slabs[i] );
      if (err_i != MB_SUCCESS)
        return err_i;
    }
    return MB_SUCCESS;
  }


  void RefineSlabs::shrink_mark_slab( Entities &slab )
  { 
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle hex = slab[i];
      set_shrink_membership( hex, SlabData::INTERNAL );
    }

    // set internal and boundary nodes of the shrink set
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle chex = slab[i];
      EntityHandle hex_nodes[8];
      get_hex_nodes( chex, hex_nodes );
      for ( size_t n = 0; n < 8; ++n )
      {
        EntityHandle hex_node = hex_nodes[n];

        // only check/set it once
        if ( get_shrink_membership(hex_node) == SlabData::EXTERNAL )
        {
          bool is_internal = true;

          // check if its internal to the refinement set
          if (get_membership(hex_node) != SlabData::INTERNAL ) 
            is_internal = false;

          // check if all its hexes are in the set
          if (is_internal)
          {
            Entities node_hexes;
            get_all_hexes( hex_node, node_hexes );
            for ( size_t k = 0; k < node_hexes.size(); ++k )
            {
              EntityHandle node_hex = node_hexes[k];
              if ( get_shrink_membership( node_hex ) == SlabData::EXTERNAL )
              {
                is_internal = false;
                break;
              }
            }
          }
          set_shrink_membership( hex_node, is_internal ? SlabData::INTERNAL : SlabData::BOUNDARY );
        }
      }
    }
  }
  void RefineSlabs::remove_shrink_mark_slab( Entities &slab )
  {
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle chex = slab[i];
      set_shrink_membership( chex, SlabData::EXTERNAL );
      EntityHandle hex_nodes[8];
      get_hex_nodes( chex, hex_nodes );
      for ( size_t n = 0; n < 8; ++n )
      {
        set_shrink_membership( hex_nodes[n], SlabData::EXTERNAL );
      }
    }
  }

  void RefineSlabs::shrink_mark_coarse_slab( Entities &slab )
  { 
    shrink_mark_slab( slab );
  }

  void RefineSlabs::shrink_mark_fine_slab( Entities &slab, Entities &shrink_set)
  {
    // todo next level of sophistication
    // the tricky part is tracking the edges and quads of the fine mesh that are on the coarse mesh quad and hex
    // for now, just put *all* the fine hexes into the set

    // choose the set of hexes

    // for every coarse hex
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle chex = slab[i];
      SlabData *coarse_data = get_slab_data( chex );
      assert( coarse_data );      
      HexRefinement *crefine = coarse_data->refinement;
      assert( crefine );
      Entities &fine_hexes = crefine->fine_hexes; // alias

      // for all its fine hexes
      for (size_t j = 0; j < fine_hexes.size(); ++j)
      {
        EntityHandle fine_hex = fine_hexes[j];

        // todo next level of sophistication
        // additional containment tests here
        // check the coarse_owner of every fine hex, see if that entity is internal or boundary to the set

        shrink_set.push_back( fine_hex );
      }
    }

    // mark the nodes for boundary or internal, as before
    shrink_mark_slab( shrink_set );
  }

  void RefineSlabs::pillow_hexes( Entities &shrink_set, Entities &new_hexes )
  {
    // shrink set and new hexes are both fine-mesh entities

    // split vertices into edges
    // gather boundary entities
    for (size_t i = 0; i < shrink_set.size(); ++i)
    {
      EntityHandle chex = shrink_set[i];
      EntityHandle hex_nodes[8];
      // these are fine hexes and nodes
      get_hex_nodes( chex, hex_nodes ); 
      for ( size_t n = 0; n < 8; ++n )
      {
        EntityHandle node = hex_nodes[n];        
        EntityHandle copy; 
        if ( get_shrink_membership( node ) == SlabData::BOUNDARY && !get_copy( node, copy ) ) 
        {
          copy = shrink_node( node );
          // todo next level of sophistication
          // pick some new geometry for the copied node, say partway towards the hex center
          replace_node( chex, n, copy);
        }
      } 

      // split edges into quads
      // not needed explicitly

      // split quads into hexes
      for (size_t q = 0; q < 6; ++q)
      {
        // get quad nodes, with outward facing normal
        EntityHandle quad_nodes[8];
        get_quad_nodes( chex, hex_nodes, q, quad_nodes );
        bool make_hex = true;
        for ( size_t n = 0; n < 4; n++ )
        {
          EntityHandle copy;
          if ( get_copy( quad_nodes[n], copy ) )
            quad_nodes[n+4] = copy;
          else 
          {
            make_hex = false;
            break;
          }
        }
        // if all nodes had a copy, make a hex out of the copy and the original
        if ( make_hex )
        {
          EntityHandle fhex;
          create_hex( quad_nodes, fhex );
          new_hexes.push_back(fhex);
        }
      }
    }
    // need to update the datastructures for traversal from one hex to another
    // either do that above during replace_node and create_hex, or below
    udpate_AHF_connectivity();
  }
  void RefineSlabs::get_quad_nodes( EntityHandle /*hex*/, const EntityHandle hex_nodes[8], int face_lid, EntityHandle* quad_nodes )
  {
    int face_node_ids[4];
    int num_face_nodes = 4;
    // see CN.hpp and MBEntityType.h
    EntityType quad_type;
    CN::SubEntityNodeIndices( moab::MBHEX, 8, 2, face_lid, quad_type, num_face_nodes, face_node_ids);
    assert( quad_type == moab::MBQUAD );
    assert( num_face_nodes >= 4 );
    for (size_t i = 0; i < 4; ++i)
      quad_nodes[i] = hex_nodes[ face_node_ids[i] ];
  }


  ErrorCode RefineSlabs::get_quad_nodes( EntityHandle quad, EntityHandle quad_nodes[4] )
  {
    int num_nodes = 0;
    const bool corners_only = true;
    const EntityHandle *const_quad_nodes;
    ErrorCode error = mbImpl->get_connectivity(quad, const_quad_nodes, num_nodes, corners_only ); MB_CHK_ERR(error);
    assert(num_nodes == 4);

    // transfer and cast away const
    for (size_t i = 0; i < 4; ++i)
    {
      quad_nodes[i] = *const_cast<EntityHandle*>(&const_quad_nodes[i]);
    }

    return MB_SUCCESS;
  }


  ErrorCode RefineSlabs::get_hex_nodes( EntityHandle hex, EntityHandle hex_nodes[8] )
  {
    int num_nodes = 0;
    const bool corners_only = true;
    const EntityHandle *const_hex_nodes;
    ErrorCode error = mbImpl->get_connectivity(hex, const_hex_nodes, num_nodes, corners_only ); MB_CHK_ERR(error);
    assert( num_nodes == 8 );

    // transfer and cast away const
    for (size_t i = 0; i < 8; ++i)
    {
      hex_nodes[i] = *const_cast<EntityHandle*>(&const_hex_nodes[i]);
    }

    return MB_SUCCESS;
  }

  void RefineSlabs::get_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge )
  {
    // CN.hpp
    EntityType quad_type;
    int edge_nodes[4];
    int num_edge_nodes;
    CN::SubEntityNodeIndices( moab::MBHEX, 8, 2, edge_lid, quad_type, num_edge_nodes, edge_nodes );
    assert( quad_type == moab::MBQUAD );
    assert( num_edge_nodes >= 2 );
    slab_edge.hex = hex;
    slab_edge.edge_lid = edge_lid;
    if (node_01 == 0)
    {
      slab_edge.head_node = edge_nodes[0];
      slab_edge.tail_node = edge_nodes[1];
    }
    else
    {
      slab_edge.head_node = edge_nodes[1];
      slab_edge.tail_node = edge_nodes[0];
    }
  }

  bool RefineSlabs::is_equal( const SlabEntity *m1, const SlabEntity *m2 ) const
  {
    // relies on these type overloading operator == 

    if ( dynamic_cast<const SlabHex*>(m1) && dynamic_cast<const SlabHex*>(m2) )
    {
      return *dynamic_cast<const SlabHex*>(m1) == *dynamic_cast<const SlabHex*>(m2);
    }
    // unused
    // if ( dynamic_cast<const SlabQuad*>(m1) && dynamic_cast<const SlabQuad*>(m2) )
    // {
    //   return *dynamic_cast<const SlabQuad*>(m1) == *dynamic_cast<const SlabQuad*>(m2);
    // }
    if ( dynamic_cast<const SlabEdge*>(m1) && dynamic_cast<const SlabEdge*>(m2) )
    {
      return *dynamic_cast<const SlabEdge*>(m1) == *dynamic_cast<const SlabEdge*>(m2);
    }
    if ( dynamic_cast<const SlabNode*>(m1) && dynamic_cast<const SlabNode*>(m2) )
    {
      return *dynamic_cast<const SlabNode*>(m1) == *dynamic_cast<const SlabNode*>(m2);
    }
    // BadEntity or different classes
    return false;
  }


  ErrorCode RefineSlabs::pillow_slab( Entities &slab )
  {

    // convert coarse set to fine set of hexes to pillow

    // mark each hex as being in the set
    shrink_mark_coarse_slab( slab );


    // find internal and boundary fine hexes
    Entities shrink_set;
    shrink_mark_fine_slab( slab, shrink_set);

    // pillow, splitting vertices, edges, then quads into hexes
    Entities new_hexes;
    pillow_hexes( shrink_set, new_hexes );

    // cleanup 
    remove_shrink_mark_slab( slab );
    remove_shrink_mark_slab( shrink_set );
    remove_shrink_mark_slab( new_hexes );

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::mark_hex_nodes(Entities &coarse_hexes)
  {
    // mark all hexes as being in the set
    for (size_t i = 0; i < coarse_hexes.size(); ++i)
    {
      EntityHandle hex = coarse_hexes[i];
      set_coarse_hex( hex );
    }

    // mark selected nodes as being on the boundary of the refinement hex set
    // a node is internal iff all its hexes are in the set.
    // geometric-surface nodes are marked as on the boundary in mark_surface_nodes
    for (size_t i = 0; i < coarse_hexes.size(); ++i)
    {
      EntityHandle hex = coarse_hexes[i];
      EntityHandle hex_nodes [8];
      get_hex_nodes( hex, hex_nodes );
      for (size_t j = 0; j < 8; ++j)
      {
        EntityHandle node = hex_nodes[j];
        // process each node only once
        if ( !get_coarse(node) )
        {
          bool is_internal(true);

          // check if internal to the geometry
          if ( get_geometry_dimension(node) < 3)
            is_internal = false;

          // check if internal to the set of hexes
          if (is_internal)
          {
            Entities hexes;
            get_all_hexes( node, hexes );
            for ( size_t k = 0; k < hexes.size(); ++k )
            {
              EntityHandle h = hexes[k];
              if ( !get_coarse(h) )
              {
                is_internal = false;
                break;
              }
            }
          }

          set_coarse_node( node );
          set_membership( node, is_internal ? SlabData::INTERNAL : SlabData::BOUNDARY );
        }
      }

    }

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::mark_surface_nodes(Entities &coarse_quads)
  {
    // mark all the quads as being on the surface == a coarse entity
    for (size_t i = 0; i < coarse_quads.size(); ++i)
    {
      EntityHandle quad = coarse_quads[i];
      set_coarse_quad( quad );
    }

    // mark selected nodes as being internal to the set
    // a node is internal iff all of quads that are owned by a geometric surface 
    // are in the coarse_quads set
    for (size_t i = 0; i < coarse_quads.size(); ++i)
    {
      EntityHandle quad = coarse_quads[i];
      EntityHandle quad_nodes [4];
      get_quad_nodes( quad, quad_nodes );
      for (size_t j = 0; j < 4; ++j)
      {
        EntityHandle node = quad_nodes[j];
        SlabData::Membership membership = get_membership( node );
        // its boundary if it was on the geometric surface
        if ( ( membership == SlabData::BOUNDARY ) && ( get_geometry_dimension(node) < 3 ) )
        {
          Entities quads;
          get_all_quads( node, quads );
          bool is_internal(true);
          for ( size_t k = 0; k < quads.size(); ++k )
          {
            EntityHandle q = quads[k];
            if ( (get_geometry_dimension(q)==2) && !get_coarse(q))
            {
              // some of its faces are not in the set, so consider it boundary
              is_internal = false;
              break;
            }
          }
          set_membership( node, is_internal ? SlabData::INTERNAL : SlabData::BOUNDARY );
        }
      }
    }
    return MB_SUCCESS;
  }

  EntityHandle RefineSlabs::copy_node( EntityHandle coarse_node )
  {
    EntityHandle fine_copy;
    create_node( coarse_node, fine_copy );
    set_copy( coarse_node, fine_copy );
    SlabData *slab_data = force_slab_data(coarse_node);    
    SlabData *copy_slab = force_slab_data(fine_copy);
    copy_slab->copy_data( slab_data );

    if ( slab_data->coarsening )
      copy_slab->coarsening = new HexCoarsening( *slab_data->coarsening );
    return fine_copy;
  }

  EntityHandle RefineSlabs::shrink_node( EntityHandle fine_node )
  {
    EntityHandle fine_copy;
    create_node( fine_node, fine_copy );
    set_fine_node( fine_node, fine_copy );
    SlabData *slab_data = force_slab_data(fine_node);
    SlabData *copy_slab = force_slab_data(fine_copy);
    copy_slab->copy_data( slab_data );

    if ( slab_data->coarsening )
      copy_slab->coarsening = new HexCoarsening( *slab_data->coarsening );
    // the coarse mesh should still point to the original fine_node, not its copy
    return fine_copy;
  }

  void RefineSlabs::reset_in_slab()
  {
    edge_data_map.clear();
  }
  bool RefineSlabs::any_in_slab( std::vector< SlabEdge > slab_edges )
  {      
    for (size_t i=0; i<slab_edges.size(); ++i )
    {
      if ( get_in_slab(slab_edges[i]) )
      {
        return true;
      }
    }
    return false;
  }
  bool RefineSlabs::get_in_slab( SlabEdge slab_edge )
  {

    EdgeDataIterator it = edge_data_map.find( std::pair<EntityHandle, EntityHandle>(slab_edge.head_node, slab_edge.tail_node));
    if ( it == edge_data_map.end() )
      return false;
    return it->second;
  }
  void RefineSlabs::set_in_slab( SlabEdge slab_edge, bool new_value )
  {
    edge_data_map[std::pair<EntityHandle, EntityHandle>(slab_edge.head_node, slab_edge.tail_node)] = new_value;
  }


}//namesapce moab

