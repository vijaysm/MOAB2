#include "MBCore.hpp"
#include "MBEdgeSizeSimpleImplicit.hpp"
#include "MBSimplexTemplateRefiner.hpp"
#include "MBMeshRefiner.hpp"
#include "MBInterface.hpp"
#include "MBParallelConventions.h"

#ifdef USE_MPI
#include "MBParallelComm.hpp"
#include "mpi.h"
#endif // USE_MPI

#include <iostream>
#include <map>

int TestMeshRefiner( int argc, char* argv[] )
{
  int nprocs, rank;
#ifdef USE_MPI
  int err = MPI_Init( &argc, &argv );
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else // USE_MPI
  nprocs = 1;
  rank = 0;
#endif // USE_MPI

  // Create the input mesh and, if -new-mesh is specified, an output mesh
  bool input_is_output = ( argc > 1 && ! strcmp( argv[1], "-new-mesh" ) ) ? false : true;
  MBInterface* imesh = new MBCore( rank, nprocs );
  MBInterface* omesh = input_is_output ? imesh : new MBCore( rank, nprocs );

  // Print out the process ranks, one at a time.
#ifdef USE_MPI
  for ( int i = 0; i < nprocs; ++ i )
    {
    MPI_Barrier( MPI_COMM_WORLD );
    if ( i == rank )
      {
      std::cout << "Rank: " << ( rank + 1 ) << " of: " << nprocs << "\n";
      }
    MPI_Barrier( MPI_COMM_WORLD );
    }
#endif // USE_MPI

  // The refiner will need an implicit function to be used as an indicator function for subdivision:
  MBEdgeSizeSimpleImplicit* eval = new MBEdgeSizeSimpleImplicit();
  eval->set_ratio( 2. );
#ifdef USE_MPI
  // Use an MBParallelComm object to help set up the input mesh
  MBParallelComm ipcomm( imesh );
#endif // USE_MPI

  // The arrays below specify the mesh we create.
  // The mesh is distributed across processes.

  // Vertex coordinates
  double coords[][6] = {
    {  0. ,  0.0,  0. ,  0. ,  0.0,  0.  }, // 0
    {  0. ,  0. ,  1. ,  0. ,  0. ,  1.  }, // 1
    { -1. ,  0.0,  0.5, -1. ,  0.0,  0.5 }, // 2
    { -0.5, -1.0,  0.5, -0.5, -1.0,  0.5 }, // 3
    {  0.5, -0.5,  0.5,  0.5, -0.5,  0.5 }, // 4
    {  0.5,  0.5,  0.5,  0.5,  0.5,  0.5 }  // 5
  };

  // Floating-point tag values for testing
  double default_floatular[2] = {
     38.7,   104. };
  double floatular_values[][2] = {
    {  38.7,   104.,    }, // 0
    {   3.141,   2.718, }, // 1
    { 123.456, 789.012, }, // 2
    {   0.,      1.,    }, // 3
    {   1.,      0.,    }, // 4
    {  -1.,      1.,    }  // 5
  };

  // Integer tag values for testing
  int default_intular[4] = {
     7, 11, 24,  7 };
  int intular_values[][4] = {
    {  7, 11, 24,  7, }, // 0
    { 10,  4, 10, 20, }, // 1
    {  1,  2,  3,  4, }, // 2
    { 13, 17, 19, 23, }, // 3
    {  3,  3,  0,  6, }, // 4
    {  5,  4,  3,  2  }  // 5
  };

  // Global IDs of all entities
  int default_gid[] = { -1 };
  int gid_values[] = {
    1, 2, 3, 4, 5, 6, // vertices
    7, 8, 9, 10,      // tetrahedra
    11, 12, 
    13, 14, 
    15, 16, 
    17, 18            // triangles
  };

  // List of vertices resident on each node, as an index into gid_values
  int proc_nodes[4][4] = {
    { 0, 1, 2, 3 },
    { 0, 1, 3, 4 },
    { 0, 1, 4, 5 },
    { 0, 1, 5, 2 },
  };

  // List of nodes with a copy of each vertex (first entry is owner)
  // Ignore this for now.
  int node_procs[6][4] = {
    { 0,  1,  2,  3 }, // 0
    { 0,  1, -1, -1 }, // 1
    { 0,  3, -1, -1 }, // 2
    { 0,  1,  2,  3 }, // 3
    { 1,  2, -1, -1 }, // 4
    { 2,  3, -1, -1 }  // 5
  };

  // List of indices into local node_handles array of two triangles that interface with neighboring processes
  int internal_bdy[4][6] = {
    { 1, 0, 2, 0, 1, 3 },
    { 1, 0, 2, 0, 1, 3 },
    { 1, 0, 2, 0, 1, 3 },
    { 1, 0, 2, 0, 1, 3 }
  };

  MBEntityHandle node_handles[4];
  MBEntityHandle tet_handle;
  MBEntityHandle tri_handles[2];
  MBEntityHandle tri_node_handles[6];

  // Create some tags so we can test the refiner's ability to copy/interpolate them.
  MBTag tag_floatular;
  imesh->tag_create( "floatular", 2 * sizeof( double ), MB_TAG_DENSE, MB_TYPE_DOUBLE, tag_floatular, default_floatular );

  MBTag tag_intular;
  imesh->tag_create( "intular", 4 * sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_intular, default_intular );

  // Get the global ID tag so we can set things up for resolve_shared_ents
  MBTag tag_gid;
  imesh->tag_create( PARALLEL_GID_TAG_NAME, sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_gid, default_gid );

#ifdef USE_MPI
  // Get tags for the various data-distributed mesh annotations
  MBTag tag_sproc;
  MBTag tag_sprocs;
  MBTag tag_shand;
  MBTag tag_shands;
  MBTag tag_pstat;
  ipcomm.get_shared_proc_tags( tag_sproc, tag_sprocs, tag_shand, tag_shands, tag_pstat );
  MBTag tag_part = ipcomm.partition_tag();
#endif // USE_MPI

  // Create vertices for our local tet and prepare tag data for each
  void const* iptrs[4];
  void const* fptrs[4];
  void const* gptrs[4];
  void const* sptrs[4];
  for ( int i = 0; i < 4; ++ i )
    {
    int pnode = proc_nodes[rank][i];
    imesh->create_vertex( coords[pnode], node_handles[i] );
    //std::cout << rank << ": global " << (pnode + 1) << " is handle " << node_handles[i] << "\n";
    iptrs[i] = (void const*) intular_values[pnode];
    fptrs[i] = (void const*) floatular_values[pnode];
    gptrs[i] = (void const*) &gid_values[pnode];
    }

  // Set tag values on vertices
  imesh->tag_set_data( tag_floatular, node_handles, 4, fptrs, 0 );
  imesh->tag_set_data( tag_intular, node_handles, 4, iptrs, 0 );
  imesh->tag_set_data( tag_gid, node_handles, 4, gptrs, 0 );

  // Create a tetrahedron from the vertices
  imesh->create_element( MBTET, node_handles, 4, tet_handle );
  imesh->tag_set_data( tag_gid, &tet_handle, 1, gid_values + 6 + rank );

  // Create two triangles on the interface of the tet with its 2 neighboring processes.
  for ( int i = 0; i < 6; ++ i )
    {
    tri_node_handles[i] = node_handles[internal_bdy[rank][i]];
    }
  imesh->create_element( MBTRI, tri_node_handles, 3, tri_handles[0] );
  imesh->create_element( MBTRI, tri_node_handles + 3, 3, tri_handles[1] );
  imesh->tag_set_data( tag_gid, tri_handles, 2, gid_values + 10 + 2 * rank );
  //imesh->tag_set_data( tag_sprocs, &tet_handle, 1, sptrs, 0 );
  //MBRange proc_ents( node_handles[0], tet_handle );

  // Create a mesh set containing the tet
  MBEntityHandle set_handle;
  imesh->create_meshset( MESHSET_SET, set_handle );
  imesh->add_entities( set_handle, &tet_handle, 1 );

#ifdef USE_MPI
  // Tag the mesh set as a parallel partition
  imesh->tag_set_data( tag_part, &set_handle, 1, &rank );

  // Resolve shared entities (we really only care about vertices but
  // ipcomm.resolve_shared_ents( 0, 0 ) doesn't mark anything up.
  ipcomm.resolve_shared_ents( 3, 0 );
  //ipcomm.resolve_shared_ents( 0, 0 );
  //ipcomm.resolve_shared_ents( 3, 3 );

  // Print out what we have so far, one process at a time
  for ( int i = 0; i < nprocs; ++ i )
    {
    MPI_Barrier( MPI_COMM_WORLD );
    if ( i == rank )
      {
      std::cout << "\n************** Rank: " << ( rank + 1 ) << " of: " << nprocs << "\n";
      imesh->list_entities( 0, 1 );
      std::cout << "**************\n\n";
      }
    MPI_Barrier( MPI_COMM_WORLD );
    }
#else // USE_MPI
  imesh->list_entities( 0, 1 );
#endif // USE_MPI

  // Refine the mesh
  MBMeshRefiner* mref = new MBMeshRefiner( imesh, omesh );
  MBSimplexTemplateRefiner* eref = new MBSimplexTemplateRefiner;
  mref->set_entity_refiner( eref );
  mref->add_vertex_tag( tag_floatular );
  mref->add_vertex_tag( tag_intular );
  // (We don't add tag_gid to the refiner's tag manager because it is special)
  eref->set_edge_size_evaluator( eval );
  MBRange ents_to_refine;
  ents_to_refine.insert( set_handle );
  mref->refine( ents_to_refine );

  // Print out the results, one process at a time
#ifdef USE_MPI
  for ( int i = 0; i < nprocs; ++ i )
    {
    MPI_Barrier( MPI_COMM_WORLD );
    if ( i == rank )
      {
      std::cout << "\n************** Rank: " << ( rank + 1 ) << " of: " << nprocs << "\n";
      omesh->list_entities( 0, 1 );
      std::cout << "**************\n\n";
      }
    MPI_Barrier( MPI_COMM_WORLD );
    }
#else // USE_MPI
  omesh->list_entities( 0, 1 );
#endif // USE_MPI

  // Clean up
  if ( omesh != imesh )
    delete omesh;
  delete imesh;
  delete mref; // mref will delete eref

#ifdef USE_MPI
  err = MPI_Barrier( MPI_COMM_WORLD );
  err = MPI_Finalize();
#endif // USE_MPI

  return 0;
}

int main( int argc, char* argv[] )
{
  return TestMeshRefiner( argc, argv );
}
