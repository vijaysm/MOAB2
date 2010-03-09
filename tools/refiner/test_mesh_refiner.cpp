#include "MBCore.hpp"
#include "MBEdgeSizeSimpleImplicit.hpp"
#include "MBSimplexTemplateRefiner.hpp"
#include "MBMeshRefiner.hpp"
#include "MBInterface.hpp"
#include "MBParallelConventions.h"

#ifdef USE_MPI
#include "MBParallelComm.hpp"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "MBmpi.h"
#endif // USE_MPI

#include <iostream>
#include <sstream>
#include <map>

#include "sys/time.h"

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
  //sleep(20);

  // Create the input mesh and, if -new-mesh is specified, an output mesh
  const char* ifname = argc > 1 ? argv[1] : "fourVolsBare.cub";
  bool input_is_output = ( argc > 2 && ! strcmp( argv[2], "-new-mesh" ) ) ? false : true;
  MBInterface* imesh = new MBCore; // ( rank, nprocs );
  MBInterface* omesh = input_is_output ? imesh : new MBCore; // ( rank, nprocs );

#ifdef USE_MPI
  // Use an MBParallelComm object to help set up the input mesh
  MBParallelComm* ipcomm = new MBParallelComm( imesh );
  //ReadParallel* readpar = new ReadParallel( imesh, ipcomm );
#endif // USE_MPI

  MBEntityHandle set_handle;
  std::ostringstream parallel_options;
  parallel_options
    << "PARALLEL=READ_DELETE" << ";" // NB: You can use BCAST_DELETE or READ_DELETE here.
    //<< "PARALLEL=BCAST_DELETE" << ";" // NB: You can use BCAST_DELETE or READ_DELETE here.
    << "PARTITION=MATERIAL_SET" << ";"
    //<< "PARTITION_DISTRIBUTE" << ";"
    << "PARTITION_VAL=" << ( rank + 1 ) << ";"
    << "PARALLEL_RESOLVE_SHARED_ENTS" << ";"
    << "CPUTIME";
#ifdef USE_MPI
  //readpar->load_file( ifname, set_handle, FileOptions( parallel_options.str().c_str() ), 0, 0 );
  imesh->load_file( ifname, set_handle, parallel_options.str().c_str() );
#  if 0
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
#  endif // 0
#else // USE_MPI
  imesh->load_file( ifname, set_handle, parallel_options.str().c_str() );
  imesh->list_entities( 0, 1 );
#endif // USE_MPI
  std::ostringstream ifs;
  ifs << "prerefiner." << nprocs << "." << rank << ".vtk";
  imesh->write_mesh( ifs.str().c_str() );

  // The refiner will need an implicit function to be used as an indicator function for subdivision:
  MBEdgeSizeSimpleImplicit* eval = new MBEdgeSizeSimpleImplicit();
  eval->set_ratio( 2. );
  // Refine the mesh
  MBMeshRefiner* mref = new MBMeshRefiner( imesh, omesh );
  MBSimplexTemplateRefiner* eref = new MBSimplexTemplateRefiner;
  mref->set_entity_refiner( eref );
  //mref->add_vertex_tag( tag_floatular );
  //mref->add_vertex_tag( tag_intular );
  // (We don't add tag_gid to the refiner's tag manager because it is special)
  eref->set_edge_size_evaluator( eval );
  MBRange ents_to_refine;
  imesh->get_entities_by_type( set_handle, MBTET, ents_to_refine ); // refine just the tets
  //ents_to_refine.insert( set_handle ); // refine everything multiple times (because subsets are not disjoint)
  struct timeval tic, toc;
  gettimeofday( &tic, 0 );
  mref->refine( ents_to_refine );
  gettimeofday( &toc, 0 );
  std::cout << "\nTime: " << ( (toc.tv_sec - tic.tv_sec) * 1000 + (toc.tv_usec - tic.tv_usec) / 1000. ) << " ms\n\n";

  std::ostringstream ofs;
  ofs << "refiner." << nprocs << "." << rank << ".vtk";
  omesh->write_mesh( ofs.str().c_str() );
  // Print out the results, one process at a time
#ifdef USE_MPI
#  if 0
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
#  endif // 0
#else // USE_MPI
  omesh->list_entities( 0, 1 );
#endif // USE_MPI

  // Clean up
#ifdef USE_MPI
  //delete readpar;
  delete ipcomm;
#endif // USE_MPI
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
