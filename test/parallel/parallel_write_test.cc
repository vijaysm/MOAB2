#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <sstream>

using namespace moab;

#define TPRINT(A) tprint( (A) )
static void tprint(const char* A) 
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  char buffer[128]; 
  sprintf(buffer,"%02d: %6.2f: %s\n", rank, (double)clock()/CLOCKS_PER_SEC, A);
  fputs( buffer, stderr ); 
}

const int DEFAULT_INTERVALS  = 2;
const char* DEFAULT_FILE_NAME = "parallel_write_test.h5m";

// Create mesh for each processor that is a cube of hexes
// with the specified interval count along each edge.  Cubes
// of mesh will be positioned and vertex IDs assigned such that
// there are shared entities.  If the cubic root of the number
// of processors is a whole number, then each processors mesh
// will be a cube in a grid with that number of processor blocks
// along each edge.  Otherwise processor blocks will be arranged
// within the subset of the grid that is the ceiling of the cubic
// root of the comm size such that there are no disjoint regions.
ErrorCode generate_mesh( Interface& moab, int intervals );

const char args[] = "[-i <intervals>] [-o <filename>] [-L <filename>] [-g <n>]";
void help() {
  std::cout << "parallel_write_test " << args << std::endl
            << "  -i <N>    Each processor owns an NxNxN cube of hex elements (default: " << DEFAULT_INTERVALS << ")" << std::endl
            << "  -o <name> Retain output file and name it as specified." << std::endl
            << "  -L <name> Write local mesh to file name prefixed with MPI rank" << std::endl
            << "  -g <n>    Specify writer debug output level" << std::endl
            << "  -R        Skip resolve of shared entities (interface ents will be duplicated in file)" << std::endl
            << std::endl
            << "This program creates a (non-strict) subset of a regular hex mesh "
               "such that the mesh is already partitioned, and then attempts to "
               "write that mesh using MOAB's parallel HDF5 writer.  The mesh size "
               "will scale with the number of processors and the number of elements "
               "per processor (the latter is a function of the value specified "
               "with the '-i' flag.)" << std::endl 
            << std::endl
            << "Let N = ceil(cbrt(P)), where P is the number of processes.  "
               "The mesh will be some subset of a cube with one corner at the "
               "origin and the other at (N,N,N).  Each processor will own a "
               "non-overlapping 1x1x1 unit block of mesh within that cube.  "
               "If P is a power of 3, then the entire NxNxN cube will be "
               "filled with hex elements.  Otherwise, some connected subset "
               "of the cube will be meshed.  Each processor is assigned a "
               "sub-block of the cube by rank where the blocks are enumerated "
               "sequentally with x increasing most rapidly and z least rapidly." << std::endl
            << std::endl
            << "The size of the mesh owned by each processor is controlled by "
               "the number of intervals along each edge of its block of mesh.  "
               "If each block has N intervals, than each processor will have "
               "N^3 hex elements." << std::endl
            << std::endl;
}
 
int main( int argc, char* argv[] )
{
  int ierr = MPI_Init( &argc, &argv );
  if (ierr) {
    std::cerr << "MPI_Init failed with error code: " << ierr << std::endl;
    return ierr;
  }
  int rank;
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (ierr) {
    std::cerr << "MPI_Comm_rank failed with error code: " << ierr << std::endl;
    return ierr;
  }
  int size;
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &size );
  if (ierr) {
    std::cerr << "MPI_Comm_size failed with error code: " << ierr << std::endl;
    return ierr;
  }

    // settings controlled by CL flags
  const char* output_file_name = 0;
  const char* indiv_file_name = 0;
  int intervals = 0, debug_level = 0;
    // state for CL flag processing
  bool expect_intervals = false;
  bool expect_file_name = false;
  bool expect_indiv_file = false;
  bool skip_resolve_shared = false;
  bool expect_debug_level = false;
    // process CL args
  for (int i = 1; i < argc; ++i) {
    if (expect_intervals) {
      char* endptr = 0;
      intervals = (int)strtol( argv[i], &endptr, 0 );
      if (*endptr || intervals < 1) {
        std::cerr << "Invalid block interval value: " << argv[i] << std::endl;
        return 1;
      }
      expect_intervals = false;
    }
    else if (expect_indiv_file) {
      indiv_file_name = argv[i];
      expect_indiv_file = false;
    }
    else if (expect_file_name) {
      output_file_name = argv[i];
      expect_file_name = false;
    }
    else if (expect_debug_level) {
      debug_level = atoi(argv[i]);
      if (debug_level < 1) {
        std::cerr << "Invalid argument following -g flag: \"" << argv[i] << '"' << std::endl;
        return 1;
      }
      expect_debug_level = false;
    }
    else if (!strcmp( "-i", argv[i]))
      expect_intervals = true;
    else if (!strcmp( "-o", argv[i]))
      expect_file_name = true;
    else if (!strcmp( "-L", argv[i]))
      expect_indiv_file = true;
    else if (!strcmp( "-R", argv[i]))
      skip_resolve_shared = true;
    else if (!strcmp( "-g", argv[i]))
      expect_debug_level = true;
    else if (!strcmp( "-h", argv[i])) {
      help();
      return 0;
    }
    else {
      std::cerr << "Unexpected argument: " << argv[i] << std::endl
                << "Usage: " << argv[0] << " " << args << std::endl
                << "       " << argv[0] << " -h" << std::endl
                << " Try '-h' for help." << std::endl;
      return 1;
    }
  }
    // Check for missing argument after last CL flag
  if (expect_file_name || expect_intervals || expect_indiv_file) {
    std::cerr << "Missing argument for '" << argv[argc-1] << "'" << std::endl;
    return 1;
  }
    // If intervals weren't specified, use default
  if (intervals == 0) {
    std::cout << "Using default interval count: " << DEFAULT_INTERVALS << std::endl;
    intervals = DEFAULT_INTERVALS;
  }
    // If no output file was specified, use default one and note that
    // we need to delete it when the test completes.
  bool keep_output_file = true;
  if (!output_file_name) {
    output_file_name = DEFAULT_FILE_NAME;
    keep_output_file = false;
  }
  
    // Create mesh
TPRINT("Generating mesh");
  double gen_time = MPI_Wtime();
  Core mb;
  Interface& moab = mb;
  ErrorCode rval = generate_mesh( moab, intervals );
  if (MB_SUCCESS != rval) {
    std::cerr << "Mesh creation failed with error code: " << rval << std::endl;
    return (int)rval;
  }
  gen_time = MPI_Wtime() - gen_time;

    // Write out local mesh on each processor if requested.
  if (indiv_file_name) {
TPRINT("Writing individual file");
    char buffer[64];
    int width = (int)ceil( log10( size ) );
    sprintf(buffer,"%0*d-", width, rank );
    std::string name(buffer);
    name += indiv_file_name;
    rval = moab.write_file( name.c_str() );
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to write file: " << name << std::endl;
      return (int)rval;
    }
  }

  double res_time = MPI_Wtime();
  Range hexes;
  moab.get_entities_by_type( 0, MBHEX, hexes );
  if (!skip_resolve_shared) {
TPRINT("Resolving shared entities");
      // Negotiate shared entities using vertex global IDs
    ParallelComm* pcomm = new ParallelComm( &moab, MPI_COMM_WORLD );
    rval = pcomm->resolve_shared_ents( 0, hexes, 3, 0 );
    if (MB_SUCCESS != rval) {
      std::cerr << "ParallelComm::resolve_shared_ents failed" << std::endl;
      return rval;
    }
  }
  res_time = MPI_Wtime() - res_time;
  
TPRINT("Beginning parallel write");
  double write_time = MPI_Wtime();
    // Do parallel write
  clock_t t = clock();
  std::ostringstream opts;
  opts << "PARALLEL=WRITE_PART";
  if (debug_level > 0)
    opts << ";DEBUG_IO=" << debug_level;
  rval = moab.write_file( output_file_name, "MOAB", opts.str().c_str() );
  t = clock() - t;
  if (MB_SUCCESS != rval) {
    std::string msg;
    moab.get_last_error( msg );
    std::cerr << "File creation failed with error code: " << moab.get_error_string( rval ) << std::endl;
    std::cerr << "\t\"" << msg << '"' << std::endl;
    return (int)rval;
  }
  write_time = MPI_Wtime() - write_time;
  
  double times[3] = { gen_time, res_time, write_time };
  double max[3] = { 0, 0, 0 };
  MPI_Reduce( times, max, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  
    // Clean up and summarize
  if (0 == rank) {
    double sec = (double)t / CLOCKS_PER_SEC;
    std::cout << "Wrote " << hexes.size()*size << " hexes in " << sec << " seconds." << std::endl;

    if (!keep_output_file) {
TPRINT("Removing written file");
      remove( output_file_name );
    }
    
    std::cout << "Wall time: generate: " << max[0] 
              << ", resovle shared: " << max[1]
              << ", write_file: " << max[2] << std::endl;
  }
  
TPRINT("Finalizing MPI");
  return MPI_Finalize();
}

#define IDX(i,j,k) ((num_interval+1)*((num_interval+1)*(k) + (j)) + (i))

ErrorCode generate_mesh( Interface& moab, int num_interval )
{
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  ErrorCode rval;
  Tag global_id;
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, global_id );
  if (MB_SUCCESS != rval)
    return rval;
  
    // Each processor will own one cube of mesh within
    // an 3D grid of cubes.  Calculate the dimensions of
    // that grid in numbers of cubes.
  int root = 1;
  while (root*root*root < size)
    ++root;
  int num_x_blocks = root;
  int num_y_blocks = root-1;
  int num_z_blocks = root-1;
  if (num_x_blocks * num_y_blocks * num_z_blocks < size)
    ++num_y_blocks;
  if (num_x_blocks * num_y_blocks * num_z_blocks < size)
    ++num_z_blocks;
  
    // calculate position of this processor in grid
  int my_z_block = rank / (num_x_blocks * num_y_blocks);
  int rem = rank % (num_x_blocks * num_y_blocks);
  int my_y_block = rem / num_x_blocks;
  int my_x_block = rem % num_x_blocks;  
  
    // Each processor's cube of mesh will be num_iterval^3 elements
    // and will be 1.0 units on a side
  
    // create vertices
  const int num_x_vtx = num_interval * num_x_blocks + 1;
  const int num_y_vtx = num_interval * num_y_blocks + 1;
  const int x_offset = my_x_block * num_interval;
  const int y_offset = my_y_block * num_interval;
  const int z_offset = my_z_block * num_interval;
  double step = 1.0 / num_interval;
  std::vector<EntityHandle> vertices( (num_interval+1)*(num_interval+1)*(num_interval+1) );
  std::vector<EntityHandle>::iterator v = vertices.begin();
  for (int k = 0; k <= num_interval; ++k) {
    for (int j = 0; j <= num_interval; ++j) {
      for (int i = 0; i <= num_interval; ++i) {
        double coords[] = { my_x_block + i*step,
                            my_y_block + j*step,
                            my_z_block + k*step };
        EntityHandle h;
        rval = moab.create_vertex( coords, h );
        if (MB_SUCCESS != rval)
          return rval;
        
        int id = 1 +  x_offset + i 
                   + (y_offset + j) * num_x_vtx
                   + (z_offset + k) * num_x_vtx * num_y_vtx;
        rval = moab.tag_set_data( global_id, &h, 1, &id );
        if (MB_SUCCESS != rval)
          return rval;
        
        assert(v != vertices.end());
        *v++ = h;
      }
    }
  }
  
    // create hexes
  for (int k = 0; k < num_interval; ++k) {
    for (int j = 0; j < num_interval; ++j) {
      for (int i = 0; i < num_interval; ++i) {
        assert( IDX(i+1,j+1,k+1) < (int)vertices.size() );
        const EntityHandle conn[] = { vertices[IDX(i,  j,  k  )],
                                        vertices[IDX(i+1,j,  k  )],
                                        vertices[IDX(i+1,j+1,k  )],
                                        vertices[IDX(i,  j+1,k  )],
                                        vertices[IDX(i,  j,  k+1)],
                                        vertices[IDX(i+1,j,  k+1)],
                                        vertices[IDX(i+1,j+1,k+1)],
                                        vertices[IDX(i,  j+1,k+1)] };
        EntityHandle elem;
        rval = moab.create_element( MBHEX, conn, 8, elem );
        if (MB_SUCCESS != rval)
          return rval;
      }
    }
  }
  /*
    // create interface quads 
  for (int j = 0; j < num_interval; ++j) {
    for (int i = 0; i < num_interval; ++i) {
      EntityHandle h;

      const EntityHandle conn1[] = { vertices[IDX(i,  j,  0)],
                                       vertices[IDX(i+1,j,  0)],
                                       vertices[IDX(i+1,j+1,0)],
                                       vertices[IDX(i,  j+1,0)] };
      rval = moab.create_element( MBQUAD, conn1, 4, h );
      if (MB_SUCCESS != rval)
        return rval;

      const EntityHandle conn2[] = { vertices[IDX(i,  j,  num_interval)],
                                       vertices[IDX(i+1,j,  num_interval)],
                                       vertices[IDX(i+1,j+1,num_interval)],
                                       vertices[IDX(i,  j+1,num_interval)] };
      rval = moab.create_element( MBQUAD, conn2, 4, h );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  for (int k = 0; k < num_interval; ++k) {
    for (int i = 0; i < num_interval; ++i) {
      EntityHandle h;

      const EntityHandle conn1[] = { vertices[IDX(i,  0,k  )],
                                       vertices[IDX(i+1,0,k  )],
                                       vertices[IDX(i+1,0,k+1)],
                                       vertices[IDX(i,  0,k+1)] };
      rval = moab.create_element( MBQUAD, conn1, 4, h );
      if (MB_SUCCESS != rval)
        return rval;

      const EntityHandle conn2[] = { vertices[IDX(i,  num_interval,k  )],
                                       vertices[IDX(i+1,num_interval,k  )],
                                       vertices[IDX(i+1,num_interval,k+1)],
                                       vertices[IDX(i,  num_interval,k+1)] };
      rval = moab.create_element( MBQUAD, conn2, 4, h );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  for (int k = 0; k < num_interval; ++k) {
    for (int j = 0; j < num_interval; ++j) {
      EntityHandle h;

      const EntityHandle conn1[] = { vertices[IDX(0,j,  k  )],
                                       vertices[IDX(0,j+1,k  )],
                                       vertices[IDX(0,j+1,k+1)],
                                       vertices[IDX(0,j,  k+1)] };
      rval = moab.create_element( MBQUAD, conn1, 4, h );
      if (MB_SUCCESS != rval)
        return rval;

      const EntityHandle conn2[] = { vertices[IDX(num_interval,j,  k  )],
                                       vertices[IDX(num_interval,j+1,k  )],
                                       vertices[IDX(num_interval,j+1,k+1)],
                                       vertices[IDX(num_interval,j,  k+1)] };
      rval = moab.create_element( MBQUAD, conn2, 4, h );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  */
  return MB_SUCCESS;
}


    
