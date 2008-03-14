/**\file varlen_ll.cpp
 * \brief Test HDF5 parallel I/O of variable-length tag data.
 * \author Jason Kraftcheck <kraftche@cae.wisc.edu>
 */

#include "MBCore.hpp"
#include "TestUtil.hpp"
#include "MBParallelConventions.h"

#include <mpi.h>

bool keep_files = false; // controllable with -k flag
bool wait_on_start = false; // start all procs and wait for input on root node

void test_var_length_parallel();

int main(int argc, char* argv[])
{
  MPI_Init( &argc, &argv );

  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-k"))
      keep_files = true;
    else if (!strcmp(argv[i], "-w"))
      wait_on_start = true;
    else {
      fprintf( stderr, "Usage: %s [-k] [-w]\n", argv[0] );
      abort();
    }
  }
  
  if (wait_on_start) {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if (0 == rank) {
      puts( "Press <enter> to begin: ");
      fflush( stdout );
      getchar();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  int err_count = 0;
  err_count += RUN_TEST( test_var_length_parallel );
  return err_count;
}


void test_var_length_parallel()
{
  MBRange::const_iterator i;
  MBErrorCode rval;
  MBCore moab;
  MBInterface &mb = moab;
  MBRange verts;
  MBTag vartag;
  const char* filename = "var-len-para.h5m";
  const char* tagname = "ParVar";
  
  // If this tag doesn't exist, writer will fail
  MBTag junk_tag;
  mb.tag_create( PARALLEL_GID_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, junk_tag, 0 );

  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  
  // Create N+1 vertices on each processor, where N is the rank 
  std::vector<double> coords( 3*rank+3, (double)rank );
  rval = mb.create_vertices( &coords[0], rank+1, verts );
  CHECK_ERR(rval);
  
  // Create a var-len tag
  rval = mb.tag_create_variable_length( tagname, MB_TAG_DENSE, MB_TYPE_INTEGER, vartag );
  CHECK_ERR(rval);
  
  // Write data on each vertex:
  // { n, rank, rank+1, ..., rank+n-1 } where n >= 1
  std::vector<int> data;
  rval = MB_SUCCESS;
  for (i = verts.begin(); i != verts.end(); ++i) {
    MBEntityHandle h = *i;
    const int n = h % 7 + 1;
    data.resize( n+1 );
    data[0] = n;
    for (int j = 0; j < n; ++j)
      data[j+1] = rank + j;
    const int s = (n + 1) * sizeof(int);
    const void* ptrarr[] = { &data[0] };
    MBErrorCode tmperr = mb.tag_set_data( vartag, &h, 1, ptrarr, &s );
    if (MB_SUCCESS != tmperr)
      rval = tmperr;
  }
  CHECK_ERR(rval);
  
  // Write file
  rval = mb.write_file( filename, "MOAB", "PARALLEL=FORMAT" );
  CHECK_ERR(rval);
  
  // Read file.  We only reset and re-read the file on the
  // root processsor.  All other processors keep the pre-write
  // mesh.  This allows many of the tests to be run on all
  // processors.  Running the tests on the pre-write mesh on
  // non-root processors allows us to verify that any problems
  // are due to the file API rather than some other bug.
  MBErrorCode rval2 = rval = MB_SUCCESS;
  if (!rank) {
    moab.~MBCore();
    new (&moab) MBCore;
    rval = mb.load_mesh( filename );
    if (!keep_files) 
      remove( filename );
    rval2 = mb.tag_get_handle( tagname, vartag );
  }
  CHECK_ERR(rval);
  CHECK_ERR(rval2);
  
  // Check that tag is correct
  int tag_size;
  rval = mb.tag_get_size( vartag, tag_size );
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, rval );
  MBTagType storage;
  rval = mb.tag_get_type( vartag, storage );
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  MBDataType type;
  rval = mb.tag_get_data_type( vartag, type);
  CHECK_EQUAL( MB_TYPE_INTEGER, type );
  
  // get vertices
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  // Check consistency of tag data on each vertex 
  // and count the number of vertices for each rank.
  std::vector<int> vtx_counts( numproc, 0 );
  for (i = verts.begin(); i != verts.end(); ++i) {
    MBEntityHandle h = *i;
    int size = -1;
    const void* ptrarr[1] = { 0 };
    rval = mb.tag_get_data( vartag, &h, 1, ptrarr, &size );
    size /= sizeof(int);
    CHECK_ERR( rval );
    const int* data = reinterpret_cast<const int*>(ptrarr[0]);
    CHECK( size >= 2 );
    CHECK( NULL != data );
    CHECK_EQUAL( size-1, data[0] );
    CHECK( data[1] >= 0 && data[1] < numproc );
    ++vtx_counts[data[1]];
    for (int j = 1; j < size-1; ++j)
      CHECK_EQUAL( data[1]+j, data[1+j] );
  }
  
  // Check number of vertices for each rank
  for (int j = 0; j < numproc; ++j) {
    // Only root should have data for other processors.
    if (rank == 0 || rank == j) 
      CHECK_EQUAL( j+1, vtx_counts[j] );
    else 
      CHECK_EQUAL( 0, vtx_counts[j] );
  }
}
