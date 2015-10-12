#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <assert.h>
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <unistd.h>
#endif


using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)


#define CHKERR(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::cerr << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl;\
    return val; \
  } \
} while (false) 

#define PCHECK(A) if (is_any_proc_error(!(A))) return report_error(__FILE__,__LINE__)

ErrorCode report_error( const char* file, int line )
{
  std::cerr << "Failure at " << file << ':' << line << std::endl;
  return MB_FAILURE;
}

ErrorCode test_read(const char *filename, const char *option);

#define RUN_TEST(A, B, C) run_test( &A, #A, B, C)

int is_any_proc_error( int is_my_error )
{
  int result = 0;
  int err = MPI_Allreduce( &is_my_error, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  return err || result;
}

int run_test( ErrorCode (*func)(const char*, const char*), 
              const char* func_name,
              const char* file_name,
              const char *option)
{
  ErrorCode result = (*func)(file_name, option);
  int is_err = is_any_proc_error( (MB_SUCCESS != result) );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (rank == 0) {
    if (is_err) 
      std::cout << func_name << " : FAILED!!" << std::endl;
    else
      std::cout << func_name << " : success" << std::endl;
  }
  
  return is_err;
}

int main( int argc, char* argv[] )
{
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  int num_errors = 0;

  const char* filename, *option;
  const char* filename2=NULL;
  if (1 < argc) 
    filename = argv[1];
  else {
#ifdef MESHDIR
    filename = STRINGIFY(MESHDIR) "/64bricks_512hex.h5m";
    filename2 = STRINGIFY(MESHDIR) "/hex_2048.vtk";
#else
    filename = "/64bricks_512hex.h5m";
#endif
  }
#ifdef MOAB_HAVE_HDF5
    //=========== read_delete, geom_dimension, resolve_shared
  option = "PARALLEL=READ_DELETE;PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== read_delete, material_set, resolve_shared
  option = "PARALLEL=READ_DELETE;PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== bcast_delete, geom_dimension, resolve_shared
  option = "PARALLEL=BCAST_DELETE;PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== bcast_delete, material_set, resolve_shared
  option = "PARALLEL=BCAST_DELETE;PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== read_delete, geom_dimension, resolve_shared, exch ghost
  option = "PARALLEL=READ_DELETE;PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== read_delete, material_set, resolve_shared, exch ghost
  option = "PARALLEL=READ_DELETE;PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== bcast_delete, geom_dimension, resolve_shared, exch ghost
  option = "PARALLEL=BCAST_DELETE;PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;";
  num_errors += RUN_TEST( test_read, filename, option );

    //=========== bcast_delete, material_set, resolve_shared, exch ghost
  option = "PARALLEL=BCAST_DELETE;PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;";
  num_errors += RUN_TEST( test_read, filename, option );
#endif
  if (filename2)
  {
    //=========== bcast_delete, trivial, resolve_shared
    option = "PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;";
    num_errors += RUN_TEST( test_read, filename2, option );
    //=========== bcast_delete, trivial, resolve_shared + ghosting
    option = "PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL;PARTITION_DISTRIBUTE;PARALLEL_RESOLVE_SHARED_ENTS;PARALLEL_GHOSTS=3.0.1;";
    num_errors += RUN_TEST( test_read, filename2, option );
  }
  MPI_Finalize();

  return num_errors;
}

ErrorCode test_read(const char *filename, const char *option) 
{
  Core mb_instance;
  Interface& moab = mb_instance;
  ErrorCode rval;

  rval = moab.load_file( filename, 0, option);
  CHKERR(rval);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab, 0);

  rval = pcomm->check_all_shared_handles();
  CHKERR(rval);

  return MB_SUCCESS;
}
