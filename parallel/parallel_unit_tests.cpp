#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "mpi.h"
#include <iostream>
#include <sstream>
#include <assert.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)


#define CHKERR(a) do { \
  MBErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::cerr << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl;\
    return val; \
  } \
} while (false) 

/**************************************************************************
                     Utility Method Declarations
 **************************************************************************/
 
// Get either all the entities in a set or only the entities of
// a specific type, depending on whether or not type == MBMAXTYPE.
MBErrorCode get_set_contents( MBInterface& moab, 
                              MBEntityHandle set, 
                              MBEntityType type,
                              MBRange& contents_out );

// Get processors an entity is shared with.
MBErrorCode get_sharing_processors( MBInterface& moab,
                                    MBEntityHandle entity,
                                    std::vector<int>& other_procs_out );

// Get mesh entities owned by a geometric entity set.
MBErrorCode get_geom_inclusive_mesh( MBInterface& moab,
                                     MBEntityHandle set,
                                     MBEntityType type,
                                     MBRange& contents_out );
                                     
// Test if is_my_error is non-zero on any processor in MPI_COMM_WORLD
int is_any_proc_error( int is_my_error );

/**************************************************************************
                           Test  Declarations
 **************************************************************************/

MBErrorCode test_elements_on_several_procs( const char* filename );


/**************************************************************************
                              Main Method
 **************************************************************************/

#define RUN_TEST(A, B) run_test( &A, #A, B )

int run_test( MBErrorCode (*func)(const char*), 
              const char* func_name,
              const char* file_name )
{
  MBErrorCode result = (*func)(file_name);
  int is_err = is_any_proc_error( (MB_SUCCESS != result) );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (rank == 0) {
    if (is_err) 
      std::cout << func_name << " FAILED!!" << std::endl;
    else
      std::cout << func_name << " success" << std::endl;
  }
  
  return is_err;
}

int main( int argc, char* argv[] )
{
  MPI_Init(&argc, &argv);


  const char* filename;
  if (argc == 1) {
#ifdef SRCDIR
    filename = STRINGIFY(SRCDIR) "/ptest.cub";
#else
    filename = "ptest.cub";
#endif
  }
  else if (argc == 2) {
    filename = argv[1];
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
    return 1;
  }
  
  int num_errors = 0;
  
  num_errors += RUN_TEST( test_elements_on_several_procs, filename );
  
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (rank == 0) {
    if (!num_errors) 
      std::cout << "All tests passed" << std::endl;
    else
      std::cout << num_errors << " TESTS FAILED!" << std::endl;
  }
  
  MPI_Finalize();
  return num_errors;
}


/**************************************************************************
                     Utility Method Implementations
 **************************************************************************/
 
MBErrorCode get_set_contents( MBInterface& moab,
                              MBEntityHandle set,
                              MBEntityType type,
                              MBRange& contents_out )
{
  if (type == MBMAXTYPE) 
    return moab.get_entities_by_handle( set, contents_out );
  else
    return moab.get_entities_by_type( set, type, contents_out );
}

MBErrorCode get_sharing_processors( MBInterface& moab,
                                    MBEntityHandle entity,
                                    std::vector<int>& other_procs_out )
{
  MBErrorCode rval;
  
    // get tags for parallel data
  MBTag sharedp_tag, sharedps_tag, pstatus_tag;
  const char* ptag_names[] = { PARALLEL_SHARED_PROC_TAG_NAME,
                               PARALLEL_SHARED_PROCS_TAG_NAME,
                               PARALLEL_STATUS_TAG_NAME };
  MBTag* tag_ptrs[] = { &sharedp_tag, &sharedps_tag, &pstatus_tag };
  const int ntags = sizeof(ptag_names)/sizeof(ptag_names[0]);
  for (int i = 0; i < ntags; ++i) {
    rval = moab.tag_get_handle( ptag_names[i], *tag_ptrs[i] ); 
    CHKERR(rval);
  }
  
  other_procs_out.clear();
  char status;
  rval = moab.tag_get_data( pstatus_tag, &entity, 1, &status ); CHKERR(rval);
  if (!(status & PSTATUS_SHARED))
    return MB_SUCCESS;
  
  int proc_id;
  rval = moab.tag_get_data( sharedp_tag, &entity, 1, &proc_id ); CHKERR(rval);
  if (proc_id >= 0) {
    other_procs_out.push_back( proc_id );
    return MB_SUCCESS;
  }
  
  int procs[MAX_SHARING_PROCS];
  rval = moab.tag_get_data( sharedps_tag, &entity, 1, procs ); CHKERR(rval);
  for (int i = 0; i < MAX_SHARING_PROCS && procs[i] >= 0; ++i)
    other_procs_out.push_back(procs[i]);
  return MB_SUCCESS;
}
  
                          

MBErrorCode get_geom_inclusive_mesh( MBInterface& moab,
                                     MBEntityHandle set,
                                     MBEntityType type,
                                     MBRange& contents_out )
{
  MBRange children, child_ents, tmp_range;
  MBErrorCode rval;
  
  rval = get_set_contents( moab, set, type, contents_out ); CHKERR(rval);
  rval = moab.get_child_meshsets( set, children );    CHKERR(rval);
    
  for (MBRange::iterator i = children.begin(); i != children.end(); ++i) {
    child_ents.clear();
    rval = get_set_contents( moab, *i, type, child_ents );  CHKERR(rval);
    tmp_range = contents_out.subtract( child_ents );
    contents_out.swap( tmp_range );
  }
  
  return MB_SUCCESS;
}

int is_any_proc_error( int is_my_error )
{
  int result;
  int err = MPI_Allreduce( &is_my_error, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  return err || result;
}


/**************************************************************************
                           Test  Implementations
 **************************************************************************/

MBErrorCode test_elements_on_several_procs( const char* filename )
{
  MBCore mb_instance;
  MBInterface& moab = mb_instance;
  MBErrorCode rval;
  MBEntityHandle set;
  const char* geom_names[] = { "vertex", "curve", "surface", "volume", "unknown" };
  
  rval = moab.load_file( filename, set, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS" );
  CHKERR(rval);
  
  MBTag geom_tag, id_tag;
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, geom_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag ); CHKERR(rval);  
  
    // search for geometric entity sets that contain a vertex
    // that is shared by more than two 
  MBRange geom_ents, several_proc_ents, invalid_proc_ents;
  rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag, 0, 1, geom_ents );
  CHKERR(rval);
  for (MBRange::iterator i = geom_ents.begin(); i != geom_ents.end(); ++i) {
    MBRange verts;
    rval = get_geom_inclusive_mesh( moab, *i, MBVERTEX, verts ); CHKERR(rval);
    if (verts.empty())
      continue;
    std::vector<int> procs;
    MBRange::iterator j = verts.begin();
    rval = get_sharing_processors( moab, *j, procs ); CHKERR(rval);
    if (procs.size() > 1)
      several_proc_ents.insert( *i );
    for (++j; j != verts.end(); ++j) {
      std::vector<int> tmp_procs;
      rval = get_sharing_processors( moab, *j, tmp_procs ); CHKERR(rval);
      if (tmp_procs != procs) 
        invalid_proc_ents.insert( *i );
    }
  }
  
    // if the vertices owned by any geometric entity do not
    // have consistent shared processor ids, list the geometric
    // entities and return failure.
  int my_error = 0;
  if (!invalid_proc_ents.empty()) {
    my_error = 1;
    std::cerr << "Vertices owned by a single geometric entity are "
              << "not shared by the same set of processors for the "
              << "following geometric entities: ";
    for (MBRange::iterator i = invalid_proc_ents.begin();  
         i != invalid_proc_ents.end(); ++i) {
      int dim;
      int id;
      rval = moab.tag_get_data( geom_tag, &*i, 1, &dim );
      if (MB_SUCCESS != rval)
        dim = 4;
      rval = moab.tag_get_data( id_tag, &*i, 1, &id );
      if (MB_SUCCESS != rval)
        id = -1;
      std::cerr << geom_names[dim] << " " << id << ", ";
    }
    std::cerr << std::endl;
    
    my_error = 1;
  }
  if (is_any_proc_error(my_error))
    return MB_FAILURE;
  
    // now scan the list of geometric entities for
    // any for which the higher-dimension entities
    // don't match the vertices.
  for (MBRange::iterator i = several_proc_ents.begin();
       i != several_proc_ents.end(); ++i) {
    MBRange ents;
    rval = get_geom_inclusive_mesh( moab, *i, MBMAXTYPE, ents ); CHKERR(rval);
    std::vector<int> exp_procs, ent_procs;
    rval = get_sharing_processors( moab, ents.front(), exp_procs ); CHKERR(rval);
    for (MBRange::iterator j = ents.upper_bound(MBVERTEX); j != ents.end(); ++j) {
      rval = get_sharing_processors( moab, *j, ent_procs ); CHKERR(rval);
      if (ent_procs != exp_procs)
        invalid_proc_ents.insert( *i );
    }
  }
  
  
    // if the elements owned by any geometric entity do not
    // have consistent shared processor ids, list the geometric
    // entities and return failure.
  my_error = 0;
  if (!invalid_proc_ents.empty()) {
    my_error = 1;
    std::cerr << "Elements owned by a single geometric entity are "
              << "not shared by the same set of processors for the "
              << "following geometric entities: ";
    for (MBRange::iterator i = invalid_proc_ents.begin();  
         i != invalid_proc_ents.end(); ++i) {
      int dim;
      int id;
      rval = moab.tag_get_data( geom_tag, &*i, 1, &dim );
      if (MB_SUCCESS != rval)
        dim = 4;
      rval = moab.tag_get_data( id_tag, &*i, 1, &id );
      if (MB_SUCCESS != rval)
        id = -1;
      std::cerr << geom_names[dim] << " " << id << ", ";
    }
    std::cerr << std::endl;
    
    my_error = 1;
  }
  if (is_any_proc_error(my_error))
    return MB_FAILURE;
  
  return MB_SUCCESS;
}
