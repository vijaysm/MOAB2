#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <sstream>
#include <assert.h>
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <unistd.h>
#endif


#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

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

/**************************************************************************
                     Utility Method Declarations
 **************************************************************************/

// Get processors an entity is shared with.
ErrorCode get_sharing_processors( Interface& moab,
                                    EntityHandle entity,
                                    std::vector<int>& other_procs_out );

// Create a parallel mesh
//
// Each processor will create four quads.
// Groups of four quads will be arranged as follows:
// +------+------+------+------+------+-----
// |             |             |
// |             |             |
// +    Proc 0   +    Proc 2   +    Proc 4
// |             |             |
// |             |             |
// +------+------+------+------+------+-----
// |             |             |
// |             |             |
// +    Proc 1   +    Proc 3   +    Proc 5
// |             |             |
// |             |             |
// +------+------+------+------+------+-----
//
// Vertices will be enumerated as follows:
// 1------6-----11-----16-----21-----26-----
// |             |             |
// |             |             |
// 2      7     12     17     22     27
// |             |             |
// |             |             |
// 3------8-----13-----18-----23-----28-----
// |             |             |
// |             |             |
// 4      9     14     19     24     29
// |             |             |
// |             |             |
// 5-----10-----15-----20-----25-----30-----
//
// Element IDs will be [4*rank+1,4*rank+5]
//
// If output_sets is non-null, it will be populated with 2 or 3 
// set handles.  A set handle is created for each group of four
// adjacent processes, such that one is created across procs 0 to 4,
// one is created for procs 2 to 5, etc.  An additional set is created
// that spans all of the processes.  The set spanning all procs is 
// returned first.  The other two sets are returned in the order of the
// largest contained process rank, where the last entry is zero if 
// there is only one additional set.
ErrorCode parallel_create_mesh( Interface& mb,
                                int output_vertx_ids[9],
                                EntityHandle output_vertex_handles[9],
                                Range& output_elements,
                                EntityHandle output_sets[3] = 0 );

// Test if is_my_error is non-zero on any processor in MPI_COMM_WORLD
int is_any_proc_error( int is_my_error );

/**************************************************************************
                           Test  Declarations
 **************************************************************************/

// Check consistancy of sharing data.  (E.g. compare sharing procs for
// vertices to that of adjacent elements, compare sharing data for 
// interfaces with that of contained entities, etc.)
ErrorCode test_elements_on_several_procs( const char* filename );
// Test correct ghosting of elements
ErrorCode test_ghost_elements_3_2_1( const char* filename );
ErrorCode test_ghost_elements_3_2_2( const char* filename );
ErrorCode test_ghost_elements_3_0_1( const char* filename );
ErrorCode test_ghost_elements_2_0_1( const char* filename );
// Test exchange of tag data on ghost elements
ErrorCode test_ghost_tag_exchange( const char* filename );
// Bug where exchange_tags fails if dense tag cannot be queried
// for all ghost entities (e.g. no default value)
ErrorCode regression_ghost_tag_exchange_no_default( const char* filename );
// Test owners for interface entities
ErrorCode test_interface_owners( const char* );
// Test data for shared interface entitites with one level of ghosting
ErrorCode regression_owners_with_ghosting( const char* );
// Verify all sharing data for vertices with one level of ghosting
ErrorCode test_ghosted_entity_shared_data( const char* );
// Test assignment of global IDs
ErrorCode test_assign_global_ids( const char* );
// Test shared sets
ErrorCode test_shared_sets( const char* );
// Test reduce_tags
ErrorCode test_reduce_tags( const char* );

/**************************************************************************
                              Main Method
 **************************************************************************/

#define RUN_TEST(A, B) run_test( &A, #A, B)

int run_test( ErrorCode (*func)(const char*), 
              const char* func_name,
              const char* file_name)
{
  ErrorCode result = (*func)(file_name);
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

  int pause_proc = -1;
  const char* filename = 0;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i],"-p")) {
      ++i;
      assert(i < argc);
      pause_proc = atoi( argv[i] );
    }
    else if (!filename) {
      filename = argv[i];
    }
    else {
      std::cerr << "Invalid arg: \"" << argv[i] << '"' << std::endl
                << "Usage: " << argv[0] << " [-p <rank>] [<filename>]" << std::endl;
      exit(1);
    }
  }

  if (!filename) {
#ifdef SRCDIR
    filename = STRINGIFY(MESHDIR) "/64bricks_512hex.h5m";
#else
    filename = "64bricks_512hex.h5m";
#endif
  }

  if (pause_proc != -1) {
#if !defined(_MSC_VER) && !defined(__MINGW32__)
    std::cout << "Processor " << rank << " of " << size << " with PID " << getpid() << std::endl;
    std::cout.flush();
#endif
      // loop forever on requested processor, giving the user time
      // to attach a debugger.  Once the debugger in attached, user
      // can change 'pause'.  E.g. on gdb do "set var pause = 0"
    if (pause_proc == rank) {
      volatile int pause = 1;
      while (pause);
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    std::cout << "Processor " << rank << " resuming" << std::endl;
  }

  
  int num_errors = 0;
  
  num_errors += RUN_TEST( test_elements_on_several_procs, filename );
  num_errors += RUN_TEST( test_ghost_elements_3_2_1, filename );
  num_errors += RUN_TEST( test_ghost_elements_3_2_2, filename );
  num_errors += RUN_TEST( test_ghost_elements_3_0_1, filename );
  num_errors += RUN_TEST( test_ghost_elements_2_0_1, filename );
  num_errors += RUN_TEST( test_ghost_tag_exchange, filename );
  num_errors += RUN_TEST( regression_ghost_tag_exchange_no_default, filename );
  num_errors += RUN_TEST( test_interface_owners, filename );
  num_errors += RUN_TEST( regression_owners_with_ghosting, filename );
  num_errors += RUN_TEST( test_ghosted_entity_shared_data, filename );
  num_errors += RUN_TEST( test_assign_global_ids, filename );
  num_errors += RUN_TEST( test_shared_sets, 0 );
  num_errors += RUN_TEST( test_reduce_tags, 0);
  
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

ErrorCode get_sharing_processors( Interface& moab,
                                    EntityHandle entity,
                                    std::vector<int>& other_procs_out )
{
  ErrorCode rval;
  
    // get tags for parallel data
  Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  rval = moab.tag_get_handle( PARALLEL_SHARED_PROC_TAG_NAME, 1, MB_TYPE_INTEGER, sharedp_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( PARALLEL_SHARED_PROCS_TAG_NAME, MAX_SHARING_PROCS, MB_TYPE_INTEGER, sharedps_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( PARALLEL_SHARED_HANDLE_TAG_NAME, 1, MB_TYPE_HANDLE, sharedh_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( PARALLEL_SHARED_HANDLES_TAG_NAME, MAX_SHARING_PROCS, MB_TYPE_HANDLE, sharedhs_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( PARALLEL_STATUS_TAG_NAME, 1, MB_TYPE_OPAQUE, pstatus_tag ); CHKERR(rval);
  
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

int is_any_proc_error( int is_my_error )
{
  int result = 0;
  int err = MPI_Allreduce( &is_my_error, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  return err || result;
}


ErrorCode parallel_create_mesh( Interface& mb,
                                int vtx_ids[9],
                                EntityHandle vtx_handles[9],
                                Range& range,
                                EntityHandle* entity_sets )
{
    // Each processor will create four quads.
    // Groups of four quads will be arranged as follows:
    // +------+------+------+------+------+-----
    // |             |             |
    // |             |             |
    // +    Proc 0   +    Proc 2   +    Proc 4
    // |             |             |
    // |             |             |
    // +------+------+------+------+------+-----
    // |             |             |
    // |             |             |
    // +    Proc 1   +    Proc 3   +    Proc 5
    // |             |             |
    // |             |             |
    // +------+------+------+------+------+-----
    //
    // Vertices will be enumerated as follows:
    // 1------6-----11-----16-----21-----26-----
    // |             |             |
    // |             |             |
    // 2      7     12     17     22     27
    // |             |             |
    // |             |             |
    // 3------8-----13-----18-----23-----28-----
    // |             |             |
    // |             |             |
    // 4      9     14     19     24     29
    // |             |             |
    // |             |             |
    // 5-----10-----15-----20-----25-----30-----
    //
    // Element IDs will be [4*rank+1,4*rank+5]
    
  int size, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  const int first_vtx_id = 10*(rank/2) + 2*(rank%2) + 1;
  const double x = 2.0*(rank/2);
  const double y = 2.0*(rank%2);

  // create vertices
  const int idoff = (size%2 && rank/2 == size/2) ? 0 : 2;
  const int idoff1 = rank ? 2 : idoff;
  const int idoff2 = idoff1+idoff;
  const int ids[9] = { first_vtx_id    , first_vtx_id + 3 + idoff1, first_vtx_id + 6 + idoff2,
                       first_vtx_id + 1, first_vtx_id + 4 + idoff1, first_vtx_id + 7 + idoff2,
                       first_vtx_id + 2, first_vtx_id + 5 + idoff1, first_vtx_id + 8 + idoff2 };
  memcpy( vtx_ids, ids, sizeof(ids) );
  const double coords[27] = {  x, y,   0,   x+1, y,   0,   x+2, y,   0,
                               x, y+1, 0,   x+1, y+1, 0,   x+2, y+1, 0,
                               x, y+2, 0,   x+1, y+2, 0,   x+2, y+2, 0 };

  ErrorCode rval;
  Tag id_tag;

  rval = mb.create_vertices( coords, 9, range ); CHKERR(rval);
  assert(range.size() == 9);
  std::copy( range.begin(), range.end(), vtx_handles );
  range.clear();
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHKERR(rval);
  rval = mb.tag_set_data( id_tag, vtx_handles, 9, &ids ); CHKERR(rval);

  const EntityHandle conn[4][4] = { 
                         { vtx_handles[0], vtx_handles[3], vtx_handles[4], vtx_handles[1] },
                         { vtx_handles[1], vtx_handles[4], vtx_handles[5], vtx_handles[2] },
                         { vtx_handles[3], vtx_handles[6], vtx_handles[7], vtx_handles[4] },
                         { vtx_handles[4], vtx_handles[7], vtx_handles[8], vtx_handles[5] } };
  for (int i = 0; i < 4; ++i) {
    const int id = 4*rank + i + 1;
    EntityHandle h;
    rval = mb.create_element( MBQUAD, conn[i], 4, h ); CHKERR(rval);
    range.insert(h);
    rval = mb.tag_set_data( id_tag, &h, 1, &id ); CHKERR(rval);
  }
  
  if (!entity_sets)
    return MB_SUCCESS;
  
  int set_ids[3] = { size+1, rank/2, rank/2+1 };
  int nsets = 0;
  rval = mb.create_meshset( MESHSET_SET, entity_sets[nsets++] ); CHKERR(rval);
  rval = mb.create_meshset( MESHSET_SET, entity_sets[nsets++] ); CHKERR(rval);
  entity_sets[nsets] = 0;
  if (rank < 2) { // if first (left-most) column
    set_ids[1] = set_ids[2];
  }
  else if (rank/2 < (size-1)/2) { // if not last (right-most) column
    rval = mb.create_meshset( MESHSET_SET, entity_sets[nsets++] ); CHKERR(rval);
  }
  for (int i = 0; i < nsets; ++i) {
    rval = mb.add_entities( entity_sets[0], range ); CHKERR(rval);
  }
  rval = mb.tag_set_data( id_tag, entity_sets, nsets, set_ids ); CHKERR(rval);
  
  return MB_SUCCESS;
}

/**************************************************************************
                           Test  Implementations
 **************************************************************************/

ErrorCode test_elements_on_several_procs( const char* filename )
{
  Core mb_instance;
  Interface& moab = mb_instance;
  ErrorCode rval;
  const char* geom_names[] = { "vertex", "curve", "surface", "volume", "unknown" };
  
  rval = moab.load_file( filename, 0, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS" );
  CHKERR(rval);
  
    // test contents of interface sets against sharedEnts structure in pcomm;
  int my_error = 0;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab, 0);
  rval = pcomm->check_all_shared_handles();
  if (MB_SUCCESS != rval) {
    my_error = 1;
    std::cerr << "check_all_shared_handles test failed on proc " 
              << pcomm->proc_config().proc_rank() << std::endl;
  }
  PCHECK(!my_error);
  
    // check adjacencies just to make sure they're consistent
  rval = mb_instance.check_adjacencies();
  if (MB_SUCCESS != rval) my_error = 1;
  PCHECK(!my_error);

  Tag geom_tag, id_tag;
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHKERR(rval);  
  
    // Because we partitioned based on geometric volume sets, then for each
    // lower-dimension geometric entity set all of the contained entities of
    // the same dimension must be shared by the same set of processors
  Range shared, invalid;
  for (int dim = 2; dim > 0; --dim) {
    Range geom_sets;
    const void* tagvals[] = { &dim };
    rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag, tagvals, 1, geom_sets );
    CHKERR(rval);
    
    for (Range::iterator j = geom_sets.begin(); j != geom_sets.end(); ++j) {
      Range ents;
      rval = moab.get_entities_by_dimension( *j, dim, ents ); CHKERR(rval);
      if (ents.empty()) continue;
      
        // get a single sharing list to compare with
      Range::iterator k = ents.begin();
      std::vector<int> procs;
      rval = get_sharing_processors( moab, *k, procs ); CHKERR(rval);
      std::sort( procs.begin(), procs.end() );
      if (procs.size() > 1)
        shared.merge( ents );
      
        // compare other elements
      for (++k; k != ents.end(); ++k) {
        std::vector<int> tmp_procs;
        rval = get_sharing_processors( moab, *k, tmp_procs ); CHKERR(rval);
        std::sort( tmp_procs.begin(), tmp_procs.end() );
        if (tmp_procs != procs)
          invalid.insert( *j );
      }
      
        // compare interior vertices
      ents.clear();
      rval = moab.get_entities_by_dimension( *j, 0, ents ); CHKERR(rval);
      for (k = ents.begin(); k != ents.end(); ++k) {
        std::vector<int> tmp_procs;
        rval = get_sharing_processors( moab, *k, tmp_procs ); CHKERR(rval);
        if (tmp_procs != procs)
          invalid.insert( *j );
      }
    }
  }
  
    // Report any geometric sets for which sharing was inconsistent
  if (!invalid.empty()) {
    std::cerr << "Elements or vertices owned by a single geometric entity are "
              << "not shared by the same set of processors for the "
              << "following geometric entities on process " << pcomm->proc_config().proc_rank()
              << ": ";
    for (Range::iterator i = invalid.begin(); i != invalid.end(); ++i) {
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
  PCHECK(!my_error);
  
    
    // for each shared element, check that its vertices are shared
    // with at least the same processes
  for (Range::iterator i = shared.begin(); i != shared.end(); ++i) {
    std::vector<int> procs;
    rval = get_sharing_processors( moab, *i, procs ); CHKERR(rval);
    std::sort( procs.begin(), procs.end() );
    
    std::vector<EntityHandle> tmp;
    const EntityHandle* conn;
    int len;
    rval = moab.get_connectivity( *i, conn, len, false, &tmp ); CHKERR(rval);
    for (int j = 0; j < len; ++j) {
      std::vector<int> vprocs;
      rval = get_sharing_processors( moab, conn[j], vprocs ); CHKERR(rval);
      std::sort( vprocs.begin(), vprocs.end() );
      std::vector<int> diff( std::max( procs.size(), vprocs.size() ) );
      std::vector<int>::iterator k =
        std::set_difference( procs.begin(), procs.end(), vprocs.begin(), vprocs.end(), diff.begin() );
      if (k != diff.begin()) // difference is not empty
        invalid.insert( conn[j] );
    }
  }
  
    // Report any vertices with insufficient sharing
  if (!invalid.empty()) {
    std::cerr << "Vertices must be shared with at least the union of the processes "
              << "sharing the elements containing the vertex.  This is NOT true for "
              << "the following vertices on process " << pcomm->proc_config().proc_rank() 
              << ": " << invalid << std::endl;
    
    my_error = 1;
  }
  PCHECK(!my_error);

  return MB_SUCCESS;
}


ErrorCode get_ghost_entities( ParallelComm& pcomm,
                                Range& ghost_ents )
{
  Range all_ents;
  ErrorCode rval;
  
  
  rval = pcomm.get_moab()->get_entities_by_handle( 0, all_ents );
  CHKERR(rval);
  std::vector<unsigned char> flags(all_ents.size());
  rval = pcomm.get_moab()->tag_get_data( pcomm.pstatus_tag(), all_ents, &flags[0] );
  CHKERR(rval);
  
  Range::iterator ins = ghost_ents.begin();
  std::vector<unsigned char>::const_iterator f = flags.begin();
  for (Range::iterator i = all_ents.begin(); i != all_ents.end(); ++i, ++f) 
    if ((*f & PSTATUS_NOT_OWNED) && !(*f & PSTATUS_INTERFACE))
      ins = ghost_ents.insert( ins, *i, *i );
  
  return MB_SUCCESS;
}


ErrorCode get_ents_from_geometric_sets( Interface& moab,
                                        const Tag tags[2],
                                        int dimension,
                                        const std::vector<int>& ids,
                                        Range& results )
{
  ErrorCode rval;
  for (size_t i = 0; i < ids.size(); ++i) {
    const void* tag_vals[2] = { &dimension, &ids[i] };
    Range sets;
    rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET,
                                              tags, tag_vals, 2, 
                                              sets ); CHKERR(rval);
    for (Range::iterator j = sets.begin(); j != sets.end(); ++j) {
      Range ents;
      rval = moab.get_entities_by_dimension( *j, dimension, ents ); CHKERR(rval);
      results.merge( ents );
    }
  }
  return MB_SUCCESS;
}

/**\brief Given entire file and IDs of geometric sets, get expected ghost entities
 *
 *\param moab               The entire mesh, loaded in serial
 *\param partition_geom_ids IDs of: geometric volumes owned by this proc and interface topology
 *\param ghost_entity_ids   output list
 */
ErrorCode get_expected_ghosts( Interface& moab,
                                 const std::vector<int> partition_geom_ids[4],
                                 std::vector<int>& ghost_entity_ids,
                                 int ghost_dimension,
                                 int bridge_dimension,
                                 int num_layers )
{
  ErrorCode rval;
  Tag tags[2];
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, tags[0] ); CHKERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, tags[1] ); CHKERR(rval);  

    // get face entities defining interface
  Range skin;
  rval = get_ents_from_geometric_sets( moab, tags, 2, partition_geom_ids[2], skin ); CHKERR(rval);
  
    // get adjacent entities
  Range iface_ghosts, iface_ents;
  if (bridge_dimension == 2) {
    iface_ents = skin;
  }
  else {
    rval = moab.get_adjacencies( skin, bridge_dimension, true, iface_ents, Interface::UNION ); 
    CHKERR(rval);
  }
  for (int n = 0; n < num_layers; ++n) {
    iface_ghosts.clear();
    rval = moab.get_adjacencies( iface_ents, ghost_dimension, false, iface_ghosts, Interface::UNION ); CHKERR(rval);
    iface_ents.clear();
    rval = moab.get_adjacencies( iface_ghosts, bridge_dimension, true, iface_ents, Interface::UNION ); CHKERR(rval);
  }
   
    // get volume entities owned by this process
  Range ents;
  rval = get_ents_from_geometric_sets( moab, tags, 3, partition_geom_ids[3], ents ); CHKERR(rval);
 
    // get entities of ghost_dimension owned by this process
  Range owned;
  if (ghost_dimension == 3)
    owned = ents;
  else {
    rval = moab.get_adjacencies( ents, ghost_dimension, false, owned, Interface::UNION ); 
    CHKERR(rval);
  }
  
    // remove owned entities from ghost list
  Range ghosts = subtract( iface_ghosts, owned );

    // get ids
  ghost_entity_ids.resize( ghosts.size() );
  rval = moab.tag_get_data( tags[1], ghosts, &ghost_entity_ids[0] );
  CHKERR(rval);
  return MB_SUCCESS;
}

ErrorCode test_ghost_elements( const char* filename,
                                 int ghost_dimension,
                                 int bridge_dimension,
                                 int num_layers )
{
  Core mb_instance;
  Interface& moab = mb_instance;
  ErrorCode rval;

  std::ostringstream file_opts;
  file_opts << "PARALLEL=READ_DELETE;"
            << "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
            << "PARTITION_DISTRIBUTE;"
            << "PARALLEL_RESOLVE_SHARED_ENTS;"
            << "PARALLEL_GHOSTS=" 
            << ghost_dimension << '.'
            << bridge_dimension << '.'
            << num_layers;
  
  rval = moab.load_file( filename, 0, file_opts.str().c_str() );
  CHKERR(rval);
  Tag geom_tag, id_tag;
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHKERR(rval);  
  
    // Get partition sets
  Range partition_geom[4];
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab, 0);
  partition_geom[3] = pcomm->partition_sets();
  PCHECK( !partition_geom[3].empty() );

  rval = pcomm->check_all_shared_handles();
  CHKERR(rval);
  
    // exchange id tags to allow comparison by id
  Range tmp_ents;
  rval = pcomm->get_shared_entities(-1, tmp_ents, -1, false, true);
  rval = pcomm->exchange_tags(id_tag, tmp_ents);
  CHKERR(rval);
  
    // Get geometric surfaces
  Range surfs, tmp;
  for (Range::iterator i = partition_geom[3].begin(); i != partition_geom[3].end(); ++i) {
    tmp.clear();
    rval = moab.get_child_meshsets( *i, tmp ); CHKERR(rval);
    surfs.merge( tmp );
  }
  
    // Get the set of geometric surfaces that represent the skin
    // of the union of the parts for this processor.  As we partition
    // based on geometric volumes, the interface must be represented
    // by some subset of the surfaces and their child geometric topology.
  
  int error = 0;
  std::ostringstream error_msg;
  Range ents, iface_surfs, iface_curves, iface_vertices;
  for (Range::iterator i = surfs.begin(); !error && i != surfs.end(); ++i) {
    ents.clear();
    rval = moab.get_entities_by_dimension( *i, ghost_dimension-1, ents ); CHKERR(rval);
    if (ents.empty())
      continue;
    
    std::vector<int> procs, tmp_procs;
    Range::iterator j = ents.begin();
    rval = get_sharing_processors( moab, *j, procs ); CHKERR(rval);
    for (++j; !error && j != ents.end(); ++j) {
      tmp_procs.clear();
      rval = get_sharing_processors( moab, *j, tmp_procs ); CHKERR(rval);
      if( tmp_procs != procs ) {
        error_msg << "Failure at " << __FILE__ << ':' << __LINE__ << std::endl
                  << "\tNot all entities in geometric surface are shared with"
                  << " same processor." << std::endl;
        error = 1;
        break;
      }
    }
    
    if (error)
      break;
    
    // if surface is not part of inter-proc interface, skip it.
    if (procs.empty())
      continue;
    if (procs.size() != 1) {
      error_msg << "Failure at " << __FILE__ << ':' << __LINE__ << std::endl
                  << "\tSurface elements shared with" << procs.size() << "processors." << std::endl;
      error = 1;
      break;
    }
    int other_rank = procs[0];
    if (other_rank == (int)pcomm->proc_config().proc_rank())
      continue;
    
    partition_geom[2].insert( *i );
    ents.clear();
    rval = moab.get_child_meshsets( *i, ents ); CHKERR(rval);
    partition_geom[1].merge( ents );
  }
  
    // Don't to global communication until outside
    // of loops.  Otherwise we will deadlock if not all
    // procs iterate the same number of times.
  if (is_any_proc_error(error)) {
    std::cerr << error_msg.str();
    return MB_FAILURE;
  }
  
  for (Range::iterator i = partition_geom[1].begin(); i != partition_geom[1].end(); ++i) {
    ents.clear();
    rval = moab.get_child_meshsets( *i, ents ); CHKERR(rval);
    partition_geom[0].merge( ents );
  }
  
  std::vector<int> partn_geom_ids[4];
  for (int dim = 0; dim <= 3; ++dim) {
    partn_geom_ids[dim].resize( partition_geom[dim].size() );
    rval = moab.tag_get_data( id_tag, partition_geom[dim], &(partn_geom_ids[dim][0]) );
    CHKERR(rval);
  }
  
    // get the global IDs of the ghosted entities
  Range ghost_ents;
  rval = get_ghost_entities( *pcomm, ghost_ents ); CHKERR(rval);
  std::pair<Range::iterator,Range::iterator> vtx = ghost_ents.equal_range(MBVERTEX);
  ghost_ents.erase( vtx.first, vtx.second );
  std::vector<int> actual_ghost_ent_ids(ghost_ents.size());
  rval = moab.tag_get_data( id_tag, ghost_ents, &actual_ghost_ent_ids[0] ); CHKERR(rval);
  
    // read file in serial
  Core moab2;
  rval = moab2.load_file( filename );  
  PCHECK(MB_SUCCESS == rval);
  
    // get the global IDs of the entities we expect to be ghosted
  std::vector<int> expected_ghost_ent_ids;
  rval = get_expected_ghosts( moab2, partn_geom_ids, expected_ghost_ent_ids,  
                             ghost_dimension, bridge_dimension, num_layers );
  PCHECK(MB_SUCCESS == rval);
  
    // check that the correct entities were ghosted
  std::sort( actual_ghost_ent_ids.begin(), actual_ghost_ent_ids.end() );
  std::sort( expected_ghost_ent_ids.begin(), expected_ghost_ent_ids.end() );
  PCHECK( expected_ghost_ent_ids == actual_ghost_ent_ids );
  
    // check we only have the partitioned and ghosted entities
    // on this processor.
  Range myents;
  for (Range::iterator i = partition_geom[3].begin(); i != partition_geom[3].end(); ++i) {
    ents.clear();
    rval = moab.get_entities_by_dimension( *i, 3, ents ); CHKERR(rval);
    myents.merge( ents );
  }
  if (ghost_dimension != 3) {
    ents.clear();
    rval = moab.get_adjacencies( myents, ghost_dimension, false, ents, Interface::UNION );
    PCHECK(MB_SUCCESS == rval);
    myents.swap(ents);
  }
  myents.merge( ghost_ents );
  ents.clear();
  rval = moab.get_entities_by_dimension( 0, ghost_dimension, ents );
  PCHECK( ents == myents );
    
  rval = pcomm->check_all_shared_handles();
  if (MB_SUCCESS != rval) error = 1;
  PCHECK(!error);
  
    // done
  return MB_SUCCESS;
}


ErrorCode test_ghost_elements_3_2_1( const char* filename )
{
  return test_ghost_elements( filename, 3, 2, 1 );
}

ErrorCode test_ghost_elements_3_2_2( const char* filename )
{
  return test_ghost_elements( filename, 3, 2, 2 );
}

ErrorCode test_ghost_elements_3_0_1( const char* filename )
{
  return test_ghost_elements( filename, 3, 0, 1 );
}

ErrorCode test_ghost_elements_2_0_1( const char* filename )
{
  return test_ghost_elements( filename, 2, 0, 1 );
}

ErrorCode get_owner_handles( ParallelComm* pcomm,
                             const Range& ents,
                             EntityHandle* handle_arr )
{
  size_t i = 0;
  int junk;
  for (Range::iterator j = ents.begin(); j != ents.end(); ++i, ++j) {
    ErrorCode rval = pcomm->get_owner_handle( *j, junk, handle_arr[i] );
    if (MB_SUCCESS != rval)
      return rval;
  }
  return MB_SUCCESS;
}

ErrorCode test_ghost_tag_exchange( const char* filename )
{
  Core mb_instance;
  Interface& moab = mb_instance;
  ErrorCode rval;

  rval = moab.load_file( filename, 0, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS;"
                         "PARALLEL_GHOSTS=3.2.1" );
  CHKERR(rval);
  
    // Get ghost elements
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab, 0);
  Range local, ghosts;
  rval = moab.get_entities_by_dimension( 0, 3, local ); CHKERR(rval);
  Range::iterator i = local.begin();
  while (i != local.end()) {
    int rank;
    rval = pcomm->get_owner( *i, rank ); CHKERR(rval);
    if (rank == (int)pcomm->proc_config().proc_rank()) {
      ++i;
    }
    else {
      ghosts.insert( *i );
      i = local.erase(i);
    }
  }
  
    // create a tag to exchange
  Tag dense_test_tag;
  EntityHandle defval = 0;
  //rval = moab.tag_get_handle( "TEST-TAG", sizeof(EntityHandle), MB_TAG_DENSE,
  //                         dense_test_tag, &defval ); CHKERR(rval);
  rval = moab.tag_get_handle( "TEST-TAG", 1, MB_TYPE_HANDLE, dense_test_tag,
                              MB_TAG_DENSE|MB_TAG_EXCL, &defval ); CHKERR(rval);
    
    // for all entities that I own, set tag to handle value
  std::vector<EntityHandle> handles(local.size()), handles2;
  std::copy( local.begin(), local.end(), handles.begin() );
  rval = moab.tag_set_data( dense_test_tag, local, &handles[0] ); CHKERR(rval);
  
    // exchange tag data
  Range tmp_range;
  rval = pcomm->exchange_tags( dense_test_tag, tmp_range ); CHKERR(rval);
  
    // make sure local values are unchanged
  handles2.resize( local.size() );
  rval = moab.tag_get_data( dense_test_tag, local, &handles2[0] ); CHKERR(rval);
  PCHECK( handles == handles2 );
  
    // compare values on ghost entities
  handles.resize( ghosts.size() );
  handles2.resize( ghosts.size() );
  rval = moab.tag_get_data( dense_test_tag, ghosts, &handles2[0] ); CHKERR(rval);
  rval = get_owner_handles( pcomm, ghosts, &handles[0] ); CHKERR(rval);
  PCHECK( handles == handles2 );


    // now do it all again for a sparse tag
  Tag sparse_test_tag;
  //rval = moab.tag_get_handle( "TEST-TAG-2", sizeof(int), MB_TAG_DENSE,
  //                         MB_TYPE_INTEGER, sparse_test_tag, 0 ); CHKERR(rval);
  rval = moab.tag_get_handle( "TEST-TAG-2", 1, MB_TYPE_INTEGER, sparse_test_tag,
                          MB_TAG_DENSE|MB_TAG_EXCL ); CHKERR(rval);
    
    // for all entiites that I own, set tag to my rank
  std::vector<int> procs1(local.size(), pcomm->proc_config().proc_rank());
  rval = moab.tag_set_data( sparse_test_tag, local, &procs1[0] ); CHKERR(rval);
  
    // exchange tag data
  tmp_range.clear();
  rval = pcomm->exchange_tags( sparse_test_tag, tmp_range ); 
  PCHECK( MB_SUCCESS == rval );
  
    // make sure local values are unchanged
  std::vector<int> procs2(local.size());
  rval = moab.tag_get_data( sparse_test_tag, local, &procs2[0] ); CHKERR(rval);
  PCHECK( procs1 == procs2 );
  
    // compare values on ghost entities
  procs1.resize( ghosts.size() );
  procs2.resize( ghosts.size() );
  rval = moab.tag_get_data( sparse_test_tag, ghosts, &procs2[0] ); CHKERR(rval);
  std::vector<int>::iterator j = procs1.begin();
  for (i = ghosts.begin(); i != ghosts.end(); ++i, ++j) {
    rval = pcomm->get_owner( *i, *j ); 
    CHKERR(rval);
  }
  PCHECK( procs1 == procs2 );
  
  return MB_SUCCESS;
}

ErrorCode regression_ghost_tag_exchange_no_default( const char* filename )
{
  Core mb_instance;
  Interface& moab = mb_instance;
  ErrorCode rval;

  rval = moab.load_file( filename, 0, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS;"
                         "PARALLEL_GHOSTS=3.2.1" );
  CHKERR(rval);
  
    // create a tag to exchange
  Tag dense_test_tag;
  //rval = moab.tag_get_handle( "TEST-TAG", sizeof(EntityHandle), MB_TAG_DENSE,
  //                         dense_test_tag, 0 ); CHKERR(rval);
  rval = moab.tag_get_handle( "TEST-TAG", 1, MB_TYPE_HANDLE, dense_test_tag, 
                              MB_TAG_DENSE|MB_TAG_EXCL); CHKERR(rval);
  
    // exchange tag data
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab, 0);
  Range tmp_range;
  rval = pcomm->exchange_tags( dense_test_tag, tmp_range ); 
  PCHECK(MB_SUCCESS == rval);
  
  return MB_SUCCESS;
}

  
// Helper for exhange_sharing_data
// Swap contens of buffer with specified processor.
int MPI_swap( void* buffer, int num_val, MPI_Datatype val_type, int other_proc )
{
  int err, rank, bytes;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Type_size( val_type, &bytes );
  bytes *= num_val;
  std::vector<unsigned char> buffer2(bytes);
  
  for (int i = 0; i < 2; ++i) {
    if (i == (rank < other_proc)) {
      err = MPI_Send( buffer, num_val, val_type, other_proc, 0, MPI_COMM_WORLD );
      if (err)
        return err;
    }
    else {
      MPI_Status status;
      err = MPI_Recv( &buffer2[0], num_val, val_type, other_proc, 0, MPI_COMM_WORLD, &status );
      if (err)
        return err;
    }
  }
  
  memcpy( buffer, &buffer2[0], bytes );
  return 0;
}

  
int valid_ghosting_owners( int comm_size, 
                           const int* ids,
                           const int* owners )
{
    // for each vertex ID, build list of {rank,owner} tuples
  std::map< int, std::vector<int> > verts;
  for (int p = 0; p < comm_size; ++p) {
    for (int i = 0; i < 9; ++i) { // nine verts per proc
      int idx = 9*p+i;
      verts[ ids[idx] ].push_back( p );
      verts[ ids[idx] ].push_back( owners[idx] );
    }
  }
  
    // now check for each vertex that the owner from
    // each processor is the same
  bool print_desc = true;
  int error_count = 0;
  std::map< int, std::vector<int> >::iterator it;
  for (it = verts.begin(); it != verts.end(); ++it) {
    int id = it->first;
    std::vector<int>& list = it->second;
    bool all_same = true;
    for (size_t i = 2; i < list.size(); i += 2)
      if (list[i+1] != list[1])
        all_same = false;
    if (all_same)
      continue;
    
    ++error_count;
    
    if (print_desc) {
      print_desc = false;
      std::cerr << "ERROR at " __FILE__ ":" << __LINE__ << std::endl
                << "  Processors have inconsistant ideas of vertex ownership:"
                << std::endl;
    }
    
    std::cerr << "  Vertex " << id << ": " << std::endl;
    for (size_t i = 0; i < list.size(); i += 2) 
      std::cerr << "    Proc " << list[i] << " thinks owner is " << list[i+1] << std::endl;
  }

  return error_count;
}

ErrorCode test_interface_owners_common( int num_ghost_layers )
{
  ErrorCode rval;  
  Core moab_instance;
  Interface& mb = moab_instance;
  ParallelComm pcomm( &mb );
  
    // build distributed quad mesh
  Range quads;
  EntityHandle verts[9];
  int ids[9];
  rval = parallel_create_mesh( mb, ids, verts, quads );  PCHECK(MB_SUCCESS == rval);
  rval = pcomm.resolve_shared_ents( 0, quads, 2, 1 ); PCHECK(MB_SUCCESS == rval);
  if (num_ghost_layers) {
    rval = pcomm.exchange_ghost_cells( 2, 0, num_ghost_layers, 0, true ); 
    PCHECK(MB_SUCCESS == rval);
  }
  
    // get vertex owners
  int owner[9];
  for (int i = 0; i < 9; ++i) {
    rval = pcomm.get_owner( verts[i], owner[i] );
    if (MB_SUCCESS != rval)
      break;
  }
  PCHECK(MB_SUCCESS == rval);
  
    // exchange vertex owners amongst all processors
  int rank, size, ierr;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  std::vector<int> all_ids(9*size), all_owner(9*size);
  ierr = MPI_Gather( ids, 9, MPI_INT, &all_ids[0], 9, MPI_INT, 0, MPI_COMM_WORLD );
  if (ierr)
    return MB_FAILURE;
  ierr = MPI_Gather( owner, 9, MPI_INT, &all_owner[0], 9, MPI_INT, 0, MPI_COMM_WORLD );
  if (ierr)
    return MB_FAILURE;
  
  int errors = rank ? 0 : valid_ghosting_owners( size, &all_ids[0], &all_owner[0] );
  MPI_Bcast( &errors, 1, MPI_INT, 0, MPI_COMM_WORLD );
  return errors ? MB_FAILURE : MB_SUCCESS;
}

// Common implementation for both:
//   test_interface
//   regression_interface_with_ghosting
ErrorCode test_interface_owners( const char* )
{
  return test_interface_owners_common(0);
}

ErrorCode regression_owners_with_ghosting( const char* )
{
  return test_interface_owners_common(1);
}


struct VtxData {
  std::vector<int> procs;
  std::vector<int> owners;
  std::vector<EntityHandle> handles;
};

ErrorCode test_ghosted_entity_shared_data( const char* )
{
  ErrorCode rval;  
  Core moab_instance;
  Interface& mb = moab_instance;
  ParallelComm pcomm( &mb );

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
   
    // build distributed quad mesh
  Range quads;
  EntityHandle verts[9];
  int ids[9];
  rval = parallel_create_mesh( mb, ids, verts, quads );  PCHECK(MB_SUCCESS == rval);
  rval = pcomm.resolve_shared_ents( 0, quads, 2, 1 ); PCHECK(MB_SUCCESS == rval);
  rval = pcomm.exchange_ghost_cells( 2, 1, 1, 0, true ); 
  PCHECK(MB_SUCCESS == rval);
  
  rval = pcomm.check_all_shared_handles();
  PCHECK(MB_SUCCESS == rval);
  
  return MB_SUCCESS;
}

ErrorCode check_consistent_ids( Interface& mb,
                                const EntityHandle* entities,
                                const int* orig_ids,
                                int num_ents,
                                const char* singular_name,
                                const char* plural_name )
{
  ErrorCode rval;  
  int rank, size, ierr;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHKERR(rval);
  std::vector<int> new_ids(num_ents);
  rval = mb.tag_get_data( id_tag, entities, num_ents, &new_ids[0] ); CHKERR(rval);
  // This test is wrong.  a) The caller can select a start ID so there's
  // no guarantee that the IDs will be in any specific range and b) There
  // is no reason to expect the global number of entities to be num_ents*size
  // if there are any shared entities. J.Kraftcheck 2010-12-22
  //rval = MB_SUCCESS;
  //for (int i = 0; i < num_ents; ++i) 
  //  if (new_ids[i] < 0 || new_ids[i] >= num_ents*size) {
  //    std::cerr << "ID out of bounds on proc " << rank 
  //              << " : " << new_ids[i] << " not in [0," << num_ents*size-1
  //              << "]" << std::endl;
  //    rval = MB_FAILURE;
  //  }
  //if (MB_SUCCESS != rval)
  //  return rval;

    // Gather up all data on root proc for consistency check
  std::vector<int> all_orig_ids(num_ents*size), all_new_ids(num_ents*size);
  ierr = MPI_Gather( (void*)orig_ids, num_ents, MPI_INT, &all_orig_ids[0], num_ents, MPI_INT, 0, MPI_COMM_WORLD );
  if (ierr)
    return MB_FAILURE;
  ierr = MPI_Gather( &new_ids[0], num_ents, MPI_INT, &all_new_ids[0], num_ents, MPI_INT, 0, MPI_COMM_WORLD );
  if (ierr)
    return MB_FAILURE;
  
    // build a local map from original ID to new ID and use it
    // to check for consistancy between all procs
  rval = MB_SUCCESS;;
  if (0 == rank) {
      // check for two processors having different global ID for same entity
    std::map<int,int> idmap; // index by original ID and contains new ID
    std::map<int,int> owner; // index by original ID and contains owning rank
    for (int i = 0; i < num_ents*size; ++i) {
      std::map<int,int>::iterator it = idmap.find(all_orig_ids[i]);
      if (it == idmap.end()) {
        idmap[all_orig_ids[i]] = all_new_ids[i];
        owner[all_orig_ids[i]] = i/num_ents;
      }
      else if (it->second != all_new_ids[i]) {
        std::cerr << "Inconsistant " << singular_name << " IDs between processors "
                  << owner[all_orig_ids[i]] << " and " << i/num_ents 
                  << " : " << it->second << " and "
                  << all_new_ids[i] << " respectively." << std::endl;
        rval = MB_FAILURE;
      }
    }
      // check for two processors having same global ID for different entities
    idmap.clear();
    owner.clear();
    for (int i = 0; i < num_ents*size; ++i) {
      std::map<int,int>::iterator it = idmap.find(all_new_ids[i]);
      if (it == idmap.end()) {
        idmap[all_new_ids[i]] = all_orig_ids[i];
        owner[all_new_ids[i]] = i/num_ents;
      }
      else if (it->second != all_orig_ids[i]) {
        std::cerr << "ID " << all_new_ids[i] 
                  << " assigned to different " << plural_name << " on processors "
                  << owner[all_new_ids[i]] << " and " << i/num_ents << std::endl;
        rval = MB_FAILURE;
      }
    }
  }
  return rval;
}


ErrorCode test_assign_global_ids( const char* )
{
  ErrorCode rval;  
  Core moab_instance;
  Interface& mb = moab_instance;
  ParallelComm pcomm( &mb );

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
    // build distributed quad mesh
  Range quad_range;
  EntityHandle verts[9];
  int vert_ids[9];
  rval = parallel_create_mesh( mb, vert_ids, verts, quad_range );  PCHECK(MB_SUCCESS == rval);
  rval = pcomm.resolve_shared_ents( 0, quad_range, 2, 1 ); PCHECK(MB_SUCCESS == rval);
  
    // get global ids for quads
  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHKERR(rval);
  assert(4u == quad_range.size());
  EntityHandle quads[4];
  std::copy( quad_range.begin(), quad_range.end(), quads );
  int quad_ids[4];
  rval = mb.tag_get_data( id_tag, quads, 4, quad_ids ); CHKERR(rval);
  
    // clear GLOBAL_ID tag
  int zero[9] = {0};
  rval = mb.tag_set_data( id_tag, verts, 9, zero ); CHKERR(rval);
  rval = mb.tag_set_data( id_tag, quads, 4, zero ); CHKERR(rval);
    
    // assign new global IDs
  rval = pcomm.assign_global_ids( 0, 2 ); PCHECK(MB_SUCCESS == rval);
  
  rval = check_consistent_ids( mb, verts, vert_ids, 9, "vertex", "vertices" );
  PCHECK(MB_SUCCESS == rval);
  rval = check_consistent_ids( mb, quads, quad_ids, 4, "quad", "quads" );
  PCHECK(MB_SUCCESS == rval);
  return rval;
}

ErrorCode test_shared_sets( const char* )
{
  ErrorCode rval;  
  Core moab_instance;
  Interface& mb = moab_instance;
  ParallelComm pcomm( &mb );

  int rank_i, size_i;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank_i );
  MPI_Comm_size( MPI_COMM_WORLD, &size_i );
  const unsigned rank = rank_i;
  const unsigned nproc = size_i;
  
    // build distributed quad mesh
  Range quads, sets;
  EntityHandle verts[9], set_arr[3];
  int ids[9];
  rval = parallel_create_mesh( mb, ids, verts, quads, set_arr );  PCHECK(MB_SUCCESS == rval);
  rval = pcomm.resolve_shared_ents( 0, quads, 2, 1 ); PCHECK(MB_SUCCESS == rval);
  sets.insert( set_arr[0] );
  sets.insert( set_arr[1] );
  if (set_arr[2]) {
    sets.insert( set_arr[2] );
  }
  else {
    set_arr[2] = set_arr[1];
  }
    
  Tag id_tag;
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHKERR(rval);
  rval = pcomm.resolve_shared_sets( sets, id_tag ); PCHECK(MB_SUCCESS == rval);
  
    // check that set data is consistent
  ErrorCode ok = MB_SUCCESS;
  for (size_t i = 0; i < 3; ++i) {
    unsigned owner;
    EntityHandle owner_handle, local;
    rval = pcomm.get_entityset_owner( set_arr[i], owner, &owner_handle );
    if (MB_SUCCESS != rval)
      ok = rval;
    else if (owner == rank) {
      if (owner_handle != set_arr[i]) {
        std::cerr << __FILE__ << ":" << __LINE__ << " rank " << rank <<
                "invalid remote handle for owned set" << std::endl;
        ok = MB_FAILURE;
      }
    }
    else if (MB_SUCCESS != pcomm.get_entityset_local_handle( owner, owner_handle, local ))
      ok = MB_FAILURE;
    else if (local != set_arr[i]) {
      std::cerr << __FILE__ << ":" << __LINE__ << " rank " << rank <<
                "invalid local handle for remote data" << std::endl;
      ok = MB_FAILURE;
    }
  }
  PCHECK(MB_SUCCESS == ok);

  const unsigned col = rank / 2; // which column is proc in
  const unsigned colrank = 2*col;// rank of first of two procs in col (rank if rank is even, rank-1 if rank is odd)
  unsigned mins[3] = { 0, 0, 0 };
  unsigned maxs[3] = { nproc-1, 0, 0 };
  if (rank < 2) { // if in first (left-most) column, then 
    mins[1] = mins[2] = 0;
    maxs[1] = maxs[2] = std::min( 3u, nproc-1 );
  }
  else if (col == (nproc-1)/2) { // else if last (right-most) column, then 
    mins[1] = mins[2] = colrank-2;
    maxs[1] = maxs[2] = std::min( colrank+1, nproc-1 );
  }
  else { // if inside column, then
    mins[1] = colrank-2;
    mins[2] = colrank;
    maxs[1] = std::min( colrank+1, nproc-1 );
    maxs[2] = std::min( colrank+3, nproc-1 );
  }
  
    // check expected sharing lists
    // set 0 is shared between all procs in the range [ mins[0], maxs[0] ]
  std::vector<unsigned> expected, list;
  for (size_t i = 0; i < 3; ++i) {
    expected.clear();
    for (unsigned r = mins[i]; r <= std::min( maxs[i], nproc-1 ); ++r) 
      if (r != rank)
        expected.push_back(r);
    list.clear();
    rval = pcomm.get_entityset_procs( set_arr[i], list );
    if (MB_SUCCESS != rval)
      ok = rval;
    else {
      std::sort( list.begin(), list.end() );
      if (expected != list) {
        std::cerr << __FILE__ << ":" << __LINE__ << " rank " << rank <<
                " incorrect sharing list for entity set" << std::endl;
        ok = MB_FAILURE;
      }
    }
  }
  PCHECK(MB_SUCCESS == ok);

    // check consistent owners 
  unsigned send_list[6], set_owners[3];
  std::vector<unsigned> recv_list(6*nproc);
  for (size_t i = 0; i < 3; ++i) {
    mb.tag_get_data( id_tag, set_arr+i, 1, &send_list[2*i] );
    pcomm.get_entityset_owner( set_arr[i], set_owners[i] );
    send_list[2*i+1] = set_owners[i];
  }
  MPI_Allgather( send_list, 6, MPI_UNSIGNED, &recv_list[0], 6, MPI_UNSIGNED, MPI_COMM_WORLD );
  std::map<unsigned,unsigned> owners;
  for (unsigned i = 0; i < 6*nproc; i+= 2) {
    unsigned id = recv_list[i];
    unsigned owner = recv_list[i+1];
    if (owners.find(id) == owners.end())
      owners[id] = owner;
    else if (owners[id] != owner) {
      std::cerr << __FILE__ << ":" << __LINE__ << " rank " << rank <<
                " mismatched owners (" << owners[id] << " and " << owner
                << ") for set with ID " << id << std::endl;
      ok = MB_FAILURE;
    }
  }
  PCHECK(MB_SUCCESS == ok);
  
    // test PComm::get_entityset_owners
  std::vector<unsigned> act_owners, exp_owners;
  for (size_t i = 0; i < 3; ++i) {
    exp_owners.push_back(set_owners[i]);
  }
  ok = pcomm.get_entityset_owners( act_owners );
  PCHECK(MB_SUCCESS == ok);
  //PCHECK(std::find(act_owners.begin(),act_owners.end(),rank) == act_owners.end());
  std::sort( act_owners.begin(), act_owners.end() );
  std::sort( exp_owners.begin(), exp_owners.end() );
  exp_owners.erase( std::unique( exp_owners.begin(), exp_owners.end() ), exp_owners.end() );
  PCHECK( exp_owners == act_owners );

    // test PComm::get_owned_sets
  ok = MB_SUCCESS;
  for (unsigned i = 0; i < nproc; ++i) {
    sets.clear();
    if (MB_SUCCESS != pcomm.get_owned_sets( i, sets )) {
      std::cerr << __FILE__ << ":" << __LINE__ << " rank " << rank <<
                " failed to get shared set list for sets owned by rank " <<
                set_owners[i] << std::endl;
      continue;
    }
    
    Range expected;
    for (size_t j = 0; j < 3; ++j)
      if (set_owners[j] == i)
        expected.insert( set_arr[j] );
    
    if (expected != sets) {
      std::cerr << __FILE__ << ":" << __LINE__ << " rank " << rank 
                << " has incorrect shared set list for sets owned by rank " 
                << set_owners[i] << std::endl << "Expected: " << expected << std::endl
                << "Actual: " << sets << std::endl;
      ok = MB_FAILURE;
    }
  }
  PCHECK( MB_SUCCESS == ok );

  return MB_SUCCESS;
}

template <typename T> ErrorCode check_shared_ents(ParallelComm &pcomm, Tag tagh, T fact, MPI_Op mpi_op) 
{
    // get the shared entities
  Range shared_ents;
  ErrorCode rval = pcomm.get_shared_entities(-1, shared_ents); CHKERR(rval);
  
  std::vector<int> shprocs(MAX_SHARING_PROCS);
  std::vector<EntityHandle> shhandles(MAX_SHARING_PROCS);
  std::vector<T> dum_vals(shared_ents.size());
  int np;
  unsigned char pstatus;

  rval = pcomm.get_moab()->tag_get_data(tagh, shared_ents, &dum_vals[0]); CHKERR(rval);
  
    // for each, compare number of sharing procs against tag value, should be 1/fact the tag value
  Range::iterator rit;
  int i = 0;
  for (rit = shared_ents.begin(); rit != shared_ents.end(); rit++, i++) {
    rval = pcomm.get_sharing_data(*rit, &shprocs[0], &shhandles[0], pstatus, np); CHKERR(rval);
    if (1 == np && shprocs[0] != (int) pcomm.proc_config().proc_rank()) np++;
    if      (mpi_op == MPI_SUM) {if (dum_vals[i] != fact*np) return MB_FAILURE;}
    else if (mpi_op == MPI_PROD) {if (dum_vals[i] != pow(fact, np)) return MB_FAILURE;}
    else if (mpi_op == MPI_MAX) {if (dum_vals[i] != fact) return MB_FAILURE;}
    else if (mpi_op == MPI_MIN) {if (dum_vals[i] != fact) return MB_FAILURE;}
    else return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}
  
template<typename T> ErrorCode test_reduce_tags( const char*, DataType tp )
{
  ErrorCode rval;  
  Core moab_instance;
  Interface& mb = moab_instance;
  ParallelComm pcomm( &mb );

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
    // build distributed quad mesh
  Range quad_range;
  EntityHandle verts[9];
  int vert_ids[9];
  rval = parallel_create_mesh( mb, vert_ids, verts, quad_range );  PCHECK(MB_SUCCESS == rval);
  rval = pcomm.resolve_shared_ents( 0, quad_range, 2, 1 ); PCHECK(MB_SUCCESS == rval);
  
  Tag test_tag;
  T def_val = 2;
  Range dum_range;
  
    // T, MPI_SUM
  rval = mb.tag_get_handle( "test_tag", 1, tp, test_tag, MB_TAG_DENSE|MB_TAG_CREAT, &def_val ); CHKERR(rval);
  rval = pcomm.reduce_tags(test_tag, MPI_SUM, dum_range); CHKERR(rval);
  rval = check_shared_ents(pcomm, test_tag, (T)2, MPI_SUM); CHKERR(rval);
  rval = mb.tag_delete(test_tag);
  
    // T, MPI_PROD
  rval = mb.tag_get_handle( "test_tag", 1, tp, test_tag, MB_TAG_DENSE|MB_TAG_CREAT, &def_val ); CHKERR(rval);
  rval = pcomm.reduce_tags(test_tag, MPI_PROD, dum_range); CHKERR(rval);
  rval = check_shared_ents(pcomm, test_tag, (T)2, MPI_PROD); CHKERR(rval);
  rval = mb.tag_delete(test_tag);
  
    // T, MPI_MAX
  rval = mb.tag_get_handle( "test_tag", 1, tp, test_tag, MB_TAG_DENSE|MB_TAG_CREAT, &def_val ); CHKERR(rval);
    // get owned entities and set a different value on them
  rval = pcomm.get_shared_entities(-1, dum_range, -1, false, false); CHKERR(rval);
  std::vector<T> dum_vals(dum_range.size(), (T)3);
  rval = mb.tag_set_data(test_tag, dum_range, &dum_vals[0]); CHKERR(rval);
  rval = pcomm.reduce_tags(test_tag, MPI_MAX, dum_range); CHKERR(rval);
  rval = check_shared_ents(pcomm, test_tag, (T)3, MPI_MAX); CHKERR(rval);
  rval = mb.tag_delete(test_tag);
  
    // T, MPI_MIN
  rval = mb.tag_get_handle( "test_tag", 1, tp, test_tag, MB_TAG_DENSE|MB_TAG_CREAT, &def_val ); CHKERR(rval);
    // get owned entities and set a different value on them
  std::fill(dum_vals.begin(), dum_vals.end(), (T)-1);
  rval = mb.tag_set_data(test_tag, dum_range, &dum_vals[0]); CHKERR(rval);
  rval = pcomm.reduce_tags(test_tag, MPI_MIN, dum_range); CHKERR(rval);
  rval = check_shared_ents(pcomm, test_tag, (T)-1, MPI_MIN); CHKERR(rval);
  rval = mb.tag_delete(test_tag);
  
  return MB_SUCCESS;
}

ErrorCode test_reduce_tags( const char*)
{
  ErrorCode rval = MB_SUCCESS, tmp_rval;
  
  tmp_rval = test_reduce_tags<int>("test_reduce_tags", MB_TYPE_INTEGER);
  if (MB_SUCCESS != tmp_rval) {
    std::cout << "test_reduce_tags failed for int data type." << std::endl;
    rval = tmp_rval;
  }
  
  tmp_rval = test_reduce_tags<double>("test_reduce_tags", MB_TYPE_DOUBLE);
  if (MB_SUCCESS != tmp_rval) {
    std::cout << "test_reduce_tags failed for double data type." << std::endl;
    rval = tmp_rval;
  }
  
  return rval;
}

