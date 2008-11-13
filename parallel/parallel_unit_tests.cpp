#include "MBParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "mpi.h"
#include <iostream>
#include <algorithm>
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

#define PCHECK(A) if (is_any_proc_error(!(A))) return report_error(__FILE__,__LINE__)

MBErrorCode report_error( const char* file, int line )
{
  std::cerr << "Failure at " << file << ':' << line << std::endl;
  return MB_FAILURE;
}

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

// Test sharing tags for mesh entities shared by more than two processors
MBErrorCode test_elements_on_several_procs( const char* filename );
// Test correct ghosting of elements
MBErrorCode test_ghost_elements_3_2_1( const char* filename );
MBErrorCode test_ghost_elements_3_2_2( const char* filename );
MBErrorCode test_ghost_elements_3_0_1( const char* filename );
MBErrorCode test_ghost_elements_2_0_1( const char* filename );
// Test exchange of tag data on ghost elements
MBErrorCode test_ghost_tag_exchange( const char* filename );
// Bug where exchange_tags fails if dense tag cannot be queried
// for all ghost entities (e.g. no default value)
MBErrorCode regression_ghost_tag_exchange_no_default( const char* filename );

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
      std::cout << func_name << " : FAILED!!" << std::endl;
    else
      std::cout << func_name << " : success" << std::endl;
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
  num_errors += RUN_TEST( test_ghost_elements_3_2_1, filename );
  num_errors += RUN_TEST( test_ghost_elements_3_2_2, filename );
  num_errors += RUN_TEST( test_ghost_elements_3_0_1, filename );
  num_errors += RUN_TEST( test_ghost_elements_2_0_1, filename );
  num_errors += RUN_TEST( test_ghost_tag_exchange, filename );
  num_errors += RUN_TEST( regression_ghost_tag_exchange_no_default, filename );
  
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
  PCHECK(!my_error);
  
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
  PCHECK(!my_error);
  
  return MB_SUCCESS;
}


MBErrorCode get_ghost_entities( MBParallelComm& pcomm,
                                MBRange& ghost_ents )
{
  MBRange all_ents;
  MBErrorCode rval;
  
  
  rval = pcomm.get_moab()->get_entities_by_handle( 0, all_ents );
  CHKERR(rval);
  std::vector<unsigned char> flags(all_ents.size());
  rval = pcomm.get_moab()->tag_get_data( pcomm.pstatus_tag(), all_ents, &flags[0] );
  CHKERR(rval);
  
  MBRange::iterator ins = ghost_ents.begin();
  std::vector<unsigned char>::const_iterator f = flags.begin();
  for (MBRange::iterator i = all_ents.begin(); i != all_ents.end(); ++i, ++f) 
    if ((*f & PSTATUS_NOT_OWNED) && !(*f & PSTATUS_INTERFACE))
      ins = ghost_ents.insert( ins, *i, *i );
  
  return MB_SUCCESS;
}

MBErrorCode get_expected_ghosts( MBInterface& moab,
                                 const std::vector<int> partition_geom_ids[4],
                                 std::vector<int>& ghost_entity_ids,
                                 int ghost_dimension,
                                 int bridge_dimension,
                                 int num_layers )
{
  MBErrorCode rval;
  MBTag tags[2];
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, tags[0] ); CHKERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, tags[1] ); CHKERR(rval);  

    // get all interface sets
  MBRange iface_sets;
  for (int d = 0; d < 3; ++d) {
    for (size_t i = 0; i < partition_geom_ids[d].size(); ++i) {
        // get the entity set
      const void* tag_vals[2] = { &d, &(partition_geom_ids[d][i]) };
      MBRange ents;
      rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                tags, tag_vals, 2, 
                                                ents ); CHKERR(rval);
      iface_sets.merge( ents );
    }
  }
    // get all interface entities
  MBRange all_iface_ents;
  for (MBRange::iterator i = iface_sets.begin(); i != iface_sets.end(); ++i) {
    rval = moab.get_entities_by_handle( *i, all_iface_ents ); 
    CHKERR(rval);
  }
  
    // for each interface set
  MBRange ghosts;
  for (MBRange::iterator i = iface_sets.begin(); i != iface_sets.end(); ++i) {
    if (num_layers < 1)
      break;
    
      // get iface dim
    int gdim = -1;
    rval = moab.tag_get_data( tags[0], &*i, 1, &gdim );
    CHKERR(rval);
     
      // get partitions adjacent to interface set
    MBRange parents; parents.insert(*i);
    for (int step = gdim; step < 3; ++step) {
      MBRange old_parents;
      old_parents.swap(parents);
      for (MBRange::iterator p = old_parents.begin(); p != old_parents.end(); ++p) {
        rval = moab.get_parent_meshsets( *p, parents ); CHKERR(rval);
      }
    }
    
      // get entities in adjacent partitions 
    MBRange adj_proc_ents;
    for (MBRange::iterator p = parents.begin(); p != parents.end(); ++p) {
      int id;
      rval = moab.tag_get_data( tags[1], &*p, 1, &id ); CHKERR(rval);
      if (std::find(partition_geom_ids[3].begin(), partition_geom_ids[3].end(), id) 
          != partition_geom_ids[3].end())
        continue;
      rval = moab.get_entities_by_dimension( *p, 3, adj_proc_ents );
      CHKERR(rval);
    }
    
      // get sub-entities implicitly within partitions
    MBRange tmprange;
    for (int d = 2; d >= 0; --d) {
      rval = moab.get_adjacencies( adj_proc_ents, d, false, tmprange, MBInterface::UNION );
      CHKERR(rval);
    }
    adj_proc_ents.merge( tmprange.subtract( all_iface_ents ) );
    
      // get adjacent entities
    MBRange iface_ghosts, iface_ents;
    rval = moab.get_entities_by_dimension( *i, bridge_dimension, iface_ents ); CHKERR(rval);
    for (int n = 0; n < num_layers; ++n) {
      iface_ghosts.clear();
      rval = moab.get_adjacencies( iface_ents, ghost_dimension, false, iface_ghosts, MBInterface::UNION ); CHKERR(rval);
      iface_ents.clear();
      rval = moab.get_adjacencies( iface_ghosts, bridge_dimension, false, iface_ents, MBInterface::UNION ); CHKERR(rval);
    }
    
      // intersect with entities in adjacent partitions
    ghosts.merge( iface_ghosts.intersect( adj_proc_ents ) );
  }
  
    // get ids
  ghost_entity_ids.resize( ghosts.size() );
  rval = moab.tag_get_data( tags[1], ghosts, &ghost_entity_ids[0] );
  CHKERR(rval);
  return MB_SUCCESS;
}

MBErrorCode test_ghost_elements( const char* filename,
                                 int ghost_dimension,
                                 int bridge_dimension,
                                 int num_layers )
{
  MBCore mb_instance;
  MBInterface& moab = mb_instance;
  MBErrorCode rval;
  MBEntityHandle set;

  std::ostringstream file_opts;
  file_opts << "PARALLEL=READ_DELETE;"
            << "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
            << "PARTITION_DISTRIBUTE;"
            << "PARALLEL_RESOLVE_SHARED_ENTS;"
            << "PARALLEL_GHOSTS=" 
            << ghost_dimension << '.'
            << bridge_dimension << '.'
            << num_layers;
  
  rval = moab.load_file( filename, set, file_opts.str().c_str() );
  CHKERR(rval);
  MBTag geom_tag, id_tag;
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, geom_tag ); CHKERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag ); CHKERR(rval);  
  
    // Get partition sets
  MBRange partition_geom[4];
  MBParallelComm* pcomm = MBParallelComm::get_pcomm(&moab, 0);
  partition_geom[3] = pcomm->partition_sets();
  PCHECK( !partition_geom[3].empty() );
  
    // Get geometric surfaces
  MBRange surfs, tmp;
  for (MBRange::iterator i = partition_geom[3].begin(); i != partition_geom[3].end(); ++i) {
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
  MBRange ents, iface_surfs, iface_curves, iface_vertices;
  for (MBRange::iterator i = surfs.begin(); !error && i != surfs.end(); ++i) {
    ents.clear();
    rval = moab.get_entities_by_dimension( *i, ghost_dimension-1, ents ); CHKERR(rval);
    if (ents.empty())
      continue;
    
    std::vector<int> procs, tmp_procs;
    MBRange::iterator j = ents.begin();
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
  
  for (MBRange::iterator i = partition_geom[1].begin(); i != partition_geom[1].end(); ++i) {
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
  MBRange ghost_ents;
  rval = get_ghost_entities( *pcomm, ghost_ents ); CHKERR(rval);
  std::pair<MBRange::iterator,MBRange::iterator> vtx = ghost_ents.equal_range(MBVERTEX);
  ghost_ents.erase( vtx.first, vtx.second );
  std::vector<int> actual_ghost_ent_ids(ghost_ents.size());
  rval = moab.tag_get_data( id_tag, ghost_ents, &actual_ghost_ent_ids[0] ); CHKERR(rval);
  
    // read file in serial
  MBCore moab2;
  MBEntityHandle set2;
  rval = moab2.load_file( filename, set2 );  
  PCHECK(MB_SUCCESS == rval);
  
    // get the global IDs of teh entities we expect to be ghosted
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
  MBRange myents;
  for (MBRange::iterator i = partition_geom[3].begin(); i != partition_geom[3].end(); ++i) {
    ents.clear();
    rval = moab.get_entities_by_dimension( *i, 3, ents ); CHKERR(rval);
    myents.merge( ents );
  }
  if (ghost_dimension != 3) {
    ents.clear();
    rval = moab.get_adjacencies( myents, ghost_dimension, false, ents, MBInterface::UNION );
    PCHECK(MB_SUCCESS == rval);
    myents.swap(ents);
  }
  myents.merge( ghost_ents );
  ents.clear();
  rval = moab.get_entities_by_dimension( 0, ghost_dimension, ents );
  PCHECK( ents == myents );
    
    // Verify correct ownership and sharing of ghost entities.
  ents.clear();
  for (MBRange::iterator i = myents.begin(); i != myents.end(); ++i) {
    int owner;
    rval = pcomm->get_owner( *i, owner ); CHKERR(rval);
    if ((unsigned)owner == pcomm->proc_config().proc_rank())
      ents.insert( *i );
  }
  myents.swap(ents);
  std::vector<int> my_ent_ids(ents.size());
  rval = moab.tag_get_data( id_tag, myents, &my_ent_ids[0] );
  PCHECK(MB_SUCCESS == rval);
  
  std::sort( my_ent_ids.begin(), my_ent_ids.end() );
  int counts[2] = { my_ent_ids.size(), actual_ghost_ent_ids.size() };
  int totals[2] = {0,0};
  error = MPI_Allreduce( counts, totals, 2, MPI_INT, MPI_SUM,
                         pcomm->proc_config().proc_comm() );
  PCHECK(!error);
  std::vector<int> all_owned(totals[0]), all_ghost(totals[1]), 
                   owned_counts(pcomm->proc_config().proc_size()),
                   owned_displs(pcomm->proc_config().proc_size()),
                   ghost_counts(pcomm->proc_config().proc_size()),
                   ghost_displs(pcomm->proc_config().proc_size());
  error = MPI_Allgather( counts, 1, MPI_INT, &owned_counts[0], 1, MPI_INT,
                         pcomm->proc_config().proc_comm() );
  PCHECK(!error);
  error = MPI_Allgather( counts+1, 1, MPI_INT, &ghost_counts[0], 1, MPI_INT,
                         pcomm->proc_config().proc_comm() );
  PCHECK(!error);
  owned_displs[0] = ghost_displs[0] = 0;
  for (unsigned i = 1; i < pcomm->proc_config().proc_size(); ++i) {
    owned_displs[i] = owned_displs[i-1] + owned_counts[i-1];
    ghost_displs[i] = ghost_displs[i-1] + ghost_counts[i-1];
  }
  
  error = MPI_Allgatherv( &my_ent_ids[0], my_ent_ids.size(), MPI_INT,
                          &all_owned[0], &owned_counts[0], &owned_displs[0],
                          MPI_INT, pcomm->proc_config().proc_comm() );
  PCHECK(!error);
  error = MPI_Allgatherv( &actual_ghost_ent_ids[0], actual_ghost_ent_ids.size(), MPI_INT,
                          &all_ghost[0], &ghost_counts[0], &ghost_displs[0],
                          MPI_INT, pcomm->proc_config().proc_comm() );
  PCHECK(!error);
 
     // for each ghost entity, check owning processor and list of
     // sharing processors.
  int k = 0;
  error = 0;
  for (MBRange::iterator i = ghost_ents.begin(); !error && i != ghost_ents.end(); ++i) {
    std::vector<int> act_procs, exp_procs;
    int act_owner, exp_owner;
    int id = actual_ghost_ent_ids[k++];
    for (unsigned j = 0; j < pcomm->proc_config().proc_size(); ++j) {
      const int* proc_owned_begin = &all_owned[0] + owned_displs[j];
      const int* proc_owned_end = proc_owned_begin + owned_counts[j];
      if (std::binary_search( proc_owned_begin, proc_owned_end, id ))
        exp_owner = j;
      
      const int* proc_ghost_begin = &all_ghost[0] + ghost_displs[j];
      const int* proc_ghost_end = proc_ghost_begin + ghost_counts[j];
      if (std::binary_search( proc_ghost_begin, proc_ghost_end, id ))
        exp_procs.push_back(j);
    }
    
    rval = pcomm->get_owner( *i, act_owner ); CHKERR(rval);
    if (exp_owner != act_owner) {
      error = 1;
      break;
    }
    
    rval = get_sharing_processors( moab, *i, act_procs ); CHKERR(rval);
    std::sort(act_procs.begin(), act_procs.end());
    if (exp_procs != act_procs) {
      error = 1;
      break;
    }
  }
    
    // done
  return MB_SUCCESS;
}


MBErrorCode test_ghost_elements_3_2_1( const char* filename )
{
  return test_ghost_elements( filename, 3, 2, 1 );
}

MBErrorCode test_ghost_elements_3_2_2( const char* filename )
{
  return test_ghost_elements( filename, 3, 2, 2 );
}

MBErrorCode test_ghost_elements_3_0_1( const char* filename )
{
  return test_ghost_elements( filename, 3, 0, 1 );
}

MBErrorCode test_ghost_elements_2_0_1( const char* filename )
{
  return test_ghost_elements( filename, 2, 0, 1 );
}


MBErrorCode test_ghost_tag_exchange( const char* filename )
{
  MBCore mb_instance;
  MBInterface& moab = mb_instance;
  MBErrorCode rval;
  MBEntityHandle set;

  rval = moab.load_file( filename, set, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS;"
                         "PARALLEL_GHOSTS=3.2.1" );
  CHKERR(rval);
  
    // Get ghost elements
  MBParallelComm* pcomm = MBParallelComm::get_pcomm(&moab, 0);
  MBRange local, ghosts;
  rval = moab.get_entities_by_dimension( 0, 3, local ); CHKERR(rval);
  MBRange::iterator i = local.begin();
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
  MBTag dense_test_tag;
  MBEntityHandle defval = 0;
  rval = moab.tag_create( "TEST-TAG", sizeof(MBEntityHandle), MB_TAG_DENSE,
                           dense_test_tag, &defval ); CHKERR(rval);
    
    // for all entiites that I own, set tag to handle value
  std::vector<MBEntityHandle> handles(local.size()), handles2;
  std::copy( local.begin(), local.end(), handles.begin() );
  rval = moab.tag_set_data( dense_test_tag, local, &handles[0] ); CHKERR(rval);
  
    // exchange tag data
  rval = pcomm->exchange_tags( dense_test_tag ); CHKERR(rval);
  
    // make sure local values are unchanged
  handles2.resize( local.size() );
  rval = moab.tag_get_data( dense_test_tag, local, &handles2[0] ); CHKERR(rval);
  PCHECK( handles == handles2 );
  
    // compare values on ghost entities
  handles.resize( ghosts.size() );
  handles2.resize( ghosts.size() );
  rval = moab.tag_get_data( dense_test_tag, ghosts, &handles2[0] ); CHKERR(rval);
  rval = moab.tag_get_data( pcomm->sharedh_tag(), ghosts, &handles[0] ); CHKERR(rval);
  PCHECK( handles == handles2 );


    // now do it all again for a sparse tag
  MBTag sparse_test_tag;
  rval = moab.tag_create( "TEST-TAG-2", sizeof(int), MB_TAG_DENSE,
                           MB_TYPE_INTEGER, sparse_test_tag, 0 ); CHKERR(rval);
    
    // for all entiites that I own, set tag to my rank
  std::vector<int> procs1(local.size(), pcomm->proc_config().proc_rank());
  rval = moab.tag_set_data( sparse_test_tag, local, &procs1[0] ); CHKERR(rval);
  
    // exchange tag data
  rval = pcomm->exchange_tags( sparse_test_tag ); 
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

MBErrorCode regression_ghost_tag_exchange_no_default( const char* filename )
{
  MBCore mb_instance;
  MBInterface& moab = mb_instance;
  MBErrorCode rval;
  MBEntityHandle set;

  rval = moab.load_file( filename, set, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS;"
                         "PARALLEL_GHOSTS=3.2.1" );
  CHKERR(rval);
  
    // create a tag to exchange
  MBTag dense_test_tag;
  rval = moab.tag_create( "TEST-TAG", sizeof(MBEntityHandle), MB_TAG_DENSE,
                           dense_test_tag, 0 ); CHKERR(rval);
  
    // exchange tag data
  MBParallelComm* pcomm = MBParallelComm::get_pcomm(&moab, 0);
  rval = pcomm->exchange_tags( dense_test_tag ); 
  PCHECK(MB_SUCCESS == rval);
  
  return MB_SUCCESS;
}

  
