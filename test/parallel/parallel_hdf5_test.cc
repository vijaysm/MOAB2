#include "moab/Range.hpp"
#include "TestUtil.hpp"

#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/Skinner.hpp"
#include "MBTagConventions.hpp"
#include "moab/CN.hpp"
#include "MBParallelConventions.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include "moab_mpi.h"
#include <unistd.h>
#include <float.h>
#include <stdio.h>
#include <ctime>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;

#ifdef MESHDIR
const char* InputFile = STRINGIFY(MESHDIR) "/ptest.cub";
const char* InputMix = STRINGIFY(MESHDIR) "/io/mix.h5m";
const char* InputOneSide = STRINGIFY(MESHDIR) "/io/oneside.h5m";
#else
#error Specify MESHDIR to compile test
#endif

void load_and_partition( Interface& moab, const char* filename, bool print_debug = false );

void save_and_load_on_root( Interface& moab, const char* tmp_filename );

void check_identical_mesh( Interface& moab1, Interface& moab2 );

void test_write_elements();
void test_write_shared_sets();
void test_var_length_parallel();

void test_read_elements_common( bool by_rank, int intervals, bool print_time,
                                const char* extra_opts = 0 );

int ReadIntervals = 0;
void test_read_elements()         { test_read_elements_common( false, ReadIntervals, false ); }
void test_read_elements_by_rank() { test_read_elements_common(  true, ReadIntervals, false ); }
void test_bcast_summary()         { test_read_elements_common( false, ReadIntervals, false, "BCAST_SUMMARY=yes" ); }
void test_read_summary()          { test_read_elements_common( false, ReadIntervals, false, "BCAST_SUMMARY=no" ); }
void test_read_time();

void test_read_tags();
void test_read_global_tags();
void test_read_sets_common( const char* extra_opts = 0 );
void test_read_sets()             { test_read_sets_common(); }
void test_read_sets_bcast_dups()  { test_read_sets_common( "BCAST_DUPLICATE_READS=yes" ); }
void test_read_sets_read_dups()   { test_read_sets_common( "BCAST_DUPLICATE_READS=no"  ); }
void test_read_bc_sets();

void test_write_different_element_types();
void test_write_different_tags();
void test_write_polygons();
void test_write_unbalanced();
void test_write_dense_tags();
void test_read_non_adjs_side();

const char PARTITION_TAG[] = "PARTITION";

bool KeepTmpFiles = false;
bool PauseOnStart = false;
bool HyperslabAppend = false;
const int DefaultReadIntervals = 2;
int ReadDebugLevel = 0;
int WriteDebugLevel = 0;
int ReadBlocks = 1;

enum Mode { READ_PART, READ_DELETE, BCAST_DELETE };
const Mode DEFAULT_MODE = READ_PART;
const bool DEFAULT_BY_RANK = false;

std::string get_read_options( bool by_rank = DEFAULT_BY_RANK, 
                              Mode mode = DEFAULT_MODE,
                              const char* extra_opts = 0 );
std::string get_read_options( const char* extra_opts )
  { return get_read_options( DEFAULT_BY_RANK, DEFAULT_MODE, extra_opts ); }
std::string get_read_options( bool by_rank, const char* extra_opts )
  { return get_read_options( by_rank, DEFAULT_MODE, extra_opts ); }

std::string get_read_options( bool by_rank, 
                              Mode mode,
                              const char* extra_opts )
{
  int numproc;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  
    // do parallel read unless only one processor
  std::ostringstream str;
  if (numproc > 1) {
    str << "PARALLEL=";
    switch (mode) {
      case READ_PART: str << "READ_PART"; break;
      case READ_DELETE: str << "READ_DELETE"; break;
      case BCAST_DELETE: str << "BCAST_DELETE"; break;
    }
    str << ";PARTITION=" << PARTITION_TAG << ";";
    if (by_rank)
      str << "PARTITION_BY_RANK;";
  }
      
  if (extra_opts)
    str << extra_opts << ";";

  if (ReadDebugLevel)
    str << "DEBUG_IO=" << ReadDebugLevel << ";";

  if (HyperslabAppend)
    str << "HYPERSLAB_APPEND;HYPERSLAB_SELECT_LIMIT=2147483647";

  return str.str();
}

int main( int argc, char* argv[] )
{
  int err = MPI_Init( &argc, &argv );
  CHECK(!err);

  for (int i = 1; i < argc; ++i) {
    if (!strcmp( argv[i], "-k"))
      KeepTmpFiles = true;
    else if (!strcmp( argv[i], "-p"))
      PauseOnStart = true;
    else if (!strcmp( argv[i], "-R")) {
      ++i;
      CHECK( i < argc );
      ReadIntervals = atoi( argv[i] );
      CHECK( ReadIntervals > 0 );
    }
    else if (!strcmp( argv[i], "-b")) {
      ++i;
      CHECK( i < argc );
      ReadBlocks = atoi( argv[i] );
      CHECK( ReadBlocks > 0 );
    }
    else if (!strcmp( argv[i], "-r")) {
      ++i;
      CHECK( i < argc );
      ReadDebugLevel = atoi( argv[i] );
      CHECK( ReadDebugLevel > 0 );
    }
    else if (!strcmp( argv[i], "-w")) {
      ++i;
      CHECK( i < argc );
      WriteDebugLevel = atoi( argv[i] );
      CHECK( WriteDebugLevel > 0 );
    }
    else if (!strcmp( argv[i], "-A")) 
      HyperslabAppend = true;
    else {
      std::cerr << "Usage: " << argv[0] << " [-k] [-p] [-R <intervals>] [-r <level>] [-w <level>] [-b <blocks>]  [-A]" << std::endl;
      return 1;
    }
  }
  
  if (PauseOnStart) {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    printf("Rank %2d PID %lu\n", rank, (unsigned long)getpid());
    sleep(30);
  }
  
  int result = 0;
  if (ReadIntervals) {
    result = RUN_TEST( test_read_time );
  }
  else {
    ReadIntervals = DefaultReadIntervals;
    result += RUN_TEST( test_write_elements );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_write_shared_sets );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_var_length_parallel );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_elements );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_elements_by_rank );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_tags );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_global_tags );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_sets );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_sets_bcast_dups );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_sets_read_dups );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_bc_sets );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_write_different_element_types );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_write_different_tags );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_write_polygons );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_write_unbalanced );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_write_dense_tags );
    MPI_Barrier(MPI_COMM_WORLD);
    result += RUN_TEST( test_read_non_adjs_side );
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Finalize();
  return result;
}

/* Assume geometric topology sets correspond to interface sets 
 * in that the entities contained inclusively in a geometric
 * topology set must be shared by the same set of processors.
 * As we partition based on geometric volumes, this assumption
 * should aways be true.  Print for each geometric topology set
 * the list of processors it is shared with.
 */
void print_partitioned_entities( Interface& moab, bool list_non_shared = false )
{
  ErrorCode rval;
  int size, rank;
  std::vector<int> ent_procs(MAX_SHARING_PROCS), tmp_ent_procs(MAX_SHARING_PROCS);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
    // expect shared entities to correspond to geometric sets
    
    // get tags for parallel data
  Tag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  rval = moab.tag_get_handle( PARALLEL_SHARED_PROC_TAG_NAME, 1, MB_TYPE_INTEGER, sharedp_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( PARALLEL_SHARED_PROCS_TAG_NAME, MAX_SHARING_PROCS, MB_TYPE_INTEGER, sharedps_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( PARALLEL_SHARED_HANDLE_TAG_NAME, 1, MB_TYPE_HANDLE, sharedh_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( PARALLEL_SHARED_HANDLES_TAG_NAME, MAX_SHARING_PROCS, MB_TYPE_HANDLE, sharedhs_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( PARALLEL_STATUS_TAG_NAME, 1, MB_TYPE_OPAQUE, pstatus_tag ); CHECK_ERR(rval);
  
    // for each geometric entity, check which processor we are sharing
    // entities with
  Tag geom_tag, id_tag;
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag ); CHECK_ERR(rval);
  const char* topo_names_s[] = { "Vertex", "Curve", "Surface", "Volume" };
//  const char* topo_names_p[] = { "Vertices", "Curves", "Surfaces", "Volumes" };
  std::ostringstream buffer; // buffer output in an attempt to prevent lines from different processsors being mixed up.
  for (int t = 0; t < 4; ++t) {
    Range geom;
    int dim = t;
    const void* ptr = &dim;
    rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag, &ptr, 1, geom );
    CHECK_ERR(rval);
  
      // for each geometric entity of dimension 't'
    for (Range::const_iterator i = geom.begin(); i != geom.end(); ++i) {
      EntityHandle set = *i;
      int id;
      rval = moab.tag_get_data( id_tag, &set, 1, &id ); CHECK_ERR(rval);

      buffer.clear();

        // get entities contained in this set but not its children
      Range entities, tmp_entities, children, diff;
      rval = moab.get_entities_by_handle( set, entities ); CHECK_ERR(rval);
      rval = moab.get_child_meshsets( set, children ); CHECK_ERR(rval);
      for (Range::const_iterator j = children.begin(); j != children.end(); ++j) {
        tmp_entities.clear();
        rval = moab.get_entities_by_handle( *j, tmp_entities ); CHECK_ERR(rval);
        diff = subtract( entities, tmp_entities );
        entities.swap( diff );
      }
      
        // for each entity, check owning processors
      std::vector<char> status_flags( entities.size(), 0 );
      std::vector<int> shared_procs( entities.size(), 0 );
      rval = moab.tag_get_data( pstatus_tag, entities, &status_flags[0] );
      if (MB_TAG_NOT_FOUND == rval) {
        // keep default values of zero (not shared)
      }
      CHECK_ERR(rval);
      unsigned num_shared = 0, num_owned = 0;
      for (size_t j = 0; j < status_flags.size(); ++j) {
        num_shared += !!(status_flags[j] & PSTATUS_SHARED);
        num_owned += !(status_flags[j] & PSTATUS_NOT_OWNED);
      }
      
      if (!num_shared) {
        if (list_non_shared)
          buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
                 << "not shared" << std::endl;
      }
      else if (num_shared != entities.size()) {
        buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
               << "ERROR: " << num_shared << " of " << entities.size() 
               << " entities marked as 'shared'" << std::endl;
      }
      else if (num_owned && num_owned != entities.size()) {
        buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
               << "ERROR: " << num_owned << " of " << entities.size() 
               << " entities owned by this processor" << std::endl;
      }
      else {
        rval = moab.tag_get_data( sharedp_tag, entities, &shared_procs[0] );
        CHECK_ERR(rval);
        int proc = shared_procs[0];
        bool all_match = true;
        for (size_t j = 1; j < shared_procs.size(); ++j)
          if (shared_procs[j] != proc)
            all_match = false;
        if (!all_match) {
          buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
                 << "ERROR: processsor IDs do not match!" << std::endl;
        }
        else if (proc != -1) {
          buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
                 << "shared with processor " << proc;
          if (num_owned)
            buffer << " (owned by this processor)";
          buffer << std::endl;
        }
        else if (entities.empty()) {
          buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
                 << "ERROR: no entities!" << std::endl;
        }
        else {
          Range::const_iterator j = entities.begin();
          rval = moab.tag_get_data( sharedps_tag, &*j, 1, &ent_procs[0] );
          CHECK_ERR(rval);
          for (++j; j != entities.end(); ++j) {
            rval = moab.tag_get_data( sharedps_tag, &*j, 1, &tmp_ent_procs[0] );
            CHECK_ERR(rval);
            if (ent_procs != tmp_ent_procs) 
              all_match = false;
          }
          if (!all_match) {
            buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
                   << "ERROR: processsor IDs do not match!" << std::endl;
          }
          else {
            buffer << rank << ":\t" << topo_names_s[t] << " " << id << ":\t"
                   << "processors ";
            for (int k = 0; k < MAX_SHARING_PROCS; ++k)
              if (ent_procs[k] != -1)
                buffer << ent_procs[k] << ", ";
            if (num_owned)
              buffer << " (owned by this processor)";
            buffer << std::endl;
          }
        }
      }
    }
  }
  for (int i = 0; i < size; ++i) {
    MPI_Barrier( MPI_COMM_WORLD );
    if (i == rank) {
      std::cout << buffer.str();
      std::cout.flush();
    }
  }
}

void load_and_partition( Interface& moab, const char* filename, bool print )
{
  ErrorCode rval;
  
  rval = moab.load_file( filename, 0, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS" );

  if (print)
    print_partitioned_entities(moab);

  CHECK_ERR(rval);
}

void save_and_load_on_root( Interface& moab, const char* tmp_filename )
{
  ErrorCode rval;
  int procnum;
  MPI_Comm_rank( MPI_COMM_WORLD, &procnum );
  
  const char* opt = "PARALLEL=WRITE_PART";
  std::string str;
  if (WriteDebugLevel) {
    std::ostringstream s;
    s << opt << ";DEBUG_IO=" << WriteDebugLevel;
    str = s.str();
    opt = str.c_str();
  }
  rval = moab.write_file( tmp_filename, 0, opt );
  if (MB_SUCCESS != rval) {
    std::cerr << "Parallel write failed on processor " << procnum << std::endl;
    if (procnum == 0 && !KeepTmpFiles)
      remove( tmp_filename );
    CHECK_ERR(rval);
  }
  
  if (procnum == 0 && KeepTmpFiles)
    std::cout << "Wrote file: \"" << tmp_filename << "\"\n";

  // All created pcomm objects should be retrieved (with the pcomm tag) and
  // deleted at this time. Otherwise, the memory used by them will be leaked
  // as the pcomm tag will be deleted by moab.delete_mesh() below.
  std::vector<ParallelComm*> pc_list;
  ParallelComm::get_all_pcomm(&moab, pc_list);
  for (std::vector<ParallelComm*>::iterator vit = pc_list.begin();
       vit != pc_list.end(); ++vit)
    delete *vit;

  moab.delete_mesh();
  std::vector<Tag> tags;
  rval = moab.tag_get_tags( tags );
  CHECK_ERR(rval);
  for (size_t i = 0; i < tags.size(); ++i) {
    rval = moab.tag_delete( tags[i] );
    CHECK_ERR(rval);
  }
  
  if (procnum == 0) {
    rval = moab.load_file( tmp_filename );
    if (!KeepTmpFiles)
      remove( tmp_filename );
    CHECK_ERR(rval);
  }
}

void count_owned_entities( Interface& moab, int counts[MBENTITYSET] )
{
  ErrorCode rval;
  ParallelComm* pcomm = ParallelComm::get_pcomm( &moab, 0 );
  CHECK(0 != pcomm);
  std::fill( counts, counts+MBENTITYSET, 0u );
  
  for (EntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
    Range range;
    rval = moab.get_entities_by_type( 0, t, range );
    CHECK_ERR(rval);
    if (!range.empty())
      rval = pcomm->filter_pstatus(range, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    CHECK_ERR(rval);
    counts[t] = range.size();
  }
}

void check_identical_mesh( Interface& mb1, Interface& mb2 )
{
  ErrorCode rval;
  std::map<EntityHandle,EntityHandle> entmap;
  
    // match vertices by coordinate
  Range r1, r2;
  Range::iterator i1, i2;
  rval = mb1.get_entities_by_type( 0, MBVERTEX, r1 );
  CHECK_ERR(rval);
  rval = mb2.get_entities_by_type( 0, MBVERTEX, r2 );
  CHECK_ERR(rval);
  CHECK_EQUAL( r1.size(), r2.size() );
  for (i1 = r1.begin(); i1 != r1.end(); ++i1) {
    double coords1[3];
    rval = mb1.get_coords( &*i1, 1, coords1 );
    CHECK_ERR(rval);
    for (i2 = r2.begin(); i2 != r2.end(); ++i2) {
      double coords2[3];
      rval = mb2.get_coords( &*i2, 1, coords2 );
      CHECK_ERR(rval);
      coords2[0] -= coords1[0];
      coords2[1] -= coords1[1];
      coords2[2] -= coords1[2];
      double lensqr = coords2[0]*coords2[0] + coords2[1]*coords2[1] + coords2[2]*coords2[2];
      if (lensqr < 1e-12)
        break;
    }
    CHECK( i2 != r2.end() );
    entmap[*i2] = *i1;
    r2.erase( i2 );
  }
  
    // match element connectivity
  std::vector<EntityHandle> conn1, conn2;
  for (EntityType t = MBEDGE; t < MBENTITYSET; ++t) {
    r1.clear();
    rval = mb1.get_entities_by_type( 0, t, r1 );
    CHECK_ERR(rval);
    r2.clear();
    rval = mb2.get_entities_by_type( 0, t, r2 );
    CHECK_ERR(rval);
    CHECK_EQUAL( r1.size(), r2.size() );
    
    for (i1 = r1.begin(); i1 != r1.end(); ++i1) {
      conn1.clear();
      rval = mb1.get_connectivity( &*i1, 1, conn1 );
      CHECK_ERR(rval);
      for (i2 = r2.begin(); i2 != r2.end(); ++i2) {
        conn2.clear();
        rval = mb2.get_connectivity( &*i2, 1, conn2 );
        CHECK_ERR(rval);
        if (conn1.size() != conn2.size())
          continue;
        for (std::vector<EntityHandle>::iterator j = conn2.begin(); j != conn2.end(); ++j)
          *j = entmap[*j];
        if (conn1 == conn2)
          break;
      }

      CHECK( i2 != r2.end() );
      entmap[*i2] = *i1;
      r2.erase( i2 );
    }
  }
}

void test_write_elements()
{
  int proc_counts[MBENTITYSET], all_counts[MBENTITYSET], file_counts[MBENTITYSET];
  int err, rank;
  ErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  
    // load and partition a .cub file
  Core moab_instance;
  Interface& moab = moab_instance;
  load_and_partition( moab, InputFile, false );
  
    // count number of owned entities of each type and sum over all procs
  count_owned_entities( moab, proc_counts );
  std::fill( all_counts, all_counts+MBENTITYSET, 0u );
  err = MPI_Allreduce( proc_counts, all_counts, MBENTITYSET, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  CHECK(!err);
  
    // do parallel write and on root proc do serial read of written file
  save_and_load_on_root( moab, "test_write_elements.h5m" );
  if (rank == 0) {
    for (EntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
      rval = moab.get_number_entities_by_type( 0, t, file_counts[t] );
      CHECK_ERR(rval);
    }
  }
  
    // check that sum of owned entities equals number of entities from serial read
    
  err = MPI_Bcast( file_counts, MBENTITYSET, MPI_INT, 0, MPI_COMM_WORLD );
  CHECK(!err);
  
  bool all_equal = true;
  for (EntityType t = MBVERTEX; t < MBENTITYSET; ++t) 
    if (file_counts[t] != all_counts[t])
      all_equal = false;
    
  if (rank == 0 && !all_equal) {
    std::cerr << "Type\tPartnd\tWritten" << std::endl;
    for (EntityType t = MBVERTEX; t < MBENTITYSET; ++t) 
      std::cerr << CN::EntityTypeName(t) << '\t' << all_counts[t] << '\t' << file_counts[t] << std::endl;
  }
  
  CHECK(all_equal);
  
    // on root processor, do serial read of original .cub file and compare
  
  if (rank == 0) {
    Core moab2;
    rval = moab2.load_file( InputFile );
    CHECK_ERR(rval);
    check_identical_mesh( moab, moab2 );
  }
}

bool check_sets_sizes( Interface& mb1, EntityHandle set1,
                       Interface& mb2, EntityHandle set2 )
{
  ErrorCode rval;
  bool result = true;
  for (EntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
    int count1, count2;
    rval = mb1.get_number_entities_by_type( set1, t, count1 );
    CHECK_ERR(rval);
    rval = mb2.get_number_entities_by_type( set2, t, count2 );
    CHECK_ERR(rval);
    if (count1 != count2) {
      std::cerr << "Sets differ in number of " << CN::EntityTypeName(t)
                << " : " << count1 << " vs. " << count2 << std::endl;
      result = false;
    }
  }
  return result;
}

void test_write_shared_sets()
{
  int err, rank, size;
  ErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  err = MPI_Comm_size( MPI_COMM_WORLD, &size );
  CHECK(!err);
  
  Core moab_instance;
  Interface& moab = moab_instance;
  load_and_partition( moab, InputFile );
  save_and_load_on_root( moab, "test_write_shared_sets.h5m" );

  if (rank != 0)
    return;
  
  Core moab2_instance;
  Interface& moab2 = moab2_instance;
  rval = moab2.load_file( InputFile );
  CHECK_ERR(rval);
  
  Tag mattag1, mattag2;
  rval = moab.tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, mattag1 );
  CHECK_ERR(rval);
  rval = moab2.tag_get_handle( MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, mattag2 );
  CHECK_ERR(rval);
  
  Range matsets;
  rval = moab2.get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag2, 0, 1, matsets );
  CHECK_ERR(rval);
  for (Range::iterator i = matsets.begin(); i != matsets.end(); ++i) {
    int block_id;
    rval = moab2.tag_get_data( mattag2, &*i, 1, &block_id );
    CHECK_ERR(rval);
    
    Range tmpents;
    void* tagdata[] = {&block_id};
    rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag1, tagdata, 1, tmpents );
    if (tmpents.size() != 1) 
      std::cerr << tmpents.size() << " sets with material set id " << block_id << std::endl;
    CHECK_EQUAL( (int)tmpents.size(), 1 );
  
    CHECK( check_sets_sizes( moab2, *i, moab, tmpents.front() ) );
  }
}


void test_var_length_parallel()
{
  Range::const_iterator i;
  ErrorCode rval;
  Core moab;
  Interface &mb = moab;
  Range verts;
  Tag vartag;
  const char* filename = "var-len-para.h5m";
  const char* tagname = "ParVar";
  
  // If this tag doesn't exist, writer will fail
  Tag junk_tag;
  mb.tag_get_handle( PARALLEL_GID_TAG_NAME, 1, MB_TYPE_INTEGER, junk_tag, MB_TAG_DENSE|MB_TAG_EXCL );

  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  
  // Create N+1 vertices on each processor, where N is the rank 
  std::vector<double> coords( 3*rank+3, (double)rank );
  rval = mb.create_vertices( &coords[0], rank+1, verts );
  CHECK_ERR(rval);
  
  // Create a var-len tag
  rval = mb.tag_get_handle( tagname, 0, MB_TYPE_INTEGER, vartag, MB_TAG_DENSE|MB_TAG_VARLEN|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
  // Write data on each vertex:
  // { n, rank, rank+1, ..., rank+n-1 } where n >= 1
  std::vector<int> data;
  rval = MB_SUCCESS;
  for (i = verts.begin(); i != verts.end(); ++i) {
    EntityHandle h = *i;
    const int n = h % 7 + 1;
    data.resize( n+1 );
    data[0] = n;
    for (int j = 0; j < n; ++j)
      data[j+1] = rank + j;
    const int s = (n + 1);
    const void* ptrarr[] = { &data[0] };
    ErrorCode tmperr = mb.tag_set_by_ptr( vartag, &h, 1, ptrarr, &s );
    if (MB_SUCCESS != tmperr)
      rval = tmperr;
  }
  CHECK_ERR(rval);
  
  // Write file
  const char* opt = "PARALLEL=WRITE_PART";
  std::string str;
  if (WriteDebugLevel) {
    std::ostringstream s;
    s << opt << ";DEBUG_IO=" << WriteDebugLevel;
    str = s.str();
    opt = str.c_str();
  }
  rval = mb.write_file( filename, "MOAB", opt );
  CHECK_ERR(rval);
  
  // Read file.  We only reset and re-read the file on the
  // root processor.  All other processors keep the pre-write
  // mesh.  This allows many of the tests to be run on all
  // processors.  Running the tests on the pre-write mesh on
  // non-root processors allows us to verify that any problems
  // are due to the file API rather than some other bug.
  ErrorCode rval2 = rval = MB_SUCCESS;
  if (!rank) {
    moab.~Core();
    new (&moab) Core;
    rval = mb.load_mesh( filename );
    if (!KeepTmpFiles) 
      remove( filename );
    rval2 = mb.tag_get_handle( tagname, 0, MB_TYPE_INTEGER, vartag );
  }
  CHECK_ERR(rval);
  CHECK_ERR(rval2);
  
  // Check that tag is correct
  int tag_size;
  rval = mb.tag_get_length( vartag, tag_size );
  CHECK_EQUAL( MB_VARIABLE_DATA_LENGTH, rval );
  TagType storage;
  rval = mb.tag_get_type( vartag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  DataType type;
  rval = mb.tag_get_data_type( vartag, type);
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TYPE_INTEGER, type );
  
  // get vertices
  verts.clear();
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  // Check consistency of tag data on each vertex 
  // and count the number of vertices for each rank.
  std::vector<int> vtx_counts( numproc, 0 );
  for (i = verts.begin(); i != verts.end(); ++i) {
    EntityHandle h = *i;
    int size = -1;
    const void* ptrarr[1] = { 0 };
    rval = mb.tag_get_by_ptr( vartag, &h, 1, ptrarr, &size );
    CHECK_ERR( rval );
    const int* tag_data = reinterpret_cast<const int*>(ptrarr[0]);
    CHECK( size >= 2 );
    CHECK( NULL != tag_data );
    CHECK_EQUAL( size - 1, tag_data[0] );
    CHECK( tag_data[1] >= 0 && tag_data[1] < numproc );
    ++vtx_counts[tag_data[1]];
    for (int j = 1; j < size - 1; ++j)
      CHECK_EQUAL( tag_data[1] + j, tag_data[1 + j] );
  }
  
  // Check number of vertices for each rank
  for (int j = 0; j < numproc; ++j) {
    // Only root should have data for other processors.
    if (rank == 0 || rank == j) 
      CHECK_EQUAL( j + 1, vtx_counts[j] );
    else 
      CHECK_EQUAL( 0, vtx_counts[j] );
  }
}

// create row of cubes of mesh
void create_input_file( const char* file_name, 
                        int intervals, 
                        int num_cpu,
                        int blocks_per_cpu = 1,
                        const char* ijk_vert_tag_name = 0,
                        const char* ij_set_tag_name = 0,
                        const char* global_tag_name = 0,
                        const int* global_mesh_value = 0, 
                        const int* global_default_value = 0,
                        bool create_bcsets = false)
{
  Core moab;
  Interface& mb = moab;
  ErrorCode rval;
  
  Tag ijk_vert_tag = 0, ij_set_tag = 0, global_tag = 0;
  if (ijk_vert_tag_name) {
    rval = mb.tag_get_handle( ijk_vert_tag_name, 3, MB_TYPE_INTEGER, 
                          ijk_vert_tag, MB_TAG_EXCL|MB_TAG_DENSE );
    CHECK_ERR(rval);
  }
  if (ij_set_tag_name) {
    rval = mb.tag_get_handle( ij_set_tag_name, 2, MB_TYPE_INTEGER, ij_set_tag,
                          MB_TAG_SPARSE|MB_TAG_EXCL );
    CHECK_ERR(rval);
  }
  if (global_tag_name) {
    rval = mb.tag_get_handle( global_tag_name, 1, MB_TYPE_INTEGER, global_tag, 
                              MB_TAG_DENSE|MB_TAG_EXCL, global_default_value );
    CHECK_ERR(rval);
    if (global_mesh_value) {
      EntityHandle root = 0;
      rval = mb.tag_set_data( global_tag, &root, 1, global_mesh_value );
      CHECK_ERR(rval);
    }
  }
  
  const int num_blk = num_cpu * blocks_per_cpu;
  int iv = intervals+1, ii = num_blk*intervals+1;
  std::vector<EntityHandle> verts(iv*iv*ii);
  int idx = 0;
  for (int i = 0; i < ii; ++i) {
    for (int j = 0; j < iv; ++j) {
      int start = idx;
      for (int k = 0; k < iv; ++k) {
        const double coords[3] = {static_cast<double>(i), static_cast<double>(j), 
                                  static_cast<double>(k)};
        rval = mb.create_vertex( coords, verts[idx] );
        CHECK_ERR(rval);
        if (ijk_vert_tag) {
          int vals[] = {i,j,k};
          rval = mb.tag_set_data( ijk_vert_tag, &verts[idx], 1, vals );
          CHECK_ERR(rval);
        }
        ++idx; 
      }
      
      if (ij_set_tag) {
        EntityHandle set;
        rval = mb.create_meshset( MESHSET_SET, set );
        CHECK_ERR(rval);
        rval = mb.add_entities( set, &verts[start], idx - start );
        CHECK_ERR(rval);
        int vals[] = { i, j };
        rval = mb.tag_set_data( ij_set_tag, &set, 1, vals );
        CHECK_ERR(rval);
      }
    }
  }
  
  const int eb = intervals*intervals*intervals;
  std::vector<EntityHandle> elems(num_blk*eb);
  idx = 0;
  for (int c = 0; c < num_blk; ++c) {
    for (int i = c*intervals; i < (c+1)*intervals; ++i) {
      for (int j = 0; j < intervals; ++j) {
        for (int k = 0; k < intervals; ++k) {
          EntityHandle conn[8] = { verts[iv*(iv* i +      j    ) + k    ],
                                   verts[iv*(iv*(i + 1) + j    ) + k    ],
                                   verts[iv*(iv*(i + 1) + j + 1) + k    ],
                                   verts[iv*(iv* i      + j + 1) + k    ],
                                   verts[iv*(iv* i      + j    ) + k + 1],
                                   verts[iv*(iv*(i + 1) + j    ) + k + 1],
                                   verts[iv*(iv*(i + 1) + j + 1) + k + 1],
                                   verts[iv*(iv* i      + j + 1) + k + 1] };

          rval = mb.create_element( MBHEX, conn, 8, elems[idx++] );
          CHECK_ERR(rval);
        }
      }
    }
  }
  
  Tag part_tag;
  rval = mb.tag_get_handle( PARTITION_TAG, 1, MB_TYPE_INTEGER, part_tag, MB_TAG_SPARSE|MB_TAG_EXCL );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> parts(num_cpu);
  for (int i = 0; i < num_cpu; ++i) {
    rval = mb.create_meshset( MESHSET_SET, parts[i] );
    CHECK_ERR(rval);
    for (int j = 0; j < blocks_per_cpu; ++j) {
      rval = mb.add_entities( parts[i], &elems[(num_cpu*j + i)*eb], eb );
      CHECK_ERR(rval);
    }
    rval = mb.tag_set_data( part_tag, &parts[i], 1, &i );
    CHECK_ERR(rval);
  }

  if (create_bcsets) {
      // neumann set
    Range skin_ents;
    rval = Skinner(&mb).find_skin(0, &elems[0], elems.size(), false, skin_ents);
    CHECK_ERR(rval);
    EntityHandle bcset;
    rval = mb.create_meshset( MESHSET_SET, bcset);
    CHECK_ERR(rval);
    Tag bcset_tag;
    rval = mb.tag_get_handle(NEUMANN_SET_TAG_NAME, 1, MB_TYPE_INTEGER, bcset_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
    CHECK_ERR(rval);
    int dum = 100;
    rval = mb.tag_set_data(bcset_tag, &bcset, 1, &dum);
    CHECK_ERR(rval);
    rval = mb.add_entities(bcset, skin_ents);
    CHECK_ERR(rval);

      // dirichlet set
    rval = mb.create_meshset( MESHSET_SET, bcset);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle(DIRICHLET_SET_TAG_NAME, 1, MB_TYPE_INTEGER, bcset_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
    CHECK_ERR(rval);
    dum = 200;
    rval = mb.tag_set_data(bcset_tag, &bcset, 1, &dum);
    CHECK_ERR(rval);
    Range nodes;
    rval = mb.get_adjacencies(skin_ents, 0, false, nodes, Interface::UNION);
    CHECK_ERR(rval);
    rval = mb.add_entities(bcset, nodes);
    CHECK_ERR(rval);

      // material set
    rval = mb.create_meshset( MESHSET_SET, bcset);
    CHECK_ERR(rval);
    rval = mb.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, bcset_tag, MB_TAG_SPARSE|MB_TAG_CREAT);
    CHECK_ERR(rval);
    dum = 300;
    rval = mb.tag_set_data(bcset_tag, &bcset, 1, &dum);
    CHECK_ERR(rval);
    rval = mb.add_entities(bcset, &elems[0], elems.size());
    CHECK_ERR(rval);
  }
  
  rval = mb.write_file( file_name, "MOAB" );
  CHECK_ERR(rval);
}

void test_read_elements_common( bool by_rank, int intervals, bool /* print_time */,
                                const char* extra_opts )
{
  const char *file_name = by_rank ? "test_read_rank.h5m" : "test_read.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;

    // if root processor, create hdf5 file for use in testing
  if (0 == rank) 
    create_input_file( file_name, intervals, numproc );
  MPI_Barrier(MPI_COMM_WORLD); // make sure root has completed writing the file
  
    // do parallel read unless only one processor
  std::string opt = get_read_options( by_rank, extra_opts );
  rval = mb.load_file( file_name, 0, opt.c_str() );
    
  MPI_Barrier(MPI_COMM_WORLD); // make sure all procs complete before removing file
  if (0 == rank && !KeepTmpFiles) remove( file_name );
  CHECK_ERR(rval);
      
  
  Tag part_tag;
  rval = mb.tag_get_handle( PARTITION_TAG, 1, MB_TYPE_INTEGER, part_tag );
  CHECK_ERR(rval);
  
  Range parts;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &part_tag, 0, 1, parts );
  CHECK_ERR(rval);
  CHECK_EQUAL( 1, (int)parts.size() );
  EntityHandle part = parts.front();
  int id;
  rval = mb.tag_get_data( part_tag, &part, 1, &id );
  CHECK_ERR(rval);
  if (by_rank) {
    CHECK_EQUAL( rank, id );
  }
  
    // check that all of the elements in the mesh are in the part
  int npart, nall;
  rval = mb.get_number_entities_by_dimension( part, 3, npart );
  CHECK_ERR(rval);
  rval = mb.get_number_entities_by_dimension( 0, 3, nall );
  CHECK_ERR(rval);
  CHECK_EQUAL( npart, nall );
  
    // check that we have the correct vertices
  const double x_min = intervals*rank;
  const double x_max = intervals*(rank+1);
  Range verts;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  std::vector<double> coords(verts.size());
  rval = mb.get_coords( verts, &coords[0], 0, 0 );
  CHECK_ERR(rval);
  const double act_x_min = *std::min_element( coords.begin(), coords.end() );
  const double act_x_max = *std::max_element( coords.begin(), coords.end() );
  CHECK_REAL_EQUAL( x_min, act_x_min, DBL_EPSILON );
  CHECK_REAL_EQUAL( x_max, act_x_max, DBL_EPSILON );
}

void test_read_time()
{
  const char file_name[] = "read_time.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  ErrorCode rval;

    // if root processor, create hdf5 file for use in testing
  if (0 == rank) 
    create_input_file( file_name, ReadIntervals, numproc, ReadBlocks );
  MPI_Barrier( MPI_COMM_WORLD );
  
  // CPU Time for true paralle, wall time for true paralle,
  //   CPU time for read and delete, wall time for read and delete
  double times[6];
  clock_t tmp_t;
  
    // Time true parallel read
  Core moab;
  Interface &mb = moab;
  times[0] = MPI_Wtime();
  tmp_t    = clock();
  std::string opt = get_read_options( true, READ_PART );
  rval = mb.load_file( file_name, 0, opt.c_str() );
  CHECK_ERR(rval);
  times[0] = MPI_Wtime() - times[0];
  times[1] = double(clock() - tmp_t) / CLOCKS_PER_SEC;
  mb.delete_mesh();
    
    // Time read and delete
  Core moab2;
  Interface& mb2 = moab2;
  std::string opt2 = get_read_options( true, READ_DELETE );
  times[2] = MPI_Wtime();
  tmp_t    = clock();
  rval = mb2.load_file( file_name, 0, opt2.c_str() );
  CHECK_ERR(rval);
  times[2] = MPI_Wtime() - times[2];
  times[3] = double(clock() - tmp_t) / CLOCKS_PER_SEC;
  mb2.delete_mesh();
    
    // Time broadcast and delete
  Core moab3;
  Interface& mb3 = moab3;
  std::string opt3 = get_read_options( true, BCAST_DELETE );
  times[4] = MPI_Wtime();
  tmp_t    = clock();
  rval = mb3.load_file( file_name, 0, opt3.c_str() );
  CHECK_ERR(rval);
  times[4] = MPI_Wtime() - times[4];
  times[5] = double(clock() - tmp_t) / CLOCKS_PER_SEC;
  mb3.delete_mesh();
  
  double max_times[6] = {0,0,0,0,0,0}, sum_times[6] = {0,0,0,0,0,0}; 
  MPI_Reduce( &times, &max_times, 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  MPI_Reduce( &times, &sum_times, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Barrier( MPI_COMM_WORLD );
  if (0 == rank) {
    printf( "%12s  %12s  %12s  %12s\n", "", "READ_PART", "READ_DELETE", "BCAST_DELETE" );
    printf( "%12s  %12g  %12g  %12g\n", "Max Wall", max_times[0], max_times[2], max_times[4] );
    printf( "%12s  %12g  %12g  %12g\n", "Total Wall", sum_times[0], sum_times[2], sum_times[4] );
    printf( "%12s  %12g  %12g  %12g\n", "Max CPU", max_times[1], max_times[3], max_times[5] );
    printf( "%12s  %12g  %12g  %12g\n", "Total CPU", sum_times[1], sum_times[3], sum_times[5] );
  }
  
  MPI_Barrier( MPI_COMM_WORLD );
  if (0 == rank && !KeepTmpFiles) remove( file_name );
}

void test_read_tags()
{
  const char tag_name[] = "test_tag_xx";
  const char file_name[] = "test_read_tags.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;

    // if root processor, create hdf5 file for use in testing
  if (0 == rank) 
    create_input_file( file_name, DefaultReadIntervals, numproc, 1, tag_name );
  MPI_Barrier(MPI_COMM_WORLD); // make sure root has completed writing the file
  
    // do parallel read unless only one processor
  std::string opt = get_read_options( );
  rval = mb.load_file( file_name, 0, opt.c_str() );
  MPI_Barrier(MPI_COMM_WORLD); // make sure all procs complete before removing file
  if (0 == rank && !KeepTmpFiles) remove( file_name );
  CHECK_ERR(rval);
      
  Tag tag;
  rval = mb.tag_get_handle( tag_name, 3, MB_TYPE_INTEGER, tag );
  CHECK_ERR(rval);
  
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  
  Range verts, tagged;
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type_and_tag( 0, MBVERTEX, &tag, 0, 1, tagged );
  CHECK_ERR(rval);
  CHECK_EQUAL( verts, tagged );
  
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    double coords[3];
    rval = mb.get_coords( &*i, 1, coords );
    CHECK_ERR(rval);
    int ijk[3];
    rval = mb.tag_get_data( tag, &*i, 1, ijk );
    CHECK_ERR(rval);
    
    CHECK( ijk[0] >= DefaultReadIntervals * rank );
    CHECK( ijk[0] <= DefaultReadIntervals * (rank+1) );
    CHECK( ijk[1] >= 0 );
    CHECK( ijk[1] <= DefaultReadIntervals );
    CHECK( ijk[2] >= 0 );
    CHECK( ijk[2] <= DefaultReadIntervals );

    CHECK_REAL_EQUAL( coords[0], (double)ijk[0], 1e-100 );
    CHECK_REAL_EQUAL( coords[1], (double)ijk[1], 1e-100 );
    CHECK_REAL_EQUAL( coords[2], (double)ijk[2], 1e-100 );
  }
}

void test_read_global_tags()
{
  const char tag_name[] = "test_tag_g";
  const char file_name[] = "test_read_global_tags.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;
  const int def_val = 0xdeadcad;
  const int global_val = -11;

    // if root processor, create hdf5 file for use in testing
  if (0 == rank) 
    create_input_file( file_name, 1, numproc, 1, 0, 0, tag_name, &global_val, &def_val );
  MPI_Barrier(MPI_COMM_WORLD); // make sure root has completed writing the file
  
    // do parallel read unless only one processor
  std::string opt = get_read_options( );
  rval = mb.load_file( file_name, 0, opt.c_str() );
  MPI_Barrier(MPI_COMM_WORLD); // make sure all procs complete before removing file
  if (0 == rank && !KeepTmpFiles) remove( file_name );
  CHECK_ERR(rval);
      
  Tag tag;
  rval = mb.tag_get_handle( tag_name, 1, MB_TYPE_INTEGER, tag );
  CHECK_ERR(rval);
  
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_DENSE, storage );
  
  int mesh_def_val, mesh_gbl_val;
  rval = mb.tag_get_default_value( tag, &mesh_def_val );
  CHECK_ERR(rval);
  CHECK_EQUAL( def_val, mesh_def_val );
  EntityHandle root = 0;
  rval = mb.tag_get_data( tag, &root, 1, &mesh_gbl_val );
  CHECK_ERR(rval);
  CHECK_EQUAL( global_val, mesh_gbl_val );
}

void test_read_sets_common( const char* extra_opts )
{
  const char tag_name[] = "test_tag_s";
  const char file_name[] = "test_read_sets.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;

    // if root processor, create hdf5 file for use in testing
  if (0 == rank) 
    create_input_file( file_name, DefaultReadIntervals, numproc, 1, 0, tag_name );
  MPI_Barrier(MPI_COMM_WORLD); // make sure root has completed writing the file
  
    // do parallel read unless only one processor
  std::string opt = get_read_options( extra_opts );
  rval = mb.load_file( file_name, 0, opt.c_str() );
  MPI_Barrier(MPI_COMM_WORLD); // make sure all procs complete before removing file
  if (0 == rank && !KeepTmpFiles) remove( file_name );
  CHECK_ERR(rval);
      
  Tag tag;
  rval = mb.tag_get_handle( tag_name, 2, MB_TYPE_INTEGER, tag );
  CHECK_ERR(rval);
  
  TagType storage;
  rval = mb.tag_get_type( tag, storage );
  CHECK_ERR(rval);
  CHECK_EQUAL( MB_TAG_SPARSE, storage );
  
  const int iv = DefaultReadIntervals + 1;
  Range sets;
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, 0, 1, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( (iv*iv), (int)sets.size() );
  
  for (Range::iterator i = sets.begin(); i != sets.end(); ++i) {
    int ij[2];
    rval = mb.tag_get_data( tag, &*i, 1, &ij );
    CHECK_ERR(rval);
    
    CHECK( ij[0] >= DefaultReadIntervals * rank );
    CHECK( ij[0] <= DefaultReadIntervals * (rank+1) );
    CHECK( ij[1] >= 0 );
    CHECK( ij[1] <= DefaultReadIntervals );
    
    Range contents;
    rval = mb.get_entities_by_handle( *i, contents );
    CHECK_ERR(rval);
    CHECK(contents.all_of_type(MBVERTEX));
    CHECK_EQUAL( iv, (int)contents.size() );
    
    for (Range::iterator v = contents.begin(); v != contents.end(); ++v) {
      double coords[3];
      rval = mb.get_coords( &*v, 1, coords );
      CHECK_ERR(rval);
      CHECK_REAL_EQUAL( coords[0], (double)ij[0], 1e-100 );
      CHECK_REAL_EQUAL( coords[1], (double)ij[1], 1e-100 );
    }
  }
}

void test_read_bc_sets()
{
  //const char tag_name[] = "test_tag_s";
  const char file_name[] = "test_read_sets.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;

    // if root processor, create hdf5 file for use in testing
  if (0 == rank) 
    create_input_file( file_name, DefaultReadIntervals, numproc, 1, NULL, NULL, NULL, NULL, NULL, true );
  MPI_Barrier(MPI_COMM_WORLD); // make sure root has completed writing the file
  
    // do parallel read unless only one processor
  std::string opt = get_read_options();
  rval = mb.load_file( file_name, 0, opt.c_str() );
  MPI_Barrier(MPI_COMM_WORLD); // make sure all procs complete before removing file
  if (0 == rank && !KeepTmpFiles) remove( file_name );
  CHECK_ERR(rval);
      
  Tag tag;
  int num_ents[3], global_num_ents[3] = {0, 0, 0};
  Range sets, contents;
  const char *names[] = {NEUMANN_SET_TAG_NAME, DIRICHLET_SET_TAG_NAME, MATERIAL_SET_TAG_NAME};
  int vints = DefaultReadIntervals+1;
  int expected_num_ents[] = {(numproc*4+2)*DefaultReadIntervals*DefaultReadIntervals, 
                             ((numproc*4+2)*(vints-2)*(vints-2) + 12*numproc*(vints-2) + 8*numproc),
                             numproc*DefaultReadIntervals*DefaultReadIntervals*DefaultReadIntervals};

  for (int i = 0; i < 3; i++) {
    rval = mb.tag_get_handle(names[i], 1, MB_TYPE_INTEGER, tag );
    CHECK_ERR(rval);
    rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, 0, 1, sets );
    CHECK_ERR(rval);
    CHECK_EQUAL(1, (int)sets.size());
    rval = mb.get_entities_by_handle(*sets.begin(), contents, true);
    CHECK_ERR(rval);
    num_ents[i] = contents.size();
    sets.clear();
    contents.clear();
  }

  MPI_Reduce(num_ents, global_num_ents, 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
  if (0 == rank) {
//    std::cout << "Global:" << global_num_ents[0] << " " << global_num_ents[1] << " " << global_num_ents[2] << " " << std::endl;
//    std::cout << "Expected:" << expected_num_ents[0] << " " << expected_num_ents[1] << " " << expected_num_ents[2] << " " << std::endl;
    
    for (int i = 0; i < 3; i++) {
      CHECK_EQUAL(global_num_ents[i], expected_num_ents[i]);
    }
  }
}


void test_write_different_element_types()
{
  const char file_name[] = "test_write_different_element_types.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;
  
  const EntityType topos[] =  
   { MBEDGE, MBEDGE, MBQUAD, MBTRI, MBQUAD, MBTRI, MBTET, MBHEX, MBPRISM, MBPYRAMID, MBHEX, MBHEX };
  const int verts[] = 
   { 2,      3,      4,      3,     9,      6,     4,     8,     6,       5,         20,    27    };
  const int ntypes = sizeof(topos)/sizeof(topos[0]);
  
  const EntityType mtype = topos[rank%ntypes];
  const int nvert = verts[rank%ntypes];
  std::vector<EntityHandle> conn(nvert);
  for (int i = 0; i < nvert; ++i) {
    const double coords[] = { static_cast<double>(rank), static_cast<double>(i), 0 };
    rval = mb.create_vertex( coords, conn[i] );
    CHECK_ERR(rval);
  }
  EntityHandle handle;
  rval = mb.create_element( mtype, &conn[0], nvert, handle );
  CHECK_ERR(rval);
  
  save_and_load_on_root( mb, file_name );
  
  if (rank)
    return;
  
  for (int i = 0; i < ntypes; ++i) {
    const int num_exp = numproc / ntypes + (i < (numproc % ntypes) ? 1 : 0);
    Range ents;
    rval = mb.get_entities_by_type( 0, topos[i], ents );
    CHECK_ERR(rval);
    int num = 0;
    for (Range::iterator e = ents.begin(); e != ents.end(); ++e) {
      const EntityHandle* junk;
      int len;
      rval = mb.get_connectivity( *e, junk, len, false );
      if (len == verts[i])
        ++num;
    }
    
    CHECK_EQUAL( num_exp, num );
  }
}


Tag get_tag( Interface& mb, int rank, bool create )
{
  DataType type = (DataType)(rank % (MB_MAX_DATA_TYPE+1));
  TagType storage = (type == MB_TYPE_BIT) ? MB_TAG_BIT :
                    (rank % 2)            ? MB_TAG_DENSE :
                                            MB_TAG_SPARSE;
  int len = rank % 3 + 1;
  TagType cbit = create ? MB_TAG_EXCL : (TagType)0;
  TagType vbit = rank% 4 == 1 && storage != MB_TAG_BIT ? MB_TAG_VARLEN : (TagType)0;
  std::ostringstream name;
  name << "TestTag" << rank;
  const void* defval = 0;
  const int defint[] = { static_cast<int>(rank), static_cast<int>(rank/2), 
                         static_cast<int>(rank+1), static_cast<int>(rank-1) };
  const double defreal[] = { 0.1*rank, 1.0/rank, 
                             static_cast<double>(-rank), static_cast<double>(rank) };
  const int defhandle[] = { 0, 0, 0, 0 };
  const unsigned char defbit = 0x1;
  const char defopq[] = "Jason";
  if (rank % 4 < 2) {
    switch (type) {
      case MB_TYPE_BIT: defval = &defbit; break;
      case MB_TYPE_INTEGER: defval = defint; break;
      case MB_TYPE_DOUBLE: defval = defreal; break;
      case MB_TYPE_HANDLE: defval = defhandle; break;
      case MB_TYPE_OPAQUE: defval = defopq; break;
    }
  }
  
  Tag result;
  ErrorCode rval = mb.tag_get_handle( name.str().c_str(),
                                      len, type,
                                      result,
                                      storage|cbit|vbit,
                                      defval );
  CHECK_ERR(rval);
  return result;
}
              

void test_write_different_tags()
{
  const char file_name[] = "test_write_different_tags.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  
  const int min_tags = 8;
  get_tag( mb, rank, true );
  for (int idx = rank+numproc; idx < min_tags; idx+=numproc)
    get_tag( mb, idx, true );
  
  save_and_load_on_root( mb, file_name );
  
  if (0 == rank) {
    for (int i = 0; i < std::max(min_tags, numproc); ++i)
      get_tag( mb, i, false );
  }
}


void test_write_polygons()
{
  const char file_name[] = "test_write_polygons.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;
  
    // create a polygon on each process
  const double r = 0.70710678118654757;
  const double points[8][3] = { { 1, 0, static_cast<double>(rank) },
                                { static_cast<double>(r), static_cast<double>(r), static_cast<double>(rank) },
                                { 0, 1, static_cast<double>(rank) },
                                {static_cast<double>(-r), static_cast<double>(r), static_cast<double>(rank) },
                                {-1, 0, static_cast<double>(rank) },
                                {static_cast<double>(-r),static_cast<double>(-r), static_cast<double>(rank) },
                                { 0,-1, static_cast<double>(rank) },
                                { static_cast<double>(r),static_cast<double>(-r), static_cast<double>(rank) } };
  const int nvtx = rank % 4 + 5;
  std::vector<EntityHandle> conn(nvtx);
  for (int i = 0; i < nvtx; ++i) {
    rval = mb.create_vertex( points[i], conn[i] );
    CHECK_ERR(rval);
  }
  
  EntityHandle h;
  rval = mb.create_element( MBPOLYGON, &conn[0], nvtx, h );
  CHECK_ERR(rval);
  
  save_and_load_on_root( mb, file_name );

  if (rank != 0)
    return;
  
    // check results on root process
    
    // determine which polygon was created by which proc by 
    // looking at the z-coordinate of the vertices
  Range range;
  rval = mb.get_entities_by_type( 0, MBPOLYGON, range );
  std::vector<EntityHandle> poly( numproc, 0 );
  CHECK_EQUAL( numproc, (int)range.size() );
  for (Range::iterator it = range.begin(); it != range.end(); ++it) {
    const EntityHandle* conn_arr;
    int len;
    rval = mb.get_connectivity( *it, conn_arr, len );
    CHECK_ERR(rval);
    double coords[3];
    rval = mb.get_coords( conn_arr, 1, coords );
    CHECK_ERR(rval);
    int proc = (int)(coords[2]);
    CHECK_EQUAL( (EntityHandle)0, poly[proc] );
    poly[proc] = *it;
  }
  
    // check that each poly has the expected number of vertices
  for (int i = 0; i < numproc; ++i) {
    const EntityHandle* conn_arr;
    int len;
    rval = mb.get_connectivity( poly[i], conn_arr, len );
    CHECK_ERR(rval);
    CHECK_EQUAL( i % 4 + 5, len );
  }
}


// Some processes have no mesh to write
void test_write_unbalanced()
{
  const char file_name[] = "test_write_unbalanced.h5m";
  int numproc, rank;
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  Core moab;
  Interface &mb = moab;
  ErrorCode rval;
  Tag idtag;
  
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, idtag, MB_TAG_CREAT );
  CHECK_ERR(rval);
  
    // create a shared set
  const int two = 2;
  Range entities, sets;

  EntityHandle set;
  rval = mb.create_meshset( MESHSET_SET, set );
  CHECK_ERR(rval);
  rval = mb.tag_set_data( idtag, &set, 1, &two );
  CHECK_ERR(rval);
  sets.insert( set );

    // create a quad on every odd processor
  if (rank % 2) {
    const double coords[4][3] = { { static_cast<double>(rank),   0, 0 },
                                  { static_cast<double>(rank+2), 0, 0 },
                                  { static_cast<double>(rank+2), 2, 0 },
                                  { static_cast<double>(rank),   2, 0 } };
    EntityHandle conn[4], quad;
    for (int i = 0; i < 4; ++i)
      mb.create_vertex( coords[i], conn[i] );
    mb.create_element( MBQUAD, conn, 4, quad );
    
    const int ids[4] = { rank, rank+2, rank+3, rank+1 };
    rval = mb.tag_set_data( idtag, conn, 4, ids );
    CHECK_ERR(rval);
    
    rval = mb.add_entities( set, &quad, 1 );
    CHECK_ERR(rval);
    
    entities.insert(quad);
  }
  
    // set up sharing data
  ParallelComm* pcomm = ParallelComm::get_pcomm( &mb, 0 );
  if (0 == pcomm)
    pcomm = new ParallelComm( &mb, MPI_COMM_WORLD );
  rval = pcomm->resolve_shared_ents( 0, entities, 2, 0, NULL, &idtag );
  CHECK_ERR(rval);
  rval = pcomm->resolve_shared_sets( sets, idtag );
  CHECK_ERR(rval);
  
    // do parallel write and serial load
  save_and_load_on_root( mb, file_name );

  if (rank != 0)
    return;

    // check that we got what we expected
  Range quads, verts;
  rval = mb.get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  rval = mb.get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  const size_t nquads = numproc/2;
  const size_t nverts = nquads ? 2+2*nquads : 0;
  CHECK_EQUAL( nquads, quads.size() );
  CHECK_EQUAL( nverts, verts.size() );
  
  rval = mb.tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, idtag );
  CHECK_ERR(rval);
  sets.clear();
  const void* vals[] = { &two };
  rval = mb.get_entities_by_type_and_tag( 0, MBENTITYSET, &idtag, vals, 1, sets );
  CHECK_ERR(rval);
  CHECK_EQUAL( (size_t)1, sets.size() );
  
  entities.clear();
  rval = mb.get_entities_by_handle( sets.front(), entities );
  CHECK_ERR(rval);
  CHECK_EQUAL( nquads, entities.size() );
  CHECK_EQUAL( quads, entities );
}

// this test will load a file that has 2 partitions (mix.h5m)
// one partition with triangles, and one with triangles and quads
// a dense tag created on this should be dense in the file
//
void test_write_dense_tags()
{
  int err, rank, size;
  ErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  err = MPI_Comm_size( MPI_COMM_WORLD, &size );
  CHECK(!err);

  Core moab;
  rval =  moab.load_file( InputMix, 0,
      "PARALLEL=READ_PART;"
      "PARTITION=PARALLEL_PARTITION;"
      "PARALLEL_RESOLVE_SHARED_ENTS" );
  CHECK_ERR(rval);

  // create an dense tag, on all elements, then write the output in parallel
  // load the file again, and test the type of the tag
  Range elems;
  moab.get_entities_by_dimension(0, 2, elems);

  const char * tagname = "element_tag";
  const double defVal = 0.;
  Tag fieldTag;
  rval = moab.tag_get_handle(tagname, 1, MB_TYPE_DOUBLE, fieldTag, MB_TAG_DENSE|MB_TAG_CREAT, &defVal);
  CHECK_ERR(rval);


  int numElems = (int)elems.size();

  for(int i=0; i<numElems; i++){
    EntityHandle elem = elems[i];
    double xyz[3];
    moab.get_coords(&elem, 1, xyz);
    // some value
    double fieldValue =  sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]+xyz[2]*xyz[2]);
    moab.tag_set_data(fieldTag, &elem, 1, &fieldValue);
  }

  // write the file in parallel
  rval = moab.write_file("newfile.h5m", 0 , "PARALLEL=WRITE_PART");
  CHECK_ERR(rval);

  // now read the new file, in a new instance, and test the tag type

  Core moab2;
  rval = moab2.load_file( "newfile.h5m",0,
      "PARALLEL=READ_PART;"
      "PARTITION=PARALLEL_PARTITION;"
      "PARALLEL_RESOLVE_SHARED_ENTS" );
  CHECK_ERR(rval);

  // find the element tag
  Tag found_tag;
  rval = moab2.tag_get_handle(tagname, 1, MB_TYPE_DOUBLE, found_tag);
  CHECK_ERR(rval);
  TagType tagt;
  rval = moab2.tag_get_type(found_tag, tagt);
  CHECK_ERR(rval);
  CHECK(tagt == MB_TAG_DENSE);

}
// this test will load a file that has 2 partitions (oneside.h5m)
// and one side set, that is adjacent to one part only
//
void test_read_non_adjs_side()
{
  int err, rank, size;
  ErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  err = MPI_Comm_size( MPI_COMM_WORLD, &size );
  CHECK(!err);

  Core moab;
  rval =  moab.load_file( InputOneSide, 0,
      "PARALLEL=READ_PART;"
      "PARTITION=PARALLEL_PARTITION;"
      "PARALLEL_RESOLVE_SHARED_ENTS" );
  CHECK_ERR(rval);

  return;

}
