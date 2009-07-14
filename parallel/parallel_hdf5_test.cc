#include "TestUtil.hpp"

#include "MBCore.hpp"
#include "MBParallelComm.hpp"
#include "MBTagConventions.hpp"
#include "MBCN.hpp"
#include "MBParallelConventions.h"

#include <iostream>
#include <sstream>
#include <mpi.h>
#include <unistd.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#ifdef SRCDIR
const char* InputFile = STRINGIFY(SRCDIR) "/ptest.cub";
#else
const char* InputFile = "ptest.cub";
#endif

void load_and_partition( MBInterface& moab, const char* filename, bool print_debug = false );

void save_and_load_on_root( MBInterface& moab, const char* tmp_filename );

void check_identical_mesh( MBInterface& moab1, MBInterface& moab2 );

void test_write_elements();
void test_write_shared_sets();
void test_var_length_parallel();

bool KeepTmpFiles = false;
bool PauseOnStart = false;

int main( int argc, char* argv[] )
{
  int err = MPI_Init( &argc, &argv );
  CHECK(!err);

  for (int i = 1; i < argc; ++i) {
    if (!strcmp( argv[i], "-k"))
      KeepTmpFiles = true;
    else if (!strcmp( argv[i], "-p"))
      PauseOnStart = true;
    else {
      std::cerr << "Usage: " << argv[0] << " [-k] [-p]" << std::endl;
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
  result += RUN_TEST( test_write_elements );
  result += RUN_TEST( test_write_shared_sets );
  result += RUN_TEST( test_var_length_parallel );
  
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
void print_partitioned_entities( MBInterface& moab, bool list_non_shared = false )
{
  MBErrorCode rval;
  int size, rank;
  std::vector<int> ent_procs(MAX_SHARING_PROCS), tmp_ent_procs(MAX_SHARING_PROCS);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
    // expect shared entities to correspond to geometric sets
    
    // get tags for parallel data
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  const char* ptag_names[] = { PARALLEL_SHARED_PROC_TAG_NAME,
                               PARALLEL_SHARED_PROCS_TAG_NAME,
                               PARALLEL_SHARED_HANDLE_TAG_NAME,
                               PARALLEL_SHARED_HANDLES_TAG_NAME,
                               PARALLEL_STATUS_TAG_NAME };
  MBTag* tag_ptrs[] = { &sharedp_tag, &sharedps_tag, &sharedh_tag, &sharedhs_tag, &pstatus_tag };
  const int ntags = sizeof(ptag_names)/sizeof(ptag_names[0]);
  for (int i = 0; i < ntags; ++i) {
    rval = moab.tag_get_handle( ptag_names[i], *tag_ptrs[i] ); CHECK_ERR(rval);
  }
  
    // for each geometric entity, check which processor we are sharing
    // entities with
  MBTag geom_tag, id_tag;
  rval = moab.tag_get_handle( GEOM_DIMENSION_TAG_NAME, geom_tag ); CHECK_ERR(rval);
  rval = moab.tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag ); CHECK_ERR(rval);
  const char* topo_names_s[] = { "Vertex", "Curve", "Surface", "Volume" };
//  const char* topo_names_p[] = { "Vertices", "Curves", "Surfaces", "Volumes" };
  std::ostringstream buffer; // buffer output in an attempt to prevent lines from different processsors being mixed up.
  for (int t = 0; t < 4; ++t) {
    MBRange geom;
    int dim = t;
    const void* ptr = &dim;
    rval = moab.get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag, &ptr, 1, geom );
    CHECK_ERR(rval);
  
      // for each geometric entity of dimension 't'
    for (MBRange::const_iterator i = geom.begin(); i != geom.end(); ++i) {
      MBEntityHandle set = *i;
      int id;
      rval = moab.tag_get_data( id_tag, &set, 1, &id ); CHECK_ERR(rval);

      buffer.clear();

        // get entities contained in this set but not its children
      MBRange entities, tmp_entities, children, diff;
      rval = moab.get_entities_by_handle( set, entities ); CHECK_ERR(rval);
      rval = moab.get_child_meshsets( set, children ); CHECK_ERR(rval);
      for (MBRange::const_iterator j = children.begin(); j != children.end(); ++j) {
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
          MBRange::const_iterator j = entities.begin();
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
            for (int j = 0; j < MAX_SHARING_PROCS; ++j)
              if (ent_procs[j] != -1)
                buffer << ent_procs[j] << ", ";
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

void load_and_partition( MBInterface& moab, const char* filename, bool print )
{
  MBErrorCode rval;
  MBEntityHandle set;
  
  rval = moab.load_file( filename, set, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS" );

  if (print)
    print_partitioned_entities(moab);

  CHECK_ERR(rval);
}

void save_and_load_on_root( MBInterface& moab, const char* tmp_filename )
{
  MBErrorCode rval;
  MBEntityHandle set;
  int procnum;
  MPI_Comm_rank( MPI_COMM_WORLD, &procnum );
  
  rval = moab.write_file( tmp_filename, 0, "PARALLEL=FORMAT" );
  if (MB_SUCCESS != rval) {
    std::cerr << "Parallel write filed on processor " << procnum << std::endl;
    if (procnum == 0 && !KeepTmpFiles)
      remove( tmp_filename );
    CHECK_ERR(rval);
  }
  
  if (procnum == 0 && KeepTmpFiles)
    std::cout << "Wrote file: \"" << tmp_filename << "\"\n";
  
  moab.delete_mesh();
  if (procnum == 0) {
    rval = moab.load_file( tmp_filename, set );
    if (!KeepTmpFiles)
      remove( tmp_filename );
    CHECK_ERR(rval);
  }
}

void count_owned_entities( MBInterface& moab, int counts[MBENTITYSET] )
{
  MBErrorCode rval;
  MBParallelComm* pcomm = MBParallelComm::get_pcomm( &moab, 0 );
  CHECK(0 != pcomm);
  std::fill( counts, counts+MBENTITYSET, 0u );
  
  for (MBEntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
    MBRange range;
    rval = moab.get_entities_by_type( 0, t, range );
    CHECK_ERR(rval);
    if (!range.empty())
      rval = pcomm->filter_pstatus(range, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    CHECK_ERR(rval);
    counts[t] = range.size();
  }
}

void check_identical_mesh( MBInterface& mb1, MBInterface& mb2 )
{
  MBErrorCode rval;
  std::map<MBEntityHandle,MBEntityHandle> entmap;
  
    // match vertices by coordinate
  MBRange r1, r2;
  MBRange::iterator i1, i2;
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
  std::vector<MBEntityHandle> conn1, conn2;
  for (MBEntityType t = MBEDGE; t < MBENTITYSET; ++t) {
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
        for (std::vector<MBEntityHandle>::iterator j = conn2.begin(); j != conn2.end(); ++j)
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
  MBErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  
    // load and partition a .cub file
  MBCore moab_instance;
  MBInterface& moab = moab_instance;
  load_and_partition( moab, InputFile, true );
  
    // count number of owned entities of each type and sum over all procs
  count_owned_entities( moab, proc_counts );
  std::fill( all_counts, all_counts+MBENTITYSET, 0u );
  err = MPI_Allreduce( proc_counts, all_counts, MBENTITYSET, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  CHECK(!err);
  
    // do parallel write and on root proc do serial read of written file
  save_and_load_on_root( moab, "test_write_elements.h5m" );
  if (rank == 0) {
    for (MBEntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
      rval = moab.get_number_entities_by_type( 0, t, file_counts[t] );
      CHECK_ERR(rval);
    }
  }
  
    // check that sum of owned entities equals number of entities from serial read
    
  err = MPI_Bcast( file_counts, MBENTITYSET, MPI_INT, 0, MPI_COMM_WORLD );
  CHECK(!err);
  
  bool all_equal = true;
  for (MBEntityType t = MBVERTEX; t < MBENTITYSET; ++t) 
    if (file_counts[t] != all_counts[t])
      all_equal = false;
    
  if (rank == 0 && !all_equal) {
    std::cerr << "Type\tPartnd\tWritten" << std::endl;
    for (MBEntityType t = MBVERTEX; t < MBENTITYSET; ++t) 
      std::cerr << MBCN::EntityTypeName(t) << '\t' << all_counts[t] << '\t' << file_counts[t] << std::endl;
  }
  
  CHECK(all_equal);
  
    // on root processor, do serial read of original .cub file and compare
  
  if (rank == 0) {
    MBCore moab2;
    MBEntityHandle set;
    rval = moab2.load_file( InputFile, set );
    CHECK_ERR(rval);
    check_identical_mesh( moab, moab2 );
  }
}

bool check_sets_sizes( MBInterface& mb1, MBEntityHandle set1,
                       MBInterface& mb2, MBEntityHandle set2 )
{
  MBErrorCode rval;
  bool result = true;
  for (MBEntityType t = MBVERTEX; t < MBMAXTYPE; ++t) {
    int count1, count2;
    rval = mb1.get_number_entities_by_type( set1, t, count1 );
    CHECK_ERR(rval);
    rval = mb2.get_number_entities_by_type( set2, t, count2 );
    CHECK_ERR(rval);
    if (count1 != count2) {
      std::cerr << "Sets differ in number of " << MBCN::EntityTypeName(t)
                << " : " << count1 << " vs. " << count2 << std::endl;
      result = false;
    }
  }
  return result;
}

void test_write_shared_sets()
{
  int err, rank, size;
  MBErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  err = MPI_Comm_size( MPI_COMM_WORLD, &size );
  CHECK(!err);
  
  MBCore moab_instance;
  MBInterface& moab = moab_instance;
  load_and_partition( moab, InputFile );
  save_and_load_on_root( moab, "test_write_shared_sets.h5m" );

  if (rank != 0)
    return;
  
  MBCore moab2_instance;
  MBInterface& moab2 = moab2_instance;
  MBEntityHandle set;
  rval = moab2.load_file( InputFile, set );
  CHECK_ERR(rval);
  
  MBTag mattag1, mattag2;
  rval = moab.tag_get_handle( MATERIAL_SET_TAG_NAME, mattag1 );
  CHECK_ERR(rval);
  rval = moab2.tag_get_handle( MATERIAL_SET_TAG_NAME, mattag2 );
  CHECK_ERR(rval);
  
  MBRange matsets;
  rval = moab2.get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag2, 0, 1, matsets );
  CHECK_ERR(rval);
  for (MBRange::iterator i = matsets.begin(); i != matsets.end(); ++i) {
    int block_id;
    rval = moab2.tag_get_data( mattag2, &*i, 1, &block_id );
    CHECK_ERR(rval);
    
    MBRange tmpents;
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
    if (!KeepTmpFiles) 
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
