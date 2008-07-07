#include "TestUtil.hpp"

#include "MBCore.hpp"
#include "MBParallelComm.hpp"
#include "MBTagConventions.hpp"
#include "MBCN.hpp"

#include <iostream>
#include <mpi.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#ifdef SRCDIR
const char* InputFile = STRINGIFY(SRCDIR) "/ptest.cub";
#else
const char* InputFile = "ptest.cub";
#endif

void load_and_partition( MBInterface& moab, const char* filename );

void save_and_load_on_root( MBInterface& moab, const char* tmp_filename );

void check_identical_mesh( MBInterface& moab1, MBInterface& moab2 );

void test_write_elements();

void test_write_shared_sets();

bool KeepTmpFiles = false;

int main( int argc, char* argv[] )
{
  int err = MPI_Init( &argc, &argv );
  CHECK(!err);

  if (argc == 2 && !strcmp(argv[1],"-k"))
    KeepTmpFiles = true;
  else if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " [-k]" << std::endl;
    return 1;
  }
  
  int result = 0;
  result += RUN_TEST( test_write_elements );
  result += RUN_TEST( test_write_shared_sets );
  
  MPI_Finalize();
  return result;
}

void load_and_partition( MBInterface& moab, const char* filename )
{
  MBErrorCode rval;
  MBEntityHandle set;
  
  rval = moab.load_file( filename, set, 
                         "PARALLEL=READ_DELETE;"
                         "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;"
                         "PARTITION_DISTRIBUTE;"
                         "PARALLEL_RESOLVE_SHARED_ENTS" );
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
  
  moab.delete_mesh();
  if (procnum == 0) {
    rval = moab.load_file( tmp_filename, set );
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
    rval = pcomm->remove_nonowned_shared( range, -1, true, false );
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
  int err, rank, size;
  MBErrorCode rval;
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  CHECK(!err);
  err = MPI_Comm_size( MPI_COMM_WORLD, &size );
  CHECK(!err);
  
  MBCore moab_instance;
  MBInterface& moab = moab_instance;
  load_and_partition( moab, InputFile );
  
  count_owned_entities( moab, proc_counts );
  std::fill( all_counts, all_counts+MBENTITYSET, 0u );
  err = MPI_Allreduce( proc_counts, all_counts, MBENTITYSET, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  CHECK(!err);
  
  save_and_load_on_root( moab, "test_write_elements.h5m" );
  if (rank == 0) {
    for (MBEntityType t = MBVERTEX; t < MBENTITYSET; ++t) {
      rval = moab.get_number_entities_by_type( 0, t, file_counts[t] );
      CHECK_ERR(rval);
    }
  }
  
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
 
  
    
  
