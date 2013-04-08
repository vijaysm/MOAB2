// A test file for Subset Normalization
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/Core.hpp"
#include <iostream>
#include <sstream>
#include <cstdlib>

using namespace moab;

bool debug = true;

// Error routines for use with MOAB API
#define CHKERR(CODE, MSG)                                 \
  do {                                                    \
    if (MB_SUCCESS != (CODE)) {                           \
      std::string errstr;  mbi->get_last_error(errstr);   \
      std::cerr << errstr << std::endl;                   \
      std::cerr << MSG << std::endl;                      \
      MPI_Finalize();                                     \
    }                                                     \
  } while(false)

// Error routines for use with MPI API
#define MPICHKERR(CODE, MSG)                              \
  do {                                                    \
    if (0 != CODE) {                                      \
      std::cerr << MSG << std::endl;                      \
      MPI_Finalize();                                     \
    }                                                     \
  } while(false)


// Function to parse input parameters
ErrorCode get_file_options(int argc, char **argv,
                           const char** filename,
                           const char** tagName,
                           double &tagValues)
{
  int npos = 1;

  assert(argc >= 2);

  // get mesh filename
  *filename = argv[npos++];

  // get tag selection options
  if (argc > 3) {
    *tagName = argv[npos++];
    tagValues = atof(argv[npos++]);
  }
  else {
    *tagName = "USERTAG";
    tagValues = 1.0;
  }

  return MB_SUCCESS;
}

#define dbgprint(MSG)                                \
  do {                                              \
      if (!rank) std::cerr << MSG << std::endl;     \
  } while(false)

//
// Start of main test program
//
int main(int argc, char **argv)
{
  ErrorCode err;
  int ierr, rank;
  const char *filename, *tagName;
  double tagValue;
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_DISTRIBUTE;PARALLEL_GHOSTS=3.0.1;PARALLEL_COMM=0";

  // Initialize MPI first
  ierr = MPI_Init(&argc, &argv);
  MPICHKERR(ierr, "MPI_Init failed");

  // Print usage if not enough arguments
  if (argc < 2) {
    std::cerr << "Usage: ";
    std::cerr << argv[0] << " <file_name> <tag_name> <tag_value>" << std::endl;
    std::cerr << "file_name    : mesh file name" << std::endl;
    std::cerr << "tag_name     : name of tag to add to mesh" << std::endl;
    std::cerr << "tag_value    : a double valued string to set for highest-dimensional entities in the mesh for the named tag" << std::endl;

    ierr = MPI_Finalize();
    MPICHKERR(ierr, "MPI_Finalize failed; Aborting");

    return 1;
  }

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPICHKERR(ierr, "MPI_Comm_rank failed");

  // Create the moab instance
  Interface *mbi = new Core();
  CHKERR(NULL == mbi, "MOAB constructor failed");

  // Get the input options
  dbgprint( "Getting options..." );
  err = get_file_options(argc, argv, &filename, &tagName, tagValue);
  CHKERR(err, "get_file_options failed");

  // Print out the input parameters
  dbgprint( "    Input Parameters - " );
  dbgprint( "      Filenames: " << filename );
  dbgprint( "      Tag: Name=" << tagName << " Value=" << tagValue );

  // Create root sets for each mesh.  Then pass these
  // to the load_file functions to be populated.
  EntityHandle rootset, partnset;
  err = mbi->create_meshset(MESHSET_SET, rootset);
  CHKERR(err, "Creating root set failed");
  err = mbi->create_meshset(MESHSET_SET, partnset);
  CHKERR(err, "Creating partition set failed");

  // Create the parallel communicator object with the partition handle associated with MOAB
  ParallelComm *parallel_communicator = ParallelComm::get_pcomm( mbi, partnset, &comm );

  // Load the file from disk with given options
  err = mbi->load_file( filename, &rootset, read_options.c_str() );
  CHKERR(err, "MOAB::load_file failed");

  // Get tag handles for passed in tags
  dbgprint( "Getting tag handle " << tagName << "..." );
  Tag tagReduce, tagExchange;
  {
    std::stringstream sstr;
    sstr << tagName << "_RED";
    err = mbi->tag_get_handle(sstr.str().c_str(), 1, MB_TYPE_DOUBLE, tagReduce, MB_TAG_CREAT|MB_TAG_DENSE, &tagValue);
    CHKERR(err, "Retrieving tag handles failed");
    sstr.str(""); sstr << tagName << "_EXC";
    err = mbi->tag_get_handle(sstr.str().c_str(), 1, MB_TYPE_INTEGER, tagExchange, MB_TAG_CREAT|MB_TAG_DENSE, &tagValue);
    CHKERR(err, "Retrieving tag handles failed");
  }

  // Set local tag data for reduction
  {
    Range partEnts;
    err = parallel_communicator->get_part_entities(partEnts);
    CHKERR(err, "ParallelComm::get_part_entities failed");
    // Output what is in current partition sets
    std::cout << "[" << rank << "]: Number of Partitioned entities: " <<  partEnts.size() << std::endl;
    MPI_Barrier(comm);

    std::vector<double> tagValues(partEnts.size(), tagValue*(rank+1));
    err = mbi->tag_set_data(tagReduce, partEnts, &tagValues[0]);
    CHKERR(err, "Setting local tag data failed during reduce phase");

    Range dummy;
    dbgprint( "Reducing tags between processors " );
    err = parallel_communicator->reduce_tags(tagReduce, MPI_SUM, dummy/*partEnts*/);
    CHKERR(err, "Reducing tags between processors failed");

    partEnts.clear();
  }

  // Set local tag data for exchange
  {
    Range partEnts, dimEnts;
    for (int dim = 0; dim <= 2; dim++) {
      err = mbi->get_entities_by_dimension(rootset, dim, dimEnts, false);
      std::vector<int> tagValues(dimEnts.size(), static_cast<int>(tagValue)*(rank+1)*(dim+1));
      err = mbi->tag_set_data(tagExchange, dimEnts, &tagValues[0]);
      CHKERR(err, "Setting local tag data failed during exchange phase");
      partEnts.merge(dimEnts);
    }

    dbgprint( "Exchanging tags between processors " );
    err = parallel_communicator->exchange_tags(tagExchange, partEnts);
    CHKERR(err, "Exchanging tags between processors failed");
  }

  // Write out to output file to visualize reduction/exchange of tag data
  mbi->write_file("test.h5m", "H5M", "PARALLEL=WRITE_PART");

  // Done, cleanup
  delete mbi;

  dbgprint( "********** reduce_exchange_tags DONE! **********" );
  MPI_Finalize();
  return 0;
}
