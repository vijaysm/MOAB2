#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/HomXform.hpp"
#include "moab/ProgOptions.hpp"
#include "MBTagConventions.hpp"
#include "TestUtil.hpp"
#include <string>
#include <iomanip>
#include <iostream>
#include <cassert>
 
using namespace moab;
 
  // Number of cells in each direction:
int NC;

/* This mesh creates a box that is NCxNCxNC in global dimension, partitioned among processors
 * automatically using ScdInterface's SQIJK algorithm.  It checks to make sure there are enough
 * procs to support this partition method.  After mesh creation, shared vertex resolution is done,
 * then ghost exchange is done.
 */
const int ITERS = 50;
 
void create_parallel_mesh();
 
int main(int argc, char *argv[]) 
{
  MPI_Init(&argc, &argv);

  ProgOptions po;
  po.addOpt<int>( "int,i", "Number of intervals on a side" );
  po.parseCommandLine( argc, argv );
  if(!po.getOpt( "int", &NC )) NC = 4;
  
  int err = RUN_TEST(create_parallel_mesh);
 
  MPI_Finalize();
  return err;
}
 
void create_parallel_mesh() 
{
  Core mbint;
  ParallelComm pc(&mbint, MPI_COMM_WORLD);
  ScdInterface *scdi;
  ErrorCode rval = mbint.query_interface(scdi);
  CHECK_ERR(rval);
    //pc.set_debug_verbosity(2);
  
    // create a structured mesh in parallel
  ScdBox *new_box;
  ScdParData par_data;
  par_data.pComm = &pc;
  par_data.gDims[0] = par_data.gDims[1] = par_data.gDims[2] = 0;
  par_data.gDims[3] = par_data.gDims[4] = par_data.gDims[5] = NC;
  if ((par_data.gDims[3]-par_data.gDims[0])*(par_data.gDims[3]-par_data.gDims[0])*(par_data.gDims[3]-par_data.gDims[0]) < (int)pc.size()) {
    std::cerr << "Too few processors for this number of elements." << std::endl;
    CHECK_ERR(MB_FAILURE);
  }
    
  par_data.partMethod = ScdParData::SQIJK;

    // timing data
  double times[5]; // tstart, tvert, tnonvert, tghost, titer;
  times[0] = MPI_Wtime();
  rval = scdi->construct_box(HomCoord(), HomCoord(), NULL, 0, // no vertex positions
                             new_box, NULL, // not locally periodic
                             &par_data, true, false); // assign global ids, don't resolve shared verts
  CHECK_ERR(rval);

    // get global id tag
  Tag tag;
  rval = mbint.tag_get_handle(GLOBAL_ID_TAG_NAME, tag);
  CHECK_ERR(rval);
  
    // resolve shared verts
  std::cout << "Resolving shared ents..." << std::endl;
  rval = pc.resolve_shared_ents(new_box->box_set(), -1, 0, &tag);
  CHECK_ERR(rval);
  times[1] = MPI_Wtime();
  
  std::cout << "Exchanging ghost cells..." << std::endl;
  rval = pc.exchange_ghost_cells(-1, -1, 0, 0, true, true);
  CHECK_ERR(rval);
  times[2] = MPI_Wtime();

//  pc.list_entities(0,-1);
  
  rval = pc.exchange_ghost_cells(-1, 0, 1, 0, true);
  if (MB_SUCCESS != rval) {
    std::string err;
    mbint.get_last_error(err);
    std::cerr << "Error: proc " << pc.rank() << ": " << err << std::endl;
  }
  CHECK_ERR(rval);
  times[3] = MPI_Wtime();

//  pc.list_entities(0,-1);
  
    // Create a tag, used in exchange_tags
  int def_val = 1.0;
  rval = mbint.tag_get_handle("test_tag", 1, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE|MB_TAG_EXCL, &def_val);
  CHECK_ERR(rval);

  Range empty_range;
  if (!pc.rank()) std::cout << "Exchanging tags: ";
  for(int i = 0; i < ITERS; i++) {
    if (!pc.rank()) std::cout << i << ";";
    pc.exchange_tags(tag, empty_range);
    CHECK_ERR(rval);
  }
  if (!pc.rank()) std::cout << std::endl;
  times[4] = MPI_Wtime();

  for (int i = 4; i >= 1; i--) times[i] -= times[i-1];

  double tottimes[5] = {0.0};
  MPI_Reduce(times+1, tottimes+1, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  if (!pc.rank()) 
    std::cout << "Times:             " << std::endl
              << "Create:            " << times[1] << std::endl
              << "Resolve verts:     " << times[2] << std::endl
              << "Resolve non-verts: " << times[3] << std::endl
              << "Exchange ghosts:   " << times[4] << std::endl;
}
