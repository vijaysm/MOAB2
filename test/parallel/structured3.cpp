#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ScdInterface.hpp"
#include "moab/HomXform.hpp"
#include "TestUtil.hpp"
#include <string>
#include <iomanip>
#include <iostream>
#include <cassert>
 
using namespace moab;
 
  // Number of cells in each direction:
const int NC = 2;
const int ITERS = 50;
 
void create_parallel_mesh();
 
int main(int argc, char *argv[]) 
{
  MPI_Init(&argc, &argv);
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
  pc.set_debug_verbosity(2);
  
    // create a structured mesh in parallel
  ScdBox *new_box;
  ScdParData par_data;
  par_data.pComm = &pc;
  par_data.gDims[0] = par_data.gDims[1] = par_data.gDims[2] = 0;
  par_data.gDims[3] = par_data.gDims[4] = par_data.gDims[5] = NC;
  par_data.gDims[5] = 1;
  par_data.partMethod = ScdParData::SQIJK;
  rval = scdi->construct_box(HomCoord(), HomCoord(), NULL, 0, // no vertex positions
                             new_box, NULL, // not locally periodic
                             &par_data, true, true); // assign global ids & resolve shared verts
  CHECK_ERR(rval);

  rval = pc.exchange_ghost_cells(-1, -1, 0, 0, true, true);
  CHECK_ERR(rval);

  pc.list_entities(0,-1);
  
  rval = pc.exchange_ghost_cells(-1, 0, 1, 0, true);
  if (MB_SUCCESS != rval) {
    std::string err;
    mbint.get_last_error(err);
    std::cerr << "Error: proc " << pc.rank() << ": " << err << std::endl;
  }
  CHECK_ERR(rval);

  pc.list_entities(0,-1);
  
  return;

    // Create a tag, used in exchange_tags
  Tag tag;
  int def_val = 1.0;
  rval = mbint.tag_get_handle("test_tag", 1, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE|MB_TAG_EXCL, &def_val);
  CHECK_ERR(rval);

  Range empty_range;
  for(int i = 0; i < ITERS; i++) {
    pc.exchange_tags(tag, empty_range);
    CHECK_ERR(rval);
  }
}
