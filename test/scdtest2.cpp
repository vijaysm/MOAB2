#include "moab/ScdInterface.hpp"
#include "moab/Core.hpp"
#include "TestUtil.hpp"

#include <iostream>

#ifdef USE_MPI
#include "moab/ParallelComm.hpp"
#endif

using namespace moab;

void test_john()
{
  Core moab;
  ScdInterface *scd;
  ErrorCode rval = moab.Interface::query_interface(scd);
  CHECK_ERR(rval);
  int gdims[] = {0, 0, 0, 48, 40, 0};
  int ldims[6];

  std::cout << "x dims are " << gdims[3] << "\n";
  std::cout << "y dims are " << gdims[4] << "\n";
  
#ifdef USE_MPI
  ParallelComm pcomm(&moab, MPI_COMM_WORLD);
  int procs = pcomm.proc_config().proc_size(),
      rank = pcomm.proc_config().proc_rank();
  ScdParData pardata;
  std::copy(gdims, gdims+6, pardata.gDims);
  pardata.gPeriodic[0] = pardata.gPeriodic[1] = 0;
  pardata.partMethod = ScdParData::SQIJ;
  int pDims[2];
  
  rval = ScdInterface::compute_partition(procs, rank, pardata, ldims, NULL, pDims); CHECK_ERR(rval);
  std::cout << "processors in x are " << pDims[0] << "\n";
  std::cout << "processors in y are " << pDims[1] << "\n";
#else
  std::copy(gdims, gdims+6, ldims);
#endif

  std::cout << "local dims are " 
            << ldims[0] << " " << ldims[1] << " " << ldims[2] << " " 
            << ldims[3] << " " << ldims[4] << " " << ldims[5] << "\n";

  ScdBox* box;

#ifdef USE_MPI
  rval = scd->construct_box(HomCoord(ldims, 3), HomCoord(ldims+3, 3), NULL, 0, box, NULL, &pardata); CHECK_ERR(rval);
  rval = pcomm.resolve_shared_ents(0); CHECK_ERR(rval);
#else
  rval = scd->construct_box(HomCoord(ldims, 3), HomCoord(ldims+3, 3), NULL, 0, box); CHECK_ERR(rval);
#endif

  moab.list_entities(0,0);
  std::vector<EntityHandle> conn;
  Range quads, verts;
  rval = moab.get_entities_by_type(0, MBQUAD, quads); CHECK_ERR(rval);
  rval = moab.get_adjacencies(quads, 0, false, verts, Interface::UNION); CHECK_ERR(rval);
  Tag gid_tag;
  rval = moab.tag_get_handle("GLOBAL_ID", gid_tag); CHECK_ERR(rval);
  int *qgids, *vgids, count;
  rval = moab.tag_iterate(gid_tag, quads.begin(), quads.end(), count, (void*&)qgids); 
  CHECK_ERR(rval); CHECK_EQUAL(count, (int)quads.size());
  rval = moab.tag_iterate(gid_tag, verts.begin(), verts.end(), count, (void*&)vgids); 
  CHECK_ERR(rval); CHECK_EQUAL(count, (int)verts.size());
  for (Range::iterator rit = quads.begin(); rit != quads.end(); rit++) {
    rval = moab.get_connectivity(&(*rit), 1, conn); CHECK_ERR(rval);
    std::cout << "Connectivity array for quad " << qgids[*rit-*quads.begin()] << " is: ";
    for (int i = 0; i < 4; i++) std::cout << vgids[conn[i]-*verts.begin()] << " ";
    std::cout << std::endl;
  }
  
}

int main(int, char**) 
{
    // test partition methods
  RUN_TEST(test_john);
}
