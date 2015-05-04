#include "moab/ScdInterface.hpp"
#include "moab/Core.hpp"
#include "TestUtil.hpp"

#include <iostream>

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

using namespace moab;

int gdims[6], np;

void test_sqijk()
{
  Core moab;
  ScdInterface *scd;
  ErrorCode rval = moab.Interface::query_interface(scd);
  CHECK_ERR(rval);
  int ldims[6];
  ScdParData par_data;
  std::copy(gdims, gdims+6, par_data.gDims);
  for (int i = 0; i < 3; i++) par_data.gPeriodic[i] = 0;
  par_data.partMethod = ScdParData::SQIJK;

  std::cout << "gDims =  (" << par_data.gDims[0] << "," << par_data.gDims[1] << "," << par_data.gDims[2] << ")--("
            << par_data.gDims[3] << "," << par_data.gDims[4] << "," << par_data.gDims[5] << ")" << std::endl;

  int lperiodic[3], pijk[3];
  
  rval = ScdInterface::compute_partition(np, 0, par_data, ldims, lperiodic, pijk); CHECK_ERR(rval);
  
  
  std::cout << "#proc in 3 directions = (" << pijk[0] << "," << pijk[1] << "," << pijk[2] << ")" << std::endl;
  std::cout << "local dims are (" 
            << ldims[0] << "," << ldims[1] << "," << ldims[2] << ")--(" 
            << ldims[3] << "," << ldims[4] << "," << ldims[5] << ")\n";
}

int main(int argc, char**argv) 
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <#proc> [<imax> [<jmax> <kmax>]]" << std::endl;
    std::cout << "Using default parameters for autotest purposes." << std::endl;
    np = 4;
    gdims[0] = gdims[1] = gdims[2] = 0;
    gdims[3] = gdims[4] = gdims[5] = 100;
  }
  else if (argc < 3) {
    np = atoi(argv[1]);
    gdims[0] = gdims[1] = gdims[2] = 0;
    gdims[3] = gdims[4] = gdims[5] = 100;
  }
  else if (argc < 4) {
    np = atoi(argv[1]);
    gdims[0] = gdims[1] = gdims[2] = 0;
    gdims[3] = gdims[4] = gdims[5] = atoi(argv[2]);
  }
  else if (argc < 6) {
    np = atoi(argv[1]);
    gdims[0] = gdims[1] = gdims[2] = 0;
    gdims[3] = atoi(argv[2]);
    gdims[4] = atoi(argv[3]);
    gdims[5] = atoi(argv[4]);
  }

    // test partition method
  RUN_TEST(test_sqijk);
}
