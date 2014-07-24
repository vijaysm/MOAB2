/** @example TestErrorHandling.cpp \n
 * Description: This example tests MOAB's trace back error handler in serial. \n
 *
 * <b>To run</b>: ./TestErrorHandling <test_case_num(1 to 3)> \n
 */

#include "moab/Core.hpp"
#ifdef USE_MPI
#include "moab_mpi.h"
#endif

#include <iostream>

using namespace moab;
using namespace std;

// In this test case, an error MB_NOT_IMPLEMENTED is returned by MOAB
ErrorCode TestErrorHandling_1()
{
  Core moab;
  Interface& mb = moab;

  // Load a CAM-FV file and read a variable on edges (not supported yet)
  string test_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, "VARIABLE=US");CHK_ERR(rval);

  return MB_SUCCESS;
}

// In this test case, an error MB_TYPE_OUT_OF_RANGE is returned by MOAB
ErrorCode TestErrorHandling_2()
{
  Core moab;
  Interface& mb = moab;

  // Load a HOMME file with an invalid GATHER_SET option
  string test_file = string(MESH_DIR) + string("/io/homme3x3458.t.3.nc");
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, "VARIABLE=T;GATHER_SET=0.1");CHK_ERR(rval);

  return MB_SUCCESS;
}

// In this test case, an error MB_FAILURE is returned by MOAB
ErrorCode TestErrorHandling_3()
{
  Core moab;
  Interface& mb = moab;

  // Load a CAM-FV file with  and a NULL file set
  string test_file = string(MESH_DIR) + string("/io/fv3x46x72.t.3.nc");
  ErrorCode rval = mb.load_file(test_file.c_str(), NULL, "NOMESH;VARIABLE=");CHK_ERR(rval);

  return MB_SUCCESS;
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <test_case_num(1 to 3)>" << endl;
    return 0;
  }

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Initialize error handler, optional for this example (using moab instances)
  MBErrorHandler_Init();

  ErrorCode rval = MB_SUCCESS;

  int test_case_num = atoi(argv[1]);
  switch (test_case_num) {
    case 1:
      rval = TestErrorHandling_1();CHK_ERR(rval);
      break;
    case 2:
      rval = TestErrorHandling_2();CHK_ERR(rval);
      break;
    case 3:
      rval = TestErrorHandling_3();CHK_ERR(rval);
      break;
    default:
      break;
  }

  // Finalize error handler, optional for this example (using moab instances)
  MBErrorHandler_Finalize();

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
