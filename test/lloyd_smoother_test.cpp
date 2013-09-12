#include "moab/Core.hpp"
#include "moab/LloydSmoother.hpp"
#include "TestUtil.hpp"

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
std::string TestDir(".");
#endif

std::string filename = TestDir + "/surfrandomtris-4part.h5m";

using namespace moab;

int main(int argc, char**argv) 
{
  if (argc > 1) filename = std::string(argv[1]);
  Core mb;
  ErrorCode rval = mb.load_file(filename.c_str());
  CHECK_ERR(rval);
  
  Range elems;
  rval = mb.get_entities_by_dimension(0, 3, elems);
  CHECK_ERR(rval);
  if (elems.empty()) {
    rval = mb.get_entities_by_dimension(0, 2, elems);
    CHECK_ERR(rval);
  }
  if (elems.empty()) {
    std::cout << "Mesh must have faces or regions for this test." << std::endl;
    CHECK_ERR(MB_FAILURE);
  }

  LloydSmoother ll(&mb, NULL, elems);
  ll.report_its(10);
  rval = ll.perform_smooth();
  CHECK_ERR(rval);
  
  std::cout << "Mesh smoothed in " << ll.num_its() << " iterations." << std::endl;
}

