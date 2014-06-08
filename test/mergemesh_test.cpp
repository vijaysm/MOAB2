#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MergeMesh.hpp"
#include <iostream>

#define STRINGIFY_(A) #A
#define STRINGIFY(A) STRINGIFY_(A)

using namespace moab;

 const char* meshfile = STRINGIFY(MESHDIR) "/16_unmerged_hex.h5m";
 const char *outfile = "mm_out.h5m";

int main( int argc, char** argv)
{
    Core moab_core;
    ErrorCode rval;
    Interface* iface = &moab_core;
    // can be generalized to load user defined input/output file

    if (argc>1)
      meshfile = argv[1];
    if (argc>2)
      outfile= argv[2];
    rval = iface->load_mesh(meshfile);
    if (MB_SUCCESS != rval) {
        std::cerr << "Error reading file: " << meshfile << std::endl;
        exit(2);
    }
    int dim = 3;
    moab::Range ents;
    iface->get_entities_by_dimension(0, dim, ents);

    MergeMesh mm(iface);
    double merge_tol = 1e-3;

    rval = mm.merge_entities(ents, merge_tol);
    if (MB_SUCCESS != rval) {
        std::cerr << "Error in MergeMesh during merging entities" << std::endl;
        exit(2);
    }

    // Fixed for now

    rval = iface->write_file( outfile);
    if (MB_SUCCESS != rval) {
        std::cerr << "Error saving file: " << outfile << std::endl;
        exit(2);
    }
    return 0;
}
