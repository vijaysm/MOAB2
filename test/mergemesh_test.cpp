#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MergeMesh.hpp"
#include <iostream>

#define STRINGIFY_(A) #A
#define STRINGIFY(A) STRINGIFY_(A)

using namespace moab;

const char* meshfile = STRINGIFY(MESHDIR) "/16_unmerged_hex.h5m";


int main( int , char** )
{
    Core moab_core;
    ErrorCode rval;
    Interface* iface = &moab_core;
    // can be generalized to load user defined input/output file
//    std::cout << "loading mesh file " << (std::string) meshfile << std::endl;
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
    bool merge_higher_dim_entities = true;

    rval = mm.merge_entities(ents, merge_tol, merge_higher_dim_entities);
    if (MB_SUCCESS != rval) {
        std::cerr << "Error in MergeMesh during merging entities" << std::endl;
        exit(2);
    }

    // Fixed for now
    const char *outfile = "mm_out.h5m";
    rval = iface->write_mesh( outfile);
    if (MB_SUCCESS != rval) {
        std::cerr << "Error saving file: " << outfile << std::endl;
        exit(2);
    }
    return 0;
}
