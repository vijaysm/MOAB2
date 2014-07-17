/*This function tests the AHF datastructures on CST meshes*/
#include <iostream>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "TestUtil.hpp"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"

using namespace moab;

#ifdef MESHDIR
std::string TestDir(STRINGIFY(MESHDIR));
#else
std::string TestDir(".");
#endif

std::string filename;

int number_tests_successful = 0;
int number_tests_failed = 0;

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
    std::cout << "Success";
    number_successful++;
  } else {
    std::cout << "Failure";
    number_failed++;
  }
}


ErrorCode ahf_test(Core *moab)
{

    Interface* mbImpl = &*moab;
    MeshTopoUtil mtu(mbImpl);

    ErrorCode error = mbImpl->load_file(filename.c_str());
    CHECK_ERR(error);

    /*Create ranges for handles of explicit elements of the mixed mesh*/
    Range verts, edges, faces, cells;
    error = mbImpl->get_entities_by_dimension( 0, 0, verts);
    error = mbImpl->get_entities_by_dimension( 0, 1, edges);
    error = mbImpl->get_entities_by_dimension( 0, 2, faces);
    error = mbImpl->get_entities_by_dimension( 0, 3, cells);

    // Create an ahf instance
    HalfFacetRep ahf(&*moab);

    // Call the initialize function which creates the maps for each dimension
    ahf.initialize();

    //ahf.print_tags();

    //Perform queries
    std::vector<EntityHandle> adjents;
    Range mbents, ahfents;

    //1D Queries //
    //IQ1: For every vertex, obtain incident edges
    for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 1, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents );
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(),mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //NQ1:  For every edge, obtain neighbor edges
    for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
        adjents.clear();
        error = ahf.get_neighbor_adjacencies( *i, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mtu.get_bridge_adjacencies( *i, 0, 1, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    // 2D Queries

    // IQ21: For every vertex, obtain incident faces
    for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 2, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //IQ22: For every edge, obtain incident faces
    for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 2, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //NQ2: For every face, obtain neighbor faces
    for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
        adjents.clear();
        error = ahf.get_neighbor_adjacencies( *i, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mtu.get_bridge_adjacencies( *i, 1, 2, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //DQ 21: For every face, obtain its edges
    for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
        adjents.clear();
        error = ahf.get_down_adjacencies( *i, 1, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    // 3D Queries
    //IQ 31: For every vertex, obtain incident cells
    for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 3, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    // IQ 32: For every edge, obtain incident cells
    for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 3, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //IQ33: For every face, obtain incident cells
    for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
        adjents.clear();
        error = ahf.get_up_adjacencies( *i, 3, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //NQ3: For every cell, obtain neighbor cells
    for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
        adjents.clear();
        error = ahf.get_neighbor_adjacencies( *i, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mtu.get_bridge_adjacencies( *i, 2, 3, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }


    //DQ 31: For every cell, obtain its edges
    for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
        adjents.clear();
        error = ahf.get_down_adjacencies( *i, 1, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }

    //DQ 32: For every cell, obtain its faces
    for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
        adjents.clear();
        error = ahf.get_down_adjacencies( *i, 2, adjents);
        CHECK_ERR(error);
        mbents.clear();
        error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
        CHECK_ERR(error);

        CHECK_EQUAL(adjents.size(), mbents.size());

        std::sort(adjents.begin(), adjents.end());
        std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
        mbents = subtract(mbents, ahfents);
        CHECK(!mbents.size());
    }



    ahf.deinitialize();

    return MB_SUCCESS;

}

int main(int argc, char *argv[])
{
    filename = TestDir + "/hexes_mixed.vtk";

    if (argc==1)
        std::cout<<"Using default input file:"<<filename<<std::endl;
    else if (argc==2)
        filename = argv[1];
    else {
            std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
            return 1;
    }

    Core moab;
    ErrorCode result;

    std::cout<<" ahf_test: ";
    result = ahf_test(&moab);
    handle_error_code(result, number_tests_failed, number_tests_successful);
    std::cout<<"\n";

    return number_tests_failed;
}

