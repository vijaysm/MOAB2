/*This function tests the AHF datastructures on CST meshes*/

#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "../TestUtil.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

std::string read_options;

int number_tests_successful = 0;
int number_tests_failed = 0;

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
    //std::cout << "Success";
    number_successful++;
  } else {
    //std::cout << "Failure";
    number_failed++;
  }
}

ErrorCode para_ahf_test(const char* filename)
{
    Core moab;
    Interface* mbImpl = &moab;
    MeshTopoUtil mtu(mbImpl);

    //std::cout<<"DB: Begin"<<std::endl;

    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

    ErrorCode error = mbImpl->load_file(filename, 0, read_options.c_str());
    CHECK_ERR(error);

    /*Create ranges for handles of explicit elements of the mixed mesh*/
    Range local_verts, local_edges, local_faces, local_cells;
    error = mbImpl->get_entities_by_dimension( 0, 0, local_verts);
    error = mbImpl->get_entities_by_dimension( 0, 1, local_edges);
    error = mbImpl->get_entities_by_dimension( 0, 2, local_faces);
    error = mbImpl->get_entities_by_dimension( 0, 3, local_cells);

    // Create an ahf instance
    HalfFacetRep ahf(&moab);

    // Call the initialize function which creates the maps for each dimension
    ahf.initialize();

    //Perform queries
    std::vector<EntityHandle> adjents;
    Range mbents, ahfents;

    //1D Queries //
    //IQ1: For every vertex, obtain incident edges
    for (Range::iterator i = local_verts.begin(); i != local_verts.end(); ++i) {
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
    for (Range::iterator i = local_edges.begin(); i != local_edges.end(); ++i) {
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
    for (Range::iterator i = local_verts.begin(); i != local_verts.end(); ++i) {
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
    for (Range::iterator i = local_edges.begin(); i != local_edges.end(); ++i) {
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
    for (Range::iterator i = local_faces.begin(); i != local_faces.end(); ++i) {
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
    for (Range::iterator i = local_faces.begin(); i != local_faces.end(); ++i) {
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
    for (Range::iterator i = local_verts.begin(); i != local_verts.end(); ++i) {
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
    for (Range::iterator i = local_edges.begin(); i != local_edges.end(); ++i) {
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
    for (Range::iterator i = local_faces.begin(); i != local_faces.end(); ++i) {
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
    for (Range::iterator i = local_cells.begin(); i != local_cells.end(); ++i) {
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
    for (Range::iterator i = local_cells.begin(); i != local_cells.end(); ++i) {
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
    for (Range::iterator i = local_cells.begin(); i != local_cells.end(); ++i) {
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
    MPI_Init(&argc, &argv);

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const char* filename = 0;

#ifdef SRCDIR
        filename = STRINGIFY(MESHDIR) "/32hex_ef.h5m";
#else
        filename = "32hex_ef.h5m";
#endif

    if (argc==1){
        if (rank==0)
            std::cout<<"Using default file:"<<filename<<std::endl;
    }
    else if (argc==2)
        filename = argv[1];
    else {
        std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
        exit(1);
    }

    ErrorCode result;
    if (rank == 0)
        std::cout<<" para_ahf_test: ";

    result = para_ahf_test(filename);
    handle_error_code(result, number_tests_failed, number_tests_successful);

    if (rank == 0){
        if (!number_tests_failed)
            std::cout << "Success" << std::endl;
        else
            std::cout << number_tests_failed << " TESTS FAILED!" << std::endl;
    }

    MPI_Finalize();

    return number_tests_failed;
}

