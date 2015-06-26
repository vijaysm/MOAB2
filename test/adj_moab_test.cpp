/*This function tests the AHF datastructures on CST meshes*/
#include <iostream>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "TestUtil.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#endif

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#ifdef MOAB_HAVE_MPI
std::string read_options;
#endif

int number_tests_successful = 0;
int number_tests_failed = 0;

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
#ifdef MOAB_HAVE_MPI
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank==0)
          std::cout << "Success";
#else
      std::cout << "Success";
#endif
    number_successful++;
  } else {
    std::cout << "Failure";
    number_failed++;
  }
}


ErrorCode ahf_test(const char* filename)
{

    Core moab;
    Interface* mbImpl = &moab;
    ParallelComm *pc=NULL;
    MeshTopoUtil mtu(mbImpl);
    ErrorCode error;
    EntityHandle fileset;
    error = mbImpl->create_meshset(moab::MESHSET_SET, fileset);CHECK_ERR(error);

#ifdef MOAB_HAVE_MPI
    int procs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs > 1){
    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";

    error = mbImpl->load_file(filename,  &fileset, read_options.c_str());CHECK_ERR(error);
    }
    else if (procs == 1) {
#endif
    error = mbImpl->load_file(filename,  &fileset);CHECK_ERR(error);
#ifdef MOAB_HAVE_MPI
    }
#endif

    /*Create ranges for handles of explicit elements of the mixed mesh*/
    Range verts, edges, faces, cells;
    error = mbImpl->get_entities_by_dimension( fileset, 0, verts);CHECK_ERR(error);
    error = mbImpl->get_entities_by_dimension( fileset, 1, edges);CHECK_ERR(error);
    error = mbImpl->get_entities_by_dimension( fileset, 2, faces);CHECK_ERR(error);
    error = mbImpl->get_entities_by_dimension( fileset, 3, cells);CHECK_ERR(error);

    // Create an ahf instance  
#ifdef MOAB_HAVE_MPI
    pc = ParallelComm::get_pcomm(mbImpl,0);
    if (!pc)
      pc = new moab::ParallelComm(&moab, MPI_COMM_WORLD);
#endif
    HalfFacetRep ahf(&moab, pc, fileset);

    // Call the initialize function which creates the maps for each dimension
    error = ahf.initialize();CHECK_ERR(error);

     std::cout<<"Finished AHF initialization"<<std::endl;

    //Perform queries
    std::vector<EntityHandle> adjents;
    Range mbents, ahfents;

    //1D Queries //
    //IQ1: For every vertex, obtain incident edges
    if (edges.size()){
        for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
            adjents.clear();
            error = ahf.get_up_adjacencies( *i, 1, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents );CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(),mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    //NQ1:  For every edge, obtain neighbor edges
    if (edges.size()){
        for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
            adjents.clear();
            error = ahf.get_neighbor_adjacencies( *i, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mtu.get_bridge_adjacencies( *i, 0, 1, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    std::cout<<"Finished 1D queries"<<std::endl;
    // 2D Queries
    // IQ21: For every vertex, obtain incident faces
    if (faces.size()){
        for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
            adjents.clear();
            error = ahf.get_up_adjacencies( *i, 2, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    std::cout<<"Finished 2D queries: v2f"<<std::endl;

    //IQ22: For every edge, obtain incident faces
    if (edges.size() && faces.size()){
        for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
            adjents.clear();

            if ((*i-*edges.begin())==159)
              std::cout<<"Here now"<<std::endl;

            error = ahf.get_up_adjacencies( *i, 2, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);CHECK_ERR(error);

            if (adjents.size() != mbents.size())
              std::cout<<"EID = "<<(*i-*edges.begin())<<std::endl;

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    //NQ2: For every face, obtain neighbor faces
    if (faces.size()){
        for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
            adjents.clear();
            error = ahf.get_neighbor_adjacencies( *i, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mtu.get_bridge_adjacencies( *i, 1, 2, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    //DQ 21: For every face, obtain its edges
    if (edges.size() && faces.size()){
        for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
            adjents.clear();
            error = ahf.get_down_adjacencies( *i, 1, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    std::cout<<"Finished 2D queries: "<<std::endl;

    // 3D Queries
    //IQ 31: For every vertex, obtain incident cells
    if (cells.size()){
        for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
            adjents.clear();
            error = ahf.get_up_adjacencies( *i, 3, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    // IQ 32: For every edge, obtain incident cells
    if (edges.size()&&cells.size()){
        for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
            adjents.clear();

            if ((*i-*edges.begin())==12)
              std::cout<<"Here now"<<std::endl;

            error = ahf.get_up_adjacencies( *i, 3, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);CHECK_ERR(error);




            if (adjents.size() != mbents.size())
              {
                std::cout<<"EID = "<<(*i-*edges.begin())<<std::endl;
                for ( int j=0; j<(int)adjents.size(); j++)
                  std::cout<<"adjents["<<j<<"] = "<<(adjents[j]-*cells.begin())<<std::endl;
                for ( int j=0; j<(int)mbents.size(); j++)
                  {
                    const EntityHandle* conn;
                    int nv;
                    error = mbImpl->get_connectivity(mbents[j], conn, nv);MB_CHK_ERR(error);

                    std::cout<<"mbents["<<j<<"] = "<<(mbents[j]-*cells.begin())<<":: conn = [ "<<conn[0]<<", "<<conn[1]<<", "<<conn[2]<<", "<<conn[3]<<" ]"<<std::endl;
                  }
              }

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    //IQ33: For every face, obtain incident cells
    if (faces.size() && cells.size()){
        for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
            adjents.clear();
            error = ahf.get_up_adjacencies( *i, 3, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);CHECK_ERR(error);

            if (adjents.size() != mbents.size())
              {
                const EntityHandle* fconn;
                int nvf;
                error = mbImpl->get_connectivity(*i, fconn, nvf);MB_CHK_ERR(error);

                std::cout<<"FID = "<<(*i-*faces.begin())<<" :: conn = [ "<<fconn[0]<<", "<<fconn[1]<<", "<<fconn[2]<<" ]"<<std::endl;

                for ( int j=0; j<(int)adjents.size(); j++)
                  std::cout<<"adjents["<<j<<"] = "<<(adjents[j]-*cells.begin())<<std::endl;
                for ( int j=0; j<(int)mbents.size(); j++)
                  {
                    const EntityHandle* conn;
                    int nv;
                    error = mbImpl->get_connectivity(mbents[j], conn, nv);MB_CHK_ERR(error);

                    std::cout<<"mbents["<<j<<"] = "<<(mbents[j]-*cells.begin())<<":: conn = [ "<<conn[0]<<", "<<conn[1]<<", "<<conn[2]<<", "<<conn[3]<<" ]"<<std::endl;
                  }
              }


            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    //NQ3: For every cell, obtain neighbor cells
    if (cells.size()){
        for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
            adjents.clear();
            error = ahf.get_neighbor_adjacencies( *i, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mtu.get_bridge_adjacencies( *i, 2, 3, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }


    //DQ 31: For every cell, obtain its edges
    if (edges.size() && cells.size()){
        for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
            adjents.clear();
            error = ahf.get_down_adjacencies( *i, 1, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    //DQ 32: For every cell, obtain its faces
    if (faces.size() && cells.size()){
        for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
            adjents.clear();
            error = ahf.get_down_adjacencies( *i, 2, adjents);CHECK_ERR(error);
            mbents.clear();
            error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);CHECK_ERR(error);

            CHECK_EQUAL(adjents.size(), mbents.size());
            std::sort(adjents.begin(), adjents.end());
            std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
            mbents = subtract(mbents, ahfents);
            CHECK(!mbents.size());
          }
      }

    std::cout<<"Finished 3D queries"<<std::endl;
   error = ahf.deinitialize();CHECK_ERR(error);

    return MB_SUCCESS;

}

int main(int argc, char *argv[])
{

#ifdef MOAB_HAVE_MPI
    MPI_Init(&argc, &argv);

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    const char* filename = 0;
#ifdef MESHDIR
 #ifdef MOAB_HAVE_HDF5
    filename = STRINGIFY(MESHDIR) "/32hex_ef.h5m";
 #else
    filename = STRINGIFY(MESHDIR) "/hexes_mixed.vtk";
 #endif
#else
 #ifdef MOAB_HAVE_HDF5
    filename = "32hex_ef.h5m";
 #else
    filename = "hexes_mixed.vtk";
 #endif
#endif


    if (argc==1)
    {
#ifdef MOAB_HAVE_MPI
        if (rank == 0)
            std::cout<<"Using default input file:"<<filename<<std::endl;
#else
        std::cout<<"Using default input file:"<<filename<<std::endl;
#endif
    }

    else if (argc==2)
        filename = argv[1];
    else {
            std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
            return 1;
    }

    ErrorCode result;

#ifdef MOAB_HAVE_MPI
    if (rank == 0)
        std::cout<<" para_ahf_test: ";
#else
    std::cout<<"ahf_test:";
#endif

    result = ahf_test(filename);
    handle_error_code(result, number_tests_failed, number_tests_successful);
    std::cout<<"\n";

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return number_tests_failed;
}

