/*This function profiles the performance of the AHF datastructure */
#include <iostream>
#include <assert.h>
#include <time.h>
#include <vector>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "RefineMesh/moab/NestedRefine.hpp"
#include "../TestUtil.hpp"
#include <sys/time.h>

#ifdef USE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "ReadParallel.hpp"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#endif

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

#ifdef USE_MPI
std::string read_options;
#endif

#ifdef MESHDIR
std::string TestDir( STRINGIFY(MESHDIR) );
#else
#error MESHDIR needs to be defined for running unit tests
#endif

int number_tests_successful = 0;
int number_tests_failed = 0;

struct mesh_mem
{
   unsigned long long total_storage;
   unsigned long long amortized_total_storage;
   unsigned long long entity_storage;
   unsigned long long amortized_entity_storage;
   unsigned long long adjacency_storage;
   unsigned long long amortized_adjacency_storage;
   unsigned long long tag_storage;
   unsigned long long amortized_tag_storage;
};

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
#ifdef USE_MPI
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

double wtime() {
  double y = -1;
  struct timeval cur_time;
  gettimeofday(&cur_time, NULL);
  y = (double)(cur_time.tv_sec) + (double)(cur_time.tv_usec)*1.e-6;
  return (y);
}

ErrorCode umr_perf_test(const char* filename, int *level_degrees, int num_levels)
{
  ErrorCode error;
  Core mb;
  Interface* mbImpl = &mb;

#ifdef USE_MPI
    int procs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs > 1){
    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;";

    error = mbImpl->load_file(filename, 0, read_options.c_str());
    CHECK_ERR(error);
    }
    else if (procs == 1) {
#endif
    error = mbImpl->load_file(filename);
    CHECK_ERR(error);
#ifdef USE_MPI
    }
#endif

    mesh_mem *umem = new mesh_mem[num_levels+2];

    //Storage Costs before mesh hierarchy generation
    std::cout<<std::endl;
    unsigned long long sTotS, sTAS, sES, sAES, sAS, sAAS, sTS, sATS;
    sTotS = sTAS = sES = sAES = sAS = sAAS = sTS = sATS = 0;
    mbImpl->estimated_memory_use(NULL, 0, &sTotS, &sTAS, &sES, &sAES, &sAS, &sAAS, NULL, 0, &sTS, &sATS);

    umem[0].total_storage = sTotS;
    umem[0].amortized_total_storage = sTAS;
    umem[0].entity_storage = sES;
    umem[0].amortized_entity_storage = sAES;
    umem[0].adjacency_storage = sAS;
    umem[0].amortized_adjacency_storage = sAAS;
    umem[0].tag_storage = sTS;
    umem[0].amortized_tag_storage = sATS;

    double time_start, time_avg, time_total;

    //Create ranges of entities in the initial mesh
    Range inverts, inedges, infaces, incells;
    error = mbImpl->get_entities_by_dimension( 0, 0, inverts);
    error = mbImpl->get_entities_by_dimension( 0, 1, inedges);
    error = mbImpl->get_entities_by_dimension( 0, 2, infaces);
    error = mbImpl->get_entities_by_dimension( 0, 3, incells);

    Range init_ents;
    int dim;
    if (inedges.size())
      {
        dim = 1;
        init_ents = inedges;
      }
    else if (infaces.size())
      {
        dim = 2;
        init_ents = infaces;
      }
    else if (incells.size())
      {
        dim = 3;
        init_ents = incells;
      }

    std::cout<<"Initial Mesh Size: "<<"NV = "<<inverts.size()<<", NE = "<<init_ents.size()<<std::endl;

    //Create an hm object and generate the hierarchy
    std::cout<<"Creating a hm object"<<std::endl;
    NestedRefine uref(&mb);
    EntityHandle *set = new EntityHandle[num_levels];

    std::cout<<"Starting hierarchy generation"<<std::endl;
    time_start = wtime();
    error = uref.generate_mesh_hierarchy(level_degrees, num_levels, set);
    CHECK_ERR(error);

    time_total = wtime() - time_start;
    std::cout<<"Total time:: generate mesh hierarchy:: L = "<<num_levels<<" :: "<<time_total<<std::endl;

    std::cout<<std::endl;
    sTotS = sTAS = sES = sAES = sAS = sAAS = sTS = sATS = 0;
    mbImpl->estimated_memory_use(NULL, 0, &sTotS, &sTAS, &sES, &sAES, &sAS, &sAAS, NULL, 0, &sTS, &sATS);

    umem[num_levels+1].total_storage = sTotS;
    umem[num_levels+1].amortized_total_storage = sTAS;
    umem[num_levels+1].entity_storage = sES;
    umem[num_levels+1].amortized_entity_storage = sAES;
    umem[num_levels+1].adjacency_storage = sAS;
    umem[num_levels+1].amortized_adjacency_storage = sAAS;
    umem[num_levels+1].tag_storage = sTS;
    umem[num_levels+1].amortized_tag_storage = sATS;


    //Loop over each mesh level and check its topological properties
    for (int l=0; l<num_levels; l++)
      {
        Range verts, ents;
        error = mbImpl->get_entities_by_type(set[l], MBVERTEX, verts);
        CHECK_ERR(error);
        error = mbImpl->get_entities_by_dimension(set[l], dim, ents);
        CHECK_ERR(error);

        std::cout<<std::endl;
        std::cout<<"Mesh size for level "<<l<<"  :: NV = "<<verts.size()<<", NE = "<<ents.size()<<std::endl;

        //Storage Costs
        std::cout<<std::endl;
        unsigned long long vTotS, vTAS, vES, vAES, vAS, vAAS, vTS, vATS;
        vTotS = vTAS = vES = vAES = vAS = vAAS = vTS = vATS = 0;
        mbImpl->estimated_memory_use(verts, &vTotS, &vTAS, &vES, &vAES, &vAS, &vAAS, NULL, 0, &vTS, &vATS);
        unsigned long long eTotS, eTAS, eES, eAES, eAS, eAAS, eTS, eATS;
        eTotS = eTAS = eES = eAES = eAS = eAAS = eTS = eATS = 0;
        mbImpl->estimated_memory_use(ents, &eTotS, &eTAS, &eES, &eAES, &eAS, &eAAS, NULL, 0, &eTS, &eATS);

        umem[l+1].total_storage = vTotS+eTotS;
        umem[l+1].amortized_total_storage = vTAS+eTAS;
        umem[l+1].entity_storage = vES+eES;
        umem[l+1].amortized_entity_storage = vAES+eAES;
        umem[l+1].adjacency_storage = vAS+eAS;
        umem[l+1].amortized_adjacency_storage = vAAS+eAAS;
        umem[l+1].tag_storage = vTS+eTS;
        umem[l+1].amortized_tag_storage = vATS+eATS;

        //Loop over all vertices and get their coordinates
        time_start = wtime();
        for (Range::iterator i = verts.begin(); i != verts.end(); ++i)
          {
            double coords[3];
            EntityHandle vid = *i;
            error = uref.get_coordinates(&vid, 1, (int)l, &coords[0]);
            CHECK_ERR(error);
          }
        time_total = wtime() - time_start;
        time_avg = time_total/(double)verts.size();

        std::cout<<"Class NR :: OPERATION:: get_coordinates"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;

        time_start = wtime();
        for (Range::iterator i = verts.begin(); i != verts.end(); ++i)
          {
            double coords[3];
            EntityHandle vid = *i;
            error = mbImpl->get_coords(&vid, 1, &coords[0]);
            CHECK_ERR(error);
          }
        time_total = wtime() - time_start;
        time_avg = time_total/(double)verts.size();

        std::cout<<"Class Core :: OPERATION:: get_coordinates"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;


        //Loop over all entities and get their connectivity
        time_start = wtime();
        for (Range::iterator i = ents.begin(); i != ents.end(); ++i)
          {
            std::vector<EntityHandle> conn;
            error = uref.get_connectivity(*i, l, conn);
            CHECK_ERR(error);
          }
        time_total = wtime() - time_start;
        time_avg = time_total/(double)ents.size();

        std::cout<<"Class NR :: OPERATION:: get_connectivity"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;

        time_start = wtime();
        for (Range::iterator i = ents.begin(); i != ents.end(); ++i)
          {
            std::vector<EntityHandle> conn;
            error = mbImpl->get_connectivity(&*i, 1, conn);
            CHECK_ERR(error);
          }
        time_total = wtime() - time_start;
        time_avg = time_total/(double)ents.size();

        std::cout<<"Class Core :: OPERATION:: get_connectivity"<<" :: Time_total = "<<time_total<<" :: Time_avg = "<<time_avg<<std::endl;
      }


    //Print out memory storage
    std::cout<<"MEMORY STORAGE"<<std::endl;
    for (int l=0; l<num_levels+2; l++)
      {
        if (l==0)
          std::cout<<"INITIAL MESH"<<std::endl;
        else if (l==num_levels+1)
          std::cout<<"MESH HIERARCHY"<<std::endl;
        else
          std::cout<<"MESH LEVEL L="<<l+1<<std::endl;

        std::cout<<std::endl;
        std::cout<<"Total storage = "<<umem[l].total_storage<<std::endl;
        std::cout<<"Total amortized storage = "<< umem[l].amortized_total_storage<<std::endl;
        std::cout<<"Entity storage = "<<umem[l+1].entity_storage<<std::endl;
        std::cout<<"Amortized entity storage = "<<umem[l].amortized_entity_storage<<std::endl;
        std::cout<<"Adjacency storage = "<< umem[l].adjacency_storage <<std::endl;
        std::cout<<"Amortized adjacency storage = "<<umem[l].amortized_adjacency_storage <<std::endl;
        std::cout<<"Tag storage = "<< umem[l].tag_storage<<std::endl;
        std::cout<<"Amortized tag storage = "<<umem[l].amortized_tag_storage <<std::endl;
        std::cout<<std::endl;
      }

    return MB_SUCCESS;
}

int main(int argc, char *argv[])
{

#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    const char* filename = 0;
#ifdef MESHDIR
 #ifdef HDF5_FILE
    filename = STRINGIFY(MESHDIR) "/hex_2048.h5m";
 #else
    filename = STRINGIFY(MESHDIR) "/hex_2048.vtk";
 #endif
#else
 #ifdef HDF5_FILE
    filename = "hex_2048.h5m";
 #else
    filename = "hex_2048.vtk";
 #endif
#endif

    if (argc==1)
    {
#ifdef USE_MPI
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

#ifdef USE_MPI
    if (rank == 0)
        std::cout<<" para_umr_perf: ";
#else
    std::cout<<"umr_perf:";
#endif

    ErrorCode result;
    int deg[3] = {2,2,2};
    int len = sizeof(deg) / sizeof(int);
    result = umr_perf_test(filename, deg, len);
    handle_error_code(result, number_tests_failed, number_tests_successful);
    std::cout<<"\n";

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return number_tests_failed;
}
