/*This function profiles the performance of the AHF datastructure */
#include <iostream>
#include <assert.h>
#include <time.h>
#include <vector>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "../TestUtil.hpp"
#include "moab/CpuTimer.hpp"

#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#include "moab/FileOptions.hpp"
#include "MBTagConventions.hpp"
#include "moab_mpi.h"
#endif

using namespace moab;

#ifdef MOAB_HAVE_MPI
std::string read_options;
#endif

int number_tests_successful = 0;
int number_tests_failed = 0;

struct query_time
{
  double ds_construction;
  double vertex_to_edges_total;
  double vertex_to_edges_avg;
  double edge_to_edges_total;
  double edge_to_edges_avg;
  double vertex_to_faces_total;
  double vertex_to_faces_avg;
  double edge_to_faces_total;
  double edge_to_faces_avg;
  double face_to_faces_total;
  double face_to_faces_avg;
  double face_to_edges_total;
  double face_to_edges_avg;
  double vertex_to_cells_total;
  double vertex_to_cells_avg;
  double edge_to_cells_total;
  double edge_to_cells_avg;
  double face_to_cells_total;
  double face_to_cells_avg;
  double cell_to_cells_total;
  double cell_to_cells_avg;
  double cell_to_edges_total;
  double cell_to_edges_avg;
  double cell_to_faces_total;
  double cell_to_faces_avg;
};

struct mesh_mem
{
   unsigned long long total_storage[2];
   unsigned long long amortized_total_storage[2];
   unsigned long long entity_storage[2];
   unsigned long long amortized_entity_storage[2];
   unsigned long long adjacency_storage[2];
   unsigned long long amortized_adjacency_storage[2];
   unsigned long long tag_storage[2];
   unsigned long long amortized_tag_storage[2];
};

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

ErrorCode adj_perf(const char* filename)
{
  ErrorCode error;
  Core moab;
  Interface* mbImpl = &moab;
  MeshTopoUtil mtu(mbImpl);

  struct query_time qtime;
  struct mesh_mem qmem;

#ifdef MOAB_HAVE_MPI
    int procs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (procs > 1){
    read_options = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

    error = mbImpl->load_file(filename, 0, read_options.c_str());
    CHECK_ERR(error);
    }
    else if (procs == 1) {
#endif
    error = mbImpl->load_file(filename);
    CHECK_ERR(error);
#ifdef MOAB_HAVE_MPI
    }
#endif

    //Storage Costs before any call to adjacencies
    unsigned long long sTotS, sTAS, sES, sAES, sAS, sAAS, sTS, sATS;
    sTotS = sTAS = sES = sAES = sAS = sAAS = sTS = sATS = 0;
    mbImpl->estimated_memory_use(NULL, 0, &sTotS, &sTAS, &sES, &sAES, &sAS, &sAAS, NULL, 0, &sTS, &sATS);

    qmem.total_storage[0] = sTotS;
    qmem.amortized_total_storage[0] = sTAS;
    qmem.entity_storage[0] = sES;
    qmem.amortized_entity_storage[0] = sAES;
    qmem.adjacency_storage[0] = sAS;
    qmem.amortized_adjacency_storage[0] = sAAS;
    qmem.tag_storage[0] = sTS;
    qmem.amortized_tag_storage[0] = sATS;


    //Create ranges for handles of explicit elements of the mixed mesh
    Range verts, edges, faces, cells;
    error = mbImpl->get_entities_by_dimension( 0, 0, verts);
    error = mbImpl->get_entities_by_dimension( 0, 1, edges);
    error = mbImpl->get_entities_by_dimension( 0, 2, faces);
    error = mbImpl->get_entities_by_dimension( 0, 3, cells);

  int nverts = verts.size();
  int nedges = edges.size();
  int nfaces = faces.size();
  int ncells = cells.size();

  std::cout<<"MESH SIZE :: "<<"NV = "<<nverts<<", NE = "<<nedges<<", NF = "<<nfaces<<", NC = "<<ncells<<std::endl;

  CpuTimer *mt = new CpuTimer;
  double time_start, time_avg, time_total;

  //Perform queries
  std::vector<EntityHandle> adjents;
  Range ngbents;

  // This call should create all the necessary ahf maps or adjacency lists
  time_start = mt->time_elapsed();
  error = mbImpl->get_adjacencies( &*verts.begin(), 1, 1, false, adjents );
  time_total = mt->time_elapsed() - time_start;
  qtime.ds_construction = time_total;

  //1D Queries

  std::cout<<"1D QUERIES Start"<<std::endl;

  //IQ1: For every vertex, obtain incident edges
  time_start = mt->time_elapsed();
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)verts.size();

  qtime.vertex_to_edges_total = time_total;
  qtime.vertex_to_edges_avg = time_avg;

  //NQ1:  For every edge, obtain neighbor edges  
#ifdef MOAB_HAVE_AHF
  time_start = mt->time_elapsed();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {    
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)edges.size();  
#else
  error = mtu.get_bridge_adjacencies( *edges.begin(), 0, 1, ngbents);
  time_start = mt->time_elapsed();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
      ngbents.clear();
      error = mtu.get_bridge_adjacencies( *i, 0, 1, ngbents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)edges.size();
#endif

  qtime.edge_to_edges_total = time_total;
  qtime.edge_to_edges_avg = time_avg;

  std::cout<<"1D QUERIES End"<<std::endl;

  // 2D Queries

  std::cout<<"2D QUERIES Start"<<std::endl;

  //IQ21: For every vertex, obtain incident faces
  time_start = mt->time_elapsed();
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)verts.size();

  qtime.vertex_to_faces_total = time_total;
  qtime.vertex_to_faces_avg = time_avg;

  //IQ22: For every edge, obtain incident faces
  time_start = mt->time_elapsed();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {   
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)edges.size();

  qtime.edge_to_faces_total = time_total;
  qtime.edge_to_faces_avg = time_avg;

  //NQ2: For every face, obtain neighbor faces 
#ifdef MOAB_HAVE_AHF
  time_start = mt->time_elapsed();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
      adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)faces.size();
#else
  error = mtu.get_bridge_adjacencies( *faces.begin(), 1, 2, ngbents);
  time_start = mt->time_elapsed();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
      ngbents.clear();
      error = mtu.get_bridge_adjacencies( *i, 1, 2, ngbents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)faces.size();
#endif

  qtime.face_to_faces_total = time_total;
  qtime.face_to_faces_avg = time_avg;

  //DQ2: For every face, obtain its edges
  time_start = mt->time_elapsed();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)faces.size();

  qtime.face_to_edges_total = time_total;
  qtime.face_to_edges_avg = time_avg;

  std::cout<<"2D QUERIES End"<<std::endl;

  // 3D Queries

  std::cout<<"3D QUERIES Start "<<std::endl;

  //IQ31: For every vertex, obtain incident cells
  time_start = mt->time_elapsed();
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
      adjents.clear();
      error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)verts.size();

  qtime.vertex_to_cells_total = time_total;
  qtime.vertex_to_cells_avg = time_avg;

  // IQ 32: For every edge, obtain incident cells
  time_start = mt->time_elapsed();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
      adjents.clear();
      error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)edges.size();

  qtime.edge_to_cells_total = time_total;
  qtime.edge_to_cells_avg = time_avg;

  //IQ32: For every face, obtain incident cells
  time_start = mt->time_elapsed();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
      adjents.clear();
      error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)faces.size();

  qtime.face_to_cells_total = time_total;
  qtime.face_to_cells_avg = time_avg;

  //NQ3: For every cell, obtain neighbor cells
#ifdef MOAB_HAVE_AHF
  time_start = mt->time_elapsed();
  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {   
      adjents.clear();
      error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)cells.size();
#else
  error = mtu.get_bridge_adjacencies( *cells.begin(), 2, 3, ngbents);
  time_start = mt->time_elapsed();
  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
    ngbents.clear();
    error = mtu.get_bridge_adjacencies( *i, 2, 3, ngbents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)cells.size();
#endif

  qtime.cell_to_cells_total = time_total;
  qtime.cell_to_cells_avg = time_avg;

  //DQ31: For every cell, obtain its edges
  time_start = mt->time_elapsed();
  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)cells.size();

  qtime.cell_to_edges_total = time_total;
  qtime.cell_to_edges_avg = time_avg;

  //DQ32: For every cell, obtain its faces
  time_start = mt->time_elapsed();
  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_total = mt->time_elapsed()-time_start;
  time_avg = time_total/(double)cells.size();

  qtime.cell_to_faces_total = time_total;
  qtime.cell_to_faces_avg = time_avg;

   std::cout<<"3D QUERIES End"<<std::endl;

  //Storage Costs after calling ahf deinitialize
  unsigned long long eTotS, eTAS, eES, eAES, eAS, eAAS, eTS, eATS;
  eTotS = eTAS = eES = eAES = eAS = eAAS = eTS = eATS = 0;
  mbImpl->estimated_memory_use(NULL, 0, &eTotS, &eTAS, &eES, &eAES, &eAS, &eAAS, NULL, 0, &eTS, &eATS);

  qmem.total_storage[1] = eTotS;
  qmem.amortized_total_storage[1] = eTAS;
  qmem.entity_storage[1] = eES;
  qmem.amortized_entity_storage[1] = eAES;
  qmem.adjacency_storage[1] = eAS;
  qmem.amortized_adjacency_storage[1] = eAAS;
  qmem.tag_storage[1] = eTS;
  qmem.amortized_tag_storage[1] = eATS;

  //Print times
  std::cout<<std::endl;
  std::cout<<" Data Structure Construction Time = "<<qtime.ds_construction<<" Secs"<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Query times in Seconds"<<std::endl;
#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Vertex -> Edges :: MOAB_AHF: Average time =  "<< qtime.vertex_to_edges_avg ;
  std::cout<<", Total time = "<<qtime.vertex_to_edges_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Vertex -> Edges :: MOAB: Average time =  "<<  qtime.vertex_to_edges_avg;
  std::cout<<", Total time = "<<qtime.vertex_to_edges_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Edge -> Edges :: MOAB_AHF: Average time = "<<qtime.edge_to_edges_avg;
  std::cout<<", Total time = "<<qtime.edge_to_edges_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Edge -> Edges :: MOAB: Average time = "<<qtime.edge_to_edges_avg;
  std::cout<<", Total time = "<<qtime.edge_to_edges_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Vertex -> Faces :: MOAB_AHF: Average time =  "<<qtime.vertex_to_faces_avg;
  std::cout<<", Total time = "<<qtime.vertex_to_faces_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Vertex -> Faces :: MOAB: Average time =  "<<qtime.vertex_to_faces_avg;
  std::cout<<", Total time = "<<qtime.vertex_to_faces_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Edge -> Faces :: MOAB_AHF: Average time =  "<<qtime.edge_to_faces_avg;
  std::cout<<", Total time = "<<qtime.edge_to_faces_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Edge -> Faces :: MOAB: Average time =  "<<qtime.edge_to_faces_avg;
  std::cout<<", Total time = "<<qtime.edge_to_faces_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Face -> Faces :: MOAB_AHF: Average time = "<< qtime.face_to_faces_avg;
  std::cout<<", Total time = "<<qtime.face_to_faces_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Face -> Faces :: MOAB: Average time = "<< qtime.face_to_faces_avg;
  std::cout<<", Total time = "<<qtime.face_to_faces_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Face -> Edges :: MOAB_AHF: Average time =  "<<qtime.face_to_edges_avg;
  std::cout<<", Total time = "<<qtime.face_to_edges_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Face -> Edges :: MOAB: Average time =  "<<qtime.face_to_edges_avg;
  std::cout<<", Total time = "<<qtime.face_to_edges_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Vertex -> Cells :: MOAB_AHF: Average time =  "<<qtime.vertex_to_cells_avg;
  std::cout<<", Total time = "<<qtime.vertex_to_cells_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Vertex -> Cells :: MOAB: Average time =  "<<qtime.vertex_to_cells_avg;
  std::cout<<", Total time = "<<qtime.vertex_to_cells_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Edge -> Cells :: MOAB_AHF: Average time =  "<<qtime.edge_to_cells_avg;
  std::cout<<", Total time = "<<qtime.edge_to_cells_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Edge -> Cells :: MOAB: Average time =  "<<qtime.edge_to_cells_avg;
  std::cout<<", Total time = "<<qtime.edge_to_cells_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Face -> Cells :: MOAB_AHF: Average time =  "<<qtime.face_to_cells_avg;
  std::cout<<", Total time = "<<qtime.face_to_cells_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Face -> Cells :: MOAB: Average time =  "<<qtime.face_to_cells_avg;
  std::cout<<", Total time = "<<qtime.face_to_cells_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Cell -> Cells :: MOAB_AHF: Average time =  "<< qtime.cell_to_cells_avg;
  std::cout<<", Total time = "<<qtime.cell_to_cells_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Cell -> Cells :: MOAB: Average time =  "<< qtime.cell_to_cells_avg;
  std::cout<<", Total time = "<<qtime.cell_to_cells_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Cell -> Edges :: MOAB_AHF: Average time =  "<<qtime.cell_to_edges_avg;
  std::cout<<", Total time = "<<qtime.cell_to_edges_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Cell -> Edges :: MOAB: Average time =  "<<qtime.cell_to_edges_avg;
  std::cout<<", Total time = "<<qtime.cell_to_edges_total<<std::endl;
  std::cout<<std::endl;
#endif

#ifdef MOAB_HAVE_AHF
  std::cout<<"QUERY: Cell -> Faces :: MOAB_AHF: Average time =  "<<qtime.cell_to_faces_avg;
  std::cout<<", Total time = "<<qtime.cell_to_faces_total<<std::endl;
  std::cout<<std::endl;
#else
  std::cout<<"QUERY: Cell -> Faces :: MOAB: Average time =  "<<qtime.cell_to_faces_avg;
  std::cout<<", Total time = "<<qtime.cell_to_faces_total<<std::endl;
  std::cout<<std::endl;
#endif


  //Print Storage
  std::cout<<std::endl;
  for (int i=0; i<2; i++){
      if (i==0)
        std::cout<<"STORAGE BEFORE CALLING ADJACENCIES"<<std::endl;
      else
        std::cout<<"STORAGE AFTER CALLING ADJACENCIES"<<std::endl;

      std::cout<<"Total storage = "<<qmem.total_storage[i]<<std::endl;
      std::cout<<"Total amortized storage = "<<qmem.amortized_total_storage[i]<<std::endl;
      std::cout<<std::endl;

      std::cout<<"Entity storage = "<<qmem.entity_storage[i]<<std::endl;
      std::cout<<"Amortized entity storage = "<<qmem.amortized_entity_storage[i]<<std::endl;
      std::cout<<std::endl;

      std::cout<<"Adjacency storage = "<<qmem.adjacency_storage[i]<<std::endl;
      std::cout<<"Amortized adjacency storage = "<<qmem.amortized_adjacency_storage[i]<<std::endl;
      std::cout<<std::endl;

      std::cout<<"Tag storage = "<<qmem.tag_storage[i]<<std::endl;
      std::cout<<"Amortized tag storage = "<<qmem.amortized_tag_storage[i]<<std::endl;
      std::cout<<std::endl;
    }

  double total_time = qtime.vertex_to_edges_total+qtime.edge_to_edges_total+qtime.vertex_to_faces_total+qtime.edge_to_faces_total+qtime.face_to_faces_total+qtime.face_to_edges_total+qtime.vertex_to_cells_total+qtime.edge_to_cells_total+qtime.face_to_cells_total+qtime.cell_to_cells_total+qtime.cell_to_edges_total+qtime.cell_to_faces_total;

  //Print values in a line to aid data copying later
  std::cout<<qtime.ds_construction<<"  "<<total_time<<"  "<<qmem.entity_storage[1]<<"  "<<qmem.adjacency_storage[1]<<"  "<<qtime.vertex_to_edges_avg<<"  "<<qtime.edge_to_edges_avg<<"  "<<qtime.vertex_to_faces_avg<<"  "<<qtime.edge_to_faces_avg<<"  "<<qtime.face_to_faces_avg<<"  "<<qtime.face_to_edges_avg<<"  "<<qtime.vertex_to_cells_avg<<"  "<<qtime.edge_to_cells_avg<<"  "<<qtime.face_to_cells_avg<<"  "<<qtime.cell_to_cells_avg<<"  "<<qtime.cell_to_edges_avg<<"  "<<qtime.cell_to_faces_avg<<std::endl;


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
        std::cout<<" para_adj_perf: ";
#else
    std::cout<<"adj_perf:";
#endif

    result = adj_perf(filename);
    handle_error_code(result, number_tests_failed, number_tests_successful);
    std::cout<<"\n";

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return number_tests_failed;
}



