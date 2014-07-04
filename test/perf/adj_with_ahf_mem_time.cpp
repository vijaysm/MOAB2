/*This function tests the AHF datastructures on CST meshes*/
#include <iostream>
#include <assert.h>
#include <time.h>
#include <vector>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include <sys/time.h>

using namespace moab;

#ifdef MESHDIR
std::string TestDir(STRINGIFY(MESHDIR));
#else
std::string TestDir(".");
#endif

std::string filename;

double wtime() {
  double y = -1;
  struct timeval cur_time;  
  gettimeofday(&cur_time, NULL);  
  y = (double)(cur_time.tv_sec) + (double)(cur_time.tv_usec)*1.e-6;  
  return (y);
}

int main(int argc, char **argv)
{
  // Read the input mesh
    filename = TestDir + "/hexes_mixed.vtk";

    if (argc==1)
        std::cout<<"Using default input file:"<<filename<<std::endl;
    else if (argc==2)
        filename = argv[1];
    else {
            std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
            return 1;
    }

  ErrorCode error;
  Core moab;
  Interface* mbImpl = &moab;
  MeshTopoUtil mtu(mbImpl);

  error = mbImpl->load_file( filename.c_str());
  if (MB_SUCCESS != error) {
    std::cerr << filename <<": failed to load file." << std::endl;
    return error;
  }

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

  std::cout<<"nverts = "<<nverts<<", nedges = "<<nedges<<", nfaces = "<<nfaces<<", ncells = "<<ncells<<std::endl;


  //Storage Costs before calling ahf functionalities
  std::cout<<std::endl;
  std::cout<<"STORAGE BEFORE CALLING ADJACENCIES"<<std::endl;
  unsigned long sTotS, sTAS, sES, sAES, sAS, sAAS, sTS, sATS;
  sTotS = sTAS = sES = sAES = sAS = sAAS = sTS = sATS = 0;
  mbImpl->estimated_memory_use(NULL, 0, &sTotS, &sTAS, &sES, &sAES, &sAS, &sAAS, NULL, 0, &sTS, &sATS);
  std::cout<<std::endl;
  std::cout<<"Total storage = "<<sTotS<<std::endl;
  std::cout<<"Total amortized storage = "<<sTAS<<std::endl;
  std::cout<<"Entity storage = "<<sES<<std::endl;
  std::cout<<"Amortized entity storage = "<<sAES<<std::endl;
  std::cout<<"Adjacency storage = "<<sAS<<std::endl;
  std::cout<<"Amortized adjacency storage = "<<sAAS<<std::endl;
  std::cout<<"Tag storage = "<<sTS<<std::endl;
  std::cout<<"Amortized tag storage = "<<sATS<<std::endl;
  std::cout<<std::endl;


  double time_start, time_avg;

  //Perform queries
  std::vector<EntityHandle> adjents;

  // This call should create all the necessary ahf maps
  error = mbImpl->get_adjacencies( &*verts.begin(), 1, 1, false, adjents );

  //1D Queries //
  //IQ1: For every vertex, obtain incident edges
  std::cout<<"1D QUERIES"<<std::endl;
  time_start = wtime();
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {    
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)verts.size();
  std::cout<<"QUERY: Vertex -> Edges :: MOAB_AHF: Average time =  "<< time_avg<<" secs" <<std::endl;

  //NQ1:  For every edge, obtain neighbor edges  
  time_start = wtime();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {    
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)edges.size();
  std::cout<<"QUERY: Edge -> Edges :: MOAB_AHF: Average time = "<<time_avg<<" secs" << std::endl;

  // 2D Queries
  std::cout<<"2D QUERIES"<<std::endl;
  //IQ21: For every vertex, obtain incident faces
  time_start = wtime();
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)edges.size();
  std::cout<<"QUERY: Vertex -> Faces :: MOAB_AHF: Average time =  "<<time_avg<<" secs" <<std::endl;


  //IQ22: For every edge, obtain incident faces
  time_start = wtime();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {   
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)edges.size();
  std::cout<<"QUERY: Edge -> Faces :: MOAB_AHF: Average time =  "<<time_avg<<" secs" <<std::endl;


  //NQ2: For every face, obtain neighbor faces 
  time_start = wtime();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {   
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)faces.size();
  std::cout<<"QUERY: Face -> Faces :: MOAB_AHF: Average time = "<< time_avg<<" secs" <<std::endl;


  // 3D Queries
  std::cout<<"3D QUERIES"<<std::endl;
  //IQ31: For every vertex, obtain incident cells
  time_start = wtime();
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)edges.size();
  std::cout<<"QUERY: Vertex -> Cells :: MOAB_AHF: Average time =  "<<time_avg <<" secs"<<std::endl;


  // IQ 32: For every edge, obtain incident cells
  time_start = wtime();
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)edges.size();
  std::cout<<"QUERY: Edge -> Cells :: MOAB_AHF: Average time =  "<<time_avg <<" secs"<<std::endl;


  //IQ32: For every face, obtain incident cells
  time_start = wtime();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)faces.size();
  std::cout<<"QUERY: Face -> Cells :: MOAB_AHF: Average time =  "<<time_avg <<" secs"<<std::endl;


  //NQ3: For every cell, obtain neighbor cells
  time_start = wtime();
  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {   
    adjents.clear();
    error = mbImpl->get_adjacencies(&*i, 1, 3, false, adjents);
  }
  time_avg = (wtime()-time_start)/(double)cells.size();
  std::cout<<"QUERY: Cell -> Cells :: MOAB_AHF: Average time =  "<< time_avg <<" secs" << std::endl;


  //Storage Costs after calling ahf deinitialize
  std::cout<<std::endl;
  std::cout<<"STORAGE AFTER CALLING ADJACENCIES"<<std::endl;
  unsigned long eTotS, eTAS, eES, eAES, eAS, eAAS, eTS, eATS;
  eTotS = eTAS = eES = eAES = eAS = eAAS = eTS = eATS = 0;
  mbImpl->estimated_memory_use(NULL, 0, &eTotS, &eTAS, &eES, &eAES, &eAS, &eAAS, NULL, 0, &eTS, &eATS);
  std::cout<<std::endl;
  std::cout<<"Total storage = "<<eTotS<<std::endl;
  std::cout<<"Total amortized storage = "<<eTAS<<std::endl;
  std::cout<<"Entity storage = "<<eES<<std::endl;
  std::cout<<"Amortized entity storage = "<<eAES<<std::endl;
  std::cout<<"Adjacency storage = "<<eAS<<std::endl;
  std::cout<<"Amortized adjacency storage = "<<eAAS<<std::endl;
  std::cout<<"Tag storage = "<<eTS<<std::endl;
  std::cout<<"Amortized tag storage = "<<eATS<<std::endl;
  std::cout<<std::endl;

  return 0;
}

