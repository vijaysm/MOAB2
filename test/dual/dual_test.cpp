// tests dual construction code
 
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MeshTopoUtil.hpp"
#include "DualTool.hpp"
#include "iostream"

MBInterface *gMB;

int main(int argc, char* argv[])
{
  if (argc <= 1)
  {
    std::cout << "Usage: dual_test <mesh_file_name>" << std::endl;
    return 0;
  }

  gMB = new MBCore();

    // read the mesh file
  MBErrorCode result = gMB->load_mesh(argv[1]);
  if (MB_SUCCESS != result) {
    std::cout << "Problems reading file " << argv[1] << "." << std::endl;
    return 0;
  }

    // make sure aentities are created
  MBRange all_verts;
  result = gMB->get_entities_by_dimension(0, 0, all_verts);
  if (MB_SUCCESS != result) 
    std::cout << "Problem getting vertices." << std::endl;
  
  MeshTopoUtil(gMB).construct_aentities(all_verts);
  
    // get counts before constructing dual
  int num_edges;
  result = gMB->get_number_entities_by_dimension(0, 1, num_edges);
  if (MB_SUCCESS != result) 
    std::cout << "Problem getting number of edges." << std::endl;
  
  
    // construct the dual for this mesh
  DualTool dt(gMB);
  result = dt.construct_dual();
  if (MB_SUCCESS != result) {
    std::cout << "Problems constructing dual." << std::endl;
    return 0;
  }
  
    // print information about the dual
  MBRange dual_cells, dual_faces;
  result = dt.get_dual_faces(dual_faces);
  if (MB_SUCCESS != result)
    std::cout << "Problem getting dual faces." << std::endl;
  else
    std::cout << "Found " << dual_faces.size() << "/" <<  num_edges << " dual faces." 
              << std::endl;
    
  result = dt.get_dual_cells(dual_cells);
  if (MB_SUCCESS != result)
    std::cout << "Problem getting dual cells." << std::endl;
  else
    std::cout << "Found " << dual_cells.size() << "/" << all_verts.size() << " dual cells." 
              << std::endl;

  delete gMB;
}
