// tests dual construction code
 
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MeshTopoUtil.hpp"
#include "DualTool.hpp"
#include "WriteGMV.hpp"
#include <iostream>
#include <string>

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

    // see if it's all-hex or all-quad, and construct hex/quad dual if so
  int num_hex, num_quad, num_3d, num_2d;
  result = gMB->get_number_entities_by_dimension(0, 2, num_2d);
  if (MB_SUCCESS != result) return 0;
  result = gMB->get_number_entities_by_dimension(0, 3, num_3d);
  if (MB_SUCCESS != result) return 0;
  result = gMB->get_number_entities_by_type(0, MBQUAD, num_quad);
  if (MB_SUCCESS != result) return 0;
  result = gMB->get_number_entities_by_type(0, MBHEX, num_hex);
  if (MB_SUCCESS != result) return 0;
  
    // construct the dual for this mesh
  DualTool dt(gMB);
  if (num_quad == num_2d && num_hex == num_3d)
    result = dt.construct_hex_dual();
  else
    result = dt.construct_dual();
  if (MB_SUCCESS != result) {
    std::cout << "Problems constructing dual." << std::endl;
    return 0;
  }
  
    // print information about the dual
  MBRange dual_cells, dual_faces;
  result = dt.get_dual_entities(2, dual_faces);
  if (MB_SUCCESS != result)
    std::cout << "Problem getting dual faces." << std::endl;
  else
    std::cout << "Found " << dual_faces.size() << "/" <<  num_edges << " dual faces." 
              << std::endl;
    
  result = dt.get_dual_entities(3, dual_cells);
  if (MB_SUCCESS != result)
    std::cout << "Problem getting dual cells." << std::endl;
  else
    std::cout << "Found " << dual_cells.size() << "/" << all_verts.size() << " dual cells." 
              << std::endl;

    // print information about dual hyperplanes, if any
  MBTag hp_tag;
  MBRange hp_sets;
  if (num_2d == num_quad) {
      // get sets with the right tag
    result = gMB->tag_get_handle(DualTool::DUAL_CURVE_TAG_NAME, hp_tag);
    if (MB_SUCCESS == result) {
      result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, &hp_tag, NULL, 1,
                                                 hp_sets);
      if (MB_SUCCESS == result && !hp_sets.empty())
        std::cout << "Found " << hp_sets.size() << " 1d dual hyperplanes (chords)." 
                  << std::endl;
    }
  }
  
  if (num_3d == num_hex) {
      // get sets with the right tag
    result = gMB->tag_get_handle(DualTool::DUAL_SURFACE_TAG_NAME, hp_tag);
    if (MB_SUCCESS == result) {
      hp_sets.clear();
      result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, &hp_tag, NULL, 1,
                                                 hp_sets);
      if (MB_SUCCESS == result && !hp_sets.empty())
        std::cout << "Found " << hp_sets.size() << " 2d dual hyperplanes (sheets)." 
                  << std::endl;
    }
  }

    // write GMV file
  WriteGMV *wgmv = new WriteGMV(gMB);
  wgmv->write_file(argv[1], 0, 3);
  delete wgmv;
  
  delete gMB;
}
