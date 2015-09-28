// simple example construct obb tree and ray-tracing the tree
// it reads triangle mesh, construct obb tree and get intersection distances

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include <iostream>
#include <math.h>

int main(int argc, char **argv) {
  if (1 == argc) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }

  // instantiate & load a mesh from a file
  moab::Core *mb = new moab::Core();
  moab::ErrorCode rval = mb->load_mesh(argv[1]);
  if (rval != moab::MB_SUCCESS) {
    std::cerr << "Couldn't load mesh." << std::endl;
    delete mb;
    return 1;
  }

  // get all triangles
  moab::EntityHandle tree_root;
  moab::Range tris;
  //moab::OrientedBoxTreeTool::Settings settings;

  rval = mb->get_entities_by_type(0, moab::MBTRI, tris);
  if (rval != moab::MB_SUCCESS) {
    std::cerr << "Couldn't get triangles." << std::endl;
    delete mb;
    return 1;
  }

  // build OBB trees for all triangles
  moab::OrientedBoxTreeTool tool(mb);
  //rval = tool.build(tris, tree_root, &settings);
  rval = tool.build(tris, tree_root);
  if (rval != moab::MB_SUCCESS) {
    std::cerr << "Could'nt build tree." << std::endl;
    delete mb;
    return 1;
  }
  
  // build box
  double box_center[3], box_axis1[3], box_axis2[3], box_axis3[3], pnt_start[3], ray_length;
  rval = tool.box(tree_root, box_center, box_axis1, box_axis2, box_axis3);
  if (rval != moab::MB_SUCCESS) {
    std::cerr << "Couldn't get box for tree root set.";
    delete mb;
    return 1;
  }

  ray_length = 2.*sqrt(box_axis1[0]*box_axis1[0] + box_axis1[1]*box_axis1[1] +
		       box_axis1[2]*box_axis1[2]);
  
  // do ray-tracing from box center side to x direction
  std::vector<double> intersections;
  std::vector<moab::EntityHandle> intersection_facets;

  for (int i = 0; i < 3; i++)
    pnt_start[i] = box_center[i] - box_axis1[i];

  if (ray_length > 0) { // normalize ray direction
    for (int j = 0; j < 3; j++)
      box_axis1[j] = 2 * box_axis1[j] / ray_length;
  }
  rval = tool.ray_intersect_triangles(intersections, intersection_facets, 
				      tree_root, 10e-12, pnt_start, box_axis1,
				      &ray_length);
  if (rval != moab::MB_SUCCESS) {
    std::cerr << "Couldn't ray tracing.";
    delete mb;
    return 1;
  }
  
  std::cout << "ray start point: " << pnt_start[0] << " "
	    << pnt_start[1] << " " << pnt_start[2] << std::endl;
  std::cout << " ray direction: " << box_axis1[0] << " " << box_axis1[1] << " " << box_axis1[2] << "\n";
  std::cout << "# of intersections : " << intersections.size() << std::endl;
  std::cout << "intersection distances are on";
  for (unsigned int i = 0; i < intersections.size(); i++)
    std::cout << " " << intersections[i];
  std::cout << " of ray length " << ray_length << std::endl;

  delete mb;

  return 0;
}
