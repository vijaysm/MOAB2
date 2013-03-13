/** \brief This test shows how to perform local point-in-element searches with MOAB's new tree searching functionality.  
 *
 * MOAB's SpatialLocator functionality performs point-in-element searches over a local or parallel mesh.
 * SpatialLocator is flexible as to what kind of tree is used and what kind of element basis functions are 
 * used to localize elements and interpolate local fields.
 */

#include <iostream>
#include <cstdlib>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

using namespace moab;

#define ERR(s) if (MB_SUCCESS != rval) \
    {std::string str;mb->get_last_error(str); std::cerr << s << str << std::endl; return 1;}

int main(int argc, char **argv) {

  int num_queries = 1000000;
  
  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << "<filename> [num_queries]" << std::endl;
    return 0;
  }
  else if (argc == 3) sscanf(argv[2], "%d", num_queries);

    // instantiate & load a file
  moab::Interface *mb = new moab::Core();

  Error err;
  ErrorCode rval = mb->query_interface(err);
  if (MB_SUCCESS != rval) return 1;
  
    // load the file
  rval = mb->load_file(argv[argc-1]); ERR("Error loading file");
  
    // get all 3d elements in the file
  Range elems;
  rval = mb->get_entities_by_dimension(0, 3, elems); ERR("Error getting 3d elements");
  
    // create a tree to use for the location service
  AdaptiveKDTree tree(mb);

    // specify an evaluator based on linear hexes
  ElemEvaluator el_eval(mb, LinearHex::get_eval_set());

    // build the SpatialLocator
  SpatialLocator sl(mb, elems, &tree, &el_eval);
  
    // get the box extents
  CartVect box_min, box_max, box_extents, pos;
  rval = sl.get_bounding_box(box_min, box_max); ERR("Problem getting tree bounding box");
  box_extents = box_max - box_min;
  
    // query at random places in the tree
  EntityHandle elem;
  CartVect params;
  bool is_inside;
  int num_inside = 0;
  for (int i = 0; i < num_queries; i++) {
    pos = box_min + 
        CartVect(box_extents[0]*.01*(rand()%100), box_extents[1]*.01*(rand()%100), box_extents[2]*.01*(rand()%100));
    
    ErrorCode tmp_rval = sl.locate_point(pos, tol, in_elem, params, &is_inside);
    if (MB_SUCCESS != tmp_rval) rval = tmp_rval;
    if (is_inside) num_inside++;
  }
  
  std::cout << "Mesh contains " << elems.size() << " elements of type " 
            << CN::EntityTypeName(mb->type_from_handle(*elems.begin())) << std::endl;
  std::cout << "Bounding box min-max = (" << box_min[0] << "," << box_min[1] << "," << box_min[2] << ")-("
            << box_max[0] << "," << box_max[1] << "," << box_max[2] << ")" << std::endl;
  std::cout << "Queries inside box = " << num_inside << "/" << num_queries << " = " 
            << (double)(num_inside/num_queries) << "%" << std::endl;
}

    
  
  
  

