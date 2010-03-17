#include "moab/MBCore.hpp"
#include "moab/MBRange.hpp"
#include <iostream>

int main(int argc, char **argv) {
  if (1 == argc) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }
  
    // instantiate & load a file
  MBInterface *mb = new MBCore();
  MBErrorCode rval = mb->load_file(argv[1]);

    // get the geometric topology tag handle
  MBTag geom_tag;
  rval = mb->tag_get_handle("GEOM_DIMENSION", geom_tag);

    // traverse the model, from dimension 3 downward
  MBRange psets, chsets;
  int dim;
  void *dim_ptr = &dim;
  for (dim = 3; dim >= 0; dim--) {
      // get parents at this dimension
    psets.clear();
    rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                            &geom_tag, &dim_ptr, 1, 
                                            psets, 1, false);

      // for each parent, get children and do something with them
    MBRange::iterator par_it;
    for (par_it = psets.begin(); par_it != psets.end(); par_it++) {
        // get the children and put in child set list
      chsets.clear();
      rval = mb ->get_child_meshsets(*par_it, chsets);

        // print # children
      std::cout << "d=" << dim << " entity has " << chsets.size() 
                << " children." << std::endl;
    }
  } 
}
