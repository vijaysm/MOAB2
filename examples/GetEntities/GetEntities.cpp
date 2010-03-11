#include "moab/MBCore.hpp"
#include "moab/MBRange.hpp"
#include <iostream>

int main(int, char **argv) {
    // instantiate & load a mesh from a file
  MBCore *mb = new MBCore();
  MBErrorCode rval = mb->load_mesh(argv[1]);

  MBRange ents;

    // iterate over dimensions
  for (int d = 0; d <= 3; d++) {
    rval = mb->get_entities_by_dimension(0, d, ents);
    for (MBRange::iterator it = ents.begin(); it != ents.end(); it++) {
      MBEntityHandle ent = *it;
      std::cout << "Found d=" << d << " entity " 
                << mb->id_from_handle(ent) << "." << std::endl;
    }
  }
  return 0;
}
