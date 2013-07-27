/** @example SetsNTags.cpp
 * Description: Get the sets representing materials and Dirichlet/Neumann boundary conditions and list their contents.\n
 * This example shows how to get entity sets, and tags on those sets.
 *
 * To run: ./SetsNTags [meshfile]\n
 * (default values can run if users don't specify a mesh file)
 */

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"

#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

string test_file_name = string(MESH_DIR) + string("/1hex.g");

// tag names for these conventional tags come from MBTagConventions.hpp
const char *tag_nms[] = {MATERIAL_SET_TAG_NAME, DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME};

int main(int argc, char **argv) {
    // get the material set tag handle
  Tag mtag;
  ErrorCode rval;
  Range sets, set_ents;

    // instantiate & load a file
  Interface *mb = new Core();

    // need option handling here for input filename
  if (argc > 1){
    //user has input a mesh file
    test_file_name = argv[1];
  }  

  rval = mb->load_file(test_file_name.c_str());
  if (MB_SUCCESS != rval) return 1;

    // loop over set types
  for (int i = 0; i < 3; i++) {
      // get the tag handle for this tag name; tag should already exist (it was created during file read)
    rval = mb->tag_get_handle(tag_nms[i], 1, MB_TYPE_INTEGER, mtag);
    if (MB_SUCCESS != rval) return 1;

      // get all the sets having that tag (with any value for that tag)
    sets.clear();
    rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &mtag, NULL, 1, sets);
    if (MB_SUCCESS != rval) return 1;

      // iterate over each set, getting the entities and printing them
    Range::iterator set_it;
    for (set_it = sets.begin(); set_it != sets.end(); set_it++)  {
        // get the id for this set
      int set_id;
      rval = mb->tag_get_data(mtag, &(*set_it), 1, &set_id);
      if (MB_SUCCESS != rval) return 1;

        // get the entities in the set, recursively
      rval = mb->get_entities_by_handle(*set_it, set_ents, true);
      if (MB_SUCCESS != rval) return 1;

      cout << tag_nms[i] << " " << set_id << " has " 
                << set_ents.size() << " entities:" << endl;
      set_ents.print("   ");
      set_ents.clear();
    }
  }

  delete mb;
}
