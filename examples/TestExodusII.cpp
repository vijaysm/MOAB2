/** @example TestExodusII.cpp
 *  TestExodusII example demonstrates how to retrieve material, dirichlet and neumann sets from an Exodus file
 * Sets are traversed to find out the number of entities contained in each set.
 * The entities contained in each set are retrieved, using a Boolean flag to indicate that any contained sets
 * should be traversed recursively to include non-set entities in the results. */

#include <iostream>

// Include header for MOAB instance and range
#include "moab/Core.hpp"

int main(int argc, char **argv) {

  // instantiate & load a file
  moab::Interface *mb = new moab::Core();

  // If no input is specified load ../MeshFiles/unittest/mbtest2.g
  const char* test_file_name =  std::string("../MeshFiles/unittest/mbtest2.g");

  // get the material set tag handle
  moab::Tag mtag;
  moab::ErrorCode rval;
  const char *tag_nms[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  moab::Range sets, set_ents;

  if (argc == 1) {
      std::cout << "Running default case, loading ../MeshFiles/unittest/mbtest2.g" << std::endl;
      std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
      rval = mb->load_file("../MeshFiles/unittest/mbtest2.g");
    }
  else{
      rval = mb->load_file(argv[argc-1]);
    }

  // loop over set types
  for (int i = 0; i < 3; i++) {
      rval = mb->tag_get_handle(tag_nms[i], 1, moab::MB_TYPE_INTEGER, mtag);
      if (moab::MB_SUCCESS != rval) return 1;

      // get all the sets of that type in the mesh
      sets.clear();
      rval = mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &mtag,
                                              NULL, 1, sets);
      if (moab::MB_SUCCESS != rval) return 1;

      // iterate over each set, getting entities
      moab::Range::iterator set_it;
      for (set_it = sets.begin(); set_it != sets.end(); set_it++)  {
          moab::EntityHandle this_set = *set_it;

          // get the id for this set
          int set_id;
          rval = mb->tag_get_data(mtag, &this_set, 1, &set_id);
          if (moab::MB_SUCCESS != rval) return 1;

          // get the entities in the set, recursively
          rval = mb->get_entities_by_handle(this_set, set_ents, true);
          if (moab::MB_SUCCESS != rval) return 1;

          std::cout << tag_nms[i] << " " << set_id << " has "
                    << set_ents.size() << " entities:" << std::endl;
          set_ents.print("   ");
          set_ents.clear();
        }
    }

  // do the same for all sets
  sets.clear();
  rval = mb->get_entities_by_type(0, moab::MBENTITYSET, sets);
  if (moab::MB_SUCCESS != rval) return 1;

  delete mb;
}
