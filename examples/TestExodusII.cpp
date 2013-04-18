/** @example TestExodusII.cpp
 * This example demonstrates how to retrieve material, dirichlet and neumann sets
 * from an ExodusII file. \n
 * Sets in MOAB contain entities and have a tag on them to give meaning to the entities.
 * Tag names: MATERIAL_SET", "DIRICHLET_SET" and "NEUMANN_SET" are reserved and
 * are associated with their corresponding entity sets.
 * Sets are traversed to find out the type number of entities contained in each set. \n
 *
 * <b>Steps in this example </b>:
 *    -# Instantiate MOAB
 *    -# Get input mesh file name and load it.
 *    -# Loop over the three sets: material, dirichlet and neumann
 *      -# Get TagHandle and EntitySet(corresponding to the TagHandle)
 *      -# Loop thru all the EntitySet's
 *        -# Get the set id and entities in this set
 *    -# Destroy the MOAB instance
 *
 *
 * <b> To compile: </b>
 *    make TestExodusII MOAB_DIR=<installdir> \n
 *
 * <b> To run: </b>
 *    -# TestExodusII <mesh-file> \n
 *    -# TestExodusII (This uses the default <mesh-file>: <MOAB_SRC_DIR>/MeshFiles/unittest/mbtest2.g)
 */
#include <iostream>

// Include header for MOAB instance and range
#include "moab/Core.hpp"

int main(int argc, char **argv) {

  // instantiate & load a file
  moab::Interface *mb = new moab::Core();

  // If no input is specified load ../MeshFiles/unittest/mbtest2.g
//  const char* test_file_name =  "../MeshFiles/unittest/mbtest2.g";

  // get the material set tag handle
  moab::Tag mtag;
  moab::ErrorCode rval;
  const char *tag_nms[] = {"MATERIAL_SET", "DIRICHLET_SET", "NEUMANN_SET"};
  moab::Range sets, set_ents;

  if (argc == 1) {
      std::cout << "Running default case, loading ../MeshFiles/unittest/mbtest2.g" << std::endl;
      std::cout << "Usage: " << argv[0] << " <filename>\n" << std::endl;
      rval = mb->load_file("../MeshFiles/unittest/mbtest2.g");
    }
  else{
      rval = mb->load_file(argv[argc-1]);
      std::cout << "Loaded mesh file: " << argv[argc-1] << std::endl;
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

          // print the entities contained in this set
          set_ents.print("   ");
          set_ents.clear();
        }
    }
  delete mb;
}
