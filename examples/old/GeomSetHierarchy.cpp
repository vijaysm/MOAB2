#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBCN.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"
#include <iostream>

const char *ent_names[] = {"Vertex", "Edge", "Face", "Region"};

int main(int argc, char **argv) {
  if (1 == argc) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }
  
    // instantiate & load a file
  moab::Interface *mb = new moab::Core();
  moab::ErrorCode rval = mb->load_file(argv[1]);

    // get the geometric topology tag handle
  moab::Tag geom_tag, gid_tag;
  rval = mb->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geom_tag);
  rval = mb->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, gid_tag);

    // traverse the model, from dimension 3 downward
  moab::Range psets, chsets;
  std::vector<moab::EntityHandle> sense_ents;
  std::vector<int> senses, pgids;
  int dim, pgid, chgid;
  void *dim_ptr = &dim;
  int sense;

  moab::GeomTopoTool gt(mb, true);
  
  for (dim = 3; dim >= 0; dim--) {
      // get parents at this dimension
    chsets.clear();
    rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                            &geom_tag, &dim_ptr, 1, 
                                            chsets, 1, false);

      // for each child, get parents and do something with them
    moab::Range::iterator ch_it, p_it;
    for (ch_it = chsets.begin(); ch_it != chsets.end(); ch_it++) {
        // get the children and put in child set list
      psets.clear();
      rval = mb ->get_parent_meshsets(*ch_it, psets);

      rval = mb->tag_get_data(gid_tag, &(*ch_it), 1, &chgid);
      
        // print # parents
      std::cout << ent_names[dim] << " " << chgid << " has " << psets.size() 
                << " parents." << std::endl;

      if (2 == dim) {
        for (p_it = psets.begin(); p_it != psets.end(); p_it++) {
          rval = mb->tag_get_data(gid_tag, &(*p_it), 1, &pgid);
          rval = gt.get_sense(*ch_it, *p_it, sense);
          if (moab::MB_SUCCESS != rval) continue;
          std::cout << ent_names[dim+1] << " "   << pgid << ", " 
                    << ent_names[dim] << " " << chgid << " sense is: ";
          if (1==sense) std::cout << "FORWARD" << std::endl;
          else std::cout << "REVERSE" << std::endl;
        }
      }
      else if (1 == dim) {
        sense_ents.clear();
        senses.clear();
        rval = gt.get_senses(*ch_it, sense_ents, senses);
        if (moab::MB_SUCCESS != rval) continue;
        for (unsigned int i = 0; i < sense_ents.size(); i++) {
          rval = mb->tag_get_data(gid_tag, &sense_ents[i], 1, &pgid);
          std::cout << ent_names[dim+1] << " "   << pgid << ", " 
                    << ent_names[dim] << " " << chgid << " sense is: ";
          if (-1 == senses[i]) std::cout << "REVERSED" << std::endl;
          else if (0 == senses[i]) std::cout << "BOTH" << std::endl;
          else if (1 == senses[i]) std::cout << "FORWARD" << std::endl;
          else std::cout << "(invalid)" << std::endl;
        }
      }
    }
  } 
}
